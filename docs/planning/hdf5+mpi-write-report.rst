==============================================================================
Parallel HDF5 from Python on buddy-2.6: working recipes for B and C
==============================================================================

This report shows how to write to a shared HDF5 file from multiple MPI ranks
(on one node or across many nodes) when the file already contains data — both
the "overwrite existing tables" case (B) and the "leave a reference subtree
alone, add new results next to it" case (C).

Tested 2026-05-29 on UCO **buddy-2.6** with:

- OpenMPI 4.1.5-GCC-12.3.0
- HDF5 1.14.0-gompi-2023a (built with parallel MPI-IO support)
- Python 3.11, h5py 3.16.0 (built from source with ``HDF5_MPI=ON``),
  mpi4py 4.1.2
- Home directory on NFS
  (``storage-01.hpc.uco.edu:/mnt/shares/users/...``)

Verified with 10 ranks both on a single node *and* across 10 different nodes.

TL;DR — the one rule
====================

    **Every metadata operation on a file opened with the** ``mpio`` **driver
    must be collective: all ranks must execute the same** ``create_group``\ **,**
    ``create_dataset``\ **, and** ``del`` **calls in the same order. Only the
    actual data writes can be per-rank.**

This is the rule the HDF5 docs state, but in subtle ways the consequences
are easy to miss. The trap we fell into: pre-creating per-rank datasets
serially on the login node, then opening the file under ``mpio`` and having
each rank only fill its own dataset. The file looks right serially, but
under parallel access HDF5's metadata cache is uncoordinated and writes
end up routed to the wrong byte ranges (we saw rank N's host string
appearing in rank M's slot, or all empty except one). It silently
"succeeds." It is not an NFS problem — pre-created datasets fail on a
single node too.

Recipe B: overwrite a pre-populated file
========================================

**Use when:** the file was produced by an earlier (serial or otherwise)
process, and you want to refill the existing tables in parallel.

**Plan:** open under ``mpio``, collectively delete every dataset that will
be rewritten, collectively recreate them, then have each rank fill its own.

.. code-block:: python

    # prepare.py — serial, runs anywhere
    import h5py, numpy as np
    with h5py.File("out.h5", "w") as f:
        f.attrs["n_ranks"] = N_RANKS
        f.attrs["n_rows"] = N_ROWS
        for r in range(N_RANKS):
            d = f.create_dataset(f"rank_{r}/data", shape=(N_ROWS,), dtype="f8")
            d[:] = SOME_INITIAL_DATA  # marker / draft / previous run
            h = f.create_dataset(f"rank_{r}/host", shape=(1,), dtype="S64")
            h[0] = b"PREP-MARKER"

.. code-block:: python

    # worker.py — every rank runs this under mpirun
    from mpi4py import MPI
    import h5py, numpy as np, socket

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    with h5py.File("out.h5", "a", driver="mpio", comm=comm) as f:
        n_ranks = int(f.attrs["n_ranks"])
        n_rows  = int(f.attrs["n_rows"])
        assert size == n_ranks

        # COLLECTIVE: delete every old dataset, then recreate.
        # Every rank must execute these calls in the same order.
        for r in range(n_ranks):
            del f[f"rank_{r}/data"]
            del f[f"rank_{r}/host"]
        for r in range(n_ranks):
            f.create_dataset(f"rank_{r}/data", shape=(n_rows,), dtype="f8")
            f.create_dataset(f"rank_{r}/host", shape=(1,), dtype="S64")

        # INDEPENDENT: each rank only touches its own slot.
        payload = np.full(n_rows, float(rank)) + np.arange(n_rows) * 1e-3
        f[f"rank_{rank}/data"][:] = payload
        f[f"rank_{rank}/host"][0] = f"{socket.gethostname()}:rank={rank}".encode()

This is the ``b`` scenario in ``worker.py``. Verified passing single-node
(10 ranks on node-302) and multi-node (10 ranks across node-104..113)
with all 10 datasets containing the new payload.

Variant: in-place update without re-creating
--------------------------------------------

You can edit the contents of an existing dataset under ``mpio`` *as long
as you do not touch any metadata*. That means: same shape, same dtype, no
attribute writes, no dataset creation. In that case skip the
``del`` / ``create_dataset`` block in ``worker.py`` and just do the
independent writes. We did not exhaustively verify this variant — if your
existing file came from a serial process and you do not control its byte
layout, prefer the explicit delete-and-recreate above.

Recipe C: add new parallel data without disturbing existing data
================================================================

**Use when:** the file holds reference / configuration / earlier-run data
that the workers should read but never modify, and the workers append new
results next to it.

**Plan:** prepare the file serially with two top-level groups: a populated
``reference`` subtree and an empty ``worker`` subtree. The workers open
the file under ``mpio``, leave ``reference`` strictly alone, and create
their datasets inside ``worker`` collectively. Then each rank writes its
own slot.

.. code-block:: python

    # prepare.py — serial
    import h5py, numpy as np
    with h5py.File("out.h5", "w") as f:
        f.attrs["n_ranks"] = N_RANKS
        f.attrs["n_rows"] = N_ROWS

        ref = f.create_group("reference")
        for r in range(N_RANKS):
            ref.create_dataset(f"rank_{r}/data", data=REFERENCE_PAYLOAD[r])
            ref.create_dataset(f"rank_{r}/host", data=np.array([b"REFERENCE"]))

        f.create_group("worker")  # empty; the workers will fill it

.. code-block:: python

    # worker.py
    from mpi4py import MPI
    import h5py, numpy as np, socket

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    with h5py.File("out.h5", "a", driver="mpio", comm=comm) as f:
        n_ranks = int(f.attrs["n_ranks"])
        n_rows  = int(f.attrs["n_rows"])
        assert size == n_ranks

        wg = f["worker"]
        # COLLECTIVE: every rank creates every worker dataset in the same order.
        for r in range(n_ranks):
            wg.create_dataset(f"rank_{r}/data", shape=(n_rows,), dtype="f8")
            wg.create_dataset(f"rank_{r}/host", shape=(1,), dtype="S64")

        # INDEPENDENT: each rank writes only its own slot.
        payload = np.full(n_rows, float(rank)) + np.arange(n_rows) * 1e-3
        f[f"worker/rank_{rank}/data"][:] = payload
        f[f"worker/rank_{rank}/host"][0] = (
            f"{socket.gethostname()}:rank={rank}".encode())

This is the ``c`` scenario. Verified passing single-node and multi-node:
all ``worker/rank_R`` slots contain the new payload, and *every*
``reference/rank_R`` dataset is bit-identical to what ``prepare.py``
wrote. There is no special locking, no exclusion, no read-only flag
needed — the rule is simply "never touch it, including no metadata
operations".

You can extend this pattern to: multiple reference subtrees,
inputs/outputs separation, run-N output groups in the same file, etc.

Launch pattern (SLURM)
======================

Two cluster-specific gotchas affect every multi-rank Python launch on
buddy-2.6 — they have nothing to do with PHDF5 but you will hit them if
you forget:

1. **Bare** ``python`` **inside an interactive** ``salloc`` / ``srun``
   **allocation fails** at ``MPI_Init_thread``, because the OpenMPI build
   looks for SLURM PMI components it cannot find. Always wrap it:
   ``mpirun -np 1 python ...`` works in interactive mode and
   ``mpirun -np N python ...`` works inside an ``sbatch`` script.

2. ``srun --mpi=pmix`` **hangs** during PMIx wireup (the OpenMPI build is
   missing the ``psec/munge`` component, and even with
   ``PMIX_MCA_psec=native`` the bootstrap stalls cross-node). Use
   ``mpirun -np N python worker.py`` inside the sbatch — it handles
   cross-node spawn through the OpenMPI ``PLM:slurm`` component and just
   works.

A working multi-node template (10 ranks across 10 nodes; for single-node
use ``--nodes=1 --ntasks-per-node=10`` instead):

.. code-block:: bash

    #!/bin/bash -l
    #SBATCH --job-name=ph5test
    #SBATCH --nodes=10
    #SBATCH --ntasks=10
    #SBATCH --ntasks-per-node=1
    #SBATCH --cpus-per-task=1
    #SBATCH --time=00:05:00
    #SBATCH --mem-per-cpu=1G
    #SBATCH --partition=general
    #SBATCH --output=slurm-%j.out

    set -euo pipefail

    # Same toolchain the python extensions were built against.
    module load OpenMPI/4.1.5-GCC-12.3.0
    module load HDF5/1.14.0-gompi-2023a

    cd /path/to/project
    source .venv/bin/activate

    # Bypass PMIx auth components the bundled OpenMPI doesn't ship.
    export PMIX_MCA_psec=native
    export PMIX_MCA_gds=hash

    # Prepare can run on the login node; here we re-do it on the compute node
    # to keep one self-contained script.
    python prepare.py

    mpirun -np 10 python worker.py

``GooseSLURM`` is fine for templating this string but its
``scripts.plain(...)`` helper emits ``--job_name`` (underscore) where
sbatch needs ``--job-name`` (hyphen), so prefer either a literal heredoc
or format string like the above, or pass the SBATCH flags through
``sbatch <args> script.sh``.

Pitfalls
========

- **Pre-created datasets** are a trap. If ``prepare.py`` calls
  ``create_dataset`` and the worker opens the file under ``mpio`` without
  recreating them, writes get cross-routed silently (this looked like an
  NFS coherency bug for a while; it isn't). Either delete and recreate
  (recipe B) or only pre-create groups (recipe C).
- **Variable-length string datasets** (``dtype=h5py.string_dtype()``)
  cannot be written under PHDF5 ("Parallel IO does not support writing VL
  or region reference datatypes yet"). Use a fixed-width byte string like
  ``"S64"`` and ``bytes`` payloads.
- **Compression** (``compression="gzip"`` or any chunked filter pipeline)
  breaks the independent-write pattern: PHDF5 raises *"Can't perform
  independent write with filters in pipeline."* You would need to switch
  to collective writes, which works in principle (h5py
  ``with dset.collective:``) but is awkward when each rank has its own
  dataset — every rank has to participate in every dataset's write,
  selecting empty slabs for the ranks that don't own that dataset, and
  we did not get a working version of that on this cluster. For now,
  leave compression off on PHDF5 outputs.
- **OpenMPI's** ``ompio`` **MPI-IO backend** is the default and behaves
  inconsistently on this NFS. We left it on for the tests above and they
  passed, but if you see strange corruption that recipe B/C can't
  explain, try ``export OMPI_MCA_io=romio321`` to force ROMIO.

Build recipe (one-time, slow)
=============================

The wheels on PyPI ship h5py linked against a serial bundled HDF5, so
parallel I/O is disabled. Build from source against the cluster's
parallel HDF5:

.. code-block:: bash

    module load OpenMPI/4.1.5-GCC-12.3.0
    module load HDF5/1.14.0-gompi-2023a
    uv venv --python 3.11 .venv
    source .venv/bin/activate
    MPICC=mpicc CC=mpicc uv pip install --no-binary mpi4py mpi4py
    HDF5_MPI=ON HDF5_DIR="$EBROOTHDF5" CC=mpicc uv pip install --no-binary h5py h5py
    uv pip install GooseSLURM   # pure python, wheel is fine

Verify with ``python -c "import h5py; print(h5py.get_config().mpi)"`` —
must print ``True``. The h5py build can take an hour on the login node;
do it once.

Files in this directory
=======================

- ``prepare.py``, ``worker.py``, ``verify.py`` — parameterized by
  ``--scenario {a,b,c,d}``; the report's B and C recipes are the ``b``
  and ``c`` cases.
- ``submit.py`` — builds a one-shot sbatch script for a given
  ``--scenario X --mode {single,multi}`` and submits it. Single-node
  uses one node with N tasks; multi-node uses N nodes with one task each.
- ``results_*.h5`` — the verified output files for each scenario/mode
  pair (B and C, both modes — they all pass).
- ``slurm_*_*.out`` — the SLURM job stdout/stderr for each run, useful
  as reference output to share with cluster maintainers if something
  regresses on the next OS / toolchain upgrade.
