==================================
Running the Microscale Simulation
==================================

Overview
--------

The ``lysis run-micro`` command (and its Python API counterpart
:class:`~lysis.execution.codeutil.FortranMicro`) executes the compiled
Fortran microscale binary for a Run (see :doc:`ontology`), then automatically
imports the results back into the same HDF5 file.

The full workflow is:

1. Read :class:`~lysis.config.parameters.MicroParameters` from an existing
   ``.h5`` file.
2. Create a temporary working directory alongside the HDF5 file.
3. Execute the Fortran binary, writing output to that directory.
4. Convert the Fortran output to the HDF5 v2.0.0 format and write it into
   the source ``.h5`` file.
5. Remove the temporary directory (unless ``--keep-tmpdir`` is set).

.. note::

   ``run-micro`` expects an HDF5 file that already contains
   ``micro_params``.  If you are starting from scratch, use
   ``lysis init-experiment`` first (see :doc:`experiment_init`).


Prerequisites
-------------

Compiled binary
~~~~~~~~~~~~~~~

The Fortran microscale binary (``micro_rates``) must be compiled before
running any simulations.  The source is at ``src/fortran/micro_rates.f90``.
A typical one-time build on the UCO HPC cluster::

    module load gcc
    gfortran -O2 -o bin/micro_rates src/fortran/micro_rates.f90

Place the resulting binary anywhere accessible; the ``--executable`` option
accepts an absolute or relative path.

HDF5 run file
~~~~~~~~~~~~~

Each Run must have an ``.h5`` file containing ``micro_params`` before
``run-micro`` can be called.  The recommended way to create these files is
``lysis init-experiment`` (see :doc:`experiment_init`).  They can also be
created programmatically — see `Python API`_ below.


CLI Usage
---------

.. code-block:: text

    lysis run-micro [OPTIONS] HDF5_PATH

``HDF5_PATH``
    Path to the ``.h5`` file for the Run.  Must contain ``micro_params``.
    Simulation results are written back into this same file on completion.

Options
~~~~~~~

.. option:: --executable PATH

    **Required.** Path to the compiled Fortran microscale binary.

.. option:: --slurm

    Dispatch via Slurm instead of executing locally.  Submits a master job
    that orchestrates the simulation and imports results on completion.
    Returns immediately after printing the Slurm job ID.

.. option:: --partition NAME

    Slurm partition to target (e.g. ``normal``, ``long``).  Only meaningful
    with ``--slurm``.

.. option:: --staging-root PATH

    Root directory for the shared staging temporary directory used by the
    Slurm job.  Defaults to the parent directory of ``HDF5_PATH``.  Only
    meaningful with ``--slurm``.

    .. note::

       The staging directory is created with a unique name under
       ``--staging-root`` using ``mktemp``-style naming.  It must be on a
       shared filesystem visible to all compute nodes.  Do **not** point this
       at a node-local path.

.. option:: --fast-tmp-root PATH

    Root directory for fast node-local scratch storage (e.g. an NVMe drive
    mounted at ``/nvme/scratch``).  When set, the Fortran binary writes to
    this location and the output is moved to the staging directory before the
    job finishes (two-tier storage).  Only meaningful with ``--slurm``.

    Leave this unset if your cluster does not expose per-node scratch storage.

.. option:: --keep-tmpdir

    Always preserve the temporary output directory after the run finishes.
    Normally the directory is deleted on success (on failure it is already
    preserved by default).  Use this flag when you need to inspect the raw
    Fortran output files for debugging.

.. option:: --file-code TEXT

    Output file code suffix appended to all Fortran output filenames
    (e.g. ``_PLG2_tPA01``).  Defaults to the empty string.  Use this when
    multiple runs write to the same directory and their output files would
    otherwise collide.

Examples
~~~~~~~~

Local execution (simplest case):

.. code-block:: bash

    lysis run-micro data/my-experiment/run-01.h5 --executable bin/micro_rates

Preserve temporary directory for debugging:

.. code-block:: bash

    lysis run-micro data/my-experiment/run-01.h5 \
        --executable bin/micro_rates \
        --keep-tmpdir

Submit via Slurm (default partition):

.. code-block:: bash

    lysis run-micro data/my-experiment/run-01.h5 \
        --executable bin/micro_rates \
        --slurm

Submit via Slurm on a specific partition with a shared staging directory:

.. code-block:: bash

    lysis run-micro data/my-experiment/run-01.h5 \
        --executable bin/micro_rates \
        --slurm \
        --partition normal \
        --staging-root /work/mygroup/staging

Submit via Slurm using fast node-local storage:

.. code-block:: bash

    lysis run-micro data/my-experiment/run-01.h5 \
        --executable bin/micro_rates \
        --slurm \
        --partition normal \
        --fast-tmp-root /nvme/scratch


Python API
----------

The same functionality is available directly from Python via
:class:`~lysis.execution.codeutil.FortranMicro`.

Simple one-call workflow
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from lysis.execution.codeutil import FortranMicro

    fm = FortranMicro.from_hdf5("data/my-experiment/run-01.h5",
                                executable="bin/micro_rates")
    fm.run_full("data/my-experiment/run-01.h5")

:meth:`~lysis.execution.codeutil.FortranMicro.run_full` creates a temporary
directory, executes the binary, imports results into the HDF5 file, and
removes the temporary directory.

Step-by-step workflow
~~~~~~~~~~~~~~~~~~~~~

Use this when you need access to intermediate products (e.g. the raw Fortran
output files) or want to handle execution and import separately:

.. code-block:: python

    from pathlib import Path
    from lysis.execution.codeutil import FortranMicro

    hdf5_path = Path("data/my-experiment/run-01.h5")
    work_dir  = Path("data/my-experiment/run-01-workdir")
    work_dir.mkdir(exist_ok=True)

    fm = FortranMicro.from_hdf5(hdf5_path, executable="bin/micro_rates")

    # Execute the binary; returns path to data/{run_code}/ inside work_dir
    data_dir = fm.exec_in_workdir(work_dir)

    # Import results into the HDF5 file
    FortranMicro.import_results(data_dir, hdf5_path, keep_tmpdir=True)

Creating a Run programmatically
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have not used ``init-experiment``, you can create a Run and its HDF5
file directly:

.. code-block:: python

    import h5py
    from lysis.config.constants import CONST, Q_
    from lysis.config.parameters import MicroParameters
    from lysis.config.run import Run
    from lysis.execution.codeutil import FortranMicro

    # Build parameters
    mp = MicroParameters(
        fiber_radius=Q_("61.5 nanometer"),
        nodes_in_micro_row=13,
        micro_simulations=50_000,
    )

    # Write a minimal HDF5 file
    hdf5_path = "data/my-run/run-01.h5"
    with h5py.File(hdf5_path, "w") as f:
        f.attrs[CONST.DATASPEC_VERSION_ATTR] = "v2.0.0"
        grp = f.require_group("micro_data")
        for k, v in mp.to_basedict().items():
            grp.attrs[k] = str(v) if not isinstance(v, (int, float, bool)) else v

    # Execute and import
    fm = FortranMicro.from_hdf5(hdf5_path, "bin/micro_rates")
    fm.run_full(hdf5_path)


Typical Workflow
----------------

The intended end-to-end workflow for a new experiment is:

1. **Define runs** — prepare a CSV file with the parameter values for each
   Run (see :doc:`experiment_init`).

2. **Initialise the experiment** — create the HDF5 files:

   .. code-block:: bash

       lysis init-experiment runs.csv /data/experiments/ --name fiber-sweep

3. **Execute each Run** — run the microscale simulation for each ``.h5``
   file.  On a cluster this is typically done for all runs in a loop:

   .. code-block:: bash

       for h5 in /data/experiments/fiber-sweep/*.h5; do
           lysis run-micro "$h5" \
               --executable bin/micro_rates \
               --slurm \
               --partition normal
       done

4. **Analyse results** — once the jobs complete, the HDF5 files each contain
   the microscale output datasets (``micro_data/sim_final_time``,
   ``micro_data/fiber_degraded``, etc.) ready for analysis.


Troubleshooting
---------------

Binary exits immediately without writing output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Fortran binary writes data to ``data/{run_code}/`` **relative to its
working directory**.  If that subdirectory does not exist the binary may fail
silently.  :meth:`~lysis.execution.codeutil.FortranMicro.exec_in_workdir`
creates it automatically; if you are invoking the binary manually, create the
directory first.

Inspect raw Fortran output
~~~~~~~~~~~~~~~~~~~~~~~~~~

Pass ``--keep-tmpdir`` to preserve the working directory.  Inside you will
find:

* ``data/{run_code}/lysis{file_code}.dat`` — binary lysis-time data
* ``data/{run_code}/micro{file_code}.txt`` — full Fortran log (stdout)
* ``data/{run_code}/params.json`` — parameter snapshot used for import

Results not appearing in the HDF5 file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The import step reads from the temporary ``data/{run_code}/`` directory.
If the binary fails (non-zero exit code), the directory is preserved for
inspection even without ``--keep-tmpdir``.  Check
``micro{file_code}.txt`` for error messages from the Fortran code.

Slurm job submitted but results never imported
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The master Slurm job is a Python script that polls ``squeue`` until the
child job finishes before importing.  If the master job itself is killed (e.g.
wall-time exceeded) the staging directory is preserved.  Re-run the import
step manually:

.. code-block:: python

    from pathlib import Path
    from lysis.execution.codeutil import FortranMicro

    # Path printed when the child job was submitted, or found in staging-root
    data_dir  = Path("/work/mygroup/staging/lysis-micro-run-01-XXXXX/data/run-01")
    hdf5_path = Path("data/my-experiment/run-01.h5")

    FortranMicro.import_results(data_dir, hdf5_path, keep_tmpdir=True)
