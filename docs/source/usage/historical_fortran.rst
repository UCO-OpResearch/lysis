===========================================
Running a Historical Fortran Implementation
===========================================

Both ``lysis run-micro`` and ``lysis run-macro`` accept a
``--fortran-commit <ref>`` flag that builds the Fortran binaries from an
arbitrary git ref before executing them.  The intended use cases are
*reproducing an earlier result* and *bisecting a behaviour change* — the
current Python wrapper drives a binary compiled from an older commit, and
the resulting HDF5 file is stamped so the deviation is auditable later.

Basic usage
-----------

::

    lysis run-micro data/run01.h5 \
        --executable micro_rates \
        --fortran-commit v0.3.1

::

    lysis run-macro data/my-experiment/ \
        --executable macro_diffuse_into_and_along__internal \
        --fortran-commit 60d200f \
        --slurm --compiler intel-compilers/2024

When ``--fortran-commit`` is set:

* ``--executable`` is interpreted as a binary **name** (with or without a
  leading ``bin/``) to pick from the historical build's ``bin/``
  directory.  The historical ``Makefile`` builds every target, so any
  binary that commit produced is selectable.

* A temporary build directory is created via ``tempfile.mkdtemp`` (in
  ``$TMPDIR`` or ``/tmp``), the source tree at that commit is extracted
  with ``git archive``, ``make`` is invoked, and the resulting binary is
  fed into the same :class:`~lysis.execution.fortran_micro.FortranMicro`
  / :class:`~lysis.execution.fortran_macro.FortranMacro` pipeline as a
  pre-built binary.

* The binary↔\ ``src/fortran/`` staleness check is bypassed.  The
  mismatch is intentional, so the usual warning banner is *not* prepended
  to the Fortran log.

* The HDF5 file is stamped with synthesised binary provenance:
  :data:`~lysis.config.constants.CONST.BINARY_COMMIT_ATTR` is the
  resolved SHA; :data:`~lysis.config.constants.CONST.BINARY_DIRTY_ATTR`
  is ``"clean"``; :data:`~lysis.config.constants.CONST.BINARY_COMPILER_ATTR`
  is the string from ``iso_fortran_env``'s ``compiler_version()`` intrinsic
  (probed by compiling a tiny F2008 stub with the same compiler the
  Makefile used); and
  :data:`~lysis.config.constants.CONST.BINARY_SOURCE_ATTR` is
  ``"historical:<sha>"``.

* The Python-side gates (``enforce_lysis_clean`` and
  ``enforce_init_commit_match``) still apply unchanged.  ``--fortran-commit``
  pins only the Fortran half of the stack; the Python wrapper is current
  code and must be coherent with the file's init stamp.

Lifecycle
---------

* The build runs once per CLI invocation, before the run loop.  Every
  HDF5 file the invocation processes shares the same binary.
* The build directory is removed when the CLI exits (success or failure),
  unless ``--keep-tmpdir`` is also set.
* In ``--slurm`` mode each submitted job receives its own copy of the
  binary in its staging directory (via the existing
  ``shutil.copy2(executable, staging_dir)`` step), so the master CLI's
  build directory is torn down immediately after the last ``sbatch``
  returns.  Parallel ``lysis run-macro`` invocations on different folders
  never collide — each gets a distinct ``mkdtemp``.

Compiler module
---------------

The ``--compiler`` flag (default ``intel-compilers/2023``) tells generated
Slurm scripts which LMod module to ``module load`` before invoking the
Fortran binary.  With ``--fortran-commit --slurm`` the build step is *also*
wrapped in ``module purge && module load <compiler>`` so the binary links
against the exact toolchain the Slurm preamble will load at runtime.

Local mode (no ``--slurm``) ignores ``--compiler``: the binary is built and
executed in the same environment the CLI was launched from, so anything you
loaded with ``module load`` before running ``lysis`` governs both steps.

Examples::

    # Local — uses whatever compiler is on PATH for both build and run.
    lysis run-micro data/run01.h5 \
        --executable micro_rates \
        --fortran-commit v0.3.1

    # Slurm — load Intel 2024 for both the build and the array tasks.
    lysis run-macro data/run01.h5 \
        --executable macro_diffuse_into_and_along__internal \
        --fortran-commit 60d200f \
        --slurm --compiler intel-compilers/2024

Limitations and caveats
-----------------------

CLI argument churn
~~~~~~~~~~~~~~~~~~

The Fortran binary's CLI has changed shape over time.  The current
wrapper passes flags like ``--runCode``, ``--outFileCode``, ``--radius``,
``--macro_seed``; older commits used different names (``--dist`` for
``--radius``, ``expCode`` for ``--runCode``, ``runs`` for
``simulations``, …) and commits before ~2024-01 had no CLI at all
(hardcoded source parameters + Matlab input files).

For the initial cut of ``--fortran-commit`` we assume the chosen commit
uses the **current** CLI.  Picking a ref that predates the current flag
names will fail with ``make`` succeeding but the binary rejecting its
arguments at runtime.  A follow-up will add a small per-SHA-range
compatibility registry to ``src/lysis/execution/fortran_cli_compat.py``
that rewrites flag names at the wrapper level; until then, use
``git log -- src/fortran`` to confirm the chosen commit is recent enough
that its CLI matches what the wrapper expects.

Required source paths
~~~~~~~~~~~~~~~~~~~~~

``git archive`` extracts ``src/fortran``, ``src/c``, and ``Makefile``
from the chosen commit.  ``src/c`` carries ``kiss.c``, the C random
number generator linked into every binary.  Refs that predate any of
these paths will produce a build error; the tail of ``.build_log`` is
shown in the resulting ``click.ClickException``.

Debugging
~~~~~~~~~

Pass ``--keep-tmpdir`` to preserve the build directory after the CLI
exits.  Its path is announced when the build starts; inside it you'll
find ``.build_log`` (full ``make`` output), the extracted ``src/`` tree,
and the ``bin/`` directory with the produced binary.
