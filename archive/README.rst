Archive
=======

Unmaintained, preserved-for-reference code. **Nothing here is imported by the
``lysis`` package, packaged for distribution, or built into the documentation.**
It is kept in version control so past experiments can be revisited, but it is
not expected to run against the current codebase without rework.

Subfolders group code by lineage: ``cupy/``, ``cpp/`` and ``c/`` are abandoned
alternative-language ports of the macroscale model; ``fortran/`` holds
superseded prior versions of the *canonical* Fortran macroscale model (the
lineage the active ``src/fortran/macro_diffuse_into_and_along__*.f90``
replaced); and ``matlab/`` holds the original MATLAB/Octave pre- and
post-processing pipeline, now superseded by the ``lysis`` Python package; and
``scripts/`` holds the pre-CLI run scripts (plus one unused package module),
superseded by the ``lysis`` command-line interface.

GPU / CuPy macroscale experiment (``cupy/``)
--------------------------------------------

An early (circa 2022-2023) proof-of-concept that ran the NumPy macroscale model
on a GPU via `CuPy <https://cupy.dev/>`_. It demonstrated that the simulation
*could* use the GPU, but was never optimized and has since diverged
significantly from ``src/lysis/np_macroscale.py``. GPU acceleration is **not**
on the near- or medium-term roadmap; if it is revisited it will likely be
redesigned rather than resumed from this code.

Files:

- ``cupy/cp_macroscale.py`` --- CuPy port of the macroscale model
  (``CudaMacroscaleRun``). Was previously exposed as ``lysis.cp_macroscale`` via
  an optional ``import cupy`` guard in ``src/lysis/__init__.py``.
- ``cupy/cp_exec.py`` --- standalone runner script for the GPU model (references
  the even older ``CudaMacroscaleSim`` name).
- ``cupy/CuPy Workbook.ipynb`` --- benchmarking/exploration notebook driving
  ``lysis.CudaMacroscaleRun``.

These require the optional GPU dependencies (``cupy``, ``nvtx``), still declared
under the ``gpu`` extra in ``pyproject.toml``.

C++ macroscale port (``cpp/``)
------------------------------

A partial conversion of the macroscale model to C++, half-completed in 2015 and
never finished. Formerly ``src/cpp/``; nothing in the active build references it.
It was the only consumer of the old ``cpp`` Makefile target (now removed).

Note on KISS RNG: the live build's shared random number generator is compiled
from ``src/c/kiss.c`` (into ``bin/kiss.o`` for the Fortran binaries and
``lib/kiss.so`` for ``np_macroscale`` via ``lysis.tools.kiss``). The
``kiss.c`` / ``kiss.h`` / ``kiss.o`` copies preserved here belong to the C++
port and were *not* part of that shared build --- the old ``cpp`` target itself
linked ``src/c``'s ``kiss.o`` and used only ``cpp/kiss.h`` as a header.

Files: ``macro_Q2.cpp`` (the port), ``fileio.h``, and the port's own
``kiss.c`` / ``kiss.h`` / ``kiss.o``.

C / OpenMPI macroscale port (``c/``)
------------------------------------

An abandoned partial conversion of the macroscale model to C, parallelised
across the lattice with OpenMPI/MPICH (built into a ``c_macro`` /
``ParallelMacro`` binary). Written by **Bryan Carroll** (misspelled "Carrol"
in the 2022-12-21 import commit) and imported into this repository in Dec 2022
-- back when it was named ``UCO-OpResearch/BloodClotting``, before the rename
to ``lysis``. His own commits are not in this history (the code was imported
wholesale by Brad Paynter), so this note is the record of his authorship. Formerly ``src/c/`` (minus the kiss files); nothing
in the active build references it (the old ``c`` / ``c-macro`` Makefile target
has been removed). ``old/`` holds earlier drafts, including a C++ variant.

KISS RNG note: the canonical ``kiss.c`` / ``kiss.h`` -- the shared RNG compiled
into ``bin/kiss.o`` (Fortran binaries) and ``lib/kiss.so`` (``np_macroscale``
via ``lysis.tools.kiss``) -- **remain in** ``src/c/`` and were deliberately not
moved. This C port ``#include``\ s ``kiss.h`` (via ``all.h``) and linked
``kiss.c``, so rebuilding it from the archive would need those paths adjusted.

Fortran macroscale variants (``fortran/``)
------------------------------------------

Superseded versions of the canonical Fortran macroscale model (``program
macrolysis``), each encoding a different tPA forced-unbinding / diffusion
hypothesis. They worked and produced results but were replaced by the active
``src/fortran/macro_diffuse_into_and_along__{internal,external}.f90``. Formerly
``src/fortran/_Archive/``.

Files:

- ``macro_Q2_diffuse_into.f90`` --- forced-unbind tPA is removed from binding
  but may diffuse *into* the clot (assumed attached to a small FDP). Formerly
  ``macro_Q2_forcedtPArebind.f90``.
- ``macro_Q2_diffuse_along.f90`` --- forced-unbind tPA may diffuse only *away
  from or along* the clot front, never into it (FDPs assumed too large).
- ``macro_Q2_always_rebind.f90`` --- forced-unbind tPA is allowed to immediately
  rebind.
- ``macro_diffuse_into_and_along_slow_micro__external.f90`` --- the
  into-and-along behavior plus the ability to slow micro-unbound tPA movement by
  an integer factor.
- ``macro_rng_array.f90`` --- the into-and-along behavior, but draws every
  random number at the start of each iteration rather than on demand (an
  RNG-ordering rule that matches how the Python model draws random numbers).

.. note::

   In the project's ontology a **Mechanism** is any set of rules governing how
   the model executes --- it need not mimic a biological process --- so each of
   these files encodes a Mechanism. Once **v2.0.0** is released, the Mechanisms
   not yet available in the Python model should be reimplemented as Python
   Mechanisms so they can be run from the current codebase:

   - "diffuse into" --- ``macro_Q2_diffuse_into.f90``
   - "diffuse along" --- ``macro_Q2_diffuse_along.f90``
   - "always rebind" --- ``macro_Q2_always_rebind.f90``
   - "slow micro movement" --- ``macro_diffuse_into_and_along_slow_micro__external.f90``
   - "all RNG up front" --- ``macro_rng_array.f90`` (draws every random number at
     the start of each iteration; the Python model already draws RNG this way,
     so the work is mainly to expose it as a selectable Mechanism rather than
     reimplement the behavior)

   The plain into-and-along baseline that these variants descend from is the
   behavior the active ``__internal`` / ``__external`` model already implements,
   so it is not listed above.

   TODO: track this porting work as a GitHub issue once v2.0.0 ships.

MATLAB / Octave pipeline (``matlab/``)
--------------------------------------

The original pre- and post-processing pipeline for the Fortran model, written
in MATLAB by Dr. Brittany Bannish, together with Octave/Jupyter notebook ports
of the same code (the ``.ipynb`` files) made by Brad Paynter for easier use
under the Octave kernel. All of this functionality is now provided by the
``lysis`` Python package. Formerly ``src/matlab/``, plus ``octave_setup.sh``
(from ``scripts/``) and ``octave_notebooks.rst`` (from ``docs/source/usage/``;
removed from the ``docs/source/index.rst`` toctree when archived).

Pipeline (each ``.m`` has a matching ``.ipynb`` Octave port):

- ``Lat_create.m`` --- generates the lattice input file read by the Fortran
  microscale code (as ``Lat``).
- ``micro_to_macro.m`` --- reads the microscale results, produces figures and
  metrics, and writes the input files required by the ``macro_Q2`` macroscale
  code.
- ``macro_scale_model_post_processing.m`` --- post-processes macroscale results
  (metrics and plots).
- ``movie.m`` / ``movie_snapshot.m`` --- animate / snapshot macroscale results.

Supporting files:

- ``octave_setup.sh`` --- one-time environment setup for running the notebooks
  under Octave on the Buddy cluster.
- ``octave_notebooks.rst`` --- the former usage guide for running these
  notebooks on Buddy.

Legacy run scripts + dead package code (``scripts/``)
-----------------------------------------------------

The pre-CLI scripts that drove the model before the ``lysis`` command-line
interface (``lysis run-micro`` / ``lysis run-macro``) replaced them, plus one
unused package module.

Run scripts (formerly ``scripts/`` and the repo-root ``exec.sh``):

- ``exec.py`` --- Python driver that ran a macroscale Run and recorded results
  through the ``DataStore``.
- ``fortran_exec.py`` --- wrapper that launched the compiled Fortran
  micro/macro executables (via ``lysis.execution.codeutil``).
- ``micro_to_macro.py`` --- Python reimplementation of the MATLAB
  ``micro_to_macro`` bridge (computed ``frac_forced`` and wrote the macroscale
  input files).
- ``exec.sh`` --- shell wrapper that invoked ``fortran_exec.py``.
- ``micro_fortran_run.sh`` / ``macro_fortran_run.sh`` --- SLURM batch scripts
  for running the Fortran micro/macro models on the cluster.

Dead package module (formerly ``src/lysis/molecule.py``):

- ``molecule.py`` --- an unused ``Molecule`` dataclass. It was re-exported by
  ``lysis/__init__.py`` but never instantiated anywhere; ``np_macroscale``
  represents molecules with arrays instead.

Note: the usage guides ``docs/source/usage/fortran_microscale.rst`` and
``fortran_macroscale.rst`` still walk through the archived shell scripts and
``micro_to_macro.py``; they should be updated to the CLI workflow.
