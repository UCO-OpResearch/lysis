Archive
=======

Unmaintained, preserved-for-reference code. **Nothing here is imported by the
``lysis`` package, packaged for distribution, or built into the documentation.**
It is kept in version control so past experiments can be revisited, but it is
not expected to run against the current codebase without rework.

Subfolders group code by lineage: ``cupy/`` and ``cpp/`` are abandoned
alternative-language ports of the macroscale model, while ``fortran/`` holds
superseded prior versions of the *canonical* Fortran macroscale model --- the
lineage that the active ``src/fortran/macro_diffuse_into_and_along__*.f90``
replaced.

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
- ``macro_brad_scratch.f90`` --- the baseline into-and-along reference build,
  heavily annotated during the Python port; its behavior is what the active
  Fortran/Python macroscale model already implements.

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

   Excluded: ``macro_brad_scratch.f90`` --- its into-and-along behavior is the
   baseline Mechanism the active Fortran/Python model already implements.

   TODO: track this porting work as a GitHub issue once v2.0.0 ships.
