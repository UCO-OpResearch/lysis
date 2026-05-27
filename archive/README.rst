Archive
=======

Unmaintained, preserved-for-reference code. **Nothing here is imported by the
``lysis`` package, packaged for distribution, or built into the documentation.**
It is kept in version control so past experiments can be revisited, but it is
not expected to run against the current codebase without rework.

GPU / CuPy macroscale experiment
--------------------------------

An early (circa 2022-2023) proof-of-concept that ran the NumPy macroscale model
on a GPU via `CuPy <https://cupy.dev/>`_. It demonstrated that the simulation
*could* use the GPU, but was never optimized and has since diverged
significantly from ``src/lysis/np_macroscale.py``. GPU acceleration is **not**
on the near- or medium-term roadmap; if it is revisited it will likely be
redesigned rather than resumed from this code.

Files:

- ``cp_macroscale.py`` --- CuPy port of the macroscale model
  (``CudaMacroscaleRun``). Was previously exposed as ``lysis.cp_macroscale`` via
  an optional ``import cupy`` guard in ``src/lysis/__init__.py``.
- ``cp_exec.py`` --- standalone runner script for the GPU model (references the
  even older ``CudaMacroscaleSim`` name).
- ``CuPy Workbook.ipynb`` --- benchmarking/exploration notebook driving
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
