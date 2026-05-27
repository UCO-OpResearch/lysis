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
