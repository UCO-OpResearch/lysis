====================================
Re-execution Environment (0.1.x)
====================================

*How the 0.1.x notebooks were made runnable again in 2026, and what breaks if
you deviate.*

The notebooks are not an accessory to this release line — they **are** the
interface. Every Dataset in :doc:`tpa_diffusion` and :doc:`internal_lysis` was
dispatched and processed from a ``F-Macro Multi-*-{Run,Process}.ipynb`` pair.
Re-running any of that analysis means first reconstructing a Python environment
old enough to load them, which is what this file records.

Reference environment
=====================

The target was the ``lysis-bak`` conda environment — the known-good stack the
0.1.x release was actually run against. It still exists at
``~/.conda/envs/lysis-bak`` ``[D]``, and every pin below was taken from it
rather than chosen. Conda is no longer used; the environment is now
``uv``-managed from ``pyproject.toml`` + ``uv.lock``.

**Python 3.10** (``requires-python = ">=3.10,<3.11"``). The scientific stack is
incompatible with numpy 2.x / pandas 2.x.

.. list-table::
   :header-rows: 1
   :widths: 26 16 58

   * - Package
     - Pin
     - Why this version
   * - ``numpy``
     - 1.24.1
     - ``lysis-bak``; 2.x breaks the model code
   * - ``pandas``
     - 1.5.3
     - ``lysis-bak``; 2.x breaks the analysis
   * - ``scipy``
     - 1.10.0
     - ``scipy.stats``, used by the analysis notebooks
   * - ``matplotlib``
     - 3.6.3
     - Every figure in both papers
   * - ``tqdm``
     - 4.64.1
     - Progress bars in the Run notebooks
   * - ``GooseSLURM``
     - 0.12.4
     - HPC job submission
   * - ``tabulate``
     - 0.9.0
     - ``DataFrame.to_markdown()``; see below
   * - ``jupyterlab``
     - 3.5.2
     - Release-era Lab; 4.x not tested here
   * - ``octave-kernel``
     - 0.35.1
     - See `Octave kernel`_
   * - ``ipywidgets``
     - 8.0.4
     - ``IProgress``, needed by ``tqdm.auto``
   * - ``matplotlib-inline``
     - 0.1.6
     - See `Transitive-dependency drift`_
   * - ``jupyter-black`` / ``black``
     - 0.3.3 / 22.12.0
     - In-notebook formatting

``matplotlib`` and ``GooseSLURM`` are **main** dependencies here, not optional
extras, because the notebooks are the primary way of driving the model on this
line. A bare ``uv sync`` provisions everything; there is no ``notebooks``
extra, unlike ``main``.

Transitive-dependency drift
===========================

The recurring failure mode on this line: ``uv`` resolves *indirect*
dependencies to their newest releases, which increasingly assume a modern
matplotlib/jupyter. Three surfaced one at a time ``[D]``:

.. list-table::
   :header-rows: 1
   :widths: 24 20 56

   * - Package
     - Floated to
     - Symptom
   * - ``black``
     - 26.x
     - Incompatible with the pinned stack
   * - ``matplotlib-inline``
     - 0.2.2
     - ``AttributeError: 'RcParams' object has no attribute '_get'`` on
       ``plt.figure()`` under the inline backend — 0.2.x calls
       ``matplotlib.rcParams._get()``, a private API absent in 3.6.3
   * - ``ipywidgets``
     - *(missing)*
     - ``TqdmWarning: IProgress not found`` from ``tqdm.auto``

Rather than keep patching reactively, the intersection of the resolved lock and
``lysis-bak`` was computed and **97 transitive dependencies pinned** in a
``[tool.uv] constraint-dependencies`` block, which caps versions without
forcing installs. This downgraded ipython 8.39 → 8.22.2, ipykernel 7.3.0 →
6.29.3, and dozens more ``[D]``.

Four packages are left deliberately floating, documented in the block's own
comment: ``metakernel``/``octave-kernel`` (absent from ``lysis-bak``),
``python-dateutil`` (bak's ``2.9.0`` is PyPI ``2.9.0.post0``) and ``tzdata``
(only modern pandas needs it).

.. warning::

   If a notebook starts failing in a way that looks like a matplotlib or
   jupyter API change, suspect an unconstrained indirect dependency before
   suspecting the notebook. This has been the cause every time so far.

Octave kernel
=============

Neither ``lysis`` nor ``lysis-bak`` still carries ``octave-kernel`` — it was
dropped at some point ``[D]``. The only release-blessed record is the original
v0.1.0 conda spec: ``octave_kernel 0.35.1`` with ``metakernel 0.29.4``.

Verified rather than assumed: 0.35.1 resolves cleanly against this env's exact
pins, and **the whole 0.35.1 → 0.39.0 range works**, because octave-kernel has
no numpy/pandas/scipy dependency at all — it shells out to the ``octave``
executable. It was confirmed end-to-end driving **Octave 10.1.0**, the module
the launch script loads ``[D]``. Pinned at 0.35.1 for fidelity; left unpinned
it pulls a modern ipykernel 7.x that conflicts with the rest of the stack.

Launching under Open OnDemand
=============================

``scripts/ood_jupyter_lab.sh`` prepares the environment. It **must be sourced,
not executed** — the ``PATH``/``VIRTUAL_ENV`` changes have to survive into the
shell that later runs ``jupyter-lab``. It deliberately does *not* start
JupyterLab: OOD appends its own ``jupyter-lab`` invocation with the flags its
websocket proxy needs.

Order matters. Modules load **first**, then ``uv sync``, then venv activation,
so ``.venv/bin`` ends up at the front of ``PATH`` and the bare ``jupyter-lab``
resolves to 3.5.2 rather than anything the modules bring.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Module
     - For
   * - ``Octave/10.1.0-foss-2023a``
     - Backs ``octave-kernel``
   * - ``intel-compilers/2023.1.0``
     - The Fortran binaries
   * - ``FFmpeg/6.0-GCCcore-12.3.0``
     - ``matplotlib``'s ``FFMpegWriter``, for the animation cells

.. note::

   **The ``uv: command not found`` trap.** ``uv`` is a user install, not a
   cluster module, so it reaches ``PATH`` via ``~/.bash_envars``. That is
   sourced from ``~/.bashrc`` *below* its interactive guard
   (``[ -z "$PS1" ] && return``), so in the non-interactive shell OOD uses,
   ``.bashrc`` returns before ever reaching it — and ``uv`` is missing ``[D]``.

   Fixed in the **personal** ``~/.bashrc`` by moving the ``.bash_envars``
   sourcing above the guard, *not* in the committed script. Hard-coding a
   ``uv`` path or sourcing a personal dotfile from a repo script would bake one
   user's install layout into everyone's launch path. The script requires only
   that ``uv`` is already on ``PATH``.

Fixes applied while re-executing
================================

Two gaps surfaced only once the notebooks actually ran, and both are now in the
committed environment:

* ``tabulate==0.9.0`` — ``DataFrame.to_markdown()`` is used in several cells of
  the ``2023-12-10-1900`` notebook and raised ``ImportError`` without it. Not
  present in ``lysis-bak``, so 0.9.0 was chosen as the release-era version
  rather than matched.
* ``FFmpeg`` module — the animation cells of the ``2023-02-02-2200`` notebook
  shell out to ``ffmpeg`` via ``FFMpegWriter``, and ``module purge`` leaves
  none on ``PATH``.

Release base
============

``dev-0.1.1`` is based on ``63a01bb``, **not** on the ``v0.1.0`` tag. The tag
points at ``368a6f2`` (2024-01-12) while the code used for the papers runs to
``63a01bb`` (2024-01-30); the two intervening commits are notebook-only. The
build commits were rebased onto ``63a01bb`` and force-pushed to correct this
for the re-release ``[D]``. See also :doc:`README`.

Known rough edges
=================

* **Restart the kernel and the Lab server after any ``uv sync``.** A freshly
  installed widget extension will not load into a running session, so e.g. the
  ``IProgress`` warning persists until restart even once ``ipywidgets`` is in.
* **Nothing timestamps the outputs.** No cell carries ``metadata.execution``,
  and the ``ExecuteTime`` entries in the ``2023-02-02-2200`` notebook are
  fossils from 2024-01-30/31 — that extension is a classic-notebook one and is
  inert under JupyterLab 3.5.2, so 2026 re-runs recorded nothing ``[D]``. To
  date future runs, enable **Notebook › Recording Timing** in JupyterLab, which
  writes real ``metadata.execution`` entries.
* **Execution counts are the only re-run fingerprint.** A fresh kernel
  renumbers ``In[]`` from 1, so comparing a notebook's execution counts and
  outputs against its committed version is what identifies which cells were
  actually re-executed.
