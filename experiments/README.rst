================================
Experiment Logs
================================

*A retrospective record of the computational Experiments run with this model.*

Each file in this directory logs one :term:`Experiment` — a collection of
:term:`Runs <Run>` compared to elucidate cause and effect. The log records, for
each Dataset in the Experiment: its code, the notebooks that executed and
processed it, where it sits in the lineage, and **why** it was superseded.

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Log
     - Span
     - Subject
   * - :doc:`tpa_diffusion`
     - 2023-02 → 2023-11
     - tPA transport through the clot under varying unbinding Mechanisms
   * - :doc:`internal_lysis`
     - 2023-04 → 2023-12
     - Lysis initiated from molecules starting *inside* the clot
   * - :doc:`reexecution_environment`
     - 2026-08
     - How to run any of the above again, and what breaks if you deviate

.. note::

   **Tidy Fortran Data** (``2024-02-02-1400``) was drafted here and then moved
   to the **v0.2.x** line, where its code and notebooks live. It derives from
   the last Internal Lysis Dataset, so it appears under "Downstream" in
   :doc:`internal_lysis` — but the log itself belongs with v0.2.x.

Provenance of these logs
========================

These logs were **reconstructed in August 2026**, long after the work was done.
They are archaeology, not a contemporaneous lab notebook. Every entry is marked
with the evidence it rests on:

``[S]``
   From ``Computational Result Data.xlsx``, the hand-maintained index of
   Dataset codes. This is a *primary* record — written at the time by the
   researcher — and is the only surviving source for most "reason for change"
   entries.

``[N]``
   Read directly out of the notebook source (``runs``/``scenarios``/
   ``mechanisms`` arrays, ``group_code``). Highly reliable.

``[D]``
   Confirmed present on disk as a data directory in one of the live checkouts.

``[A]``
   Confirmed present as a ``.tar.gz`` in ``~/Archive/lysis_data/``.

``[O]``
   Confirmed present in the OneDrive archive, from a directory listing taken
   2026-08-21. Seen in a listing, not opened.

``[?]``
   Inferred. Treat as a hypothesis, not a record.

Where the data lives
====================

Data is **not committed to this repository** — ``data`` is in ``.gitignore``.
These logs are therefore deliberately kept outside the data tree, so that they
can be version-controlled. A log entry names its Dataset; it never contains it.

Raw data is spread across at least four places, and no single one is complete:

.. list-table::
   :header-rows: 1
   :widths: 30 18 52

   * - Location
     - Status
     - Notes
   * - ``~/git/UCO-OpResearch/lysis/data/``
     - **canonical**
     - The reference copy. Prefer this when a Dataset exists in more than one
       place.
   * - ``~/git/UCO-OpResearch/lysis-v0.1.0/data/``
     - working copy
     - Selected Datasets copied into this worktree so the 0.1.x analysis could
       be re-executed without overwriting the canonical set.
   * - ``~/Archive/lysis_data/``
     - cold archive
     - 173 per-Run ``.tar.gz`` (a few 2024 sets are ``.tar.xz``). Holds several
       Datasets that exist nowhere else.
   * - OneDrive — ``Lysis/Archive``
     - offsite
     - The deepest store. Holds every 2023 Internal Lysis Dataset, including
       two that exist nowhere else. Mostly ``.tar.gz``; some sets also
       unpacked.
   * - OneDrive — ``Lysis/Current``
     - offsite
     - The 2024-era working set (``2024-01-26`` onward) plus unpacked copies of
       the headline 2023 Datasets.

.. note::

   OneDrive **has** now been surveyed, from ``ls`` listings of both folders
   taken 2026-08-21. With those included, **every Dataset in both Experiments
   logged here is accounted for.** No Dataset is known to be lost.

   The listings record names, sizes and dates only — nothing in OneDrive has
   been opened — so any ``[O]``-only claim about file *contents* remains
   unverified.

The 0.1.x notebooks read from a ``data_root`` variable, so pointing them at a
different store is a one-cell edit. ``2023-02-02-2200 - F-Macro
Multi-Process.ipynb`` still carries a commented-out ``data_root`` line aimed at
``/home/bpaynter/Archive/lysis_data`` from an earlier such switch.

Data specification alignment
============================

Every Dataset logged here was surveyed for which degradation-output file it
carries, and the answer partitions them cleanly by date ``[D]``:

.. list-table::
   :header-rows: 1
   :widths: 26 16 20 38

   * - Datasets
     - Spec
     - Degradation file
     - Distinguishing layout
   * - ``2023-02-02-22xx``
     - **v1.85.0**
     - ``f_deg_time``
     - Flat — all Simulations combined in one top-level file
   * - ``2023-05-17`` … ``2023-12-10``
     - **v1.90.0**
     - ``f_deg_time``
     - One directory per Simulation (``00``–``09``)
   * - ``2024-02-02`` onward
     - **v1.95.0**
     - ``f_deg_list``
     - Per-Simulation dirs; ``params.json`` holds bare floats

**No v0.1.x Dataset contains ``f_deg_list``.** The changeover is a hard cutover
at ``2024-02-02-1400`` — the Tidy Fortran Data Dataset, whose entire purpose was
to make it (logged with v0.2.x). Confirmed across the canonical checkout, this worktree,
and four spot-checked archive tarballs ``[D]`` ``[A]``.

Note that the spec numbers run **opposite** to intuition when reading the
package source: ``v1.85.0`` is the oldest and ``v1.99.0`` the current Fortran
format, and each older spec is defined in ``dataspec.py`` as a modified copy of
the next *newer* one.

.. note::

   The ``v1.85.0`` test fixture (``tests/fixtures/fortran_v185_sample/``) is
   file-for-file a sample of ``2023-02-02-2200`` — it carries the same
   ``PLG2_tPA01_along_Q2`` file code ``[D]``. The tPA Diffusion baseline is
   therefore the Dataset that the modern v1.85.0 reader was written against.

Prior investigations
====================

These logs were written in one session, but they rest on archaeology done
across several earlier ones. What those established, so it is not re-derived:

**2026-08-04 — microscale/macroscale provenance.** Determined that the
``2023-02-01-22xx`` Datasets are the microscale inputs to the
``2023-02-02-22xx`` macroscale Runs, mapped three-to-twelve by Kd Scenario, and
that the micro log copied into each macroscale directory is a *verbose variant*
of the microscale original rather than the same file. Recorded in
:doc:`tpa_diffusion`.

**2026-08-04 — the v0.1.0 tag is cut from the wrong commit.** ``v0.1.0`` points
at ``368a6f2`` (2024-01-12), but the code used for the 2023 papers runs to
``63a01bb`` (2024-01-30). The two intervening commits are notebook-only — no
package code, Fortran, docs or Makefile changed. ``dev-0.1.1`` was rebased onto
``63a01bb`` to correct this for the re-release. Consequence for these logs: the
**tag** is not a reliable marker of the paper-era code; ``63a01bb`` is. See
:doc:`reexecution_environment`.

**2026-08-04 — paper figure lineage.** The B&W figures for the tPA Diffusion
paper come from ``63a01bb``; the earlier ``42263b5`` is a scratch intermediate
that renders a 2x2 into the top-left quadrant of a 4x4 grid and crops it out
with a hand-built ``Bbox``. Cell 31 of the 2200 notebook has been unchanged
since ``63a01bb``.

**2026-08-04 — the two checkouts render different figures.** Every one of the
59 PNGs common to ``lysis-v0.1.0/data/2023-02-02-220*`` and
``lysis/data/2023-02-02-220*`` differs, and not merely in embedded metadata: 13
differ in pixel *dimensions*, with the v0.1.0 renders generally far higher
resolution (e.g. ``microscale_measure_plots`` at 6190x3751 against 845x555).
Ten figures exist on only one side. 69 side-by-side diff composites were
rendered to ``lysis/data/2023-02-02-2200/diffs/`` (35 MB, gitignored).

   **Do not treat a figure as identified by filename alone.** The same name
   denotes different renders in the two trees.

**2026-08-05 — ``lenlysisvect`` has four conflicting definitions.** The
quantity is the position of the first ``6000`` in a sorted column of
``lysismat``, which reads as ``count + 1`` 1-indexed (MATLAB/Fortran) and
``count`` 0-indexed (Python). Verified empirically:
``lenlysisvect == argmax(lysismat, axis=0) + 1 == count_degraded + 1``, always.
The MATLAB generator also writes a magic ``999`` when a column contains no
``6000`` at all. The modern conversion path handles the index-base split
correctly; the shared all-degraded edge case became issue #122.

   Relevant to any Dataset here whose microscale summary tables are read
   directly — ``lenlysisvect`` counts are off by one from fibre counts if the
   index base is assumed rather than checked.

Conventions
===========

**Dataset code**
   A ``YYYY-MM-DD-HHMM`` timestamp identifying a group of Runs. The last two
   digits are the Run index within the group, so Dataset ``2023-02-02-2200``
   with ``group_code = "2023-02-02-22"`` contains Runs ``2023-02-02-2200``
   through ``2023-02-02-2214``. The group code and the first Run code are
   therefore easy to confuse — they are the same string.

**Category** (from the spreadsheet)
   ``Main Sequence``
      A step along the principal line of development. Supersedes its parent.
   ``Comparative``
      A branch run alongside its parent to isolate the effect of one change.
   ``HEAD``
      The furthest point reached on that line of work.

**Notebook pairs**
   Most Datasets have two notebooks: ``... - F-Macro Multi-Array-Run.ipynb``
   (dispatches the Simulations) and ``... - F-Macro Multi-Array-Process.ipynb``
   (reads the results and produces figures). Earlier ones use ``Multi-Run`` /
   ``Multi-Process`` without ``Array``. Notebooks under ``_Archive/`` are
   superseded but retained.
