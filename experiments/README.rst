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
   * - :doc:`tidy_fortran_data`
     - 2024-02
     - Degradation output as an event log instead of a snapshot array

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
   * - OneDrive
     - offsite
     - Further archives. **Not surveyed** — contents unknown to these logs.

.. note::

   Because OneDrive has not been searched, a Dataset marked *not located* below
   means exactly that: not found in the three places that *were* searched. It is
   **not** a claim that the data is gone. Only two Datasets across both
   Experiments fall into this category.

The 0.1.x notebooks read from a ``data_root`` variable, so pointing them at a
different store is a one-cell edit. ``2023-02-02-2200 - F-Macro
Multi-Process.ipynb`` still carries a commented-out ``data_root`` line aimed at
``/home/bpaynter/Archive/lysis_data`` from an earlier such switch.

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
