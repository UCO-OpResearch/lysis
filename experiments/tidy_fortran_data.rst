================================
Tidy Fortran Data
================================

*Can the macroscale model record fibre degradation as an event log instead of a
dense snapshot array, without changing the science?*

A single-Dataset Experiment. It re-executes the ten :doc:`internal_lysis`
Scenarios **unchanged** — same parameters, same seeds — and varies only the
output format. The Runs are therefore a controlled comparison of file formats,
not of model behaviour.

.. important::

   **This Dataset sits on the far side of the v0.1.x boundary.** It is logged
   here because its data lineage is continuous with :doc:`internal_lysis`, but
   the code and notebooks that produced it are **v0.2.x**. See
   `Version seam`_ below before relying on this log's placement.

Lineage
=======

.. code-block:: text

    2023-12-10-1900        (Internal Lysis HEAD of that branch)
        |
    2024-02-02-1400  [Comparative]   f_deg_time array -> f_deg_list event log

Dataset
=======

2024-02-02-1400
---------------

:Category: Comparative ``[S]``
:Derived from: ``2023-12-10-1900`` ``[S]``
:Reason for change: **"No change to input data, new output format"** ``[S]``
:Notebooks: ``2024-02-02-1400 - F-Macro Multi-Array-{Run,Process}.ipynb`` —
   **not in this worktree**; live in the ``lysis-v0.2.0`` worktree and archived
   in ``main`` under ``notebooks/_Archive/`` ``[N]``
:Microscale: ``micro_rates.f90`` ``[S]``
:Macroscale: ``macro_diffuse_into_and_along__internal.f90`` ``[S]`` ``[N]``
:Data: present in the ``main`` checkout, 10 Run directories ``[D]``; **absent
   from this worktree and from the cold archive** ``[A]``

**Inputs are byte-identical to the parent.** The ``scenarios`` array, the
``mechanisms`` array and all ten seeds match ``2023-12-10-1900`` exactly,
including the seed-sequence entropy comment
(``3881821051554698152964433817123076384``) ``[N]``. This is what makes the
Dataset a usable format control.

**The change, concretely** ``[D]``. Per-Run output is identical in every file
except one:

.. list-table::
   :header-rows: 1
   :widths: 26 18 18 38

   * - File
     - Dataset
     - Size
     - Form
   * - ``f_deg_time_*.dat``
     - ``2023-12-10-1900``
     - 19.4 MB
     - Binary snapshot array — one degradation time per fibre per save
   * - ``f_deg_list_*.dat``
     - ``2024-02-02-1400``
     - 1.5 MB
     - CSV event log — ``time, fibre_index, degrade_time`` per event

A **12x reduction**, which is evidently the point: commit ``4996830`` of the
same day adds a "Data Size" notebook ``[?]``.

.. note::

   **Run ``02`` (TN-D_684) is commented out** of the ``runs`` array, so only
   nine Runs are processed ``[N]`` — but all ten data directories
   ``1400``–``1409`` exist ``[D]``. The data for ``1402`` was produced and then
   excluded from processing; no reason is recorded ``[?]``.

Version seam
============

The placement of this log is the one open question about it. The evidence:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Criterion
     - Points to
   * - Data lineage
     - **v0.1.x** — parent ``2023-12-10-1900`` is an Internal Lysis Dataset
   * - Producing commits
     - **v0.2.x** — ``4b7e129`` "Made f_deg_list changes to Fortran" and
       ``f676e01`` "Completed runs for Fortran data redesign" (both
       2024-02-02) are **not** ancestors of ``63a01bb``, the v0.1.x tip, but
       **are** ancestors of ``v0.2.0``
   * - Notebooks
     - **v0.2.x** — live in the ``lysis-v0.2.0`` worktree; absent from this one
   * - Branch
     - **v0.2.x** — developed on ``tidy-fortran-data``, merged into
       ``lysis-front-study`` on 2024-02-05 (``892d714``)

Three of four criteria put this Dataset in v0.2.x, alongside the Lysis Front
Experiment's ``2024-01-26-1000`` — whose notebooks are likewise live only in the
v0.2.0 worktree. If the v0.2.x logs are organised by *code line*, this file
should move there and leave a cross-reference behind. If they are organised by
*data lineage*, it belongs here.

.. note::

   The v0.1.x tip is ``63a01bb`` (2024-01-30); this Dataset was produced three
   days later. The ``v0.1.0`` **tag** itself is older still — ``368a6f2``
   (2024-01-12) — and is known to have been cut from the wrong commit.

Afterlife
=========

``f_deg_list`` **is the format that stuck.** It is still the current Fortran
representation: the base data specification ``v1.99.0`` defines
``f_deg_list{file_code}_{sim:02}.dat`` as a three-field text event log
(``Simulation Time Elapsed``, ``Grid Location Index``,
``Fiber New Degrade Time``) ``[D]``.

The older specs are built by deriving *backwards* from ``v1.99.0``, and
``_create_v1_90()`` is the one that removes it:

.. code-block:: text

    v1.99.0  f_deg_list   <- base (current Fortran spec)
    v1.95.0  f_deg_list   <- copy of v1.99.0, pre-Pint parameters
    v1.90.0  f_deg_time   <- copy of v1.95.0, swaps the event log for the array
    v1.85.0  f_deg_time   <- copy of v1.90.0, combined-simulation layout

So this Experiment marks the **permanent** transition from ``f_deg_time`` to
``f_deg_list``, and the spec boundary between ``v1.90.0`` and ``v1.95.0`` falls
exactly here.

.. note::

   An earlier draft of this log recorded the opposite — that ``f_deg_list`` had
   been superseded by a return to ``f_deg_time`` in v1.90.0. That was a
   misreading of the version numbers: v1.90.0 is **older** than v1.95.0, not
   newer. The modern dataset ``2026-02-28-1907`` does contain ``f_deg_time``,
   but it was produced in 2026 specifically to exercise the legacy v1.90.0
   reader, and is the source of ``tests/fixtures/fortran_v190_sample/`` ``[D]``.

Data specification
==================

``2024-02-02-1400`` aligns with **v1.95.0** ``[D]``:

* ``f_deg_list`` rather than ``f_deg_time`` — rules out v1.90.0 and v1.85.0.
* One directory per Simulation (``00``–``09``) — rules out v1.85.0's
  combined-simulation layout.
* ``params.json`` holds **bare floats** (``"binding_rate": 0.1``) and
  ``"micro_params": null``, so microscale parameters must be parsed out of
  ``micro_PLG2_tPA01.txt`` — this is precisely what ``_create_v1_95()``
  describes, and it rules out v1.99.0, where parameters are Pint Quantities.

This makes ``2024-02-02-1400`` the **earliest v1.95.0 Dataset** in the
collection, and its parent ``2023-12-10-1900`` the latest v1.90.0 one.
