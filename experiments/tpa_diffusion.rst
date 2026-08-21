================================
tPA Diffusion
================================

*How far does tPA travel into a clot before it binds, and how does that depend
on the rule governing where a molecule goes after it unbinds?*

Published as the 2023 "tPA diffusion" paper. Figures for that paper were
produced by the ``2023-02-02-2200`` Process notebook; the B&W-compatible
revisions are commits ``42263b5`` and ``63a01bb`` (2024-01-30).

Lineage
=======

.. code-block:: text

    2023-02-02-2200   [HEAD]        4 Mechanisms x 3 Kd Scenarios
        |
        +-- 2023-11-24-1000  [Comparative]   slower microscale diffusion

Two Datasets only. Unlike :doc:`internal_lysis`, this Experiment did not go
through a long Main Sequence — the second Dataset is a Comparative branch that
re-uses the baseline and adds slowed-diffusion variants.

Datasets
========

2023-02-02-2200
---------------

:Category: HEAD ``[S]``
:Derived from: — (origin of both this Experiment and :doc:`internal_lysis`) ``[S]``
:Notebooks: ``2023-02-02-2200 - F-Macro Multi-Process.ipynb`` ``[N]``
:Microscale: ``micro_rates.f90`` ``[S]``
:Macroscale: recorded only as ``macro_XXXXX.f90`` — **not identified** ``[S]`` ``[?]``
:Data: present in both checkouts, 13 Run directories ``[D]``; **not** in the
   cold archive ``[A]``
:Data specification: **v1.85.0** — ``f_deg_time``, and the combined-simulation
   flat layout (one ``f_deg_time_PLG2_tPA01_along_Q2.dat`` per Run directory,
   no ``00``–``09`` subdirectories) ``[D]``. This Dataset is the source of
   ``tests/fixtures/fortran_v185_sample/``.

The baseline. A 3 x 4 grid of Scenarios and Mechanisms.

**Scenarios** (dissociation constant of tPA) ``[N]``

.. list-table::
   :header-rows: 1
   :widths: 22 18 18 20

   * - Descriptor
     - File code
     - ``diss_const``
     - ``k_off``
   * - 10x smaller Kd
     - ``_Kd00020036``
     - 0.5143
     - 277.8
   * - Physiological Kd
     - *(none)*
     - 8.52e-2
     - 27.8
   * - 10x bigger Kd
     - ``_Kd0236``
     - 5.4e-3
     - 2.78

**Mechanisms** (where a molecule may go after unbinding) ``[N]``

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Descriptor
     - File code
     - Executable
   * - Always bind
     - ``_always``
     - ``macro_Q2_always_rebind``
   * - Diffuse along clot
     - ``_along``
     - ``macro_Q2_diffuse_along``
   * - Diffuse into clot
     - ``_into``
     - ``macro_Q2_diffuse_into``
   * - Diffuse into and along clot
     - ``_into_and_along_fixed``
     - ``macro_Q2_diffuse_into_and_along_fixed``

.. important::

   **Runs 02, 07 and 12 are excluded.** They used an earlier
   ``macro_Q2_diffuse_into_and_along`` executable, labelled in the notebook as
   ``"Diffuse into and along clot - BUGGED"``, and were superseded by the
   ``_fixed`` variant at Runs 03, 08 and 13. Both the Mechanism entry and the
   three Run entries are commented out rather than deleted ``[N]``.

   Only ``2023-02-02-2202`` still exists on disk; ``2207`` and ``2212`` are
   gone ``[D]``. Do not process ``2202`` — it is bugged output.

The 12 valid Runs: ``00`` ``01`` ``03`` ``04`` (Physiological), ``05`` ``06``
``08`` ``09`` (10x smaller), ``10`` ``11`` ``13`` ``14`` (10x bigger) ``[N]``.
Each carries an explicit seed and a wall-clock budget of 900–1800 s ``[N]``.

2023-11-24-1000
---------------

:Category: Comparative ``[S]``
:Derived from: ``2023-02-02-2200`` ``[S]``
:Reason for change: **"Experiments with slower micro diffusion"** ``[S]``
:Notebooks: ``2023-11-24-1000 - F-Macro Multi-Array-{Run,Process}.ipynb`` ``[N]``
:Microscale: ``micro_rates.f90`` ``[S]``
:Macroscale: three executables, one per Mechanism ``[N]`` —
   ``macro_diffuse_into_and_along__external`` (baseline) and
   ``macro_diffuse_into_and_along_slow_micro__external`` (2x and 4x).
   The spreadsheet records only the ``slow_micro`` one ``[S]``.
:Data: present in both checkouts, 9 Run directories ``[D]``; also archived,
   ``03 08 13 15 16 17 18 19 20`` ``[A]``
:Data specification: **v1.90.0** — ``f_deg_time``, one directory per
   Simulation ``[D]``

Revisits the paper's conclusions by slowing the diffusion of fibrin degradation
products at the microscale. Keeps the three Kd Scenarios and replaces the four
unbinding Mechanisms with three diffusion speeds:

.. list-table::
   :header-rows: 1
   :widths: 30 20 20 20

   * - Mechanism
     - Physiological
     - 10x smaller
     - 10x bigger
   * - Baseline diffusion
     - ``1003``
     - ``1008``
     - ``1013``
   * - 2x slower FDP diffusion
     - ``1015``
     - ``1016``
     - ``1017``
   * - 4x slower FDP diffusion
     - ``1018``
     - ``1019``
     - ``1020``

All nine exist on disk in both checkouts, and all nine are additionally held
as tarballs in ``~/Archive/lysis_data/`` ``[D]`` ``[A]``.

.. note::

   **"Baseline diffusion" is a rename, not a new Mechanism** ``[?]``. Runs
   ``03``/``08``/``13`` occupy the same indices as the ``_into_and_along_fixed``
   Runs of ``2023-02-02-2200`` and the notebook carries the same commented-out
   ladder of superseded Mechanisms at ``00``–``02``, ``04``–``07``, ``09``–``12``,
   ``14`` ``[N]``. Confirmed by the ``mechanisms`` array: the entry reads
   ``("_into_and_along", "Baseline diffusion",
   "macro_diffuse_into_and_along__external")`` — same ``_into_and_along`` file
   code as the parent's fixed Mechanism, new descriptor ``[N]``. The Dataset
   was laid out on the parent's index scheme and only the into-and-along arm
   was re-executed.

Open questions
==============

* The macroscale executable for ``2023-02-02-2200`` is recorded only as the
  placeholder ``macro_XXXXX.f90``. The Mechanism table gives per-Run executable
  names (``macro_Q2_*``), so the source file may be recoverable from
  ``archive/fortran/`` in the modern repository.
* Whether the ``2023-11-24-1000`` "Baseline diffusion" Runs were re-executed or
  copied from the parent Dataset is unresolved. Comparing seeds would settle it:
  the parent's Run ``03`` seed is ``-2137354075`` while ``1003``'s is
  ``1896691809``, which suggests **re-execution with fresh seeds** ``[N]`` ``[?]``.

Processing status (August 2026)
===============================

Both Process notebooks were re-executed under the ``dev-0.1.1`` uv environment.

* ``2023-02-02-2200`` — ``In[1]``–``In[52]``, **no errors**. Animations
  re-enabled; 12 ``combined_animation_*.mp4`` rendered.
* ``2023-11-24-1000`` — ``In[16]``–``In[28]``, **2 failures**: ``In[23]`` and
  ``In[26]`` both raise ``IndexError: index 0 is out of bounds for axis 0 with
  size 0``, from a lookup of the Mechanism descriptor
  ``"Diffuse into and along clot"``, which does not exist in this Dataset's
  ``mechanisms`` array (it was renamed to ``"Baseline diffusion"``). These cells
  appear to have been copied from the parent notebook without updating the name.
