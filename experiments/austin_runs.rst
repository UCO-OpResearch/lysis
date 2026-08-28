================================
austin-runs
================================

*How do the microscale outputs respond to order-of-magnitude changes in the tPA
and plasminogen binding kinetics — and can the archived answers be reproduced?*

Austin Segrest's microscale parameter sweep: **175 Runs of 50 000 Simulations
each** (8 750 000 Simulations), executed 2024-10 → 2026-03. In 2026 the archived
results were found not to reproduce, and the investigation that followed turned
this Experiment into a second, methodological one about the Fortran ``Lat``
allocation bug and the Intel Fortran compiler version. Both are logged here,
because the artefacts are inseparable — the corrected data set *is* an output of
the reproducibility investigation.

Companion document
==================

A much longer orientation document was written for this Experiment alongside
this log, aimed at students doing statistical analysis on the data:

   ``/shared/lysis-group/austin-runs.rst``

It sits **next to the data** and covers the full swept-parameter grid, the Run
naming quirks, per-folder inventories, worked ``lysis`` CLI commands and an
explicit exclusions section. **This log does not repeat it.** The division is:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Document
     - Covers
   * - ``/shared/lysis-group/austin-runs.rst``
     - What the data *is* and how to analyse it. Audience: students.
   * - This log
     - What was executed, when, with what code, and why each artefact
       supersedes the last. Audience: maintainers.

Note the filenames differ by one character and are **not** the same file:
``austin_runs.rst`` (this log, snake_case, matching the ``experiments/``
convention) and ``austin-runs.rst`` (the data-side document, hyphenated,
matching the Experiment name as it appears in folder and file codes).

Evidence tags
=============

Every claim below carries the evidence it rests on, using the vocabulary defined
in ``experiments/README.rst``. Six of those tags appear in this log; the
definitions are repeated here because that index file does not reach ``main``
until the v0.1.1 merge (see `Filing note`_).

``[S]``
   From a contemporaneous primary record — written at the time by the person who
   did the work. For this Experiment that means one of the two narrative
   ``.rst`` documents stored beside the data in ``/shared/lysis-group/`` (see
   `Sources`_), which the vocabulary counts as ``[S]`` in the same way as the
   spreadsheet index used by the v0.1.x logs.

``[H]``
   Read directly out of HDF5 attributes with ``h5py`` in mode ``"r"``. The
   strongest evidence available for this Experiment — the modern files carry
   full provenance stamps.

``[D]``
   Confirmed present on disk.

``[G]``
   From the ``lysis`` git history — a commit, branch or diff.

``[X]``
   Computed by cross-checking one artefact against another — comparing Run-code
   lists, seeds, hashes, or array contents across data sets. Reproducible:
   wherever this tag appears, the comparison that produced the claim is stated
   with it.

``[M]``
   **Maintainer testimony** — a recollection, not a contemporaneous record.
   Weaker than ``[S]``, and *not* to be conflated with ``[?]``: ``[M]`` says
   someone remembers this happening, ``[?]`` says this is inferred, neither
   remembered nor recorded. Several facts in this log rest on ``[M]`` alone and
   are marked as such; they are the answers the maintainer gave on 2026-08-28 to
   the questions this log originally left open.

``[?]``
   Inferred. Treat as a hypothesis, not a record.

``[N]``, ``[A]`` and ``[O]`` do not appear below: no notebook drove this
Experiment, and none of it is in the cold archive or on OneDrive.

Lineage
=======

Unlike :doc:`tpa_diffusion` and :doc:`internal_lysis`, the lineage here is not a
chain of Datasets superseding one another by *Scenario*. The Scenarios are fixed
from the start; what changes down the chain is the **format** and then the
**code and compiler** the same Scenarios were executed under.

.. code-block:: text

    old_data/                          original scattered tree — verified, then DELETED
        |  archive_script/archive.sh       flatten + de-duplicate, 199 units
        v
    old_data_compiled/                 [HEAD, raw]  199 run folders, 1.7 T
        |  verify_script/verify.sh         XXH128, 8783/8784 matched
        |  import to HDF5
        v
    austin-old-data-imported/          [Main Sequence]  145 .h5 — the diff baseline
        |
        |  re-execute 4 ways:
        |     {pre-lat, post-lat} x {ifort-2021.2.0, ifort-2021.9.0}
        v
    bpaynter/data/austin-runs__*/      [Comparative x4]  175 .h5 each
        |                                  |
        |                                  +--> lysis diff --> bpaynter/austin-runs__*.txt
        |                                                      bpaynter/data/austin-lat-bug-diffs/
        v
    wpumphrey/austin-runs-corrected/   [HEAD]  175 .h5 — authoritative
                                       == austin-runs__post-lat__ifort-2021.9.0  [X]

    (off to one side)
    wpumphrey/austin-runs-no-seedA/    input-only deck, no data — never executed [H]
    austin_segrest/austin-pre-lat-rerun/  input-only deck — DELETED 2026-08-21 [S]

The chain is reproduced from ``austin_segrest/README.rst`` ``[S]`` and confirmed
step by step against the artefacts ``[D]`` ``[H]`` ``[X]``.

Scenarios and Mechanisms
========================

**One Mechanism throughout.** All 175 Runs use the Fortran microscale model
``bin/micro_rates`` — ``micro_version = 'micro_rates'``,
``backend_type = 'fortran'``, constant across all 175 files ``[H]``. **No
Mechanism is varied anywhere in this Experiment.** That makes it the mirror
image of :doc:`tpa_diffusion`, which varies four Mechanisms over three
Scenarios.

**Microscale only.** The ``.h5`` files hold a ``micro_data`` Data Collection and
a ``log_files`` group and nothing else ``[H]``. The ``experiment.json`` files do
carry a 26-key ``macro_params`` block per Run, but those are unexecuted defaults
— there is no macroscale Dataset anywhere in the Experiment ``[H]``.

**Scenarios: a non-factorial sweep.** Seven kinetic parameters over a five-point
:math:`10^{-2}` … :math:`10^{2}` multiplier grid (suffixes ``_M0`` … ``_M4``,
so ``_M2`` is always the baseline), plus ``nodes_in_micro_row`` (5 or 7) and
``snap_proportion``. Derived by reading all 175 files and grouping Runs by which
parameters deviate from the modal value ``[H]`` ``[X]``:

.. list-table::
   :header-rows: 1
   :widths: 40 10 10 40

   * - Family
     - Runs
     - Nodes
     - Note
   * - ``<param>_M0`` … ``_M4``, 7 parameters
     - 45
     - 7
     - One-at-a-time. ``KdtPAplg_M*`` appears twice (``__v1``/``__v2``)
   * - ``ktPAon_M<i>_kplgon_M<j>``
     - 25
     - 7
     - Full 5 × 5
   * - ``ktPAon_M<i>_kplioff_M<j>``
     - 25
     - 7
     - Full 5 × 5
   * - ``ktPAon_M<i>_KdtPAyesplg_M<j>``
     - 25
     - 7
     - Full 5 × 5
   * - ``KdtPAyesplg_M<i>_KdtPAnoplg_M<i>``
     - 5
     - 7
     - Matched multipliers only — the diagonal
   * - Absolute-valued names
     - 29
     - 5
     - The same sweeps executed earlier at 5 nodes
   * - ``2131``–``2135`` / ``2141``–``2145``
     - 10
     - 5
     - 5-node ``ktPAon_M0`` × ``kplioff`` / × ``kplgon``
   * - ``sp_f12/13/14/18`` + ``sp_f23``
     - 5
     - 5
     - ``snap_proportion``; ``sp_f23`` is the 5-node baseline
   * - ``kdeg_kncat005/05/50/500``
     - 4
     - 5
     - ``deg_rate_fibrin`` and ``exposure_rate_binding_site`` together
   * - ``P0_M100``
     - 1
     - 7
     - Baseline Scenario; an unintentional replicant. **Part of the sweep** ``[M]``
   * - ``testing``
     - 1
     - 5
     - Baseline Scenario; a historic artefact, **not part of the sweep** ``[M]``

The full grid, with every value and unit, is in the companion document. Two
traps worth surfacing here because they bite maintainers, not just students:

.. warning::

   **``unbind_rate_*`` is derived, not swept.** It is the product of the
   corresponding dissociation constant and binding rate
   (``diss_const_PLG_intact`` 38 µM × ``bind_rate_PLG`` 0.1 µM⁻¹s⁻¹ =
   ``unbind_rate_PLG_intact`` 3.8 s⁻¹) ``[H]``. Likewise ``binding_sites``,
   ``fibrin_conc_per_fiber`` and ``protein_per_fiber`` each take exactly two
   values, splitting precisely on the 7-node / 5-node boundary ``[H]`` ``[X]``.
   None of these four are independent sweep axes.

   **14 Runs share one Scenario.** The ``_M2`` centre point of each family is a
   distinct Run at the identical baseline Scenario, differing only by
   ``micro_seed`` ``[H]`` ``[X]``. Pooled naively across families, the baseline
   is counted 14 times. ``P0_M100`` is one of the 14.

The two odd Run codes
---------------------

``P0_M100`` and ``testing`` both sit at the baseline Scenario and were the two
Runs this log could not originally account for. The maintainer has since
identified both ``[M]``.

``P0_M100`` — **an unintentional replicant. Keep it.**
   Its parameters are the plain 7-node baseline; it differs from the other
   thirteen baseline Runs only by ``micro_seed`` ``[H]``. The name implies a
   100× multiplier that no parameter carries. It was not meant to be a separate
   Scenario — it duplicates one — but it is a valid Run of 50 000 Simulations
   at a known Scenario with a known seed, and it **stays in the data set**.
   Treat it as one more baseline replicate, subject to the 14-Run
   double-counting caveat above.

``testing`` — **a historic artefact. Not part of the sweep Experiment.**
   Believed to be either the first Run executed on Buddy-2.5 or the last on
   Buddy-2.0 ``[M]``. It is a 5-node baseline Run, identical in parameters to
   ``sp_f23`` bar the seed ``[H]``.

   **This explains the one reproduction anomaly this log otherwise could not
   account for.** In the pre-``Lat`` / ifort-2021.2.0 configuration — the
   combination that reproduces the archive almost perfectly — exactly one of the
   145 Runs fails to reproduce, and that Run is ``testing`` ``[D]``. If it
   straddles a cluster rebuild, it is the one Run in the collection whose
   archived output was produced under a *different* platform state from every
   other, so a mismatch there is expected rather than anomalous. The 144/145
   result is therefore better read as **144 of 144** for the sweep proper, with
   ``testing`` correctly excluded.

   It is retained as a historic artefact. It should be **excluded from pooled
   analysis** of the sweep.

Datasets
========

old_data_compiled
-----------------

:Category: HEAD (raw) ``[S]``
:Derived from: — (origin) ``[S]``
:Executed: 2024-10-10 → 2026-03-31 (``.dat`` mtime span) ``[D]``
:Assembled: 2026-05-08 by ``archive_script/archive.sh``; ``archive.log`` records
   ``ok=199 skipped=0 failed=0`` ``[D]`` ``[S]``
:Extent: 199 run folders, 4 649 files, **1.7 T** ``[D]``
:Microscale: ``micro_rates.f90``, mixed code and compiler versions — see
   `The Lat bug and the compiler`_
:Macroscale: 24 of the 199 folders are ``macro_*`` units with ``00``–``10``
   sub-run folders ``[D]``. **Not** imported anywhere; see `Resolved questions`_
:Data specification: none — predates the HDF5 specification entirely. Flat
   ``.dat`` plus a ``micro_*.txt`` log per folder ``[D]``

The canonical, de-duplicated copy of the original raw Fortran output, and now
the **only** copy: the original ``old_data/`` tree was verified byte-for-byte by
``verify_script/verify.sh`` (XXH128; 8 783 of 8 784 files matched, the single
exception a documented expected drop) and then deleted ``[S]``.

Also holds ``experiment_results.csv`` (per-Run summary — seed, nodes, rate
constants, ``Ratio_LysComplete``, ``Median_Lysis``) and a copy of
``research.ipynb`` ``[D]``.

``__v1``/``__v2`` suffixes mark six basename collisions that were genuinely
divergent Runs; ``__v1`` is the older by mtime ``[S]``.

.. important::

   **30 of the 199 folders are log-only, and the reason is now known.**

   Thirty folders contain their ``micro_*.txt`` log and no ``.dat`` output at
   all ``[D]`` ``[X]``. They are the two-parameter cross sweeps
   ``ktPAon_M{0..4}_KdtPAyesplg_M{0..4}`` (25) and
   ``KdtPAyesplg_M{i}_KdtPAnoplg_M{i}`` (5) — exactly the set of Runs absent
   from ``austin-old-data-imported/`` ``[X]``.

   ``austin_segrest/README.rst`` records this as unresolved: *"Whether they were
   never written, were deleted before the archive was assembled, or were lost
   earlier is not recorded."* ``[S]``

   **The logs themselves settle it: the Simulations never executed.** Every one
   of the 30 logs terminates immediately after the line naming its first output
   file — ``data/<run_code>/lysis_<run_code>.dat`` — and **none** contains a
   single ``stats=`` progress line, where a complete Run's log carries ``stats=``
   at every thousandth Simulation up to ``stats=50000`` followed by a closing
   rate-constant dump ``[D]``. The parameter echo and the seed are present; the
   Simulation loop never produced a line.

   The mtimes confirm two aborted dispatch batches, each finishing in seconds
   ``[D]``:

   * the 25 ``ktPAon`` × ``KdtPAyesplg`` logs, **2025-09-12 23:14:54–57**
     (a 3-second window);
   * the 5 ``KdtPAyesplg`` × ``KdtPAnoplg`` logs, **2026-02-09 15:52:16–21**
     (a 5-second window).

   50 000 Simulations do not complete in three seconds. So the output was never
   written — not deleted later, and not lost by ``archive_script`` or
   ``verify_script``, both of which behaved correctly (they copied and verified
   faithfully; neither checks that a run folder is *complete*) ``[S]``.

   **They were recorded at the time as cancelled.**
   ``bpaynter/data/austin-lat-bug-diffs/Austin Diffs Report.xlsx`` indexes all
   199 units of ``old_data_compiled`` and marks exactly these 30
   ``Not included: Cancelled Run`` ``[S]``. That is a contemporaneous
   classification, and it agrees with what the logs show — a dispatch that
   stopped before doing any work.

   **What it does not say is why, or by whom, and no reason was recorded or has
   been determined since** ``[M]``. "Cancelled" is consistent with a scheduler
   kill, a manual cancellation, or a job that aborted and was written off;
   nothing distinguishes them. An earlier draft of this log proposed that the
   sweep drivers' ``mkdir -p data/$RUN_CODE`` step might have been skipped,
   making the output open fail ``[?]`` — that is a guess from where the logs
   stop, it is *less* well supported than the recorded "cancelled", and it
   should not be repeated as the explanation.

   The maintainer's recollection is that **the same Scenarios were executed to
   completion elsewhere** ``[M]`` — so these 30 aborted dispatches may be a
   duplicate attempt rather than the only one. Where that complete output went
   has not been established; it is not in ``old_data_compiled``.

   What survived here is the parameters and the seed, which is why these 30 Runs
   could later be recreated exactly rather than merely re-parameterised ``[S]``.

austin-old-data-imported
------------------------

:Category: Main Sequence ``[S]``
:Derived from: ``old_data_compiled`` ``[S]``
:Reason for change: import the raw Fortran output into HDF5 so it could be
   diffed against reruns ``[S]``
:Created: 2026-05-14, 11:33–11:43 — a single ten-minute batch ``[D]``
:Extent: **145** ``.h5``, 162 M ``[D]``. 145 of 175 Runs; the other 30 had no
   output to import ``[X]``
:Data specification: ``dataspec_version = 'v2.0.0'``,
   ``converted_from = 'v1.95.0'`` ``[H]``
:Also at: ``bpaynter/data/austin-old-data-imported`` — a **symlink** to this
   path, created 2026-05-25, not a second copy ``[D]``

The baseline every rerun was diffed against. The 145 seeds here match
``austin-runs-corrected`` exactly, Run for Run, with zero mismatches ``[X]`` —
that agreement is what makes the diffs meaningful.

.. note::

   **Two fingerprints of the legacy import path** ``[H]``, both worth knowing
   before writing code against this set:

   *Numeric parameters are floats.* ``micro_simulations`` is ``50000.0`` and
   ``nodes_in_micro_row`` is ``7.0``, against ``50000`` and ``7`` in
   ``austin-runs-corrected``.

   *The attribute set is thinner and irregular.* No ``fiber_radius``,
   ``binding_sites``, ``fibrin_conc_per_fiber``, ``protein_per_fiber``,
   ``fibrinogen_*``, ``protofibril_radius``, ``micro_version`` or
   ``micro_log_lvl``; ``snap_proportion`` appears on only 4 of 145 files. Three
   stray legacy attributes scraped out of the old logs — ``runcode``,
   ``percent_degrade`` and ``t`` — appear on a handful of files each. And there
   are **no provenance attributes** at all (``backend_*``, ``init_*``,
   ``pipeline_*``): these numbers were imported, not executed, so there was
   nothing to stamp.

bpaynter/data/austin-runs__{pre,post}-lat__ifort-2021.{2,9}.0
--------------------------------------------------------------

:Category: Comparative × 4 ``[S]``
:Derived from: ``austin-old-data-imported`` (parameters and seeds) ``[S]`` ``[X]``
:Reason for change: isolate the ``Lat`` bug from the compiler version ``[S]``
:Prepared: ``experiment.json`` ``created`` 2026-05-19T20:17 – 20:59 ``[D]``
:Extent: 175 ``.h5`` + ``experiment.json`` each, ≈192 M per folder ``[D]``
:Microscale: ``bin/micro_rates``. Post-``lat`` legs at ``3e34eb2c`` ``[H]``;
   pre-``lat`` legs from branch ``archive/hdf-pre-lat-fix`` ``[S]`` ``[G]``
:Data specification: ``v2.0.0``, ``converted_from = 'v1.99.0'`` ``[H]``

The four-way rerun matrix — {pre-``Lat`` fix, post-``Lat`` fix} ×
{ifort-2021.2.0, ifort-2021.9.0}. Every Run re-executed from the archived seed
and parameters, then diffed against the 145-Run baseline.

The exact ``uv run lysis run-micro`` command line for each of the four legs is
recorded in ``wpumphrey/austin-runs-corrected/austin-runs-corrected.rst``
``[S]`` — including the full ``--compiler`` module string needed to reach
ifort-2021.2.0. Do not reconstruct those by hand; copy them.

austin-runs-corrected
---------------------

:Category: **HEAD** — authoritative ``[S]``
:Derived from: ``austin-old-data-imported`` (seeds and parameters), regenerated
   under the post-``Lat`` codebase ``[S]`` ``[H]``
:Reason for change: the archived results were affected by the ``Lat`` bug in
   every 5-node configuration; this set replaces them ``[S]`` ``[X]``
:Executed / written: 2026-05-28 → 2026-06-02 ``[D]``
:Extent: **175** ``.h5``, 192 M, plus ``experiment.json``,
   ``austin-runs-corrected.rst``, ``Austin Diff Report - bpaynter.xlsx`` and
   ``.slurm/`` (350 job logs) ``[D]``
:Microscale: ``bin/micro_rates`` at ``backend_commit = 3e34eb2c``,
   ``backend_compiler`` Intel Fortran **2021.9.0** Build 20230302,
   ``backend_dirty = 'clean'``; pipeline at ``1dccf4c6`` ``[H]``
:Macroscale: none executed ``[H]``
:Data specification: ``dataspec_version = 'v2.0.0'``,
   ``converted_from = 'v1.99.0'`` ``[H]``

**Use this one.** Both ``austin_segrest/README.rst`` and
``austin-runs-corrected.rst`` name it the authoritative replacement ``[S]``.

Three things distinguish it:

* **It is the only complete set.** 175 Runs, against 145 in the baseline. For
  the 30 cross-sweep Runs whose Simulations never executed in 2025/2026, this
  is not a *corrected* replacement — it is the **only copy of those results that
  has ever existed** ``[S]`` ``[X]``.
* **It is internally consistent.** One code version and one compiler across all
  175 Runs ``[H]``, which ``old_data_compiled`` is not: that archive accumulated
  across two compiler versions and both sides of the ``Lat`` fix.
* **Clean modern attributes** — integers where integers belong, the full
  geometry set, full provenance stamps ``[H]``.

.. important::

   **It is the same data as** ``bpaynter/data/austin-runs__post-lat__ifort-2021.9.0``.

   Neither existing document says so, and the two folders look like separate
   Runs. Established three ways ``[X]``:

   #. Both ``experiment.json`` files carry
      ``name = "austin-runs__post-lat__ifort-2021.9.0"`` and the same
      ``created = "2026-05-19T20:37:49"``.
   #. Every provenance attribute matches — ``backend_commit``,
      ``backend_compiler``, ``init_timestamp``, ``init_hostname``,
      ``pipeline_timestamp``, ``pipeline_hostname`` ``[H]``.
   #. A 20-Run random sample compared array-by-array and attribute-by-attribute
      with ``h5py``: **zero** differences. ``uv run lysis diff`` reports all 8
      tables ``OK``.

   The ``.h5`` files are **not** byte-identical — different md5 — but that is
   HDF5 container layout after a copy, not content. ``austin-runs-corrected/``
   is the published copy of one leg of the four-way matrix.

austin-runs-no-seedA
--------------------

:Category: input-only deck — **holds no data** ``[H]``
:Derived from: the same source CSV family as the corrected set ``[S]``
:Reason: test whether results depend on the archived seeds ``[S]``
:Created: 2026-05-27 23:38 — all 176 files within one second ``[D]``
:Extent: 175 ``.h5`` at exactly 18 368 bytes each + ``experiment.json``,
   3.7 M ``[D]``
:Status: **never executed** — confirmed by the maintainer ``[M]``, and visible
   in the files: every Dataset is present but empty, shape ``(0,)`` ``[H]``
   ``[S]``

**An abandoned idea, recorded here because that is what an experiment log is
for.** The plan was to re-execute the sweep without the archived seeds, to see
whether the results depended on them. The decks were built and never run. The
maintainer deliberately omits this artefact from the Experiment's folder listing
in the companion document ``[M]``: it holds no data, so it has nothing to offer
a student, and its name invites the mistake of counting it as a fourth data set.

It is already documented in detail — including the full parameter comparison
against ``austin-runs-corrected`` and the verification that all 175 files are
empty — in:

   ``/shared/lysis-group/austin_segrest/README.rst``

That account is not repeated here. In brief, its parameters match
``austin-runs-corrected`` with two exceptions ``[S]``:

* ``micro_seed`` is ``0`` in all 175 files, against 175 distinct archived seeds.
* ``snap_proportion`` for ``sp_f13`` is truncated to ``0.333333333``, inherited
  from ``wpumphrey/data/austin-runs-no-seed.csv``. The archived Fortran log
  records the full value, so **this deck is the deviant copy**.

It exists at two paths — ``wpumphrey/austin-runs-no-seedA`` and
``wpumphrey/data/austin-runs-no-seedA`` — which are two real directories, not
symlinks, holding the same 175 files ``[D]`` ``[X]``.

A sibling deck, ``austin_segrest/austin-pre-lat-rerun/`` (the *seeded* variant),
was **deleted 2026-08-21** as redundant after a node-by-node comparison found
nothing in it that was not also in ``austin-runs-corrected`` ``[S]``. Its
``experiment.json`` metadata was preserved in ``austin_segrest/README.rst``
rather than in a folder of empty files. ``wpumphrey/data/austin-pre-lat-rerun``
is now a **dangling symlink** to the deleted target, left in place because
``wpumphrey/data/`` is not group-writable ``[S]``.

The Lat bug and the compiler
============================

**The bug.** Commit ``b44355a`` (2026-04-23, *"Fixing fortran microscale 'Lat'
allocation bug"*) adds a **single line**, ``Lat = 0``, to
``src/fortran/micro_rates.f90``, immediately before the loop that writes 1s
along the diagonal of ``Lat`` ``[G]``. Before that line existed the off-diagonal
entries held whatever was in the freshly allocated memory — which is why the bug
and the compiler version are entangled, and why the effect is not deterministic
across builds.

Neither ``austin-runs-corrected.rst`` nor ``austin_segrest/README.rst`` states
the mechanism; both name the bug without explaining it.

**Compilers.** Cluster nicknames map to compiler versions, and both appear in
filenames ``[S]``:

* **Sooner** → ``ifort-2021.2.0``
* **Buddy** → ``ifort-2021.9.0``

Buddy was upgraded around July 2025; ``ifort-2021.2.0`` is no longer available
there and must be reached through Sooner ``[S]``.

**Results.** Counted directly from the four Markdown tables at
``bpaynter/austin-runs__*.txt`` (145 Runs each; columns
``Run | Result | Max % Diff | Worst Table``) ``[D]``:

.. list-table::
   :header-rows: 1
   :widths: 46 18 18 18

   * - Configuration
     - Exact Match
     - Differ
     - Note
   * - pre-``Lat`` × ifort-2021.2.0
     - **144**
     - 1
     - the one is the ``testing`` Run
   * - pre-``Lat`` × ifort-2021.9.0
     - 97
     - 48
     -
   * - post-``Lat`` × ifort-2021.9.0
     - 96
     - 49
     -
   * - post-``Lat`` × ifort-2021.2.0
     - 96
     - 49
     - report file **byte-identical** to the row above

Two readings follow. Pre-``Lat`` code built with ifort-2021.2.0 reproduces the
archive almost completely — so that is the combination the bulk of the archive
was originally produced with. And once the bug is fixed, the compiler version
stops mattering: the two post-``Lat`` reports are the same file ``[D]``.

.. important::

   **The node-count finding, confirmed independently.**
   ``austin-runs-corrected.rst`` attributes the non-reproducing Runs to 5-node
   configurations, and explains the apparent "September 2025 transition" as the
   adoption of 7-node Runs rather than any software change ``[S]``.

   That checks out exactly. The 49 Runs marked "differ" in the post-``Lat``
   reports are, as a set, **precisely** the 49 Runs whose
   ``nodes_in_micro_row`` attribute is ``5``. The two sorted lists are
   identical, not merely similar ``[H]`` ``[X]``.

   Consequence for anyone reading ``austin-old-data-imported``: its 49 five-node
   Runs carry the bug and are not comparable with its 96 seven-node Runs.

Related artefacts
=================

Scripts and reports supporting the Experiment, all under ``/shared/lysis-group/``:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Path
     - What it is
   * - ``austin_segrest/README.rst``
     - The folder map and full provenance chain. Read first ``[S]``
   * - ``austin_segrest/archive_script/``
     - ``archive.sh``, ``dest_map.tsv`` (199 rows), ``archive.log``,
       ``README.txt`` — how ``old_data/`` was flattened
   * - ``austin_segrest/verify_script/``
     - ``verify.sh`` (three-phase resumable XXH128) + hash TSVs — the proof the
       archive was faithful before ``old_data/`` was deleted
   * - ``austin_segrest/bash_scripts/``
     - 7 SLURM sweep drivers, one per parameter. **Load
       ``intel-compilers/2021.4.0``** — a *third* compiler, distinct from either
       leg of the investigation ``[D]``. Post-hoc reconstruction ``[S]``
   * - ``austin_segrest/experiment_csv/``
     - The same 7 sweeps as ``lysis`` Experiment CSVs. Also written after the
       fact; the folder README warns they may be unreliable ``[S]``
   * - ``austin_segrest/code/research.ipynb``
     - Parses the legacy ``.txt``/``.dat`` output into dataframes — the tool for
       reading ``old_data_compiled`` directly
   * - ``wpumphrey/data/austin-runs.csv``
     - **The authoritative sweep definition**: 16 parameter rows × 175 Run
       columns; the file the Run decks were built from ``[D]``
   * - ``bpaynter/austin-runs__*.txt`` (4)
     - Per-Run diff results as Markdown tables ``[D]``
   * - ``bpaynter/Austin Diff Report - bpaynter.xlsx``
     - Spreadsheet roll-up of those four; copied into
       ``austin-runs-corrected/``
   * - ``bpaynter/data/austin-lat-bug-diffs/``
     - ``generate_file_list.sh`` (builds the 145-entry baseline list) and
       ``split_matching_diffs.sh`` (executes ``uv run lysis diff`` per file).
       Duplicated at ``wpumphrey/data/austin-lat-bug-diffs/`` — see below

The two copies of ``austin-lat-bug-diffs``
------------------------------------------

The diff bookkeeping exists at two paths. **They are two independent copies, and
their contents are byte-identical.** Established by hashing both trees ``[X]``:

* Neither is a symlink; both are real directories with different inodes, so they
  are not hardlinks either.
* Both hold the same 10 files. Every one matches on ``sha256sum``, and
  ``diff -rq`` across the trees reports no differences.
* ``wpumphrey``'s copy is the **earlier** (files 2026-05-24 19:41);
  ``bpaynter``'s is the later (2026-05-25 10:40), so the copy went
  wpumphrey → bpaynter. This matches the maintainer's account that Tyler began
  the diff analysis and ``bpaynter`` completed it ``[M]``.

The ``matched_`` / ``remaining_`` lists are a **cascade**, not four independent
runs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the thing to understand before reading those files, and it is easy to
get wrong. ``split_matching_diffs.sh`` takes a *file list* as input, and the
workflow feeds each configuration the **previous configuration's**
``remaining_`` output — the plan set out in
`issue #27 <https://github.com/UCO-OpResearch/lysis/issues/27>`_ ("Record the
matches and remove from list for the next run"). Verified by set comparison
``[X]``:

==================================  =======  =======  =========  ==========
Step / configuration                Input    Matched  Remaining  Output kept
==================================  =======  =======  =========  ==========
1. post-``Lat`` × ifort-2021.9.0    145          96         49   yes
2. post-``Lat`` × ifort-2021.2.0     49           0         49   yes (empty)
3. pre-``Lat`` × ifort-2021.9.0      49           1         48   yes
4. pre-``Lat`` × ifort-2021.2.0      48          48          0   **no files**
==================================  =======  =======  =========  ==========

Each step's input set is exactly the previous step's ``remaining_`` set ``[X]``.

So the empty ``matched_..._post-lat__ifort-2021.2.0.txt`` does **not** mean that
configuration reproduced nothing — its standalone report shows 96 exact matches.
It means none of the 49 runs *still outstanding at that point* matched. Read the
two side by side without knowing this and they look contradictory.

What these folders hold that the reports do not
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The per-Run *verdicts* are not unique: every run in every ``matched_`` /
``remaining_`` file appears in the corresponding ``bpaynter/austin-runs__*.txt``
report with a consistent verdict, and the reports are strictly richer — they
carry ``Max % Diff`` and ``Worst Table`` columns and cover all 145 runs for all
four configurations ``[X]``. Three things, however, exist **only** here:

**1.** ``Austin Diffs Report.xlsx`` **— the exclusion register.** This is the
most valuable file in either folder. It is indexed on **199** run codes — every
unit in ``old_data_compiled``, not the 145 the reports cover — and records, for
each of the 54 units the reports omit, *why* it was omitted ``[S]``:

* **24 × "Not included: Macroscale Sim"** — the ``macro_*`` units.
* **30 × "Not included: Cancelled Run"** — the log-only cross sweeps.

The remaining 145 rows each carry a single checkmark, in the column for the
first configuration in which that run reproduced: 96 / 0 / 1 / 48 across the
four steps, totalling 145 ``[X]``. That is the cascade above, summarised.

.. note::

   This is a **different file** from ``bpaynter/Austin Diff Report - bpaynter.xlsx``
   (note the singular "Diff" and the suffix). The top-level one is later
   (13:48 vs 10:40) and much richer: four sheets, adding Kolmogorov–Smirnov
   statistics and p-values on the lysis-time and tPA-leaving-time distributions,
   plus percentage differences on five summary measures, for the
   Pre→Post (Sooner), Pre→Post (Buddy) and Sooner→Buddy (Pre-Lat) comparisons —
   49 runs each. It also adds ``nodes_in_micro_row`` and a modified-date column
   to its Sheet1, which is what makes the node-count conclusion visible at a
   glance. Neither spreadsheet contains the other: only the diffs-folder copy
   carries the 54 exclusion reasons ``[D]``.

**2. The method.** ``generate_file_list.sh`` (builds the 145-name baseline list)
and ``split_matching_diffs.sh``. The reports record results, not how they were
produced. In particular the match criterion is defined here and nowhere else: a
run counts as reproduced when ``uv run lysis diff --markdown`` returns ``OK`` for
**all 8** tables ``[D]``.

**3. The cascade order itself**, which is not recoverable from the reports —
they are standalone per-configuration tables over all 145 runs.

Step 4 is the one genuine gap: ``pre-lat__ifort-2021.2.0`` has no ``matched_`` or
``remaining_`` file at all. But the spreadsheet records its 48 matches, so the
*analysis* was completed even though that step's script output was not kept
``[X]``.

Not part of this Experiment
===========================

The companion document carries the full exclusions list for students — the
separate lysis-front and binding-site Experiments, ``polymerization/``,
``papers/``, and the 2023 Python-vs-Fortran macroscale tarball. Two exclusions
need recording here as well, because both are *adjacent* to this Experiment and
a maintainer could reasonably mistake either for part of it.

``austin_segrest/lysis/data/`` (144 G) — Austin's incomplete work
-----------------------------------------------------------------

**Not part of this Experiment** ``[M]``. Austin's own clone of the repository
holds 45 folders of raw Fortran output named ``micro_<param>_scale_<scale>`` and
``macro_<param>_scale_<scale>`` ``[D]``. The parameters are the same seven, and
the scales the same :math:`10^{-2}` … :math:`10^{2}`, which is exactly why it
invites confusion.

The distinction is one of completeness. **Austin compiled the data he considered
complete, and** ``old_data_compiled`` **is what was derived from that** ``[M]``.
What remains in ``lysis/data/`` is the residue of the work he did *not* consider
finished. It was never imported to HDF5, never diffed against the archive, and
was built under ``intel-compilers/2021.4.0`` ``[D]`` — a third compiler, distinct
from both legs of the ``Lat`` investigation, so it is not even comparable with
the rerun matrix.

.. warning::

   **Do not build analysis on** ``austin_segrest/lysis/data/``. It is expected
   to be archived or erased ``[M]``. Anything derived from it will not be
   reproducible once it goes, and it has no provenance stamps to reconstruct it
   from.

The 24 ``macro_*`` units in ``old_data_compiled`` — pending, not discarded
--------------------------------------------------------------------------

Out of scope **for now**, but explicitly *not* abandoned ``[M]``. Twenty-four of
the 199 folders in ``old_data_compiled`` are named ``macro_*`` and hold
macroscale output in ``00``–``10`` sub-run folders ``[D]``. They correspond to
macroscale Runs derived from some of the austin-runs microscale Scenarios.

They were never imported into any of the three HDF5 data sets — all three are
microscale-only ``[H]`` — and the diff baseline list has 145 entries, every one
microscale ``[D]``. That is why this log treats the Experiment as
microscale-only.

They were, however, **consciously excluded rather than overlooked**.
``austin-lat-bug-diffs/Austin Diffs Report.xlsx`` indexes all 199 units of
``old_data_compiled`` and marks each of these 24 ``Not included: Macroscale
Sim`` ``[S]``. Someone looked at them, classified them, and set them aside — a
better starting point for the eventual import than a gap in the record would
have been.

**They still need to be imported and analysed.** The work is captured against:

* `issue #25 — Import Austin's Data <https://github.com/UCO-OpResearch/lysis/issues/25>`_
* `issue #27 — Recreate Austin's Data Exactly <https://github.com/UCO-OpResearch/lysis/issues/27>`_

Both are open against the **v1.0.0** milestone, and #27 is a sub-issue of #25.
A note recording the 24 folders and their status was added to #25 on 2026-08-28.

.. note::

   Issue #27 independently corroborates a finding in this log. Its stated task
   is to *"create an 'austin-runs-corrected' folder with
   'post_lat__ifort_2021.9.0' contents, spreadsheet, and rst"* — which is
   precisely what ``austin-runs-corrected/`` turned out to be when compared
   against ``bpaynter/data/austin-runs__post-lat__ifort-2021.9.0`` ``[X]``. The
   equivalence was designed in, not accidental.

Data specification alignment
============================

This Experiment sits **entirely past** the v0.1.x span logged in
:doc:`tpa_diffusion` and :doc:`internal_lysis`, which end at v1.95.0. Both HDF5
sets here are ``v2.0.0`` ``[H]``:

.. list-table::
   :header-rows: 1
   :widths: 34 18 22 26

   * - Data set
     - ``dataspec_version``
     - ``converted_from``
     - Path into the format
   * - ``old_data_compiled``
     - *(none)*
     - —
     - Raw Fortran ``.dat``; pre-specification
   * - ``austin-old-data-imported``
     - ``v2.0.0``
     - ``v1.95.0``
     - Legacy import, then converted forward
   * - ``austin-runs__*`` (4)
     - ``v2.0.0``
     - ``v1.99.0``
     - Current execution pipeline
   * - ``austin-runs-corrected``
     - ``v2.0.0``
     - ``v1.99.0``
     - Current execution pipeline

The two ``converted_from`` values are the cleanest single discriminator between
the imported baseline and everything re-executed, and they are readable without
opening a Dataset — ``v1.95.0`` means *imported from the archive*, ``v1.99.0``
means *executed by the modern pipeline* ``[H]``.

Each Run's ``micro_data`` group holds the eight standard microscale Datasets,
each of length 50 000 ``[H]``: ``fiber_degraded``, ``sim_final_time``,
``pli_first_time``, ``pli_generated_num``, ``tpa_final_num``,
``tpa_leaving_time``, ``tpa_unbound_by_pli``, ``tpa_unbound_kinetic``. Note the
HDF5 group is named ``micro_data`` while ``lysis diff`` and the data
specification call the same Data Collection ``microscale_out``.

Filing note
===========

This log lives in ``experiments/`` on ``main``. It is currently the only file in
that directory: the index ``README.rst``, and the sibling logs
:doc:`tpa_diffusion`, :doc:`internal_lysis` and :doc:`reexecution_environment`
that this file cross-references, all arrive when the **v0.1.1** line is merged.
Until then those four cross-references do not resolve. Nothing builds this
directory — ``experiments/`` is not wired into any Sphinx toctree — so nothing
is broken in the meantime.

**The tag vocabulary this log needs is already in place.** ``[H]``, ``[G]``,
``[X]`` and ``[M]`` were added to ``experiments/README.rst`` on the
``dev-0.1.1`` line on 2026-08-28, and ``[S]`` was generalised there from "the
spreadsheet index" to "a contemporaneous primary record" — which absorbs the
narrative ``.rst`` documents this log cites, so no separate tag for them is
needed. The definitions in `Evidence tags`_ above track that file; once it
reaches ``main`` they can be reduced to a pointer. **Do not modify**
``experiments/README.rst`` — it is maintained on ``dev-0.1.1``.

One mechanical edit falls due at the merge, and is deliberately *not* applied on
``dev-0.1.1`` because this log lands on ``main``. Add this row to the log
list-table, after the ``reexecution_environment`` row:

.. code-block:: rst

   * - :doc:`austin_runs`
     - 2024-10 → 2026-06
     - Microscale rate-constant sweep; the ``Lat`` bug and compiler reproducibility

Resolved questions
==================

This log originally left seven questions open. The maintainer answered all of
them on 2026-08-28 ``[M]``; two were also settled from the artefacts ``[X]``.
Recorded here so the reasoning is not rediscovered.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Question
     - Answer
   * - What is ``P0_M100``?
     - An unintentional replicant of the baseline Scenario. **Kept** in the data
       set ``[M]``. See `The two odd Run codes`_.
   * - Is ``testing`` a throwaway?
     - A historic artefact — probably the first Run on Buddy-2.5 or the last on
       Buddy-2.0. Retained, but **not part of the sweep Experiment**, and it
       explains the single pre-``Lat``/2021.2.0 reproduction anomaly ``[M]``.
   * - Is ``austin_segrest/lysis/data/`` (144 G) part of this Experiment?
     - **No.** It is the remnant of Austin's *incomplete* work ``[M]``. See
       `Not part of this Experiment`_.
   * - Were the no-seed Runs ever executed?
     - **No** ``[M]``, matching the files themselves ``[H]``. See
       `austin-runs-no-seedA`_.
   * - Are the 24 ``macro_*`` units in scope?
     - Out of scope **for now**, but pending import — not discarded. See
       `Not part of this Experiment`_.
   * - What caused the 30 aborted dispatches?
     - Recorded at the time as ``Not included: Cancelled Run`` ``[S]``, but
       **no reason was given or determined** ``[M]``. That the Simulations never
       executed is established ``[D]``; *why* is not. See the warning under
       `old_data_compiled`_.
   * - Are the two ``austin-lat-bug-diffs`` folders one artefact or two?
     - **Two independent copies, byte-identical**, wpumphrey → bpaynter ``[X]``.
       Their ``matched_``/``remaining_`` lists are a cascade, and their
       spreadsheet carries the 54 exclusion reasons found nowhere else. See
       `The two copies of austin-lat-bug-diffs`_.

Still open
==========

* **Where did the complete counterparts of the 30 aborted Runs go?** The
  maintainer recalls the same Scenarios being executed to completion somewhere
  ``[M]``, but that output is not in ``old_data_compiled`` and has not been
  located. If it is found, it would give a genuine archival baseline for the 30
  cross-sweep Runs, which currently exist only in ``austin-runs-corrected``.
* **When will the 24 ``macro_*`` units be imported?** Tracked against
  `issue #25 <https://github.com/UCO-OpResearch/lysis/issues/25>`_.

Sources
=======

Everything above is either quoted from a document that already existed, read out
of the data, or computed by cross-checking artefacts. In descending order of
authority:

``[S]`` ``/shared/lysis-group/austin_segrest/README.rst``
   Folder map, provenance chain, the 30 log-only Runs, the deleted decks.

``[S]`` ``/shared/lysis-group/wpumphrey/austin-runs-corrected/austin-runs-corrected.rst``
   The investigation narrative, the node-count conclusion, and the four
   ``run-micro`` command lines.

``[H]`` The HDF5 files themselves
   Opened with ``h5py`` in mode ``"r"``. All 175 files in
   ``austin-runs-corrected/`` and all 145 in ``austin-old-data-imported/`` were
   scanned attribute by attribute.

``[D]`` ``bpaynter/austin-runs__*.txt``, ``old_data_compiled/*/micro_*.txt``, mtimes, ``experiment.json``
   The diff tables, the aborted Run logs, and file metadata.

``[D]`` ``austin-runs-analysis-inventory.rst``
   The companion inventory of every analysis artefact — what each contains,
   what is unique to it, which are byte-identical duplicates, and the plan for
   compiling one canonical set. The duplication map and keep/drop decisions
   summarised in this log are derived there.

``[S]`` ``austin-lat-bug-diffs/Austin Diffs Report.xlsx``
   The exclusion register: 199 units, with the 24 ``Macroscale Sim`` and 30
   ``Cancelled Run`` classifications, and the cascade of first-match
   configurations.

``[G]`` The ``lysis`` git history
   Commit ``b44355a`` (the ``Lat`` fix), ``3e34eb2c`` and ``1dccf4c6`` (stamped
   in the data), branch ``archive/hdf-pre-lat-fix``.

``[X]`` Cross-checks performed 2026-08-28
   Comparisons of Run-code lists, seeds and array contents across the data sets.
   Each is stated inline with the claim it supports, so any of them can be
   re-executed from this file alone.

.. note::

   All of the above was gathered read-only. **Nothing under
   ``/shared/lysis-group/`` was created, modified, moved or deleted** in the
   course of writing this log — the data sets it describes are unchanged by it.
