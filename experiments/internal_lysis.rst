================================
Internal Lysis
================================

*What happens when lysis begins from tPA molecules that start* **inside** *the
clot, rather than arriving at its face?*

Published as the 2023 "Internal Lysis" paper; a Zenodo dataset was generated
for it (commit ``c2993ae``). This Experiment has by far the longest lineage in
the project — nine Datasets over eight months. Six of the nine have **no data in
any live checkout**, but all but two of those were recovered from the cold
archive at ``~/Archive/lysis_data/``.

Lineage
=======

.. code-block:: text

    2023-02-02-2200                  (origin — see tpa_diffusion)
        |
    2023-04-13-2100  [Comparative]   remaining-bind-time unbinding + internal start
        |
    2023-04-15-1400  [Main]          first Internal Lysis parameters
        |
    2023-04-17-1400  [Main]          'degrade' -> 't_degrade' tracking
        |
    2023-04-20-1400  [Main]          fixed the degradation-tracking bug
        |
    2023-04-23-2000  [Main]          reordered loops, separated Simulations
        |
        +---------------------+---------------------+
        |                     |                     |
    2023-05-23-1000       2023-05-23-1100       2023-05-23-2000
      [Main]                [Comparative]         [Comparative]
      height + empty rows   empty rows only       height only
        |
    2023-06-09-1000  [HEAD]          pore size & transit time data
        |
    2023-12-04-1000  [Comparative]   molecule locations at each bind/unbind
        |
    2023-12-10-1900  [Comparative]   binding locations on thick fibers
        |
        `-- 2024-02-02-1400           ("Tidy Fortran Data" — separate Experiment)

.. note::

   The three-way split at ``2023-05-23`` is a **factorial isolation**, and the
   only one in the project ``[?]``. ``1000`` applies both changes (adjusted clot
   height *and* eliminated empty rows), while ``1100`` and ``2000`` each apply
   one. All three derive from ``2023-04-23-2000``, not from each other ``[S]``.
   The combined Dataset is the one that carries the Main Sequence forward.

   Only the combined arm's data has been located. Both single-factor arms are
   missing, which is the single biggest gap in this Experiment's record.

Scenario naming
===============

From ``2023-04-15-1400`` onward the Scenarios use codes of the form
``T{N|K}-{L|D}_{molecules}``, e.g. ``TK-L_307``. Decoding them ``[?]``, from the
parameter values in the notebooks ``[N]``:

.. list-table::
   :header-rows: 1
   :widths: 14 20 30 36

   * - Token
     - Reading
     - Parameter
     - Values
   * - ``TN``
     - **T**\ hi\ **n** fibres
     - ``fiber_diameter``
     - 72.7 nm
   * - ``TK``
     - **T**\ hic\ **k** fibres
     - ``fiber_diameter``
     - 145.4 nm (exactly 2x thin)
   * - ``L``
     - **L**\ oose
     - ``pore_size``
     - 1.0135 micron
   * - ``D``
     - **D**\ ense
     - ``pore_size``
     - 0.22 micron
   * - ``_NNNN``
     - molecule count
     - ``total_molecules``
     - 307 / 684 / 3042 / 9350

The thin/thick and loose/dense readings are inferences from the numbers; the
letters are never expanded in any notebook. Grid dimensions track the pore size
(loose thin = 93x93, dense thin = 342x116) and ``forced_unbind`` tracks the
fibre thickness (0.0852 thin, 0.0729129 thick) ``[N]``.

Ten Scenarios were used from ``2023-04-23-2000`` onward, and the same ten,
with the same seeds, recur unchanged through ``2023-12-10-1900`` ``[N]``. The
Mechanism is constant across the whole Experiment:
``Into and along - Internal`` / ``macro_diffuse_into_and_along__internal``
``[N]`` ``[S]``.

Datasets
========

2023-04-13-2100
---------------

:Category: Comparative ``[S]``
:Derived from: ``2023-02-02-2200`` ``[S]``
:Reason for change: **"Changed macro-unbinding wait time to remaining bind time, added internal molecule start"** ``[S]``
:Notebooks: ``_Archive/2023-04-13-2100 - F-Macro Multi-{Run,Process}.ipynb`` ``[N]``
:Macroscale: recorded as ``macro_diffuse_into_and_along__XXXX.f90`` — arm not identified ``[S]``
:Data: archived — ``02 03 07 08 12 13`` ``[A]``; absent from all checkouts ``[D]``

The branch point. Still on the tPA Diffusion Scenario set (Physiological /
10x smaller / 10x bigger Kd) with Runs at ``02``, ``07`` and ``12``, under a new
Mechanism ``Internal - Remaining Unbind Time`` ``[N]``. Three further Run entries
are present but commented out ``[N]``.

2023-04-15-1400
---------------

:Category: Main Sequence ``[S]``
:Derived from: ``2023-04-13-2100`` ``[S]``
:Reason for change: **"Initial runs with Internal Lysis parameters"** ``[S]``
:Notebooks: ``_Archive/2023-04-15-1400 - F-Macro Multi-{Run,Process}.ipynb`` ``[N]``
:Data: archived — all ten, ``00``–``09`` ``[A]``

First appearance of the ``T*-*`` Scenarios — 5 Runs: ``TN-L_9350``, ``TN-L_307``,
``TK-L_3042``, ``TK-L_9350``, ``TK-L_307`` (loose only) ``[N]``.

2023-04-17-1400
---------------

:Category: Main Sequence ``[S]``
:Derived from: ``2023-04-15-1400`` ``[S]``
:Reason for change: **"Changed tracking of fiber degradation from 'degrade' to 't_degrade'"** ``[S]``
:Notebooks: ``_Archive/2023-04-17-1400 - F-Macro Multi-{Run,Process}.ipynb`` ``[N]``
:Data: archived — ``02 03 04 08 09`` ``[A]``

4 Runs, all **dense** Scenarios: ``TN-D_684``, ``TN-D_307``, ``TK-D_307``,
``TK-D_9350`` ``[N]``. Complements its parent's loose-only set rather than
repeating it.

2023-04-20-1400
---------------

:Category: Main Sequence ``[S]``
:Derived from: ``2023-04-17-1400`` ``[S]``
:Reason for change: **"Fixed bug with fiber degradation tracking"** ``[S]``
:Notebooks: ``_Archive/2023-04-20-1400 - F-Macro Multi-{Run,Process}.ipynb`` ``[N]``
:Data: archived — ``03`` and ``05`` ``[A]``

A single Run — ``TK-L_3042`` at index ``05`` ``[N]``. Consistent with a
spot-check of the fix rather than a full re-execution ``[?]``.

2023-04-23-2000
---------------

:Category: Main Sequence ``[S]``
:Derived from: ``2023-04-20-1400`` ``[S]``
:Reason for change: **"Reordered loops, separated simulations"** ``[S]``
:Notebooks: ``_Archive/2023-04-23-2000 - F-Macro Multi-Array-{Run,Process}.ipynb`` ``[N]``
:Data: archived — ``00``–``09`` plus ``12`` and ``14`` ``[A]``

First **Array** notebook, and the first full 10-Scenario sweep (``00``–``09``)
``[N]``. The switch from ``Multi-Process`` to ``Multi-Array-Process`` marks the
move to SLURM array jobs ``[?]``. The ten seeds fixed here persist unchanged
through the rest of the Experiment ``[N]``.

2023-05-23-1000
---------------

:Category: Main Sequence ``[S]``
:Derived from: ``2023-04-23-2000`` ``[S]``
:Reason for change: **"Adjusted clot height, eliminated empty rows, selectable save interval"** ``[S]``
:Notebooks: ``_Archive/2023-05-23-1000 - F-Macro Multi-Array-{Run,Process}.ipynb`` ``[N]``
:Data: archived — all ten, ``00``–``09`` ``[A]``

Both structural changes together, plus a new per-Run save interval — the Run
tuples gain a sixth field (10 or 100) absent from ``2023-04-23-2000`` ``[N]``.

2023-05-23-1100
---------------

:Category: Comparative ``[S]``
:Derived from: ``2023-04-23-2000`` ``[S]``
:Reason for change: **"Eliminated empty rows"** ``[S]``
:Notebooks: ``_Archive/2023-05-23-1100 - F-Macro Multi-Array-{Run,Process}.ipynb`` ``[N]``
:Data: **not located** — absent from all checkouts and from the archive; OneDrive unsearched ``[D]`` ``[A]``

2023-05-23-2000
---------------

:Category: Comparative ``[S]``
:Derived from: ``2023-04-23-2000`` ``[S]``
:Reason for change: **"Adjusted clot height"** ``[S]``
:Notebooks: ``_Archive/2023-05-23-2000 - F-Macro Multi-Array-{Run,Process}.ipynb`` ``[N]``
:Data: **not located** — absent from all checkouts and from the archive; OneDrive unsearched ``[D]`` ``[A]``

.. note::

   ``_Archive/2023-05-23-2200 - Compare.ipynb`` has no Dataset of its own and no
   entry in the spreadsheet ``[N]`` ``[S]``. It is presumably the notebook that
   compared the three ``2023-05-23`` arms ``[?]``. It is also the only notebook
   in the project that actually *uses* ``networkx``, for a connected-components
   analysis of degraded fibres.

2023-06-09-1000
---------------

:Category: HEAD ``[S]``
:Derived from: ``2023-05-23-1000`` ``[S]``
:Reason for change: **"Pore size & Transit Time data"** ``[S]``
:Notebooks: ``2023-06-09-1000 - F-Macro Multi-Array-{Run,Process}.ipynb`` (live) ``[N]``
:Data: present in the ``main`` checkout, 10 Run directories; **absent from this worktree** ``[D]``

The tip of the Main Sequence and the reference Dataset for the paper. Adds pore
size and transit time collection. Full 10-Scenario sweep ``[N]``. Related commit:
``a0d7dbb`` "Added notebook code to dump pore size and transit time data".

2023-12-04-1000
---------------

:Category: Comparative ``[S]``
:Derived from: ``2023-06-09-1000`` ``[S]``
:Reason for change: **"Output molecule locations at each binding/unbinding event"** ``[S]``
:Notebooks: ``2023-12-04-1000 - F-Macro Multi-Array-{Run,Process}.ipynb`` (live) ``[N]``
:Data: archived — all ten, ``00``–``09`` ``[A]``; absent from all checkouts ``[D]``

Same ten Scenarios and seeds as its parent; the change is in what the macroscale
model *records*, not what it simulates ``[N]``. Related commit: ``8795732``
"Added molecule tracking".

.. note::

   The notebook is live and un-archived — it looks current — but its data is in
   no checkout, only in the cold archive. Re-executing this notebook therefore
   requires unpacking ``2023-12-04-10{00..09}.tar.gz`` first and repointing
   ``data_root``. Nothing in the notebook says so.

2023-12-10-1900
---------------

:Category: Comparative ``[S]``
:Derived from: ``2023-12-04-1000`` ``[S]``
:Reason for change: **"Changed number of binding locations for thick fibers"** ``[S]``
:Notebooks: ``2023-12-10-1900 - F-Macro Multi-Array-{Run,Process}.ipynb`` (live) ``[N]``
:Data: present, 10 Run directories ``[D]``

The last Dataset in the Experiment. The change is directly visible in the
notebook: the ``scenario_type`` dtype gains a field between ``fiber_diameter``
and ``cols``, set to **427** for thin-fibre Scenarios and **213** for thick
``[N]`` — the binding-site count per fibre.

Downstream
==========

``2024-02-02-1400`` derives from ``2023-12-10-1900`` but belongs to a different
Experiment Sequence, **"Tidy Fortran Data"** — *"No change to input data, new
output format"* ``[S]``. It is a format migration, not an Internal Lysis result,
and its data is present in the ``main`` checkout ``[D]``. Log it separately if
that line of work is ever written up.

Open questions
==============

* **The archive holds Runs the notebooks never processed.** Four Datasets have
  more archived Runs than their Process notebook lists: ``2023-04-15-1400``
  (10 archived vs 5 processed), ``2023-04-17-1400`` (5 vs 4),
  ``2023-04-20-1400`` (2 vs 1) and ``2023-04-23-2000`` (12 archived, including
  indices ``12`` and ``14``, vs 10 processed) ``[A]`` ``[N]``. Whether these are
  abandoned Runs, failed Runs, or results that were simply never written up is
  unknown ``[?]``.

* The macroscale arm for ``2023-04-13-2100`` is recorded only as
  ``macro_diffuse_into_and_along__XXXX.f90``. Given the Mechanism descriptor
  (``Internal - Remaining Unbind Time``) it is presumably the ``__internal``
  arm ``[?]``, but the spreadsheet's placeholder suggests uncertainty at the
  time of writing.
* ``2023-05-23-1100`` and ``2023-05-23-2000`` are the only Datasets in this
  Experiment not located anywhere. They are the two single-factor arms of the
  factorial split, so losing both leaves the combined ``1000`` Dataset without
  its controls. **OneDrive has not been searched** and is the obvious next
  place to look; the Zenodo deposit (commit ``c2993ae``) is a second lead.
* ``2023-05-17-1400`` ("Model Schematic") sits between ``2023-05-23`` Datasets
  in time and shares the Internal Lysis machinery, but the spreadsheet assigns
  it its own Sequence with no parent ``[S]``. Its data survives in ``main``
  ``[D]``. Treated here as out of scope.

Processing status (August 2026)
===============================

Only ``2023-12-10-1900`` has been re-executed under the ``dev-0.1.1`` uv
environment — ``In[1]``–``In[35]``, with 4 failures:

* ``In[22]``, ``In[33]`` — ``ImportError: Missing optional dependency
  'tabulate'``. **Fixed** by commit ``ba57528``; not yet re-executed.
* ``In[31]`` — ``NameError: name 'macro_unbound_times' is not defined``.
* ``In[32]`` — ``TypeError: find_degradation_marker_frames() missing 1 required
  positional argument: 'degrade_percent_markers'``.

The last two are stale calls against a changed helper API and remain open.
``2023-06-09-1000`` has not been re-executed here (its data lives in the ``main``
checkout); ``2023-12-04-1000`` cannot be, as its data is lost.
