=============================================
Statistics Commands: A Reference
=============================================

.. note::

   Anything written in the form ``<___>`` — for example ``<_data_path_>`` — is a
   placeholder.  Replace the whole token, angle brackets and underscores
   included, with your own value.  Do not leave the brackets in the command.

Overview
========

This page documents the seven ``lysis`` commands that report numbers rather
than produce data:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Command
     - What it reports
   * - ``lysis micro-stats``
     - Fiber-degradation counts and lysis times, plus tPA leaving times, from
       the microscale output of one or more Runs.
   * - ``lysis macro-stats``
     - Six macroscale summary metrics (degradation rate, lag time, time to full
       degradation, tPA back-row percentage, first-passage time, lysis-front
       velocity) for one or more Runs.
   * - ``lysis deg-rate``
     - Macroscale degradation rates over **configurable** degradation
       intervals (default 20 %→80 %, 20 %→50 %, 50 %→80 %).
   * - ``lysis deg-time``
     - Macroscale times to reach **configurable** degradation milestones
       (default 5 %, 20 %, 50 %, 80 %, 100 %).
   * - ``lysis parameters``
     - The Scenario parameters stored in each Run's file, with units, and a
       re-feedable CSV export.
   * - ``lysis compare micro-stats``
     - Two-sample Kolmogorov–Smirnov tests plus percent differences between two
       sets of microscale Runs.
   * - ``lysis compare macro-stats``
     - The same, for macroscale Runs.

All seven are **read-only**.  They open each ``.h5`` file in HDF5 read mode,
compute in memory, and print.  They never write to the data file.  The only
files they can create are the ones you explicitly name with ``--markdown`` or
``--csv``.

Who should use them: anyone doing statistical analysis on existing simulation
output.  They are the fastest way to get a table of numbers out of a directory
of Runs and into R, Python, or a spreadsheet, without writing any HDF5 code.

.. note::

   ``lysis diff`` also compares two Runs, but element by element over the raw
   HDF5 tables.  It is a diagnostic tool for checking whether two Runs produced
   identical data, not an analysis tool, and it is not documented here — see
   ``lysis diff --help``.  For *statistical* comparison of two Runs, use
   :ref:`statistics_cli_compare`.

.. _statistics_cli_common:

Common concepts
===============

Everything in this section applies to all five commands.

The ``PATH`` argument
---------------------

Every command takes one or more paths.  A path may be either:

**A single** ``.h5`` **file**
    The Run code is the filename with the ``.h5`` extension removed.  A file
    named ``TF-x__9_951.h5`` is the Run whose run code is ``TF-x__9_951``.
    Output is a single-Run report.

**A directory**
    Every entry in that directory whose name ends in ``.h5``
    (case-insensitively) is treated as one Run.  The scan is **not**
    recursive — only the directory's own entries are considered, and
    subdirectories are ignored.  Output is a table with one row per Run.

The path must already exist; Click rejects a missing path before the command
runs (exit code 2).

For ``lysis compare``, both paths must be of the *same* kind: two directories
or two ``.h5`` files.  Mixing a directory with a file is rejected.

How a Run is identified
-----------------------

A **Run** (see :doc:`ontology`) is identified by its **run code**, which is
always the HDF5 file's basename without the extension.  There is no separate
identifier stored inside the file that the CLI consults.  This has one
important consequence for ``lysis compare`` in directory mode: the two
directories are matched by *run code*, so ``dirA/2131.h5`` is compared against
``dirB/2131.h5``.  Run codes that appear in only one of the two directories are
silently skipped.

.. _statistics_cli_sort:

``--sort``
----------

``--sort`` chooses the row order in directory mode.  Default: ``smart``.

``alpha``
    Plain lexicographic sort of the run codes (Python's ``sorted``).

``smart``
    Splits each run code on ``__`` (double underscore) or ``-`` (hyphen), then
    sorts position by position.  A position where *every* run code carries a
    valid Roman numeral is sorted by Roman value (so ``v`` < ``vii`` < ``ix`` <
    ``xi`` < ``xiii``); a position where every code carries a Python-style
    integer, optionally with underscore thousands separators, is sorted
    numerically (so ``9_951`` < ``21_105`` < ``1_582_867``).  Positions that
    are not uniformly numeric are compared as text.  If the codes cannot be
    tokenised into a consistent number of positions, ``smart`` falls back to a
    plain lexicographic sort.

    Source of truth: ``src/lysis/tools/runcode_sort.py``.

``--markdown``
--------------

``--markdown <_output_file_>`` emits the result as a GitHub-flavoured Markdown
table instead of a Rich terminal table.  Use ``-`` in place of the filename to
print to standard output.  Passing ``--markdown`` implies ``--no-progress``.

Markdown output is what you want when you are going to read the table into
another tool.  A Markdown pipe table is trivially parsed in R with
``readr::read_delim(..., delim = "|")`` or in Python with
``pandas.read_csv(..., sep="|")`` after stripping the separator row.

.. warning::

   Error and warning messages are printed to **standard output**, the same
   stream ``--markdown -`` writes the table to.  If a single Run in a
   directory fails to load, its error message lands in the middle of your
   piped table.  When scripting, prefer ``--markdown <_output_file_>`` over
   ``--markdown -``: the table goes to the file and only the messages remain on
   the terminal.

``--no-progress``
-----------------

``--no-progress`` suppresses the Rich progress bar and status spinner.  Set it
when executing under Slurm or redirecting output to a log file.  It is implied
by ``--markdown``.

Global options
--------------

These belong to the ``lysis`` group itself and must appear *before* the
subcommand name, e.g. ``lysis -v micro-stats <_data_path_>``.

.. program:: lysis

.. option:: -v, --verbose

   Increase verbosity.  Repeatable (``-v``, ``-vv``).

.. option:: --version

   Print the ``lysis`` version and exit.

Exit codes and failure modes
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Code
     - Meaning
   * - ``0``
     - At least one Run was reported successfully.
   * - ``1``
     - Nothing could be reported: the single named file failed to load or
       compute; or the directory contained no ``.h5`` files; or every Run in
       the directory failed; or (``compare`` only) the two paths were of
       different kinds, were not both ``.h5`` files, or shared no run codes.
   * - ``2``
     - Click usage error: the path does not exist, an option value is invalid,
       or the measure-set name is not one of ``micro-stats`` / ``macro-stats``.

In **directory mode a single failing Run is not fatal**.  Its error is printed
and the Run is dropped from the table; the command still exits ``0`` as long as
one other Run succeeded.  Check the row count, not just the exit code.

Two failures you are likely to meet with the shared data:

``Error computing statistics for <_run_code_>: 'NoneType' object has no attribute 'macro_simulations'``
    The file has no macroscale Data Collection.  ``lysis macro-stats`` and
    ``lysis compare macro-stats`` need one.  Microscale-only files —
    which includes every file under
    ``/shared/lysis-group/austin_segrest/austin-old-data-imported`` and
    ``/shared/lysis-group/wpumphrey/austin-runs-corrected`` — support
    ``micro-stats``, ``parameters``, and ``compare micro-stats`` only.

``Error loading <_run_code_>: seed must be an int or str, got float64``
    The file stores its RNG seed as a float.  Some older converted files do.
    The Run cannot be opened by the current reader at all, so every command
    fails on it.

.. _statistics_cli_notation:

Notation
========

Every formula below uses the following symbols.  All symbols are defined here
and reused unchanged.

Microscale (from the ``microscale_out`` Data Collection, see
:doc:`data_specification`)
--------------------------------------------------------------------------

:math:`N_\mu`
    Number of microscale Simulations in the Run
    (``micro_params.micro_simulations``); the length of every microscale
    dataset.

:math:`d_k \in \{0, 1\}`
    ``fiber_degraded[k]`` — whether microscale Simulation :math:`k` ended with
    the fiber fully degraded.  Stored as a boolean.

:math:`T_k`
    ``sim_final_time[k]`` — elapsed time of microscale Simulation :math:`k`, in
    **seconds**.

:math:`L_k`
    ``tpa_leaving_time[k]`` — the simulation time at which the tPA molecule
    left the system in Simulation :math:`k`, in **seconds**.

:math:`\mathcal{D} = \{\, k : d_k = 1 \,\}`
    The set of Simulations whose fiber degraded.  :math:`N_D = |\mathcal{D}|`.

Macroscale (from the ``macroscale_out`` Data Collection)
--------------------------------------------------------

:math:`S`
    Number of macroscale Simulations (``macro_params.macro_simulations``).
    Simulations are indexed :math:`s = 0, \dots, S-1` and stored as
    ``macro_data/sim_00``, ``macro_data/sim_01``, ...

:math:`K_s` and :math:`t^{(s)}_0 < \dots < t^{(s)}_{K_s - 1}`
    The number of save points in Simulation :math:`s`, and the save-point times
    ``snapshot_time``, in **seconds**.

:math:`M`
    Number of tPA molecules (``macro_params.total_molecules``); the length of
    ``tpa_transit_time`` in every Simulation.

:math:`\tau^{(s)}_m`
    ``tpa_transit_time[m]`` in Simulation :math:`s` — the time at which
    molecule :math:`m` first reached the back row of the fiber grid, in
    **seconds**.  A molecule that never reached the back row is stored as
    exactly :math:`0`.

:math:`R,\ C,\ E`
    ``rows``, ``cols`` and ``empty_rows`` of the edge grid.

:math:`W = 3C - 1`, :math:`Z = 2C - 1`
    ``full_row`` (edges in a normal row) and ``xz_row`` (edges in the final
    row, which has no y-edges).

:math:`F = W\,(R - E - 1) + Z`
    ``total_fibers`` — the number of fibrin fibers in the grid.

:math:`WE`
    ``empty_edges`` — the number of edges in the fibrin-free region.

:math:`\Sigma^{(s)}(t) \in \mathbb{R}^{R \times W}`
    The **replayed degrade-time state** of Simulation :math:`s` at time
    :math:`t`.  Entry :math:`\Sigma^{(s)}(t)_{ij}` is the currently scheduled
    degradation time (seconds) of the fiber at grid row :math:`i`, rank
    :math:`j`.  It is built by
    ``FiberReplayCursor``, which starts from

    .. math::

       \Sigma^{(s)}(0)_{ij} =
       \begin{cases}
         0        & i < E                        \quad \text{(fibrin-free rows)}\\
         \text{NaN} & i = R-1 \text{ and } j \ge Z \quad \text{(padding)}\\
         \infty   & \text{otherwise}             \quad \text{(not yet scheduled)}
       \end{cases}

    and then applies, in time order, every record of the ``fiber_degrade_time``
    event log whose ``Simulation Time Elapsed`` is :math:`\le t`, writing that
    record's ``Fiber New Degrade Time`` into the named ``(row, rank)`` cell.

    Source of truth: ``src/lysis/analysis/fiber_replay.py``.

:math:`\Delta`
    ``grid_node_distance`` in **microns**, computed as
    :math:`\texttt{pore\_size} + 2 \cdot \texttt{fiber\_radius}`.

:math:`y_i = i\,\Delta`
    The y-distance (microns) of grid row :math:`i` from the front of the clot.

Degradation milestones
----------------------

Two commands on this page — :ref:`statistics_cli_deg_rate` and
:ref:`statistics_cli_deg_time` — locate the moment a Simulation crosses a given
degradation threshold, so the crossing is given a symbol here.

:math:`k_s(\theta)`
    The **milestone frame**: the index of the first save point of Simulation
    :math:`s` at which the degraded fraction reaches the threshold
    :math:`\theta \in [0, 1]`.

    .. math::

       k_s(\theta) = \min\bigl\{\, k : \phi^{(s)}_k \ge \theta \,\bigr\}

    Implemented as ``numpy.argmax`` over the boolean array
    :math:`\phi^{(s)} \ge \theta`, which returns ``0`` when the threshold is
    never reached — see
    `issue #126 <https://github.com/UCO-OpResearch/lysis/issues/126>`_ and the
    warnings that cite it.  Source of truth:
    ``find_degradation_marker_frames`` in
    ``src/lysis/analysis/degradation.py``.

:math:`T_s(\theta)`
    The **milestone time**, in **minutes**:

    .. math::

       T_s(\theta) = \frac{t^{(s)}_{k_s(\theta)}}{60}

    Source of truth: ``find_degradation_marker_times`` in the same module.

Each :math:`k_s(\theta)` depends only on its own threshold, so asking for more
milestones or more intervals never changes the value of the ones you already
had.

Statistical conventions
-----------------------

Mean and standard deviation
    :math:`\operatorname{mean}(x) = \frac{1}{n}\sum_{i=1}^{n} x_i`.

    .. math::

       \operatorname{sd}(x) = \sqrt{\frac{1}{n}\sum_{i=1}^{n}(x_i - \operatorname{mean}(x))^2}

    Every standard deviation reported by these commands is the **population**
    form, :math:`\text{ddof} = 0` — the default of ``numpy.std``.  None of the
    code passes ``ddof=1``.  If you need the sample standard deviation, rescale
    by :math:`\sqrt{n / (n-1)}`, or recompute from the raw datasets.

Median
    ``numpy.median``, i.e. the 50th percentile with **linear interpolation**:
    for even :math:`n` the median is the arithmetic mean of the two middle
    order statistics.  No other quantile is reported by these commands.

Line fits
    Every "rate" or "velocity" below is the slope of a degree-1 fit by
    ``numpy.polynomial.polynomial.polyfit``, i.e. an ordinary unweighted
    least-squares straight line through the selected points.

Censoring and sentinels
-----------------------

Three "missing value" conventions appear in the data and it matters which
statistics see them.

**Fibers that never degrade.**  In the replayed state
:math:`\Sigma^{(s)}` these hold :math:`\infty`.  Since every use is a
comparison of the form :math:`\Sigma \le t` or :math:`X < t_{\text{end}}`, they
are excluded, never counted, and never contribute a numeric value.

**The Fortran sentinel** :math:`9.9\times10^{100}`.  This is the on-disk marker
for "fiber not yet scheduled for degradation" in the pre-HDF5 specifications
(≤ v1.90.0, see :doc:`data_specification`); it is defined as
``CONST.UNSCHEDULED_DEGRADE_TIME`` in ``src/lysis/config/constants.py``.  The
v2.0.0 files these commands read store an *event log*, not a snapshot array, so
the sentinel is normally absent.  Where a converted file does carry it as a
degrade time, it behaves exactly like :math:`\infty` for every statistic on
this page — it is larger than any save-point time, so it fails every
:math:`\le t` test.

**Molecules that never reach the back row.**  These hold exactly :math:`0` in
``tpa_transit_time``.  Both statistics derived from that dataset filter on
:math:`\tau > 0`, so they are excluded — see the two formulas below for exactly
what that does to each denominator.

.. warning::

   ``tpa_leaving_time`` is specified as :math:`\infty` when the tPA molecule is
   still bound at the end of a microscale Simulation.  The microscale
   statistics apply **no filter** to that dataset.  A Run containing even one
   infinite value will therefore report an infinite mean and standard
   deviation for tPA leaving time.  (No such value occurs in the austin-runs
   data; the median would remain finite in any case.)

Units
-----

Everything printed by ``micro-stats``, ``macro-stats`` and ``compare`` is a
plain Python ``float`` or ``int``, **not** a Pint ``Quantity``.  The units are
carried only in the column heading.  ``lysis parameters`` is the exception: it
reads Pint ``Quantity`` values from the file and strips ``.magnitude`` at the
last moment, after converting to the display unit.

The model mixes seconds and minutes.  Times read from HDF5 are in **seconds**;
several statistics are divided by 60 before reporting.  Each formula below
states the unit of every quantity.

.. _statistics_cli_micro:

``lysis micro-stats``
=====================

Synopsis
--------

.. code-block:: text

    lysis micro-stats [OPTIONS] PATH

Reports fiber-degradation counts and lysis times together with tPA leaving
times, aggregated over all microscale Simulations in each Run.  Requires only
a ``microscale_out`` Data Collection, so it works on microscale-only files.

Options
-------

.. program:: lysis micro-stats

.. option:: --sort <smart|alpha>

   Row order in directory mode.  Default ``smart``.  See
   :ref:`statistics_cli_sort`.

.. option:: --no-progress

   Suppress progress indicators.

.. option:: --markdown <FILE>

   Write a Markdown table to ``FILE``, or to standard output when ``FILE`` is
   ``-``.  Implies ``--no-progress``.

.. option:: --help

   Show the command's help and exit.

Statistics
----------

The five columns, in the order printed.  They are named by
``MICRO_STATS_COLUMNS`` in ``src/lysis/analysis/summary.py``; the numbers come
from ``compute_micro_statistics`` in ``src/lysis/analysis/microscale.py``.

.. _stat-fibers-degraded:

``Fibers Degraded``
~~~~~~~~~~~~~~~~~~~

The number of microscale Simulations that ended with the fiber fully degraded.
Dimensionless count, printed with thousands separators.

.. math::

   N_D = \sum_{k=0}^{N_\mu - 1} d_k

Aggregated over: **all** :math:`N_\mu` microscale Simulations.  Nothing is
excluded.  Note that this is a count of *Simulations*, and each microscale
Simulation models one fiber — so "fibers degraded" and "simulations in which
the fiber degraded" are the same number.

``Mean Lysis Time (min)``
~~~~~~~~~~~~~~~~~~~~~~~~~

Printed as ``mean ± sd``, both in **minutes**, both to three decimals.

.. math::

   \overline{\ell} = \frac{1}{60\,N_D} \sum_{k \in \mathcal{D}} T_k
   \qquad
   \operatorname{sd}(\ell)
     = \frac{1}{60}\sqrt{\frac{1}{N_D}\sum_{k \in \mathcal{D}}
       \bigl(T_k - \textstyle\frac{1}{N_D}\sum_{k' \in \mathcal{D}} T_{k'}\bigr)^2}

Aggregated over: **only the degraded Simulations** :math:`\mathcal{D}`.  This
is a *conditional* mean — the mean lysis time given that lysis occurred.
Simulations whose fiber did not degrade contribute nothing, neither a value nor
a denominator.  Comparing this number across Runs with very different
``Fibers Degraded`` compares two conditional means over different
sub-populations; that is a real confound, not a subtlety.

If :math:`N_D = 0` all three lysis-time statistics are ``nan``.

Standard deviation is the population form (:math:`\text{ddof} = 0`).

``Median Lysis Time (min)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \tilde{\ell} = \frac{1}{60}\,\operatorname{median}_{k \in \mathcal{D}} T_k

Same domain as the mean: degraded Simulations only.  ``numpy.median``, linear
interpolation for even :math:`N_D`.

``Mean tPA Leaving Time (sec)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Printed as ``mean ± sd``, both in **seconds**.

.. math::

   \overline{L} = \frac{1}{N_\mu} \sum_{k=0}^{N_\mu - 1} L_k
   \qquad
   \operatorname{sd}(L) = \sqrt{\frac{1}{N_\mu}\sum_{k=0}^{N_\mu-1}(L_k - \overline{L})^2}

Aggregated over: **all** :math:`N_\mu` Simulations, degraded or not.  This is
the one place where the domain differs from the lysis-time columns — do not
assume the two ``Mean ...`` columns are means over the same set.  No unit
conversion is applied: the value is reported in seconds exactly as stored.

``Median tPA Leaving Time (sec)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \tilde{L} = \operatorname{median}_{0 \le k < N_\mu} L_k

All Simulations, in seconds.

Example
-------

.. code-block:: console

    $ lysis micro-stats /shared/lysis-group/wpumphrey/austin-runs-corrected/2131.h5
    2131

      Fibers Degraded: 46,685
      Mean Lysis Time (min): 8.891 ± 4.039
      Median Lysis Time (min): 8.156
      Mean tPA Leaving Time (sec): 192.342 ± 190.910
      Median tPA Leaving Time (sec): 134.165

.. _statistics_cli_macro:

``lysis macro-stats``
=====================

Synopsis
--------

.. code-block:: text

    lysis macro-stats [OPTIONS] PATH

Reports six macroscale summary metrics.  Requires a ``macroscale_out`` Data
Collection; on a microscale-only file it fails with
``'NoneType' object has no attribute 'macro_simulations'``.

Options
-------

.. program:: lysis macro-stats

.. option:: --sort <smart|alpha>

   Row order in directory mode.  Default ``smart``.

.. option:: --no-progress

   Suppress progress indicators.

.. option:: --markdown <FILE>

   Write a Markdown table to ``FILE``, or ``-`` for standard output.  In
   single-file mode the Markdown table has two columns, ``Metric`` and
   ``Value``, with the mean and standard deviation combined into one
   ``mean ± std`` cell.

.. option:: --help

   Show the command's help and exit.

Derived quantities
------------------

Two intermediate quantities are used by several of the six metrics.  Both are
defined in ``src/lysis/analysis/degradation.py``.

.. _macro-degraded-fraction:

Degraded fraction :math:`\phi^{(s)}_k`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fraction of fibrin fibers degraded in Simulation :math:`s` at save point
:math:`k`.  Dimensionless, in :math:`[0, 1]`.

.. math::

   \phi^{(s)}_k =
   \frac{\Bigl|\bigl\{(i,j) : \Sigma^{(s)}\!\bigl(t^{(s)}_k\bigr)_{ij}
          \le t^{(s)}_k \bigr\}\Bigr| \;-\; WE}{F}

The count runs over the whole :math:`R \times W` state array.  The
:math:`WE` fibrin-free edges are initialised to :math:`0` and therefore always
satisfy the test; subtracting :math:`WE` removes them, leaving only genuine
fibrin fibers in the numerator.  Padding cells hold NaN and fail the test.
Fibers that are not yet scheduled hold :math:`\infty` and fail the test.

Source of truth: ``find_degraded_fraction`` in
``src/lysis/analysis/degradation.py``.

.. _macro-exposure-time:

Row exposure time :math:`X^{(s)}_{i,j}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The time in **minutes** at which grid row :math:`i` in column :math:`j` first
becomes exposed — that is, the moment the fiber above it has degraded — using
the final replayed state :math:`\Sigma^{(s)}_{\text{end}} =
\Sigma^{(s)}(t^{(s)}_{K_s-1})`.

.. math::

   X^{(s)}_{0,j} = 0,
   \qquad
   X^{(s)}_{i,j} = \frac{1}{60}\max\Bigl(
       60\,X^{(s)}_{i-1,j},\;
       \bigl(\Sigma^{(s)}_{\text{end}}\bigr)_{i,j}
   \Bigr),
   \quad 1 \le i \le R-2,\ 0 \le j < C

(The code accumulates the maximum in seconds and divides the whole array by 60
at the end; the nesting above says the same thing.)

.. note::

   The column index :math:`j` here runs over the **first** :math:`C` **ranks**
   of each grid row, not over :math:`C` edges of one orientation.  Within a row
   the ranks interleave orientations — y-edges at :math:`j \bmod 3 = 0`,
   z-edges at :math:`j \bmod 3 = 1`, x-edges at :math:`j \bmod 3 = 2`
   (``src/lysis/geometry/edge_grid.py``).  So the "columns" over which the
   front velocity is averaged are the leftmost :math:`C` edges of each row and
   are of mixed orientation.  This is what the code does; whether it is what
   the model intends has not been established.  The question is filed for
   analysis as
   `issue #127 <https://github.com/UCO-OpResearch/lysis/issues/127>`_ — it is
   an open question about intent, not a confirmed defect.

Source of truth: ``calculate_time_row_exposed`` in
``src/lysis/analysis/degradation.py``.

Statistics
----------

Each of the six metrics is reported as a ``Mean ± Standard Deviation`` pair,
formatted to three decimals with thousands separators.  The pair is built by
``compute_run_statistics`` in ``src/lysis/analysis/degradation.py`` and
rendered by ``macro_stats_table`` in ``src/lysis/analysis/summary.py``.

``Degradation rate (%/min)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The slope of a straight line fitted through the **rapid-degradation phase** of
each Simulation, in **percent of the clot per minute**.

First define the per-save-point increment and the rapid-phase index set for
Simulation :math:`s`:

.. math::

   \delta^{(s)}_0 = \phi^{(s)}_0, \qquad
   \delta^{(s)}_k = \phi^{(s)}_k - \phi^{(s)}_{k-1} \ \ (k \ge 1)

.. math::

   \mathcal{R}^{(s)} = \Bigl\{\, k \;:\; \delta^{(s)}_k >
       \tfrac{1}{2}\max_{0 \le j < K_s} \delta^{(s)}_j \,\Bigr\}

Then :math:`m_s` is the least-squares slope of :math:`\phi^{(s)}_k` regressed
on :math:`t^{(s)}_k / 60` over :math:`k \in \mathcal{R}^{(s)}`, in
fraction per minute:

.. math::

   m_s = \operatorname*{arg\,min}_{m}\ \min_{b}
     \sum_{k \in \mathcal{R}^{(s)}}
     \Bigl(\phi^{(s)}_k - b - m\,\tfrac{t^{(s)}_k}{60}\Bigr)^{2}

.. math::

   \text{Mean} = 100 \cdot \operatorname{mean}_{s}(m_s),
   \qquad
   \text{Standard Deviation} = 100 \cdot \operatorname{sd}_{s}(m_s)

Aggregated over: the :math:`S` macroscale Simulations.  The rapid-phase
threshold is re-derived independently for each Simulation, so the fitted window
may cover a different number of save points in each.

.. _stat-lysis-lag:

``Lysis lag time (min)``
~~~~~~~~~~~~~~~~~~~~~~~~

The time of the **first** save point in the rapid-degradation phase, in
**minutes**.

.. math::

   g_s = \frac{t^{(s)}_{k^*_s}}{60},
   \qquad k^*_s = \min \mathcal{R}^{(s)}

.. math::

   \text{Mean} = \operatorname{mean}_{s}(g_s),
   \qquad
   \text{Standard Deviation} = \operatorname{sd}_{s}(g_s)

Aggregated over: the :math:`S` Simulations.  Because :math:`g_s` is a
save-point time and not an interpolated crossing, its resolution is the save
interval.

``Time to full clot degradation (min)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The time of the first save point at which the degraded fraction reaches
:math:`1`, in **minutes**.

.. math::

   u_s = \frac{t^{(s)}_{c_s}}{60},
   \qquad c_s = \min\bigl\{\, k : \phi^{(s)}_k \ge 1 \,\bigr\}

.. math::

   \text{Mean} = \operatorname{mean}_{s}(u_s),
   \qquad
   \text{Standard Deviation} = \operatorname{sd}_{s}(u_s)

.. warning::

   **This statistic is not censored — it is silently wrong when the clot never
   fully degrades.**  The index :math:`c_s` is computed with ``numpy.argmax``
   over the boolean array :math:`\phi^{(s)} \ge 1`.  When no element is true,
   ``argmax`` returns ``0``, so :math:`c_s = 0` and the reported time becomes
   :math:`t^{(s)}_0 / 60`, which is normally ``0.0`` minutes — indistinguishable
   from instantaneous total lysis.  Before trusting this column, confirm that
   the clot did fully degrade in every Simulation.  A near-zero mean, or a mean
   far below the ``Lysis lag time``, is the tell.

   This is a confirmed defect and is tracked as
   `issue #126 <https://github.com/UCO-OpResearch/lysis/issues/126>`_, which
   the maintainer intends to fix.  Until it is fixed, check for full
   degradation yourself.

The same ``argmax`` rule governs all degradation milestones, not only 100 %,
so the same defect reaches every milestone of :ref:`statistics_cli_deg_time`
and both endpoints of every interval of :ref:`statistics_cli_deg_rate`.  See
``find_degradation_marker_frames`` in ``src/lysis/analysis/degradation.py``.

``Percent of molecules that reached the back row``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The percentage of tPA molecules whose first-passage time is strictly positive,
computed **per Simulation** and then aggregated.  Dimensionless percent.

.. math::

   p_s = \frac{100}{M}\;\Bigl|\bigl\{\, m : \tau^{(s)}_m > 0 \,\bigr\}\Bigr|

.. math::

   \text{Mean} = \operatorname{mean}_{s}(p_s),
   \qquad
   \text{Standard Deviation} = \operatorname{sd}_{s}(p_s)

Aggregated over: the :math:`S` Simulations.  The denominator is the constant
``total_molecules`` :math:`M` taken from the Run's parameters, not the length
of the stored array, so a Run whose ``total_molecules`` disagrees with its
stored ``tpa_transit_time`` length will report a distorted percentage.

``First passage time (min)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The mean and standard deviation of the transit times of molecules that *did*
reach the back row, in **minutes**.

.. math::

   \mathcal{F} = \Bigl\{\, \tfrac{\tau^{(s)}_m}{60}
      \;:\; 0 \le s < S,\ 0 \le m < M,\ \tau^{(s)}_m > 0 \,\Bigr\}

.. math::

   \text{Mean} = \operatorname{mean}\bigl(\mathcal{F}\bigr),
   \qquad
   \text{Standard Deviation} = \operatorname{sd}\bigl(\mathcal{F}\bigr)

Aggregated over: **every molecule of every Simulation, pooled**, not
Simulation-by-Simulation.  This is the only macroscale metric whose mean is
taken over molecules rather than over Simulations, so its :math:`n` is
:math:`|\mathcal{F}| \approx S \cdot M \cdot p/100`, not :math:`S`.  Its
standard deviation therefore measures molecule-to-molecule spread, and is not
comparable with the across-Simulation standard deviations of the other five
metrics.  Molecules that never reached the back row (:math:`\tau = 0`) are
excluded from both the numerator and the count — this is a **conditional**
mean, given arrival.

``Front Velocity (microns/min)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The speed at which the lysis front advances into the clot, in **microns per
minute**.

For Simulation :math:`s` and column :math:`j`, collect the front events — the
rows that became exposed strictly later than the row above them and strictly
before the end of the Simulation:

.. math::

   \mathcal{E}^{(s)}_j = \Bigl\{\, \bigl(X^{(s)}_{i,j},\, y_i\bigr)
     \;:\; 1 \le i \le R-2,\;\;
     X^{(s)}_{i-1,j} < X^{(s)}_{i,j} < \tfrac{t^{(s)}_{K_s-1}}{60} + 1 \,\Bigr\}

The upper bound is the censoring rule: a row whose exposure time exceeds the
Simulation's end time by more than one minute — including every row whose fiber
never degraded, whose exposure time is :math:`\infty` — is dropped from the
fit.  The ``+ 1`` is a literal one-minute slack in the code.

:math:`v_{s,j}` is the least-squares slope of :math:`y` on exposure time over
:math:`\mathcal{E}^{(s)}_j`, in microns per minute.  Then

.. math::

   \mu_s = \frac{1}{C}\sum_{j=0}^{C-1} v_{s,j},
   \qquad
   \sigma_s = \operatorname{sd}_{j}\bigl(v_{s,j}\bigr)

.. math::

   \text{Mean} = \operatorname{mean}_{s}(\mu_s),
   \qquad
   \text{Standard Deviation} = \operatorname{mean}_{s}(\sigma_s)

.. warning::

   **The reported "Standard Deviation" for front velocity is not a standard
   deviation across Simulations.**  It is the *mean over Simulations of the
   within-Simulation standard deviation across columns* — a measure of how much
   the front velocity varies from column to column inside a typical
   Simulation, not of how much the Run-level velocity varies between
   Simulations.  The two differ substantially in practice.  If you want the
   across-Simulation spread, compute :math:`\operatorname{sd}_s(\mu_s)`
   yourself from
   ``per_sim_front_velocity``
   (``src/lysis/analysis/degradation.py``), which returns the :math:`(\mu_s)` and :math:`(\sigma_s)` vectors directly.

   This aggregation is **deliberate for now** — it is not scheduled for change,
   and the formula and the numbers above are correct as the command stands.  It
   is tracked as
   `issue #124 <https://github.com/UCO-OpResearch/lysis/issues/124>`_, together
   with the ``.. todo::`` in ``mean_front_velocity``
   (``src/lysis/analysis/degradation.py``) proposing a single mean and standard
   deviation pooled over all columns of all Simulations.  Read that issue
   before quoting this column's ``±`` as a Run-level uncertainty.

Example
-------

.. code-block:: console

    $ lysis macro-stats /shared/lysis-group/bpaynter/data/lysis-front-pre-lat/TF-x__9_951.h5
    TF-x__9_951

      Degradation rate (%/min): 2.946 ± 0.065
      Lysis lag time (min): 3.433 ± 0.978
      Time to full clot degradation (min): 42.133 ± 0.521
      Percent of molecules that reached the back row: 21.321 ± 0.485
      First passage time (min): 33.426 ± 8.369
      Front Velocity (microns/min): 3.471 ± 0.183

.. _statistics_cli_deg_rate:

``lysis deg-rate``
==================

Synopsis
--------

.. code-block:: text

    lysis deg-rate [OPTIONS] PATH

Reports the macroscale degradation rate over one or more **configurable**
degradation intervals — by default 20 %→80 %, 20 %→50 % and 50 %→80 %.  Like
:ref:`statistics_cli_macro` it requires a ``macroscale_out`` Data Collection and
fails on a microscale-only file.

Options
-------

.. program:: lysis deg-rate

.. option:: --sort <smart|alpha>

   Row order in directory mode.  Default ``smart``.  See
   :ref:`statistics_cli_sort`.

.. option:: --no-progress

   Suppress progress indicators.

.. option:: --markdown <FILE>

   Write a Markdown table to ``FILE``, or ``-`` for standard output.  Implies
   ``--no-progress``.  Markdown column headers carry the units inline, as
   ``20% to 80% (%/min)``.

.. option:: --add <START-END>

   Add one degradation interval.  Repeatable.  ``START`` and ``END`` are
   integer percentages written with a hyphen between them, e.g. ``--add 0-100``.
   They must satisfy :math:`0 \le \text{START} < \text{END} \le 100`.

.. option:: --drop <START-END>

   Remove one interval from the default set.  Repeatable, same ``START-END``
   syntax, e.g. ``--drop 20-80``.

.. option:: --help

   Show the command's help and exit.

Choosing intervals
------------------

The effective interval list is built as **drop first, then add**
(``src/lysis/cli/deg_rate.py``):

1. start from the defaults ``[(20, 80), (20, 50), (50, 80)]``;
2. remove every interval named by ``--drop``;
3. append every interval named by ``--add`` that is not already present.

Three consequences, all confirmed by executing the command:

* **Column order is not sorted.**  Surviving defaults keep their default order,
  and added intervals follow in the order you wrote them.  ``--add 0-100 --add
  0-50 --drop 20-80`` yields the columns ``20% to 50%``, ``50% to 80%``,
  ``0% to 100%``, ``0% to 50%`` — in that order.
* **``--add`` beats ``--drop`` for the same interval.**  Dropping is applied to
  the defaults before adding, so ``--add 20-80 --drop 20-80`` leaves
  ``20% to 80%`` in the table.
* **Duplicates are ignored.**  ``--add 20-80`` when it is already present is a
  no-op.

An interval that fails to parse, or that violates
:math:`0 \le \text{START} < \text{END} \le 100`, is reported and the command
exits ``1`` — not Click's usual ``2``, because the value is validated inside the
command body rather than by a Click parameter type.  Dropping every default
without adding anything also exits ``1``
(``Error: No intervals remain after applying --drop.``).

Statistics
----------

One column per interval.  The DataFrame column name is exactly
``"<START>% to <END>%"``; the Rich table appends ``(%/min)`` on a second header
line and the Markdown table appends it inline.  Each cell is
``mean ± std`` formatted to **four** decimal places (``deg_rate_table`` in
``src/lysis/analysis/summary.py``).

``<START>% to <END>%``
~~~~~~~~~~~~~~~~~~~~~~

The average degradation rate between two milestones, in **percent of the clot
per minute**.  Plain Python floats, not Pint ``Quantity`` values.

For an interval given as integer percentages :math:`(a, b)`, write
:math:`\theta_a = a/100` and :math:`\theta_b = b/100`.  The per-Simulation rate
is the **two-point secant slope** of the degraded fraction between the two
milestone frames:

.. math::

   r_s(a, b) = 100 \cdot 60 \cdot
     \frac{\phi^{(s)}_{k_s(\theta_b)} - \phi^{(s)}_{k_s(\theta_a)}}
          {t^{(s)}_{k_s(\theta_b)} - t^{(s)}_{k_s(\theta_a)}}

The factor :math:`60` converts the per-second denominator to per-minute; the
factor :math:`100` converts the dimensionless fraction to percent.  The
reported pair is

.. math::

   \text{Mean} = \operatorname{mean}_{s}\bigl(r_s(a, b)\bigr),
   \qquad
   \text{Standard Deviation} = \operatorname{sd}_{s}\bigl(r_s(a, b)\bigr)

Aggregated over: the :math:`S` macroscale Simulations, all of them, with no
exclusions.  Standard deviation is the population form
(:math:`\text{ddof} = 0`); no quantile is involved.

Because the endpoints are save-point frames rather than interpolated crossings,
:math:`r_s` is a secant between two points that both lie *on or after* their
thresholds — its time resolution is the save interval, and it is systematically
a slight underestimate of the instantaneous rate at the interval's midpoint.

Source of truth: ``degradation_rates`` and ``compute_degradation_rate_stats``
in ``src/lysis/analysis/degradation.py``; interval parsing and assembly in
``src/lysis/cli/deg_rate.py``.

.. warning::

   **The** `#126 <https://github.com/UCO-OpResearch/lysis/issues/126>`_
   **milestone defect reaches this command through both endpoints.**  When a
   milestone is never reached, :math:`k_s(\theta)` is ``0`` rather than
   undefined, and the secant is silently taken from the wrong frames:

   * **End milestone never reached** — :math:`k_s(\theta_b) = 0`, so both the
     numerator and the denominator change sign and cancel.  The result is a
     **positive, entirely plausible-looking rate** that is in fact the
     :math:`0\% \rightarrow a\%` rate, not the :math:`a\% \rightarrow b\%` rate
     you asked for.  Nothing in the output marks it.
   * **Both milestones never reached** — :math:`k_s(\theta_a) = k_s(\theta_b) =
     0`, the denominator is exactly zero, and the cell reports ``nan``.

   Only the second case is visible in the table.  Before trusting a
   ``deg-rate`` column, confirm with :ref:`statistics_cli_deg_time` that every
   Simulation actually reached **both** endpoints of the interval — a
   milestone time of ``0.00`` there is the tell.

Relationship to ``macro-stats``
-------------------------------

``lysis deg-rate`` and the ``Degradation rate (%/min)`` metric of
:ref:`statistics_cli_macro` both report a percent-per-minute degradation rate,
and they are **not the same quantity**.  Neither reproduces the other for any
choice of interval:

.. list-table::
   :header-rows: 1
   :widths: 22 39 39

   * -
     - ``macro-stats`` ``Degradation rate``
     - ``deg-rate`` ``<a>% to <b>%``
   * - Estimator
     - Least-squares slope of a straight line fitted to *all* save points in
       the rapid-degradation phase.
     - Two-point secant between the :math:`a\%` and :math:`b\%` milestone
       frames.
   * - Window
     - Chosen automatically per Simulation: the save points where the
       incremental change exceeds half that Simulation's maximum increment.
     - Chosen by you, as a pair of degradation percentages.
   * - Points used
     - Every save point in the window.
     - Exactly two.
   * - Configurable
     - No.
     - Yes.

On ``TF-x__9_951`` the two give ``2.946 ± 0.065`` (``macro-stats``) and
``3.0033 ± 0.0691`` (``deg-rate``, 20 %→80 %).  Both are correct; they measure
different things.  Use ``macro-stats`` when you want one robust
whole-Simulation rate, and ``deg-rate`` when you want to compare the early and
late phases of lysis against each other — which is exactly what the default
20 %→50 % and 50 %→80 % pair is for.

Example
-------

.. code-block:: console

    $ lysis deg-rate /shared/lysis-group/bpaynter/data/lysis-front-pre-lat/TF-x__9_951.h5
    TF-x__9_951

      20% to 80%: 3.0033 ± 0.0691 %/min
      20% to 50%: 2.9692 ± 0.0726 %/min
      50% to 80%: 3.0391 ± 0.0811 %/min

.. _statistics_cli_deg_time:

``lysis deg-time``
==================

Synopsis
--------

.. code-block:: text

    lysis deg-time [OPTIONS] PATH

Reports how long the clot took to reach each of one or more **configurable**
degradation milestones — by default 5 %, 20 %, 50 %, 80 % and 100 %.  Requires a
``macroscale_out`` Data Collection.

Options
-------

.. program:: lysis deg-time

.. option:: --sort <smart|alpha>

   Row order in directory mode.  Default ``smart``.

.. option:: --no-progress

   Suppress progress indicators.

.. option:: --markdown <FILE>

   Write a Markdown table to ``FILE``, or ``-`` for standard output.  Implies
   ``--no-progress``.  Markdown column headers carry the units inline, as
   ``50% (min)``.

.. option:: --add <PCT>

   Add one degradation milestone.  Repeatable.  ``PCT`` is an integer
   percentage in :math:`[0, 100]`, with an optional trailing ``%`` — both
   ``--add 90`` and ``--add 90%`` work.

.. option:: --drop <PCT>

   Remove one milestone from the default set.  Repeatable, same syntax.

.. option:: --help

   Show the command's help and exit.

Choosing milestones
-------------------

Assembly follows the same **drop first, then add** rule as ``deg-rate``
(``src/lysis/cli/deg_time.py``), with one difference that matters:

* **The final list is sorted ascending**, unlike ``deg-rate``'s.  ``--add 90%
  --add 10 --drop 5`` yields the columns ``10%``, ``20%``, ``50%``, ``80%``,
  ``90%``, ``100%`` — your ``--add`` order is not preserved.
* ``--add`` still beats ``--drop`` for the same value, and duplicates are still
  ignored.
* A milestone that is not an integer in :math:`[0, 100]` is reported and the
  command exits ``1``.  Dropping every default without adding anything also
  exits ``1`` (``Error: No milestones remain after applying --drop.``).

.. note::

   ``--add 0`` is accepted and produces a column that is always
   ``0.00 ± 0.00``: every Simulation satisfies :math:`\phi^{(s)}_0 \ge 0` at its
   first save point, so :math:`k_s(0) = 0` and :math:`T_s(0) = t^{(s)}_0 / 60`,
   which is zero.  That zero is legitimate arithmetic, **not** the ``#126``
   defect — but it is also the same value ``#126`` produces for an unreached
   milestone, which is precisely why an unreached milestone is hard to spot.

Statistics
----------

One column per milestone.  The DataFrame column name is exactly ``"<PCT>%"``;
the Rich table appends ``(min)`` on a second header line and the Markdown table
appends it inline.  Each cell is ``mean ± std`` formatted to **two** decimal
places (``deg_time_table`` in ``src/lysis/analysis/summary.py``).

``<PCT>%``
~~~~~~~~~~

The time at which the clot first reached the given degradation milestone, in
**minutes**.  Plain Python floats, not Pint ``Quantity`` values.

For a milestone given as an integer percentage :math:`m`, with
:math:`\theta = m/100`:

.. math::

   \text{Mean} = \operatorname{mean}_{s}\bigl(T_s(\theta)\bigr),
   \qquad
   \text{Standard Deviation} = \operatorname{sd}_{s}\bigl(T_s(\theta)\bigr)

where :math:`T_s(\theta) = t^{(s)}_{k_s(\theta)} / 60` is the milestone time of
:ref:`statistics_cli_notation`.

Aggregated over: the :math:`S` macroscale Simulations, all of them, with no
exclusions.  Standard deviation is the population form
(:math:`\text{ddof} = 0`); no quantile is involved.

Because :math:`T_s` is a save-point time and not an interpolated crossing, its
resolution is the save interval and it is systematically a slight
**over**\ estimate of the true crossing time.

Source of truth: ``compute_degradation_marker_stats`` in
``src/lysis/analysis/degradation.py``; milestone parsing and assembly in
``src/lysis/cli/deg_time.py``.

.. warning::

   **Every column of this command is exposed to the**
   `#126 <https://github.com/UCO-OpResearch/lysis/issues/126>`_ **milestone
   defect.**  A milestone that a Simulation never reaches is reported as having
   been reached at save point 0 — normally ``0.00`` minutes — rather than as
   missing.  A single such Simulation drags that column's mean down and inflates
   its standard deviation, with nothing in the output to say so.

   The high milestones are the ones at risk: a Run that stalls at 90 % lysis
   will show a plausible ``80%`` column and a badly wrong ``100%`` one.  Look
   for a column whose mean is *lower* than the column to its left, or whose
   standard deviation is large relative to its neighbours — milestone times must
   increase monotonically across a row, so any inversion is a positive
   indication of the defect.

Relationship to ``macro-stats``
-------------------------------

Unlike ``deg-rate``, this command **does** reproduce one of the
:ref:`statistics_cli_macro` metrics exactly.  ``macro-stats``'
``Time to full clot degradation (min)`` is defined as
:math:`u_s = t^{(s)}_{c_s}/60` with :math:`c_s = \min\{k : \phi^{(s)}_k \ge 1\}`
— which is precisely :math:`T_s(1)`, the ``100%`` column here.  The two go
through the same ``find_degradation_marker_times`` call and differ only in
display precision: ``macro-stats`` prints three decimals, ``deg-time`` two.

On ``TF-x__9_951``, ``macro-stats`` reports ``42.133 ± 0.521`` and ``deg-time``
reports ``42.13 ± 0.52`` for ``100%``.  They also share the ``#126`` exposure,
for the same reason.

The other four default milestones have no ``macro-stats`` equivalent — 5 %,
20 %, 50 % and 80 % are only available here.

Example
-------

.. code-block:: console

    $ lysis deg-time /shared/lysis-group/bpaynter/data/lysis-front-pre-lat/TF-x__9_951.h5
    TF-x__9_951

      5%: 4.25 ± 0.08 min
      20%: 10.75 ± 0.23 min
      50%: 20.87 ± 0.32 min
      80%: 30.72 ± 0.52 min
      100%: 42.13 ± 0.52 min

.. _statistics_cli_parameters:

``lysis parameters``
====================

Synopsis
--------

.. code-block:: text

    lysis parameters [OPTIONS] PATH

``lysis parameters`` is a **single command, not a command group** — it has
options, not subcommands.  It prints the Scenario (see :doc:`ontology`) behind
each Run: the numeric parameters stored in the file's ``micro_params`` and
``macro_params``.

Options
-------

.. program:: lysis parameters

.. option:: --sort <smart|alpha>

   Column order in directory mode (Runs are columns here, parameters are rows).
   Default ``smart``.

.. option:: --no-progress

   Suppress progress indicators.

.. option:: --add <PARAM>

   Add one parameter row to the table.  Repeatable.  ``PARAM`` must be an
   attribute name on
   ``MacroParameters`` or ``MicroParameters``
   (``src/lysis/config/parameters.py``).  Added rows appear in an
   ``Additional`` section below the curated table.

.. option:: --drop <PARAM>

   Remove one parameter row from the curated default table.  Repeatable.

.. option:: --markdown <FILE>

   Write the table as Markdown to ``FILE``, or ``-`` for standard output.
   Markdown output is unstyled: the non-default highlighting described below is
   **not** applied.  Implies ``--no-progress``.

.. option:: --csv <FILE>

   Export a re-feedable parameter CSV to ``FILE``, or ``-`` for standard
   output.  This takes a completely different path through the code — see
   `The --csv export`_ below.

.. option:: --help

   Show the command's help and exit.

The curated table
-----------------

Rows are not "every attribute of the parameter objects".  They come from a
hand-maintained list, ``_DEFAULT_PARAMS`` in ``src/lysis/cli/parameters.py``,
whose entries are 4-tuples ``(attr_name, source, display_units, fmt)``:

``attr_name``
    The attribute to read.

``source``
    ``"macro"``, ``"micro"`` or ``"computed"``.  Informational only, except
    that ``"computed"`` rows are dropped alongside the macroscale rows when no
    Run has macroscale data.  There is exactly one computed parameter,
    ``grid_node_distance`` :math:`= \texttt{pore\_size} + 2 \cdot
    \texttt{fiber\_radius}`, which is not stored in the file.

``display_units``
    A Pint unit string to convert to before the magnitude is taken (e.g.
    ``pore_size`` is shown in microns, ``fiber_radius`` in nanometers,
    ``total_time`` in minutes), or ``None`` to use the value's natural units.

``fmt``
    A Python format string applied to the final number, or ``None`` to use
    ``str()``.  Seeds use ``None`` so they print as exact integers.

Row labels combine the attribute name with the display units, or with the
natural units from
``MacroParameters.units()`` when no display unit is given — hence labels like ``pore_size (microns)`` and ``deg_rate_fibrin
(sec^-1)``.  A parameter with no units, such as ``cols``, is labelled by name
alone.

Rows are grouped into a ``Macroscale`` section and a ``Microscale`` section,
plus an ``Additional`` section for ``--add`` rows.  Source of truth:
``parameters_table`` in ``src/lysis/analysis/summary.py``.

Value lookup and formatting
---------------------------

For each row, ``_get_raw`` resolves the value in this order: a computed
parameter, then ``macro_params``, then ``micro_params``; ``None`` if the
attribute exists nowhere or the relevant parameter object is absent.

``_format_value`` then renders it.  A Pint ``Quantity`` is converted with
``.to(display_units)`` when display units are given and its ``.magnitude`` is
taken; a plain number is used as is; the format string is applied last.  A
``None`` value renders as ``N/A``.

``--add`` rows bypass the curated formats and go through ``_auto_format``
instead: a Pint ``Quantity`` prints as ``{:.4g~P}`` (four significant figures
with its own pretty-printed units), an ``int`` with thousands separators, a
``float`` as ``{:.4g}``, anything else via ``str()``.

Microscale-only files
---------------------

If a file has no macroscale Data Collection, ``macro_params`` is ``None`` and
every macroscale row would be ``N/A``.  Rather than print a wall of ``N/A``,
``_drop_macro_specs`` removes the macroscale and computed rows entirely.  The
set of macroscale-only names is derived from the dataclass definitions —
fields on ``MacroParameters`` that are not on ``MicroParameters``, excluding
the nested ``micro_params`` container — so a *microscale* parameter that
happens to be missing still shows as ``N/A`` and is never silently dropped.

In directory mode the rows are dropped only when **no** Run in the directory
has macroscale data, so columns stay aligned across a mixed set.

Non-default highlighting
------------------------

In Rich (terminal) output only, a cell whose value differs from the model
default is printed in a distinct colour, so a scan down a column shows at a
glance what makes a Run different.

The comparison is made on the **formatted string**, not the underlying number.
``_default_formatted`` constructs a default ``MicroParameters`` and a default
``MacroParameters`` and pushes them through the identical formatting path, then
``_nondefault_flags`` marks a cell non-default when its string differs from the
default's string *and* is not ``"N/A"``.  Two consequences worth knowing:

* a difference smaller than the row's display precision is **not** flagged;
* an ``N/A`` is never flagged, so absent macroscale data does not light up.

The flags are attached to the DataFrame as ``df.attrs["nondefault"]``.
Markdown output ignores them, which is why ``--markdown`` skips computing them
at all.

The ``--csv`` export
--------------------

``--csv`` short-circuits everything above.  It does **not** use the curated
table, and ``--add``, ``--drop`` and ``--markdown`` have no effect on it.
Instead it loads each Run's raw ``MicroParameters`` / ``MacroParameters``
objects and hands them to
``write_params_csv``
(``src/lysis/config/parameters.py``), which writes the same
transposed CSV format that ``lysis init-experiment`` and
``lysis init-macroscale`` read back in: column 0 is the shared parameter-name
column, and there is one value column per Run.

Two rules make the output sparse:

* only **independent** parameters are written — dataclass fields with
  ``init=True``, excluding the nested ``micro_params`` container.  Derived
  quantities such as ``total_edges`` and ``full_row`` are omitted because they
  are recomputed on read;
* a cell is left **blank** when the Run's value equals that field's effective
  default, and a row that is blank for every Run is omitted entirely.

So the CSV shows you exactly what distinguishes these Runs from a default
Scenario — which is usually what you want for an Experiment table — and it
round-trips: feeding it back re-derives the same parameters.  Note that values
are written with their Pint units attached, e.g. ``0.001 / micromolar /
second``.

Worked example
--------------

The parameters of one Run, as Markdown:

.. code-block:: console

    $ lysis parameters /shared/lysis-group/bpaynter/data/lysis-front-pre-lat/TF-x__9_951.h5 --markdown -
    ## TF-x__9_951

    ### Macroscale Parameters

    | Parameter | TF-x__9_951 |
    | --- | --- |
    | pore_size (microns) | 2.5000 |
    | diffusion_coeff (cm^2/s) | 5.00000e-07 |
    | forced_unbind | 0.079 |
    | average_bound_time (seconds) | 27.78 |
    | grid_node_distance (microns) | 2.605 |
    | cols | 39 |
    | rows | 385 |
    ...
    | total_edges | 44,621 |
    | total_fibers | 6,805 |
    | total_molecules | 9,951 |
    ...

    ### Microscale Parameters

    | Parameter | TF-x__9_951 |
    | --- | --- |
    | bind_rate_tPA ((micromolar*sec)^-1) | 0.100 |
    ...
    | micro_seed | 4209415086 |

The same Run's non-default parameters as a re-feedable CSV — here for a
microscale-only Run, which is why only four rows survive:

.. code-block:: console

    $ lysis parameters /shared/lysis-group/wpumphrey/austin-runs-corrected/2131.h5 --csv -
    parameter,2131
    bind_rate_tPA,0.001 / micromolar / second
    unbind_rate_PLi,0.576 / second
    nodes_in_micro_row,5
    micro_seed,2616089003

.. _statistics_cli_compare:

``lysis compare micro-stats`` and ``lysis compare macro-stats``
===============================================================

Synopsis
--------

.. code-block:: text

    lysis compare [OPTIONS] {macro-stats|micro-stats} PATH1 PATH2

.. important::

   ``micro-stats`` and ``macro-stats`` here are **not subcommands**.  ``lysis
   compare`` is a single Click command whose first positional argument,
   ``WHICH``, is a choice of *measure set*.  ``lysis compare micro-stats
   --help`` therefore does not exist; use ``lysis compare --help``.  Options
   may be given before or after ``WHICH``.

Both paths must be of the same kind:

**Folder-pair mode** (``PATH1`` and ``PATH2`` are directories)
    Run codes present in both directories are compared pairwise, one table row
    per common run code.  Codes present in only one directory are skipped.

**File-pair mode** (``PATH1`` and ``PATH2`` are ``.h5`` files)
    A single row is produced, labelled ``<basename1> vs <basename2>``.  The two
    files may have the same run code.

Options
-------

.. program:: lysis compare

.. option:: --sort <smart|alpha>

   Order of the common run codes.  Default ``smart``.

.. option:: --no-progress

   Suppress progress indicators.

.. option:: --markdown <FILE>

   Write the table as Markdown to ``FILE``, or ``-`` for standard output.

.. option:: --help

   Show the command's help and exit.

What is compared
----------------

For each Run pair the command produces two kinds of column, in this order:
first every KS column, then every percent-difference column.  Which quantities
those are is fixed by the measure set.

Source of truth: ``MEASURE_EXTRACTORS`` and ``STATS_COMPUTERS`` in
``src/lysis/analysis/compare.py``; column assembly in ``compare_stats_table``
in ``src/lysis/analysis/summary.py``.

Kolmogorov–Smirnov columns
~~~~~~~~~~~~~~~~~~~~~~~~~~

Each entry produces two columns: ``<label> KS`` (the statistic, four decimals)
and ``<label> p-value`` (three significant figures).

For ``micro-stats``, the two samples are **per-Simulation vectors from the
microscale output**:

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Label
     - Sample vector
   * - ``Lysis Time (min)``
     - :math:`\bigl(T_k / 60 : k \in \mathcal{D}\bigr)` — the degraded
       Simulations only, in minutes.  **Length** :math:`N_D` **differs between
       the two Runs.**
   * - ``tPA Leaving Time (sec)``
     - :math:`\bigl(L_k : 0 \le k < N_\mu\bigr)` — all Simulations, in seconds.
       Length :math:`N_\mu`, typically 50 000.

For ``macro-stats``, the two samples are **per-Simulation vectors of length**
:math:`S` (one value per macroscale Simulation, typically 10 — not per fiber
and not per molecule):

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Label
     - Sample vector
   * - ``Degradation rate (%/min)``
     - :math:`(100\,m_s)_{s=0}^{S-1}`, the per-Simulation rapid-phase slopes.
   * - ``Time to full clot degradation (min)``
     - :math:`(u_s)_{s=0}^{S-1}`.
   * - ``Front Velocity (microns/min)``
     - :math:`(\mu_s)_{s=0}^{S-1}`, the per-Simulation **means over columns**
       from ``per_sim_front_velocity``.  The
       within-Simulation spread :math:`\sigma_s` is not tested.

The test
~~~~~~~~

Let :math:`x_1, \dots, x_n` be the sample from ``PATH1`` and
:math:`y_1, \dots, y_m` the sample from ``PATH2``, with empirical
distribution functions :math:`F_n` and :math:`G_m`.  The reported statistic is
the two-sided two-sample Kolmogorov–Smirnov statistic

.. math::

   D_{n,m} = \sup_{x \in \mathbb{R}} \bigl| F_n(x) - G_m(x) \bigr|

testing

.. math::

   H_0: F = G
   \qquad\text{against}\qquad
   H_1: F \neq G

The call site is ``ks_2samp(arr1, measures_2[label])`` in ``compare_runs``,
``src/lysis/analysis/compare.py``.  **No keyword arguments are passed**, so
SciPy's defaults apply:

``alternative="two-sided"``
    The test is two-sided: it detects a difference in either direction, not a
    stochastic ordering.

``method="auto"``
    SciPy uses the **exact** null distribution when
    :math:`\max(n, m) \le 10\,000`, and Smirnov's **asymptotic** approximation
    otherwise (``scipy/stats/_stats_py.py``).  With this project's data that
    means a ``macro-stats`` comparison (:math:`n = m = S \approx 10`) is exact,
    while a ``micro-stats`` comparison (:math:`n, m \approx 50\,000`) is
    asymptotic.  Quote the method when you report a p-value.

Two properties of the KS test are worth restating for this data:

* it is sensitive to any difference in distribution — location, scale, or
  shape — so a significant result does not tell you that the means differ;
* with :math:`n \approx m \approx 50\,000` (the ``micro-stats`` case) it has
  enormous power, and will return a very small p-value for a difference far
  too small to matter scientifically.  Read :math:`D_{n,m}` as the effect
  size and treat the p-value as secondary.

.. warning::

   **No multiple-comparison correction is applied — anywhere.**  A single
   ``lysis compare macro-stats`` over a directory of :math:`n` Runs performs
   :math:`3n` independent KS tests and prints :math:`3n` raw p-values; the
   ``micro-stats`` measure set performs :math:`2n`.  Testing a directory of 20
   Runs on the macroscale set is 60 tests: at :math:`\alpha = 0.05` you should
   expect about three "significant" results even if every pair of Runs is
   identically distributed.

   The command will not correct for this and does not warn you.  Apply a
   correction yourself over the full set of tests you actually ran — e.g.
   Holm–Bonferroni, or Benjamini–Hochberg if you want false-discovery-rate
   control — using ``p.adjust()`` in R or
   ``statsmodels.stats.multitest.multipletests`` in Python.  Decide the family
   of tests *before* you look at the table.

Percent-difference columns
~~~~~~~~~~~~~~~~~~~~~~~~~~

Each entry produces one column, ``<label> % diff``.  These compare **scalar
summary statistics**, not distributions, using the symmetric relative percent
difference

.. math::

   \operatorname{pct}(a, b) = 100 \cdot \frac{b - a}{(a + b)/2}
                            = \frac{200\,(b - a)}{a + b}

where :math:`a` is the value from ``PATH1`` and :math:`b` the value from
``PATH2``.  The sign is preserved, so **positive means PATH2 > PATH1**.  The
symmetric denominator bounds the result to :math:`\pm 200` for same-sign
inputs.  Special cases: :math:`0` when both inputs are exactly zero; ``NaN``
when either input is ``NaN`` or when :math:`a + b = 0`.

Source of truth: ``percent_difference`` in ``src/lysis/analysis/compare.py``.

The compared scalars are:

``micro-stats``
    ``Fibers Degraded``, ``Mean Lysis Time (min)``,
    ``Median Lysis Time (min)``, ``Mean tPA Leaving Time (sec)``,
    ``Median tPA Leaving Time (sec)`` — the five values of
    :ref:`statistics_cli_micro`.  Standard deviations are deliberately not
    compared.

``macro-stats``
    ``Mean Degradation rate (%/min)``,
    ``Mean Time to full clot degradation (min)``,
    ``Mean First passage time (min)``,
    ``Mean Front Velocity (microns/min)`` — the ``Mean`` half of four of the
    six metrics of :ref:`statistics_cli_macro`.  ``Lysis lag time`` and
    ``Percent of molecules that reached the back row`` are not compared.

Percent-difference cells are printed without a trailing ``%`` and without a
leading ``+``, and are aligned on their decimal points down each column,
switching to scientific notation for values below ``0.005``.

Examples
--------

Two individual files:

.. code-block:: console

    $ lysis compare micro-stats \
        /shared/lysis-group/austin_segrest/austin-old-data-imported/2131.h5 \
        /shared/lysis-group/wpumphrey/austin-runs-corrected/2131.h5 --markdown -
    ## 2131.h5 vs 2131.h5

    | Metric | Value |
    | --- | --- |
    | Lysis Time (min) KS | 0.0343 |
    | Lysis Time (min) p-value | 2.95e-24 |
    | tPA Leaving Time (sec) KS | 0.0040 |
    | tPA Leaving Time (sec) p-value | 0.808 |
    | Fibers Degraded % diff | -0.03 |
    | Mean Lysis Time (min) % diff | 3.30 |
    | Median Lysis Time (min) % diff | 3.75 |
    | Mean tPA Leaving Time (sec) % diff | -0.27 |
    | Median tPA Leaving Time (sec) % diff | -0.17 |

Read that as: the lysis-time distributions differ by :math:`D = 0.034` — a
three-percentage-point maximum gap between the two empirical CDFs — and with
:math:`n \approx 46\,000` per sample that is overwhelmingly "significant".
Whether a 3 % shift in mean lysis time matters is a scientific question, not a
statistical one.

Two directories, macroscale:

.. code-block:: console

    $ lysis compare macro-stats \
        /shared/lysis-group/bpaynter/data/lysis-front-pre-lat/ \
        /shared/lysis-group/bpaynter/data/lysis-front-post-lat-single/ --markdown -
    | Run | Degradation rate (%/min) KS | Degradation rate (%/min) p-value | ... |
    | --- | --- | --- | --- |
    | Q1 | 0.6000 | 0.0524 | ... |
    | TB-ix__21_105 | 0.0000 | 1 | ... |
    ...

With :math:`S = 10` Simulations per Run the smallest attainable p-value is far
from zero, so a table like this one distinguishes "identical" (``KS = 0``,
``p = 1``) from "different" but has little power to resolve small effects.

.. _statistics_cli_recipes:

Recipes
=======

All of these were executed against the shared read-only data.  Replace
``<_output_file_>`` with a path in **your own** scratch space — never write
anything under ``/shared/lysis-group/``.

Microscale statistics for one Run
---------------------------------

.. code-block:: console

    $ lysis micro-stats /shared/lysis-group/wpumphrey/austin-runs-corrected/2131.h5

Every Run in a directory, to a Markdown file
--------------------------------------------

.. code-block:: console

    $ lysis micro-stats /shared/lysis-group/wpumphrey/austin-runs-corrected/ \
        --markdown <_output_file_>

Sending the table to a file rather than to ``-`` keeps per-Run error messages
out of the table.

Macroscale statistics for a whole Experiment
--------------------------------------------

.. code-block:: console

    $ lysis macro-stats /shared/lysis-group/bpaynter/data/lysis-front-pre-lat/ \
        --markdown -
    | Run | Degradation rate (%/min) | Lysis lag time (min) | ... |
    | --- | --- | --- | --- |
    | Q1 | 0.907 ± 0.009 | 1.117 ± 0.707 | ... |
    | TB-ix__21_105 | 1.403 ± 0.019 | 1.233 ± 0.692 | ... |
    | TB-xi__21_105 | 1.793 ± 0.024 | 2.483 ± 0.751 | ... |
    | TB-xiii__21_105 | 2.358 ± 0.050 | 3.967 ± 0.862 | ... |
    | TF-v__9_951 | 1.036 ± 0.009 | 1.067 ± 0.260 | ... |
    | TF-vii__9_951 | 1.782 ± 0.014 | 2.917 ± 1.158 | ... |
    | TF-x__9_951 | 2.946 ± 0.065 | 3.433 ± 0.978 | ... |

Note that ``smart`` sorting has put ``TB-ix`` before ``TB-xi`` before
``TB-xiii`` — Roman-numeral order, not alphabetical order.

Compare two versions of the same Experiment
-------------------------------------------

.. code-block:: console

    $ lysis compare micro-stats \
        /shared/lysis-group/austin_segrest/austin-old-data-imported/ \
        /shared/lysis-group/wpumphrey/austin-runs-corrected/ \
        --markdown <_output_file_>

Remember to correct the resulting p-values for multiple comparisons.

Check that every Run actually finished lysing
---------------------------------------------

Do this **before** trusting any milestone-based column, because of
`issue #126 <https://github.com/UCO-OpResearch/lysis/issues/126>`_:

.. code-block:: console

    $ lysis deg-time /shared/lysis-group/bpaynter/data/lysis-front-pre-lat/ \
        --markdown -
    | Run | 5% (min) | 20% (min) | 50% (min) | 80% (min) | 100% (min) |
    | --- | --- | --- | --- | --- | --- |
    | Q1 | 7.28 ± 0.08 | 22.50 ± 0.30 | 54.00 ± 0.60 | 88.45 ± 0.75 | 117.58 ± 0.45 |
    | TB-ix__21_105 | 3.90 ± 0.08 | 11.88 ± 0.15 | 32.35 ± 0.29 | 56.08 ± 0.53 | 76.63 ± 0.72 |
    | TB-xi__21_105 | 4.38 ± 0.11 | 11.35 ± 0.12 | 27.32 ± 0.30 | 45.15 ± 0.44 | 61.30 ± 0.85 |
    | TB-xiii__21_105 | 4.45 ± 0.11 | 10.18 ± 0.16 | 22.37 ± 0.30 | 35.68 ± 0.51 | 48.78 ± 0.49 |
    | TF-v__9_951 | 5.92 ± 0.08 | 18.92 ± 0.19 | 46.42 ± 0.43 | 76.82 ± 0.63 | 102.60 ± 0.93 |
    | TF-vii__9_951 | 5.08 ± 0.08 | 14.38 ± 0.18 | 30.57 ± 0.41 | 47.60 ± 0.28 | 63.72 ± 0.58 |
    | TF-x__9_951 | 4.25 ± 0.08 | 10.75 ± 0.23 | 20.87 ± 0.32 | 30.72 ± 0.52 | 42.13 ± 0.52 |

Read across each row: the times must increase monotonically.  They do here, and
no cell is ``0.00``, so every Simulation of every Run reached 100 % and the
milestone columns — and the ``deg-rate`` and ``macro-stats`` numbers derived
from them — are trustworthy for this folder.

Compare the early and late phases of lysis
------------------------------------------

The default ``deg-rate`` intervals are chosen for exactly this: 20 %→50 % is the
early half, 50 %→80 % the late half, and 20 %→80 % the whole rapid phase.

.. code-block:: console

    $ lysis deg-rate /shared/lysis-group/bpaynter/data/lysis-front-pre-lat/ \
        --markdown -
    | Run | 20% to 80% (%/min) | 20% to 50% (%/min) | 50% to 80% (%/min) |
    | --- | --- | --- | --- |
    | Q1 | 0.9096 ± 0.0094 | 0.9520 ± 0.0152 | 0.8711 ± 0.0125 |
    | TB-ix__21_105 | 1.3564 ± 0.0163 | 1.4665 ± 0.0187 | 1.2617 ± 0.0221 |
    | TB-xi__21_105 | 1.7721 ± 0.0220 | 1.8754 ± 0.0308 | 1.6799 ± 0.0274 |
    | TB-xiii__21_105 | 2.3517 ± 0.0478 | 2.4650 ± 0.0616 | 2.2499 ± 0.0668 |
    | TF-v__9_951 | 1.0362 ± 0.0102 | 1.0915 ± 0.0131 | 0.9863 ± 0.0124 |
    | TF-vii__9_951 | 1.8061 ± 0.0144 | 1.8529 ± 0.0326 | 1.7624 ± 0.0277 |
    | TF-x__9_951 | 3.0033 ± 0.0691 | 2.9692 ± 0.0726 | 3.0391 ± 0.0811 |

To widen the window, add your own interval — the defaults are only defaults:

.. code-block:: console

    $ lysis deg-rate /shared/lysis-group/bpaynter/data/lysis-front-pre-lat/TF-x__9_951.h5 \
        --add 0-100 --drop 20-80 --markdown <_output_file_>

Export parameters for analysis in R
-----------------------------------

.. code-block:: console

    $ lysis parameters /shared/lysis-group/wpumphrey/austin-runs-corrected/ \
        --csv <_output_file_>

Then, in R:

.. code-block:: r

    params <- read.csv("<_output_file_>", check.names = FALSE)
    # column 1 is the parameter name; one further column per Run.
    # A blank cell means "this Run uses the model default".

Because blank means *default*, do not treat a blank as missing data: fill it
from a default Run, or read the full value table with
``lysis parameters ... --markdown`` instead, which prints every curated
parameter for every Run.

Recomputing a statistic yourself
--------------------------------

Every formula on this page can be checked directly against the raw datasets.
For example, to confirm the ``micro-stats`` numbers for one Run:

.. code-block:: python

    import h5py
    import numpy as np

    path = "<_data_path_>"           # e.g. ".../austin-runs-corrected/2131.h5"
    with h5py.File(path, "r") as f:  # always open read-only
        degraded = f["micro_data/fiber_degraded"][:]
        final_time = f["micro_data/sim_final_time"][:]
        leaving = f["micro_data/tpa_leaving_time"][:]

    lysis_times = final_time[degraded] / 60          # minutes, degraded only
    print(int(degraded.sum()))                        # Fibers Degraded
    print(lysis_times.mean(), lysis_times.std())      # Mean Lysis Time ± sd
    print(np.median(lysis_times))                     # Median Lysis Time
    print(leaving.mean(), leaving.std())              # Mean tPA Leaving ± sd

Note ``lysis_times.std()`` uses ``ddof=0``, matching the CLI.

See also
========

* :doc:`ontology` — what a Run, Simulation, Scenario and Mechanism are.
* :doc:`data_specification` — the datasets, shapes, dtypes and units the
  formulas on this page read.
* :doc:`fiber_size_conventions` — how ``pore_size``, ``fiber_radius`` and the
  grid geometry relate.
* :doc:`rstudio_ondemand` — opening this data from RStudio on Open OnDemand.
