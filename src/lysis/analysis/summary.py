"""Summary table functions that return pandas DataFrames for CLI display.

Each function accepts already-computed statistics (from the lower-level analysis
functions) organised by run code and returns a :class:`pandas.DataFrame` with
pre-formatted string values suitable for rendering in Rich tables or Markdown.

Available functions:

- :func:`micro_stats_table` — microscale fiber and tPA statistics
- :func:`macro_stats_table` — macroscale run summary statistics
- :func:`deg_rate_table` — degradation-rate table for configurable intervals
- :func:`deg_time_table` — degradation-time table for configurable milestones
- :func:`parameters_table` — Scenario/Run parameter comparison table
- :func:`compare_stats_table` — 2-sample KS test statistics across Runs
"""

import math

import pandas as pd


# ---------------------------------------------------------------------------
# micro_stats_table
# ---------------------------------------------------------------------------

#: Ordered display column names for :func:`micro_stats_table`.
MICRO_STATS_COLUMNS = [
    "Fibers Degraded",
    "Mean Lysis Time (min)",
    "Median Lysis Time (min)",
    "Mean tPA Leaving Time (sec)",
    "Median tPA Leaving Time (sec)",
]


def _fmt_micro(col: str, stats: dict) -> str:
    """Format one microscale metric from *stats*.

    :param col: Column name from :data:`MICRO_STATS_COLUMNS`.
    :type col: str
    :param stats: Stats dict from
        :func:`~lysis.analysis.microscale.compute_micro_statistics`.
    :type stats: dict
    :return: Formatted string.
    :rtype: str
    :raises KeyError: If *col* is not a recognised metric name.
    """
    if col == "Fibers Degraded":
        return f"{stats['fibers_degraded']:,}"
    if col == "Mean Lysis Time (min)":
        return f"{stats['lysis_time_mean']:.3f} \u00b1 {stats['lysis_time_std']:.3f}"
    if col == "Median Lysis Time (min)":
        return f"{stats['lysis_time_median']:.3f}"
    if col == "Mean tPA Leaving Time (sec)":
        return f"{stats['tpa_leaving_mean']:.3f} \u00b1 {stats['tpa_leaving_std']:.3f}"
    if col == "Median tPA Leaving Time (sec)":
        return f"{stats['tpa_leaving_median']:.3f}"
    raise KeyError(col)


def micro_stats_table(stats_by_run: dict) -> pd.DataFrame:
    """Build a summary DataFrame for microscale statistics.

    :param stats_by_run: Maps run code to a stats dict as returned by
        :func:`~lysis.analysis.microscale.compute_micro_statistics`.
    :type stats_by_run: dict[str, dict]
    :return: DataFrame with run codes as index and pre-formatted string values.
        Columns are defined by :data:`MICRO_STATS_COLUMNS`.
    :rtype: pandas.DataFrame
    """
    rows = {
        rc: {col: _fmt_micro(col, stats) for col in MICRO_STATS_COLUMNS}
        for rc, stats in stats_by_run.items()
    }
    return pd.DataFrame.from_dict(rows, orient="index", columns=MICRO_STATS_COLUMNS)


# ---------------------------------------------------------------------------
# macro_stats_table
# ---------------------------------------------------------------------------


def macro_stats_table(stats_by_run: dict) -> pd.DataFrame:
    """Build a summary DataFrame for macroscale run statistics.

    :param stats_by_run: Maps run code to a :class:`pandas.Series` with a
        two-level MultiIndex ``(metric_name, "Mean"/"Standard Deviation")``
        as returned by
        :func:`~lysis.analysis.degradation.compute_run_statistics`.
    :type stats_by_run: dict[str, pandas.Series]
    :return: DataFrame with run codes as index and one column per metric.
        Values are ``"mean ± std"`` formatted strings (``:.3f`` with commas).
    :rtype: pandas.DataFrame
    """
    if not stats_by_run:
        return pd.DataFrame()
    metrics = (
        next(iter(stats_by_run.values()))
        .index.get_level_values(0)
        .unique()
        .tolist()
    )
    rows = {
        rc: {
            m: f"{series[m]['Mean']:,.3f} \u00b1 {series[m]['Standard Deviation']:,.3f}"
            for m in metrics
        }
        for rc, series in stats_by_run.items()
    }
    return pd.DataFrame.from_dict(rows, orient="index")[metrics]


# ---------------------------------------------------------------------------
# deg_rate_table
# ---------------------------------------------------------------------------


def deg_rate_table(stats_by_run: dict, intervals: list) -> pd.DataFrame:
    """Build a summary DataFrame for degradation rates.

    Column names use the plain ``"s% to e%"`` form (without units), so the
    DataFrame is suitable for programmatic use.  Display layers add units to
    the rendered headers.

    :param stats_by_run: Maps run code to ``{(start%, end%): (mean, std)}`` as
        returned by
        :func:`~lysis.analysis.degradation.compute_degradation_rate_stats`.
    :type stats_by_run: dict[str, dict[tuple[int, int], tuple[float, float]]]
    :param intervals: Ordered list of ``(start, end)`` integer percent pairs.
    :type intervals: list[tuple[int, int]]
    :return: DataFrame with run codes as index and one column per interval.
        Values are ``"mean ± std"`` formatted strings (``:.4f``).
    :rtype: pandas.DataFrame
    """
    columns = [f"{s}% to {e}%" for s, e in intervals]
    rows = {
        rc: {
            f"{s}% to {e}%": f"{stats[(s, e)][0]:.4f} \u00b1 {stats[(s, e)][1]:.4f}"
            for s, e in intervals
        }
        for rc, stats in stats_by_run.items()
    }
    return pd.DataFrame.from_dict(rows, orient="index")[columns]


# ---------------------------------------------------------------------------
# deg_time_table
# ---------------------------------------------------------------------------


def deg_time_table(stats_by_run: dict, markers: list) -> pd.DataFrame:
    """Build a summary DataFrame for degradation times.

    Column names use the plain ``"X%"`` form (without units).

    :param stats_by_run: Maps run code to ``{marker_pct: (mean_min, std_min)}``
        as returned by
        :func:`~lysis.analysis.degradation.compute_degradation_marker_stats`.
    :type stats_by_run: dict[str, dict[int, tuple[float, float]]]
    :param markers: Ordered list of integer percentage milestones.
    :type markers: list[int]
    :return: DataFrame with run codes as index and one column per milestone.
        Values are ``"mean ± std"`` formatted strings (``:.2f``).
    :rtype: pandas.DataFrame
    """
    columns = [f"{m}%" for m in markers]
    rows = {
        rc: {
            f"{m}%": f"{stats[m][0]:.2f} \u00b1 {stats[m][1]:.2f}"
            for m in markers
        }
        for rc, stats in stats_by_run.items()
    }
    return pd.DataFrame.from_dict(rows, orient="index")[columns]


# ---------------------------------------------------------------------------
# parameters_table
# ---------------------------------------------------------------------------


def _make_param_label(attr_name: str, display_units, natural_units: dict) -> str:
    """Build the display label for a parameter row.

    :param attr_name: Attribute name on the parameters object.
    :type attr_name: str
    :param display_units: Explicit display unit string, or ``None`` to use the
        parameter's natural units.
    :param natural_units: Mapping of attribute name → natural unit string.
    :type natural_units: dict[str, str]
    :return: Label string, e.g. ``"pore_size (microns)"`` or ``"cols"``.
    :rtype: str
    """
    units = display_units or natural_units.get(attr_name)
    return f"{attr_name} ({units})" if units else attr_name


def parameters_table(
    values_by_run: dict,
    param_specs: list,
    add_names: list,
    natural_units: dict,
    run_codes: list,
) -> pd.DataFrame:
    """Build a summary DataFrame for Run parameters.

    The DataFrame index is a :class:`pandas.MultiIndex` with two levels:

    * **section** — ``"Macroscale"``, ``"Microscale"``, or ``"Additional"``
    * **parameter** — display label (e.g. ``"pore_size (microns)"``)

    Display functions can use the section level to add section separators or
    heading rows when rendering.

    :param values_by_run: Maps run code to ``{attr_name: formatted_str}`` as
        returned by ``_load_run_params`` in
        :mod:`lysis.cli.parameters`.
    :type values_by_run: dict[str, dict[str, str]]
    :param param_specs: Effective parameter spec list; each entry is a
        ``(attr_name, source, display_units, fmt)`` tuple where *source* is
        ``"macro"``, ``"micro"``, or ``"computed"``.
    :type param_specs: list[tuple]
    :param add_names: Extra attribute names added via ``--add``.
    :type add_names: list[str]
    :param natural_units: Mapping of attribute name → natural unit string
        (from :meth:`~lysis.config.parameters.MacroParameters.units`).
    :type natural_units: dict[str, str]
    :param run_codes: Ordered list of run codes (determines column order).
    :type run_codes: list[str]
    :return: DataFrame with MultiIndex rows and run codes as columns.
    :rtype: pandas.DataFrame
    """
    index_tuples = []
    attr_names = []

    for attr_name, source, display_units, _fmt in param_specs:
        section = "Microscale" if source == "micro" else "Macroscale"
        label = _make_param_label(attr_name, display_units, natural_units)
        index_tuples.append((section, label))
        attr_names.append(attr_name)

    for attr_name in add_names:
        label = _make_param_label(attr_name, None, natural_units)
        index_tuples.append(("Additional", label))
        attr_names.append(attr_name)

    present = [rc for rc in run_codes if rc in values_by_run]
    data = {
        rc: [values_by_run[rc].get(a, "N/A") for a in attr_names]
        for rc in present
    }

    index = pd.MultiIndex.from_tuples(index_tuples, names=["section", "parameter"])
    return pd.DataFrame(data, index=index, columns=present)


# ---------------------------------------------------------------------------
# compare_stats_table
# ---------------------------------------------------------------------------


def _fmt_data_diff_cell(entry: dict) -> str:
    """Format one :func:`~lysis.analysis.compare.compare_data_tables` cell.

    :param entry: Result dict with ``"status"`` plus optional
        ``"max_pct_diff"``, ``"location"``, ``"detail"`` keys.
    :type entry: dict
    :return: Display string (e.g. ``"OK"``, ``"+1.50% @ (3, 7)"``).
    :rtype: str
    """
    status = entry.get("status")
    if status == "match":
        return "OK"
    if status == "diff":
        pct = entry.get("max_pct_diff")
        loc = entry.get("location")
        if pct is None or pct != pct:
            return f"NaN @ {loc}"
        return f"{pct:+.2f}% @ {loc}"
    if status == "shape_mismatch":
        return f"shapes {entry.get('detail', '')}"
    if status == "missing":
        return f"missing ({entry.get('detail', '')})"
    return str(entry)


def compare_stats_table(results_by_run: dict) -> pd.DataFrame:
    """Build a summary DataFrame for :func:`compare_runs` results across Runs.

    For each entry in the ``"ks"`` sub-dict, emits two columns:
    ``"{label} KS"`` (test statistic, ``:.4f``) and ``"{label} p-value"``
    (``:.3g``).  For each entry in the ``"pct_diff"`` sub-dict, emits one
    column: ``"{label} % diff"`` (signed, ``:+.2f`` with a ``%`` suffix).
    For each entry in the ``"data_diff"`` sub-dict, emits one column using
    the dataset label as the header and a formatted cell (``OK`` / percent
    + location / shape-mismatch / missing).  Column groups appear in the
    order KS → pct-diff → data-diff; orderings within each group are
    taken from the first Run's result dict.

    :param results_by_run: Maps run code to a compare-runs result dict
        (see :func:`~lysis.analysis.compare.compare_runs`) or, for
        backwards compatibility, to a plain ``{label: KstestResult}`` mapping.
    :type results_by_run: dict[str, dict]
    :return: DataFrame with run codes as index and pre-formatted string values.
        Returns an empty DataFrame if *results_by_run* is empty.
    :rtype: pandas.DataFrame
    """
    if not results_by_run:
        return pd.DataFrame()

    first = next(iter(results_by_run.values()))
    # Support both the new {"ks": ..., "pct_diff": ..., "data_diff": ...}
    # format and the plain {label: KstestResult} format.
    new_format = isinstance(first, dict) and (
        "ks" in first or "pct_diff" in first or "data_diff" in first
    )

    if new_format:
        ks_labels = list(first.get("ks", {}).keys())
        pct_labels = list(first.get("pct_diff", {}).keys())
        data_labels = list(first.get("data_diff", {}).keys())
    else:
        ks_labels = list(first.keys())
        pct_labels = []
        data_labels = []

    columns = []
    for label in ks_labels:
        columns.append(f"{label} KS")
        columns.append(f"{label} p-value")
    for label in pct_labels:
        columns.append(f"{label} % diff")
    for label in data_labels:
        columns.append(label)

    rows = {}
    for rc, entry in results_by_run.items():
        if new_format:
            ks_dict = entry.get("ks", {})
            pct_dict = entry.get("pct_diff", {})
            data_dict = entry.get("data_diff", {})
        else:
            ks_dict = entry
            pct_dict = {}
            data_dict = {}

        row = {}
        for label in ks_labels:
            result = ks_dict[label]
            row[f"{label} KS"] = f"{result.statistic:.4f}"
            row[f"{label} p-value"] = f"{result.pvalue:.3g}"
        for label in pct_labels:
            pct = pct_dict[label]
            row[f"{label} % diff"] = "N/A" if pct != pct else f"{pct:+.2f}%"
        for label in data_labels:
            row[label] = _fmt_data_diff_cell(data_dict[label])
        rows[rc] = row
    return pd.DataFrame.from_dict(rows, orient="index", columns=columns)


# ---------------------------------------------------------------------------
# compare_data_diff_summary_table
# ---------------------------------------------------------------------------


#: Column names for :func:`compare_data_diff_summary_table`.
DATA_DIFF_SUMMARY_COLUMNS = ["Result", "Max % Diff", "Worst Table"]


def _worst_mismatch(mismatches: dict):
    """Pick the mismatch entry with the largest absolute percent difference.

    Non-finite percent differences (NaN / Inf from symmetric-formula
    cancellation) rank highest.  If no entry has a numeric percent
    difference (e.g. all mismatches are ``shape_mismatch`` or
    ``missing``), returns the first entry unchanged.

    :param mismatches: Mapping ``{label: entry}`` of non-match entries.
    :type mismatches: dict
    :return: ``(label, entry)`` of the worst mismatch.
    :rtype: tuple[str, dict]
    """
    def _key(item):
        pct = item[1].get("max_pct_diff")
        if pct is None:
            return -1.0  # rank structural issues below numeric diffs
        if not math.isfinite(pct):
            return math.inf
        return abs(pct)

    return max(mismatches.items(), key=_key)


def compare_data_diff_summary_table(results_by_run: dict) -> pd.DataFrame:
    """Build a compact summary of :func:`compare_data_tables` results.

    For each Run, reports whether every data table matched exactly or —
    if not — how many tables differ and which one has the largest
    symmetric percent difference.  Designed for console display where
    the full per-table breakdown from :func:`compare_stats_table` is
    too wide.

    Columns (:data:`DATA_DIFF_SUMMARY_COLUMNS`):

    - ``"Result"`` — ``"Exact Match"`` when every table matched, or
      ``"<n> of <total> differ"`` when some did not.
    - ``"Max % Diff"`` — signed ``:+.2f`` percent of the worst numeric
      mismatch, or ``"—"`` when no numeric difference exists.
    - ``"Worst Table"`` — dataset label (plus index tuple for numeric
      diffs, or structural detail for shape/missing issues) pointing
      at the mismatch most worth investigating.  ``"—"`` on exact
      match.

    :param results_by_run: Maps run code to a ``compare_runs`` result
        dict with a ``"data_diff"`` sub-dict.  Entries without
        ``"data_diff"`` are treated as exact matches over zero tables.
    :type results_by_run: dict[str, dict]
    :return: DataFrame with run codes as index and three string columns.
        Returns an empty DataFrame if *results_by_run* is empty.
    :rtype: pandas.DataFrame
    """
    if not results_by_run:
        return pd.DataFrame()

    rows = {}
    for rc, entry in results_by_run.items():
        data_diff = (
            entry.get("data_diff", {}) if isinstance(entry, dict) else {}
        )
        total = len(data_diff)
        mismatches = {
            label: info
            for label, info in data_diff.items()
            if info.get("status") != "match"
        }

        if not mismatches:
            rows[rc] = {
                "Result": "Exact Match",
                "Max % Diff": "—",
                "Worst Table": "—",
            }
            continue

        worst_label, worst_entry = _worst_mismatch(mismatches)
        status = worst_entry.get("status")
        if status == "diff":
            pct = worst_entry.get("max_pct_diff")
            loc = worst_entry.get("location")
            pct_str = (
                "NaN" if pct is None or pct != pct else f"{pct:+.2f}%"
            )
            worst_str = f"{worst_label} @ {loc}" if loc is not None else worst_label
        elif status == "shape_mismatch":
            pct_str = "—"
            worst_str = f"{worst_label} (shapes {worst_entry.get('detail', '')})"
        elif status == "missing":
            pct_str = "—"
            worst_str = f"{worst_label} (missing: {worst_entry.get('detail', '')})"
        else:
            pct_str = "—"
            worst_str = worst_label

        rows[rc] = {
            "Result": f"{len(mismatches)} of {total} differ",
            "Max % Diff": pct_str,
            "Worst Table": worst_str,
        }

    return pd.DataFrame.from_dict(
        rows, orient="index", columns=DATA_DIFF_SUMMARY_COLUMNS
    )
