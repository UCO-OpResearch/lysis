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
"""

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
