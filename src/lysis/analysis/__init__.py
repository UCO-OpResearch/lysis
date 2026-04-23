"""Analysis functions for fibrinolysis simulation Runs and Experiments.

This package provides tools for quantifying and visualising simulation output.
Currently available submodules:

- :mod:`~lysis.analysis.degradation`: Clot-degradation analysis and plots.
- :mod:`~lysis.analysis.microscale`: Microscale simulation summary statistics.
- :mod:`~lysis.analysis.summary`: Summary table functions returning
  :class:`pandas.DataFrame` objects for CLI display and notebook use.

Example::

    from lysis.analysis.degradation import (
        find_degraded_fraction,
        find_degradation_marker_frames,
        find_degradation_marker_times,
        mean_degradation_rate,
    )
    from lysis.analysis.microscale import compute_micro_statistics
    from lysis.analysis.summary import micro_stats_table, macro_stats_table
"""

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2026, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"

from .degradation import (
    calculate_time_row_exposed,
    compute_degradation_marker_stats,
    compute_degradation_rate_stats,
    compute_run_statistics,
    degradation_rates,
    fiber_degradation_linear_extrapolation,
    find_degradation_fronts,
    find_degradation_marker_frames,
    find_degradation_marker_times,
    find_degraded_fraction,
    find_front,
    find_row_deg_fraction,
    get_processing_time,
    get_total_binds,
    get_unbind_amounts,
    mean_degradation_rate,
    mean_front_velocity,
    plot_degradation_percent,
    plot_front_degradation,
)
from .compare import (
    DATA_TABLE_EXTRACTORS,
    MEASURE_EXTRACTORS,
    STATS_COMPUTERS,
    available_measure_sets,
    compare_data_tables,
    compare_runs,
    compute_stats,
    extract_measures,
    percent_difference,
)
from .fiber_replay import FiberReplayCursor
from .microscale import compute_micro_statistics
from .summary import (
    compare_data_diff_summary_table,
    compare_stats_table,
    deg_rate_table,
    deg_time_table,
    macro_stats_table,
    micro_stats_table,
    parameters_table,
)

__all__ = [
    "DATA_TABLE_EXTRACTORS",
    "FiberReplayCursor",
    "MEASURE_EXTRACTORS",
    "STATS_COMPUTERS",
    "available_measure_sets",
    "calculate_time_row_exposed",
    "compare_data_diff_summary_table",
    "compare_data_tables",
    "compare_runs",
    "compare_stats_table",
    "compute_stats",
    "deg_rate_table",
    "deg_time_table",
    "extract_measures",
    "percent_difference",
    "macro_stats_table",
    "micro_stats_table",
    "parameters_table",
    "compute_degradation_marker_stats",
    "compute_degradation_rate_stats",
    "compute_micro_statistics",
    "compute_run_statistics",
    "degradation_rates",
    "fiber_degradation_linear_extrapolation",
    "find_degradation_fronts",
    "find_degradation_marker_frames",
    "find_degradation_marker_times",
    "find_degraded_fraction",
    "find_front",
    "find_row_deg_fraction",
    "get_processing_time",
    "get_total_binds",
    "get_unbind_amounts",
    "mean_degradation_rate",
    "mean_front_velocity",
    "plot_degradation_percent",
    "plot_front_degradation",
]
