"""Analysis functions for fibrinolysis simulation Runs and Experiments.

This package provides tools for quantifying and visualising simulation output.
Currently available submodules:

- :mod:`~lysis.analysis.degradation`: Clot-degradation analysis and plots.

Example::

    from lysis.analysis.degradation import (
        find_degraded_fraction,
        find_degradation_marker_frames,
        find_degradation_marker_times,
        mean_degradation_rate,
    )
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

__all__ = [
    "calculate_time_row_exposed",
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
