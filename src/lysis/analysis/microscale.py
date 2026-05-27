"""Analysis functions for microscale fibrinolysis simulation output.

Provides functions for computing summary statistics from microscale simulation
output datasets.  Mirrors the "## Microscale Measures" notebook section.

Typical workflow::

    from lysis.analysis.microscale import compute_micro_statistics

    stats = compute_micro_statistics(run)
    print(stats["fibers_degraded"])
    print(stats["lysis_time_mean"], "±", stats["lysis_time_std"], "min")
"""

__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from lysis.config.run import Run

__all__ = ["compute_micro_statistics"]


def compute_micro_statistics(run: "Run") -> dict:
    """Compute summary statistics from microscale simulation output.

    Reads ``fiber_degraded``, ``sim_final_time``, and ``tpa_leaving_time``
    from ``microscale_out`` and returns summary statistics matching the
    notebook's ``microscale_df`` table.

    :param run: Run object with open microscale data.
    :type run: Run
    :return: Dict with keys:

        - ``fibers_degraded`` — total number of degraded fibers (int)
        - ``lysis_time_mean`` — mean lysis time across degraded fibers (min)
        - ``lysis_time_std`` — standard deviation of lysis times (min)
        - ``lysis_time_median`` — median lysis time (min)
        - ``tpa_leaving_mean`` — mean tPA leaving time across all simulations (sec)
        - ``tpa_leaving_std`` — standard deviation of tPA leaving times (sec)
        - ``tpa_leaving_median`` — median tPA leaving time (sec)

    :rtype: dict
    """
    micro = run.data.microscale_out
    fiber_degraded = micro.fiber_degraded[:]
    sim_final_time = micro.sim_final_time[:]
    tpa_leaving_time = micro.tpa_leaving_time[:]

    lysis_times = sim_final_time[fiber_degraded]

    if lysis_times.size == 0:
        lysis_mean = lysis_std = lysis_median = float("nan")
    else:
        lysis_mean = float(np.mean(lysis_times) / 60)
        lysis_std = float(np.std(lysis_times) / 60)
        lysis_median = float(np.median(lysis_times) / 60)

    return {
        "fibers_degraded": int(np.sum(fiber_degraded)),
        "lysis_time_mean": lysis_mean,
        "lysis_time_std": lysis_std,
        "lysis_time_median": lysis_median,
        "tpa_leaving_mean": float(np.mean(tpa_leaving_time)),
        "tpa_leaving_std": float(np.std(tpa_leaving_time)),
        "tpa_leaving_median": float(np.median(tpa_leaving_time)),
    }
