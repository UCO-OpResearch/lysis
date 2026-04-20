"""Two-sample comparison of simulation Runs.

Provides functions for comparing two Runs using:

- 2-sample Kolmogorov-Smirnov tests (:func:`scipy.stats.ks_2samp`) on
  per-simulation arrays from :data:`MEASURE_EXTRACTORS`;
- symmetric percent differences on scalar summary statistics from
  :data:`STATS_COMPUTERS`.

Both dispatch tables are keyed by a short measure-set name (e.g.
``"micro-stats"``) and are designed to be extended by adding a new key to
either or both.

Typical workflow::

    from lysis.analysis.compare import compare_runs

    result = compare_runs(run1, run2, "micro-stats")
    for label, r in result["ks"].items():
        print(label, r.statistic, r.pvalue)
    for label, pct in result["pct_diff"].items():
        print(label, f"{pct:+.2f}%")
"""

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2026, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"

from typing import TYPE_CHECKING, Callable, Dict

import numpy as np
from scipy.stats import ks_2samp

if TYPE_CHECKING:
    from lysis.config.run import Run

__all__ = [
    "MEASURE_EXTRACTORS",
    "STATS_COMPUTERS",
    "compare_runs",
    "compute_stats",
    "extract_measures",
    "percent_difference",
]


# ---------------------------------------------------------------------------
# Measure extractors (arrays, for KS tests)
# ---------------------------------------------------------------------------


def _extract_micro_stats_measures(run: "Run") -> Dict[str, np.ndarray]:
    """Return per-simulation arrays for the ``micro-stats`` measure set.

    Mirrors the arrays used by
    :func:`lysis.analysis.microscale.compute_micro_statistics`: lysis times
    are restricted to degraded fibers and converted to minutes; tPA leaving
    times cover all simulations in seconds.

    :param run: Run with open microscale data.
    :type run: Run
    :return: Mapping from human-readable measure label to 1-D array.
    :rtype: dict[str, numpy.ndarray]
    """
    micro = run.data.microscale_out
    fiber_degraded = micro.fiber_degraded[:]
    sim_final_time = micro.sim_final_time[:]
    tpa_leaving_time = micro.tpa_leaving_time[:]

    lysis_times = sim_final_time[fiber_degraded] / 60.0
    return {
        "Lysis Time (min)": lysis_times,
        "tPA Leaving Time (sec)": tpa_leaving_time,
    }


#: Dispatch table mapping a measure-set name to an extractor function.
#:
#: An extractor takes an opened :class:`~lysis.config.run.Run` and returns
#: a dict ``{measure_label: 1-D array}`` used for KS testing.
MEASURE_EXTRACTORS: Dict[str, Callable[["Run"], Dict[str, np.ndarray]]] = {
    "micro-stats": _extract_micro_stats_measures,
}


# ---------------------------------------------------------------------------
# Scalar stat computers (scalars, for percent-difference comparison)
# ---------------------------------------------------------------------------


def _compute_micro_stats_scalars(run: "Run") -> Dict[str, float]:
    """Return scalar micro-stats for percent-difference comparison.

    Emits the five scalar quantities shown in the ``lysis micro-stats``
    table (std values are omitted since percent difference of a standard
    deviation is not usually meaningful here).

    :param run: Run with open microscale data.
    :type run: Run
    :return: Mapping from display label to scalar value.
    :rtype: dict[str, float]
    """
    from lysis.analysis.microscale import compute_micro_statistics

    stats = compute_micro_statistics(run)
    return {
        "Fibers Degraded": stats["fibers_degraded"],
        "Mean Lysis Time (min)": stats["lysis_time_mean"],
        "Median Lysis Time (min)": stats["lysis_time_median"],
        "Mean tPA Leaving Time (sec)": stats["tpa_leaving_mean"],
        "Median tPA Leaving Time (sec)": stats["tpa_leaving_median"],
    }


#: Dispatch table mapping a measure-set name to a scalar-stats computer.
#:
#: A computer takes an opened :class:`~lysis.config.run.Run` and returns a
#: dict ``{stat_label: float}`` used for percent-difference comparison.
STATS_COMPUTERS: Dict[str, Callable[["Run"], Dict[str, float]]] = {
    "micro-stats": _compute_micro_stats_scalars,
}


# ---------------------------------------------------------------------------
# Percent difference
# ---------------------------------------------------------------------------


def percent_difference(val1: float, val2: float) -> float:
    """Symmetric percent difference between two scalar values.

    Uses the symmetric (relative percent difference) formula
    ``200 * (val2 - val1) / (val1 + val2)``, with sign preserved so positive
    means *val2 > val1*.  Returns ``0.0`` when both inputs are exactly 0,
    and ``NaN`` when the inputs cancel to a zero denominator or either input
    is ``NaN``.

    :param val1: First value (e.g. from ``run1``).
    :type val1: float
    :param val2: Second value (e.g. from ``run2``).
    :type val2: float
    :return: Symmetric percent difference (``-200`` to ``+200`` for
        same-sign inputs), ``NaN`` on zero denominator, ``0.0`` when both
        inputs are 0.
    :rtype: float
    """
    if np.isnan(val1) or np.isnan(val2):
        return float("nan")
    if val1 == 0 and val2 == 0:
        return 0.0
    denom = (val1 + val2) / 2.0
    if denom == 0:
        return float("nan")
    return 100.0 * (val2 - val1) / denom


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def extract_measures(run: "Run", which: str) -> Dict[str, np.ndarray]:
    """Extract the per-simulation measure arrays named by ``which``.

    :param run: Run with the required data collections open.
    :type run: Run
    :param which: Key into :data:`MEASURE_EXTRACTORS`.
    :type which: str
    :return: Mapping from measure label to 1-D array.
    :rtype: dict[str, numpy.ndarray]
    :raises KeyError: If ``which`` is not a known measure set.
    """
    if which not in MEASURE_EXTRACTORS:
        raise KeyError(
            f"Unknown measure set {which!r}. "
            f"Available: {sorted(MEASURE_EXTRACTORS.keys())}"
        )
    return MEASURE_EXTRACTORS[which](run)


def compute_stats(run: "Run", which: str) -> Dict[str, float]:
    """Compute the scalar summary stats named by ``which``.

    :param run: Run with the required data collections open.
    :type run: Run
    :param which: Key into :data:`STATS_COMPUTERS`.
    :type which: str
    :return: Mapping from stat label to scalar value.  Returns ``{}`` if
        ``which`` has no registered stats computer.
    :rtype: dict[str, float]
    """
    computer = STATS_COMPUTERS.get(which)
    if computer is None:
        return {}
    return computer(run)


def compare_runs(run1: "Run", run2: "Run", which: str) -> dict:
    """Compare two Runs with KS tests and scalar percent-differences.

    For the measure set selected by ``which``:

    - apply :func:`scipy.stats.ks_2samp` to each per-simulation array pair
      from :data:`MEASURE_EXTRACTORS`;
    - compute :func:`percent_difference` between each scalar stat pair
      from :data:`STATS_COMPUTERS`.

    :param run1: First Run (with data open).
    :type run1: Run
    :param run2: Second Run (with data open).
    :type run2: Run
    :param which: Key into :data:`MEASURE_EXTRACTORS` / :data:`STATS_COMPUTERS`.
    :type which: str
    :return: Dict with two sub-dicts:

        - ``"ks"`` — ``{measure_label: KstestResult}`` with ``statistic`` and
          ``pvalue`` attributes;
        - ``"pct_diff"`` — ``{stat_label: float}`` in percent, signed so
          positive means *run2 > run1*.

        Either sub-dict may be empty if no entries are registered for
        ``which``.
    :rtype: dict
    :raises KeyError: If ``which`` is not a known measure set, or if the
        two extractors disagree on labels.
    """
    measures_1 = extract_measures(run1, which)
    measures_2 = extract_measures(run2, which)
    ks_results = {}
    for label, arr1 in measures_1.items():
        if label not in measures_2:
            raise KeyError(
                f"Measure {label!r} missing from second Run's extractor output."
            )
        ks_results[label] = ks_2samp(arr1, measures_2[label])

    stats_1 = compute_stats(run1, which)
    stats_2 = compute_stats(run2, which)
    pct_diffs = {}
    for label, v1 in stats_1.items():
        if label not in stats_2:
            raise KeyError(
                f"Stat {label!r} missing from second Run's stats output."
            )
        pct_diffs[label] = percent_difference(float(v1), float(stats_2[label]))

    return {"ks": ks_results, "pct_diff": pct_diffs}
