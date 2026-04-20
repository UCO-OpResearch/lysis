"""Two-sample comparison of simulation Runs.

Provides functions for comparing per-simulation measure arrays between two
Runs using the 2-sample Kolmogorov-Smirnov test
(:func:`scipy.stats.ks_2samp`).

Measure sets are keyed by a short label (e.g. ``"micro-stats"``) and map to
extractor functions that take an opened Run and return a dict of
``{measure_label: 1-D array}``.  The dispatch table
:data:`MEASURE_EXTRACTORS` is designed to be extended: add a new key mapping
to an extractor function to make a new measure set available.

Typical workflow::

    from lysis.analysis.compare import compare_runs_ks

    results = compare_runs_ks(run1, run2, "micro-stats")
    for label, res in results.items():
        print(label, res.statistic, res.pvalue)
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
    "extract_measures",
    "compare_runs_ks",
]


# ---------------------------------------------------------------------------
# Measure extractors
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
#: a dict ``{measure_label: 1-D array}``.  Add new entries to extend the set
#: of comparisons available to :func:`compare_runs_ks` and the
#: ``lysis compare`` CLI command.
MEASURE_EXTRACTORS: Dict[str, Callable[["Run"], Dict[str, np.ndarray]]] = {
    "micro-stats": _extract_micro_stats_measures,
}


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


def compare_runs_ks(run1: "Run", run2: "Run", which: str) -> dict:
    """Compare two Runs using the 2-sample Kolmogorov-Smirnov test.

    For each measure in the set selected by ``which``, extract the
    corresponding 1-D array from each Run and apply
    :func:`scipy.stats.ks_2samp`.

    :param run1: First Run (with data open).
    :type run1: Run
    :param run2: Second Run (with data open).
    :type run2: Run
    :param which: Key into :data:`MEASURE_EXTRACTORS`.
    :type which: str
    :return: Ordered mapping ``{measure_label: KstestResult}``, where each
        value has ``statistic`` and ``pvalue`` attributes (see
        :func:`scipy.stats.ks_2samp`).
    :rtype: dict
    :raises KeyError: If ``which`` is not a known measure set, or if the
        two extractors disagree on measure labels.
    """
    measures_1 = extract_measures(run1, which)
    measures_2 = extract_measures(run2, which)

    results = {}
    for label, arr1 in measures_1.items():
        if label not in measures_2:
            raise KeyError(
                f"Measure {label!r} missing from second Run's extractor output."
            )
        results[label] = ks_2samp(arr1, measures_2[label])
    return results
