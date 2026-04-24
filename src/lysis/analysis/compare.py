"""Two-sample comparison of simulation Runs.

Provides functions for comparing two Runs using:

- 2-sample Kolmogorov-Smirnov tests (:func:`scipy.stats.ks_2samp`) on
  per-simulation arrays from :data:`MEASURE_EXTRACTORS`;
- symmetric percent differences on scalar summary statistics from
  :data:`STATS_COMPUTERS`.

Both dispatch tables are keyed by a short measure-set name (e.g.
``"micro-stats"``, ``"macro-stats"``) and are designed to be extended by
adding a new key to either or both.

Element-wise data-table comparisons (:func:`compare_data_tables`) apply a
:data:`_MAX_ULPS`-ULP tolerance on floating-point dtypes: two float
values are treated as equivalent when they are within that many units
in the last place of their IEEE-754 representation.  This suppresses
not just the single-bit rounding drift between different Fortran
compilations but also the small multi-ULP accumulation that shows up
downstream (``count*tstep`` compounding across time steps).  Integer,
boolean, and string dtypes still use exact equality.

Typical workflow::

    from lysis.analysis.compare import compare_runs
    from lysis.tools.display import _format_pct_base

    result = compare_runs(run1, run2, "macro-stats")
    for label, r in result["ks"].items():
        print(label, r.statistic, r.pvalue)
    for label, pct in result["pct_diff"].items():
        print(label, _format_pct_base(pct))
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
    "DATA_TABLE_EXTRACTORS",
    "MEASURE_EXTRACTORS",
    "STATS_COMPUTERS",
    "available_measure_sets",
    "compare_data_tables",
    "compare_runs",
    "compute_stats",
    "extract_measures",
    "percent_difference",
]

#: Dataset names (log tables) that must never be compared.
_LOG_DATASETS = frozenset({"micro_log", "macro_log"})

#: Floating-point equality tolerance, in units in the last place.
#:
#: Two float values are treated as equivalent by :func:`_values_match`
#: (and therefore by :func:`compare_data_tables` and every display path
#: that routes through it) when they differ by at most this many ULPs in
#: their IEEE-754 representation.  Raise this number to suppress more
#: compiler-rounding noise at the cost of masking small real drifts;
#: lower it to surface more bit-level disagreement.
_MAX_ULPS = 2


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


def _extract_macro_stats_measures(run: "Run") -> Dict[str, np.ndarray]:
    """Return per-simulation arrays for the ``macro-stats`` measure set.

    Uses the per-simulation outputs of
    :func:`lysis.analysis.degradation.mean_degradation_rate`,
    :func:`lysis.analysis.degradation.find_degradation_marker_times`, and
    :func:`lysis.analysis.degradation.per_sim_front_velocity`.  The
    degradation rate is converted from fraction/min to %/min so the K-S
    column units match the paired ``% diff`` column.

    :param run: Run with open macroscale data.
    :type run: Run
    :return: Mapping from human-readable measure label to 1-D array of
        length ``n_sims``.
    :rtype: dict[str, numpy.ndarray]
    """
    from lysis.analysis.degradation import (
        find_degradation_marker_times,
        mean_degradation_rate,
        per_sim_front_velocity,
    )

    deg_rate, _offset, _deg_start = mean_degradation_rate(run)
    marker_times = find_degradation_marker_times(run, [0.0, 0.25, 0.5, 0.75, 1.0])
    front_mean, _front_std = per_sim_front_velocity(run)
    return {
        "Degradation rate (%/min)": deg_rate * 100,
        "Time to full clot degradation (min)": marker_times[:, -1],
        "Front Velocity (microns/min)": front_mean,
    }


#: Dispatch table mapping a measure-set name to an extractor function.
#:
#: An extractor takes an opened :class:`~lysis.config.run.Run` and returns
#: a dict ``{measure_label: 1-D array}`` used for KS testing.
MEASURE_EXTRACTORS: Dict[str, Callable[["Run"], Dict[str, np.ndarray]]] = {
    "micro-stats": _extract_micro_stats_measures,
    "macro-stats": _extract_macro_stats_measures,
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


def _compute_macro_stats_scalars(run: "Run") -> Dict[str, float]:
    """Return scalar macro-stats for percent-difference comparison.

    Emits the four scalar summaries most commonly tracked for macroscale
    Runs: mean degradation rate, mean time to full clot degradation, mean
    first-passage time, and mean front velocity.  Delegates to
    :func:`lysis.analysis.degradation.compute_run_statistics` and pulls
    the ``"Mean"`` row of each metric.

    :param run: Run with open macroscale data.
    :type run: Run
    :return: Mapping from display label to scalar value.
    :rtype: dict[str, float]
    """
    from lysis.analysis.degradation import compute_run_statistics

    stats = compute_run_statistics(run)
    return {
        "Mean Degradation rate (%/min)": float(
            stats[("Degradation rate (%/min)", "Mean")]
        ),
        "Mean Time to full clot degradation (min)": float(
            stats[("Time to full clot degradation (min)", "Mean")]
        ),
        "Mean First passage time (min)": float(
            stats[("First passage time (min)", "Mean")]
        ),
        "Mean Front Velocity (microns/min)": float(
            stats[("Front Velocity (microns/min)", "Mean")]
        ),
    }


#: Dispatch table mapping a measure-set name to a scalar-stats computer.
#:
#: A computer takes an opened :class:`~lysis.config.run.Run` and returns a
#: dict ``{stat_label: float}`` used for percent-difference comparison.
STATS_COMPUTERS: Dict[str, Callable[["Run"], Dict[str, float]]] = {
    "micro-stats": _compute_micro_stats_scalars,
    "macro-stats": _compute_macro_stats_scalars,
}


# ---------------------------------------------------------------------------
# Data-table extractors (for exact-match comparison)
# ---------------------------------------------------------------------------


def _extract_data_arrays(run: "Run", include_macro: bool) -> Dict[str, np.ndarray]:
    """Extract every on-disk (non-log) data table from a Run's DataStore.

    Dataset contents are read into numpy arrays keyed by a display label of
    the form ``microscale_out/<name>`` or ``macroscale_out[<sim>]/<name>``.
    Log tables (:data:`_LOG_DATASETS`) are skipped.  Derived datasets with
    ``data_location=None`` are already excluded by the :attr:`datasets`
    property on :class:`DataCollection` / :class:`SimulationView`.

    :param run: Run with data open.
    :type run: Run
    :param include_macro: If ``True``, also include per-simulation
        ``macroscale_out`` datasets.
    :type include_macro: bool
    :return: Mapping from label to numpy array copy.
    :rtype: dict[str, numpy.ndarray]
    """
    tables: Dict[str, np.ndarray] = {}
    data = run.data

    collections = data.collections
    if "microscale_out" in collections:
        micro = data.microscale_out
        for name in micro.datasets:
            if name in _LOG_DATASETS:
                continue
            tables[f"microscale_out/{name}"] = getattr(micro, name)[:]

    if include_macro and "macroscale_out" in collections:
        macro = data.macroscale_out
        n_sims = run.macro_params.macro_simulations
        for sim in range(n_sims):
            view = macro[sim]
            for name in view.datasets:
                if name in _LOG_DATASETS:
                    continue
                tables[f"macroscale_out[{sim:02}]/{name}"] = getattr(view, name)[:]

    return tables


#: Dispatch table mapping a data-comparison name to a table extractor.
#:
#: An extractor takes an opened :class:`~lysis.config.run.Run` and returns
#: a dict ``{label: ndarray}`` of on-disk data tables (log tables and
#: HDF5 attributes are always excluded).
DATA_TABLE_EXTRACTORS: Dict[str, Callable[["Run"], Dict[str, np.ndarray]]] = {
    "micro-data": lambda run: _extract_data_arrays(run, include_macro=False),
    "data": lambda run: _extract_data_arrays(run, include_macro=True),
}


def available_measure_sets() -> list:
    """Return the full list of measure-set names accepted by :func:`compare_runs`.

    Combines :data:`MEASURE_EXTRACTORS` (KS + scalar pct-diff comparisons)
    and :data:`DATA_TABLE_EXTRACTORS` (exact data-table comparisons).

    :return: Sorted list of measure-set keys.
    :rtype: list[str]
    """
    return sorted(set(MEASURE_EXTRACTORS) | set(DATA_TABLE_EXTRACTORS))


# ---------------------------------------------------------------------------
# Data-table comparison helpers
# ---------------------------------------------------------------------------


def _values_match(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Element-wise equivalence with an N-ULP tolerance on floating-point arrays.

    For float dtypes a position is considered equivalent when *b* lies in
    the closed interval ``[prev_n(a), next_n(a)]``, where ``prev_n`` /
    ``next_n`` are ``np.nextafter`` applied :data:`_MAX_ULPS` times
    toward ``-inf`` / ``+inf`` respectively.  In other words, *a* and
    *b* are either exactly equal or differ by at most :data:`_MAX_ULPS`
    units in the last place.  This absorbs the last-bit rounding
    variation produced by different Fortran compilations (operation
    reordering, FMA) plus the small multi-ULP accumulation that shows
    up downstream of a drifting time step.

    NaN positions compare as False (matching NumPy's ``==`` semantics).
    Signed zeros and infinities compare equal to themselves.  Integer,
    boolean, string, and other non-floating dtypes fall back to exact
    ``==`` equality regardless of :data:`_MAX_ULPS`.

    :param a: First array.
    :type a: numpy.ndarray
    :param b: Second array (broadcastable to *a*).
    :type b: numpy.ndarray
    :return: Boolean ndarray, same broadcast shape.
    :rtype: numpy.ndarray
    """
    if np.issubdtype(a.dtype, np.floating) and np.issubdtype(b.dtype, np.floating):
        af = np.asarray(a, dtype=np.float64)
        bf = np.asarray(b, dtype=np.float64)
        lo = af
        hi = af
        for _ in range(_MAX_ULPS):
            lo = np.nextafter(lo, -np.inf)
            hi = np.nextafter(hi, np.inf)
        return (lo <= bf) & (bf <= hi)
    return a == b


def _arraywise_max_pct(a: np.ndarray, b: np.ndarray):
    """Find the element with the largest absolute percent difference.

    Uses the same symmetric formula as :func:`percent_difference`
    (``200 * (b - a) / (a + b)``) element-wise.  Positions where both
    inputs are zero contribute 0.  Positions that produce NaN or Inf
    (e.g. opposite-sign cancellation) are ranked highest so true
    mismatches are surfaced even when the symmetric denominator breaks
    down.  Positions that :func:`_values_match` treats as equivalent
    (within :data:`_MAX_ULPS` ULPs for floats, exact otherwise) are
    excluded from the argmax so they cannot be reported as the worst
    element.

    :param a: First array.
    :type a: numpy.ndarray
    :param b: Second array (same shape as *a*).
    :type b: numpy.ndarray
    :return: ``(signed_pct, location)`` where *location* is a tuple of
        ``int`` indices into the original array shape.
    :rtype: tuple[float, tuple[int, ...]]
    """
    a_float = np.asarray(a, dtype=np.float64)
    b_float = np.asarray(b, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = (a_float + b_float) / 2.0
        pct = 100.0 * (b_float - a_float) / denom
    both_zero = (a_float == 0) & (b_float == 0)
    pct = np.where(both_zero, 0.0, pct)
    abs_pct = np.abs(pct)
    # Rank non-finite (NaN from 0/0 cancellation, Inf from a+b=0) highest
    abs_pct_ranked = np.where(np.isfinite(abs_pct), abs_pct, np.inf)
    # Exclude positions that match under the 1-ULP-tolerant rule.
    match_mask = np.asarray(_values_match(a, b))
    abs_pct_ranked = np.where(match_mask, -1.0, abs_pct_ranked)
    idx_flat = int(np.argmax(abs_pct_ranked))
    location = np.unravel_index(idx_flat, pct.shape)
    return float(pct.flat[idx_flat]), tuple(int(i) for i in location)


def _count_element_mismatches(a: np.ndarray, b: np.ndarray) -> int:
    """Count positions where *a* and *b* differ element-wise.

    Floating-point positions that differ by at most :data:`_MAX_ULPS`
    ULPs are treated as equivalent and not counted (see
    :func:`_values_match`).  For structured arrays, an element counts as
    a mismatch if *any* of its fields differs (so the count is the
    number of differing records, not the sum across fields).

    :param a: First array.
    :type a: numpy.ndarray
    :param b: Second array (same shape as *a*).
    :type b: numpy.ndarray
    :rtype: int
    """
    if a.dtype.names is not None:
        any_diff = np.zeros(a.shape, dtype=bool)
        for field in a.dtype.names:
            any_diff |= ~_values_match(a[field], b[field])
        return int(np.sum(any_diff))
    return int(np.sum(~_values_match(a, b)))


def _compare_arrays(a: np.ndarray, b: np.ndarray) -> dict:
    """Compare two arrays with N-ULP float tolerance, reporting max percent diff.

    Floating-point positions that differ by at most :data:`_MAX_ULPS`
    units in the last place (see :func:`_values_match`) are treated as
    equivalent; all other dtypes use exact equality.  This distinguishes
    true numerical divergence from the last-bit rounding variation that
    different Fortran compilations inject into quantities derived from
    ``q*delx**2/(12*Diff)`` (e.g. ``snapshot_time``), plus the small
    multi-ULP drift such values accumulate downstream.

    Returns a dict with these keys:

    - ``"status"`` — one of ``"match"``, ``"diff"``, ``"shape_mismatch"``.
    - ``"max_pct_diff"`` — signed percent difference at the worst
      element (``0.0`` on match, ``None`` on shape mismatch).
    - ``"location"`` — tuple of indices identifying the worst element
      (``None`` on match or shape mismatch).  For structured arrays the
      first entry is the field name.
    - ``"mismatches"`` — count of elements that differ under the
      :data:`_MAX_ULPS`-ULP-tolerant rule.  ``0`` on match, ``None`` on
      shape mismatch.  For structured arrays this counts records where
      any field differs, not the sum across fields.
    - ``"total"`` — total element count of the (shape-matched) arrays.
      ``None`` on shape mismatch.
    - ``"detail"`` — present only on ``"shape_mismatch"`` with a
      human-readable ``shape1 vs shape2`` string.

    :param a: First array (from ``run1``).
    :type a: numpy.ndarray
    :param b: Second array (from ``run2``).
    :type b: numpy.ndarray
    :rtype: dict
    """
    if a.shape != b.shape:
        return {
            "status": "shape_mismatch",
            "max_pct_diff": None,
            "location": None,
            "mismatches": None,
            "total": None,
            "detail": f"{a.shape} vs {b.shape}",
        }
    total = int(a.size)
    mismatches = _count_element_mismatches(a, b)
    if mismatches == 0:
        return {
            "status": "match",
            "max_pct_diff": 0.0,
            "location": None,
            "mismatches": 0,
            "total": total,
        }

    if a.dtype.names is not None:
        # Structured array: scan each field, keep the worst.
        best_pct = 0.0
        best_loc = None
        for field in a.dtype.names:
            field_a = a[field]
            field_b = b[field]
            if bool(np.all(_values_match(field_a, field_b))):
                continue
            pct, loc = _arraywise_max_pct(field_a, field_b)
            worse = best_loc is None or (
                not np.isfinite(best_pct) or abs(pct) > abs(best_pct)
            )
            if worse:
                best_pct = pct
                best_loc = (field,) + loc
        if best_loc is None:
            # Defensive: mismatches > 0 but no field differs under the
            # tolerant rule.  Treat as a match.
            return {
                "status": "match",
                "max_pct_diff": 0.0,
                "location": None,
                "mismatches": 0,
                "total": total,
            }
        return {
            "status": "diff",
            "max_pct_diff": best_pct,
            "location": best_loc,
            "mismatches": mismatches,
            "total": total,
        }

    pct, loc = _arraywise_max_pct(a, b)
    return {
        "status": "diff",
        "max_pct_diff": pct,
        "location": loc,
        "mismatches": mismatches,
        "total": total,
    }


def compare_data_tables(run1: "Run", run2: "Run", which: str) -> dict:
    """Compare on-disk data tables between two Runs element-wise.

    For the data-table set selected by ``which`` (see
    :data:`DATA_TABLE_EXTRACTORS`), reads every non-log dataset from
    each Run's HDF5 file and compares the two arrays element-wise.
    Floating-point positions that differ by at most :data:`_MAX_ULPS`
    ULPs are treated as equivalent (see :func:`_values_match`); every
    other dtype uses exact equality.  When any position still differs,
    reports the maximum symmetric percent difference
    (:func:`percent_difference`) across the dataset and the location of
    that worst element.

    HDF5 attributes and log tables (``micro_log``, ``macro_log``) are
    never compared.  Datasets that exist in one Run but not the other
    are reported with ``status="missing"``.

    :param run1: First Run with data open.
    :type run1: Run
    :param run2: Second Run with data open.
    :type run2: Run
    :param which: Key into :data:`DATA_TABLE_EXTRACTORS`.
    :type which: str
    :return: Dict with a single ``"data_diff"`` sub-dict mapping each
        dataset label to a result dict.  The result dict always has
        ``"status"``; ``"match"`` results additionally carry
        ``max_pct_diff=0.0`` and ``location=None``; ``"diff"`` results
        carry a signed ``"max_pct_diff"`` and index tuple
        ``"location"``; ``"shape_mismatch"`` results carry a
        ``"detail"`` string; ``"missing"`` results carry a
        ``"detail"`` describing which Run is missing the dataset.
    :rtype: dict
    :raises KeyError: If ``which`` is not a known data-table set.
    """
    if which not in DATA_TABLE_EXTRACTORS:
        raise KeyError(
            f"Unknown data-table set {which!r}. "
            f"Available: {sorted(DATA_TABLE_EXTRACTORS.keys())}"
        )
    tables1 = DATA_TABLE_EXTRACTORS[which](run1)
    tables2 = DATA_TABLE_EXTRACTORS[which](run2)

    all_labels = sorted(set(tables1) | set(tables2))
    results: Dict[str, dict] = {}
    for label in all_labels:
        if label not in tables1:
            results[label] = {
                "status": "missing",
                "max_pct_diff": None,
                "location": None,
                "detail": "not in run1",
            }
        elif label not in tables2:
            results[label] = {
                "status": "missing",
                "max_pct_diff": None,
                "location": None,
                "detail": "not in run2",
            }
        else:
            results[label] = _compare_arrays(tables1[label], tables2[label])
    return {"data_diff": results}


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

        When ``which`` selects a data-table comparison (see
        :data:`DATA_TABLE_EXTRACTORS`) the call is routed to
        :func:`compare_data_tables`, which returns a single ``"data_diff"``
        sub-dict instead.
    :rtype: dict
    :raises KeyError: If ``which`` is not a known measure set, or if the
        two extractors disagree on labels.
    """
    if which in DATA_TABLE_EXTRACTORS:
        return compare_data_tables(run1, run2, which)

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
