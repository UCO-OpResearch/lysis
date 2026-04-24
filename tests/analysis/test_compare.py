"""Unit tests for :mod:`lysis.analysis.compare`.

Focus on the dispatch-table wiring for both ``micro-stats`` and
``macro-stats``:

* MEASURE_EXTRACTORS registers both keys with callable extractors.
* STATS_COMPUTERS registers both keys with callable scalar computers.
* ``macro-stats`` extractor returns the three expected labels with
  length-``n_sims`` arrays.
* ``macro-stats`` scalar computer returns the four expected labels and
  delegates to :func:`compute_run_statistics`.
* :func:`compare_runs` on two identical Runs yields KS statistic == 0
  and pct_diff == 0 across every registered label.
"""

from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from lysis.analysis.compare import (
    MEASURE_EXTRACTORS,
    STATS_COMPUTERS,
    _compare_arrays,
    _values_match,
    available_measure_sets,
    compare_runs,
    compute_stats,
    extract_measures,
    percent_difference,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


_MACRO_STATS_SERIES = pd.Series(
    {
        ("Degradation rate (%/min)", "Mean"): 12.5,
        ("Degradation rate (%/min)", "Standard Deviation"): 0.3,
        ("Lysis lag time (min)", "Mean"): 1.1,
        ("Lysis lag time (min)", "Standard Deviation"): 0.05,
        ("Time to full clot degradation (min)", "Mean"): 5.0,
        ("Time to full clot degradation (min)", "Standard Deviation"): 0.2,
        ("Percent of molecules that reached the back row", "Mean"): 48.0,
        ("Percent of molecules that reached the back row", "Standard Deviation"): 2.0,
        ("First passage time (min)", "Mean"): 3.7,
        ("First passage time (min)", "Standard Deviation"): 0.4,
        ("Front Velocity (microns/min)", "Mean"): 2.8,
        ("Front Velocity (microns/min)", "Standard Deviation"): 0.1,
    }
)


# ---------------------------------------------------------------------------
# Dispatch-table registration
# ---------------------------------------------------------------------------


class TestDispatchTables:
    """``macro-stats`` is registered in both dispatch tables."""

    def test_measure_extractors_has_macro_stats(self):
        assert "macro-stats" in MEASURE_EXTRACTORS
        assert callable(MEASURE_EXTRACTORS["macro-stats"])

    def test_stats_computers_has_macro_stats(self):
        assert "macro-stats" in STATS_COMPUTERS
        assert callable(STATS_COMPUTERS["macro-stats"])

    def test_micro_stats_still_registered(self):
        assert "micro-stats" in MEASURE_EXTRACTORS
        assert "micro-stats" in STATS_COMPUTERS


# ---------------------------------------------------------------------------
# Extractor
# ---------------------------------------------------------------------------


class TestExtractMacroStatsMeasures:
    """Tests for the ``macro-stats`` entry in :data:`MEASURE_EXTRACTORS`."""

    def test_returns_three_labels(self, stub_run):
        result = extract_measures(stub_run, "macro-stats")
        assert set(result.keys()) == {
            "Degradation rate (%/min)",
            "Time to full clot degradation (min)",
            "Front Velocity (microns/min)",
        }

    def test_arrays_have_n_sims_length(self, stub_run):
        result = extract_measures(stub_run, "macro-stats")
        n_sims = stub_run.macro_params.macro_simulations
        for label, arr in result.items():
            assert arr.shape == (n_sims,), f"{label} has wrong shape"

    def test_degradation_rate_in_percent_per_min(self, stub_run):
        """Degradation rate is multiplied by 100 to match the scalar units."""
        from lysis.analysis.degradation import mean_degradation_rate

        result = extract_measures(stub_run, "macro-stats")
        raw_rate, _, _ = mean_degradation_rate(stub_run)
        np.testing.assert_allclose(
            result["Degradation rate (%/min)"], raw_rate * 100
        )


# ---------------------------------------------------------------------------
# Scalar computer
# ---------------------------------------------------------------------------


class TestComputeMacroStatsScalars:
    """Tests for the ``macro-stats`` entry in :data:`STATS_COMPUTERS`."""

    def test_returns_four_labels(self, stub_run):
        with patch(
            "lysis.analysis.degradation.compute_run_statistics",
            return_value=_MACRO_STATS_SERIES,
        ):
            result = compute_stats(stub_run, "macro-stats")
        assert set(result.keys()) == {
            "Mean Degradation rate (%/min)",
            "Mean Time to full clot degradation (min)",
            "Mean First passage time (min)",
            "Mean Front Velocity (microns/min)",
        }

    def test_values_come_from_run_statistics_means(self, stub_run):
        with patch(
            "lysis.analysis.degradation.compute_run_statistics",
            return_value=_MACRO_STATS_SERIES,
        ):
            result = compute_stats(stub_run, "macro-stats")
        assert result["Mean Degradation rate (%/min)"] == pytest.approx(12.5)
        assert result["Mean Time to full clot degradation (min)"] == pytest.approx(5.0)
        assert result["Mean First passage time (min)"] == pytest.approx(3.7)
        assert result["Mean Front Velocity (microns/min)"] == pytest.approx(2.8)

    def test_values_are_python_floats(self, stub_run):
        with patch(
            "lysis.analysis.degradation.compute_run_statistics",
            return_value=_MACRO_STATS_SERIES,
        ):
            result = compute_stats(stub_run, "macro-stats")
        for v in result.values():
            assert isinstance(v, float)


# ---------------------------------------------------------------------------
# compare_runs
# ---------------------------------------------------------------------------


class TestCompareRunsMacroStats:
    """End-to-end dispatch for ``macro-stats`` through :func:`compare_runs`."""

    def test_identical_runs_give_zero_ks_and_zero_pct_diff(self, stub_run):
        """Comparing a Run to itself yields all-zero KS stats and pct diffs."""
        with patch(
            "lysis.analysis.degradation.compute_run_statistics",
            return_value=_MACRO_STATS_SERIES,
        ):
            result = compare_runs(stub_run, stub_run, "macro-stats")

        assert set(result.keys()) == {"ks", "pct_diff"}

        for label, ks in result["ks"].items():
            assert ks.statistic == pytest.approx(0.0), f"{label} KS != 0"

        for label, pct in result["pct_diff"].items():
            assert pct == pytest.approx(0.0), f"{label} pct_diff != 0"


# ---------------------------------------------------------------------------
# percent_difference
# ---------------------------------------------------------------------------


class TestPercentDifference:
    """Smoke tests for :func:`percent_difference` edge cases."""

    def test_equal_values_give_zero(self):
        assert percent_difference(5.0, 5.0) == pytest.approx(0.0)

    def test_both_zero_gives_zero(self):
        assert percent_difference(0.0, 0.0) == 0.0

    def test_opposite_sign_same_magnitude_is_nan(self):
        """v1 + v2 == 0 denominator → NaN."""
        assert np.isnan(percent_difference(1.0, -1.0))

    def test_nan_input_propagates(self):
        assert np.isnan(percent_difference(float("nan"), 1.0))
        assert np.isnan(percent_difference(1.0, float("nan")))

    def test_sign_follows_v2_greater_than_v1(self):
        assert percent_difference(1.0, 2.0) > 0
        assert percent_difference(2.0, 1.0) < 0


# ---------------------------------------------------------------------------
# _compare_arrays
# ---------------------------------------------------------------------------


class TestCompareArrays:
    def test_match_reports_zero(self):
        a = np.array([1.0, 2.0, 3.0])
        result = _compare_arrays(a, a.copy())
        assert result["status"] == "match"
        assert result["max_pct_diff"] == 0.0
        assert result["location"] is None
        assert result["mismatches"] == 0
        assert result["total"] == 3

    def test_shape_mismatch_reports_shapes(self):
        a = np.zeros((3,))
        b = np.zeros((4,))
        result = _compare_arrays(a, b)
        assert result["status"] == "shape_mismatch"
        assert result["detail"] == "(3,) vs (4,)"
        assert result["max_pct_diff"] is None
        assert result["location"] is None
        assert result["mismatches"] is None
        assert result["total"] is None

    def test_diff_reports_max_location(self):
        a = np.array([10.0, 10.0, 10.0])
        b = np.array([10.0, 10.0, 12.0])
        result = _compare_arrays(a, b)
        assert result["status"] == "diff"
        assert result["location"] == (2,)
        # symmetric pct diff of 10 vs 12: 200*(2)/22 ≈ 18.18
        assert result["max_pct_diff"] == pytest.approx(200 * 2 / 22)
        assert result["mismatches"] == 1
        assert result["total"] == 3

    def test_diff_counts_multiple_mismatches(self):
        a = np.array([1.0, 1.0, 1.0, 1.0])
        b = np.array([1.0, 2.0, 3.0, 1.0])
        result = _compare_arrays(a, b)
        assert result["mismatches"] == 2
        assert result["total"] == 4

    def test_structured_mismatch_count_is_records_not_fields(self):
        dtype = np.dtype([("t", np.float64), ("loc", np.int32)])
        a = np.array([(1.0, 5), (2.0, 7), (3.0, 9)], dtype=dtype)
        b = np.array([(1.0, 5), (2.5, 7), (3.0, 99)], dtype=dtype)
        result = _compare_arrays(a, b)
        # record 1 differs in 't', record 2 differs in 'loc' → 2 records
        assert result["mismatches"] == 2
        assert result["total"] == 3

    def test_2d_mismatch_count(self):
        a = np.ones((3, 4))
        b = a.copy()
        b[0, 0] = 2.0
        b[2, 3] = 2.0
        result = _compare_arrays(a, b)
        assert result["mismatches"] == 2
        assert result["total"] == 12

    def test_diff_sign_positive_when_b_larger(self):
        a = np.array([1.0])
        b = np.array([2.0])
        result = _compare_arrays(a, b)
        assert result["max_pct_diff"] > 0

    def test_diff_sign_negative_when_b_smaller(self):
        a = np.array([2.0])
        b = np.array([1.0])
        result = _compare_arrays(a, b)
        assert result["max_pct_diff"] < 0

    def test_structured_array_reports_field_in_location(self):
        dtype = np.dtype([("t", np.float64), ("loc", np.int32)])
        a = np.array([(1.0, 5), (2.0, 7), (3.0, 9)], dtype=dtype)
        b = np.array([(1.0, 5), (2.0, 7), (3.0, 99)], dtype=dtype)
        result = _compare_arrays(a, b)
        assert result["status"] == "diff"
        assert result["location"][0] == "loc"
        assert result["location"][1] == 2

    def test_2d_array_reports_multi_index_location(self):
        a = np.ones((3, 4))
        b = a.copy()
        b[1, 2] = 2.0
        result = _compare_arrays(a, b)
        assert result["status"] == "diff"
        assert result["location"] == (1, 2)


# ---------------------------------------------------------------------------
# N-ULP float tolerance
# ---------------------------------------------------------------------------


def _shift_ulps(a: np.ndarray, n: int) -> np.ndarray:
    """Return *a* shifted by *n* ULPs (toward +inf if *n* > 0, -inf if < 0)."""
    direction = np.inf if n >= 0 else -np.inf
    out = np.asarray(a, dtype=np.float64).copy()
    for _ in range(abs(n)):
        out = np.nextafter(out, direction)
    return out


class TestValuesMatchUlpTolerance:
    """:func:`_values_match` uses 2-ULP tolerance on float dtypes only."""

    def test_equal_floats_match(self):
        a = np.array([1.0, 2.0, 3.0])
        assert _values_match(a, a.copy()).all()

    def test_one_ulp_above_matches(self):
        a = np.array([1.0, 2.0, 3.0])
        assert _values_match(a, _shift_ulps(a, 1)).all()

    def test_one_ulp_below_matches(self):
        a = np.array([1.0, 2.0, 3.0])
        assert _values_match(a, _shift_ulps(a, -1)).all()

    def test_two_ulps_above_matches(self):
        a = np.array([1.0, 2.0, 3.0])
        assert _values_match(a, _shift_ulps(a, 2)).all()

    def test_two_ulps_below_matches(self):
        a = np.array([1.0, 2.0, 3.0])
        assert _values_match(a, _shift_ulps(a, -2)).all()

    def test_three_ulps_does_not_match(self):
        a = np.array([1.0, 2.0, 3.0])
        assert not _values_match(a, _shift_ulps(a, 3)).any()

    def test_nan_does_not_match_nan(self):
        a = np.array([np.nan, 1.0])
        b = np.array([np.nan, 1.0])
        assert _values_match(a, b).tolist() == [False, True]

    def test_signed_zero_matches(self):
        a = np.array([0.0])
        b = np.array([-0.0])
        assert _values_match(a, b).all()

    def test_inf_matches_itself(self):
        a = np.array([np.inf, -np.inf])
        assert _values_match(a, a.copy()).all()

    def test_integer_dtype_uses_exact_equality(self):
        # With float semantics, 2^53 and 2^53 + 2 are within 2 ULPs;
        # integers must NOT get that leniency.
        a = np.array([2**53], dtype=np.int64)
        b = np.array([2**53 + 2], dtype=np.int64)
        assert _values_match(a, b).tolist() == [False]

    def test_boolean_dtype_uses_exact_equality(self):
        a = np.array([True, False, True])
        b = np.array([True, True, True])
        assert _values_match(a, b).tolist() == [True, False, True]


class TestCompareArraysUlpTolerance:
    """:func:`_compare_arrays` treats up to 2 ULPs of float drift as a match."""

    def test_one_ulp_drift_reports_match(self):
        a = np.array([1.0, 2.0, 3.0])
        result = _compare_arrays(a, _shift_ulps(a, 1))
        assert result["status"] == "match"
        assert result["mismatches"] == 0
        assert result["max_pct_diff"] == 0.0
        assert result["location"] is None

    def test_two_ulp_drift_reports_match(self):
        a = np.array([1.0, 2.0, 3.0])
        result = _compare_arrays(a, _shift_ulps(a, 2))
        assert result["status"] == "match"
        assert result["mismatches"] == 0

    def test_mixed_tolerated_and_real_diff_reports_only_real(self):
        a = np.array([1.0, 2.0, 3.0])
        b = a.copy()
        b[0] = _shift_ulps(np.array([a[0]]), 2)[0]  # 2 ULPs — tolerated
        b[2] = 99.0  # real difference
        result = _compare_arrays(a, b)
        assert result["status"] == "diff"
        assert result["mismatches"] == 1
        assert result["location"] == (2,)

    def test_three_ulp_drift_still_reports_diff(self):
        a = np.array([10.0])
        result = _compare_arrays(a, _shift_ulps(a, 3))
        assert result["status"] == "diff"
        assert result["mismatches"] == 1
        assert result["location"] == (0,)

    def test_structured_array_float_field_tolerates_two_ulps(self):
        dtype = np.dtype([("t", np.float64), ("loc", np.int32)])
        a = np.array([(1.0, 5), (2.0, 7)], dtype=dtype)
        t_shifted = _shift_ulps(np.array([1.0]), 2)[0]
        b = np.array([(t_shifted, 5), (2.0, 7)], dtype=dtype)
        result = _compare_arrays(a, b)
        assert result["status"] == "match"
        assert result["mismatches"] == 0

    def test_structured_array_integer_field_remains_exact(self):
        dtype = np.dtype([("t", np.float64), ("loc", np.int32)])
        a = np.array([(1.0, 5)], dtype=dtype)
        b = np.array([(1.0, 6)], dtype=dtype)  # loc differs by 1 -> still a diff
        result = _compare_arrays(a, b)
        assert result["status"] == "diff"
        assert result["mismatches"] == 1
        assert result["location"][0] == "loc"


# ---------------------------------------------------------------------------
# available_measure_sets
# ---------------------------------------------------------------------------


class TestAvailableMeasureSets:
    def test_returns_only_stats_keys(self):
        assert set(available_measure_sets()) == {"micro-stats", "macro-stats"}
