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
