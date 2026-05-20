"""Tests for ``lysis.analysis.microscale`` — microscale statistics computation."""

import numpy as np
import pytest
from unittest.mock import MagicMock

from lysis.analysis.microscale import compute_micro_statistics


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_run(fiber_degraded, sim_final_time, tpa_leaving_time):
    """Build a mock Run with the required microscale_out datasets.

    The function accesses each dataset via ``dataset[:]``, so we configure
    ``__getitem__`` to return the numpy array regardless of the slice argument.
    """
    run = MagicMock()
    micro = run.data.microscale_out
    micro.fiber_degraded.__getitem__.return_value = fiber_degraded
    micro.sim_final_time.__getitem__.return_value = sim_final_time
    micro.tpa_leaving_time.__getitem__.return_value = tpa_leaving_time
    return run


# Reference data: 5 fibers, 3 degraded (indices 0, 1, 3)
_FD = np.array([True, True, False, True, False])
_SFT = np.array([120.0, 180.0, 300.0, 240.0, 360.0])  # seconds; degraded: 120, 180, 240
_TPA = np.array([10.0, 20.0, 30.0, 40.0, 50.0])       # seconds


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestComputeMicroStatistics:
    def test_returns_dict(self):
        run = _make_run(_FD, _SFT, _TPA)
        assert isinstance(compute_micro_statistics(run), dict)

    def test_result_has_all_keys(self):
        run = _make_run(_FD, _SFT, _TPA)
        result = compute_micro_statistics(run)
        expected = {
            "fibers_degraded",
            "lysis_time_mean",
            "lysis_time_std",
            "lysis_time_median",
            "tpa_leaving_mean",
            "tpa_leaving_std",
            "tpa_leaving_median",
        }
        assert set(result.keys()) == expected

    def test_fibers_degraded_count(self):
        run = _make_run(_FD, _SFT, _TPA)
        assert compute_micro_statistics(run)["fibers_degraded"] == 3

    def test_fibers_degraded_is_int(self):
        run = _make_run(_FD, _SFT, _TPA)
        assert isinstance(compute_micro_statistics(run)["fibers_degraded"], int)

    def test_lysis_time_mean_in_minutes(self):
        # Degraded sim_final_time: 120, 180, 240 s → 2, 3, 4 min → mean = 3.0
        run = _make_run(_FD, _SFT, _TPA)
        assert compute_micro_statistics(run)["lysis_time_mean"] == pytest.approx(3.0)

    def test_lysis_time_std_in_minutes(self):
        run = _make_run(_FD, _SFT, _TPA)
        expected = float(np.std(np.array([120.0, 180.0, 240.0]) / 60))
        assert compute_micro_statistics(run)["lysis_time_std"] == pytest.approx(expected)

    def test_lysis_time_median_in_minutes(self):
        # Median of [2, 3, 4] = 3.0
        run = _make_run(_FD, _SFT, _TPA)
        assert compute_micro_statistics(run)["lysis_time_median"] == pytest.approx(3.0)

    def test_tpa_leaving_mean_in_seconds(self):
        # Mean of [10, 20, 30, 40, 50] = 30.0
        run = _make_run(_FD, _SFT, _TPA)
        assert compute_micro_statistics(run)["tpa_leaving_mean"] == pytest.approx(30.0)

    def test_tpa_leaving_std_in_seconds(self):
        run = _make_run(_FD, _SFT, _TPA)
        expected = float(np.std(_TPA))
        assert compute_micro_statistics(run)["tpa_leaving_std"] == pytest.approx(expected)

    def test_tpa_leaving_median_in_seconds(self):
        # Median of [10, 20, 30, 40, 50] = 30.0
        run = _make_run(_FD, _SFT, _TPA)
        assert compute_micro_statistics(run)["tpa_leaving_median"] == pytest.approx(30.0)

    def test_all_fibers_degraded(self):
        fd = np.array([True, True, True])
        sft = np.array([60.0, 120.0, 180.0])
        tpa = np.array([5.0, 10.0, 15.0])
        run = _make_run(fd, sft, tpa)
        result = compute_micro_statistics(run)
        assert result["fibers_degraded"] == 3
        assert result["lysis_time_mean"] == pytest.approx(2.0)  # (1+2+3)/3 min

    def test_single_fiber_degraded(self):
        fd = np.array([True, False, False])
        sft = np.array([300.0, 600.0, 900.0])
        tpa = np.array([5.0, 10.0, 15.0])
        run = _make_run(fd, sft, tpa)
        result = compute_micro_statistics(run)
        assert result["fibers_degraded"] == 1
        assert result["lysis_time_mean"] == pytest.approx(5.0)  # 300 / 60
        assert result["lysis_time_std"] == pytest.approx(0.0)
        assert result["lysis_time_median"] == pytest.approx(5.0)

    def test_no_fibers_degraded_count(self):
        fd = np.array([False, False, False])
        sft = np.array([60.0, 120.0, 180.0])
        tpa = np.array([5.0, 10.0, 15.0])
        run = _make_run(fd, sft, tpa)
        result = compute_micro_statistics(run)
        assert result["fibers_degraded"] == 0

    def test_no_fibers_degraded_lysis_times_nan(self):
        fd = np.array([False, False, False])
        sft = np.array([60.0, 120.0, 180.0])
        tpa = np.array([5.0, 10.0, 15.0])
        run = _make_run(fd, sft, tpa)
        result = compute_micro_statistics(run)
        # numpy mean/std/median of empty array returns nan
        assert np.isnan(result["lysis_time_mean"])
        assert np.isnan(result["lysis_time_median"])

    def test_tpa_leaving_uses_all_simulations(self):
        # tpa_leaving_time is computed over all simulations, not just degraded ones
        fd = np.array([True, False])
        sft = np.array([60.0, 120.0])
        tpa = np.array([10.0, 30.0])  # both included
        run = _make_run(fd, sft, tpa)
        result = compute_micro_statistics(run)
        assert result["tpa_leaving_mean"] == pytest.approx(20.0)
