"""Shared fixtures for lysis.analysis unit tests.

Grid geometry (rows=5, cols=3, empty_rows=1, macro_simulations=2):

    full_row    = 3*3-1  = 8
    xz_row      = 2*3-1  = 5
    total_edges = 8*4+5  = 37
    empty_edges = 8*1    = 8
    total_fibers= 8*3+5  = 29
    fiber_rows  = 5-1    = 4

Synthetic degrade-time schedule (same for both simulations):

    Fiber indices 0-7   (edges 8-15)  degrade at t=50s
    Fiber indices 8-14  (edges 16-22) degrade at t=100s
    Fiber indices 15-21 (edges 23-29) degrade at t=175s
    Fiber indices 22-28 (edges 30-36) degrade at t=250s

With tsave = [0, 75, 150, 225, 300] seconds this produces
degraded_fraction values of [0, 8/29, 15/29, 22/29, 1].

Marker frames for percent_markers=[0.0, 0.25, 0.50, 0.75, 1.0]
are [0, 1, 2, 3, 4] for both simulations.
"""

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from lysis.config.constants import Q_
from lysis.config.parameters import MacroParameters, MicroParameters

# ---------------------------------------------------------------------------
# Analysis configuration constants used across multiple test classes
# ---------------------------------------------------------------------------

PERCENT_MARKERS = [0.0, 0.25, 0.50, 0.75, 1.0]
SLOPE_PAIRS = [(0.25, 0.75)]
FRONT_THRESHOLD = 0.5

# Degrade-time sentinel for edges not yet degraded
SENTINEL = 9.9e100


# ---------------------------------------------------------------------------
# Parameters fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def small_macro():
    """MacroParameters with a minimal 5×3 grid and 2 simulations."""
    micro = MicroParameters()
    return MacroParameters(
        micro_params=micro,
        rows=5,
        cols=3,
        empty_rows=1,
        macro_simulations=2,
        total_time=Q_("300 sec"),
        save_interval=Q_("75 sec"),
    )


# ---------------------------------------------------------------------------
# Stub Run
# ---------------------------------------------------------------------------


class _StubRun:
    """Minimal object exposing the Run attributes used by analysis functions."""

    def __init__(self, macro_params, os_path="."):
        self.macro_params = macro_params
        self.os_path = os_path


@pytest.fixture
def stub_run(small_macro):
    """Minimal Run-like stub carrying the small_macro parameters."""
    return _StubRun(small_macro)


# ---------------------------------------------------------------------------
# Synthetic simulation data
# ---------------------------------------------------------------------------


@pytest.fixture
def tsave_arrays(small_macro):
    """Per-simulation save-point arrays in seconds: [0, 75, 150, 225, 300]."""
    arr = np.array([0.0, 75.0, 150.0, 225.0, 300.0])
    return [arr.copy() for _ in range(small_macro.macro_simulations)]


@pytest.fixture
def deg_arrays(small_macro, tsave_arrays):
    """Synthetic fiber degrade-time arrays for a controlled degradation scenario.

    ``deg[sim]`` has shape ``(5, 37)`` with the following schedule:

    - Fiber indices  0-7  (edges  8-15): degrade at t=50s
    - Fiber indices  8-14 (edges 16-22): degrade at t=100s
    - Fiber indices 15-21 (edges 23-29): degrade at t=175s
    - Fiber indices 22-28 (edges 30-36): degrade at t=250s
    """
    total_edges = small_macro.total_edges  # 37
    empty_edges = small_macro.empty_edges  # 8
    n_sims = small_macro.macro_simulations  # 2
    n_saves = len(tsave_arrays[0])  # 5

    fiber_degrade_times = np.empty(total_edges - empty_edges, dtype=np.float64)
    fiber_degrade_times[0:8] = 50.0
    fiber_degrade_times[8:15] = 100.0
    fiber_degrade_times[15:22] = 175.0
    fiber_degrade_times[22:29] = 250.0

    result = []
    for _ in range(n_sims):
        d = np.full((n_saves, total_edges), SENTINEL, dtype=np.float64)
        # Empty edges are permanently "degraded" (degrade time = 0)
        d[:, :empty_edges] = 0.0
        # Fiber edges: record actual degrade time once tsave >= degrade_time
        for t_idx, t in enumerate(tsave_arrays[0]):
            mask = fiber_degrade_times <= t
            d[t_idx, empty_edges:][mask] = fiber_degrade_times[mask]
        result.append(d)
    return result


# ---------------------------------------------------------------------------
# Pre-computed derived fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def degraded_fraction(stub_run, deg_arrays, tsave_arrays):
    """Result of find_degraded_fraction for the synthetic data."""
    from lysis.analysis.degradation import find_degraded_fraction

    return find_degraded_fraction(stub_run, deg_arrays, tsave_arrays)


@pytest.fixture
def marker_frames(degraded_fraction):
    """Result of find_degradation_marker_frames for PERCENT_MARKERS."""
    from lysis.analysis.degradation import find_degradation_marker_frames

    return find_degradation_marker_frames(degraded_fraction, PERCENT_MARKERS)


@pytest.fixture
def marker_times(marker_frames, tsave_arrays):
    """Result of find_degradation_marker_times."""
    from lysis.analysis.degradation import find_degradation_marker_times

    return find_degradation_marker_times(marker_frames, tsave_arrays)


@pytest.fixture
def rates(stub_run, marker_frames, degraded_fraction, tsave_arrays):
    """Result of degradation_rates for SLOPE_PAIRS."""
    from lysis.analysis.degradation import degradation_rates

    return degradation_rates(
        stub_run, marker_frames, degraded_fraction, tsave_arrays,
        SLOPE_PAIRS, PERCENT_MARKERS,
    )


@pytest.fixture
def exposed_time(stub_run, deg_arrays):
    """Result of calculate_time_row_exposed for the synthetic data."""
    from lysis.analysis.degradation import calculate_time_row_exposed

    return calculate_time_row_exposed(stub_run, deg_arrays)


@pytest.fixture
def deg_fronts(stub_run, exposed_time, tsave_arrays):
    """Result of find_degradation_fronts for the synthetic data."""
    from lysis.analysis.degradation import find_degradation_fronts

    return find_degradation_fronts(stub_run, exposed_time, tsave_arrays)


@pytest.fixture
def row_deg(stub_run, deg_arrays, tsave_arrays):
    """Result of find_row_deg_fraction for the synthetic data."""
    from lysis.analysis.degradation import find_row_deg_fraction

    return find_row_deg_fraction(stub_run, deg_arrays, tsave_arrays)
