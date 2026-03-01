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
from lysis.geometry.edge_grid import from_fortran_edge_index

# ---------------------------------------------------------------------------
# Analysis configuration constants used across multiple test classes
# ---------------------------------------------------------------------------

PERCENT_MARKERS = [0.0, 0.25, 0.50, 0.75, 1.0]
SLOPE_PAIRS = [(0.25, 0.75)]
FRONT_THRESHOLD = 0.5

# Event log structured-array dtype (matches v2.0.0 DataSpec)
_EVENT_DTYPE = np.dtype(
    [
        ("Simulation Time Elapsed", np.float64),
        ("Grid Location Row", np.uint32),
        ("Grid Location Rank", np.uint32),
        ("Fiber New Degrade Time", np.float64),
    ]
)


# ---------------------------------------------------------------------------
# Mock DataStore infrastructure
# ---------------------------------------------------------------------------


class _StubDataset:
    """Wrap a numpy array to support h5py-like slice access."""

    def __init__(self, data):
        self._data = data

    def __getitem__(self, key):
        return self._data[key]


class _StubSimView:
    """Per-simulation view exposing snapshot_time and fiber_degrade_time."""

    def __init__(self, snapshot_time, fiber_degrade_time):
        self.snapshot_time = _StubDataset(snapshot_time)
        self.fiber_degrade_time = _StubDataset(fiber_degrade_time)


class _StubMacroOut:
    """Indexable container of _StubSimView objects."""

    def __init__(self, views):
        self._views = views

    def __getitem__(self, idx):
        return self._views[idx]


class _StubData:
    """Minimal DataStore-like object with macroscale_out attribute."""

    def __init__(self, macroscale_out):
        self.macroscale_out = macroscale_out


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

    def __init__(self, macro_params, os_path=".", data=None):
        self.macro_params = macro_params
        self.os_path = os_path
        self.data = data


# ---------------------------------------------------------------------------
# Synthetic event log construction
# ---------------------------------------------------------------------------


def _build_event_log(rows, cols):
    """Build a synthetic fiber_degrade_time event log from the test schedule.

    Converts each Fortran edge index (8-36) to (row, rank) coordinates
    and creates one event per fiber, sorted by time.

    Returns a structured array with _EVENT_DTYPE.
    """
    # Degrade schedule: Fortran edge index ranges → degrade times
    schedule = [
        (range(8, 16), 50.0),   # fiber indices 0-7
        (range(16, 23), 100.0),  # fiber indices 8-14
        (range(23, 30), 175.0),  # fiber indices 15-21
        (range(30, 37), 250.0),  # fiber indices 22-28
    ]

    events = []
    for edge_range, degrade_time in schedule:
        for k in edge_range:
            row, rank = from_fortran_edge_index(k, rows, cols)
            events.append((degrade_time, row, rank, degrade_time))

    event_log = np.array(events, dtype=_EVENT_DTYPE)
    event_log.sort(order="Simulation Time Elapsed")
    return event_log


@pytest.fixture
def stub_run(small_macro):
    """Minimal Run-like stub with mock DataStore wired to synthetic event logs."""
    rows = small_macro.rows
    cols = small_macro.cols
    n_sims = small_macro.macro_simulations

    snapshot_time = np.array([0.0, 75.0, 150.0, 225.0, 300.0])
    event_log = _build_event_log(rows, cols)

    views = [
        _StubSimView(snapshot_time.copy(), event_log.copy())
        for _ in range(n_sims)
    ]
    data = _StubData(_StubMacroOut(views))
    return _StubRun(small_macro, data=data)


# ---------------------------------------------------------------------------
# Pre-computed derived fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def degraded_fraction(stub_run):
    """Result of find_degraded_fraction for the synthetic data."""
    from lysis.analysis.degradation import find_degraded_fraction

    return find_degraded_fraction(stub_run)


@pytest.fixture
def marker_frames(degraded_fraction):
    """Result of find_degradation_marker_frames for PERCENT_MARKERS."""
    from lysis.analysis.degradation import find_degradation_marker_frames

    return find_degradation_marker_frames(degraded_fraction, PERCENT_MARKERS)


@pytest.fixture
def marker_times(stub_run, marker_frames):
    """Result of find_degradation_marker_times."""
    from lysis.analysis.degradation import find_degradation_marker_times

    return find_degradation_marker_times(stub_run, marker_frames)


@pytest.fixture
def rates(stub_run, marker_frames, degraded_fraction):
    """Result of degradation_rates for SLOPE_PAIRS."""
    from lysis.analysis.degradation import degradation_rates

    return degradation_rates(
        stub_run, marker_frames, degraded_fraction,
        SLOPE_PAIRS, PERCENT_MARKERS,
    )


@pytest.fixture
def exposed_time(stub_run):
    """Result of calculate_time_row_exposed for the synthetic data."""
    from lysis.analysis.degradation import calculate_time_row_exposed

    return calculate_time_row_exposed(stub_run)


@pytest.fixture
def deg_fronts(stub_run, exposed_time):
    """Result of find_degradation_fronts for the synthetic data."""
    from lysis.analysis.degradation import find_degradation_fronts

    return find_degradation_fronts(stub_run, exposed_time)


@pytest.fixture
def row_deg(stub_run):
    """Result of find_row_deg_fraction for the synthetic data."""
    from lysis.analysis.degradation import find_row_deg_fraction

    return find_row_deg_fraction(stub_run)
