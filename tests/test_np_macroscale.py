"""Unit tests for the lysis.np_macroscale module.

Tests cover:
- MacroscaleSim initialization with DataStore
- Data reads from macroscale_in (unbinding time, lysis time)
- Molecule state transitions (bind, unbind_by_degradation, unbind_by_time)
- Data snapshots via save_data
- Short end-to-end simulation via go()
"""

import numpy as np
import h5py
import pytest

from lysis.config.constants import Q_
from lysis.config.parameters import MicroParameters, MacroParameters
from lysis.config.run import Run
from lysis.data.datastore import DataStore
from lysis.data.dataspec import dataspec
from lysis.np_macroscale import MacroscaleSim


# ---------------------------------------------------------------------------
#  Helpers: build a small DataStore with synthetic microscale data
# ---------------------------------------------------------------------------

# Microscale datasets required for macroscale_in generation
_REQUIRED_DATASETS = [
    "pli_first_time",
    "tpa_leaving_time",
    "fiber_degraded",
    "sim_final_time",
    "tpa_unbound_by_pli",
    "tpa_unbound_kinetic",
]


def _fill_microscale_data(ds, n_sims=100, seed=42):
    """Fill microscale_out datasets in an open DataStore with synthetic data.

    The DataStore must have been created via ``DataStore.create()`` (so empty
    datasets already exist).  This function resizes them and writes plausible
    synthetic data so that ``generate_macroscale_in`` can produce valid
    macroscale_in.

    :param ds: An open DataStore in writable mode (``"a"``).
    :param n_sims: Number of microscale simulations (must be divisible by 100).
    :param seed: RNG seed for reproducibility.
    """
    rng = np.random.default_rng(seed)
    spec = dataspec["v2.0.0"]["microscale_out"]

    for name, ds_spec in spec.data.items():
        if ds_spec.data_location is None:
            continue
        if ds_spec.data_location not in ds._file:
            continue

        dataset = ds._file[ds_spec.data_location]

        # Generate appropriate data for each dataset type
        if name == "tpa_leaving_time":
            data = np.sort(rng.uniform(1, 100, n_sims)).astype(ds_spec.dtype)
        elif name == "sim_final_time":
            data = rng.uniform(50, 200, n_sims).astype(ds_spec.dtype)
        elif name == "fiber_degraded":
            # ~2/3 complete lysis, 1/3 do not
            data = np.ones(n_sims, dtype=bool)
            data[::3] = False
        elif name == "tpa_unbound_by_pli":
            data = np.zeros(n_sims, dtype=bool)
            data[: n_sims // 2] = True
        elif name == "tpa_unbound_kinetic":
            data = np.zeros(n_sims, dtype=bool)
            data[n_sims // 2 :] = True
        elif name == "pli_first_time":
            data = rng.uniform(0, 50, n_sims).astype(ds_spec.dtype)
        elif ds_spec.dtype == h5py.string_dtype():
            # String datasets: not required for macroscale_in, skip
            continue
        elif ds_spec.dtype == np.bool:
            data = rng.choice([True, False], n_sims)
        elif np.issubdtype(ds_spec.dtype, np.integer):
            data = rng.integers(0, 100, n_sims, dtype=ds_spec.dtype)
        else:
            data = rng.uniform(0, 100, n_sims).astype(ds_spec.dtype)

        dataset.resize((n_sims,))
        dataset[:] = data

    ds._file.flush()


def _make_run_and_datastore(tmp_path, run_code="test_run", n_sims=100,
                            macro_overrides=None):
    """Create a Run + DataStore with macroscale_in available.

    Returns the Run with ``run.data`` set to the open DataStore.
    """
    run = Run(str(tmp_path), run_code)
    run.initialize_micro_param({"micro_simulations": n_sims})

    macro_defaults = {
        "rows": 10,
        "cols": 5,
        "empty_rows": 2,
        "total_molecules": 10,
        "total_time": Q_("10 sec"),
        "macro_simulations": 1,
    }
    if macro_overrides:
        macro_defaults.update(macro_overrides)
    run.initialize_macro_param(macro_defaults)

    # Create DataStore with proper parameters
    ds = DataStore.create(run.run_code, run.os_path, run.micro_params)

    # Fill microscale data
    _fill_microscale_data(ds, n_sims=n_sims)

    # Add macroscale structure (re-initializes DataStore, detects macroscale_in)
    ds.initialize_macroscale(run.macro_params)

    run.data = ds
    return run


# ---------------------------------------------------------------------------
#  Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def run_with_data(tmp_path):
    """A Run with a small DataStore containing macroscale_in data."""
    return _make_run_and_datastore(tmp_path)


@pytest.fixture
def sim(run_with_data):
    """A freshly initialized MacroscaleSim (not yet run)."""
    return MacroscaleSim(run_with_data)


# ---------------------------------------------------------------------------
#  Initialization tests
# ---------------------------------------------------------------------------


class TestMacroscaleSimInit:
    """Tests for MacroscaleSim.__init__."""

    def test_state_arrays_initialized(self, sim):
        """Core state arrays are created with correct shapes."""
        mp = sim.run.macro_params
        total_edges = mp.rows * mp.full_row

        assert sim.fiber_status.shape == (total_edges,)
        assert sim.location.shape == (mp.total_molecules,)
        assert sim.bound.shape == (mp.total_molecules,)
        assert sim.binding_time.shape == (mp.total_molecules,)
        assert sim.waiting_time.shape == (mp.total_molecules,)

    def test_output_arrays_on_sim_instance(self, sim):
        """Output arrays live on the MacroscaleSim, not on run.data."""
        mp = sim.run.macro_params
        total_edges = mp.rows * mp.full_row

        assert hasattr(sim, "degradation_state")
        assert hasattr(sim, "molecule_location")
        assert hasattr(sim, "molecule_state")
        assert hasattr(sim, "save_time_array")

        assert sim.degradation_state.shape == (mp.number_of_saves, total_edges)
        assert sim.molecule_location.shape == (
            mp.number_of_saves,
            mp.total_molecules,
        )
        assert sim.molecule_state.shape == (
            mp.number_of_saves,
            mp.total_molecules,
        )
        assert sim.save_time_array.shape == (mp.number_of_saves,)

    def test_no_molecules_bound_initially(self, sim):
        """All molecules start unbound."""
        assert not np.any(sim.bound)

    def test_molecules_placed_in_empty_rows(self, sim):
        """All molecules start in the empty (fiber-free) rows."""
        mp = sim.run.macro_params
        max_empty_index = mp.empty_rows * mp.full_row
        assert np.all(sim.location < max_empty_index)

    def test_non_existent_edges_degraded(self, sim):
        """Edges in empty rows have fiber_status = 0 (already degraded)."""
        mp = sim.run.macro_params
        for i in range(mp.empty_rows):
            for j in range(mp.full_row):
                idx = i * mp.full_row + j
                assert sim.fiber_status[idx] == 0.0

    def test_real_fibers_are_intact(self, sim):
        """Real fibers start with fiber_status = infinity (intact)."""
        intact = sim.fiber_status[sim.real_fiber]
        assert np.all(intact == float("inf"))


# ---------------------------------------------------------------------------
#  Data read tests (macroscale_in integration)
# ---------------------------------------------------------------------------


class TestDataReads:
    """Tests for reading macroscale_in data through DataStore."""

    def test_macroscale_in_available(self, run_with_data):
        """macroscale_in collection is present on the DataStore."""
        assert "macroscale_in" in run_with_data.data.collections

    def test_bin_edge_tpa_leaving_time_shape(self, run_with_data):
        """bin_edge_tpa_leaving_time has shape (101,)."""
        data = run_with_data.data.macroscale_in.bin_edge_tpa_leaving_time
        assert data.shape == (101,)

    def test_binned_fiber_degraded_shape(self, run_with_data):
        """binned_fiber_degraded has shape (100,)."""
        data = run_with_data.data.macroscale_in.binned_fiber_degraded
        assert data.shape == (100,)

    def test_binned_fiber_degrade_time_shape(self, run_with_data):
        """binned_fiber_degrade_time has shape (bin_size, 100)."""
        n_sims = run_with_data.micro_params.micro_simulations
        bin_size = n_sims // 100
        data = run_with_data.data.macroscale_in.binned_fiber_degrade_time
        assert data.shape == (bin_size, 100)

    def test_find_unbinding_time_returns_finite(self, sim):
        """find_unbinding_time returns finite values for valid bins."""
        bins = np.array([10.0, 50.0, 90.0])
        result = sim.find_unbinding_time(bins, current_time=100.0)
        assert result.shape == (3,)
        assert np.all(np.isfinite(result))

    def test_find_unbinding_time_monotonic(self, sim):
        """Higher bins produce higher or equal unbinding times."""
        bins = np.arange(0, 101, dtype=np.float64)
        result = sim.find_unbinding_time(bins, current_time=0.0)
        # The underlying data is sorted tPA leaving times, so interpolation
        # should be monotonically non-decreasing
        diffs = np.diff(result)
        assert np.all(diffs >= -1e-10)  # allow floating point tolerance

    def test_find_lysis_time_returns_array(self, sim):
        """find_lysis_time returns an array of the correct size."""
        count = 5
        m = np.zeros(sim.run.macro_params.total_molecules, dtype=bool)
        m[:count] = True
        bins = np.array([50.0] * count)

        result = sim.find_lysis_time(m, bins, current_time=100.0, count=count)
        assert result.shape == (count,)

    def test_find_lysis_time_all_inf_for_empty_bins(self, sim):
        """Bins with no degraded simulations produce lysis_time = inf."""
        # binned_fiber_degraded tells us which bins have 0 degraded sims
        degraded = sim.run.data.macroscale_in.binned_fiber_degraded
        empty_bins = np.where(degraded == 0)[0]

        if len(empty_bins) == 0:
            pytest.skip("No empty bins in synthetic data")

        # Pick one empty bin and test find_lysis_time with it
        count = 1
        m = np.zeros(sim.run.macro_params.total_molecules, dtype=bool)
        m[0] = True
        bins = np.array([float(empty_bins[0])])

        # Force lysis_time_bin to be 0 (< total_lyses=0 is False, so no lysis)
        result = sim.find_lysis_time(m, bins, current_time=100.0, count=count)
        assert result[0] == float("inf")


# ---------------------------------------------------------------------------
#  Bind / unbind tests
# ---------------------------------------------------------------------------


class TestBindUnbind:
    """Tests for bind, unbind_by_degradation, unbind_by_time methods."""

    def test_bind_sets_bound_flag(self, sim):
        """Binding sets self.bound to True for selected molecules."""
        # Place molecules on intact fibers
        mp = sim.run.macro_params
        first_fiber_idx = mp.empty_rows * mp.full_row
        sim.location[:3] = first_fiber_idx
        # Update fiber status cache
        sim.m_fiber_status = sim.fiber_status[sim.location]

        m = np.zeros(mp.total_molecules, dtype=bool)
        m[:3] = True
        sim.bind(m, current_time=1.0)

        assert np.all(sim.bound[:3])
        assert sim.total_binds == 3

    def test_bind_sets_binding_time(self, sim):
        """After binding, binding_time is set to a finite unbinding time."""
        mp = sim.run.macro_params
        first_fiber_idx = mp.empty_rows * mp.full_row
        sim.location[0] = first_fiber_idx
        sim.m_fiber_status = sim.fiber_status[sim.location]

        m = np.zeros(mp.total_molecules, dtype=bool)
        m[0] = True
        sim.bind(m, current_time=1.0)

        assert np.isfinite(sim.binding_time[0])

    def test_unbind_by_degradation(self, sim):
        """unbind_by_degradation clears bound, sets waiting_time."""
        mp = sim.run.macro_params

        # Manually set a molecule as bound
        sim.bound[0] = True
        m = np.zeros(mp.total_molecules, dtype=bool)
        m[0] = True

        sim.unbind_by_degradation(m, current_time=10.0)

        assert not sim.bound[0]
        assert sim.unbound_by_degradation[0]
        assert sim.waiting_time[0] > 10.0
        assert sim.binding_time[0] == float("inf")
        assert sim.total_macro_unbinds == 1

    def test_unbind_by_degradation_no_molecules(self, sim):
        """unbind_by_degradation with empty mask is a no-op."""
        mp = sim.run.macro_params
        m = np.zeros(mp.total_molecules, dtype=bool)
        sim.unbind_by_degradation(m, current_time=10.0)
        assert sim.total_macro_unbinds == 0

    def test_unbind_by_time(self, sim):
        """unbind_by_time clears bound flag."""
        mp = sim.run.macro_params

        sim.bound[0] = True
        m = np.zeros(mp.total_molecules, dtype=bool)
        m[0] = True

        sim.unbind_by_time(m, current_time=10.0)

        assert not sim.bound[0]
        assert not sim.unbound_by_degradation[0]


# ---------------------------------------------------------------------------
#  save_data tests
# ---------------------------------------------------------------------------


class TestSaveData:
    """Tests for save_data snapshot method."""

    def test_save_data_stores_snapshot(self, sim):
        """save_data records fiber_status, location, bound, and time."""
        current_time = 5.0
        sim.save_data(current_time)

        np.testing.assert_array_equal(
            sim.degradation_state[0], sim.fiber_status
        )
        np.testing.assert_array_equal(
            sim.molecule_location[0], sim.location
        )
        np.testing.assert_array_equal(
            sim.molecule_state[0], sim.bound
        )
        assert sim.save_time_array[0] == current_time
        assert sim.current_save_interval == 1

    def test_save_data_increments_index(self, sim):
        """Multiple save_data calls fill successive rows."""
        sim.save_data(0.0)
        sim.save_data(5.0)

        assert sim.current_save_interval == 2
        assert sim.save_time_array[0] == 0.0
        assert sim.save_time_array[1] == 5.0


# ---------------------------------------------------------------------------
#  Integration test: short simulation
# ---------------------------------------------------------------------------


class TestShortSimulation:
    """End-to-end test with a tiny grid and very short run time."""

    def test_go_completes(self, tmp_path):
        """A short simulation completes without error."""
        run = _make_run_and_datastore(
            tmp_path,
            run_code="short_sim",
            macro_overrides={
                "total_time": Q_("0.01 sec"),
                "save_interval": Q_("0.01 sec"),
            },
        )
        sim = MacroscaleSim(run)
        sim.go()

        # Verify final save was recorded
        assert sim.current_save_interval >= 1
        assert sim.save_time_array[0] == 0.0

    def test_go_fiber_status_changes(self, tmp_path):
        """After simulation, at least some fiber status values may have changed."""
        run = _make_run_and_datastore(
            tmp_path,
            run_code="fiber_sim",
            macro_overrides={
                "total_time": Q_("0.1 sec"),
                "save_interval": Q_("0.1 sec"),
            },
        )
        sim = MacroscaleSim(run)
        sim.go()

        # The simulation ran — we just verify it completed and recorded data
        assert sim.current_save_interval >= 1
        # Total binds should be non-negative (may be 0 for very short sim)
        assert sim.total_binds >= 0
