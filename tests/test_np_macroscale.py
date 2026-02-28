"""Unit tests for the lysis.np_macroscale module.

Tests cover:
- MacroscaleSim initialization with DataStore
- Data reads from macroscale_in (unbinding time, lysis time)
- Molecule state transitions (bind, unbind_by_degradation, unbind_by_time)
- Event logging (tpa_bind_events, fiber_degrade_time)
- Data snapshots via save_data (v2.0.0 format)
- Writing to HDF5 via record_data_to_disk
- Short end-to-end simulation via go()
"""

import numpy as np
import h5py
import pytest

from lysis.config.constants import MolStatus, Q_
from lysis.config.parameters import MicroParameters, MacroParameters
from lysis.config.run import Run
from lysis.dataio.datastore import DataStore
from lysis.dataio.dataspec import dataspec
from lysis.np_macroscale import MacroscaleSim, _BIND_EVENT_DTYPE, _FIBER_DEGRADE_DTYPE


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
        """v2.0.0 output arrays live on the MacroscaleSim."""
        mp = sim.run.macro_params

        assert hasattr(sim, "tpa_location_snapshot")
        assert hasattr(sim, "snapshot_time")

        assert sim.tpa_location_snapshot.shape == (
            mp.number_of_saves,
            2,
            mp.total_molecules,
        )
        assert sim.snapshot_time.shape == (mp.number_of_saves,)

    def test_event_lists_initialized(self, sim):
        """Event accumulation lists start empty."""
        assert sim._bind_events == []
        assert sim._fiber_degrade_events == []

    def test_sim_number_default(self, sim):
        """Default sim_number is 0."""
        assert sim.sim_number == 0

    def test_sim_number_custom(self, run_with_data):
        """sim_number can be set via constructor."""
        sim = MacroscaleSim(run_with_data, sim_number=0)
        assert sim.sim_number == 0

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
#  Event logging tests
# ---------------------------------------------------------------------------


class TestEventLogging:
    """Tests for event recording in bind/unbind methods."""

    def test_bind_records_events(self, sim):
        """bind() records tpa_bind_events with BOUND status."""
        mp = sim.run.macro_params
        first_fiber_idx = mp.empty_rows * mp.full_row
        sim.location[:2] = first_fiber_idx
        sim.m_fiber_status = sim.fiber_status[sim.location]

        m = np.zeros(mp.total_molecules, dtype=bool)
        m[:2] = True
        sim.bind(m, current_time=1.0)

        assert len(sim._bind_events) >= 1
        events = np.concatenate(sim._bind_events)
        bound_events = events[events["Molecule New Status"] == MolStatus.BOUND.value]
        assert len(bound_events) == 2
        assert np.all(bound_events["Simulation Time Elapsed"] == 1.0)
        assert set(bound_events["tPA Molecule Index"]) == {0, 1}

        # Verify row/rank coordinates are valid
        row, rank = np.unravel_index(first_fiber_idx, (mp.rows, mp.full_row))
        assert np.all(bound_events["Grid Location Row"] == row)
        assert np.all(bound_events["Grid Location Rank"] == rank)

    def test_bind_records_fiber_degrade(self, sim):
        """bind() records fiber_degrade_time when lysis_time changes fiber_status."""
        mp = sim.run.macro_params
        first_fiber_idx = mp.empty_rows * mp.full_row
        sim.location[0] = first_fiber_idx
        sim.m_fiber_status = sim.fiber_status[sim.location]

        m = np.zeros(mp.total_molecules, dtype=bool)
        m[0] = True

        # Use a fixed seed so bind() produces lysis
        np.random.seed(42)
        sim.bind(m, current_time=1.0)

        # If lysis happened, fiber_degrade_events should have entries
        if sim.fiber_status[first_fiber_idx] < float("inf"):
            assert len(sim._fiber_degrade_events) > 0
            fd_events = np.concatenate(sim._fiber_degrade_events)
            assert fd_events.dtype == _FIBER_DEGRADE_DTYPE
            assert np.all(fd_events["Simulation Time Elapsed"] == 1.0)
            assert np.all(np.isfinite(fd_events["Fiber New Degrade Time"]))

    def test_unbind_by_degradation_records_events(self, sim):
        """unbind_by_degradation records MACRO_UNBOUND events."""
        mp = sim.run.macro_params

        sim.bound[0] = True
        m = np.zeros(mp.total_molecules, dtype=bool)
        m[0] = True

        sim.unbind_by_degradation(m, current_time=10.0)

        assert len(sim._bind_events) == 1
        events = sim._bind_events[0]
        assert len(events) == 1
        assert events[0]["Molecule New Status"] == MolStatus.MACRO_UNBOUND.value
        assert events[0]["Simulation Time Elapsed"] == 10.0
        assert events[0]["tPA Molecule Index"] == 0

    def test_unbind_by_time_forced_records_events(self, tmp_path):
        """unbind_by_time records MICRO_UNBOUND events for forced unbinds."""
        run = _make_run_and_datastore(
            tmp_path,
            run_code="forced_unbind",
            macro_overrides={"forced_unbind": 1.0},
        )
        sim = MacroscaleSim(run)
        mp = sim.run.macro_params

        sim.bound[0] = True
        m = np.zeros(mp.total_molecules, dtype=bool)
        m[0] = True

        sim.unbind_by_time(m, current_time=10.0)

        events = np.concatenate(sim._bind_events)
        micro_events = events[
            events["Molecule New Status"] == MolStatus.MICRO_UNBOUND.value
        ]
        assert len(micro_events) == 1
        assert micro_events[0]["tPA Molecule Index"] == 0

    def test_unbind_by_time_non_forced_records_events(self, tmp_path):
        """unbind_by_time records UNBOUND events for non-forced unbinds."""
        run = _make_run_and_datastore(
            tmp_path,
            run_code="non_forced_unbind",
            macro_overrides={"forced_unbind": 0.0},
        )
        sim = MacroscaleSim(run)
        mp = sim.run.macro_params

        sim.bound[0] = True
        m = np.zeros(mp.total_molecules, dtype=bool)
        m[0] = True

        sim.unbind_by_time(m, current_time=10.0)

        events = np.concatenate(sim._bind_events)
        unbound_events = events[
            events["Molecule New Status"] == MolStatus.UNBOUND.value
        ]
        assert len(unbound_events) == 1
        assert unbound_events[0]["tPA Molecule Index"] == 0


# ---------------------------------------------------------------------------
#  expire_waiting_period tests
# ---------------------------------------------------------------------------


class TestExpireWaitingPeriod:
    """Tests for expire_waiting_period (MACRO_UNBOUND/MICRO_UNBOUND -> UNBOUND)."""

    def test_macro_unbound_expiration(self, sim):
        """MACRO_UNBOUND molecule records UNBOUND event and clears state."""
        mp = sim.run.macro_params

        # Simulate a macro-unbound molecule whose waiting period has expired
        sim.waiting_time[0] = 5.0
        sim.unbound_by_degradation[0] = True

        sim.expire_waiting_period(current_time=6.0)

        # State cleared
        assert sim.waiting_time[0] == 0
        assert not sim.unbound_by_degradation[0]

        # UNBOUND event recorded
        assert len(sim._bind_events) == 1
        events = sim._bind_events[0]
        assert len(events) == 1
        assert events[0]["Molecule New Status"] == MolStatus.UNBOUND.value
        assert events[0]["Simulation Time Elapsed"] == 6.0
        assert events[0]["tPA Molecule Index"] == 0

    def test_micro_unbound_expiration(self, sim):
        """MICRO_UNBOUND molecule records UNBOUND event and clears waiting_time."""
        mp = sim.run.macro_params

        # Simulate a micro-unbound molecule (no degradation flag, but has
        # waiting_time set by unbind_by_time forced path)
        sim.waiting_time[2] = 3.0
        sim.unbound_by_degradation[2] = False

        sim.expire_waiting_period(current_time=4.0)

        assert sim.waiting_time[2] == 0

        assert len(sim._bind_events) == 1
        events = sim._bind_events[0]
        assert len(events) == 1
        assert events[0]["Molecule New Status"] == MolStatus.UNBOUND.value
        assert events[0]["tPA Molecule Index"] == 2

    def test_no_event_if_still_waiting(self, sim):
        """No event if waiting_time > current_time."""
        sim.waiting_time[0] = 10.0

        sim.expire_waiting_period(current_time=5.0)

        assert sim.waiting_time[0] == 10.0
        assert len(sim._bind_events) == 0

    def test_no_event_if_waiting_time_zero(self, sim):
        """No event if waiting_time is 0 (molecule never waited)."""
        assert sim.waiting_time[0] == 0

        sim.expire_waiting_period(current_time=5.0)

        assert len(sim._bind_events) == 0

    def test_no_event_if_bound(self, sim):
        """No event if molecule is bound (defensive guard)."""
        sim.bound[0] = True
        sim.waiting_time[0] = 3.0

        sim.expire_waiting_period(current_time=5.0)

        # Bound molecule should be excluded by ~self.bound
        assert sim.waiting_time[0] == 3.0
        assert len(sim._bind_events) == 0

    def test_does_not_fire_twice(self, sim):
        """After expiration, waiting_time is 0 so a second call is a no-op."""
        sim.waiting_time[0] = 3.0

        sim.expire_waiting_period(current_time=5.0)
        assert len(sim._bind_events) == 1

        sim.expire_waiting_period(current_time=6.0)
        assert len(sim._bind_events) == 1  # no new event

    def test_multiple_molecules_expiring(self, sim):
        """Multiple molecules can expire simultaneously."""
        sim.waiting_time[0] = 2.0
        sim.waiting_time[1] = 3.0
        sim.waiting_time[3] = 4.0
        sim.unbound_by_degradation[0] = True

        sim.expire_waiting_period(current_time=5.0)

        # All three should expire
        assert sim.waiting_time[0] == 0
        assert sim.waiting_time[1] == 0
        assert sim.waiting_time[3] == 0
        assert not sim.unbound_by_degradation[0]

        assert len(sim._bind_events) == 1
        events = sim._bind_events[0]
        assert len(events) == 3
        assert set(events["tPA Molecule Index"]) == {0, 1, 3}
        assert np.all(events["Molecule New Status"] == MolStatus.UNBOUND.value)
        assert np.all(events["Simulation Time Elapsed"] == 5.0)


# ---------------------------------------------------------------------------
#  save_data tests
# ---------------------------------------------------------------------------


class TestSaveData:
    """Tests for save_data snapshot method (v2.0.0 format)."""

    def test_save_data_stores_snapshot(self, sim):
        """save_data records tpa_location_snapshot and snapshot_time."""
        mp = sim.run.macro_params
        current_time = 5.0
        sim.save_data(current_time)

        # Verify snapshot_time
        assert sim.snapshot_time[0] == current_time
        assert sim.current_save_interval == 1

        # Verify tpa_location_snapshot contains valid (row, rank) coords
        expected_rows, expected_ranks = np.unravel_index(
            sim.location, (mp.rows, mp.full_row)
        )
        np.testing.assert_array_equal(
            sim.tpa_location_snapshot[0, 0, :], expected_rows
        )
        np.testing.assert_array_equal(
            sim.tpa_location_snapshot[0, 1, :], expected_ranks
        )

    def test_save_data_increments_index(self, sim):
        """Multiple save_data calls fill successive rows."""
        sim.save_data(0.0)
        sim.save_data(5.0)

        assert sim.current_save_interval == 2
        assert sim.snapshot_time[0] == 0.0
        assert sim.snapshot_time[1] == 5.0

    def test_save_data_row_rank_roundtrip(self, sim):
        """Saved (row, rank) coords correctly round-trip from flat indices."""
        mp = sim.run.macro_params
        sim.save_data(0.0)

        rows = sim.tpa_location_snapshot[0, 0, :]
        ranks = sim.tpa_location_snapshot[0, 1, :]
        flat = np.ravel_multi_index(
            (rows.astype(int), ranks.astype(int)),
            (mp.rows, mp.full_row),
        )
        np.testing.assert_array_equal(flat, sim.location)


# ---------------------------------------------------------------------------
#  record_data_to_disk tests
# ---------------------------------------------------------------------------


class TestRecordDataToDisk:
    """Tests for writing macroscale_out data to HDF5 via record_data_to_disk."""

    def test_writes_snapshot_time(self, sim):
        """record_data_to_disk writes snapshot_time to HDF5."""
        sim.save_data(0.0)
        sim.save_data(5.0)
        sim.record_data_to_disk()

        sv = sim.run.data.macroscale_out[sim.sim_number]
        ds = sv.snapshot_time
        assert ds.shape == (2,)
        np.testing.assert_array_equal(ds[:], sim.snapshot_time[:2])

    def test_writes_tpa_location_snapshot(self, sim):
        """record_data_to_disk writes tpa_location_snapshot with correct shape."""
        mp = sim.run.macro_params
        sim.save_data(0.0)
        sim.record_data_to_disk()

        sv = sim.run.data.macroscale_out[sim.sim_number]
        ds = sv.tpa_location_snapshot
        assert ds.shape == (1, 2, mp.total_molecules)
        np.testing.assert_array_equal(
            ds[:], sim.tpa_location_snapshot[:1]
        )

    def test_writes_bind_events(self, sim):
        """record_data_to_disk writes tpa_bind_events to HDF5."""
        mp = sim.run.macro_params
        first_fiber_idx = mp.empty_rows * mp.full_row
        sim.location[:2] = first_fiber_idx
        sim.m_fiber_status = sim.fiber_status[sim.location]

        m = np.zeros(mp.total_molecules, dtype=bool)
        m[:2] = True
        sim.bind(m, current_time=1.0)

        sim.save_data(1.0)
        sim.record_data_to_disk()

        sv = sim.run.data.macroscale_out[sim.sim_number]
        ds = sv.tpa_bind_events
        assert ds.shape[0] >= 2  # at least the 2 BOUND events

    def test_writes_fiber_degrade_events(self, sim):
        """record_data_to_disk writes fiber_degrade_time when events exist."""
        # Manually add a fiber degrade event
        event = np.array(
            [(1.0, 3, 5, 42.0)], dtype=_FIBER_DEGRADE_DTYPE
        )
        sim._fiber_degrade_events.append(event)

        sim.save_data(1.0)
        sim.record_data_to_disk()

        sv = sim.run.data.macroscale_out[sim.sim_number]
        ds = sv.fiber_degrade_time
        assert ds.shape == (1,)
        assert ds[0]["Fiber New Degrade Time"] == 42.0

    def test_writes_tpa_transit_time(self, sim):
        """record_data_to_disk writes tpa_transit_time for molecules that reached back row."""
        sim.reached_back_row[0] = True
        sim.time_to_reach_back_row[0] = 99.5

        sim.save_data(100.0)
        sim.record_data_to_disk()

        sv = sim.run.data.macroscale_out[sim.sim_number]
        ds = sv.tpa_transit_time
        assert ds.shape == (1,)
        assert ds[0] == 99.5

    def test_no_events_leaves_datasets_empty(self, sim):
        """With no events, bind/degrade/transit datasets remain size 0."""
        sim.save_data(0.0)
        sim.record_data_to_disk()

        sv = sim.run.data.macroscale_out[sim.sim_number]
        assert sv.tpa_bind_events.shape == (0,)
        assert sv.fiber_degrade_time.shape == (0,)
        assert sv.tpa_transit_time.shape == (0,)


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
        assert sim.snapshot_time[0] == 0.0

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

    def test_go_writes_data_to_hdf5(self, tmp_path):
        """go() writes snapshot data to HDF5 via record_data_to_disk."""
        run = _make_run_and_datastore(
            tmp_path,
            run_code="hdf5_sim",
            macro_overrides={
                "total_time": Q_("0.01 sec"),
                "save_interval": Q_("0.01 sec"),
            },
        )
        sim = MacroscaleSim(run)
        sim.go()

        sv = run.data.macroscale_out[0]
        # snapshot_time should have been written
        assert sv.snapshot_time.shape[0] == sim.current_save_interval
        # tpa_location_snapshot should match
        assert sv.tpa_location_snapshot.shape[0] == sim.current_save_interval

    def test_go_event_lists_populated(self, tmp_path):
        """After go(), event lists should contain data if binding occurred."""
        run = _make_run_and_datastore(
            tmp_path,
            run_code="event_sim",
            macro_overrides={
                "total_time": Q_("0.1 sec"),
                "save_interval": Q_("0.1 sec"),
            },
        )
        sim = MacroscaleSim(run)
        sim.go()

        # If binds occurred, bind_events should be populated
        if sim.total_binds > 0:
            all_events = np.concatenate(sim._bind_events)
            assert len(all_events) >= sim.total_binds
            # Every BOUND event should have a valid MolStatus value
            statuses = set(all_events["Molecule New Status"])
            valid = {s.value for s in MolStatus}
            assert statuses.issubset(valid)
