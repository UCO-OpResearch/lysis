"""Comprehensive pytest unit tests for :mod:`lysis.data.dataconvert`.

Tests cover:

* safe_np_int_conversion — bounds checking and float-to-int
* safe_np_bool_conversion — strict 0/1 validation
* safe_np_string_conversion — string-to-string only
* replay_event_log_to_snapshot — event log to time-series
* convert_structured_grid_fields — 1D/2D grid index conversion
* convert_location_snapshot — molecule location shape/index transform
* convert_bind_events_to_bound — event log to bound status snapshots
* convert_data — main entry point, tag resolution, params copy
"""

import numpy as np
import pytest

from lysis.config.constants import CONST
from lysis.data.dataspec import DataSetSpec, DataCollectionSpec, dataspec, tags
from lysis.data.dataconvert import (
    safe_np_int_conversion,
    safe_np_bool_conversion,
    safe_np_string_conversion,
    replay_event_log_to_snapshot,
    convert_structured_grid_fields,
    convert_location_snapshot,
    convert_bind_events_to_bound,
    convert_data,
)
from lysis.geometry.edge_grid import (
    from_fortran_edge_index_array,
    to_fortran_edge_index_array,
)


# ---------------------------------------------------------------------------
# safe_np_int_conversion
# ---------------------------------------------------------------------------


class TestSafeNpIntConversion:
    """Tests for :func:`safe_np_int_conversion`."""

    @pytest.mark.parametrize(
        "values, in_dtype, target, expected_values",
        [
            pytest.param([1, 2, 3], np.int64, np.int32, [1, 2, 3], id="int64_to_int32"),
            pytest.param([0, 127, 255], np.int64, np.uint8, [0, 127, 255], id="int64_to_uint8"),
            pytest.param([0, 100], np.int32, np.uint8, [0, 100], id="int32_to_uint8_signed_unsigned"),
            pytest.param([1, 2], np.int32, np.int32, [1, 2], id="identity_same_type"),
        ],
    )
    def test_successful_conversion(self, values, in_dtype, target, expected_values):
        """Values within target range convert successfully."""
        arr = np.array(values, dtype=in_dtype)
        result = safe_np_int_conversion(arr, dtype=target)
        np.testing.assert_array_equal(result, expected_values)
        assert result.dtype == np.dtype(target)

    def test_float_with_integer_values(self):
        """Float array containing whole numbers converts to int."""
        arr = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        result = safe_np_int_conversion(arr, dtype=np.int32)
        np.testing.assert_array_equal(result, [1, 2, 3])
        assert result.dtype == np.int32

    def test_float_with_non_integer_values_raises(self):
        """Float array with fractional values raises TypeError."""
        arr = np.array([1.5, 2.0], dtype=np.float64)
        with pytest.raises(TypeError, match="non-integer"):
            safe_np_int_conversion(arr, dtype=np.int32)

    @pytest.mark.parametrize(
        "values, in_dtype, target",
        [
            pytest.param([2**31], np.int64, np.int32, id="int64_overflow_int32"),
            pytest.param([-1], np.int64, np.uint8, id="negative_to_uint8"),
            pytest.param([256], np.int64, np.uint8, id="256_to_uint8"),
        ],
    )
    def test_overflow_raises(self, values, in_dtype, target):
        """Values outside target range raise OverflowError."""
        arr = np.array(values, dtype=in_dtype)
        with pytest.raises(OverflowError):
            safe_np_int_conversion(arr, dtype=target)

    def test_empty_array(self):
        """Empty array converts successfully."""
        arr = np.array([], dtype=np.int64)
        result = safe_np_int_conversion(arr, dtype=np.uint8)
        assert result.size == 0
        assert result.dtype == np.uint8


# ---------------------------------------------------------------------------
# safe_np_bool_conversion
# ---------------------------------------------------------------------------


class TestSafeNpBoolConversion:
    """Tests for :func:`safe_np_bool_conversion`."""

    def test_valid_01_values(self):
        """Array of 0s and 1s converts to boolean."""
        arr = np.array([0, 1, 0, 1], dtype=np.int32)
        result = safe_np_bool_conversion(arr)
        np.testing.assert_array_equal(result, [False, True, False, True])

    def test_all_zeros(self):
        """All zeros convert to all False."""
        arr = np.array([0, 0, 0], dtype=np.int32)
        result = safe_np_bool_conversion(arr)
        np.testing.assert_array_equal(result, [False, False, False])

    def test_value_2_raises(self):
        """Value 2 (not 0 or 1) raises OverflowError."""
        arr = np.array([0, 1, 2], dtype=np.int32)
        with pytest.raises(OverflowError):
            safe_np_bool_conversion(arr)

    def test_negative_raises(self):
        """Negative values raise OverflowError."""
        arr = np.array([-1, 0, 1], dtype=np.int32)
        with pytest.raises(OverflowError):
            safe_np_bool_conversion(arr)

    def test_empty_array(self):
        """Empty array converts successfully."""
        arr = np.array([], dtype=np.int32)
        result = safe_np_bool_conversion(arr)
        assert result.size == 0

    def test_already_bool(self):
        """Already-boolean array passes through."""
        arr = np.array([True, False, True])
        result = safe_np_bool_conversion(arr)
        np.testing.assert_array_equal(result, arr)


# ---------------------------------------------------------------------------
# safe_np_string_conversion
# ---------------------------------------------------------------------------


class TestSafeNpStringConversion:
    """Tests for :func:`safe_np_string_conversion`."""

    def test_object_to_unicode(self):
        """Object array of strings converts to fixed-width Unicode."""
        arr = np.array(["hello", "world"], dtype=object)
        result = safe_np_string_conversion(arr, dtype="<U75")
        assert result.dtype.kind == "U"
        np.testing.assert_array_equal(result, ["hello", "world"])

    def test_unicode_to_object(self):
        """Fixed-width Unicode converts to object dtype."""
        arr = np.array(["hello", "world"], dtype="<U75")
        result = safe_np_string_conversion(arr, dtype=object)
        assert result.dtype == np.dtype(object)

    def test_int_to_string_raises(self):
        """Non-string to string conversion raises TypeError."""
        arr = np.array([1, 2, 3], dtype=np.int32)
        with pytest.raises(TypeError):
            safe_np_string_conversion(arr, dtype="<U10")

    def test_object_to_int_raises(self):
        """String to non-string conversion raises TypeError."""
        arr = np.array(["hello"], dtype=object)
        with pytest.raises(TypeError):
            safe_np_string_conversion(arr, dtype=np.int32)

    def test_empty_object_array(self):
        """Empty object array converts successfully."""
        arr = np.array([], dtype=object)
        result = safe_np_string_conversion(arr, dtype="<U10")
        assert result.size == 0


# ---------------------------------------------------------------------------
# replay_event_log_to_snapshot
# ---------------------------------------------------------------------------


class TestReplayEventLogToSnapshot:
    """Tests for :func:`replay_event_log_to_snapshot`."""

    @pytest.fixture
    def event_dtype(self):
        return np.dtype([
            ("time", np.float64),
            ("entity", np.int32),
            ("value", np.int32),
        ])

    def _make_mapper(self):
        """Simple state_mapper that returns the 'value' field."""
        return lambda events: events["value"]

    def test_empty_log(self, event_dtype):
        """Empty event log produces all initial_value."""
        log = np.array([], dtype=event_dtype)
        snapshots = np.array([1.0, 2.0])
        result = replay_event_log_to_snapshot(
            log, snapshots, n_entities=3,
            entity_field="entity", time_field="time",
            state_mapper=self._make_mapper(),
            initial_value=0,
        )
        assert result.shape == (2, 3)
        np.testing.assert_array_equal(result, 0)

    def test_single_event(self, event_dtype):
        """Single event before second snapshot updates entity state."""
        log = np.array([(0.5, 2, 1)], dtype=event_dtype)
        snapshots = np.array([0.0, 1.0])
        result = replay_event_log_to_snapshot(
            log, snapshots, n_entities=4,
            entity_field="entity", time_field="time",
            state_mapper=self._make_mapper(),
            initial_value=0,
        )
        # At t=0.0: no events yet → all 0
        np.testing.assert_array_equal(result[0], [0, 0, 0, 0])
        # At t=1.0: event at t=0.5 set entity 2 to 1
        np.testing.assert_array_equal(result[1], [0, 0, 1, 0])

    def test_multiple_events(self, event_dtype):
        """Multiple events at different times produce correct per-snapshot state."""
        log = np.array([
            (0.5, 0, 10),
            (1.5, 1, 20),
            (2.5, 2, 30),
        ], dtype=event_dtype)
        snapshots = np.array([1.0, 2.0, 3.0])
        result = replay_event_log_to_snapshot(
            log, snapshots, n_entities=3,
            entity_field="entity", time_field="time",
            state_mapper=self._make_mapper(),
            initial_value=0,
        )
        # t=1.0: only event at t=0.5 applied
        np.testing.assert_array_equal(result[0], [10, 0, 0])
        # t=2.0: events at t=0.5 and t=1.5 applied
        np.testing.assert_array_equal(result[1], [10, 20, 0])
        # t=3.0: all events applied
        np.testing.assert_array_equal(result[2], [10, 20, 30])

    def test_overwrite_last_wins(self, event_dtype):
        """Two events on same entity: last write wins."""
        log = np.array([
            (0.5, 1, 10),
            (0.8, 1, 99),
        ], dtype=event_dtype)
        snapshots = np.array([2.0])
        result = replay_event_log_to_snapshot(
            log, snapshots, n_entities=3,
            entity_field="entity", time_field="time",
            state_mapper=self._make_mapper(),
            initial_value=0,
        )
        np.testing.assert_array_equal(result[0], [0, 99, 0])

    def test_no_events_before_first_snapshot(self, event_dtype):
        """Events after all snapshots leave output as initial_value."""
        log = np.array([(5.0, 0, 1)], dtype=event_dtype)
        snapshots = np.array([1.0, 2.0])
        result = replay_event_log_to_snapshot(
            log, snapshots, n_entities=2,
            entity_field="entity", time_field="time",
            state_mapper=self._make_mapper(),
            initial_value=0,
        )
        np.testing.assert_array_equal(result, 0)


# ---------------------------------------------------------------------------
# convert_structured_grid_fields
# ---------------------------------------------------------------------------


class TestConvertStructuredGridFields:
    """Tests for :func:`convert_structured_grid_fields`."""

    @pytest.fixture
    def grid_params(self):
        """Params for a 5x3 grid."""
        return {
            "micro_params": {"micro_simulations": 100},
            "macro_params": {
                "rows": 5,
                "cols": 3,
                "total_edges": 57,
                "total_molecules": 10,
            },
        }

    def _v199_dtype(self):
        """v1.99.0 dtype for f_deg_list."""
        return dataspec["v1.99.0"]["macroscale_out"].data["f_deg_list"].dtype

    def _v200_dtype(self):
        """v2.0.0 dtype for fiber_degrade_time."""
        return dataspec["v2.0.0"]["macroscale_out"].data["fiber_degrade_time"].dtype

    def test_1d_to_2d_conversion(self, grid_params):
        """v1.99.0 1D Fortran index converts to v2.0.0 2D row/rank."""
        v199_dt = self._v199_dtype()
        # Create a record with a known Fortran index
        # Fortran index 1 (1-based) → 0-based index 0 → row=0, rank=0
        table = np.array([(1.0, 1, 2.5)], dtype=v199_dt)
        input_data = {
            "params": grid_params,
            "f_deg_list": [table],
        }
        result = convert_structured_grid_fields(
            input_data,
            input_dataset="f_deg_list",
            output_dataset="fiber_degrade_time",
            input_spec="v1.99.0",
            output_spec="v2.0.0",
        )
        assert len(result) == 1
        out = result[0]
        # Common field preserved
        assert out["Simulation Time Elapsed"][0] == 1.0
        assert out["Fiber New Degrade Time"][0] == 2.5
        # Grid index was converted from 1-based 1D to 0-based 2D
        expected_coords = from_fortran_edge_index_array(
            np.array([0]), 5, 3
        )
        assert out["Grid Location Row"][0] == expected_coords[0, 0]
        assert out["Grid Location Rank"][0] == expected_coords[0, 1]

    def test_2d_to_1d_conversion(self, grid_params):
        """v2.0.0 2D row/rank converts to v1.99.0 1D Fortran index."""
        v200_dt = self._v200_dtype()
        # row=0, rank=0 → Fortran index 0 + 1 = 1
        table = np.array([(1.0, 0, 0, 2.5)], dtype=v200_dt)
        input_data = {
            "params": grid_params,
            "fiber_degrade_time": [table],
        }
        result = convert_structured_grid_fields(
            input_data,
            input_dataset="fiber_degrade_time",
            output_dataset="f_deg_list",
            input_spec="v2.0.0",
            output_spec="v1.99.0",
        )
        out = result[0]
        expected_fortran = to_fortran_edge_index_array(
            np.array([[0, 0]], dtype=np.uint32), 5, 3
        )
        assert out["Grid Location Index"][0] == expected_fortran[0] + 1

    def test_field_offset_applied(self, grid_params):
        """field_offsets increments tPA Molecule Index during conversion."""
        v200_bind_dt = dataspec["v2.0.0"]["macroscale_out"].data["tpa_bind_events"].dtype
        # tPA Molecule Index = 5 (0-based in v2.0.0)
        table = np.array([(1.0, 5, 0, 0, 0)], dtype=v200_bind_dt)
        input_data = {
            "params": grid_params,
            "tpa_bind_events": [table],
        }
        result = convert_structured_grid_fields(
            input_data,
            input_dataset="tpa_bind_events",
            output_dataset="m_bind_t",
            input_spec="v2.0.0",
            output_spec="v1.99.0",
            field_offsets={"tPA Molecule Index": 1},
        )
        # 5 + 1 = 6 (1-based in v1.99.0)
        assert result[0]["tPA Molecule Index"][0] == 6

    def test_common_fields_preserved(self, grid_params):
        """Non-grid common fields are copied without modification."""
        v199_dt = self._v199_dtype()
        table = np.array([(42.0, 1, 99.5)], dtype=v199_dt)
        input_data = {
            "params": grid_params,
            "f_deg_list": [table],
        }
        result = convert_structured_grid_fields(
            input_data,
            input_dataset="f_deg_list",
            output_dataset="fiber_degrade_time",
            input_spec="v1.99.0",
            output_spec="v2.0.0",
        )
        assert result[0]["Simulation Time Elapsed"][0] == 42.0
        assert result[0]["Fiber New Degrade Time"][0] == 99.5

    def test_round_trip(self, grid_params):
        """v1.99.0 → v2.0.0 → v1.99.0 recovers original data."""
        v199_dt = self._v199_dtype()
        original = np.array([(1.0, 5, 2.5), (3.0, 10, 4.5)], dtype=v199_dt)
        input_data = {
            "params": grid_params,
            "f_deg_list": [original],
        }
        # Forward: v1.99.0 → v2.0.0
        v200_result = convert_structured_grid_fields(
            input_data,
            input_dataset="f_deg_list",
            output_dataset="fiber_degrade_time",
            input_spec="v1.99.0",
            output_spec="v2.0.0",
        )
        # Reverse: v2.0.0 → v1.99.0
        mid_data = {
            "params": grid_params,
            "fiber_degrade_time": v200_result,
        }
        v199_result = convert_structured_grid_fields(
            mid_data,
            input_dataset="fiber_degrade_time",
            output_dataset="f_deg_list",
            input_spec="v2.0.0",
            output_spec="v1.99.0",
        )
        out = v199_result[0]
        np.testing.assert_array_almost_equal(
            out["Simulation Time Elapsed"], original["Simulation Time Elapsed"]
        )
        np.testing.assert_array_equal(
            out["Grid Location Index"], original["Grid Location Index"]
        )
        np.testing.assert_array_almost_equal(
            out["Fiber New Degrade Time"], original["Fiber New Degrade Time"]
        )


# ---------------------------------------------------------------------------
# convert_location_snapshot
# ---------------------------------------------------------------------------


class TestConvertLocationSnapshot:
    """Tests for :func:`convert_location_snapshot`."""

    @pytest.fixture
    def grid_params(self):
        """Params for a 5x3 grid with 10 molecules."""
        return {
            "micro_params": {"micro_simulations": 100},
            "macro_params": {
                "rows": 5,
                "cols": 3,
                "total_edges": 57,
                "total_molecules": 10,
            },
        }

    def _make_v199_m_loc(self, rows, cols, n_mol, n_snapshots):
        """Create a synthetic v1.99.0 m_loc array with valid 1-based Fortran indices."""
        full_row = 3 * cols - 1
        xz_row = 2 * cols - 1
        total_edges = full_row * (rows - 1) + xz_row
        rng = np.random.RandomState(42)
        # 1-based Fortran indices
        return rng.randint(1, total_edges + 1, size=(n_snapshots, n_mol)).astype(np.int32)

    def test_v199_to_v200_shape(self, grid_params):
        """v1.99.0 (3,10) → v2.0.0 (10,2,3) shape transform."""
        m_loc = self._make_v199_m_loc(5, 3, 10, 3)
        input_data = {"params": grid_params, "m_loc": [m_loc]}
        result = convert_location_snapshot(
            input_data, input_dataset="m_loc",
            input_spec="v1.99.0", output_spec="v2.0.0",
        )
        assert result[0].shape == (10, 2, 3)

    def test_v200_to_v199_shape(self, grid_params):
        """v2.0.0 (10,2,3) → v1.99.0 (3,10) shape transform."""
        # Create v2.0.0 format: (n_mol, 2, n_snapshots) 0-based coords
        rng = np.random.RandomState(42)
        # Generate valid 2D coords — top row (row 4) has no y-edges
        # (rank % 3 == 0), so restrict to non-y ranks there.
        coords_2d = np.zeros((10, 2, 3), dtype=np.int32)
        for i in range(10):
            for j in range(3):
                row = rng.randint(0, 5)
                if row == 4:
                    # Top row: only x-edges (rank%3==1) and z-edges (rank%3==2)
                    # Valid ranks: 1, 2, 4, 5, 7
                    valid_ranks = [1, 2, 4, 5, 7]
                    rank = valid_ranks[rng.randint(0, len(valid_ranks))]
                else:
                    rank = rng.randint(0, 8)
                coords_2d[i, 0, j] = row
                coords_2d[i, 1, j] = rank
        input_data = {"params": grid_params, "tpa_location_snapshot": [coords_2d]}
        result = convert_location_snapshot(
            input_data, input_dataset="tpa_location_snapshot",
            input_spec="v2.0.0", output_spec="v1.99.0",
        )
        assert result[0].shape == (3, 10)

    def test_v199_to_v200_index_conversion(self, grid_params):
        """v1.99.0 1-based indices convert to v2.0.0 0-based 2D coords."""
        m_loc = self._make_v199_m_loc(5, 3, 10, 3)
        input_data = {"params": grid_params, "m_loc": [m_loc]}
        result = convert_location_snapshot(
            input_data, input_dataset="m_loc",
            input_spec="v1.99.0", output_spec="v2.0.0",
        )
        # Verify specific molecule at specific snapshot
        out = result[0]
        for snap_idx in range(3):
            for mol_idx in range(10):
                fortran_1based = m_loc[snap_idx, mol_idx]
                expected_coords = from_fortran_edge_index_array(
                    np.array([fortran_1based - 1]), 5, 3
                )
                assert out[mol_idx, 0, snap_idx] == expected_coords[0, 0]
                assert out[mol_idx, 1, snap_idx] == expected_coords[0, 1]

    def test_round_trip(self, grid_params):
        """v1.99.0 → v2.0.0 → v1.99.0 recovers original data."""
        m_loc = self._make_v199_m_loc(5, 3, 10, 3)
        input_data = {"params": grid_params, "m_loc": [m_loc]}
        # Forward
        v200 = convert_location_snapshot(
            input_data, input_dataset="m_loc",
            input_spec="v1.99.0", output_spec="v2.0.0",
        )
        # Reverse
        mid_data = {"params": grid_params, "tpa_location_snapshot": v200}
        v199 = convert_location_snapshot(
            mid_data, input_dataset="tpa_location_snapshot",
            input_spec="v2.0.0", output_spec="v1.99.0",
        )
        np.testing.assert_array_equal(v199[0], m_loc)

    def test_unsupported_spec_raises(self, grid_params):
        """Unsupported spec version raises NotImplementedError."""
        m_loc = self._make_v199_m_loc(5, 3, 10, 3)
        input_data = {"params": grid_params, "m_loc": [m_loc]}
        with pytest.raises(NotImplementedError):
            convert_location_snapshot(
                input_data, input_dataset="m_loc",
                input_spec="v3.0.0", output_spec="v2.0.0",
            )


# ---------------------------------------------------------------------------
# convert_bind_events_to_bound
# ---------------------------------------------------------------------------


class TestConvertBindEventsToBound:
    """Tests for :func:`convert_bind_events_to_bound`."""

    @pytest.fixture
    def bind_event_dtype(self):
        """dtype matching v2.0.0 tpa_bind_events."""
        return dataspec["v2.0.0"]["macroscale_out"].data["tpa_bind_events"].dtype

    @pytest.fixture
    def base_params(self):
        return {
            "micro_params": {},
            "macro_params": {"total_molecules": 5},
        }

    def test_no_events_all_zeros(self, bind_event_dtype, base_params):
        """Empty event log produces all-zero bound status."""
        events = np.array([], dtype=bind_event_dtype)
        snap_times = np.array([1.0, 2.0, 3.0])
        input_data = {
            "params": base_params,
            "tpa_bind_events": [events],
            "snapshot_time": [snap_times],
        }
        result = convert_bind_events_to_bound(input_data)
        assert len(result) == 1
        assert result[0].shape == (3, 5)
        np.testing.assert_array_equal(result[0], 0)

    def test_bind_then_unbind(self, bind_event_dtype, base_params):
        """Bind at t=1, unbind at t=3 → bound between, unbound after."""
        events = np.array([
            (1.0, 2, CONST.MOL_STATUS.BOUND, 0, 0),
            (3.0, 2, CONST.MOL_STATUS.UNBOUND, 0, 0),
        ], dtype=bind_event_dtype)
        snap_times = np.array([0.5, 2.0, 4.0])
        input_data = {
            "params": base_params,
            "tpa_bind_events": [events],
            "snapshot_time": [snap_times],
        }
        result = convert_bind_events_to_bound(input_data)
        m_bound = result[0]
        # t=0.5: no events yet
        assert m_bound[0, 2] == 0
        # t=2.0: bound at t=1
        assert m_bound[1, 2] == 1
        # t=4.0: unbound at t=3
        assert m_bound[2, 2] == 0

    def test_multiple_molecules(self, bind_event_dtype, base_params):
        """Different molecules tracked independently."""
        events = np.array([
            (0.5, 0, CONST.MOL_STATUS.BOUND, 0, 0),
            (1.5, 3, CONST.MOL_STATUS.BOUND, 0, 0),
        ], dtype=bind_event_dtype)
        snap_times = np.array([1.0, 2.0])
        input_data = {
            "params": base_params,
            "tpa_bind_events": [events],
            "snapshot_time": [snap_times],
        }
        result = convert_bind_events_to_bound(input_data)
        m_bound = result[0]
        # t=1.0: mol 0 bound, mol 3 not yet
        assert m_bound[0, 0] == 1
        assert m_bound[0, 3] == 0
        # t=2.0: both bound
        assert m_bound[1, 0] == 1
        assert m_bound[1, 3] == 1


# ---------------------------------------------------------------------------
# convert_data (main entry point)
# ---------------------------------------------------------------------------


class TestConvertData:
    """Tests for :func:`convert_data`."""

    def test_params_preserved(self):
        """Parameters are copied directly without conversion."""
        params = {
            "micro_params": {"micro_simulations": 100},
            "macro_params": {"rows": 5, "cols": 3, "total_edges": 57, "total_molecules": 10},
        }
        # Create minimal microscale data for v1.99.0
        n = 100
        input_data = {
            "params": params,
            "micro_log": np.array(["sim1", "sim2"], dtype="<U75"),
            "firstPLi": np.array([1.0] * n, dtype=np.float64),
            "lasttPA": np.array([0] * n, dtype=np.int32),
            "lyscomplete": np.array([1] * n, dtype=np.int32),
            "lysis": np.array([10.0] * n, dtype=np.float64),
            "PLi": np.array([5] * n, dtype=np.int32),
            "tPA_time": np.array([2.0] * n, dtype=np.float64),
            "tPAPLiunbd": np.array([0] * n, dtype=np.int32),
            "tPAunbind": np.array([1] * n, dtype=np.int32),
        }
        result = convert_data(input_data, "v1.99.0", "v2.0.0")
        assert result["params"] is input_data["params"]

    def test_tag_resolution(self):
        """Tag aliases are resolved before conversion."""
        params = {
            "micro_params": {"micro_simulations": 100},
            "macro_params": {"rows": 5, "cols": 3, "total_edges": 57, "total_molecules": 10},
        }
        n = 100
        input_data = {
            "params": params,
            "micro_log": np.array(["sim1"], dtype="<U75"),
            "firstPLi": np.array([1.0] * n, dtype=np.float64),
            "lasttPA": np.array([0] * n, dtype=np.int32),
            "lyscomplete": np.array([1] * n, dtype=np.int32),
            "lysis": np.array([10.0] * n, dtype=np.float64),
            "PLi": np.array([5] * n, dtype=np.int32),
            "tPA_time": np.array([2.0] * n, dtype=np.float64),
            "tPAPLiunbd": np.array([0] * n, dtype=np.int32),
            "tPAunbind": np.array([1] * n, dtype=np.int32),
        }
        # "fortran" should resolve to "v1.99.0", "hdf5" to "v2.0.0"
        result = convert_data(input_data, "fortran", "hdf5")
        # If it didn't raise, tags resolved correctly
        assert "pli_first_time" in result  # v2.0.0 dataset name

    def test_partial_data_skips_missing_collections(self):
        """Collections not in input data are skipped without error."""
        params = {
            "micro_params": {"micro_simulations": 100},
            "macro_params": {"rows": 5, "cols": 3, "total_edges": 57, "total_molecules": 10},
        }
        n = 100
        # Only microscale data, no macroscale
        input_data = {
            "params": params,
            "micro_log": np.array(["sim1"], dtype="<U75"),
            "firstPLi": np.array([1.0] * n, dtype=np.float64),
            "lasttPA": np.array([0] * n, dtype=np.int32),
            "lyscomplete": np.array([1] * n, dtype=np.int32),
            "lysis": np.array([10.0] * n, dtype=np.float64),
            "PLi": np.array([5] * n, dtype=np.int32),
            "tPA_time": np.array([2.0] * n, dtype=np.float64),
            "tPAPLiunbd": np.array([0] * n, dtype=np.int32),
            "tPAunbind": np.array([1] * n, dtype=np.int32),
        }
        # Should not raise even though macroscale data is missing
        result = convert_data(input_data, "v1.99.0", "v2.0.0")
        assert "pli_first_time" in result
        # Macroscale datasets should not be in the output
        assert "snapshot_time" not in result
