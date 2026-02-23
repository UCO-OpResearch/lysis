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
* _build_conversion_graph — undirected graph from data_converters
* _build_spanning_tree — BFS spanning tree
* _extract_tree_path — LCA-based path extraction
* _build_conversion_paths — orchestrator for precomputed paths
* conversion_paths — module-level precomputed paths
* _convert_single_step — single-hop conversion
* convert_data multi-step — chaining, short-circuit, ValueError
"""

import warnings

import numpy as np
import pytest

from lysis.config.constants import CONST, Q_
from lysis.config.parameters import Parameters
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
    conversion_paths,
    data_converters,
    params_converters,
    _build_conversion_graph,
    _build_spanning_tree,
    _extract_tree_path,
    _build_conversion_paths,
    _convert_single_step,
    _convert_params_add_units,
    _convert_params_strip_units,
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
    """Tests for :func:`convert_data`.

    Synthetic v1.99.0 data: dimensioned params are stored as unit strings
    (e.g., ``"0.000534 centimeter"``); unitless params are plain types.
    """

    @pytest.fixture
    def v199_microscale_data(self):
        """Minimal v1.99.0 microscale data with unit-string params."""
        n = 100
        return {
            "params": {
                "micro_params": {
                    "pore_size": "0.000534 centimeter",
                    "total_time": "3600.0 second",
                    "micro_simulations": 100,
                },
                "macro_params": {
                    "diffusion_coeff": "5e-07 centimeter ** 2 / second",
                    "rows": 5,
                    "cols": 3,
                    "total_edges": 57,
                    "total_molecules": 10,
                },
            },
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

    def test_params_preserved(self, v199_microscale_data):
        """v1.99.0 → v2.0.0: params (including unit strings) are copied by reference."""
        result = convert_data(v199_microscale_data, "v1.99.0", "v2.0.0")
        assert result["params"] is v199_microscale_data["params"]
        # Unit strings pass through unchanged
        assert result["params"]["micro_params"]["pore_size"] == "0.000534 centimeter"

    def test_tag_resolution(self, v199_microscale_data):
        """Tag aliases are resolved before conversion."""
        # "fortran" should resolve to "v1.99.0", "hdf5" to "v2.0.0"
        result = convert_data(v199_microscale_data, "fortran", "hdf5")
        # If it didn't raise, tags resolved correctly
        assert "pli_first_time" in result  # v2.0.0 dataset name

    def test_partial_data_skips_missing_collections(self, v199_microscale_data):
        """Collections not in input data are skipped without error."""
        # Should not raise even though macroscale data is missing
        result = convert_data(v199_microscale_data, "v1.99.0", "v2.0.0")
        assert "pli_first_time" in result
        # Macroscale datasets should not be in the output
        assert "snapshot_time" not in result

    def test_missing_params_raises(self):
        """Data without parameters raises ValueError."""
        input_data = {
            "micro_log": np.array(["sim1"], dtype="<U75"),
        }
        with pytest.raises(ValueError, match="parameters"):
            convert_data(input_data, "v1.99.0", "v2.0.0")

    def test_empty_params_raises(self):
        """Data with empty parameters dict raises ValueError."""
        input_data = {
            "params": {},
            "micro_log": np.array(["sim1"], dtype="<U75"),
        }
        with pytest.raises(ValueError, match="parameters"):
            convert_data(input_data, "v1.99.0", "v2.0.0")


# ---------------------------------------------------------------------------
# _convert_params_add_units  (v1.95.0 → v1.99.0)
# ---------------------------------------------------------------------------


class TestConvertParamsAddUnits:
    """Tests for :func:`_convert_params_add_units`."""

    def test_bare_magnitude_gets_units(self):
        """Numeric value for a known unit key becomes a Quantity string."""
        params = {
            "micro_params": {"pore_size": 0.000534},
        }
        result = _convert_params_add_units(params)
        assert isinstance(result["micro_params"]["pore_size"], str)
        q = Q_(result["micro_params"]["pore_size"])
        assert q.magnitude == pytest.approx(0.000534)
        assert str(q.units) == "centimeter"

    def test_unitless_int_passes_through(self):
        """Integer value for a key NOT in Parameters.units() is unchanged."""
        params = {
            "macro_params": {"rows": 19, "cols": 3},
        }
        result = _convert_params_add_units(params)
        assert result["macro_params"]["rows"] == 19
        assert result["macro_params"]["cols"] == 3

    def test_non_dict_section_passes_through(self):
        """Non-dict sections (None, strings, lists) are copied verbatim."""
        params = {
            "micro_params": {"pore_size": 0.000534},
            "metadata": "some string",
            "extra_list": [1, 2, 3],
        }
        result = _convert_params_add_units(params)
        assert result["metadata"] == "some string"
        assert result["extra_list"] == [1, 2, 3]

    def test_string_value_for_unit_key_passes_through(self):
        """A string value (already has units) is not double-wrapped."""
        params = {
            "micro_params": {"pore_size": "0.000534 centimeter"},
        }
        result = _convert_params_add_units(params)
        assert result["micro_params"]["pore_size"] == "0.000534 centimeter"

    def test_unknown_key_passes_through(self):
        """Keys not in Parameters.units() are copied as-is."""
        params = {
            "micro_params": {"custom_key": 42.0},
        }
        result = _convert_params_add_units(params)
        assert result["micro_params"]["custom_key"] == 42.0

    def test_multiple_sections(self):
        """Both micro_params and macro_params are converted."""
        units = Parameters.units()
        params = {
            "micro_params": {"pore_size": 0.000534, "micro_simulations": 100},
            "macro_params": {"total_time": 3600.0, "rows": 19},
        }
        result = _convert_params_add_units(params)
        # Dimensioned values become strings
        assert isinstance(result["micro_params"]["pore_size"], str)
        assert isinstance(result["macro_params"]["total_time"], str)
        # Unitless values unchanged
        assert result["micro_params"]["micro_simulations"] == 100
        assert result["macro_params"]["rows"] == 19


# ---------------------------------------------------------------------------
# _convert_params_strip_units  (v1.99.0 → v1.95.0)
# ---------------------------------------------------------------------------


class TestConvertParamsStripUnits:
    """Tests for :func:`_convert_params_strip_units`."""

    def test_unit_string_becomes_magnitude(self):
        """String value for a known unit key becomes a bare magnitude."""
        params = {
            "micro_params": {"pore_size": "0.000534 centimeter"},
        }
        result = _convert_params_strip_units(params)
        assert isinstance(result["micro_params"]["pore_size"], float)
        assert result["micro_params"]["pore_size"] == pytest.approx(0.000534)

    def test_unit_conversion_to_canonical(self):
        """Value in compatible units is converted to the canonical unit."""
        # pore_size canonical unit is "centimeters"
        params = {
            "micro_params": {"pore_size": "5.34 micrometer"},
        }
        result = _convert_params_strip_units(params)
        expected = Q_("5.34 micrometer").to("centimeters").magnitude
        assert result["micro_params"]["pore_size"] == pytest.approx(expected)

    def test_unitless_int_passes_through(self):
        """Integer value for a key NOT in Parameters.units() is unchanged."""
        params = {
            "macro_params": {"rows": 19},
        }
        result = _convert_params_strip_units(params)
        assert result["macro_params"]["rows"] == 19

    def test_non_dict_section_passes_through(self):
        """Non-dict sections are copied verbatim."""
        params = {
            "micro_params": {"pore_size": "0.000534 centimeter"},
            "metadata": "some string",
        }
        result = _convert_params_strip_units(params)
        assert result["metadata"] == "some string"

    def test_numeric_value_for_unit_key_passes_through(self):
        """A numeric value (already bare) is not re-stripped."""
        params = {
            "micro_params": {"pore_size": 0.000534},
        }
        result = _convert_params_strip_units(params)
        assert result["micro_params"]["pore_size"] == pytest.approx(0.000534)

    def test_string_value_for_non_unit_key_passes_through(self):
        """A string value for a key not in Parameters.units() passes through."""
        params = {
            "micro_params": {"micro_version": "micro_rates"},
        }
        result = _convert_params_strip_units(params)
        assert result["micro_params"]["micro_version"] == "micro_rates"


# ---------------------------------------------------------------------------
# params_converters registry
# ---------------------------------------------------------------------------


class TestParamsConvertersRegistry:
    """Tests for the :data:`params_converters` registry."""

    def test_both_directions_registered(self):
        """Both v1.95.0↔v1.99.0 entries exist."""
        assert ("v1.95.0", "v1.99.0") in params_converters
        assert ("v1.99.0", "v1.95.0") in params_converters

    def test_functions_match(self):
        """Registry entries point to the correct functions."""
        assert params_converters["v1.95.0", "v1.99.0"] is _convert_params_add_units
        assert params_converters["v1.99.0", "v1.95.0"] is _convert_params_strip_units


# ---------------------------------------------------------------------------
# Round-trip params conversion
# ---------------------------------------------------------------------------


class TestParamsRoundTrip:
    """Round-trip tests for params conversion."""

    def test_add_then_strip_recovers_magnitudes(self):
        """v1.95.0 → v1.99.0 → v1.95.0 recovers original bare magnitudes."""
        original = {
            "micro_params": {
                "pore_size": 0.000534,
                "total_time": 3600.0,
                "micro_simulations": 100,
                "micro_version": "micro_rates",
            },
            "macro_params": {
                "rows": 19,
                "cols": 3,
                "diffusion_coeff": 5e-7,
            },
        }
        with_units = _convert_params_add_units(original)
        recovered = _convert_params_strip_units(with_units)

        # Dimensioned values recovered
        assert recovered["micro_params"]["pore_size"] == pytest.approx(0.000534)
        assert recovered["micro_params"]["total_time"] == pytest.approx(3600.0)
        assert recovered["macro_params"]["diffusion_coeff"] == pytest.approx(5e-7)

        # Unitless values unchanged
        assert recovered["micro_params"]["micro_simulations"] == 100
        assert recovered["micro_params"]["micro_version"] == "micro_rates"
        assert recovered["macro_params"]["rows"] == 19
        assert recovered["macro_params"]["cols"] == 3

    def test_strip_then_add_recovers_strings(self):
        """v1.99.0 → v1.95.0 → v1.99.0 recovers original unit strings."""
        original = {
            "micro_params": {
                "pore_size": "0.000534 centimeter",
                "total_time": "3600.0 second",
                "micro_simulations": 100,
            },
        }
        stripped = _convert_params_strip_units(original)
        recovered = _convert_params_add_units(stripped)

        # Parse both to Quantity and compare magnitudes and units
        orig_ps = Q_(original["micro_params"]["pore_size"])
        recov_ps = Q_(recovered["micro_params"]["pore_size"])
        assert recov_ps.magnitude == pytest.approx(orig_ps.magnitude)
        assert recov_ps.units == orig_ps.units

        assert recovered["micro_params"]["micro_simulations"] == 100


# ---------------------------------------------------------------------------
# convert_data with v1.95.0 spec
# ---------------------------------------------------------------------------


class TestConvertDataV195:
    """Tests for :func:`convert_data` with v1.95.0 ↔ v1.99.0 conversion."""

    @pytest.fixture
    def v195_microscale_data(self):
        """Minimal v1.95.0 microscale data with bare-magnitude params."""
        n = 100
        return {
            "params": {
                "micro_params": {
                    "pore_size": 0.000534,
                    "total_time": 3600.0,
                    "micro_simulations": 100,
                },
                "macro_params": {
                    "rows": 5,
                    "cols": 3,
                    "total_edges": 57,
                    "total_molecules": 10,
                    "diffusion_coeff": 5e-7,
                },
            },
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

    def test_v195_to_v199_params_converted(self, v195_microscale_data):
        """v1.95.0 → v1.99.0 wraps bare magnitudes as unit strings."""
        result = convert_data(v195_microscale_data, "v1.95.0", "v1.99.0")
        # Params should be a new dict (not same object)
        assert result["params"] is not v195_microscale_data["params"]
        # Dimensioned param now has units
        assert isinstance(result["params"]["micro_params"]["pore_size"], str)
        q = Q_(result["params"]["micro_params"]["pore_size"])
        assert q.magnitude == pytest.approx(0.000534)
        # Unitless param unchanged
        assert result["params"]["micro_params"]["micro_simulations"] == 100

    def test_v195_to_v199_data_identity(self, v195_microscale_data):
        """v1.95.0 → v1.99.0 data datasets are passed through unchanged."""
        result = convert_data(v195_microscale_data, "v1.95.0", "v1.99.0")
        np.testing.assert_array_equal(
            result["firstPLi"], v195_microscale_data["firstPLi"]
        )
        np.testing.assert_array_equal(
            result["lysis"], v195_microscale_data["lysis"]
        )

    def test_v199_to_v195_params_converted(self):
        """v1.99.0 → v1.95.0 strips unit strings to bare magnitudes."""
        n = 100
        input_data = {
            "params": {
                "micro_params": {
                    "pore_size": "0.000534 centimeter",
                    "micro_simulations": 100,
                },
                "macro_params": {
                    "rows": 5,
                    "cols": 3,
                    "total_edges": 57,
                    "total_molecules": 10,
                },
            },
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
        result = convert_data(input_data, "v1.99.0", "v1.95.0")
        assert isinstance(result["params"]["micro_params"]["pore_size"], float)
        assert result["params"]["micro_params"]["pore_size"] == pytest.approx(0.000534)

    def test_v195_to_v199_to_v200_chained(self, v195_microscale_data):
        """v1.95.0 → v1.99.0 → v2.0.0 chained conversion produces valid output."""
        intermediate = convert_data(v195_microscale_data, "v1.95.0", "v1.99.0")
        result = convert_data(intermediate, "v1.99.0", "v2.0.0")
        # v2.0.0 dataset names should be present
        assert "pli_first_time" in result
        assert "sim_final_time" in result
        # Params carried through
        assert result["params"]["micro_params"]["micro_simulations"] == 100

    def test_partial_data_skips_missing_collections(self, v195_microscale_data):
        """Collections not in input data are skipped without error."""
        result = convert_data(v195_microscale_data, "v1.95.0", "v1.99.0")
        assert "firstPLi" in result
        # macroscale datasets should not be in the output
        assert "tsave" not in result


# ---------------------------------------------------------------------------
# _build_conversion_graph
# ---------------------------------------------------------------------------


class TestBuildConversionGraph:
    """Tests for :func:`_build_conversion_graph`."""

    def test_returns_undirected_graph(self):
        """Graph has symmetric edges for all bidirectional converter pairs."""
        graph = _build_conversion_graph()
        for node, neighbors in graph.items():
            for neighbor in neighbors:
                assert node in graph[neighbor], (
                    f"{node!r} has neighbor {neighbor!r} but not vice-versa"
                )

    def test_all_converter_versions_present(self):
        """Every version mentioned in data_converters appears in the graph."""
        graph = _build_conversion_graph()
        for a, b in data_converters:
            assert a in graph
            assert b in graph

    def test_bidirectional_edges_included(self):
        """Pairs that exist in both directions appear as edges."""
        graph = _build_conversion_graph()
        # v1.99.0 ↔ v2.0.0 is bidirectional
        assert "v2.0.0" in graph["v1.99.0"]
        assert "v1.99.0" in graph["v2.0.0"]

    def test_known_edges(self):
        """Known bidirectional pairs are present in the graph."""
        graph = _build_conversion_graph()
        # v1.95.0 ↔ v1.99.0
        assert "v1.99.0" in graph["v1.95.0"]
        assert "v1.95.0" in graph["v1.99.0"]
        # v1.99.0 ↔ v2.0.0
        assert "v2.0.0" in graph["v1.99.0"]
        assert "v1.99.0" in graph["v2.0.0"]

    def test_one_way_edge_warns(self, monkeypatch):
        """A converter pair with only one direction emits a warning."""
        # Temporarily add a one-way edge
        fake_converters = dict(data_converters)
        fake_converters[("v99.0.0", "v2.0.0")] = {}
        monkeypatch.setattr(
            "lysis.data.dataconvert.data_converters", fake_converters
        )
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            graph = _build_conversion_graph()

        warning_messages = [str(x.message) for x in w]
        assert any("v99.0.0" in msg and "v2.0.0" in msg for msg in warning_messages)
        # The one-way edge should NOT appear in the undirected graph
        assert "v2.0.0" not in graph.get("v99.0.0", set())


# ---------------------------------------------------------------------------
# _build_spanning_tree
# ---------------------------------------------------------------------------


class TestBuildSpanningTree:
    """Tests for :func:`_build_spanning_tree`."""

    def test_root_parent_is_none(self):
        """Root node has None as parent."""
        graph = {"A": {"B"}, "B": {"A", "C"}, "C": {"B"}}
        parent = _build_spanning_tree(graph, "A")
        assert parent["A"] is None

    def test_all_reachable_nodes_included(self):
        """All reachable nodes appear in the parent dict."""
        graph = {"A": {"B"}, "B": {"A", "C"}, "C": {"B"}}
        parent = _build_spanning_tree(graph, "A")
        assert set(parent.keys()) == {"A", "B", "C"}

    def test_parent_pointers_form_valid_tree(self):
        """Every non-root node's parent is in the tree."""
        graph = {"A": {"B", "C"}, "B": {"A"}, "C": {"A", "D"}, "D": {"C"}}
        parent = _build_spanning_tree(graph, "A")
        for node, par in parent.items():
            if par is not None:
                assert par in parent

    def test_disconnected_component_excluded(self):
        """Nodes unreachable from root are not in the parent dict."""
        graph = {
            "A": {"B"},
            "B": {"A"},
            "X": {"Y"},
            "Y": {"X"},
        }
        parent = _build_spanning_tree(graph, "A")
        assert "A" in parent
        assert "B" in parent
        assert "X" not in parent
        assert "Y" not in parent

    def test_deterministic_with_sorted_neighbors(self):
        """Sorted neighbor iteration produces deterministic tree."""
        graph = {"A": {"C", "B"}, "B": {"A"}, "C": {"A"}}
        parent1 = _build_spanning_tree(graph, "A")
        parent2 = _build_spanning_tree(graph, "A")
        assert parent1 == parent2

    def test_single_node(self):
        """Single-node graph produces tree with just the root."""
        graph = {"A": set()}
        parent = _build_spanning_tree(graph, "A")
        assert parent == {"A": None}

    def test_linear_chain(self):
        """A → B → C → D produces linear parent chain."""
        graph = {
            "A": {"B"},
            "B": {"A", "C"},
            "C": {"B", "D"},
            "D": {"C"},
        }
        parent = _build_spanning_tree(graph, "A")
        assert parent["B"] == "A"
        assert parent["C"] == "B"
        assert parent["D"] == "C"


# ---------------------------------------------------------------------------
# _extract_tree_path
# ---------------------------------------------------------------------------


class TestExtractTreePath:
    """Tests for :func:`_extract_tree_path`."""

    @pytest.fixture
    def linear_tree(self):
        """Linear tree: A → B → C → D."""
        return {"A": None, "B": "A", "C": "B", "D": "C"}

    @pytest.fixture
    def branching_tree(self):
        """Branching tree:
              A
             / \\
            B   C
           /     \\
          D       E
        """
        return {"A": None, "B": "A", "C": "A", "D": "B", "E": "C"}

    def test_adjacent_nodes(self, linear_tree):
        """Path between adjacent nodes is just the two nodes."""
        path = _extract_tree_path(linear_tree, "A", "B")
        assert path == ["A", "B"]

    def test_same_direction_as_tree(self, linear_tree):
        """Path following tree direction: A → B → C."""
        path = _extract_tree_path(linear_tree, "A", "C")
        assert path == ["A", "B", "C"]

    def test_reverse_direction(self, linear_tree):
        """Path against tree direction: D → C → B → A."""
        path = _extract_tree_path(linear_tree, "D", "A")
        assert path == ["D", "C", "B", "A"]

    def test_middle_nodes(self, linear_tree):
        """Path between two non-root nodes: B → C → D."""
        path = _extract_tree_path(linear_tree, "B", "D")
        assert path == ["B", "C", "D"]

    def test_cross_branch_path(self, branching_tree):
        """Path across branches goes through LCA: D → B → A → C → E."""
        path = _extract_tree_path(branching_tree, "D", "E")
        assert path == ["D", "B", "A", "C", "E"]

    def test_node_to_root(self, branching_tree):
        """Path from leaf to root: E → C → A."""
        path = _extract_tree_path(branching_tree, "E", "A")
        assert path == ["E", "C", "A"]

    def test_root_to_leaf(self, branching_tree):
        """Path from root to leaf: A → C → E."""
        path = _extract_tree_path(branching_tree, "A", "E")
        assert path == ["A", "C", "E"]

    def test_sibling_path(self, branching_tree):
        """Path between siblings (same parent): B → A → C."""
        path = _extract_tree_path(branching_tree, "B", "C")
        assert path == ["B", "A", "C"]

    def test_path_start_equals_end(self, linear_tree):
        """Path from a node to itself is just that node."""
        path = _extract_tree_path(linear_tree, "B", "B")
        assert path == ["B"]


# ---------------------------------------------------------------------------
# _build_conversion_paths
# ---------------------------------------------------------------------------


class TestBuildConversionPaths:
    """Tests for :func:`_build_conversion_paths`."""

    def test_returns_dict(self):
        """Return type is a dict."""
        paths = _build_conversion_paths()
        assert isinstance(paths, dict)

    def test_all_pairs_present(self):
        """Every (src, dst) pair where src != dst is covered."""
        paths = _build_conversion_paths()
        graph = _build_conversion_graph()
        parent = _build_spanning_tree(graph, min(graph.keys()))
        reachable = set(parent.keys())
        for src in reachable:
            for dst in reachable:
                if src != dst:
                    assert (src, dst) in paths

    def test_no_self_loops(self):
        """No (x, x) entries exist."""
        paths = _build_conversion_paths()
        for src, dst in paths:
            assert src != dst

    def test_path_endpoints_match_key(self):
        """Each path starts with src and ends with dst."""
        paths = _build_conversion_paths()
        for (src, dst), path in paths.items():
            assert path[0] == src, f"Path for ({src}, {dst}) doesn't start with {src}"
            assert path[-1] == dst, f"Path for ({src}, {dst}) doesn't end with {dst}"

    def test_consecutive_pairs_have_converters(self):
        """Each consecutive pair in a path has a direct converter."""
        paths = _build_conversion_paths()
        for (src, dst), path in paths.items():
            for step_in, step_out in zip(path[:-1], path[1:]):
                assert (step_in, step_out) in data_converters, (
                    f"Path ({src} → {dst}): no converter for "
                    f"({step_in}, {step_out})"
                )

    def test_known_multi_step_path(self):
        """v1.95.0 → v2.0.0 routes through v1.99.0."""
        paths = _build_conversion_paths()
        assert paths[("v1.95.0", "v2.0.0")] == ["v1.95.0", "v1.99.0", "v2.0.0"]

    def test_known_reverse_multi_step_path(self):
        """v2.0.0 → v1.95.0 routes through v1.99.0."""
        paths = _build_conversion_paths()
        assert paths[("v2.0.0", "v1.95.0")] == ["v2.0.0", "v1.99.0", "v1.95.0"]

    def test_direct_paths_have_length_two(self):
        """Adjacent versions produce paths of length 2."""
        paths = _build_conversion_paths()
        assert len(paths[("v1.99.0", "v2.0.0")]) == 2
        assert len(paths[("v1.95.0", "v1.99.0")]) == 2

    def test_warns_on_disconnected_graph(self, monkeypatch):
        """Unreachable versions produce a warning."""
        # Add a disconnected island
        fake_converters = dict(data_converters)
        fake_converters[("v99.0.0", "v98.0.0")] = {}
        fake_converters[("v98.0.0", "v99.0.0")] = {}
        monkeypatch.setattr(
            "lysis.data.dataconvert.data_converters", fake_converters
        )
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            paths = _build_conversion_paths()

        warning_messages = [str(x.message) for x in w]
        assert any("not reachable" in msg for msg in warning_messages)

    def test_empty_converters(self, monkeypatch):
        """Empty data_converters produces empty paths."""
        monkeypatch.setattr("lysis.data.dataconvert.data_converters", {})
        paths = _build_conversion_paths()
        assert paths == {}


# ---------------------------------------------------------------------------
# conversion_paths (module-level precomputed)
# ---------------------------------------------------------------------------


class TestConversionPaths:
    """Tests for the module-level :data:`conversion_paths` variable."""

    def test_is_dict(self):
        """Module-level conversion_paths is a dict."""
        assert isinstance(conversion_paths, dict)

    def test_contains_direct_pairs(self):
        """Direct conversion pairs are present."""
        assert ("v1.99.0", "v2.0.0") in conversion_paths
        assert ("v2.0.0", "v1.99.0") in conversion_paths
        assert ("v1.95.0", "v1.99.0") in conversion_paths
        assert ("v1.99.0", "v1.95.0") in conversion_paths

    def test_contains_multi_step_pairs(self):
        """Multi-step pairs are present."""
        assert ("v1.95.0", "v2.0.0") in conversion_paths
        assert ("v2.0.0", "v1.95.0") in conversion_paths

    def test_v195_to_v200_path(self):
        """v1.95.0 → v2.0.0 path goes through v1.99.0."""
        assert conversion_paths[("v1.95.0", "v2.0.0")] == [
            "v1.95.0",
            "v1.99.0",
            "v2.0.0",
        ]

    def test_v200_to_v195_path(self):
        """v2.0.0 → v1.95.0 path goes through v1.99.0."""
        assert conversion_paths[("v2.0.0", "v1.95.0")] == [
            "v2.0.0",
            "v1.99.0",
            "v1.95.0",
        ]


# ---------------------------------------------------------------------------
# _convert_single_step
# ---------------------------------------------------------------------------


class TestConvertSingleStep:
    """Tests for :func:`_convert_single_step`."""

    @pytest.fixture
    def v199_microscale_data(self):
        """Minimal v1.99.0 microscale data."""
        n = 100
        return {
            "params": {
                "micro_params": {
                    "pore_size": "0.000534 centimeter",
                    "total_time": "3600.0 second",
                    "micro_simulations": 100,
                },
                "macro_params": {
                    "diffusion_coeff": "5e-07 centimeter ** 2 / second",
                    "rows": 5,
                    "cols": 3,
                    "total_edges": 57,
                    "total_molecules": 10,
                },
            },
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

    def test_single_step_produces_correct_datasets(self, v199_microscale_data):
        """Single step v1.99.0 → v2.0.0 produces v2.0.0 dataset names."""
        result = _convert_single_step(
            v199_microscale_data, "v1.99.0", "v2.0.0"
        )
        assert "pli_first_time" in result
        assert "sim_final_time" in result

    def test_single_step_converts_params(self):
        """Single step v1.95.0 → v1.99.0 converts params."""
        n = 100
        input_data = {
            "params": {
                "micro_params": {
                    "pore_size": 0.000534,
                    "micro_simulations": 100,
                },
                "macro_params": {
                    "rows": 5,
                    "cols": 3,
                    "total_edges": 57,
                    "total_molecules": 10,
                },
            },
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
        result = _convert_single_step(input_data, "v1.95.0", "v1.99.0")
        assert isinstance(result["params"]["micro_params"]["pore_size"], str)

    def test_single_step_missing_params_raises(self):
        """Missing params raises ValueError."""
        input_data = {"micro_log": np.array(["sim1"], dtype="<U75")}
        with pytest.raises(ValueError, match="parameters"):
            _convert_single_step(input_data, "v1.99.0", "v2.0.0")

    def test_single_step_preserves_params_without_converter(self, v199_microscale_data):
        """When no params_converter exists, params are copied by reference."""
        result = _convert_single_step(
            v199_microscale_data, "v1.99.0", "v2.0.0"
        )
        assert result["params"] is v199_microscale_data["params"]


# ---------------------------------------------------------------------------
# convert_data — multi-step conversion
# ---------------------------------------------------------------------------


class TestConvertDataMultiStep:
    """Tests for multi-step conversion via :func:`convert_data`."""

    @pytest.fixture
    def v195_microscale_data(self):
        """Minimal v1.95.0 microscale data with bare-magnitude params."""
        n = 100
        return {
            "params": {
                "micro_params": {
                    "pore_size": 0.000534,
                    "total_time": 3600.0,
                    "micro_simulations": 100,
                },
                "macro_params": {
                    "rows": 5,
                    "cols": 3,
                    "total_edges": 57,
                    "total_molecules": 10,
                    "diffusion_coeff": 5e-7,
                },
            },
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

    def test_v195_to_v200_auto_chains(self, v195_microscale_data):
        """v1.95.0 → v2.0.0 chains through v1.99.0 automatically."""
        result = convert_data(v195_microscale_data, "v1.95.0", "v2.0.0")
        # v2.0.0 dataset names present
        assert "pli_first_time" in result
        assert "sim_final_time" in result
        assert "tpa_leaving_time" in result

    def test_v195_to_v200_params_converted_through_chain(self, v195_microscale_data):
        """Multi-step conversion applies params converter at each step."""
        result = convert_data(v195_microscale_data, "v1.95.0", "v2.0.0")
        # v1.95.0 → v1.99.0 adds units, v1.99.0 → v2.0.0 has no params converter
        # so params should have unit strings from the v1.95.0 → v1.99.0 step
        assert isinstance(result["params"]["micro_params"]["pore_size"], str)
        q = Q_(result["params"]["micro_params"]["pore_size"])
        assert q.magnitude == pytest.approx(0.000534)

    def test_v200_to_v195_auto_chains(self):
        """v2.0.0 → v1.95.0 chains through v1.99.0 automatically."""
        n = 100
        input_data = {
            "params": {
                "micro_params": {
                    "pore_size": "0.000534 centimeter",
                    "total_time": "3600.0 second",
                    "micro_simulations": 100,
                },
                "macro_params": {
                    "diffusion_coeff": "5e-07 centimeter ** 2 / second",
                    "rows": 5,
                    "cols": 3,
                    "total_edges": 57,
                    "total_molecules": 10,
                },
            },
            "micro_log": np.array(["sim1", "sim2"], dtype="<U75"),
            "pli_first_time": np.array([1.0] * n, dtype=np.float64),
            "tpa_final_num": np.array([0] * n, dtype=np.int32),
            "fiber_degraded": np.array([1] * n, dtype=np.int32),
            "sim_final_time": np.array([10.0] * n, dtype=np.float64),
            "pli_generated_num": np.array([5] * n, dtype=np.int32),
            "tpa_leaving_time": np.array([2.0] * n, dtype=np.float64),
            "tpa_unbound_by_pli": np.array([0] * n, dtype=np.int32),
            "tpa_unbound_kinetic": np.array([1] * n, dtype=np.int32),
        }
        result = convert_data(input_data, "v2.0.0", "v1.95.0")
        # v1.95.0 dataset names present
        assert "firstPLi" in result
        assert "lysis" in result
        # Params should have bare magnitudes (units stripped)
        assert isinstance(result["params"]["micro_params"]["pore_size"], float)
        assert result["params"]["micro_params"]["pore_size"] == pytest.approx(0.000534)

    def test_v195_to_v200_data_values_correct(self, v195_microscale_data):
        """Multi-step conversion preserves data values through chain."""
        result = convert_data(v195_microscale_data, "v1.95.0", "v2.0.0")
        np.testing.assert_array_equal(
            result["pli_first_time"],
            v195_microscale_data["firstPLi"],
        )

    def test_same_spec_short_circuit(self, v195_microscale_data):
        """Same input and output spec returns data as-is."""
        result = convert_data(
            v195_microscale_data, "v1.95.0", "v1.95.0"
        )
        assert result is v195_microscale_data

    def test_same_spec_via_tags_short_circuit(self):
        """Tags that resolve to the same version short-circuit."""
        data = {
            "params": {"micro_params": {}, "macro_params": {}},
            "micro_log": np.array(["sim1"], dtype="<U75"),
        }
        result = convert_data(data, "hdf5", "current")
        # "hdf5" and "current" both resolve to "v2.0.0"
        assert result is data

    def test_no_path_raises_valueerror(self):
        """Requesting a conversion with no path raises ValueError."""
        data = {
            "params": {"micro_params": {}},
            "micro_log": np.array(["sim1"], dtype="<U75"),
        }
        with pytest.raises(ValueError, match="No conversion path"):
            convert_data(data, "v1.95.0", "v999.0.0")

    def test_tag_resolution_with_multi_step(self, v195_microscale_data):
        """Tag aliases work with multi-step conversion."""
        # "current" resolves to "hdf5" → "v2.0.0"
        result = convert_data(v195_microscale_data, "v1.95.0", "current")
        assert "pli_first_time" in result

    def test_partial_data_through_chain(self, v195_microscale_data):
        """Partial data (microscale only) converts through multi-step chain."""
        result = convert_data(v195_microscale_data, "v1.95.0", "v2.0.0")
        assert "pli_first_time" in result
        # Macroscale datasets should not be present
        assert "snapshot_time" not in result
