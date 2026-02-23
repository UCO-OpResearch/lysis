"""Integration tests using real Fortran simulation data.

Tests exercise the complete data pipeline with real output from the
Fortran fibrinolysis simulator:

1. **Fixture tests** (always run in CI):
   Use the truncated fixture at ``tests/fixtures/fortran_sample/``.
   Validates reading, conversion, writing, and macroscale_in generation.

2. **Full dataset tests** (marked ``@pytest.mark.real_data``, local only):
   Use the complete ~665 MB dataset at ``data/2026-02-18-1723/``.
   Skipped when the dataset is not present.
"""

import os

import h5py
import numpy as np
import pytest
from click.testing import CliRunner

from lysis.cli import cli
from lysis.data.dataconvert import convert_data, generate_macroscale_in
from lysis.data.dataspec import dataspec
from lysis.data.fileops import read_data_collection, write_data_collection


# ─── Constants matching the real data ────────────────────────────────────────

MICRO_FILE_CODE = "_PLG2_tPA01_TB-xiii"
MACRO_FILE_CODE = "_TB-xiii__21_105"
TOTAL_MOLECULES = 21105

PARAM_OVERRIDES = {
    "fibrinogen_length": "45nm",
    "fibrinogen_radius": "1.2nm",
    "micro_log_lvl": 40,
    "micro_version": "micro_rates",
    "snap_proportion": 0.66666667,
}


# ─── Helpers ─────────────────────────────────────────────────────────────────

def _read_micro(path):
    """Read microscale_out from *path* (v1.95.0)."""
    return read_data_collection(
        path,
        collections=[dataspec["v1.95.0"]["microscale_out"]],
        file_codes=[MICRO_FILE_CODE],
        param_overrides=PARAM_OVERRIDES,
    )


def _read_micro_and_macro_in(path):
    """Read microscale_out + macroscale_in from *path* (v1.95.0)."""
    return read_data_collection(
        path,
        collections=[
            dataspec["v1.95.0"]["microscale_out"],
            dataspec["v1.95.0"]["macroscale_in"],
        ],
        file_codes=[MICRO_FILE_CODE, MICRO_FILE_CODE],
        param_overrides=PARAM_OVERRIDES,
    )


def _read_all(path):
    """Read all three collections from *path* (v1.95.0)."""
    return read_data_collection(
        path,
        collections=[
            dataspec["v1.95.0"]["microscale_out"],
            dataspec["v1.95.0"]["macroscale_in"],
            dataspec["v1.95.0"]["macroscale_out"],
        ],
        file_codes=[MICRO_FILE_CODE, MICRO_FILE_CODE, MACRO_FILE_CODE],
        param_overrides=PARAM_OVERRIDES,
    )


# ═════════════════════════════════════════════════════════════════════════════
# Fixture tests (always run in CI)
# ═════════════════════════════════════════════════════════════════════════════


class TestReadMicroscale:
    """Read microscale_out from the truncated fixture and verify."""

    @pytest.fixture
    def micro_data(self, fortran_sample_path):
        return _read_micro(fortran_sample_path)

    def test_params_loaded(self, micro_data):
        """Parameters were parsed from the micro log file."""
        micro = micro_data["params"]["micro_params"]
        assert micro["micro_simulations"] == pytest.approx(50000)
        assert micro["nodes_in_micro_row"] == pytest.approx(13)

    def test_overridden_params_present(self, micro_data):
        """Parameters supplied via param_overrides appear in the output."""
        micro = micro_data["params"]["micro_params"]
        assert "fibrinogen_length" in micro
        assert "fibrinogen_radius" in micro
        assert micro["micro_log_lvl"] == 40
        assert micro["micro_version"] == "micro_rates"
        assert micro["snap_proportion"] == pytest.approx(0.66666667)

    @pytest.mark.parametrize("dataset,dtype", [
        ("firstPLi", np.float64),
        ("lysis", np.float64),
        ("tPA_time", np.float64),
        ("lasttPA", np.int32),
        ("lyscomplete", np.int32),
        ("PLi", np.int32),
        ("tPAPLiunbd", np.int32),
        ("tPAunbind", np.int32),
    ])
    def test_binary_dataset_shape_and_dtype(self, micro_data, dataset, dtype):
        arr = micro_data[dataset]
        assert arr.dtype == dtype
        assert arr.shape == (50000,)

    def test_lyscomplete_values_are_0_or_1(self, micro_data):
        assert set(np.unique(micro_data["lyscomplete"])).issubset({0, 1})

    def test_firstPLi_nonnegative(self, micro_data):
        """firstPLi times are >= 0 (zero if no PLi was generated)."""
        assert np.all(micro_data["firstPLi"] >= 0)

    def test_micro_log_is_string_array(self, micro_data):
        log = micro_data["micro_log"]
        assert isinstance(log, np.ndarray)
        assert log.dtype.kind in ("U", "S", "O")


class TestReadMacroscaleIn:
    """Read macroscale_in from the truncated fixture and verify."""

    @pytest.fixture
    def data(self, fortran_sample_path):
        return _read_micro_and_macro_in(fortran_sample_path)

    def test_tPAleave_shape(self, data):
        assert data["tPAleave"].shape == (101,)
        assert data["tPAleave"].dtype == np.float64

    def test_tsectPA_shape(self, data):
        assert data["tsectPA"].shape == (101,)
        assert data["tsectPA"].dtype == np.float64

    def test_lenlysisvect_shape(self, data):
        assert data["lenlysisvect"].shape == (100,)
        assert data["lenlysisvect"].dtype == np.float64

    def test_lysismat_shape(self, data):
        assert data["lysismat"].shape[1] == 100
        assert data["lysismat"].dtype == np.float64

    def test_neighbors_dtype(self, data):
        assert data["neighbors"].dtype == np.int32

    def test_neighbors_values_are_positive(self, data):
        """Fortran 1-based indices should all be >= 1."""
        assert np.all(data["neighbors"] >= 1)


class TestReadMacroscaleOut:
    """Read macroscale_out from the truncated fixture and verify."""

    @pytest.fixture
    def data(self, fortran_sample_path):
        return _read_all(fortran_sample_path)

    def test_macro_params_loaded(self, data):
        macro = data["params"]["macro_params"]
        assert macro["total_molecules"] == TOTAL_MOLECULES
        assert macro["rows"] == 184
        assert macro["cols"] == 19

    def test_one_simulation_loaded(self, data):
        """Only simulation 00 is present in the truncated data."""
        assert len(data["tsave"]) == 1

    def test_tsave_shape(self, data):
        tsave = data["tsave"][0]
        assert tsave.shape == (3,)
        assert tsave.dtype == np.float64

    def test_nsave_is_scalar(self, data):
        nsave = data["Nsave"][0]
        assert nsave.dtype == np.int32
        assert nsave == 2

    def test_m_loc_shape(self, data):
        m_loc = data["m_loc"][0]
        assert m_loc.shape == (3, TOTAL_MOLECULES)
        assert m_loc.dtype == np.int32

    def test_m_loc_values_are_valid_fortran_indices(self, data):
        assert np.all(data["m_loc"][0] >= 1)

    def test_m_bound_shape(self, data):
        m_bound = data["m_bound"][0]
        assert m_bound.shape == (3, TOTAL_MOLECULES)
        assert m_bound.dtype == np.int32

    def test_m_bound_values_are_0_or_1(self, data):
        assert set(np.unique(data["m_bound"][0])).issubset({0, 1})

    def test_mfpt_shape(self, data):
        mfpt = data["mfpt"][0]
        assert mfpt.shape == (TOTAL_MOLECULES,)
        assert mfpt.dtype == np.float64

    def test_f_deg_list_is_structured(self, data):
        f_deg = data["f_deg_list"][0]
        assert f_deg.dtype.names is not None
        assert len(f_deg.dtype.names) == 3

    def test_f_deg_list_timestamps_below_cutoff(self, data):
        """All f_deg_list timestamps are below the last tsave value."""
        cutoff = data["tsave"][0][-1]
        times = data["f_deg_list"][0]["Simulation Time Elapsed"]
        assert np.all(times < cutoff)

    def test_m_bind_t_is_structured(self, data):
        m_bind = data["m_bind_t"][0]
        assert m_bind.dtype.names is not None
        assert len(m_bind.dtype.names) == 4

    def test_m_bind_t_timestamps_below_cutoff(self, data):
        """All m_bind_t timestamps are below the last tsave value."""
        cutoff = data["tsave"][0][-1]
        times = data["m_bind_t"][0]["Simulation Time Elapsed"]
        assert np.all(times < cutoff)

    def test_macro_log_is_string_array(self, data):
        assert len(data["macro_log"]) == 1
        log = data["macro_log"][0]
        assert isinstance(log, np.ndarray)


class TestConvertTruncatedData:
    """Convert truncated Fortran data to HDF5 format and verify."""

    @pytest.fixture
    def converted(self, fortran_sample_path):
        raw = _read_all(fortran_sample_path)
        v199 = convert_data(raw, "v1.95.0", "v1.99.0")
        return convert_data(v199, "v1.99.0", "v2.0.0")

    # ── Parameters ───────────────────────────────────────────────────

    def test_overridden_params_survive_conversion(self, converted):
        """Parameters supplied via param_overrides are in the converted output."""
        micro = converted["params"]["micro_params"]
        assert "fibrinogen_length" in micro
        assert "fibrinogen_radius" in micro
        assert micro["micro_log_lvl"] == 40
        assert micro["micro_version"] == "micro_rates"
        assert micro["snap_proportion"] == pytest.approx(0.66666667)

    # ── Microscale datasets ──────────────────────────────────────────

    def test_microscale_datasets_present(self, converted):
        for name in dataspec["v2.0.0"]["microscale_out"].data:
            assert name in converted, f"Missing: {name}"

    def test_pli_first_time(self, converted):
        """pli_first_time is a single array (simulations_combined)."""
        arr = converted["pli_first_time"]
        assert arr.shape == (50000,)
        assert arr.dtype == np.float64

    def test_fiber_degraded_is_bool(self, converted):
        assert converted["fiber_degraded"].dtype == np.bool_

    # ── Macroscale input datasets ────────────────────────────────────

    def test_macroscale_in_datasets_present(self, converted):
        for name in dataspec["v2.0.0"]["macroscale_in"].data:
            assert name in converted, f"Missing: {name}"

    def test_edge_grid_neighbors_shape(self, converted):
        """edge_grid_neighbors is a single 2D array (simulations_combined)."""
        arr = converted["edge_grid_neighbors"]
        assert arr.ndim == 2
        assert arr.shape[1] == 8
        assert arr.dtype == np.uint32

    def test_edge_grid_neighbors_0_based(self, converted):
        """v2.0.0 uses 0-based grid indices."""
        assert np.min(converted["edge_grid_neighbors"]) >= 0

    # ── Macroscale output datasets ───────────────────────────────────

    def test_macroscale_out_datasets_present(self, converted):
        for name in dataspec["v2.0.0"]["macroscale_out"].data:
            assert name in converted, f"Missing: {name}"

    def test_tpa_location_snapshot(self, converted):
        arr = converted["tpa_location_snapshot"]
        assert len(arr) == 1
        assert arr[0].shape == (TOTAL_MOLECULES, 2, 3)

    def test_fiber_degrade_time_structured(self, converted):
        arr = converted["fiber_degrade_time"]
        assert len(arr) == 1
        assert arr[0].dtype.names is not None
        assert len(arr[0].dtype.names) == 4

    def test_tpa_bind_events_structured(self, converted):
        arr = converted["tpa_bind_events"]
        assert len(arr) == 1
        assert arr[0].dtype.names is not None
        assert len(arr[0].dtype.names) == 5

    def test_snapshot_time(self, converted):
        assert len(converted["snapshot_time"]) == 1
        assert converted["snapshot_time"][0].shape == (3,)


class TestGenerateMacroscaleIn:
    """Validate generate_macroscale_in() against the fixture files.

    Reads the macroscale_in data from the fixture (originally produced by
    the Fortran code), converts to v2.0.0, then independently generates
    macroscale_in from the microscale_out data.  The generated values
    should match the fixture values.
    """

    @pytest.fixture
    def reference_and_generated(self, fortran_sample_path):
        """Return (reference, generated) macroscale_in data in v2.0.0 format."""
        raw = _read_all(fortran_sample_path)
        v199 = convert_data(raw, "v1.95.0", "v1.99.0")
        converted = convert_data(v199, "v1.99.0", "v2.0.0")

        # Save reference macroscale_in before generate modifies data in place.
        # macroscale_in is simulations_combined=True, so convert_data returns
        # single arrays (not lists).
        reference = {
            name: converted[name].copy()
            for name in dataspec["v2.0.0"]["macroscale_in"].data
        }

        generated = generate_macroscale_in(converted)
        return reference, generated

    def test_bin_edge_proportions(self, reference_and_generated):
        ref, gen = reference_and_generated
        np.testing.assert_array_almost_equal(
            gen["bin_edge_proportions"], ref["bin_edge_proportions"]
        )

    def test_bin_edge_tpa_leaving_time(self, reference_and_generated):
        ref, gen = reference_and_generated
        np.testing.assert_array_almost_equal(
            gen["bin_edge_tpa_leaving_time"], ref["bin_edge_tpa_leaving_time"]
        )

    def test_binned_fiber_degrade_time(self, reference_and_generated):
        ref, gen = reference_and_generated
        np.testing.assert_array_almost_equal(
            gen["binned_fiber_degrade_time"], ref["binned_fiber_degrade_time"]
        )

    def test_binned_fiber_degraded(self, reference_and_generated):
        ref, gen = reference_and_generated
        np.testing.assert_array_equal(
            gen["binned_fiber_degraded"], ref["binned_fiber_degraded"]
        )

    def test_edge_grid_neighbors(self, reference_and_generated):
        ref, gen = reference_and_generated
        np.testing.assert_array_equal(
            gen["edge_grid_neighbors"], ref["edge_grid_neighbors"]
        )


class TestWriteConvertedData:
    """Write converted HDF5 data and verify the output file."""

    @pytest.fixture
    def converted_h5(self, fortran_sample_path, tmp_path):
        raw = _read_all(fortran_sample_path)
        v199 = convert_data(raw, "v1.95.0", "v1.99.0")
        converted = convert_data(v199, "v1.99.0", "v2.0.0")
        h5_path = str(tmp_path / "converted.h5")
        collections_out = [
            dataspec["v2.0.0"]["microscale_out"],
            dataspec["v2.0.0"]["macroscale_out"],
        ]
        write_data_collection(converted, h5_path, collections_out, ["", ""])
        return h5_path

    def test_hdf5_file_created(self, converted_h5):
        assert os.path.exists(converted_h5)

    def test_hdf5_has_version_attr(self, converted_h5):
        with h5py.File(converted_h5, "r") as f:
            assert "dataspec_version" in f.attrs
            assert f.attrs["dataspec_version"] == "v2.0.0"

    def test_hdf5_has_micro_data_group(self, converted_h5):
        with h5py.File(converted_h5, "r") as f:
            assert "micro_data" in f

    def test_hdf5_has_macro_data_group(self, converted_h5):
        with h5py.File(converted_h5, "r") as f:
            assert "macro_data" in f
            assert "macro_data/sim_00" in f

    def test_hdf5_pli_first_time_shape(self, converted_h5):
        with h5py.File(converted_h5, "r") as f:
            ds = f["micro_data/pli_first_time"]
            assert ds.shape == (50000,)
            assert ds.dtype == np.float64


class TestCLIConvertTruncated:
    """End-to-end test of ``lysis convert`` with the truncated fixture."""

    @pytest.mark.skip(
        reason="CLI does a single convert_data() call; v1.95.0 → v2.0.0 "
        "requires multi-hop routing (v1.95.0 → v1.99.0 → v2.0.0) "
        "which the CLI does not yet support."
    )
    def test_cli_convert_fortran_to_hdf5(self, fortran_sample_path, tmp_path):
        output_path = str(tmp_path / "cli_output.h5")
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "convert",
                fortran_sample_path,
                output_path,
                "-f", "v1.95.0",
                "-t", "hdf5",
                "--file-code",
                f"{MICRO_FILE_CODE},{MACRO_FILE_CODE}",
                "--collections", "microscale_out,macroscale_out",
                "--param-override", "fibrinogen_length=45nm",
                "--param-override", "fibrinogen_radius=1.2nm",
                "--param-override", "micro_log_lvl=40",
                "--param-override", "micro_version=micro_rates",
                "--param-override", "snap_proportion=0.66666667",
            ],
        )
        assert result.exit_code == 0, f"CLI failed:\n{result.output}"
        assert os.path.exists(output_path)

        with h5py.File(output_path, "r") as f:
            assert "micro_data" in f
            assert "macro_data" in f


# ═════════════════════════════════════════════════════════════════════════════
# Full dataset tests (local only)
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.real_data
class TestReadFullMicroscale:
    """Read microscale_out from the full dataset and verify."""

    @pytest.fixture
    def micro_data(self, full_data_path):
        return _read_micro(full_data_path)

    def test_all_50000_simulations(self, micro_data):
        assert micro_data["lysis"].shape == (50000,)

    @pytest.mark.parametrize("dataset,dtype,expected_shape", [
        ("firstPLi", np.float64, (50000,)),
        ("lysis", np.float64, (50000,)),
        ("tPA_time", np.float64, (50000,)),
        ("lasttPA", np.int32, (50000,)),
        ("lyscomplete", np.int32, (50000,)),
        ("PLi", np.int32, (50000,)),
        ("tPAPLiunbd", np.int32, (50000,)),
        ("tPAunbind", np.int32, (50000,)),
    ])
    def test_binary_dataset_shape(self, micro_data, dataset, dtype, expected_shape):
        arr = micro_data[dataset]
        assert arr.dtype == dtype
        assert arr.shape == expected_shape


@pytest.mark.real_data
class TestReadFullMacroscale:
    """Read macroscale_out from the full dataset and verify."""

    @pytest.fixture
    def data(self, full_data_path):
        return _read_all(full_data_path)

    def test_ten_simulations_loaded(self, data):
        assert len(data["tsave"]) == 10

    def test_m_loc_full_shape(self, data):
        m_loc = data["m_loc"][0]
        assert m_loc.shape[1] == TOTAL_MOLECULES
        assert m_loc.shape[0] > 100  # hundreds of snapshots


@pytest.mark.real_data
class TestConvertFullData:
    """Convert the full dataset and verify."""

    def test_full_conversion_completes(self, full_data_path):
        data = _read_all(full_data_path)
        v199 = convert_data(data, "v1.95.0", "v1.99.0")
        converted = convert_data(v199, "v1.99.0", "v2.0.0")
        assert "pli_first_time" in converted
        assert "tpa_location_snapshot" in converted
        assert len(converted["tpa_location_snapshot"]) == 10
