"""Comprehensive pytest unit tests for :mod:`lysis.dataio.fileops`.

Tests cover:

* Low-level readers (``_read_file_text``, ``_read_file_binary``,
  ``_read_file_json``, ``_read_file_parsed``, ``_read_hdf5_dataset``,
  ``_read_hdf5_attr``)
* Low-level writers (``_write_file_text``, ``_write_file_binary``,
  ``_write_file_json``, ``_write_hdf5_dataset``, ``_write_hdf5_attr``)
* Mid-level dispatch (``read_dataset``, ``write_dataset``)
* High-level collection I/O (``read_data_collection``,
  ``write_data_collection``)
"""

import json
import textwrap

import h5py
import numpy as np
import pytest

from pint import Quantity

from lysis.config.constants import CONST
from lysis.dataio.dataspec import (
    DataCollectionSpec,
    DataSetSpec,
    check_dataset_spec,
    dataspec,
    parse_shape,
)
from lysis.dataio.fileops import (
    _read_file_binary,
    _read_file_json,
    _read_file_parsed,
    _read_file_text,
    _read_hdf5_attr,
    _read_hdf5_dataset,
    _validate_hdf5_version,
    _write_file_binary,
    _write_file_json,
    _write_file_text,
    _write_hdf5_attr,
    _write_hdf5_dataset,
    data_readers,
    data_writers,
    ensure_hdf5_version,
    parse_macro_log,
    parse_micro_file_code,
    parse_micro_log,
    read_data_collection,
    read_dataset,
    write_data_collection,
    write_dataset,
)


# ---------------------------------------------------------------------------
# Low-level readers
# ---------------------------------------------------------------------------


class TestReadFileText:
    """Tests for :func:`_read_file_text`."""

    def test_float64_array(self, tmp_path, text_spec):
        """Read a (5,) float64 array written with np.savetxt."""
        original = np.array([1.1, 2.2, 3.3, 4.4, 5.5], dtype=np.float64)
        np.savetxt(tmp_path / "data.dat", original)
        result = _read_file_text(str(tmp_path), text_spec)
        np.testing.assert_array_almost_equal(result, original)

    def test_structured_array_with_delimiter(self, tmp_path):
        """Read a structured array from a delimited text file."""
        dt = np.dtype([("time", np.float64), ("idx", np.int32), ("val", np.float64)])
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
            dtype=dt,
            data_location="struct{file_code}.dat",
            shape=(-1,),
            delimiter=",",
        )
        rows = np.array(
            [(1.0, 2, 3.0), (4.0, 5, 6.0)],
            dtype=dt,
        )
        np.savetxt(
            tmp_path / "struct.dat",
            np.column_stack([rows["time"], rows["idx"], rows["val"]]),
            delimiter=",",
            fmt=["%.18e", "%d", "%.18e"],
        )
        result = _read_file_text(str(tmp_path), spec)
        # loadtxt with structured dtype returns 1D structured array
        np.testing.assert_array_almost_equal(result["time"], rows["time"])
        np.testing.assert_array_almost_equal(result["idx"], rows["idx"])
        np.testing.assert_array_almost_equal(result["val"], rows["val"])

    def test_comma_delimiter(self, tmp_path):
        """Ensure that a comma delimiter is respected during read."""
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
            dtype=np.float64,
            data_location="csv{file_code}.dat",
            shape=(-1, 3),
            delimiter=",",
        )
        original = np.arange(12, dtype=np.float64).reshape(4, 3)
        np.savetxt(tmp_path / "csv.dat", original, delimiter=",")
        result = _read_file_text(str(tmp_path), spec)
        np.testing.assert_array_almost_equal(result, original)

    def test_missing_text_file(self, tmp_path, text_spec):
        """Raise FileNotFoundError when the text file does not exist."""
        with pytest.raises((FileNotFoundError, OSError)):
            _read_file_text(str(tmp_path), text_spec)


class TestReadFileBinary:
    """Tests for :func:`_read_file_binary`."""

    def test_int32_array_reshape(self, tmp_path, binary_spec):
        """Read a (3,4) int32 binary array and verify reshape."""
        original = np.arange(12, dtype=np.int32).reshape(3, 4)
        original.tofile(tmp_path / "data.dat")
        result = _read_file_binary(str(tmp_path), binary_spec)
        np.testing.assert_array_equal(result, original)

    def test_reshape_from_params(self, tmp_path, sample_params):
        """Resolve shape from params using (-1, 'macro_params.total_molecules')."""
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
            dtype=np.float64,
            data_location="param_shape{file_code}.dat",
            shape=(-1, "macro_params.total_molecules"),
        )
        # total_molecules = 10 in sample_params
        original = np.arange(30, dtype=np.float64).reshape(3, 10)
        original.tofile(tmp_path / "param_shape.dat")
        result = _read_file_binary(str(tmp_path), spec, params=sample_params)
        np.testing.assert_array_equal(result, original)
        assert result.shape == (3, 10)

    def test_missing_binary_file(self, tmp_path, binary_spec):
        """Raise FileNotFoundError when the binary file does not exist."""
        with pytest.raises((FileNotFoundError, OSError)):
            _read_file_binary(str(tmp_path), binary_spec)


class TestReadFileJson:
    """Tests for :func:`_read_file_json`."""

    def test_round_trip(self, tmp_path, json_spec):
        """Write a params dict with json.dump, then read it back."""
        params = {
            "micro_params": {"micro_simulations": 100},
            "macro_params": {"total_molecules": 10},
        }
        with open(tmp_path / "params.json", "w") as fh:
            json.dump(params, fh)
        result = _read_file_json(str(tmp_path), json_spec)
        assert result == params


class TestReadHdf5Dataset:
    """Tests for :func:`_read_hdf5_dataset`."""

    def test_read_dataset(self, tmp_path, hdf5_dataset_spec):
        """Create an h5 file with a dataset and verify read matches."""
        original = np.linspace(0, 1, 20, dtype=np.float64)
        h5_path = str(tmp_path / "test.h5")
        with h5py.File(h5_path, "w") as fh:
            fh.create_dataset("data/arr", data=original)
        result = _read_hdf5_dataset(h5_path, hdf5_dataset_spec)
        np.testing.assert_array_equal(result, original)


class TestReadHdf5Attr:
    """Tests for :func:`_read_hdf5_attr`."""

    def test_read_attrs(self, tmp_path, hdf5_attr_spec):
        """Create an h5 file with group attributes and verify read matches."""
        h5_path = str(tmp_path / "test.h5")
        with h5py.File(h5_path, "w") as fh:
            grp = fh.require_group("test_data")
            grp.attrs["alpha"] = 1.5
            grp.attrs["beta"] = 42
        result = _read_hdf5_attr(h5_path, hdf5_attr_spec)
        assert "test_data" in result
        assert result["test_data"]["alpha"] == 1.5
        assert result["test_data"]["beta"] == 42


# ---------------------------------------------------------------------------
# Low-level writers
# ---------------------------------------------------------------------------


class TestWriteFileText:
    """Tests for :func:`_write_file_text`."""

    def test_write_float64(self, tmp_path, text_spec):
        """Write a float64 array to text and verify with np.loadtxt."""
        original = np.array([1.1, 2.2, 3.3], dtype=np.float64)
        _write_file_text(original, str(tmp_path), text_spec)
        outfile = tmp_path / "data.dat"
        assert outfile.exists()
        loaded = np.loadtxt(str(outfile), dtype=np.float64)
        np.testing.assert_array_almost_equal(loaded, original)

    def test_spec_mismatch_raises(self, tmp_path):
        """Raise TypeError when data does not match the spec."""
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
            dtype=np.float64,
            data_location="bad{file_code}.dat",
            shape=(3, 4),
        )
        wrong_data = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float64)
        with pytest.raises(TypeError):
            _write_file_text(wrong_data, str(tmp_path), spec)


class TestWriteFileBinary:
    """Tests for :func:`_write_file_binary`."""

    def test_write_int32(self, tmp_path, binary_spec):
        """Write an int32 array to binary and verify with np.fromfile."""
        original = np.arange(12, dtype=np.int32).reshape(3, 4)
        _write_file_binary(original, str(tmp_path), binary_spec)
        outfile = tmp_path / "data.dat"
        assert outfile.exists()
        loaded = np.fromfile(str(outfile), dtype=np.int32).reshape(3, 4)
        np.testing.assert_array_equal(loaded, original)

    def test_spec_mismatch_raises(self, tmp_path, binary_spec):
        """Raise TypeError when shape does not match the spec."""
        wrong_data = np.array([1, 2, 3, 4, 5], dtype=np.int32)
        with pytest.raises(TypeError):
            _write_file_binary(wrong_data, str(tmp_path), binary_spec)


class TestWriteFileJson:
    """Tests for :func:`_write_file_json`."""

    def test_write_params(self, tmp_path, json_spec):
        """Write a params dict and verify with json.load."""
        params = {"micro_params": {"micro_simulations": 50}}
        _write_file_json(params, str(tmp_path), json_spec)
        outfile = tmp_path / "params.json"
        assert outfile.exists()
        with open(outfile, "r") as fh:
            loaded = json.load(fh)
        assert loaded == params

    def test_non_serializable_converted_to_str(self, tmp_path, json_spec):
        """Non-JSON-serializable values (e.g. Pint Quantity) must be stored as strings."""
        q = Quantity("5 nanometer")
        params = {"micro_params": {"fiber_radius": q}}
        _write_file_json(params, str(tmp_path), json_spec)
        with open(tmp_path / "params.json", "r") as fh:
            loaded = json.load(fh)
        assert isinstance(loaded["micro_params"]["fiber_radius"], str)
        assert "5" in loaded["micro_params"]["fiber_radius"]

    def test_pure_serializable_unchanged(self, tmp_path, json_spec):
        """Plain int/float/str values must survive the round-trip unchanged."""
        params = {"micro_params": {"micro_simulations": 42, "label": "test"}}
        _write_file_json(params, str(tmp_path), json_spec)
        with open(tmp_path / "params.json", "r") as fh:
            loaded = json.load(fh)
        assert loaded["micro_params"]["micro_simulations"] == 42
        assert loaded["micro_params"]["label"] == "test"


class TestWriteHdf5Dataset:
    """Tests for :func:`_write_hdf5_dataset`."""

    def test_write_float64(self, tmp_path, hdf5_dataset_spec):
        """Write a float64 array to HDF5 and verify with h5py read."""
        original = np.linspace(0, 1, 15, dtype=np.float64)
        h5_path = str(tmp_path / "out.h5")
        _write_hdf5_dataset(original, h5_path, hdf5_dataset_spec)
        with h5py.File(h5_path, "r") as fh:
            loaded = fh["data/arr"][:]
        np.testing.assert_array_almost_equal(loaded, original)

    def test_overwrite_true_replaces_existing_dataset(self, tmp_path, hdf5_dataset_spec):
        """overwrite=True must delete a pre-existing dataset and write the new data."""
        h5_path = str(tmp_path / "out.h5")
        original = np.zeros(15, dtype=np.float64)
        replacement = np.linspace(0, 1, 15, dtype=np.float64)

        _write_hdf5_dataset(original, h5_path, hdf5_dataset_spec)
        _write_hdf5_dataset(replacement, h5_path, hdf5_dataset_spec, overwrite=True)

        with h5py.File(h5_path, "r") as fh:
            loaded = fh["data/arr"][:]
        np.testing.assert_array_almost_equal(loaded, replacement)

    def test_overwrite_false_raises_on_existing_dataset(self, tmp_path, hdf5_dataset_spec):
        """Without overwrite=True, writing to an existing dataset must raise."""
        import h5py as _h5py

        h5_path = str(tmp_path / "out.h5")
        data = np.zeros(15, dtype=np.float64)
        _write_hdf5_dataset(data, h5_path, hdf5_dataset_spec)
        with pytest.raises(Exception):
            _write_hdf5_dataset(data, h5_path, hdf5_dataset_spec, overwrite=False)

    def test_overwrite_true_empty_dataset_correct_shape(self, tmp_path, hdf5_dataset_spec):
        """overwrite=True on a zero-length pre-allocated dataset must produce the right shape."""
        h5_path = str(tmp_path / "out.h5")
        replacement = np.linspace(0, 1, 15, dtype=np.float64)

        # Pre-allocate an empty (zero-length) resizable dataset like DataStore.create() does
        with h5py.File(h5_path, "a") as fh:
            fh.attrs[CONST.DATASPEC_VERSION_ATTR] = hdf5_dataset_spec.version
            fh.create_dataset("data/arr", shape=(0,), maxshape=(None,), dtype=np.float64)

        _write_hdf5_dataset(replacement, h5_path, hdf5_dataset_spec, overwrite=True)

        with h5py.File(h5_path, "r") as fh:
            loaded = fh["data/arr"][:]
        assert loaded.shape == replacement.shape
        np.testing.assert_array_almost_equal(loaded, replacement)


class TestWriteHdf5Attr:
    """Tests for :func:`_write_hdf5_attr`."""

    def test_write_attrs(self, tmp_path, hdf5_attr_spec):
        """Write params as HDF5 group attributes and verify with h5py."""
        # data_location is "test_data", param group name is "test_params"
        data = {"test_params": {"key1": 3.14, "key2": 99}}
        h5_path = str(tmp_path / "out.h5")
        _write_hdf5_attr(data, h5_path, hdf5_attr_spec)
        with h5py.File(h5_path, "r") as fh:
            grp = fh["test_data"]
            assert grp.attrs["key1"] == 3.14
            assert grp.attrs["key2"] == 99


# ---------------------------------------------------------------------------
# Mid-level dispatch: read_dataset / write_dataset
# ---------------------------------------------------------------------------


class TestReadWriteDatasetDispatch:
    """Round-trip tests for :func:`read_dataset` and :func:`write_dataset`."""

    def test_text_round_trip(self, tmp_path, text_spec):
        """Write then read a text dataset and verify equality."""
        original = np.array([10.0, 20.0, 30.0], dtype=np.float64)
        write_dataset(original, str(tmp_path), text_spec)
        result = read_dataset(str(tmp_path), text_spec)
        np.testing.assert_array_almost_equal(result, original)

    def test_binary_round_trip(self, tmp_path, binary_spec):
        """Write then read a binary dataset and verify equality."""
        original = np.arange(12, dtype=np.int32).reshape(3, 4)
        write_dataset(original, str(tmp_path), binary_spec)
        result = read_dataset(str(tmp_path), binary_spec)
        np.testing.assert_array_equal(result, original)

    def test_json_round_trip(self, tmp_path, json_spec):
        """Write then read a JSON params dataset and verify equality."""
        params = {"macro_params": {"total_molecules": 42}}
        write_dataset(params, str(tmp_path), json_spec)
        result = read_dataset(str(tmp_path), json_spec)
        assert result == params

    def test_hdf5_dataset_round_trip(self, tmp_path, hdf5_dataset_spec):
        """Write then read an HDF5 dataset and verify equality."""
        original = np.arange(8, dtype=np.float64)
        h5_path = str(tmp_path / "rt.h5")
        write_dataset(original, h5_path, hdf5_dataset_spec)
        result = read_dataset(h5_path, hdf5_dataset_spec)
        np.testing.assert_array_equal(result, original)

    def test_hdf5_attr_round_trip(self, tmp_path, hdf5_attr_spec):
        """Write then read HDF5 attributes and verify equality."""
        data = {"test_params": {"x": 1.0, "y": 2.0}}
        h5_path = str(tmp_path / "rt.h5")
        write_dataset(data, h5_path, hdf5_attr_spec)
        result = read_dataset(h5_path, hdf5_attr_spec)
        assert result["test_data"]["x"] == 1.0
        assert result["test_data"]["y"] == 2.0


# ---------------------------------------------------------------------------
# High-level: read_data_collection / write_data_collection
# ---------------------------------------------------------------------------


class TestReadDataCollectionCombined:
    """Tests for :func:`read_data_collection` with simulations_combined=True."""

    def test_combined_single_array(self, tmp_path, combined_collection_spec, json_spec):
        """Combined collection returns a single array (not a list)."""
        # Write params file
        params = {"micro_params": {"micro_simulations": 10}}
        with open(tmp_path / "params.json", "w") as fh:
            json.dump(params, fh)
        # Write data file
        arr = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        np.savetxt(tmp_path / "data.dat", arr)

        data = read_data_collection(
            str(tmp_path),
            collections=[combined_collection_spec],
            file_codes=[""],
        )
        assert isinstance(data["arr"], np.ndarray)
        np.testing.assert_array_almost_equal(data["arr"], arr)

    def test_with_params(self, tmp_path, combined_collection_spec):
        """Collection with params spec populates data['params']."""
        params = {"micro_params": {"micro_simulations": 10}}
        with open(tmp_path / "params.json", "w") as fh:
            json.dump(params, fh)
        arr = np.array([1.0], dtype=np.float64)
        np.savetxt(tmp_path / "data.dat", arr)

        data = read_data_collection(
            str(tmp_path),
            collections=[combined_collection_spec],
            file_codes=[""],
        )
        assert data["params"]["micro_params"]["micro_simulations"] == 10

    def test_no_params(self, tmp_path, text_spec):
        """Collection with params=None leaves data['params'] as empty dict."""
        no_params_collection = DataCollectionSpec(
            simulations_combined=True,
            params=None,
            data={"arr": text_spec},
        )
        arr = np.array([1.0, 2.0], dtype=np.float64)
        np.savetxt(tmp_path / "data.dat", arr)

        data = read_data_collection(
            str(tmp_path),
            collections=[no_params_collection],
            file_codes=[""],
        )
        assert data["params"] == {}


class TestReadDataCollectionPerSim:
    """Tests for :func:`read_data_collection` with simulations_combined=False."""

    def _setup_per_sim_files(self, tmp_path, n_sims=3):
        """Helper: create params.json and per-sim binary files in sub-dirs."""
        params = {"micro_params": {"micro_simulations": 10}}
        with open(tmp_path / "params.json", "w") as fh:
            json.dump(params, fh)

        arrays = []
        for i in range(n_sims):
            sim_dir = tmp_path / f"{i:02}"
            sim_dir.mkdir(exist_ok=True)
            arr = np.arange(5, dtype=np.float64) + i * 10.0
            arr.tofile(sim_dir / "data.dat")
            arrays.append(arr)
        return arrays

    def test_per_sim_list_of_arrays(self, tmp_path, per_sim_collection_spec):
        """Per-sim read with 3 files returns list of 3 arrays."""
        expected = self._setup_per_sim_files(tmp_path, n_sims=3)
        data = read_data_collection(
            str(tmp_path),
            collections=[per_sim_collection_spec],
            file_codes=[""],
        )
        assert isinstance(data["arr"], list)
        assert len(data["arr"]) == 3
        for i in range(3):
            np.testing.assert_array_equal(data["arr"][i], expected[i])

    def test_per_sim_graceful_stop(self, tmp_path, per_sim_collection_spec):
        """Reading stops gracefully when sim file does not exist (after sim 0)."""
        self._setup_per_sim_files(tmp_path, n_sims=3)
        # Only 3 directories exist (00, 01, 02). Reading stops at 03.
        data = read_data_collection(
            str(tmp_path),
            collections=[per_sim_collection_spec],
            file_codes=[""],
        )
        assert len(data["arr"]) == 3

    def test_per_sim_sim0_missing_raises(self, tmp_path, per_sim_collection_spec):
        """Raise FileNotFoundError when sim 0 file is missing."""
        # Only write params, no sim directories
        params = {"micro_params": {"micro_simulations": 10}}
        with open(tmp_path / "params.json", "w") as fh:
            json.dump(params, fh)

        with pytest.raises(FileNotFoundError):
            read_data_collection(
                str(tmp_path),
                collections=[per_sim_collection_spec],
                file_codes=[""],
            )


class TestDataCollectionOptional:
    """Tests that the ``optional`` DataSetSpec flag is honored end-to-end.

    Covers both a combined (``simulations_combined=True``, like the existing
    ``neighbors`` dataset) and a per-simulation optional dataset, on both the
    read and write paths. See GitHub issue #83.
    """

    def _optional_combined_spec(self, json_spec):
        """Combined collection with a required ``arr`` and an optional ``opt``."""
        required = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
            dtype=np.float64,
            data_location="arr{file_code}.dat",
            shape=(-1,),
        )
        optional = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
            dtype=np.float64,
            data_location="opt{file_code}.dat",
            shape=(-1,),
            optional=True,
        )
        return DataCollectionSpec(
            simulations_combined=True,
            params=json_spec,
            data={"arr": required, "opt": optional},
        )

    def _optional_per_sim_spec(self, json_spec):
        """Per-sim collection with a required ``arr`` and an optional ``opt``."""
        required = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
            dtype=np.float64,
            data_location="{sim:02}/arr{file_code}.dat",
            shape=(-1,),
        )
        optional = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
            dtype=np.float64,
            data_location="{sim:02}/opt{file_code}.dat",
            shape=(-1,),
            optional=True,
        )
        return DataCollectionSpec(
            simulations_combined=False,
            params=json_spec,
            data={"arr": required, "opt": optional},
        )

    def test_combined_optional_absent_no_raise(self, tmp_path, json_spec):
        """Combined: an absent optional dataset does not raise (unchanged)."""
        spec = self._optional_combined_spec(json_spec)
        with open(tmp_path / "params.json", "w") as fh:
            json.dump({"micro_params": {"micro_simulations": 10}}, fh)
        # Write only the required dataset; the optional one is absent on disk.
        np.savetxt(tmp_path / "arr.dat", np.array([1.0, 2.0, 3.0]))

        data = read_data_collection(
            str(tmp_path),
            collections=[spec],
            file_codes=[""],
        )
        np.testing.assert_array_almost_equal(data["arr"], [1.0, 2.0, 3.0])
        assert "opt" not in data

    def test_per_sim_optional_absent_returns_empty(self, tmp_path, json_spec):
        """Per-sim: an entirely absent optional dataset reads as an empty list."""
        spec = self._optional_per_sim_spec(json_spec)
        with open(tmp_path / "params.json", "w") as fh:
            json.dump({"micro_params": {"micro_simulations": 10}}, fh)
        # Write the required dataset for 2 sims; never write the optional one.
        for i in range(2):
            sim_dir = tmp_path / f"{i:02}"
            sim_dir.mkdir(exist_ok=True)
            (np.arange(5, dtype=np.float64) + i * 10.0).tofile(sim_dir / "arr.dat")

        data = read_data_collection(
            str(tmp_path),
            collections=[spec],
            file_codes=[""],
        )
        assert len(data["arr"]) == 2
        assert data["opt"] == []

    def test_write_optional_absent_skipped(self, tmp_path, json_spec):
        """write_data_collection skips an optional key absent from the data dict."""
        spec = self._optional_per_sim_spec(json_spec)
        arrays = [np.arange(5, dtype=np.float64) + i * 100.0 for i in range(2)]
        # Note: no "opt" key in the data dict at all.
        data = {
            "params": {"micro_params": {"micro_simulations": 10}},
            "arr": arrays,
        }
        for i in range(2):
            (tmp_path / f"{i:02}").mkdir(exist_ok=True)

        # Must not raise (previously a KeyError on the absent "opt" key).
        write_data_collection(
            data=data,
            path=str(tmp_path),
            collections=[spec],
            file_codes=[""],
        )

        # The optional dataset was skipped, the required one written.
        for i in range(2):
            assert (tmp_path / f"{i:02}" / "arr.dat").exists()
            assert not (tmp_path / f"{i:02}" / "opt.dat").exists()

    def test_per_sim_optional_round_trip(self, tmp_path, json_spec):
        """Per-sim: an optional dataset that IS present round-trips normally."""
        spec = self._optional_per_sim_spec(json_spec)
        arrs = [np.arange(5, dtype=np.float64) + i * 100.0 for i in range(2)]
        opts = [np.arange(3, dtype=np.float64) + i for i in range(2)]
        data = {
            "params": {"micro_params": {"micro_simulations": 10}},
            "arr": arrs,
            "opt": opts,
        }
        for i in range(2):
            (tmp_path / f"{i:02}").mkdir(exist_ok=True)

        write_data_collection(
            data=data,
            path=str(tmp_path),
            collections=[spec],
            file_codes=[""],
        )
        loaded = read_data_collection(
            str(tmp_path),
            collections=[spec],
            file_codes=[""],
        )
        assert len(loaded["opt"]) == 2
        for i in range(2):
            np.testing.assert_array_equal(loaded["opt"][i], opts[i])


class TestWriteCollectionParamsRequired:
    """Tests that write_data_collection requires parameters."""

    def test_missing_params_raises(self, tmp_path, combined_collection_spec):
        """Writing without params key raises ValueError."""
        data = {"arr": np.array([1.0, 2.0], dtype=np.float64)}
        with pytest.raises(ValueError, match="parameters"):
            write_data_collection(
                data=data,
                path=str(tmp_path),
                collections=[combined_collection_spec],
                file_codes=[""],
            )

    def test_empty_params_raises(self, tmp_path, combined_collection_spec):
        """Writing with empty params dict raises ValueError."""
        data = {
            "params": {},
            "arr": np.array([1.0, 2.0], dtype=np.float64),
        }
        with pytest.raises(ValueError, match="parameters"):
            write_data_collection(
                data=data,
                path=str(tmp_path),
                collections=[combined_collection_spec],
                file_codes=[""],
            )


class TestWriteReadCollectionRoundTrip:
    """Round-trip tests for write_data_collection + read_data_collection."""

    def test_combined_round_trip(self, tmp_path, combined_collection_spec):
        """Write then read a combined collection; data matches."""
        original_data = {
            "params": {"micro_params": {"micro_simulations": 10}},
            "arr": np.array([5.0, 6.0, 7.0, 8.0], dtype=np.float64),
        }
        write_data_collection(
            data=original_data,
            path=str(tmp_path),
            collections=[combined_collection_spec],
            file_codes=[""],
        )
        loaded = read_data_collection(
            str(tmp_path),
            collections=[combined_collection_spec],
            file_codes=[""],
        )
        np.testing.assert_array_almost_equal(loaded["arr"], original_data["arr"])
        assert loaded["params"]["micro_params"]["micro_simulations"] == 10

    def test_per_sim_round_trip(self, tmp_path, per_sim_collection_spec):
        """Write 3 per-sim arrays then read back; all 3 match."""
        arrays = [np.arange(5, dtype=np.float64) + i * 100.0 for i in range(3)]
        original_data = {
            "params": {"micro_params": {"micro_simulations": 10}},
            "arr": arrays,
        }

        # Create sim directories (write_data_collection expects them)
        for i in range(3):
            (tmp_path / f"{i:02}").mkdir(exist_ok=True)

        write_data_collection(
            data=original_data,
            path=str(tmp_path),
            collections=[per_sim_collection_spec],
            file_codes=[""],
        )
        loaded = read_data_collection(
            str(tmp_path),
            collections=[per_sim_collection_spec],
            file_codes=[""],
        )
        assert len(loaded["arr"]) == 3
        for i in range(3):
            np.testing.assert_array_equal(loaded["arr"][i], arrays[i])


# ---------------------------------------------------------------------------
# ensure_hdf5_version tests
# ---------------------------------------------------------------------------


class TestEnsureHdf5Version:
    """Tests for :func:`ensure_hdf5_version`."""

    def test_creates_new_file_with_version(self, tmp_path):
        """When file doesn't exist, creates it with the version attribute."""
        h5_path = str(tmp_path / "new.h5")
        ensure_hdf5_version(h5_path, "v2.0.0")
        with h5py.File(h5_path, "r") as f:
            assert f.attrs[CONST.DATASPEC_VERSION_ATTR] == "v2.0.0"

    def test_existing_file_correct_version(self, tmp_path):
        """When file exists with correct version, no error raised."""
        h5_path = str(tmp_path / "ok.h5")
        with h5py.File(h5_path, "w") as f:
            f.attrs[CONST.DATASPEC_VERSION_ATTR] = "v2.0.0"
        ensure_hdf5_version(h5_path, "v2.0.0")  # Should not raise

    def test_existing_file_wrong_version_raises(self, tmp_path):
        """When file exists with wrong version, raises ValueError."""
        h5_path = str(tmp_path / "wrong.h5")
        with h5py.File(h5_path, "w") as f:
            f.attrs[CONST.DATASPEC_VERSION_ATTR] = "v1.0.0"
        with pytest.raises(ValueError, match="mismatch"):
            ensure_hdf5_version(h5_path, "v2.0.0")

    def test_existing_file_missing_version_raises(self, tmp_path):
        """When file exists with no version attr, raises ValueError."""
        h5_path = str(tmp_path / "nover.h5")
        with h5py.File(h5_path, "w"):
            pass
        with pytest.raises(ValueError, match=CONST.DATASPEC_VERSION_ATTR):
            ensure_hdf5_version(h5_path, "v2.0.0")

    def test_skips_when_no_version(self, tmp_path):
        """When version is empty string, skips all checks."""
        h5_path = str(tmp_path / "skip.h5")
        # File doesn't exist, empty version => no file created
        ensure_hdf5_version(h5_path, "")
        assert not (tmp_path / "skip.h5").exists()


# ---------------------------------------------------------------------------
# HDF5 version validation in readers
# ---------------------------------------------------------------------------


class TestHdf5VersionValidation:
    """Tests for version validation in HDF5 read functions."""

    def test_read_hdf5_dataset_with_correct_version(self, tmp_path):
        """Reading with matching version succeeds."""
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
            dtype=np.float64,
            data_location="data/arr",
        )
        object.__setattr__(spec, "version", "v2.0.0")
        h5_path = str(tmp_path / "ok.h5")
        with h5py.File(h5_path, "w") as f:
            f.attrs[CONST.DATASPEC_VERSION_ATTR] = "v2.0.0"
            f.create_dataset("data/arr", data=np.array([1.0, 2.0]))
        result = _read_hdf5_dataset(h5_path, spec)
        np.testing.assert_array_equal(result, [1.0, 2.0])

    def test_read_hdf5_dataset_wrong_version_raises(self, tmp_path):
        """Reading with mismatched version raises ValueError."""
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
            dtype=np.float64,
            data_location="data/arr",
        )
        object.__setattr__(spec, "version", "v2.0.0")
        h5_path = str(tmp_path / "wrong.h5")
        with h5py.File(h5_path, "w") as f:
            f.attrs[CONST.DATASPEC_VERSION_ATTR] = "v1.0.0"
            f.create_dataset("data/arr", data=np.array([1.0, 2.0]))
        with pytest.raises(ValueError, match="mismatch"):
            _read_hdf5_dataset(h5_path, spec)

    def test_read_hdf5_attr_wrong_version_raises(self, tmp_path):
        """Reading attrs with mismatched version raises ValueError."""
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_ATTR,
            dtype=dict,
            data_location="test_data",
        )
        object.__setattr__(spec, "version", "v2.0.0")
        h5_path = str(tmp_path / "wrong.h5")
        with h5py.File(h5_path, "w") as f:
            f.attrs[CONST.DATASPEC_VERSION_ATTR] = "v1.0.0"
            grp = f.require_group("test_data")
            grp.attrs["key"] = 42
        with pytest.raises(ValueError, match="mismatch"):
            _read_hdf5_attr(h5_path, spec)

    def test_read_skips_validation_for_versionless_spec(self, tmp_path):
        """When spec.version is empty, validation is skipped."""
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
            dtype=np.float64,
            data_location="data/arr",
        )
        # spec.version defaults to "" — no validation
        h5_path = str(tmp_path / "nover.h5")
        with h5py.File(h5_path, "w") as f:
            f.create_dataset("data/arr", data=np.array([1.0]))
        result = _read_hdf5_dataset(h5_path, spec)
        assert result[0] == 1.0


# ---------------------------------------------------------------------------
# HDF5 version in writers
# ---------------------------------------------------------------------------


class TestWriteCreatesVersion:
    """Tests for version attribute creation by HDF5 writers."""

    def test_write_hdf5_dataset_creates_version_attr(self, tmp_path):
        """Writing a dataset to a new file creates the version attribute."""
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
            dtype=np.float64,
            data_location="data/arr",
        )
        object.__setattr__(spec, "version", "v2.0.0")
        h5_path = str(tmp_path / "new.h5")
        data = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        _write_hdf5_dataset(data, h5_path, spec)
        with h5py.File(h5_path, "r") as f:
            assert f.attrs[CONST.DATASPEC_VERSION_ATTR] == "v2.0.0"

    def test_write_hdf5_attr_creates_version_attr(self, tmp_path):
        """Writing attrs to a new file creates the version attribute."""
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_ATTR,
            dtype=dict,
            data_location="test_data",
        )
        object.__setattr__(spec, "version", "v2.0.0")
        h5_path = str(tmp_path / "new.h5")
        data = {"test_params": {"key1": 3.14}}
        _write_hdf5_attr(data, h5_path, spec)
        with h5py.File(h5_path, "r") as f:
            assert f.attrs[CONST.DATASPEC_VERSION_ATTR] == "v2.0.0"

    def test_write_to_wrong_version_file_raises(self, tmp_path):
        """Writing to a file with wrong version raises ValueError."""
        spec = DataSetSpec(
            dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
            dtype=np.float64,
            data_location="data/arr",
        )
        object.__setattr__(spec, "version", "v2.0.0")
        h5_path = str(tmp_path / "wrong.h5")
        with h5py.File(h5_path, "w") as f:
            f.attrs[CONST.DATASPEC_VERSION_ATTR] = "v1.0.0"
        data = np.array([1.0, 2.0], dtype=np.float64)
        with pytest.raises(ValueError, match="mismatch"):
            _write_hdf5_dataset(data, h5_path, spec)


# ---------------------------------------------------------------------------
# _read_file_parsed tests
# ---------------------------------------------------------------------------

# Minimal micro log content for _read_file_parsed tests.
# Uses a subset of parameters that parse_micro_log can handle.
_MICRO_LOG_FOR_FILEOPS = textwrap.dedent(
    """\
 seed=  2133256963
 nodes=          13
 KdtPAnoplg=  0.360000000000000
 simulations=       50000
 KdtPAyesplg=  2.000000000000000E-002
 KdPLGnicked=   2.20000000000000
 KdPLGintact=   38.0000000000000
  kncat=   5.00000000000000
  kapcat=  0.100000000000000
  ktPAon=  0.100000000000000
  kaoff10=  3.600000000000000E-002
  kaoff12=  2.000000000000000E-003
  kplioff=   57.6000000000000
  kplgoffnick=  0.220000000000000
  kplgon=  0.100000000000000
  kplgoff=   3.80000000000000
  freeplg=   2.00000000000000
  kdeg=   5.00000000000000
  stats=        1000
"""
)


def _make_parsed_spec(version="v1.95.0", collection="microscale_out"):
    """Build a FILE_PARSED DataSetSpec with the given version and collection name."""
    spec = DataSetSpec(
        dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_PARSED,
        dtype=Quantity,
        data_location="micro{file_code}.txt",
    )
    object.__setattr__(spec, "version", version)
    object.__setattr__(spec, "collection", collection)
    return spec


class TestReadFileParsed:
    """Tests for :func:`_read_file_parsed`."""

    def test_parses_micro_log_to_params_dict(self, tmp_path):
        """Parsed log returns dict with 'micro_params' key containing raw Fortran names."""
        (tmp_path / "micro.txt").write_text(_MICRO_LOG_FOR_FILEOPS)
        spec = _make_parsed_spec()
        result = _read_file_parsed(str(tmp_path), spec)

        assert "micro_params" in result
        micro = result["micro_params"]
        # Raw Fortran names are preserved (no resolution)
        assert micro["nodes"] == pytest.approx(13)
        assert micro["simulations"] == pytest.approx(50000)

    def test_values_are_floats(self, tmp_path):
        """Values are raw floats (no Quantity wrapping or string conversion)."""
        (tmp_path / "micro.txt").write_text(_MICRO_LOG_FOR_FILEOPS)
        spec = _make_parsed_spec()
        result = _read_file_parsed(str(tmp_path), spec)

        micro = result["micro_params"]
        assert "kdtpanoplg" in micro
        assert isinstance(micro["kdtpanoplg"], float)
        assert micro["kdtpanoplg"] == pytest.approx(0.36)

    def test_unsupported_collection_raises(self, tmp_path):
        """NotImplementedError raised for unknown collection."""
        spec = _make_parsed_spec(collection="unknown_collection")
        with pytest.raises(NotImplementedError, match="unknown_collection"):
            _read_file_parsed(str(tmp_path), spec)

    def test_v190_stop_pattern(self, tmp_path):
        """v1.90.0 uses init_entry stop pattern instead of stats."""
        content = textwrap.dedent(
            """\
         KdtPAnoplg=  0.360000000000000
          stats=           1
          ktPAon=  0.100000000000000
         init_entry=          13
         p_rebind=  2.616532229748891E-005
        """
        )
        (tmp_path / "micro.txt").write_text(content)
        spec = _make_parsed_spec(version="v1.90.0")
        result = _read_file_parsed(str(tmp_path), spec)

        micro = result["micro_params"]
        # ktPAon appears after stats, before init_entry — should be captured
        assert "ktpaon" in micro
        assert micro["ktpaon"] == pytest.approx(0.1)
        # stats should also be captured (raw parser doesn't filter it)
        assert "stats" in micro
        # init_entry and p_rebind are after the stop pattern
        assert "init_entry" not in micro


class TestReadDataCollectionFileParsed:
    """Integration test: read_data_collection with FILE_PARSED params."""

    def test_read_collection_file_parsed(self, tmp_path):
        """End-to-end: read v1.95.0 microscale_out with parsed params and binary data."""
        # Write micro log file
        (tmp_path / "micro.txt").write_text(_MICRO_LOG_FOR_FILEOPS)

        # Write a simple binary data file for one of the datasets (lysis)
        lysis_data = np.array([100.0, 200.0, 300.0], dtype=np.float64)
        lysis_data.tofile(tmp_path / "lysis.dat")

        # Use the real v1.95.0 microscale_out spec but only read 'lysis' dataset
        # to keep the test simple. Build a trimmed collection spec.
        real_spec = dataspec["v1.95.0"]["microscale_out"]
        trimmed = DataCollectionSpec(
            simulations_combined=True,
            params=real_spec.params,
            data={"lysis": real_spec.data["lysis"]},
        )
        object.__setattr__(trimmed, "version", "")
        object.__setattr__(trimmed, "collection", "microscale_out")

        data = read_data_collection(
            str(tmp_path),
            collections=[trimmed],
            file_codes=[""],
        )

        # Verify params were parsed — raw Fortran names
        assert "micro_params" in data["params"]
        assert data["params"]["micro_params"]["simulations"] == pytest.approx(50000)

        # Verify data was read
        np.testing.assert_array_equal(data["lysis"], lysis_data)


# ---------------------------------------------------------------------------
# Fortran log parsing tests (moved from test_paramcheck.py)
# ---------------------------------------------------------------------------


def _write_log(content, tmp_path, filename="test_log.txt"):
    """Write inline content to a temp file and return its Path."""
    path = tmp_path / filename
    path.write_text(content)
    return path


class TestParseMicroLog:
    """Tests for parse_micro_log(): Fortran micro log → raw KV dict."""

    def test_parses_kv_pairs(self, tmp_path):
        """Numeric KV pairs are extracted with lowercase Fortran names."""
        content = textwrap.dedent("""\
             KdtPAnoplg=  0.360000000000000
              stats=        1000
        """)
        path = _write_log(content, tmp_path)
        result = parse_micro_log(path)

        assert "kdtpanoplg" in result
        assert result["kdtpanoplg"] == pytest.approx(0.36)

    def test_stops_at_stats_line(self, tmp_path):
        """Lines at and after the first stats= line are not parsed."""
        content = textwrap.dedent("""\
             KdtPAnoplg=  0.360000000000000
              stats=        1000
              KdtPAyesplg=  0.999
        """)
        path = _write_log(content, tmp_path)
        result = parse_micro_log(path)

        assert "kdtpanoplg" in result
        # KdtPAyesplg appears AFTER stats=; must not be parsed
        assert "kdtpayesplg" not in result

    def test_unknown_keys_pass_through(self, tmp_path):
        """Unknown numeric keys pass through silently (no ValueError)."""
        content = " bogusparam=       50000\n stats=1\n"
        path = _write_log(content, tmp_path)
        result = parse_micro_log(path)

        assert "bogusparam" in result
        assert result["bogusparam"] == pytest.approx(50000)

    def test_setting_lines_parsed(self, tmp_path):
        """'Setting key = value' lines from command-line parsing are extracted."""
        content = textwrap.dedent("""\
             Setting nodes =           13
             Setting simulations =        50000
             Setting seed =   2133256963
             Setting outFileCode = _PLG2_tPA01_TB-xiii
              stats=        1000
        """)
        path = _write_log(content, tmp_path)
        result = parse_micro_log(path)

        assert result["nodes"] == pytest.approx(13)
        assert result["simulations"] == pytest.approx(50000)
        assert result["seed"] == pytest.approx(2133256963)
        # String-valued Setting lines (outFileCode) are silently skipped
        assert "outfilecode" not in result

    def test_setting_overwritten_by_later_kv(self, tmp_path):
        """A later key=value line overwrites an earlier Setting line."""
        content = textwrap.dedent("""\
             Setting nodes =           99
             nodes=          13
              stats=        1000
        """)
        path = _write_log(content, tmp_path)
        result = parse_micro_log(path)

        assert result["nodes"] == pytest.approx(13)

    def test_custom_stop_pattern(self, tmp_path):
        """Custom stop_pattern overrides the default stats= pattern."""
        content = textwrap.dedent("""\
             KdtPAnoplg=  0.360000000000000
              stats=           1
              ktPAon=  0.100000000000000
             init_entry=          13
        """)
        path = _write_log(content, tmp_path)
        result = parse_micro_log(path, stop_pattern=r"^\s*init_entry\s*=")

        # ktPAon appears after stats but before init_entry — should be captured
        assert "ktpaon" in result
        assert result["ktpaon"] == pytest.approx(0.1)
        # stats is also captured (it's before init_entry)
        assert "stats" in result
        # init_entry is at/after stop — not captured
        assert "init_entry" not in result

    def test_non_numeric_values_skipped(self, tmp_path):
        """Non-numeric values like filetype=binary are silently skipped."""
        content = " filetype=binary\n nodes=13\n stats=1\n"
        path = _write_log(content, tmp_path)
        result = parse_micro_log(path)

        assert "filetype" not in result
        assert "nodes" in result


class TestParseMacroLog:
    """Tests for parse_macro_log(): Fortran macro log → raw KV dict."""

    def test_parses_grid_params(self, tmp_path):
        """N=, F=, M= are extracted with lowercase names and float values."""
        content = textwrap.dedent("""\
             N=          19
              F=         184
              M=                 21105
            After     10. sec, 0 fibers are degraded
        """)
        path = _write_log(content, tmp_path)
        result = parse_macro_log(path)

        assert result["n"] == pytest.approx(19)
        assert result["f"] == pytest.approx(184)
        assert result["m"] == pytest.approx(21105)

    def test_stops_at_after_line(self, tmp_path):
        """Lines beginning with 'After ' are not parsed."""
        content = textwrap.dedent("""\
             N=          19
             After     10. sec, 0 fibers are degraded
             F=         999
        """)
        path = _write_log(content, tmp_path)
        result = parse_macro_log(path)

        assert result["n"] == pytest.approx(19)
        # F= appears after the After line; must not be parsed
        assert "f" not in result


class TestParseMicroFileCode:
    """Tests for parse_micro_file_code(): extracting params from file codes."""

    def test_q4_code(self):
        """Q4 maps to fiber_radius=72.7 nm and nodes_in_micro_row=13."""
        result = parse_micro_file_code("_PLG2_tPA01_Q4")
        assert "fiber_radius" in result
        assert result["fiber_radius"].to("nm").magnitude == pytest.approx(72.7)
        assert result["nodes_in_micro_row"] == 13

    def test_q2_code(self):
        """Q2 maps to fiber_radius=36.35 nm and nodes_in_micro_row=7."""
        result = parse_micro_file_code("_PLG2_tPA01_Q2")
        assert result["fiber_radius"].to("nm").magnitude == pytest.approx(36.35)
        assert result["nodes_in_micro_row"] == 7

    def test_tb_xiii_code(self):
        """TB-xiii maps to fiber_radius=61.5 nm and nodes_in_micro_row=13."""
        result = parse_micro_file_code("_PLG2_tPA01_TB-xiii")
        assert result["fiber_radius"].to("nm").magnitude == pytest.approx(61.5)
        assert result["nodes_in_micro_row"] == 13

    def test_no_fiber_code(self):
        """A file code with no recognized fiber type returns empty dict."""
        result = parse_micro_file_code("_TK-L_307")
        assert result == {}

    def test_empty_code(self):
        """An empty file code returns empty dict."""
        result = parse_micro_file_code("")
        assert result == {}
