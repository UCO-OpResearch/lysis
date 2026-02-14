"""Comprehensive pytest unit tests for :mod:`lysis.data.fileops`.

Tests cover:

* Low-level readers (``_read_file_text``, ``_read_file_binary``,
  ``_read_file_json``, ``_read_hdf5_dataset``, ``_read_hdf5_attr``)
* Low-level writers (``_write_file_text``, ``_write_file_binary``,
  ``_write_file_json``, ``_write_hdf5_dataset``, ``_write_hdf5_attr``)
* Mid-level dispatch (``read_dataset``, ``write_dataset``)
* High-level collection I/O (``read_data_collection``,
  ``write_data_collection``)
"""

import json

import h5py
import numpy as np
import pytest

from lysis.config.constants import CONST
from lysis.data.dataspec import (
    DataCollectionSpec,
    DataSetSpec,
    check_dataset_spec,
    parse_shape,
)
from lysis.data.fileops import (
    _read_file_binary,
    _read_file_json,
    _read_file_text,
    _read_hdf5_attr,
    _read_hdf5_dataset,
    _write_file_binary,
    _write_file_json,
    _write_file_text,
    _write_hdf5_attr,
    _write_hdf5_dataset,
    data_readers,
    data_writers,
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
        dt = np.dtype(
            [("time", np.float64), ("idx", np.int32), ("val", np.float64)]
        )
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
            np.column_stack(
                [rows["time"], rows["idx"], rows["val"]]
            ),
            delimiter=",",
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
        result = _read_file_binary(
            str(tmp_path), spec, params=sample_params
        )
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

    def test_combined_single_array(
        self, tmp_path, combined_collection_spec, json_spec
    ):
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

    def test_per_sim_list_of_arrays(
        self, tmp_path, per_sim_collection_spec
    ):
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

    def test_per_sim_graceful_stop(
        self, tmp_path, per_sim_collection_spec
    ):
        """Reading stops gracefully when sim file does not exist (after sim 0)."""
        self._setup_per_sim_files(tmp_path, n_sims=3)
        # Only 3 directories exist (00, 01, 02). Reading stops at 03.
        data = read_data_collection(
            str(tmp_path),
            collections=[per_sim_collection_spec],
            file_codes=[""],
        )
        assert len(data["arr"]) == 3

    def test_per_sim_sim0_missing_raises(
        self, tmp_path, per_sim_collection_spec
    ):
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

    def test_combined_round_trip(
        self, tmp_path, combined_collection_spec
    ):
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
        np.testing.assert_array_almost_equal(
            loaded["arr"], original_data["arr"]
        )
        assert (
            loaded["params"]["micro_params"]["micro_simulations"] == 10
        )

    def test_per_sim_round_trip(self, tmp_path, per_sim_collection_spec):
        """Write 3 per-sim arrays then read back; all 3 match."""
        arrays = [
            np.arange(5, dtype=np.float64) + i * 100.0
            for i in range(3)
        ]
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
