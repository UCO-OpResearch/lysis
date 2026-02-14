"""Comprehensive pytest unit tests for the lysis.data.datastore module.

Tests cover:
- DataStatus enum values and flag operations
- h5_tree output formatting
- DataStore: collection detection, parameter loading, dot access, context manager
- DataCollection: combined vs per-simulation, dataset access, indexing
- SimulationView: per-sim dataset dot access
"""

import os
import warnings

import numpy as np
import h5py
import pytest

from lysis.data.datastore import (
    DataStore,
    DataCollection,
    SimulationView,
    DataStatus,
    h5_tree,
    _group_path_from_spec,
)
from lysis.data.dataspec import DataCollectionSpec, DataSetSpec, dataspec
from lysis.config.parameters import MicroParameters, MacroParameters


# ---------------------------------------------------------------------------
#  Helpers: create HDF5 files matching v2.0.0 structure
# ---------------------------------------------------------------------------


def _write_micro_attrs(h5file, **overrides):
    """Write MicroParameters-compatible attributes to micro_data group.

    Uses a minimal set of attributes. Override any with keyword args.
    """
    grp = h5file.require_group("micro_data")
    defaults = {
        "micro_simulations": 100,
        "fiber_radius": "36.35 nanometer",
    }
    defaults.update(overrides)
    for k, v in defaults.items():
        grp.attrs[k] = v


def _write_macro_attrs(h5file, n_sims=3, **overrides):
    """Write MacroParameters-compatible attributes to macro_data group.

    Uses a minimal set of attributes. Override any with keyword args.

    :param n_sims: Number of simulations (written as ``macro_simulations``).
    """
    grp = h5file.require_group("macro_data")
    defaults = {
        "rows": 5,
        "macro_simulations": n_sims,
    }
    defaults.update(overrides)
    for k, v in defaults.items():
        grp.attrs[k] = v


def _write_micro_datasets(h5file, n_sims=10):
    """Create microscale_out datasets in the HDF5 file."""
    spec = dataspec["v2.0.0"]["microscale_out"]
    for name, ds_spec in spec.data.items():
        if ds_spec.data_location is None:
            continue
        # Create parent groups as needed
        path = ds_spec.data_location
        if ds_spec.dtype == h5py.string_dtype():
            data = np.array(["log line"] * n_sims, dtype=object)
            h5file.create_dataset(path, data=data, dtype=h5py.string_dtype())
        elif ds_spec.dtype == np.bool:
            h5file.create_dataset(path, data=np.ones(n_sims, dtype=bool))
        elif np.issubdtype(ds_spec.dtype, np.integer):
            h5file.create_dataset(
                path, data=np.arange(n_sims, dtype=ds_spec.dtype)
            )
        else:
            h5file.create_dataset(
                path, data=np.random.rand(n_sims).astype(ds_spec.dtype)
            )


def _write_macro_datasets(h5file, n_sims=3, n_snapshots=5):
    """Create macroscale_out datasets in the HDF5 file."""
    spec = dataspec["v2.0.0"]["macroscale_out"]
    for sim in range(n_sims):
        for name, ds_spec in spec.data.items():
            if ds_spec.data_location is None:
                continue
            path = ds_spec.data_location.format(sim=sim)
            if ds_spec.dtype == h5py.string_dtype():
                data = np.array([f"sim {sim} log"] * 5, dtype=object)
                h5file.create_dataset(path, data=data, dtype=h5py.string_dtype())
            elif hasattr(ds_spec.dtype, "names") and ds_spec.dtype.names:
                # Structured dtype
                data = np.zeros(n_snapshots, dtype=ds_spec.dtype)
                h5file.create_dataset(path, data=data)
            elif ds_spec.dtype == np.int32:
                # tpa_location_snapshot: shape (n_tpa, 2, n_snapshots)
                data = np.zeros((2, 2, n_snapshots), dtype=np.int32)
                h5file.create_dataset(path, data=data)
            else:
                data = np.random.rand(n_snapshots).astype(ds_spec.dtype)
                h5file.create_dataset(path, data=data)


# ---------------------------------------------------------------------------
#  DataStatus enum tests
# ---------------------------------------------------------------------------


class TestDataStatus:
    """Tests for the DataStatus Flag enum."""

    def test_all_values_exist(self):
        """All expected members are present on the enum."""
        assert hasattr(DataStatus, "NONE")
        assert hasattr(DataStatus, "INITIALIZED")
        assert hasattr(DataStatus, "LOADED")
        assert hasattr(DataStatus, "SAVED")
        assert hasattr(DataStatus, "FILLED")

    def test_values_are_distinct(self):
        """Every member has a unique underlying value."""
        values = [
            DataStatus.NONE.value,
            DataStatus.INITIALIZED.value,
            DataStatus.LOADED.value,
            DataStatus.SAVED.value,
            DataStatus.FILLED.value,
        ]
        assert len(values) == len(set(values))

    def test_none_is_falsy(self):
        """DataStatus.NONE should evaluate as falsy."""
        assert not DataStatus.NONE

    def test_non_none_are_truthy(self):
        """All members other than NONE should evaluate as truthy."""
        assert DataStatus.INITIALIZED
        assert DataStatus.LOADED
        assert DataStatus.SAVED
        assert DataStatus.FILLED

    def test_flag_combination(self):
        """Flag members can be combined with the | operator."""
        combined = DataStatus.INITIALIZED | DataStatus.LOADED
        assert DataStatus.INITIALIZED in combined
        assert DataStatus.LOADED in combined
        assert DataStatus.SAVED not in combined

    def test_none_value_is_zero(self):
        """NONE must have the integer value 0."""
        assert DataStatus.NONE.value == 0


# ---------------------------------------------------------------------------
#  h5_tree tests
# ---------------------------------------------------------------------------


class TestH5Tree:
    """Tests for the h5_tree recursive tree-printing utility."""

    def test_group_with_dataset(self, tmp_path):
        """A group containing a dataset shows the dataset name and shape."""
        filepath = tmp_path / "tree1.h5"
        with h5py.File(filepath, "w") as f:
            grp = f.create_group("sensors")
            grp.create_dataset("temperature", data=np.zeros((10, 3)))

        with h5py.File(filepath, "r") as f:
            output = h5_tree(f)

        assert "sensors" in output
        assert "temperature" in output
        assert "(10, 3)" in output

    def test_scalar_dataset(self, tmp_path):
        """A scalar (0-d) dataset shows its empty-tuple shape in the output."""
        filepath = tmp_path / "tree_scalar.h5"
        with h5py.File(filepath, "w") as f:
            f.create_dataset("version", data=42)

        with h5py.File(filepath, "r") as f:
            output = h5_tree(f)

        assert "version" in output
        assert "()" in output

    def test_scalar_fallback_label(self):
        """When accessing .shape raises TypeError, the item is labelled (scalar)."""
        from unittest.mock import MagicMock

        mock_item = MagicMock(spec=[])
        type(mock_item).shape = property(
            lambda self: (_ for _ in ()).throw(TypeError)
        )

        mock_group = MagicMock()
        mock_group.__len__ = MagicMock(return_value=1)
        mock_group.items = MagicMock(return_value=[("bad_shape", mock_item)])

        output = h5_tree(mock_group)
        assert "bad_shape" in output
        assert "(scalar)" in output

    def test_nested_groups(self, tmp_path):
        """Multi-level nesting produces the correct tree drawing characters."""
        filepath = tmp_path / "tree_nested.h5"
        with h5py.File(filepath, "w") as f:
            g1 = f.create_group("level1")
            g2 = g1.create_group("level2")
            g2.create_dataset("values", data=np.arange(5))

        with h5py.File(filepath, "r") as f:
            output = h5_tree(f)

        assert "level1" in output
        assert "level2" in output
        assert "values" in output
        assert "└──" in output

    def test_multiple_items_tree_characters(self, tmp_path):
        """When a group has multiple children, non-last items use the tee character."""
        filepath = tmp_path / "tree_multi.h5"
        with h5py.File(filepath, "w") as f:
            f.create_dataset("alpha", data=np.array([1, 2]))
            f.create_dataset("beta", data=np.array([3, 4]))

        with h5py.File(filepath, "r") as f:
            output = h5_tree(f)

        assert "├──" in output
        assert "└──" in output

    def test_empty_file(self, tmp_path):
        """An empty HDF5 file produces an empty string."""
        filepath = tmp_path / "tree_empty.h5"
        with h5py.File(filepath, "w"):
            pass

        with h5py.File(filepath, "r") as f:
            output = h5_tree(f)

        assert output == ""


# ---------------------------------------------------------------------------
#  _group_path_from_spec tests
# ---------------------------------------------------------------------------


class TestGroupPathFromSpec:
    """Tests for the helper that extracts HDF5 group paths from collection specs."""

    def test_microscale_out_returns_micro_data(self):
        spec = dataspec["v2.0.0"]["microscale_out"]
        assert _group_path_from_spec(spec) == "micro_data"

    def test_macroscale_out_returns_macro_data(self):
        spec = dataspec["v2.0.0"]["macroscale_out"]
        assert _group_path_from_spec(spec) == "macro_data"

    def test_macroscale_in_returns_none(self):
        spec = dataspec["v2.0.0"]["macroscale_in"]
        assert _group_path_from_spec(spec) is None


# ---------------------------------------------------------------------------
#  DataStore tests
# ---------------------------------------------------------------------------


class TestDataStoreInit:
    """Tests for DataStore construction and collection detection."""

    def test_open_empty_hdf5(self, tmp_path):
        """Opening an empty HDF5 file produces no collections and no params."""
        filepath = tmp_path / "empty.h5"
        with h5py.File(filepath, "w"):
            pass

        with DataStore("empty", str(tmp_path)) as ds:
            assert ds.collections == {}
            assert ds.micro_params is None
            assert ds.macro_params is None

    def test_open_with_micro_data(self, tmp_path):
        """micro_data group detected as microscale_out collection."""
        filepath = tmp_path / "micro.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)

        with DataStore("micro", str(tmp_path)) as ds:
            assert "microscale_out" in ds.collections
            assert "macroscale_out" not in ds.collections

    def test_macro_without_micro_raises(self, tmp_path):
        """Opening a file with macro_data but no micro_data raises ValueError."""
        filepath = tmp_path / "macro.h5"
        with h5py.File(filepath, "w") as f:
            _write_macro_attrs(f, n_sims=2)
            _write_macro_datasets(f, n_sims=2)

        with pytest.raises(ValueError, match="microscale_out"):
            DataStore("macro", str(tmp_path))

    def test_open_with_both_collections(self, tmp_path):
        """Both collections detected when both groups present."""
        filepath = tmp_path / "both.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f, n_sims=2)
            _write_macro_datasets(f, n_sims=2)

        with DataStore("both", str(tmp_path)) as ds:
            assert "microscale_out" in ds.collections
            assert "macroscale_out" in ds.collections

    def test_macroscale_in_not_detected(self, tmp_path):
        """macroscale_in (derived collection) is never in collections."""
        filepath = tmp_path / "both.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f, n_sims=2)
            _write_macro_datasets(f, n_sims=2)

        with DataStore("both", str(tmp_path)) as ds:
            assert "macroscale_in" not in ds.collections

    def test_context_manager_closes_file(self, tmp_path):
        """The HDF5 file is closed after exiting the context manager."""
        filepath = tmp_path / "ctx.h5"
        with h5py.File(filepath, "w"):
            pass

        ds = DataStore("ctx", str(tmp_path))
        ds.__enter__()
        assert ds._file.id.valid
        ds.__exit__(None, None, None)
        assert not ds._file.id.valid

    def test_close_method(self, tmp_path):
        """Calling close() closes the HDF5 file."""
        filepath = tmp_path / "close.h5"
        with h5py.File(filepath, "w"):
            pass

        ds = DataStore("close", str(tmp_path))
        assert ds._file.id.valid
        ds.close()
        assert not ds._file.id.valid


class TestDataStoreDotAccess:
    """Tests for DataStore attribute access to collections."""

    def test_dot_access_to_microscale_out(self, tmp_path):
        """datastore.microscale_out returns a DataCollection."""
        filepath = tmp_path / "dot.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)

        with DataStore("dot", str(tmp_path)) as ds:
            coll = ds.microscale_out
            assert isinstance(coll, DataCollection)

    def test_dot_access_to_macroscale_out(self, tmp_path):
        """datastore.macroscale_out returns a DataCollection."""
        filepath = tmp_path / "dot2.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f, n_sims=2)
            _write_macro_datasets(f, n_sims=2)

        with DataStore("dot2", str(tmp_path)) as ds:
            coll = ds.macroscale_out
            assert isinstance(coll, DataCollection)

    def test_nonexistent_attribute_raises(self, tmp_path):
        """Accessing a nonexistent attribute raises AttributeError."""
        filepath = tmp_path / "noattr.h5"
        with h5py.File(filepath, "w"):
            pass

        with DataStore("noattr", str(tmp_path)) as ds:
            with pytest.raises(AttributeError, match="no attribute"):
                _ = ds.nonexistent_collection

    def test_repr(self, tmp_path):
        """__repr__ includes run code and collection names."""
        filepath = tmp_path / "repr.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)

        with DataStore("repr", str(tmp_path)) as ds:
            r = repr(ds)
            assert "repr" in r
            assert "microscale_out" in r

    def test_str_returns_tree(self, tmp_path):
        """__str__ returns h5_tree output."""
        filepath = tmp_path / "strt.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)

        with DataStore("strt", str(tmp_path)) as ds:
            s = str(ds)
            assert "micro_data" in s


# ---------------------------------------------------------------------------
#  Parameter loading tests
# ---------------------------------------------------------------------------


class TestParameterLoading:
    """Tests for loading MicroParameters and MacroParameters from HDF5 attrs."""

    def test_micro_params_loaded(self, tmp_path):
        """micro_params is a MicroParameters instance when attrs are present."""
        filepath = tmp_path / "params.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            with DataStore("params", str(tmp_path)) as ds:
                assert isinstance(ds.micro_params, MicroParameters)

    def test_micro_params_value_loaded(self, tmp_path):
        """micro_params attributes match what was written to HDF5."""
        filepath = tmp_path / "pval.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f, micro_simulations=42000)
            _write_micro_datasets(f)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            with DataStore("pval", str(tmp_path)) as ds:
                assert ds.micro_params.micro_simulations == 42000

    def test_micro_data_without_params_raises(self, tmp_path):
        """micro_data group without parameter attributes raises ValueError."""
        filepath = tmp_path / "noattr.h5"
        with h5py.File(filepath, "w") as f:
            f.create_group("micro_data")
            _write_micro_datasets(f)

        with pytest.raises(ValueError, match="no parameter attributes"):
            DataStore("noattr", str(tmp_path))

    def test_macro_params_loaded(self, tmp_path):
        """macro_params is a MacroParameters instance when attrs are present."""
        filepath = tmp_path / "mparams.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f, n_sims=2)
            _write_macro_datasets(f, n_sims=2)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            with DataStore("mparams", str(tmp_path)) as ds:
                assert isinstance(ds.macro_params, MacroParameters)

    def test_macro_params_has_micro_params(self, tmp_path):
        """macro_params.micro_params is the same MicroParameters loaded."""
        filepath = tmp_path / "mpmicro.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f, n_sims=2)
            _write_macro_datasets(f, n_sims=2)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            with DataStore("mpmicro", str(tmp_path)) as ds:
                assert ds.macro_params.micro_params is ds.micro_params

    def test_macro_params_value_loaded(self, tmp_path):
        """macro_params attributes match what was written to HDF5."""
        filepath = tmp_path / "mpval.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f, n_sims=2, rows=7)
            _write_macro_datasets(f, n_sims=2)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            with DataStore("mpval", str(tmp_path)) as ds:
                assert ds.macro_params.rows == 7

    def test_macro_data_without_params_raises(self, tmp_path):
        """macro_data group without parameter attributes raises ValueError."""
        filepath = tmp_path / "mnoattr.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            f.create_group("macro_data")
            _write_macro_datasets(f, n_sims=2)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            with pytest.raises(ValueError, match="no parameter attributes"):
                DataStore("mnoattr", str(tmp_path))

    def test_no_params_for_empty_file(self, tmp_path):
        """Both params are None when no groups exist."""
        filepath = tmp_path / "empty.h5"
        with h5py.File(filepath, "w"):
            pass

        with DataStore("empty", str(tmp_path)) as ds:
            assert ds.micro_params is None
            assert ds.macro_params is None


# ---------------------------------------------------------------------------
#  DataCollection (combined) tests
# ---------------------------------------------------------------------------


class TestDataCollectionCombined:
    """Tests for DataCollection with simulations_combined=True."""

    def test_dot_access_returns_h5py_dataset(self, tmp_path):
        """Accessing a dataset name returns an h5py.Dataset."""
        filepath = tmp_path / "comb.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f, n_sims=10)

        with DataStore("comb", str(tmp_path)) as ds:
            dataset = ds.microscale_out.pli_first_time
            assert isinstance(dataset, h5py.Dataset)

    def test_dataset_data_is_lazy(self, tmp_path):
        """The returned h5py.Dataset is not loaded until sliced."""
        filepath = tmp_path / "lazy.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f, n_sims=10)

        with DataStore("lazy", str(tmp_path)) as ds:
            dataset = ds.microscale_out.pli_first_time
            # It's an h5py.Dataset, not a numpy array
            assert not isinstance(dataset, np.ndarray)
            # Slicing loads the data
            arr = dataset[:]
            assert isinstance(arr, np.ndarray)

    def test_dataset_shape_matches(self, tmp_path):
        """Dataset shape matches what was written."""
        filepath = tmp_path / "shape.h5"
        n_sims = 15
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f, n_sims=n_sims)

        with DataStore("shape", str(tmp_path)) as ds:
            dataset = ds.microscale_out.pli_first_time
            assert dataset.shape == (n_sims,)

    def test_datasets_property(self, tmp_path):
        """The datasets property lists all available dataset names."""
        filepath = tmp_path / "dsets.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)

        with DataStore("dsets", str(tmp_path)) as ds:
            names = ds.microscale_out.datasets
            assert "pli_first_time" in names
            assert "tpa_leaving_time" in names
            assert "sim_final_time" in names

    def test_nonexistent_dataset_raises(self, tmp_path):
        """Accessing a nonexistent dataset raises AttributeError."""
        filepath = tmp_path / "noset.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)

        with DataStore("noset", str(tmp_path)) as ds:
            with pytest.raises(AttributeError, match="no dataset"):
                _ = ds.microscale_out.totally_fake_dataset

    def test_contains_known_dataset(self, tmp_path):
        """__contains__ returns True for known dataset names."""
        filepath = tmp_path / "cont.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)

        with DataStore("cont", str(tmp_path)) as ds:
            assert "pli_first_time" in ds.microscale_out
            assert "fake_dataset" not in ds.microscale_out

    def test_indexing_combined_raises_type_error(self, tmp_path):
        """Indexing a combined collection raises TypeError."""
        filepath = tmp_path / "noindex.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)

        with DataStore("noindex", str(tmp_path)) as ds:
            with pytest.raises(TypeError, match="combined simulations"):
                _ = ds.microscale_out[0]

    def test_len_combined_raises_type_error(self, tmp_path):
        """len() on a combined collection raises TypeError."""
        filepath = tmp_path / "nolen.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)

        with DataStore("nolen", str(tmp_path)) as ds:
            with pytest.raises(TypeError, match="combined simulations"):
                len(ds.microscale_out)

    def test_repr_combined(self, tmp_path):
        """__repr__ for combined collection shows 'combined'."""
        filepath = tmp_path / "repr.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)

        with DataStore("repr", str(tmp_path)) as ds:
            r = repr(ds.microscale_out)
            assert "combined" in r
            assert "microscale_out" in r


# ---------------------------------------------------------------------------
#  DataCollection (per-simulation) tests
# ---------------------------------------------------------------------------


class TestDataCollectionPerSim:
    """Tests for DataCollection with simulations_combined=False."""

    def test_getitem_returns_simulation_view(self, tmp_path):
        """collection[i] returns a SimulationView."""
        filepath = tmp_path / "persim.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f)
            _write_macro_datasets(f, n_sims=3)

        with DataStore("persim", str(tmp_path)) as ds:
            view = ds.macroscale_out[0]
            assert isinstance(view, SimulationView)

    def test_len_returns_simulation_count(self, tmp_path):
        """len(collection) returns the number of simulations."""
        filepath = tmp_path / "simlen.h5"
        n_sims = 5
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f, n_sims=n_sims)
            _write_macro_datasets(f, n_sims=n_sims)

        with DataStore("simlen", str(tmp_path)) as ds:
            assert len(ds.macroscale_out) == n_sims

    def test_index_out_of_range_raises(self, tmp_path):
        """Accessing sim beyond range raises IndexError."""
        filepath = tmp_path / "oor.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f)
            _write_macro_datasets(f, n_sims=3)

        with DataStore("oor", str(tmp_path)) as ds:
            with pytest.raises(IndexError, match="out of range"):
                _ = ds.macroscale_out[5]

    def test_negative_index_raises(self, tmp_path):
        """Negative simulation index raises IndexError."""
        filepath = tmp_path / "neg.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f)
            _write_macro_datasets(f, n_sims=3)

        with DataStore("neg", str(tmp_path)) as ds:
            with pytest.raises(IndexError, match="out of range"):
                _ = ds.macroscale_out[-1]

    def test_non_int_index_raises_type_error(self, tmp_path):
        """Non-integer index raises TypeError."""
        filepath = tmp_path / "noint.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f)
            _write_macro_datasets(f, n_sims=3)

        with DataStore("noint", str(tmp_path)) as ds:
            with pytest.raises(TypeError, match="integer"):
                _ = ds.macroscale_out["sim_00"]

    def test_direct_dataset_access_raises_type_error(self, tmp_path):
        """Dot-access to dataset on per-sim collection raises TypeError."""
        filepath = tmp_path / "directds.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f)
            _write_macro_datasets(f, n_sims=3)

        with DataStore("directds", str(tmp_path)) as ds:
            with pytest.raises(TypeError, match="per-simulation"):
                _ = ds.macroscale_out.snapshot_time

    def test_datasets_property(self, tmp_path):
        """The datasets property lists available dataset names."""
        filepath = tmp_path / "persimds.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f, n_sims=2)
            _write_macro_datasets(f, n_sims=2)

        with DataStore("persimds", str(tmp_path)) as ds:
            names = ds.macroscale_out.datasets
            assert "snapshot_time" in names
            assert "fiber_degrade_time" in names

    def test_repr_per_sim(self, tmp_path):
        """__repr__ shows 'per-simulation' and simulation count."""
        filepath = tmp_path / "persimrepr.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f, n_sims=4)
            _write_macro_datasets(f, n_sims=4)

        with DataStore("persimrepr", str(tmp_path)) as ds:
            r = repr(ds.macroscale_out)
            assert "per-simulation" in r
            assert "simulations=4" in r


# ---------------------------------------------------------------------------
#  SimulationView tests
# ---------------------------------------------------------------------------


class TestSimulationView:
    """Tests for SimulationView per-sim dataset access."""

    def test_dot_access_returns_h5py_dataset(self, tmp_path):
        """Accessing a dataset on SimulationView returns h5py.Dataset."""
        filepath = tmp_path / "simview.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f)
            _write_macro_datasets(f, n_sims=3)

        with DataStore("simview", str(tmp_path)) as ds:
            dataset = ds.macroscale_out[1].snapshot_time
            assert isinstance(dataset, h5py.Dataset)

    def test_dataset_from_correct_simulation(self, tmp_path):
        """Each simulation's dataset points to the right HDF5 path."""
        filepath = tmp_path / "correctsim.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f)
            # Write with distinctive data per sim
            for sim in range(3):
                path = "macro_data/sim_{sim:02}/snapshot_time".format(sim=sim)
                data = np.full(5, fill_value=float(sim * 100), dtype=np.float64)
                f.create_dataset(path, data=data)

        with DataStore("correctsim", str(tmp_path)) as ds:
            for sim in range(3):
                arr = ds.macroscale_out[sim].snapshot_time[:]
                assert np.all(arr == float(sim * 100))

    def test_datasets_property(self, tmp_path):
        """SimulationView.datasets lists available dataset names."""
        filepath = tmp_path / "svdsets.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f, n_sims=2)
            _write_macro_datasets(f, n_sims=2)

        with DataStore("svdsets", str(tmp_path)) as ds:
            names = ds.macroscale_out[0].datasets
            assert "snapshot_time" in names
            assert "fiber_degrade_time" in names

    def test_nonexistent_dataset_raises(self, tmp_path):
        """Accessing a nonexistent dataset on SimulationView raises AttributeError."""
        filepath = tmp_path / "svnoset.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f, n_sims=2)
            _write_macro_datasets(f, n_sims=2)

        with DataStore("svnoset", str(tmp_path)) as ds:
            with pytest.raises(AttributeError, match="no dataset"):
                _ = ds.macroscale_out[0].totally_fake

    def test_repr(self, tmp_path):
        """SimulationView __repr__ shows sim index."""
        filepath = tmp_path / "svrepr.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f)
            _write_macro_datasets(f, n_sims=3)

        with DataStore("svrepr", str(tmp_path)) as ds:
            r = repr(ds.macroscale_out[2])
            assert "sim=2" in r

    def test_full_access_chain(self, tmp_path):
        """Full access chain: datastore.macroscale_out[i].dataset works end-to-end."""
        filepath = tmp_path / "chain.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f)
            _write_macro_datasets(f, n_sims=2, n_snapshots=7)

        with DataStore("chain", str(tmp_path)) as ds:
            arr = ds.macroscale_out[1].snapshot_time[:]
            assert isinstance(arr, np.ndarray)
            assert arr.shape == (7,)
