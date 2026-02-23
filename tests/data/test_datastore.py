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
    COMPATIBLE_DATASPEC_VERSION,
    DataStore,
    DataCollection,
    DerivedDataCollection,
    SimulationView,
    DataStatus,
    h5_tree,
    _group_path_from_spec,
)
from lysis.config.constants import CONST
from lysis.data.dataspec import DataCollectionSpec, DataSetSpec, dataspec
from lysis.config.parameters import MicroParameters, MacroParameters


# ---------------------------------------------------------------------------
#  Helpers: create HDF5 files matching v2.0.0 structure
# ---------------------------------------------------------------------------


def _write_version_attr(h5file, version="v2.0.0"):
    """Write the dataspec_version root attribute to an HDF5 file."""
    h5file.attrs[CONST.DATASPEC_VERSION_ATTR] = version


def _write_micro_attrs(h5file, **overrides):
    """Write MicroParameters-compatible attributes to micro_data group.

    Uses a minimal set of attributes. Override any with keyword args.
    Also writes the dataspec_version root attribute if not already present.
    """
    _write_version_attr(h5file)
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
        with h5py.File(filepath, "w") as f:
            _write_version_attr(f)

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
            _write_version_attr(f)
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

    def test_macroscale_in_detected_when_conditions_met(self, tmp_path):
        """macroscale_in (derived collection) is advertised when microscale_out
        has non-empty required datasets and macro_params is present."""
        filepath = tmp_path / "both.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)
            _write_macro_attrs(f, n_sims=2)
            _write_macro_datasets(f, n_sims=2)

        with DataStore("both", str(tmp_path)) as ds:
            assert "macroscale_in" in ds.collections
            assert isinstance(ds.collections["macroscale_in"], DerivedDataCollection)

    def test_macroscale_in_absent_without_macro_params(self, tmp_path):
        """macroscale_in is NOT present when only microscale_out exists
        (no macro_params)."""
        filepath = tmp_path / "micro_only.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f)
            _write_micro_datasets(f)

        with DataStore("micro_only", str(tmp_path)) as ds:
            assert "macroscale_in" not in ds.collections

    def test_context_manager_closes_file(self, tmp_path):
        """The HDF5 file is closed after exiting the context manager."""
        filepath = tmp_path / "ctx.h5"
        with h5py.File(filepath, "w") as f:
            _write_version_attr(f)

        ds = DataStore("ctx", str(tmp_path))
        ds.__enter__()
        assert ds._file.id.valid
        ds.__exit__(None, None, None)
        assert not ds._file.id.valid

    def test_close_method(self, tmp_path):
        """Calling close() closes the HDF5 file."""
        filepath = tmp_path / "close.h5"
        with h5py.File(filepath, "w") as f:
            _write_version_attr(f)

        ds = DataStore("close", str(tmp_path))
        assert ds._file.id.valid
        ds.close()
        assert not ds._file.id.valid

    def test_default_mode_is_read(self, tmp_path):
        """Default mode is 'r' (read-only)."""
        filepath = tmp_path / "defmode.h5"
        with h5py.File(filepath, "w") as f:
            _write_version_attr(f)

        with DataStore("defmode", str(tmp_path)) as ds:
            assert ds.mode == "r"

    def test_open_append_mode(self, tmp_path):
        """Opening with mode='a' sets mode property to 'a'."""
        filepath = tmp_path / "append.h5"
        with h5py.File(filepath, "w") as f:
            _write_version_attr(f)

        with DataStore("append", str(tmp_path), mode="a") as ds:
            assert ds.mode == "a"
            assert ds._file.id.valid


class TestDataStoreVersioning:
    """Tests for dataspec_version attribute validation on DataStore."""

    def test_missing_version_raises(self, tmp_path):
        """HDF5 file with no dataspec_version attribute raises ValueError."""
        filepath = tmp_path / "nover.h5"
        with h5py.File(filepath, "w"):
            pass  # No version attr

        with pytest.raises(ValueError, match=CONST.DATASPEC_VERSION_ATTR):
            DataStore("nover", str(tmp_path))

    def test_wrong_version_raises(self, tmp_path):
        """HDF5 file with wrong version raises ValueError."""
        filepath = tmp_path / "wrongver.h5"
        with h5py.File(filepath, "w") as f:
            f.attrs[CONST.DATASPEC_VERSION_ATTR] = "v1.0.0"

        with pytest.raises(ValueError, match="mismatch"):
            DataStore("wrongver", str(tmp_path))

    def test_correct_version_accepted(self, tmp_path):
        """HDF5 file with correct version opens successfully."""
        filepath = tmp_path / "goodver.h5"
        with h5py.File(filepath, "w") as f:
            _write_version_attr(f)

        with DataStore("goodver", str(tmp_path)) as ds:
            assert ds.collections == {}

    def test_dataspec_version_property(self, tmp_path):
        """dataspec_version property returns the version from the file."""
        filepath = tmp_path / "verprop.h5"
        with h5py.File(filepath, "w") as f:
            _write_version_attr(f)

        with DataStore("verprop", str(tmp_path)) as ds:
            assert ds.dataspec_version == "v2.0.0"


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
        with h5py.File(filepath, "w") as f:
            _write_version_attr(f)

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
            _write_version_attr(f)
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
            _write_version_attr(f)
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
        with h5py.File(filepath, "w") as f:
            _write_version_attr(f)

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


# ---------------------------------------------------------------------------
#  DataStore.create() tests
# ---------------------------------------------------------------------------


class TestDataStoreCreate:
    """Tests for the DataStore.create() factory class method."""

    def test_create_produces_valid_file(self, tmp_path):
        """create() produces an HDF5 file with the dataspec_version attribute."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()

        ds = DataStore.create("test_run", str(tmp_path), micro)
        try:
            filepath = tmp_path / "test_run.h5"
            assert filepath.exists()
            with h5py.File(filepath, "r") as f:
                assert f.attrs[CONST.DATASPEC_VERSION_ATTR] == COMPATIBLE_DATASPEC_VERSION
        finally:
            ds.close()

    def test_create_has_micro_params(self, tmp_path):
        """create() stores MicroParameters that survive round-trip."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters(micro_simulations=42)

        ds = DataStore.create("micro_rt", str(tmp_path), micro)
        try:
            assert isinstance(ds.micro_params, MicroParameters)
            assert ds.micro_params.micro_simulations == 42
        finally:
            ds.close()

    def test_create_has_microscale_out(self, tmp_path):
        """create() produces a DataStore with microscale_out collection."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()

        ds = DataStore.create("coll", str(tmp_path), micro)
        try:
            assert "microscale_out" in ds.collections
            assert isinstance(ds.microscale_out, DataCollection)
        finally:
            ds.close()

    def test_create_empty_datasets_exist(self, tmp_path):
        """create() creates zero-length datasets for all microscale_out specs."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()

        ds = DataStore.create("empty_ds", str(tmp_path), micro)
        try:
            spec = dataspec["v2.0.0"]["microscale_out"]
            for name, ds_spec in spec.data.items():
                if ds_spec.data_location is None:
                    continue
                # Dataset should be accessible via DataCollection
                dataset = getattr(ds.microscale_out, name)
                assert isinstance(dataset, h5py.Dataset)
                # First dimension should be 0 (empty)
                assert dataset.shape[0] == 0
        finally:
            ds.close()

    def test_create_datasets_have_correct_dtype(self, tmp_path):
        """Empty datasets have the correct dtype from the spec."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()

        ds = DataStore.create("dtype_check", str(tmp_path), micro)
        try:
            # Check a float64 dataset
            dataset = ds.microscale_out.pli_first_time
            assert dataset.dtype == np.float64
            # Check an integer dataset
            dataset = ds.microscale_out.tpa_final_num
            assert dataset.dtype == np.uint8
        finally:
            ds.close()

    def test_create_datasets_are_resizable(self, tmp_path):
        """Empty datasets have maxshape=None for variable dimensions."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()

        ds = DataStore.create("resize", str(tmp_path), micro)
        try:
            dataset = ds.microscale_out.pli_first_time
            # maxshape should have None for variable dims
            assert dataset.maxshape == (None,)
        finally:
            ds.close()

    def test_create_no_macroscale(self, tmp_path):
        """create() produces no macroscale_out collection."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()

        ds = DataStore.create("no_macro", str(tmp_path), micro)
        try:
            assert "macroscale_out" not in ds.collections
            assert ds.macro_params is None
        finally:
            ds.close()

    def test_create_status_initialized(self, tmp_path):
        """create() sets status to INITIALIZED."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()

        ds = DataStore.create("status", str(tmp_path), micro)
        try:
            assert ds.status == DataStatus.INITIALIZED
        finally:
            ds.close()

    def test_create_file_exists_error(self, tmp_path):
        """create() raises FileExistsError if file already exists."""
        filepath = tmp_path / "exists.h5"
        with h5py.File(filepath, "w") as f:
            _write_version_attr(f)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()

        with pytest.raises(FileExistsError, match="already exists"):
            DataStore.create("exists", str(tmp_path), micro)

    def test_create_type_error(self, tmp_path):
        """create() raises TypeError for non-MicroParameters."""
        with pytest.raises(TypeError, match="MicroParameters"):
            DataStore.create("bad", str(tmp_path), {"not": "params"})

    def test_create_context_manager(self, tmp_path):
        """DataStore from create() works with context manager."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()

        with DataStore.create("ctx", str(tmp_path), micro) as ds:
            assert "microscale_out" in ds.collections
        assert not ds._file.id.valid

    def test_create_mode_is_append(self, tmp_path):
        """DataStore from create() is opened in 'a' (read/write) mode."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()

        with DataStore.create("mode", str(tmp_path), micro) as ds:
            assert ds.mode == "a"

    def test_create_dataspec_version(self, tmp_path):
        """DataStore from create() reports correct dataspec_version."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()

        with DataStore.create("ver", str(tmp_path), micro) as ds:
            assert ds.dataspec_version == COMPATIBLE_DATASPEC_VERSION


# ---------------------------------------------------------------------------
#  DataStore.initialize_macroscale() tests
# ---------------------------------------------------------------------------


class TestDataStoreInitializeMacroscale:
    """Tests for the DataStore.initialize_macroscale() instance method."""

    def _create_micro_store(self, tmp_path, run_code="init_macro"):
        """Helper: create a microscale-only DataStore."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()
        return DataStore.create(run_code, str(tmp_path), micro)

    def test_initialize_has_both_collections(self, tmp_path):
        """initialize_macroscale() produces DataStore with both collections."""
        ds = self._create_micro_store(tmp_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            macro = MacroParameters(micro_params=ds.micro_params, macro_simulations=2)
        ds.initialize_macroscale(macro)
        try:
            assert "microscale_out" in ds.collections
            assert "macroscale_out" in ds.collections
        finally:
            ds.close()

    def test_initialize_has_macro_params(self, tmp_path):
        """initialize_macroscale() stores MacroParameters that survive round-trip."""
        ds = self._create_micro_store(tmp_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            macro = MacroParameters(
                micro_params=ds.micro_params, macro_simulations=3, rows=7
            )
        ds.initialize_macroscale(macro)
        try:
            assert isinstance(ds.macro_params, MacroParameters)
            assert ds.macro_params.rows == 7
            assert ds.macro_params.macro_simulations == 3
        finally:
            ds.close()

    def test_initialize_macro_params_has_micro_params(self, tmp_path):
        """macro_params.micro_params links to the loaded MicroParameters."""
        ds = self._create_micro_store(tmp_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            macro = MacroParameters(micro_params=ds.micro_params, macro_simulations=2)
        ds.initialize_macroscale(macro)
        try:
            assert ds.macro_params.micro_params is ds.micro_params
        finally:
            ds.close()

    def test_initialize_sim_groups_exist(self, tmp_path):
        """initialize_macroscale() creates per-simulation groups in HDF5."""
        n_sims = 3
        ds = self._create_micro_store(tmp_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            macro = MacroParameters(
                micro_params=ds.micro_params, macro_simulations=n_sims
            )
        ds.initialize_macroscale(macro)
        try:
            filepath = tmp_path / "init_macro.h5"
            with h5py.File(filepath, "r") as f:
                for sim in range(n_sims):
                    assert f"macro_data/sim_{sim:02}" in f
        finally:
            ds.close()

    def test_initialize_empty_datasets_exist(self, tmp_path):
        """initialize_macroscale() creates zero-length per-sim datasets."""
        n_sims = 2
        ds = self._create_micro_store(tmp_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            macro = MacroParameters(
                micro_params=ds.micro_params, macro_simulations=n_sims
            )
        ds.initialize_macroscale(macro)
        try:
            spec = dataspec["v2.0.0"]["macroscale_out"]
            for sim in range(n_sims):
                view = ds.macroscale_out[sim]
                for name, ds_spec in spec.data.items():
                    if ds_spec.data_location is None:
                        continue
                    dataset = getattr(view, name)
                    assert isinstance(dataset, h5py.Dataset)
                    assert dataset.shape[0] == 0
        finally:
            ds.close()

    def test_initialize_per_sim_datasets_resizable(self, tmp_path):
        """Per-sim empty datasets have maxshape=None for variable dims."""
        ds = self._create_micro_store(tmp_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            macro = MacroParameters(
                micro_params=ds.micro_params, macro_simulations=2
            )
        ds.initialize_macroscale(macro)
        try:
            dataset = ds.macroscale_out[0].snapshot_time
            assert dataset.maxshape == (None,)
        finally:
            ds.close()

    def test_initialize_len_matches(self, tmp_path):
        """len(macroscale_out) matches macro_simulations."""
        n_sims = 4
        ds = self._create_micro_store(tmp_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            macro = MacroParameters(
                micro_params=ds.micro_params, macro_simulations=n_sims
            )
        ds.initialize_macroscale(macro)
        try:
            assert len(ds.macroscale_out) == n_sims
        finally:
            ds.close()

    def test_initialize_status_initialized(self, tmp_path):
        """initialize_macroscale() sets status to INITIALIZED."""
        ds = self._create_micro_store(tmp_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            macro = MacroParameters(
                micro_params=ds.micro_params, macro_simulations=2
            )
        ds.initialize_macroscale(macro)
        try:
            assert ds.status == DataStatus.INITIALIZED
        finally:
            ds.close()

    def test_initialize_returns_none(self, tmp_path):
        """initialize_macroscale() returns None (modifies in place)."""
        ds = self._create_micro_store(tmp_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            macro = MacroParameters(
                micro_params=ds.micro_params, macro_simulations=2
            )
        result = ds.initialize_macroscale(macro)
        try:
            assert result is None
        finally:
            ds.close()

    def test_initialize_no_micro_raises(self, tmp_path):
        """initialize_macroscale() raises ValueError if no microscale_out."""
        filepath = tmp_path / "empty.h5"
        with h5py.File(filepath, "w") as f:
            _write_version_attr(f)

        ds = DataStore("empty", str(tmp_path), mode="a")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()
            macro = MacroParameters(micro_params=micro, macro_simulations=2)
        try:
            with pytest.raises(ValueError, match="microscale_out is not present"):
                ds.initialize_macroscale(macro)
        finally:
            ds.close()

    def test_initialize_already_macro_raises(self, tmp_path):
        """initialize_macroscale() raises ValueError if macroscale_out exists."""
        ds = self._create_micro_store(tmp_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            macro = MacroParameters(
                micro_params=ds.micro_params, macro_simulations=2
            )
        ds.initialize_macroscale(macro)
        try:
            with pytest.raises(ValueError, match="already present"):
                ds.initialize_macroscale(macro)
        finally:
            ds.close()

    def test_initialize_type_error(self, tmp_path):
        """initialize_macroscale() raises TypeError for non-MacroParameters."""
        ds = self._create_micro_store(tmp_path)
        try:
            with pytest.raises(TypeError, match="MacroParameters"):
                ds.initialize_macroscale({"not": "params"})
        finally:
            ds.close()

    def test_initialize_preserves_micro_data(self, tmp_path):
        """initialize_macroscale() preserves existing microscale empty datasets."""
        ds = self._create_micro_store(tmp_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            macro = MacroParameters(
                micro_params=ds.micro_params, macro_simulations=2
            )
        ds.initialize_macroscale(macro)
        try:
            # Microscale datasets should still exist and be empty
            dataset = ds.microscale_out.pli_first_time
            assert isinstance(dataset, h5py.Dataset)
            assert dataset.shape[0] == 0
        finally:
            ds.close()

    def test_initialize_in_context_manager(self, tmp_path):
        """initialize_macroscale() works inside a context manager."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()
        with DataStore.create("ctx_init", str(tmp_path), micro) as ds:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                macro = MacroParameters(
                    micro_params=ds.micro_params, macro_simulations=2
                )
            ds.initialize_macroscale(macro)
            assert "macroscale_out" in ds.collections
        assert not ds._file.id.valid

    def test_initialize_read_only_raises(self, tmp_path):
        """initialize_macroscale() raises IOError on a read-only DataStore."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()
        # Create the file first
        ds = DataStore.create("ro_init", str(tmp_path), micro)
        ds.close()
        # Re-open read-only
        ds = DataStore("ro_init", str(tmp_path))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            macro = MacroParameters(micro_params=ds.micro_params, macro_simulations=2)
        try:
            with pytest.raises(IOError, match="read-only"):
                ds.initialize_macroscale(macro)
        finally:
            ds.close()

    def test_initialize_preserves_mode(self, tmp_path):
        """initialize_macroscale() preserves the DataStore's mode."""
        ds = self._create_micro_store(tmp_path)
        assert ds.mode == "a"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            macro = MacroParameters(
                micro_params=ds.micro_params, macro_simulations=2
            )
        ds.initialize_macroscale(macro)
        try:
            assert ds.mode == "a"
        finally:
            ds.close()


# ---------------------------------------------------------------------------
#  Dataset write-through tests (append mode)
# ---------------------------------------------------------------------------


class TestDataStoreWriteThrough:
    """Tests that HDF5 datasets are writable through DataStore in 'a' mode.

    DataStore.create() returns a store in ``"a"`` mode.  The datasets
    exposed via DataCollection and SimulationView are live ``h5py.Dataset``
    objects, so external code should be able to resize them and write data
    directly.
    """

    # -- helpers --

    @staticmethod
    def _create_micro_store(tmp_path, run_code="write_test"):
        """Create a microscale-only DataStore in append mode."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()
        return DataStore.create(run_code, str(tmp_path), micro)

    @staticmethod
    def _create_full_store(tmp_path, run_code="write_full", n_sims=2):
        """Create a DataStore with both microscale and macroscale."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()
        ds = DataStore.create(run_code, str(tmp_path), micro)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            macro = MacroParameters(
                micro_params=ds.micro_params, macro_simulations=n_sims
            )
        ds.initialize_macroscale(macro)
        return ds

    # -- combined collection (microscale_out) --

    def test_resize_and_write_combined_float64(self, tmp_path):
        """Resize a combined float64 dataset and write data through DataStore."""
        with self._create_micro_store(tmp_path) as ds:
            dataset = ds.microscale_out.pli_first_time
            assert dataset.shape == (0,)

            dataset.resize((5,))
            data = np.array([1.1, 2.2, 3.3, 4.4, 5.5], dtype=np.float64)
            dataset[:] = data

            assert dataset.shape == (5,)
            np.testing.assert_array_equal(dataset[:], data)

    def test_resize_and_write_combined_uint8(self, tmp_path):
        """Resize a combined uint8 dataset and write data through DataStore."""
        with self._create_micro_store(tmp_path) as ds:
            dataset = ds.microscale_out.tpa_final_num
            assert dataset.shape == (0,)

            dataset.resize((3,))
            data = np.array([10, 20, 30], dtype=np.uint8)
            dataset[:] = data

            np.testing.assert_array_equal(dataset[:], data)

    def test_resize_and_write_combined_bool(self, tmp_path):
        """Resize a combined bool dataset and write data through DataStore."""
        with self._create_micro_store(tmp_path) as ds:
            dataset = ds.microscale_out.fiber_degraded
            assert dataset.shape == (0,)

            dataset.resize((4,))
            data = np.array([True, False, True, False])
            dataset[:] = data

            np.testing.assert_array_equal(dataset[:], data)

    def test_incremental_resize_combined(self, tmp_path):
        """Datasets can be resized incrementally, appending rows."""
        with self._create_micro_store(tmp_path) as ds:
            dataset = ds.microscale_out.sim_final_time

            # First batch
            dataset.resize((3,))
            dataset[:3] = [1.0, 2.0, 3.0]

            # Append more
            dataset.resize((6,))
            dataset[3:6] = [4.0, 5.0, 6.0]

            expected = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
            np.testing.assert_array_equal(dataset[:], expected)

    def test_slice_write_combined(self, tmp_path):
        """Writing to a slice of a combined dataset works."""
        with self._create_micro_store(tmp_path) as ds:
            dataset = ds.microscale_out.pli_first_time
            dataset.resize((5,))
            dataset[:] = np.zeros(5, dtype=np.float64)

            # Overwrite middle slice
            dataset[1:4] = [10.0, 20.0, 30.0]

            expected = np.array([0.0, 10.0, 20.0, 30.0, 0.0])
            np.testing.assert_array_equal(dataset[:], expected)

    def test_write_persists_after_close_and_reopen(self, tmp_path):
        """Data written through DataStore persists after close and reopen."""
        data = np.array([7.7, 8.8, 9.9], dtype=np.float64)

        ds = self._create_micro_store(tmp_path)
        ds.microscale_out.pli_first_time.resize((3,))
        ds.microscale_out.pli_first_time[:] = data
        ds.close()

        # Reopen read-only and verify
        with DataStore("write_test", str(tmp_path)) as ds2:
            np.testing.assert_array_equal(
                ds2.microscale_out.pli_first_time[:], data
            )

    # -- per-simulation collection (macroscale_out) --

    def test_resize_and_write_per_sim_float64(self, tmp_path):
        """Resize a per-sim float64 dataset and write data through DataStore."""
        with self._create_full_store(tmp_path) as ds:
            dataset = ds.macroscale_out[0].snapshot_time
            assert dataset.shape == (0,)

            dataset.resize((4,))
            data = np.array([0.1, 0.5, 1.0, 2.0], dtype=np.float64)
            dataset[:] = data

            np.testing.assert_array_equal(dataset[:], data)

    def test_write_per_sim_structured_dtype(self, tmp_path):
        """Write to a structured-dtype per-sim dataset through DataStore."""
        with self._create_full_store(tmp_path) as ds:
            dataset = ds.macroscale_out[0].fiber_degrade_time
            assert dataset.shape == (0,)

            record = np.array(
                [(1.5, 2, 3, 4.5), (6.0, 0, 1, 7.0)],
                dtype=np.dtype(
                    [
                        ("Simulation Time Elapsed", np.float64),
                        ("Grid Location Row", np.uint32),
                        ("Grid Location Rank", np.uint32),
                        ("Fiber New Degrade Time", np.float64),
                    ]
                ),
            )
            dataset.resize((2,))
            dataset[:] = record

            result = dataset[:]
            np.testing.assert_array_equal(
                result["Simulation Time Elapsed"], [1.5, 6.0]
            )
            np.testing.assert_array_equal(
                result["Grid Location Row"], [2, 0]
            )

    def test_write_different_simulations_independently(self, tmp_path):
        """Each simulation's datasets are independent; writing one doesn't affect another."""
        with self._create_full_store(tmp_path, n_sims=2) as ds:
            ds0 = ds.macroscale_out[0].snapshot_time
            ds1 = ds.macroscale_out[1].snapshot_time

            ds0.resize((3,))
            ds0[:] = [10.0, 20.0, 30.0]

            ds1.resize((2,))
            ds1[:] = [100.0, 200.0]

            np.testing.assert_array_equal(ds0[:], [10.0, 20.0, 30.0])
            np.testing.assert_array_equal(ds1[:], [100.0, 200.0])

    def test_write_per_sim_persists_after_reopen(self, tmp_path):
        """Per-sim data written through DataStore persists after close and reopen."""
        data = np.array([3.14, 2.72], dtype=np.float64)

        ds = self._create_full_store(tmp_path)
        ds.macroscale_out[1].snapshot_time.resize((2,))
        ds.macroscale_out[1].snapshot_time[:] = data
        ds.close()

        with DataStore("write_full", str(tmp_path)) as ds2:
            np.testing.assert_array_equal(
                ds2.macroscale_out[1].snapshot_time[:], data
            )

    # -- read-only mode blocks writes --

    def test_read_only_dataset_write_raises(self, tmp_path):
        """Writing to a dataset through a read-only DataStore raises an error."""
        # Create and populate a dataset, then close
        ds = self._create_micro_store(tmp_path)
        ds.microscale_out.pli_first_time.resize((3,))
        ds.microscale_out.pli_first_time[:] = [1.0, 2.0, 3.0]
        ds.close()

        # Reopen read-only
        with DataStore("write_test", str(tmp_path)) as ds2:
            dataset = ds2.microscale_out.pli_first_time
            with pytest.raises(OSError):
                dataset[:] = [9.0, 9.0, 9.0]

    def test_read_only_dataset_resize_raises(self, tmp_path):
        """Resizing a dataset through a read-only DataStore raises an error."""
        ds = self._create_micro_store(tmp_path)
        ds.close()

        with DataStore("write_test", str(tmp_path)) as ds2:
            dataset = ds2.microscale_out.pli_first_time
            with pytest.raises((OSError, RuntimeError)):
                dataset.resize((10,))


# ---------------------------------------------------------------------------
#  DerivedDataCollection unit tests
# ---------------------------------------------------------------------------


class TestDerivedDataCollection:
    """Unit tests for DerivedDataCollection using a mock generator."""

    @staticmethod
    def _make_spec():
        """Return the macroscale_in spec from the v2.0.0 dataspec."""
        return dataspec["v2.0.0"]["macroscale_in"]

    @staticmethod
    def _make_generator(call_counter=None):
        """Return a mock generator that produces fake macroscale_in data.

        :param call_counter: If provided, a list that gets an item appended
            on each call (to count invocations).
        """

        def _gen():
            if call_counter is not None:
                call_counter.append(1)
            return {
                "bin_edge_proportions": np.arange(101, dtype=np.float64),
                "bin_edge_tpa_leaving_time": np.arange(101, dtype=np.float64),
                "binned_fiber_degrade_time": np.zeros((50, 100), dtype=np.float64),
                "binned_fiber_degraded": np.zeros(100, dtype=np.uint16),
                "edge_grid_neighbors": np.zeros((10, 8), dtype=np.uint32),
                "params": {"micro_params": {}, "macro_params": {"forced_unbind": 0.5}},
            }

        return _gen

    def test_datasets_lists_all_spec_datasets(self):
        """datasets property returns all dataset names from the spec."""
        spec = self._make_spec()
        dc = DerivedDataCollection("macroscale_in", spec, self._make_generator())
        assert set(dc.datasets) == set(spec.data.keys())

    def test_repr_shows_pending_before_access(self):
        """__repr__ shows 'pending' before any dataset is accessed."""
        spec = self._make_spec()
        dc = DerivedDataCollection("macroscale_in", spec, self._make_generator())
        r = repr(dc)
        assert "pending" in r
        assert "macroscale_in" in r

    def test_repr_shows_generated_after_access(self):
        """__repr__ shows 'generated' after a dataset is accessed."""
        spec = self._make_spec()
        dc = DerivedDataCollection("macroscale_in", spec, self._make_generator())
        _ = dc.bin_edge_proportions
        r = repr(dc)
        assert "generated" in r

    def test_contains_for_spec_datasets(self):
        """__contains__ returns True for spec dataset names."""
        spec = self._make_spec()
        dc = DerivedDataCollection("macroscale_in", spec, self._make_generator())
        assert "bin_edge_proportions" in dc
        assert "binned_fiber_degraded" in dc
        assert "totally_fake" not in dc

    def test_getattr_triggers_generator_once(self):
        """Accessing a dataset triggers the generator exactly once."""
        counter = []
        spec = self._make_spec()
        dc = DerivedDataCollection(
            "macroscale_in", spec, self._make_generator(call_counter=counter)
        )
        assert len(counter) == 0

        _ = dc.bin_edge_proportions
        assert len(counter) == 1

        # Second access should not call again
        _ = dc.binned_fiber_degraded
        assert len(counter) == 1

    def test_multiple_accesses_reuse_cached_data(self):
        """Multiple dataset accesses return the same cached arrays."""
        spec = self._make_spec()
        dc = DerivedDataCollection("macroscale_in", spec, self._make_generator())
        a = dc.bin_edge_proportions
        b = dc.bin_edge_proportions
        assert a is b

    def test_returns_numpy_arrays(self):
        """Accessed datasets are numpy arrays."""
        spec = self._make_spec()
        dc = DerivedDataCollection("macroscale_in", spec, self._make_generator())
        assert isinstance(dc.bin_edge_proportions, np.ndarray)
        assert isinstance(dc.binned_fiber_degrade_time, np.ndarray)

    def test_unknown_dataset_raises_attribute_error(self):
        """Accessing a name not in the spec raises AttributeError."""
        spec = self._make_spec()
        dc = DerivedDataCollection("macroscale_in", spec, self._make_generator())
        with pytest.raises(AttributeError, match="no dataset"):
            _ = dc.totally_fake

    def test_params_property(self):
        """params property returns the params dict from the generator."""
        spec = self._make_spec()
        dc = DerivedDataCollection("macroscale_in", spec, self._make_generator())
        params = dc.params
        assert isinstance(params, dict)
        assert "macro_params" in params
        assert params["macro_params"]["forced_unbind"] == 0.5


# ---------------------------------------------------------------------------
#  DataStore macroscale_in integration tests
# ---------------------------------------------------------------------------


def _write_micro_datasets_for_macro_in(h5file, n_sims=100):
    """Create microscale_out datasets suitable for generate_macroscale_in.

    Requires n_sims divisible by 100.  Creates non-empty datasets with
    at least one True in both tpa_unbound_by_pli and tpa_unbound_kinetic
    to avoid division by zero in forced_unbind calculation.
    """
    spec = dataspec["v2.0.0"]["microscale_out"]
    rng = np.random.default_rng(42)
    for name, ds_spec in spec.data.items():
        if ds_spec.data_location is None:
            continue
        path = ds_spec.data_location
        if name == "pli_first_time":
            data = rng.uniform(0, 100, n_sims).astype(ds_spec.dtype)
        elif name == "tpa_leaving_time":
            data = rng.uniform(0, 100, n_sims).astype(ds_spec.dtype)
        elif name == "fiber_degraded":
            data = rng.choice([True, False], n_sims)
        elif name == "sim_final_time":
            data = rng.uniform(50, 200, n_sims).astype(ds_spec.dtype)
        elif name == "tpa_unbound_by_pli":
            # Ensure at least one True
            data = np.zeros(n_sims, dtype=bool)
            data[:n_sims // 2] = True
        elif name == "tpa_unbound_kinetic":
            # Ensure at least one True
            data = np.zeros(n_sims, dtype=bool)
            data[n_sims // 2:] = True
        elif ds_spec.dtype == h5py.string_dtype():
            data = np.array(["log line"] * n_sims, dtype=object)
            h5file.create_dataset(path, data=data, dtype=h5py.string_dtype())
            continue
        elif ds_spec.dtype == np.bool:
            data = rng.choice([True, False], n_sims)
        elif np.issubdtype(ds_spec.dtype, np.integer):
            data = np.arange(n_sims, dtype=ds_spec.dtype)
        else:
            data = rng.random(n_sims).astype(ds_spec.dtype)
        h5file.create_dataset(path, data=data)


class TestDataStoreMacroscaleIn:
    """Integration tests for macroscale_in via DataStore."""

    def test_macroscale_in_in_collections(self, tmp_path):
        """macroscale_in appears in ds.collections when conditions are met."""
        filepath = tmp_path / "macro_in.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f, micro_simulations=100)
            _write_micro_datasets_for_macro_in(f, n_sims=100)
            _write_macro_attrs(f, n_sims=2)
            _write_macro_datasets(f, n_sims=2)

        with DataStore("macro_in", str(tmp_path)) as ds:
            assert "macroscale_in" in ds.collections

    def test_macroscale_in_absent_without_macro_params(self, tmp_path):
        """macroscale_in absent when only microscale_out present (no macro_params)."""
        filepath = tmp_path / "no_macro.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f, micro_simulations=100)
            _write_micro_datasets_for_macro_in(f, n_sims=100)

        with DataStore("no_macro", str(tmp_path)) as ds:
            assert "macroscale_in" not in ds.collections

    def test_macroscale_in_absent_when_datasets_empty(self, tmp_path):
        """macroscale_in absent when microscale_out datasets are empty."""
        filepath = tmp_path / "empty_ds.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f, micro_simulations=0)
            # Write empty (zero-length) datasets
            spec = dataspec["v2.0.0"]["microscale_out"]
            for name, ds_spec in spec.data.items():
                if ds_spec.data_location is None:
                    continue
                if ds_spec.dtype == h5py.string_dtype():
                    f.create_dataset(
                        ds_spec.data_location,
                        shape=(0,),
                        maxshape=(None,),
                        dtype=h5py.string_dtype(),
                    )
                else:
                    f.create_dataset(
                        ds_spec.data_location,
                        shape=(0,),
                        maxshape=(None,),
                        dtype=ds_spec.dtype,
                    )
            _write_macro_attrs(f, n_sims=2)
            _write_macro_datasets(f, n_sims=2)

        with DataStore("empty_ds", str(tmp_path)) as ds:
            assert "macroscale_in" not in ds.collections

    def test_dot_access_returns_derived_collection(self, tmp_path):
        """ds.macroscale_in returns a DerivedDataCollection."""
        filepath = tmp_path / "dot_acc.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f, micro_simulations=100)
            _write_micro_datasets_for_macro_in(f, n_sims=100)
            _write_macro_attrs(f, n_sims=2)
            _write_macro_datasets(f, n_sims=2)

        with DataStore("dot_acc", str(tmp_path)) as ds:
            assert isinstance(ds.macroscale_in, DerivedDataCollection)

    def test_bin_edge_proportions_shape(self, tmp_path):
        """ds.macroscale_in.bin_edge_proportions has shape (101,)."""
        filepath = tmp_path / "bep.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f, micro_simulations=100)
            _write_micro_datasets_for_macro_in(f, n_sims=100)
            _write_macro_attrs(f, n_sims=2)
            _write_macro_datasets(f, n_sims=2)

        with DataStore("bep", str(tmp_path)) as ds:
            arr = ds.macroscale_in.bin_edge_proportions
            assert isinstance(arr, np.ndarray)
            assert arr.shape == (101,)

    def test_datasets_lists_all_five(self, tmp_path):
        """ds.macroscale_in.datasets returns all five dataset names."""
        filepath = tmp_path / "dslist.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f, micro_simulations=100)
            _write_micro_datasets_for_macro_in(f, n_sims=100)
            _write_macro_attrs(f, n_sims=2)
            _write_macro_datasets(f, n_sims=2)

        with DataStore("dslist", str(tmp_path)) as ds:
            names = ds.macroscale_in.datasets
            expected = {
                "bin_edge_proportions",
                "bin_edge_tpa_leaving_time",
                "binned_fiber_degrade_time",
                "binned_fiber_degraded",
                "edge_grid_neighbors",
            }
            assert set(names) == expected

    def test_lazy_generation(self, tmp_path):
        """Data is lazily generated — _data is None before first access."""
        filepath = tmp_path / "lazy.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f, micro_simulations=100)
            _write_micro_datasets_for_macro_in(f, n_sims=100)
            _write_macro_attrs(f, n_sims=2)
            _write_macro_datasets(f, n_sims=2)

        with DataStore("lazy", str(tmp_path)) as ds:
            coll = ds.macroscale_in
            assert coll._data is None

            # Trigger generation
            _ = coll.bin_edge_proportions
            assert coll._data is not None

    def test_params_has_forced_unbind(self, tmp_path):
        """Generated params contain the forced_unbind calculation."""
        filepath = tmp_path / "params.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f, micro_simulations=100)
            _write_micro_datasets_for_macro_in(f, n_sims=100)
            _write_macro_attrs(f, n_sims=2)
            _write_macro_datasets(f, n_sims=2)

        with DataStore("params", str(tmp_path)) as ds:
            params = ds.macroscale_in.params
            assert "macro_params" in params
            assert "forced_unbind" in params["macro_params"]
            fu = params["macro_params"]["forced_unbind"]
            assert 0.0 < fu < 1.0

    def test_binned_fiber_degrade_time_shape(self, tmp_path):
        """binned_fiber_degrade_time has shape (n_per_bin, 100)."""
        filepath = tmp_path / "bfdt.h5"
        n_sims = 100
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f, micro_simulations=n_sims)
            _write_micro_datasets_for_macro_in(f, n_sims=n_sims)
            _write_macro_attrs(f, n_sims=2)
            _write_macro_datasets(f, n_sims=2)

        with DataStore("bfdt", str(tmp_path)) as ds:
            arr = ds.macroscale_in.binned_fiber_degrade_time
            assert arr.shape == (n_sims // 100, 100)
