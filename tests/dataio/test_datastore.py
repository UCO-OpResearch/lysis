"""Comprehensive pytest unit tests for the lysis.dataio.datastore module.

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

from lysis.dataio.datastore import (
    COMPATIBLE_DATASPEC_VERSION,
    DataStore,
    DataCollection,
    DerivedDataCollection,
    ImportCollectionError,
    SimulationView,
    DataStatus,
    HDF5State,
    h5_tree,
    _group_path_from_spec,
)
from lysis.config.constants import CONST
from lysis.dataio.dataspec import DataCollectionSpec, DataSetSpec, dataspec
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
            elif len(ds_spec.shape) == 3:
                # Multi-dimensional dataset (tpa_location_snapshot: n × 2 × n_snapshots)
                data = np.zeros((2, 2, n_snapshots), dtype=ds_spec.dtype)
                h5file.create_dataset(path, data=data)
            else:
                data = np.random.rand(n_snapshots).astype(ds_spec.dtype)
                h5file.create_dataset(path, data=data)


def _fill_microscale_unbinding_data(filepath, n_sims=10):
    """Resize and fill microscale datasets in an HDF5 created by DataStore.create().

    DataStore.create() produces empty (shape 0) resizable datasets.  This helper
    fills them with dummy values so that initialize_macroscale() can compute
    forced_unbind.  Half the simulations have tpa_unbound_by_pli=True, the other
    half have tpa_unbound_kinetic=True, so forced_unbind will be exactly 0.5.
    """
    spec = dataspec["v2.0.0"]["microscale_out"]
    with h5py.File(filepath, "a") as f:
        for name, ds_spec in spec.data.items():
            if ds_spec.data_location is None:
                continue
            if ds_spec.dtype == h5py.string_dtype():
                continue
            path = ds_spec.data_location
            if path not in f:
                continue
            dataset = f[path]
            dataset.resize((n_sims,) + dataset.shape[1:])
            if name == "tpa_unbound_by_pli":
                dataset[:] = np.array(
                    [True] * (n_sims // 2) + [False] * (n_sims - n_sims // 2)
                )
            elif name == "tpa_unbound_kinetic":
                dataset[:] = np.array(
                    [False] * (n_sims // 2) + [True] * (n_sims - n_sims // 2)
                )
            elif ds_spec.dtype == np.bool_:
                dataset[:] = np.ones(n_sims, dtype=bool)
            elif np.issubdtype(ds_spec.dtype, np.integer):
                dataset[:] = np.arange(n_sims, dtype=ds_spec.dtype)
            else:
                dataset[:] = np.arange(n_sims, dtype=ds_spec.dtype)


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
#  HDF5State enum and hdf5_state property tests
# ---------------------------------------------------------------------------


class TestHDF5State:
    """Tests for the HDF5State enum and DataStore.hdf5_state property."""

    def test_all_values_exist(self):
        """All expected members are present on the enum."""
        assert hasattr(HDF5State, "MICRO_EMPTY")
        assert hasattr(HDF5State, "MICRO_FILLED")
        assert hasattr(HDF5State, "MACRO_EMPTY")
        assert hasattr(HDF5State, "MACRO_FILLED")

    def test_values_are_unique(self):
        """Every member has a unique underlying value."""
        values = [m.value for m in HDF5State]
        assert len(values) == len(set(values))

    def test_micro_empty_state(self, tmp_path):
        """micro_data with params but no datasets → MICRO_EMPTY."""
        h5_path = tmp_path / "micro-empty.h5"
        with h5py.File(str(h5_path), "w") as f:
            _write_micro_attrs(f)
        with DataStore("micro-empty", str(tmp_path)) as ds:
            assert ds.hdf5_state == HDF5State.MICRO_EMPTY

    def test_micro_filled_state(self, tmp_path):
        """micro_data with non-empty tpa_leaving_time → MICRO_FILLED."""
        h5_path = tmp_path / "micro-filled.h5"
        with h5py.File(str(h5_path), "w") as f:
            _write_micro_attrs(f)
            f.create_dataset("micro_data/tpa_leaving_time", data=np.array([1.0]))
        with DataStore("micro-filled", str(tmp_path)) as ds:
            assert ds.hdf5_state == HDF5State.MICRO_FILLED

    def test_macro_empty_state(self, tmp_path):
        """micro filled + macro_data with params but no datasets → MACRO_EMPTY."""
        h5_path = tmp_path / "macro-empty.h5"
        with h5py.File(str(h5_path), "w") as f:
            _write_micro_attrs(f)
            f.create_dataset("micro_data/tpa_leaving_time", data=np.array([1.0]))
            _write_macro_attrs(f)
        with DataStore("macro-empty", str(tmp_path)) as ds:
            assert ds.hdf5_state == HDF5State.MACRO_EMPTY

    def test_macro_filled_state(self, tmp_path):
        """micro filled + macro filled → MACRO_FILLED."""
        h5_path = tmp_path / "macro-filled.h5"
        with h5py.File(str(h5_path), "w") as f:
            _write_micro_attrs(f)
            f.create_dataset("micro_data/tpa_leaving_time", data=np.array([1.0]))
            _write_macro_attrs(f)
            f.create_dataset(
                "macro_data/sim_00/snapshot_time",
                data=np.array([1.0], dtype=np.float64),
            )
        with DataStore("macro-filled", str(tmp_path)) as ds:
            assert ds.hdf5_state == HDF5State.MACRO_FILLED

    def test_inconsistent_state(self, tmp_path):
        """macro_data present but micro not filled → INCONSISTENT."""
        h5_path = tmp_path / "inconsistent.h5"
        with h5py.File(str(h5_path), "w") as f:
            _write_micro_attrs(f)
            _write_macro_attrs(f)
        with DataStore("inconsistent", str(tmp_path)) as ds:
            assert ds.hdf5_state == HDF5State.INCONSISTENT


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

    _N_MICRO_SIMS = 10  # number of dummy microscale simulations to write

    def _create_micro_store(self, tmp_path, run_code="init_macro"):
        """Helper: create a microscale-only DataStore with dummy simulation data.

        initialize_macroscale() now requires non-empty tpa_unbound_by_pli and
        tpa_unbound_kinetic datasets to compute forced_unbind, so we resize and
        fill all microscale datasets after DataStore.create().
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()
        ds = DataStore.create(run_code, str(tmp_path), micro)
        # Datasets created by DataStore.create() are empty (shape 0) and
        # resizable. Fill them with dummy data so forced_unbind can be computed.
        filepath = tmp_path / f"{run_code}.h5"
        _fill_microscale_unbinding_data(filepath, n_sims=self._N_MICRO_SIMS)
        # Close and reopen so the DataStore sees the newly written data.
        ds.close()
        return DataStore(run_code, str(tmp_path), mode="a")

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
        """initialize_macroscale() does not erase existing microscale datasets."""
        ds = self._create_micro_store(tmp_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            macro = MacroParameters(
                micro_params=ds.micro_params, macro_simulations=2
            )
        ds.initialize_macroscale(macro)
        try:
            # Microscale datasets must still exist and retain their pre-existing data.
            dataset = ds.microscale_out.pli_first_time
            assert isinstance(dataset, h5py.Dataset)
            assert dataset.shape[0] == self._N_MICRO_SIMS
        finally:
            ds.close()

    def test_initialize_in_context_manager(self, tmp_path):
        """initialize_macroscale() works inside a context manager."""
        # Use _create_micro_store (which writes microscale data) then reopen
        # as a context manager so we can verify file handle closure on exit.
        ds_setup = self._create_micro_store(tmp_path, run_code="ctx_init")
        ds_setup.close()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()
        with DataStore("ctx_init", str(tmp_path), mode="a") as ds:
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

    def test_initialize_computes_forced_unbind(self, tmp_path):
        """initialize_macroscale() stores forced_unbind computed from microscale data.

        _create_micro_store sets tpa_unbound_by_pli[:n//2] = True and
        tpa_unbound_kinetic[n//2:] = True, so forced_unbind = 0.5.
        The value passed in MacroParameters (nan by default) must be replaced.
        """
        ds = self._create_micro_store(tmp_path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            # Pass the default (nan) to confirm it gets overwritten.
            macro = MacroParameters(micro_params=ds.micro_params, macro_simulations=2)
        ds.initialize_macroscale(macro)
        try:
            stored = ds.macro_params.forced_unbind
            assert 0.0 < stored < 1.0, f"expected a valid fraction, got {stored}"
            # n//2 forced out of n total → 0.5
            assert stored == pytest.approx(0.5)
        finally:
            ds.close()

    def test_initialize_raises_on_empty_microscale(self, tmp_path):
        """initialize_macroscale() raises ValueError when microscale datasets are empty."""
        # Create a store with empty datasets (do NOT fill them).
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters()
        ds = DataStore.create("empty_micro", str(tmp_path), micro)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            macro = MacroParameters(micro_params=ds.micro_params, macro_simulations=2)
        try:
            with pytest.raises(ValueError, match="microscale datasets are empty"):
                ds.initialize_macroscale(macro)
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
        # initialize_macroscale() requires non-empty tpa_unbound_* datasets.
        # Write minimal dummy microscale data before proceeding.
        filepath = tmp_path / f"{run_code}.h5"
        _fill_microscale_unbinding_data(filepath)
        ds.close()
        ds = DataStore(run_code, str(tmp_path), mode="a")
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
#  import_collection tests
# ---------------------------------------------------------------------------


class TestDataStoreImportCollection:
    """Tests for DataStore.import_collection().

    Covers importing microscale_out and macroscale_out from both Fortran
    (v1.95.0) and HDF5 (v2.0.0) sources, plus all precondition error paths.
    """

    # File codes and param overrides matching tests/fixtures/fortran_sample/
    _MICRO_FILE_CODE = "_PLG2_tPA01_TB-xiii"
    _MACRO_FILE_CODE = "_TB-xiii__21_105"
    _PARAM_OVERRIDES = {
        "fibrinogen_length": "45nm",
        "fibrinogen_radius": "1.2nm",
        "micro_log_lvl": 40,
        "micro_version": "micro_rates",
        "snap_proportion": 0.66666667,
    }
    # The truncated fixture has 50 000 microscale simulations and 1 macro sim.
    _N_MICRO_SIMS = 50000
    _N_MACRO_SIMS = 1  # only simulation 00 in the fixture

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _create_empty_store(tmp_path, run_code="import_test", n_sims=10):
        """Create a DataStore via DataStore.create() with empty datasets."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            micro = MicroParameters(micro_simulations=n_sims)
        return DataStore.create(run_code, str(tmp_path), micro)

    @staticmethod
    def _create_filled_micro_hdf5(tmp_path, n_micro=10):
        """Write a filled v2.0.0 HDF5 containing only microscale_out datasets.

        Returns the HDF5 file path as a string.
        """
        filepath = tmp_path / "micro_source.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f, **{"micro_simulations": n_micro})
            _write_micro_datasets(f, n_sims=n_micro)
        return str(filepath)

    @staticmethod
    def _create_filled_macro_hdf5(tmp_path, n_micro=10, n_macro=3,
                                  n_snapshots=5):
        """Write a filled v2.0.0 HDF5 with both micro and macro datasets.

        Returns the HDF5 file path as a string.
        """
        filepath = tmp_path / "macro_source.h5"
        with h5py.File(filepath, "w") as f:
            _write_micro_attrs(f, **{"micro_simulations": n_micro})
            _write_micro_datasets(f, n_sims=n_micro)
            _write_macro_attrs(f, n_sims=n_macro)
            _write_macro_datasets(f, n_sims=n_macro, n_snapshots=n_snapshots)
        return str(filepath)

    @staticmethod
    def _prepare_store_for_macro_import(tmp_path, run_code, n_micro, n_macro):
        """Create a DataStore ready to receive macroscale data.

        Flow: create → fill microscale → initialize_macroscale → return.
        The returned DataStore is in 'a' mode with empty macroscale datasets.
        """
        ds = TestDataStoreImportCollection._create_empty_store(
            tmp_path, run_code=run_code, n_sims=n_micro
        )
        filepath = tmp_path / f"{run_code}.h5"
        _fill_microscale_unbinding_data(filepath, n_sims=n_micro)
        ds.close()

        ds = DataStore(run_code, str(tmp_path), mode="a")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            macro = MacroParameters(
                micro_params=ds.micro_params, macro_simulations=n_macro
            )
        ds.initialize_macroscale(macro)
        return ds

    # ------------------------------------------------------------------
    # Microscale import from Fortran fixture (v1.95.0)
    # ------------------------------------------------------------------

    def test_import_microscale_from_fortran(self, fortran_sample_path, tmp_path):
        """Importing microscale_out from v1.95.0 Fortran data fills all datasets."""
        ds = self._create_empty_store(
            tmp_path, n_sims=self._N_MICRO_SIMS
        )
        try:
            ds.import_collection(
                "microscale_out",
                "v1.95.0",
                fortran_sample_path,
                [self._MICRO_FILE_CODE],
                param_overrides=self._PARAM_OVERRIDES,
            )
            assert "microscale_out" in ds.collections
            # At least one numeric dataset must be non-empty.
            pli = ds.microscale_out.pli_first_time
            assert isinstance(pli, h5py.Dataset)
            assert pli.shape[0] == self._N_MICRO_SIMS
        finally:
            ds.close()

    def test_import_microscale_all_datasets_filled(self, fortran_sample_path,
                                                    tmp_path):
        """All numeric microscale_out datasets have shape[0] == N after import."""
        ds = self._create_empty_store(
            tmp_path, n_sims=self._N_MICRO_SIMS
        )
        try:
            ds.import_collection(
                "microscale_out",
                "v1.95.0",
                fortran_sample_path,
                [self._MICRO_FILE_CODE],
                param_overrides=self._PARAM_OVERRIDES,
            )
            spec = dataspec["v2.0.0"]["microscale_out"]
            for ds_name, ds_spec in spec.data.items():
                if ds_spec.data_location is None:
                    continue
                if ds_spec.dtype == h5py.string_dtype():
                    continue
                dataset = ds.microscale_out.__getattr__(ds_name)
                assert dataset.shape[0] == self._N_MICRO_SIMS, (
                    f"Dataset '{ds_name}' has shape {dataset.shape}, "
                    f"expected first dim {self._N_MICRO_SIMS}"
                )
        finally:
            ds.close()

    def test_import_microscale_data_integrity(self, fortran_sample_path,
                                              tmp_path):
        """Imported data matches an independent read + convert of the same source."""
        from lysis.dataio.fileops import read_data_collection
        from lysis.dataio.dataconvert import convert_data
        from lysis.dataio.dataspec import dataspec as _ds

        # Independent pipeline: read → override → convert
        ref_data = read_data_collection(
            fortran_sample_path,
            [_ds["v1.95.0"]["microscale_out"]],
            [self._MICRO_FILE_CODE],
        )
        for key, value in self._PARAM_OVERRIDES.items():
            for section in ref_data["params"].values():
                if isinstance(section, dict):
                    section[key] = value
        ref_converted = convert_data(ref_data, "v1.95.0", "v2.0.0")

        ds = self._create_empty_store(
            tmp_path, n_sims=self._N_MICRO_SIMS
        )
        try:
            ds.import_collection(
                "microscale_out",
                "v1.95.0",
                fortran_sample_path,
                [self._MICRO_FILE_CODE],
                param_overrides=self._PARAM_OVERRIDES,
            )
            # Compare a numeric dataset.
            np.testing.assert_array_equal(
                ds.microscale_out.pli_first_time[:],
                ref_converted["pli_first_time"],
            )
            np.testing.assert_array_equal(
                ds.microscale_out.tpa_leaving_time[:],
                ref_converted["tpa_leaving_time"],
            )
        finally:
            ds.close()

    def test_import_microscale_preserves_params(self, fortran_sample_path,
                                                tmp_path):
        """micro_params loaded from HDF5 attributes are unchanged after import."""
        ds = self._create_empty_store(
            tmp_path, n_sims=self._N_MICRO_SIMS
        )
        original_n = ds.micro_params.micro_simulations
        try:
            ds.import_collection(
                "microscale_out",
                "v1.95.0",
                fortran_sample_path,
                [self._MICRO_FILE_CODE],
                param_overrides=self._PARAM_OVERRIDES,
            )
            assert ds.micro_params.micro_simulations == original_n
        finally:
            ds.close()

    def test_import_microscale_preserves_mode(self, fortran_sample_path,
                                               tmp_path):
        """DataStore mode is unchanged after import."""
        ds = self._create_empty_store(
            tmp_path, n_sims=self._N_MICRO_SIMS
        )
        try:
            ds.import_collection(
                "microscale_out",
                "v1.95.0",
                fortran_sample_path,
                [self._MICRO_FILE_CODE],
                param_overrides=self._PARAM_OVERRIDES,
            )
            assert ds.mode == "a"
        finally:
            ds.close()

    def test_import_microscale_returns_none(self, fortran_sample_path,
                                            tmp_path):
        """import_collection() returns None (in-place modification)."""
        ds = self._create_empty_store(
            tmp_path, n_sims=self._N_MICRO_SIMS
        )
        try:
            result = ds.import_collection(
                "microscale_out",
                "v1.95.0",
                fortran_sample_path,
                [self._MICRO_FILE_CODE],
                param_overrides=self._PARAM_OVERRIDES,
            )
            assert result is None
        finally:
            ds.close()

    # ------------------------------------------------------------------
    # Microscale import from HDF5 source (v2.0.0)
    # ------------------------------------------------------------------

    def test_import_microscale_from_hdf5(self, tmp_path):
        """Importing microscale_out from a v2.0.0 HDF5 source works."""
        n_micro = 10
        source_path = self._create_filled_micro_hdf5(tmp_path, n_micro=n_micro)

        ds = self._create_empty_store(tmp_path, run_code="hdf5_import",
                                      n_sims=n_micro)
        try:
            ds.import_collection(
                "microscale_out", "v2.0.0", source_path, [""]
            )
            assert "microscale_out" in ds.collections
            assert ds.microscale_out.pli_first_time.shape[0] == n_micro
        finally:
            ds.close()

    def test_import_microscale_hdf5_tag_alias(self, tmp_path):
        """The 'hdf5' tag alias is resolved correctly as a source spec."""
        n_micro = 5
        source_path = self._create_filled_micro_hdf5(tmp_path, n_micro=n_micro)

        ds = self._create_empty_store(tmp_path, run_code="tag_import",
                                      n_sims=n_micro)
        try:
            ds.import_collection(
                "microscale_out", "hdf5", source_path, [""]
            )
            assert ds.microscale_out.pli_first_time.shape[0] == n_micro
        finally:
            ds.close()

    def test_import_microscale_hdf5_data_integrity(self, tmp_path):
        """Data round-trips correctly through HDF5 → DataStore import."""
        n_micro = 8
        source_path = self._create_filled_micro_hdf5(tmp_path, n_micro=n_micro)

        # Read the reference data directly from the source HDF5
        with h5py.File(source_path, "r") as f:
            ref = f["micro_data/pli_first_time"][:]

        ds = self._create_empty_store(tmp_path, run_code="integrity_test",
                                      n_sims=n_micro)
        try:
            ds.import_collection(
                "microscale_out", "v2.0.0", source_path, [""]
            )
            np.testing.assert_array_equal(
                ds.microscale_out.pli_first_time[:], ref
            )
        finally:
            ds.close()

    # ------------------------------------------------------------------
    # Macroscale import from HDF5 source (v2.0.0)
    # ------------------------------------------------------------------

    def test_import_macroscale_from_hdf5(self, tmp_path):
        """Importing macroscale_out from a v2.0.0 HDF5 source works."""
        n_micro, n_macro = 10, 3
        source_path = self._create_filled_macro_hdf5(
            tmp_path, n_micro=n_micro, n_macro=n_macro
        )

        ds = self._prepare_store_for_macro_import(
            tmp_path, "macro_import", n_micro=n_micro, n_macro=n_macro
        )
        try:
            ds.import_collection(
                "macroscale_out", "v2.0.0", source_path, [""]
            )
            assert "macroscale_out" in ds.collections
            # Verify at least one per-sim dataset is non-empty.
            snap = ds.macroscale_out[0].snapshot_time
            assert snap.shape[0] > 0
        finally:
            ds.close()

    def test_import_macroscale_all_sims_filled(self, tmp_path):
        """All macro simulations have non-empty datasets after import."""
        n_micro, n_macro, n_snapshots = 10, 3, 5
        source_path = self._create_filled_macro_hdf5(
            tmp_path, n_micro=n_micro, n_macro=n_macro, n_snapshots=n_snapshots
        )

        ds = self._prepare_store_for_macro_import(
            tmp_path, "all_sims", n_micro=n_micro, n_macro=n_macro
        )
        try:
            ds.import_collection(
                "macroscale_out", "v2.0.0", source_path, [""]
            )
            for sim_idx in range(n_macro):
                snap = ds.macroscale_out[sim_idx].snapshot_time
                assert snap.shape[0] == n_snapshots, (
                    f"sim {sim_idx}: expected {n_snapshots} snapshots, "
                    f"got {snap.shape[0]}"
                )
        finally:
            ds.close()

    def test_import_macroscale_per_sim_data_integrity(self, tmp_path):
        """snapshot_time round-trips correctly through HDF5 → DataStore import."""
        n_micro, n_macro, n_snapshots = 10, 2, 4
        source_path = self._create_filled_macro_hdf5(
            tmp_path, n_micro=n_micro, n_macro=n_macro, n_snapshots=n_snapshots
        )

        # Capture reference data from the source HDF5
        with h5py.File(source_path, "r") as f:
            ref_0 = f["macro_data/sim_00/snapshot_time"][:]

        ds = self._prepare_store_for_macro_import(
            tmp_path, "integrity_macro", n_micro=n_micro, n_macro=n_macro
        )
        try:
            ds.import_collection(
                "macroscale_out", "v2.0.0", source_path, [""]
            )
            np.testing.assert_array_equal(
                ds.macroscale_out[0].snapshot_time[:], ref_0
            )
        finally:
            ds.close()

    # ------------------------------------------------------------------
    # Macroscale import from Fortran fixture (v1.95.0)
    # ------------------------------------------------------------------

    def test_import_macroscale_from_fortran(self, fortran_sample_path, tmp_path):
        """Importing macroscale_out from v1.95.0 Fortran data fills datasets."""
        # Create store, import microscale, then initialize and import macroscale
        ds = self._create_empty_store(
            tmp_path, run_code="fort_macro", n_sims=self._N_MICRO_SIMS
        )
        ds.import_collection(
            "microscale_out",
            "v1.95.0",
            fortran_sample_path,
            [self._MICRO_FILE_CODE],
            param_overrides=self._PARAM_OVERRIDES,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            macro = MacroParameters(
                micro_params=ds.micro_params,
                macro_simulations=self._N_MACRO_SIMS,
            )
        ds.initialize_macroscale(macro)

        try:
            ds.import_collection(
                "macroscale_out",
                "v1.95.0",
                fortran_sample_path,
                [self._MACRO_FILE_CODE],
                param_overrides=self._PARAM_OVERRIDES,
            )
            assert "macroscale_out" in ds.collections
            snap = ds.macroscale_out[0].snapshot_time
            assert snap.shape[0] > 0
        finally:
            ds.close()

    # ------------------------------------------------------------------
    # Error cases
    # ------------------------------------------------------------------

    def test_import_read_only_raises(self, tmp_path):
        """import_collection() raises IOError on a read-only DataStore."""
        ds = self._create_empty_store(tmp_path)
        ds.close()

        ds = DataStore("import_test", str(tmp_path))  # read-only by default
        try:
            with pytest.raises(IOError, match="read-only"):
                ds.import_collection(
                    "microscale_out", "v2.0.0", str(tmp_path), [""]
                )
        finally:
            ds.close()

    def test_import_invalid_collection_raises(self, tmp_path):
        """import_collection() raises ValueError for unsupported collection names."""
        ds = self._create_empty_store(tmp_path)
        try:
            with pytest.raises(ValueError, match="Only 'microscale_out'"):
                ds.import_collection(
                    "macroscale_in", "v2.0.0", str(tmp_path), [""]
                )
        finally:
            ds.close()

    def test_import_microscale_not_created_raises(self, tmp_path):
        """Importing microscale_out raises if create() was not called first."""
        # Create a bare HDF5 with only the version attribute.
        bare_path = tmp_path / "bare.h5"
        with h5py.File(bare_path, "w") as f:
            _write_version_attr(f)

        ds = DataStore("bare", str(tmp_path), mode="a")
        try:
            with pytest.raises(ValueError, match="does not exist"):
                ds.import_collection(
                    "microscale_out", "v2.0.0", str(bare_path), [""]
                )
        finally:
            ds.close()

    def test_import_macroscale_before_initialize_raises(self, tmp_path):
        """Importing macroscale_out raises if initialize_macroscale() was not called."""
        ds = self._create_empty_store(tmp_path)
        try:
            with pytest.raises(ValueError, match="does not exist"):
                ds.import_collection(
                    "macroscale_out", "v2.0.0", str(tmp_path), [""]
                )
        finally:
            ds.close()

    def test_import_microscale_already_filled_raises(self, fortran_sample_path,
                                                      tmp_path):
        """Second import_collection() call raises when datasets already contain data."""
        ds = self._create_empty_store(
            tmp_path, n_sims=self._N_MICRO_SIMS
        )
        ds.import_collection(
            "microscale_out",
            "v1.95.0",
            fortran_sample_path,
            [self._MICRO_FILE_CODE],
            param_overrides=self._PARAM_OVERRIDES,
        )
        try:
            with pytest.raises(ValueError, match="already contains data"):
                ds.import_collection(
                    "microscale_out",
                    "v1.95.0",
                    fortran_sample_path,
                    [self._MICRO_FILE_CODE],
                    param_overrides=self._PARAM_OVERRIDES,
                )
        finally:
            ds.close()

    def test_import_macroscale_sim_count_mismatch_raises(self, tmp_path):
        """Importing macroscale_out raises if simulation count doesn't match."""
        n_micro, n_macro = 10, 3
        # Source has 3 macro sims, DataStore initialized for 2.
        source_path = self._create_filled_macro_hdf5(
            tmp_path, n_micro=n_micro, n_macro=n_macro
        )
        ds = self._prepare_store_for_macro_import(
            tmp_path, "mismatch", n_micro=n_micro, n_macro=2  # 2, not 3
        )
        try:
            with pytest.raises(ValueError, match="simulation"):
                ds.import_collection(
                    "macroscale_out", "v2.0.0", source_path, [""]
                )
        finally:
            ds.close()

    # ------------------------------------------------------------------
    # Revert-to-empty + loud failure on write-phase errors (#85)
    # ------------------------------------------------------------------

    @staticmethod
    def _fail_after_first_write(monkeypatch, error):
        """Patch ``_write_array_to_dataset`` to succeed once, then raise.

        Leaves one partial write on disk so the rollback has something to undo.
        """
        original = DataStore._write_array_to_dataset
        calls = {"n": 0}

        def flaky(self, ds_spec, arr, sim=None):
            calls["n"] += 1
            if calls["n"] >= 2:
                raise error
            return original(self, ds_spec, arr, sim=sim)

        monkeypatch.setattr(DataStore, "_write_array_to_dataset", flaky)

    def test_import_microscale_reverts_to_empty_on_failure(
        self, tmp_path, monkeypatch
    ):
        """A write-phase failure rolls microscale_out back to MICRO_EMPTY."""
        n_micro = 10
        source_path = self._create_filled_micro_hdf5(tmp_path, n_micro=n_micro)
        ds = self._create_empty_store(
            tmp_path, run_code="revert_micro", n_sims=n_micro
        )
        boom = RuntimeError("disk full mid-write")
        self._fail_after_first_write(monkeypatch, boom)

        try:
            with pytest.raises(ImportCollectionError) as excinfo:
                ds.import_collection(
                    "microscale_out", "v2.0.0", source_path, [""]
                )
            # Loud failure chains the original error.
            assert excinfo.value.__cause__ is boom

            # File is back to a clean, freshly-initialized empty state.
            assert ds.hdf5_state == HDF5State.MICRO_EMPTY
            micro_spec = dataspec[COMPATIBLE_DATASPEC_VERSION]["microscale_out"]
            for ds_spec in micro_spec.data.values():
                if ds_spec.data_location is None:
                    continue
                assert ds._file[ds_spec.data_location].shape[0] == 0
            assert CONST.CONVERTED_FROM_ATTR not in ds._file.attrs
            # Init provenance was re-stamped by the create() rollback.
            assert ds.read_init_provenance("micro") is not None
        finally:
            ds.close()

    def test_import_macroscale_reverts_to_empty_on_failure(
        self, tmp_path, monkeypatch
    ):
        """A write-phase failure rolls macroscale_out back to MACRO_EMPTY,
        leaving the filled microscale data intact."""
        n_micro, n_macro, n_snapshots = 10, 3, 5
        source_path = self._create_filled_macro_hdf5(
            tmp_path, n_micro=n_micro, n_macro=n_macro, n_snapshots=n_snapshots
        )
        ds = self._prepare_store_for_macro_import(
            tmp_path, "revert_macro", n_micro=n_micro, n_macro=n_macro
        )
        boom = RuntimeError("write blew up mid-sim")
        self._fail_after_first_write(monkeypatch, boom)

        try:
            with pytest.raises(ImportCollectionError) as excinfo:
                ds.import_collection(
                    "macroscale_out", "v2.0.0", source_path, [""]
                )
            assert excinfo.value.__cause__ is boom

            # macroscale reverted to empty; microscale untouched.
            assert ds.hdf5_state == HDF5State.MACRO_EMPTY
            macro_spec = dataspec[COMPATIBLE_DATASPEC_VERSION]["macroscale_out"]
            snap_loc = macro_spec.data["snapshot_time"].data_location
            for sim_idx in range(n_macro):
                assert ds._file[snap_loc.format(sim=sim_idx)].shape[0] == 0
            micro_spec = dataspec[COMPATIBLE_DATASPEC_VERSION]["microscale_out"]
            micro_check = micro_spec.data["tpa_leaving_time"].data_location
            assert ds._file[micro_check].shape[0] > 0
            # Macro init provenance was re-stamped by the rollback.
            assert ds.read_init_provenance("macro") is not None
        finally:
            ds.close()

    def test_import_succeeds_after_revert(self, tmp_path, monkeypatch):
        """After a reverted failure the collection is empty again, so a fresh
        import passes the 'must be empty' precondition and fills the data."""
        n_micro = 10
        source_path = self._create_filled_micro_hdf5(tmp_path, n_micro=n_micro)
        ds = self._create_empty_store(
            tmp_path, run_code="reimport", n_sims=n_micro
        )

        def always_boom(self, *args, **kwargs):
            raise RuntimeError("nope")

        monkeypatch.setattr(
            DataStore, "_write_array_to_dataset", always_boom
        )
        try:
            with pytest.raises(ImportCollectionError):
                ds.import_collection(
                    "microscale_out", "v2.0.0", source_path, [""]
                )
            assert ds.hdf5_state == HDF5State.MICRO_EMPTY

            # Restore the real writer and re-import: must succeed.
            monkeypatch.undo()
            ds.import_collection(
                "microscale_out", "v2.0.0", source_path, [""]
            )
            assert ds.hdf5_state == HDF5State.MICRO_FILLED
            assert ds.microscale_out.pli_first_time.shape[0] == n_micro
        finally:
            ds.close()

    # ------------------------------------------------------------------
    # Rollback-also-fails: still ImportCollectionError + recovery dump (#85)
    # ------------------------------------------------------------------

    def test_rollback_failure_still_raises_and_dumps_recovery(
        self, tmp_path, monkeypatch
    ):
        """If the rollback itself fails, the caller still gets an
        ImportCollectionError (not the raw rollback error) and the in-memory
        params are dumped to a re-feedable recovery CSV."""
        n_micro = 10
        source_path = self._create_filled_micro_hdf5(tmp_path, n_micro=n_micro)
        ds = self._create_empty_store(
            tmp_path, run_code="rollback_boom", n_sims=n_micro
        )

        write_boom = RuntimeError("write failed")
        rollback_boom = OSError("disk full during rollback")

        def fail_write(self, *args, **kwargs):
            raise write_boom

        def fail_rollback(self, collection_name):
            raise rollback_boom

        monkeypatch.setattr(DataStore, "_write_array_to_dataset", fail_write)
        monkeypatch.setattr(
            DataStore, "_revert_collection_to_empty", fail_rollback
        )

        try:
            with pytest.raises(ImportCollectionError) as excinfo:
                ds.import_collection(
                    "microscale_out", "v2.0.0", source_path, [""]
                )
            # Caller sees ImportCollectionError; rollback error is the cause,
            # original write error is the cause's context.
            assert excinfo.value.__cause__ is rollback_boom
            assert excinfo.value.__cause__.__context__ is write_boom
            assert "param_recovery" in str(excinfo.value)

            recovery = tmp_path / "rollback_boom_param_recovery.csv"
            assert recovery.exists()
        finally:
            ds.close()

        # The recovery CSV is re-feedable to init-experiment.
        from lysis.config.experiment import Experiment

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            exp = Experiment.from_csv(
                str(recovery), str(tmp_path), name="recover", dry_run=True
            )
        assert exp.runs[0].micro_params.micro_simulations == n_micro

    def test_rollback_failure_json_fallback(self, tmp_path, monkeypatch):
        """When the CSV dump fails too, recovery falls back to JSON."""
        import lysis.config.parameters as params_mod

        n_micro = 10
        source_path = self._create_filled_micro_hdf5(tmp_path, n_micro=n_micro)
        ds = self._create_empty_store(
            tmp_path, run_code="json_fallback", n_sims=n_micro
        )

        def fail_write(self, *args, **kwargs):
            raise RuntimeError("write failed")

        def fail_rollback(self, collection_name):
            raise OSError("rollback failed")

        def fail_csv(*args, **kwargs):
            raise OSError("cannot write csv")

        monkeypatch.setattr(DataStore, "_write_array_to_dataset", fail_write)
        monkeypatch.setattr(
            DataStore, "_revert_collection_to_empty", fail_rollback
        )
        monkeypatch.setattr(params_mod, "write_params_csv", fail_csv)

        try:
            with pytest.raises(ImportCollectionError):
                ds.import_collection(
                    "microscale_out", "v2.0.0", source_path, [""]
                )
        finally:
            ds.close()

        assert not (tmp_path / "json_fallback_param_recovery.csv").exists()
        assert (tmp_path / "json_fallback_param_recovery.json").exists()


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
