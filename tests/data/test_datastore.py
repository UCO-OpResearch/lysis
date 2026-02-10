"""Comprehensive pytest unit tests for the lysis.data.datastore module.

Tests cover:
- DataStatus enum values and flag operations
- DataStore.__init__ with various modes and pre-existing data
- h5_tree output formatting
- DataStore.__repr__ and __str__
"""

import os

import numpy as np
import h5py
import pytest

from lysis.data.datastore import DataStore, DataStatus, h5_tree


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
#  DataStore.__init__ tests
# ---------------------------------------------------------------------------


class TestDataStoreInit:
    """Tests for DataStore construction with various modes and HDF5 states."""

    def test_create_new_file_write_mode(self, tmp_path):
        """Creating a DataStore in write mode produces a new HDF5 file."""
        ds = DataStore("run1", str(tmp_path), mode="w")
        try:
            expected_file = tmp_path / "run1.h5"
            assert expected_file.exists()
            assert ds._status["self"] == DataStatus.INITIALIZED
            assert ds._status["micro"] == DataStatus.NONE
            assert ds._status["macro"] == DataStatus.NONE
        finally:
            ds._data.close()

    def test_open_existing_empty_read_mode(self, tmp_path):
        """Opening an existing empty HDF5 file in read mode succeeds."""
        # Pre-create an empty HDF5 file
        filepath = tmp_path / "empty.h5"
        with h5py.File(filepath, "w"):
            pass

        ds = DataStore("empty", str(tmp_path), mode="r")
        try:
            assert ds._status["self"] == DataStatus.INITIALIZED
            assert ds._status["micro"] == DataStatus.NONE
            assert ds._status["macro"] == DataStatus.NONE
        finally:
            ds._data.close()

    def test_open_with_micro_data(self, tmp_path):
        """micro status is INITIALIZED when the file contains a micro_data group."""
        filepath = tmp_path / "micro.h5"
        with h5py.File(filepath, "w") as f:
            f.create_group("micro_data")

        ds = DataStore("micro", str(tmp_path), mode="r")
        try:
            assert ds._status["micro"] == DataStatus.INITIALIZED
            assert ds._status["macro"] == DataStatus.NONE
        finally:
            ds._data.close()

    def test_open_with_macro_data(self, tmp_path):
        """macro status is INITIALIZED when the file contains a macro_data group."""
        filepath = tmp_path / "macro.h5"
        with h5py.File(filepath, "w") as f:
            f.create_group("macro_data")

        ds = DataStore("macro", str(tmp_path), mode="r")
        try:
            assert ds._status["macro"] == DataStatus.INITIALIZED
            assert ds._status["micro"] == DataStatus.NONE
        finally:
            ds._data.close()

    def test_open_with_both_groups(self, tmp_path):
        """Both micro and macro are INITIALIZED when both groups are present."""
        filepath = tmp_path / "both.h5"
        with h5py.File(filepath, "w") as f:
            f.create_group("micro_data")
            f.create_group("macro_data")

        ds = DataStore("both", str(tmp_path), mode="r")
        try:
            assert ds._status["micro"] == DataStatus.INITIALIZED
            assert ds._status["macro"] == DataStatus.INITIALIZED
        finally:
            ds._data.close()

    def test_append_mode(self, tmp_path):
        """Opening in append mode succeeds and reports correct status."""
        # Pre-create a file so append has something to open
        filepath = tmp_path / "app.h5"
        with h5py.File(filepath, "w"):
            pass

        ds = DataStore("app", str(tmp_path), mode="a")
        try:
            assert ds._status["self"] == DataStatus.INITIALIZED
            assert ds._status["micro"] == DataStatus.NONE
            assert ds._status["macro"] == DataStatus.NONE
            # Verify the file is actually open in append mode
            assert ds._mode == "a"
        finally:
            ds._data.close()


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
        # h5py scalar datasets have shape (), which formats as "()"
        assert "()" in output

    def test_scalar_fallback_label(self):
        """When accessing .shape raises TypeError, the item is labelled (scalar).

        This exercises the except-TypeError branch in h5_tree by using a mock
        object whose .shape property raises TypeError.
        """
        from unittest.mock import MagicMock

        mock_item = MagicMock(spec=[])
        # Give mock_item a .shape that raises TypeError
        type(mock_item).shape = property(lambda self: (_ for _ in ()).throw(TypeError))

        # Build a fake group-like container that yields one item
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
        # The last item at each level uses the corner connector
        assert "└──" in output

    def test_multiple_items_tree_characters(self, tmp_path):
        """When a group has multiple children, non-last items use the tee character."""
        filepath = tmp_path / "tree_multi.h5"
        with h5py.File(filepath, "w") as f:
            f.create_dataset("alpha", data=np.array([1, 2]))
            f.create_dataset("beta", data=np.array([3, 4]))

        with h5py.File(filepath, "r") as f:
            output = h5_tree(f)

        # One item should be preceded by tee and the last by corner
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
#  DataStore.__repr__ and __str__ tests
# ---------------------------------------------------------------------------


class TestDataStoreReprStr:
    """Tests for the string representations of DataStore."""

    def test_repr_contains_path_and_status(self, tmp_path):
        """__repr__ includes the file path and the status dictionary."""
        ds = DataStore("reptest", str(tmp_path), mode="w")
        try:
            r = repr(ds)
            assert str(tmp_path) in r
            assert "reptest.h5" in r
            assert "self" in r
            assert "INITIALIZED" in r
        finally:
            ds._data.close()

    def test_str_returns_tree_output(self, tmp_path):
        """__str__ delegates to h5_tree and returns matching output."""
        # Create a file with some content so the tree is non-empty
        filepath = tmp_path / "strtest.h5"
        with h5py.File(filepath, "w") as f:
            f.create_dataset("data", data=np.ones(4))

        ds = DataStore("strtest", str(tmp_path), mode="r")
        try:
            s = str(ds)
            assert "data" in s
            assert "(4,)" in s
        finally:
            ds._data.close()

    def test_str_empty_datastore(self, tmp_path):
        """__str__ on an empty DataStore returns an empty string."""
        ds = DataStore("emptystr", str(tmp_path), mode="w")
        try:
            assert str(ds) == ""
        finally:
            ds._data.close()
