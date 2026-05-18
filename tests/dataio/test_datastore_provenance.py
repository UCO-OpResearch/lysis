"""Tests for :meth:`DataStore.stamp_provenance` and
:meth:`DataStore.read_init_provenance`.

The stamp method centralises provenance writing for three "kinds" of
stamp (``init``, ``execution``, ``binary``) and replaces the inline
logic that used to live in :meth:`DataStore.import_collection`.
"""

import h5py
import numpy as np
import pytest

from lysis.config.constants import CONST
from lysis.config.parameters import MicroParameters
from lysis.dataio.datastore import DataStore


@pytest.fixture
def micro_only_ds(tmp_path):
    """A DataStore with only the microscale collection created."""
    mp = MicroParameters()
    ds = DataStore.create("run-01", str(tmp_path), mp)
    yield ds
    ds.close()


@pytest.fixture
def micro_and_macro_ds(tmp_path):
    """A DataStore with both micro and macro params groups."""
    mp = MicroParameters()
    ds = DataStore.create("run-02", str(tmp_path), mp)
    # Manually create the macro_data group with an empty attrs dict so
    # stamp_provenance has a target.  (initialize_macroscale also creates
    # this group, but it requires populated microscale datasets — too
    # much setup for these unit tests.)
    ds._file.require_group("macro_data")
    ds._file.flush()
    yield ds
    ds.close()


# ----------------------------------------------------------------------
# stamp_provenance: validation
# ----------------------------------------------------------------------


class TestStampValidation:
    def test_bad_scale_raises_value_error(self, micro_only_ds):
        with pytest.raises(ValueError, match="scale must be"):
            micro_only_ds.stamp_provenance("oops", "init")

    def test_bad_kind_raises_value_error(self, micro_only_ds):
        with pytest.raises(ValueError, match="kind must be"):
            micro_only_ds.stamp_provenance("micro", "oops")

    def test_binary_without_executable_raises(self, micro_only_ds):
        with pytest.raises(ValueError, match="requires an 'executable'"):
            micro_only_ds.stamp_provenance("micro", "binary")

    def test_missing_macro_group_raises(self, micro_only_ds):
        # Only micro_data exists in this fixture.
        with pytest.raises(ValueError, match="macro_data"):
            micro_only_ds.stamp_provenance("macro", "init")

    def test_read_only_mode_raises(self, tmp_path):
        mp = MicroParameters()
        ds = DataStore.create("ro", str(tmp_path), mp)
        ds.close()
        with DataStore("ro", str(tmp_path), mode="r") as ro_ds:
            with pytest.raises(IOError, match="read-only"):
                ro_ds.stamp_provenance("micro", "init")


# ----------------------------------------------------------------------
# stamp_provenance: init
# ----------------------------------------------------------------------


class TestStampInit:
    def test_init_writes_expected_attrs(self, micro_only_ds, tmp_path):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            micro_only_ds.stamp_provenance("micro", "init")
        # Verify on disk.
        h5_path = tmp_path / "run-01.h5"
        micro_only_ds.close()
        with h5py.File(str(h5_path), "r") as f:
            attrs = f["micro_data"].attrs
            assert CONST.INIT_VERSION_ATTR in attrs
            assert CONST.INIT_DIRTY_ATTR in attrs
            assert CONST.INIT_TIMESTAMP_ATTR in attrs
            assert CONST.INIT_HOSTNAME_ATTR in attrs
            # Dirty stored as string state.
            dirty = attrs[CONST.INIT_DIRTY_ATTR]
            if isinstance(dirty, bytes):
                dirty = dirty.decode()
            assert dirty in {"clean", "dirty", "unknown"}

    def test_macro_init_writes_to_macro_data(self, micro_and_macro_ds, tmp_path):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            micro_and_macro_ds.stamp_provenance("macro", "init")
        h5_path = tmp_path / "run-02.h5"
        micro_and_macro_ds.close()
        with h5py.File(str(h5_path), "r") as f:
            assert CONST.INIT_VERSION_ATTR in f["macro_data"].attrs
            assert CONST.INIT_VERSION_ATTR not in f["micro_data"].attrs


# ----------------------------------------------------------------------
# stamp_provenance: execution
# ----------------------------------------------------------------------


class TestStampExecution:
    def test_execution_writes_expected_attrs(self, micro_only_ds, tmp_path):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            micro_only_ds.stamp_provenance("micro", "execution")
        h5_path = tmp_path / "run-01.h5"
        micro_only_ds.close()
        with h5py.File(str(h5_path), "r") as f:
            attrs = f["micro_data"].attrs
            assert CONST.EXECUTION_VERSION_ATTR in attrs
            assert CONST.EXECUTION_DIRTY_ATTR in attrs
            assert CONST.EXECUTION_TIMESTAMP_ATTR in attrs
            assert CONST.EXECUTION_HOSTNAME_ATTR in attrs
            # Execution dirty stays a bool for backward compat.
            assert isinstance(attrs[CONST.EXECUTION_DIRTY_ATTR], (bool, np.bool_))


# ----------------------------------------------------------------------
# stamp_provenance: binary
# ----------------------------------------------------------------------


class TestStampBinary:
    def test_binary_writes_three_attrs(self, micro_only_ds, tmp_path, monkeypatch):
        # Stub gather_binary_provenance via the binary module so we don't
        # actually need a working executable on disk.
        from lysis.tools.provenance import binary as binary_mod
        monkeypatch.setattr(
            binary_mod,
            "query_binary_version",
            lambda exe, **kw: ("abc123", "clean", "Intel"),
        )
        micro_only_ds.stamp_provenance(
            "micro", "binary", executable="/fake/bin"
        )
        h5_path = tmp_path / "run-01.h5"
        micro_only_ds.close()
        with h5py.File(str(h5_path), "r") as f:
            attrs = f["micro_data"].attrs
            assert attrs[CONST.BINARY_COMMIT_ATTR].decode() == "abc123" \
                if isinstance(attrs[CONST.BINARY_COMMIT_ATTR], bytes) \
                else attrs[CONST.BINARY_COMMIT_ATTR] == "abc123"
            assert attrs[CONST.BINARY_DIRTY_ATTR] in (b"clean", "clean")
            assert attrs[CONST.BINARY_COMPILER_ATTR] in (b"Intel", "Intel")
            # No override flag when not provided.
            assert CONST.STALE_BINARY_OVERRIDE_ATTR not in attrs

    def test_binary_merges_override_attrs(
        self, micro_only_ds, tmp_path, monkeypatch
    ):
        from lysis.tools.provenance import binary as binary_mod
        monkeypatch.setattr(
            binary_mod,
            "query_binary_version",
            lambda exe, **kw: ("abc123", "clean", "Intel"),
        )
        micro_only_ds.stamp_provenance(
            "micro",
            "binary",
            executable="/fake/bin",
            binary_override={CONST.STALE_BINARY_OVERRIDE_ATTR: True},
        )
        h5_path = tmp_path / "run-01.h5"
        micro_only_ds.close()
        with h5py.File(str(h5_path), "r") as f:
            attrs = f["micro_data"].attrs
            assert bool(attrs[CONST.STALE_BINARY_OVERRIDE_ATTR]) is True


# ----------------------------------------------------------------------
# read_init_provenance
# ----------------------------------------------------------------------


class TestReadInitProvenance:
    def test_returns_none_when_not_stamped(self, micro_only_ds):
        # Newly-created DataStore.create() does not auto-stamp init.
        assert micro_only_ds.read_init_provenance("micro") is None

    def test_returns_dict_after_stamping(self, micro_only_ds):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            micro_only_ds.stamp_provenance("micro", "init")
        result = micro_only_ds.read_init_provenance("micro")
        assert result is not None
        assert CONST.INIT_VERSION_ATTR in result
        assert CONST.INIT_DIRTY_ATTR in result

    def test_bad_scale_raises(self, micro_only_ds):
        with pytest.raises(ValueError, match="scale must be"):
            micro_only_ds.read_init_provenance("oops")

    def test_missing_group_raises(self, micro_only_ds):
        with pytest.raises(ValueError, match="macro_data"):
            micro_only_ds.read_init_provenance("macro")
