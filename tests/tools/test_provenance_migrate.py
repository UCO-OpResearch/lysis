"""Tests for the provenance-attr migration: pure transform + file round-trip (#61)."""

import importlib.util
from pathlib import Path

import h5py
import numpy as np
import pytest

from lysis.config.constants import CONST
from lysis.tools.provenance.migrate import migrate_provenance_group


# --------------------------------------------------------------------------- #
# Load the standalone migration script as a module (it lives under scripts/).
# --------------------------------------------------------------------------- #
_SCRIPT = (
    Path(__file__).resolve().parents[2] / "scripts" / "migrate_provenance_attrs.py"
)
_spec = importlib.util.spec_from_file_location("migrate_provenance_attrs", _SCRIPT)
migrate_script = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(migrate_script)


# --------------------------------------------------------------------------- #
# Pure transform
# --------------------------------------------------------------------------- #
class TestMigrateProvenanceGroup:
    def test_full_old_group_maps_to_new(self):
        old = {
            "execution_version": "abc123",
            "execution_dirty": np.False_,
            "execution_timestamp": "2025-01-01T00:00:00",
            "execution_hostname": "node01",
            "execution_backend": "fortran",
            "binary_commit": "d" * 40,
            "binary_dirty": "clean",
            "binary_compiler": "ifort 2023",
            "binary_source": "historical:" + "d" * 40,
            "init_version": "zzz",          # untouched
            "init_dirty": "clean",          # untouched
            "dataspec_version": "v2.0.0",   # untouched
        }
        sets, deletes = migrate_provenance_group(old)

        assert sets[CONST.PIPELINE_VERSION_ATTR] == "abc123"
        assert sets[CONST.PIPELINE_DIRTY_ATTR] == "clean"  # bool False -> "clean"
        assert sets[CONST.PIPELINE_TIMESTAMP_ATTR] == "2025-01-01T00:00:00"
        assert sets[CONST.PIPELINE_HOSTNAME_ATTR] == "node01"
        assert sets[CONST.BACKEND_TYPE_ATTR] == "fortran"
        assert sets[CONST.BACKEND_COMMIT_ATTR] == "d" * 40
        assert sets[CONST.BACKEND_DIRTY_ATTR] == "clean"
        assert sets[CONST.BACKEND_COMPILER_ATTR] == "ifort 2023"
        assert sets[CONST.BACKEND_HISTORICAL_ATTR] is True

        # Old keys removed; init_*/dataspec_version are not in deletes.
        assert deletes == {
            "execution_version", "execution_dirty", "execution_timestamp",
            "execution_hostname", "execution_backend", "binary_commit",
            "binary_dirty", "binary_compiler", "binary_source",
        }
        assert "init_version" not in deletes
        assert "dataspec_version" not in deletes
        assert CONST.INIT_VERSION_ATTR not in sets

    def test_dirty_true_maps_to_dirty_string(self):
        sets, _ = migrate_provenance_group({"execution_dirty": True})
        assert sets[CONST.PIPELINE_DIRTY_ATTR] == "dirty"

    def test_backend_type_defaults_to_fortran_without_execution_backend(self):
        sets, _ = migrate_provenance_group(
            {"binary_commit": "a" * 40, "binary_dirty": "clean"}
        )
        assert sets[CONST.BACKEND_TYPE_ATTR] == "fortran"

    def test_backend_type_preserves_python(self):
        sets, _ = migrate_provenance_group(
            {"execution_backend": "python", "binary_commit": "a" * 40}
        )
        assert sets[CONST.BACKEND_TYPE_ATTR] == "python"

    def test_no_old_attrs_is_noop(self):
        # Already-migrated group → empty result (idempotent).
        already = {
            CONST.PIPELINE_VERSION_ATTR: "x",
            CONST.BACKEND_COMMIT_ATTR: "y",
            CONST.INIT_VERSION_ATTR: "z",
        }
        sets, deletes = migrate_provenance_group(already)
        assert sets == {}
        assert deletes == set()

    def test_binary_source_absent_no_historical_flag(self):
        sets, _ = migrate_provenance_group({"binary_commit": "a" * 40})
        assert CONST.BACKEND_HISTORICAL_ATTR not in sets

    def test_stale_override_renamed(self):
        sets, deletes = migrate_provenance_group(
            {"binary_commit": "a" * 40, "stale_binary_override": True}
        )
        assert sets["stale_backend_override"] is True
        assert "stale_binary_override" in deletes


# --------------------------------------------------------------------------- #
# File round-trip via the script's migrate_file()
# --------------------------------------------------------------------------- #
def _make_old_file(path):
    """Write a minimal v2.0.0 file with OLD provenance attrs + a couple datasets."""
    with h5py.File(path, "w") as f:
        f.attrs["dataspec_version"] = "v2.0.0"
        g = f.create_group("micro_data")
        g.attrs["execution_version"] = "abc123"
        g.attrs["execution_dirty"] = np.False_
        g.attrs["execution_timestamp"] = "2025-01-01T00:00:00"
        g.attrs["execution_hostname"] = "node01"
        g.attrs["binary_commit"] = "d" * 40
        g.attrs["binary_dirty"] = "clean"
        g.attrs["binary_compiler"] = "ifort 2023"
        g.attrs["binary_source"] = "historical:" + "d" * 40
        g.attrs["init_version"] = "init-sha"   # must survive untouched
        g.attrs["init_dirty"] = "clean"
        g.create_dataset("ints", data=np.arange(50, dtype=np.int64))
        f.create_dataset(
            "micro_data/strings",
            data=np.array(["a", "bb", "ccc"], dtype=h5py.string_dtype()),
        )


def test_migrate_file_round_trip(tmp_path):
    path = str(tmp_path / "old.h5")
    _make_old_file(path)

    migrate_script.migrate_file(path, backup_dir=str(tmp_path / "bak"))

    with h5py.File(path, "r") as f:
        attrs = f["micro_data"].attrs
        # New names present, correctly typed.
        assert attrs[CONST.PIPELINE_VERSION_ATTR] == "abc123"
        assert attrs[CONST.PIPELINE_DIRTY_ATTR] == "clean"
        assert attrs[CONST.BACKEND_COMMIT_ATTR] == "d" * 40
        assert attrs[CONST.BACKEND_TYPE_ATTR] == "fortran"
        assert bool(attrs[CONST.BACKEND_HISTORICAL_ATTR]) is True
        # Old names gone.
        for old in ("execution_version", "execution_dirty", "binary_commit",
                    "binary_source", "binary_dirty", "binary_compiler"):
            assert old not in attrs
        # init_* untouched; datasets intact.
        assert attrs["init_version"] == "init-sha"
        assert attrs["init_dirty"] == "clean"
        assert list(f["micro_data/ints"][...]) == list(range(50))
        assert f.attrs["dataspec_version"] == "v2.0.0"

    # Backup was created.
    assert any((tmp_path / "bak").iterdir())


def test_migrate_file_idempotent(tmp_path):
    path = str(tmp_path / "old.h5")
    _make_old_file(path)
    migrate_script.migrate_file(path, backup_dir=str(tmp_path / "bak"))
    # Second pass: classifier says already-migrated, and migrate_file is a no-op.
    assert migrate_script.classify(path) == "already-migrated"
    assert migrate_script.migrate_file(path, backup_dir=str(tmp_path / "bak2")) == {}


def test_migrate_file_writes_sha_log(tmp_path):
    path = str(tmp_path / "old.h5")
    _make_old_file(path)
    bak = str(tmp_path / "bak")
    migrate_script.migrate_file(path, backup_dir=bak)

    # A per-file .sha256.txt sits next to the backup copy.
    flat = migrate_script._backup_basename(path)
    log = Path(bak) / (flat + ".sha256.txt")
    assert log.is_file()
    text = log.read_text()

    assert "[pre-edit]" in text and "[post-edit]" in text
    pre_block, post_block = text.split("[post-edit]")
    # Every dataset appears in both sections, and (since datasets are untouched)
    # each dataset's pre and post sha are identical.
    for ds in ("micro_data/ints", "micro_data/strings"):
        assert ds in pre_block and ds in post_block
        pre_sha = [l for l in pre_block.splitlines() if l.startswith(ds + "  ")][0].split()[1]
        post_sha = [l for l in post_block.splitlines() if l.startswith(ds + "  ")][0].split()[1]
        assert pre_sha == post_sha
        assert len(pre_sha) == 64  # sha256 hexdigest


def test_migrate_file_halts_on_dataset_change(tmp_path, monkeypatch):
    """If a dataset checksum differs post-edit, migrate_file raises (halt)."""
    path = str(tmp_path / "old.h5")
    _make_old_file(path)

    # Force every checksum call to return a unique value so pre != post.
    counter = {"n": 0}

    def _drifting_checksum(ds):
        counter["n"] += 1
        return f"sha-{counter['n']}"

    monkeypatch.setattr(migrate_script, "_dataset_checksum", _drifting_checksum)

    with pytest.raises(migrate_script.IntegrityError, match="dataset changed"):
        migrate_script.migrate_file(path, backup_dir=str(tmp_path / "bak"))
