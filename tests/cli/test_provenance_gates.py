"""End-to-end tests for the src/lysis/ provenance gates.

Exercises the dirty-tree warning, the ``--allow-dirty`` override, the
init-time provenance stamp written to each HDF5 file, and the
init→run commit-match check.

The session-scoped fixture in ``tests/conftest.py`` sets
``LYSIS_ALLOW_DIRTY=1`` and ``LYSIS_ALLOW_COMMIT_MISMATCH=1`` so the
rest of the test suite isn't blocked by an in-progress working tree.
These tests clear those env vars locally to exercise the gates.
"""

import csv
import os
import warnings
from unittest.mock import MagicMock, patch

import h5py
import pytest
from click.testing import CliRunner

from lysis.cli import cli
from lysis.config.constants import CONST
from lysis.config.parameters import MacroParameters, MicroParameters
from lysis.dataio.datastore import DataStore
from lysis.tools.provenance import execution as execution_mod

# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------


def _write_minimal_csv(path, rows):
    param_names = list(rows[0].keys())
    run_codes = [r.get("run_code", f"run-{i:02d}") for i, r in enumerate(rows)]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["parameter"] + run_codes)
        for param in param_names:
            if param == "run_code":
                continue
            writer.writerow([param] + [str(r.get(param, "")) for r in rows])


def _one_run_csv(tmp_path):
    import inspect

    micro = MicroParameters().to_basedict()
    micro_ind = set(inspect.signature(MicroParameters).parameters.keys())
    micro = {k: v for k, v in micro.items() if k in micro_ind}

    macro = MacroParameters(micro_params=MicroParameters()).to_basedict()
    macro_ind = set(inspect.signature(MacroParameters).parameters.keys())
    macro_ind.discard("micro_params")
    macro = {k: v for k, v in macro.items() if k in macro_ind}

    row = {**micro, **macro, "run_code": "gate-test", "run_description": "x"}
    csv_path = tmp_path / "gate.csv"
    _write_minimal_csv(csv_path, [row])
    return csv_path


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def force_dirty(monkeypatch):
    """Pretend src/lysis/ is dirty and clear the allow-dirty env override."""
    monkeypatch.delenv(CONST.LYSIS_ALLOW_DIRTY_ENV, raising=False)
    monkeypatch.delenv(CONST.LYSIS_ALLOW_COMMIT_MISMATCH_ENV, raising=False)
    monkeypatch.setattr(execution_mod, "_resolve_dirty", lambda: "dirty")
    monkeypatch.setattr(execution_mod, "_resolve_version", lambda: "deadbeefdeadbeef")
    # Reset the once-per-process dedup flag so each test sees the warning.
    execution_mod._dirty_warning_emitted = False
    yield


@pytest.fixture
def force_clean(monkeypatch):
    """Pretend src/lysis/ is clean at a stable known commit."""
    monkeypatch.delenv(CONST.LYSIS_ALLOW_DIRTY_ENV, raising=False)
    monkeypatch.delenv(CONST.LYSIS_ALLOW_COMMIT_MISMATCH_ENV, raising=False)
    monkeypatch.setattr(execution_mod, "_resolve_dirty", lambda: "clean")
    monkeypatch.setattr(execution_mod, "_resolve_version", lambda: "commit-A" * 4)
    execution_mod._dirty_warning_emitted = False
    yield


# ----------------------------------------------------------------------
# Top-level dirty warning
# ----------------------------------------------------------------------


class TestTopLevelDirtyWarning:
    def test_warning_appears_on_dirty_tree(self, runner, force_dirty, tmp_path):
        # Use a real subcommand so the group callback fires.  We invoke
        # init-experiment with --allow-dirty so it would otherwise
        # succeed; the test only asserts that the warning text appears.
        csv_path = _one_run_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()
        result = runner.invoke(
            cli,
            [
                "init-experiment",
                str(csv_path),
                str(data_root),
                "--name",
                "warn",
                "--allow-dirty",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "src/lysis/ has uncommitted changes" in result.output


# ----------------------------------------------------------------------
# Dirty gate on init-experiment
# ----------------------------------------------------------------------


class TestInitExperimentDirtyGate:
    def test_dirty_blocks_init_experiment(self, runner, force_dirty, tmp_path):
        csv_path = _one_run_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()
        result = runner.invoke(
            cli,
            ["init-experiment", str(csv_path), str(data_root), "--name", "g"],
        )
        assert result.exit_code != 0
        assert "uncommitted changes" in result.output

    def test_allow_dirty_flag_proceeds(self, runner, force_dirty, tmp_path):
        csv_path = _one_run_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()
        result = runner.invoke(
            cli,
            [
                "init-experiment",
                str(csv_path),
                str(data_root),
                "--name",
                "g",
                "--allow-dirty",
            ],
        )
        assert result.exit_code == 0, result.output
        # init_dirty="dirty" stamped on the resulting HDF5.
        with h5py.File(str(data_root / "g" / "gate-test.h5"), "r") as f:
            val = f["micro_data"].attrs[CONST.INIT_DIRTY_ATTR]
            if isinstance(val, bytes):
                val = val.decode()
            assert val == "dirty"

    def test_env_var_proceeds(self, runner, force_dirty, tmp_path, monkeypatch):
        monkeypatch.setenv(CONST.LYSIS_ALLOW_DIRTY_ENV, "1")
        csv_path = _one_run_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()
        result = runner.invoke(
            cli,
            ["init-experiment", str(csv_path), str(data_root), "--name", "g"],
        )
        assert result.exit_code == 0, result.output


# ----------------------------------------------------------------------
# Clean init writes init_version
# ----------------------------------------------------------------------


class TestInitStampsCommit:
    def test_init_experiment_stamps_version(self, runner, force_clean, tmp_path):
        csv_path = _one_run_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()
        result = runner.invoke(
            cli,
            ["init-experiment", str(csv_path), str(data_root), "--name", "g"],
        )
        assert result.exit_code == 0, result.output
        with h5py.File(str(data_root / "g" / "gate-test.h5"), "r") as f:
            recorded = f["micro_data"].attrs[CONST.INIT_VERSION_ATTR]
            if isinstance(recorded, bytes):
                recorded = recorded.decode()
            assert recorded == "commit-A" * 4


# ----------------------------------------------------------------------
# run-micro commit-match check
# ----------------------------------------------------------------------


class TestRunMicroCommitMatch:
    @patch("lysis.execution.fortran_micro.FortranMicro")
    def test_legacy_file_warns_and_proceeds(
        self, mock_cls, runner, force_clean, tmp_path
    ):
        """HDF5 without init_version → warning, but run continues."""
        # Build a minimal valid v2.0.0 HDF5, then strip the init_* attrs that
        # DataStore.create() auto-stamps to simulate a file that pre-dates the
        # init-stamp feature.
        mp = MicroParameters()
        ds = DataStore.create("legacy", str(tmp_path), mp)
        attrs = ds._file["micro_data"].attrs
        for key in (
            CONST.INIT_VERSION_ATTR,
            CONST.INIT_DIRTY_ATTR,
            CONST.INIT_TIMESTAMP_ATTR,
            CONST.INIT_HOSTNAME_ATTR,
        ):
            if key in attrs:
                del attrs[key]
        ds.close()

        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", UserWarning)
            result = runner.invoke(
                cli,
                [
                    "run-micro",
                    str(tmp_path / "legacy.h5"),
                    "--executable",
                    "/fake/bin",
                ],
            )

        assert result.exit_code == 0, result.output
        # The "legacy file" warning surfaces (from Click runner or the env).
        assert (
            any("pre-dates" in str(w.message) for w in caught)
            or "pre-dates" in result.output
        )

    @patch("lysis.execution.fortran_micro.FortranMicro")
    def test_mismatch_blocks_run_micro(
        self, mock_cls, runner, force_clean, tmp_path, monkeypatch
    ):
        # First stamp a DIFFERENT commit as the init_version.
        mp = MicroParameters()
        ds = DataStore.create("mm", str(tmp_path), mp)
        monkeypatch.setattr(execution_mod, "_resolve_dirty", lambda: "clean")
        monkeypatch.setattr(execution_mod, "_resolve_version", lambda: "old-commit-zzz")
        ds.stamp_provenance("micro", "init")
        ds.close()

        # Now switch the "current" commit back to commit-A.
        monkeypatch.setattr(execution_mod, "_resolve_version", lambda: "commit-A" * 4)

        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            [
                "run-micro",
                str(tmp_path / "mm.h5"),
                "--executable",
                "/fake/bin",
            ],
        )
        assert result.exit_code != 0
        assert "commit-match check failed" in result.output

    @patch("lysis.execution.fortran_micro.FortranMicro")
    def test_mismatch_override_warns_and_continues(
        self, mock_cls, runner, force_clean, tmp_path, monkeypatch
    ):
        mp = MicroParameters()
        ds = DataStore.create("mm2", str(tmp_path), mp)
        monkeypatch.setattr(execution_mod, "_resolve_dirty", lambda: "clean")
        monkeypatch.setattr(execution_mod, "_resolve_version", lambda: "old-commit")
        ds.stamp_provenance("micro", "init")
        ds.close()
        monkeypatch.setattr(execution_mod, "_resolve_version", lambda: "commit-A" * 4)

        mock_cls.from_hdf5.return_value = MagicMock()

        # Capture the warning via pytest.warns so it is asserted *and* consumed
        # — otherwise it escapes to pytest's run-wide warnings summary even
        # though firing here is the expected --allow-commit-mismatch behaviour.
        with pytest.warns(UserWarning, match="commit-match check failed") as record:
            result = runner.invoke(
                cli,
                [
                    "run-micro",
                    str(tmp_path / "mm2.h5"),
                    "--executable",
                    "/fake/bin",
                    "--allow-commit-mismatch",
                ],
            )
        assert result.exit_code == 0, result.output
        assert any("old-commit" in str(w.message) for w in record)
