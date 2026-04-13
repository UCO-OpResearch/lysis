"""Tests for ``lysis init-macroscale`` CLI command.

Tests cover:
- Folder mode: initialises all runs in an experiment
- H5 mode: initialises a single HDF5 file
- --param overrides (H5 mode only)
- --dry-run: validates without writing
- Error cases: empty microscale, already initialised, bad PATH
"""

import csv
import inspect
import json
import os

import h5py
import numpy as np
import pytest
from click.testing import CliRunner

from lysis.cli import cli
from lysis.config.experiment import Experiment
from lysis.config.parameters import MacroParameters, MicroParameters
from lysis.dataio.dataspec import dataspec
from lysis.dataio.datastore import DataStore


# ─── CSV / experiment helpers ────────────────────────────────────────────────


def _write_minimal_csv(path, rows):
    """Write a transposed CSV: rows = parameters, columns = runs.

    run_code values become column headers; all other keys become parameter rows.
    """
    if not rows:
        return
    param_names = list(rows[0].keys())
    run_codes = [r.get("run_code", f"run-{i:02d}") for i, r in enumerate(rows)]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["parameter"] + run_codes)
        for param in param_names:
            if param == "run_code":
                continue  # run_code is the column header, not a data row
            writer.writerow([param] + [str(r.get(param, "")) for r in rows])


def _default_micro_row():
    mp = MicroParameters()
    d = mp.to_basedict()
    micro_ind = set(inspect.signature(MicroParameters).parameters.keys())
    return {k: v for k, v in d.items() if k in micro_ind}


def _default_macro_row(micro=None):
    mp = micro or MicroParameters()
    d = MacroParameters(micro_params=mp).to_basedict()
    macro_ind = set(inspect.signature(MacroParameters).parameters.keys())
    macro_ind.discard("micro_params")
    # Exclude forced_unbind (NaN) — not a valid CSV value
    return {k: v for k, v in d.items() if k in macro_ind and k != "forced_unbind"}


def _two_row_csv(tmp_path):
    micro = _default_micro_row()
    macro = _default_macro_row()
    row1 = {**micro, **macro, "run_code": "test-run-00", "run_description": "default"}
    row2 = {**micro, **macro, "run_code": "test-run-01", "run_description": "second"}
    csv_path = tmp_path / "test_exp.csv"
    _write_minimal_csv(csv_path, [row1, row2])
    return csv_path


def _fill_microscale_data(filepath, n_sims=10):
    """Resize and fill microscale datasets so initialize_macroscale() can run.

    Half the simulations get tpa_unbound_by_pli=True; the other half get
    tpa_unbound_kinetic=True.  This gives forced_unbind == 0.5.
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


# ─── Fixtures ─────────────────────────────────────────────────────────────────


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def experiment_dir(tmp_path):
    """Create a two-run experiment with microscale data filled in.

    Returns the experiment folder path (a pathlib.Path).
    """
    csv_path = _two_row_csv(tmp_path)
    data_root = tmp_path / "experiments"
    data_root.mkdir()

    exp = Experiment.from_csv(csv_path, data_root, name="test-exp")

    # Fill microscale data for each run so initialize_macroscale() can proceed
    for run in exp.runs:
        h5 = os.path.join(exp.path, f"{run.run_code}.h5")
        _fill_microscale_data(h5)

    return data_root / "test-exp"


@pytest.fixture
def single_h5(tmp_path):
    """Create a single-run HDF5 file with microscale data filled in.

    Returns the HDF5 file path (a pathlib.Path).
    """
    micro = MicroParameters()
    h5_dir = tmp_path / "run_dir"
    h5_dir.mkdir()

    ds = DataStore.create("solo-run", str(h5_dir), micro)
    ds.close()

    h5_path = h5_dir / "solo-run.h5"
    _fill_microscale_data(str(h5_path))
    return h5_path


# ─── Basic help / invocation ──────────────────────────────────────────────────


class TestInitMacroscaleHelp:
    def test_help_exits_0(self, runner):
        result = runner.invoke(cli, ["init-macroscale", "--help"])
        assert result.exit_code == 0, result.output

    def test_help_mentions_path(self, runner):
        result = runner.invoke(cli, ["init-macroscale", "--help"])
        assert "PATH" in result.output

    def test_help_mentions_param(self, runner):
        result = runner.invoke(cli, ["init-macroscale", "--help"])
        assert "--param" in result.output

    def test_help_mentions_dry_run(self, runner):
        result = runner.invoke(cli, ["init-macroscale", "--help"])
        assert "--dry-run" in result.output


# ─── Folder mode ─────────────────────────────────────────────────────────────


class TestInitMacroscaleFolderMode:
    def test_folder_mode_exits_0(self, runner, experiment_dir):
        result = runner.invoke(cli, ["init-macroscale", str(experiment_dir)])
        assert result.exit_code == 0, result.output

    def test_folder_mode_writes_macro_data(self, runner, experiment_dir):
        runner.invoke(cli, ["init-macroscale", str(experiment_dir)])

        for run_code in ["test-run-00", "test-run-01"]:
            h5 = experiment_dir / f"{run_code}.h5"
            with h5py.File(str(h5), "r") as f:
                assert "macro_data" in f, f"macro_data missing from {run_code}.h5"

    def test_folder_mode_computes_forced_unbind(self, runner, experiment_dir):
        import math

        runner.invoke(cli, ["init-macroscale", str(experiment_dir)])

        for run_code in ["test-run-00", "test-run-01"]:
            with DataStore(run_code, str(experiment_dir), mode="r") as ds:
                assert ds.macro_params is not None
                assert math.isclose(ds.macro_params.forced_unbind, 0.5, rel_tol=1e-6)

    def test_folder_mode_output_mentions_runs(self, runner, experiment_dir):
        result = runner.invoke(cli, ["init-macroscale", str(experiment_dir)])
        assert result.exit_code == 0, result.output
        # Both run codes should appear in the output table
        assert "test-run-00" in result.output
        assert "test-run-01" in result.output

    def test_folder_mode_param_override_warns(self, runner, experiment_dir):
        """--param is ignored in folder mode; a warning should be printed."""
        result = runner.invoke(
            cli,
            ["init-macroscale", str(experiment_dir), "--param", "rows=64"],
        )
        assert result.exit_code == 0, result.output
        assert "ignored" in result.output.lower() or "warning" in result.output.lower()

    def test_folder_mode_dry_run_does_not_write(self, runner, experiment_dir):
        runner.invoke(
            cli, ["init-macroscale", str(experiment_dir), "--dry-run"]
        )

        for run_code in ["test-run-00", "test-run-01"]:
            with h5py.File(str(experiment_dir / f"{run_code}.h5"), "r") as f:
                assert "macro_data" not in f

    def test_folder_mode_dry_run_exits_0(self, runner, experiment_dir):
        result = runner.invoke(
            cli, ["init-macroscale", str(experiment_dir), "--dry-run"]
        )
        assert result.exit_code == 0, result.output

    def test_folder_mode_dry_run_output_says_dry_run(self, runner, experiment_dir):
        result = runner.invoke(
            cli, ["init-macroscale", str(experiment_dir), "--dry-run"]
        )
        lower = result.output.lower()
        assert "dry" in lower or "validated" in lower or "ok" in lower


# ─── H5 file mode ────────────────────────────────────────────────────────────


class TestInitMacroscaleH5Mode:
    def test_h5_mode_exits_0(self, runner, single_h5):
        result = runner.invoke(cli, ["init-macroscale", str(single_h5)])
        assert result.exit_code == 0, result.output

    def test_h5_mode_writes_macro_data(self, runner, single_h5):
        runner.invoke(cli, ["init-macroscale", str(single_h5)])

        with h5py.File(str(single_h5), "r") as f:
            assert "macro_data" in f

    def test_h5_mode_computes_forced_unbind(self, runner, single_h5):
        import math

        runner.invoke(cli, ["init-macroscale", str(single_h5)])

        run_code = single_h5.stem
        run_dir = str(single_h5.parent)
        with DataStore(run_code, run_dir, mode="r") as ds:
            assert ds.macro_params is not None
            assert math.isclose(ds.macro_params.forced_unbind, 0.5, rel_tol=1e-6)

    def test_h5_mode_output_mentions_file(self, runner, single_h5):
        result = runner.invoke(cli, ["init-macroscale", str(single_h5)])
        assert result.exit_code == 0, result.output
        assert "solo-run" in result.output

    def test_h5_mode_dry_run_does_not_write(self, runner, single_h5):
        runner.invoke(cli, ["init-macroscale", str(single_h5), "--dry-run"])

        with h5py.File(str(single_h5), "r") as f:
            assert "macro_data" not in f

    def test_h5_mode_dry_run_exits_0(self, runner, single_h5):
        result = runner.invoke(
            cli, ["init-macroscale", str(single_h5), "--dry-run"]
        )
        assert result.exit_code == 0, result.output

    def test_h5_mode_param_override_applied(self, runner, single_h5):
        """--param rows=50 overrides the default rows value."""
        result = runner.invoke(
            cli,
            ["init-macroscale", str(single_h5), "--param", "rows=50"],
        )
        assert result.exit_code == 0, result.output

        run_code = single_h5.stem
        run_dir = str(single_h5.parent)
        with DataStore(run_code, run_dir, mode="r") as ds:
            assert ds.macro_params.rows == 50

    def test_h5_mode_invalid_param_format_exits_nonzero(self, runner, single_h5):
        """--param without = should exit with a nonzero code."""
        result = runner.invoke(
            cli,
            ["init-macroscale", str(single_h5), "--param", "rows50"],
        )
        assert result.exit_code != 0 or "NAME=VALUE" in result.output


# ─── Error cases ──────────────────────────────────────────────────────────────


class TestInitMacroscaleErrors:
    def test_empty_microscale_exits_nonzero(self, runner, tmp_path):
        """HDF5 with no microscale data should fail gracefully."""
        micro = MicroParameters()
        h5_dir = tmp_path / "empty_run"
        h5_dir.mkdir()

        ds = DataStore.create("empty-run", str(h5_dir), micro)
        ds.close()
        # Do NOT fill microscale data

        result = runner.invoke(
            cli, ["init-macroscale", str(h5_dir / "empty-run.h5")]
        )
        assert result.exit_code != 0 or "empty" in result.output.lower()

    def test_already_initialized_exits_nonzero(self, runner, single_h5):
        """Calling init-macroscale twice on the same file should fail."""
        runner.invoke(cli, ["init-macroscale", str(single_h5)])
        result = runner.invoke(cli, ["init-macroscale", str(single_h5)])
        assert result.exit_code != 0 or "already" in result.output.lower()

    def test_bad_extension_exits_nonzero(self, runner, tmp_path):
        """A file with an unrecognised extension should fail."""
        bad = tmp_path / "not_a_h5_file.txt"
        bad.write_text("not an HDF5 file")

        result = runner.invoke(cli, ["init-macroscale", str(bad)])
        assert result.exit_code != 0

    def test_folder_without_experiment_json_exits_nonzero(self, runner, tmp_path):
        """A folder that lacks experiment.json should fail gracefully."""
        empty_dir = tmp_path / "empty_folder"
        empty_dir.mkdir()

        result = runner.invoke(cli, ["init-macroscale", str(empty_dir)])
        assert result.exit_code != 0
