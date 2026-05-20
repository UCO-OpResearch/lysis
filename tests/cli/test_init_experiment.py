"""Tests for ``lysis init-experiment`` CLI command."""

import csv
import json
import os

import pytest
from click.testing import CliRunner

from lysis.cli import cli
from lysis.config.parameters import MacroParameters, MicroParameters


# ─── CSV helpers (shared with test_experiment.py) ────────────────────────────


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
    import inspect
    mp = MicroParameters()
    d = mp.to_basedict()
    micro_ind = set(inspect.signature(MicroParameters).parameters.keys())
    return {k: v for k, v in d.items() if k in micro_ind}


def _default_macro_row(micro=None):
    import inspect
    mp = micro or MicroParameters()
    d = MacroParameters(micro_params=mp).to_basedict()
    macro_ind = set(inspect.signature(MacroParameters).parameters.keys())
    macro_ind.discard("micro_params")
    return {k: v for k, v in d.items() if k in macro_ind}


def _two_row_csv(tmp_path):
    micro = _default_micro_row()
    macro1 = _default_macro_row()
    macro2 = dict(macro1)
    macro2["total_molecules"] = 86148
    row1 = {**micro, **macro1, "run_code": "test-run-00", "run_description": "default"}
    row2 = {**micro, **macro2, "run_code": "test-run-01", "run_description": "double tPA"}
    csv_path = tmp_path / "test_exp.csv"
    _write_minimal_csv(csv_path, [row1, row2])
    return csv_path


@pytest.fixture
def runner():
    return CliRunner()


# ─── Basic invocation ─────────────────────────────────────────────────────────


class TestInitExperimentBasic:
    def test_help_exits_0(self, runner):
        result = runner.invoke(cli, ["init-experiment", "--help"])
        assert result.exit_code == 0
        assert "CSV_PATH" in result.output

    def test_creates_experiment_folder(self, runner, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        result = runner.invoke(
            cli,
            ["init-experiment", str(csv_path), str(data_root), "--name", "cli-exp"],
        )
        assert result.exit_code == 0, result.output
        assert os.path.isdir(str(data_root / "cli-exp"))

    def test_creates_hdf5_files(self, runner, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        result = runner.invoke(
            cli,
            ["init-experiment", str(csv_path), str(data_root), "--name", "cli-exp"],
        )
        assert result.exit_code == 0, result.output
        exp_dir = data_root / "cli-exp"
        assert os.path.isfile(str(exp_dir / "test-run-00.h5"))
        assert os.path.isfile(str(exp_dir / "test-run-01.h5"))

    def test_creates_experiment_json(self, runner, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        result = runner.invoke(
            cli,
            ["init-experiment", str(csv_path), str(data_root), "--name", "cli-exp"],
        )
        assert result.exit_code == 0, result.output
        json_path = data_root / "cli-exp" / "experiment.json"
        assert json_path.exists()
        with open(json_path) as fh:
            meta = json.load(fh)
        assert meta["name"] == "cli-exp"
        assert len(meta["runs"]) == 2

    def test_output_includes_experiment_path(self, runner, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        result = runner.invoke(
            cli,
            ["init-experiment", str(csv_path), str(data_root), "--name", "cli-exp"],
        )
        assert result.exit_code == 0, result.output
        # Rich may wrap long paths; check that "Experiment created" and "cli-exp" both appear
        # anywhere in the output (possibly across lines)
        assert "Experiment created" in result.output
        assert "cli-exp" in result.output.replace("\n", "")

    def test_no_progress_flag(self, runner, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        result = runner.invoke(
            cli,
            [
                "init-experiment",
                str(csv_path),
                str(data_root),
                "--name",
                "cli-exp",
                "--no-progress",
            ],
        )
        assert result.exit_code == 0, result.output

    def test_description_stored(self, runner, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        result = runner.invoke(
            cli,
            [
                "init-experiment",
                str(csv_path),
                str(data_root),
                "--name",
                "cli-exp",
                "--description",
                "my description",
            ],
        )
        assert result.exit_code == 0, result.output
        json_path = data_root / "cli-exp" / "experiment.json"
        with open(json_path) as fh:
            meta = json.load(fh)
        assert meta["description"] == "my description"


# ─── Dry run ──────────────────────────────────────────────────────────────────


class TestInitExperimentDryRun:
    def test_dry_run_creates_no_folder(self, runner, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        result = runner.invoke(
            cli,
            ["init-experiment", str(csv_path), str(data_root), "--name", "cli-exp", "--dry-run"],
        )
        assert result.exit_code == 0, result.output
        assert not os.path.exists(str(data_root / "cli-exp"))

    def test_dry_run_output_says_dry_run(self, runner, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        result = runner.invoke(
            cli,
            ["init-experiment", str(csv_path), str(data_root), "--dry-run"],
        )
        assert result.exit_code == 0, result.output
        assert "Dry run" in result.output or "dry run" in result.output.lower()

    def test_dry_run_shows_run_codes(self, runner, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        result = runner.invoke(
            cli,
            ["init-experiment", str(csv_path), str(data_root), "--dry-run"],
        )
        assert result.exit_code == 0, result.output
        assert "test-run-00" in result.output
        assert "test-run-01" in result.output


# ─── Error cases ─────────────────────────────────────────────────────────────


class TestInitExperimentErrors:
    def test_folder_already_exists_exits_nonzero(self, runner, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()
        (data_root / "cli-exp").mkdir()

        result = runner.invoke(
            cli,
            ["init-experiment", str(csv_path), str(data_root), "--name", "cli-exp"],
        )
        assert result.exit_code != 0

    def test_folder_exists_suggests_name_option(self, runner, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()
        (data_root / "cli-exp").mkdir()

        result = runner.invoke(
            cli,
            ["init-experiment", str(csv_path), str(data_root), "--name", "cli-exp"],
        )
        assert "--name" in result.output

    def test_bad_csv_column_exits_nonzero(self, runner, tmp_path):
        micro = _default_micro_row()
        macro = _default_macro_row()
        row = {**micro, **macro, "not_a_param": "42"}
        csv_path = tmp_path / "bad.csv"
        _write_minimal_csv(csv_path, [row])
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        result = runner.invoke(
            cli,
            ["init-experiment", str(csv_path), str(data_root)],
        )
        assert result.exit_code != 0

    def test_inconsistent_row_exits_nonzero(self, runner, tmp_path):
        micro = _default_micro_row()
        macro = _default_macro_row()
        bad_micro = dict(micro)
        bad_micro["unbind_rate_tPA_woPLG"] = "9.99 1 / second"
        row = {**bad_micro, **macro}
        csv_path = tmp_path / "bad.csv"
        _write_minimal_csv(csv_path, [row])
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        result = runner.invoke(
            cli,
            ["init-experiment", str(csv_path), str(data_root)],
        )
        assert result.exit_code != 0
