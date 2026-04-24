"""Tests for ``lysis rename`` CLI command."""

import csv
import json
import os

import h5py
import pytest
from click.testing import CliRunner

from lysis.cli import cli
from lysis.config.constants import CONST
from lysis.config.experiment import Experiment
from lysis.config.parameters import MacroParameters, MicroParameters


def _write_minimal_csv(path, rows):
    if not rows:
        return
    param_names = list(rows[0].keys())
    run_codes = [r.get("run_code", f"run-{i:02d}") for i, r in enumerate(rows)]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["parameter"] + run_codes)
        for param in param_names:
            if param == "run_code":
                continue
            writer.writerow([param] + [str(r.get(param, "")) for r in rows])


def _default_micro_row():
    import inspect

    mp = MicroParameters()
    d = mp.to_basedict()
    micro_ind = set(inspect.signature(MicroParameters).parameters.keys())
    return {k: v for k, v in d.items() if k in micro_ind}


def _default_macro_row():
    import inspect

    mp = MicroParameters()
    d = MacroParameters(micro_params=mp).to_basedict()
    macro_ind = set(inspect.signature(MacroParameters).parameters.keys())
    macro_ind.discard("micro_params")
    return {k: v for k, v in d.items() if k in macro_ind}


def _two_row_csv(tmp_path):
    micro = _default_micro_row()
    macro1 = _default_macro_row()
    macro2 = dict(macro1)
    macro2["total_molecules"] = 86148
    row1 = {**micro, **macro1, "run_code": "test-run-00"}
    row2 = {**micro, **macro2, "run_code": "test-run-01"}
    csv_path = tmp_path / "test_exp.csv"
    _write_minimal_csv(csv_path, [row1, row2])
    return csv_path


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def experiment(tmp_path):
    csv_path = _two_row_csv(tmp_path)
    data_root = tmp_path / "experiments"
    data_root.mkdir()
    return Experiment.from_csv(csv_path, data_root, name="cli-exp")


class TestRenameCli:
    def test_help_exits_0(self, runner):
        result = runner.invoke(cli, ["rename", "--help"])
        assert result.exit_code == 0
        assert "PATH" in result.output
        assert "NEW_NAME" in result.output

    def test_rename_folder(self, runner, experiment):
        old_path = experiment.path
        parent = os.path.dirname(old_path)

        result = runner.invoke(cli, ["rename", old_path, "renamed"])

        assert result.exit_code == 0, result.output
        assert not os.path.exists(old_path)
        new_path = os.path.join(parent, "renamed")
        assert os.path.isdir(new_path)
        with open(os.path.join(new_path, "experiment.json")) as fh:
            assert json.load(fh)["name"] == "renamed"

    def test_rename_h5_updates_experiment_json(self, runner, experiment):
        old_code = experiment.runs[0].run_code
        h5_path = os.path.join(experiment.path, f"{old_code}.h5")

        result = runner.invoke(cli, ["rename", h5_path, "new-run"])

        assert result.exit_code == 0, result.output
        new_h5 = os.path.join(experiment.path, "new-run.h5")
        assert os.path.isfile(new_h5)
        assert not os.path.exists(h5_path)
        with h5py.File(new_h5, "r") as f:
            assert f.attrs[CONST.RENAMED_FROM_ATTR] == old_code
        with open(os.path.join(experiment.path, "experiment.json")) as fh:
            codes = [r["run_code"] for r in json.load(fh)["runs"]]
        assert "new-run" in codes
        assert old_code not in codes

    def test_rename_h5_strips_extension_in_new_name(self, runner, experiment):
        old_code = experiment.runs[0].run_code
        h5_path = os.path.join(experiment.path, f"{old_code}.h5")

        result = runner.invoke(cli, ["rename", h5_path, "new-run.h5"])

        assert result.exit_code == 0, result.output
        assert os.path.isfile(os.path.join(experiment.path, "new-run.h5"))

    def test_rename_h5_without_experiment_json(self, runner, experiment):
        # Remove experiment.json so we exercise the standalone-run path.
        os.remove(os.path.join(experiment.path, "experiment.json"))
        old_code = experiment.runs[0].run_code
        h5_path = os.path.join(experiment.path, f"{old_code}.h5")

        result = runner.invoke(cli, ["rename", h5_path, "standalone"])

        assert result.exit_code == 0, result.output
        new_h5 = os.path.join(experiment.path, "standalone.h5")
        assert os.path.isfile(new_h5)
        with h5py.File(new_h5, "r") as f:
            assert f.attrs[CONST.RENAMED_FROM_ATTR] == old_code

    def test_rename_rejects_path_separator_in_new_name(self, runner, experiment):
        result = runner.invoke(cli, ["rename", experiment.path, "a/b"])
        assert result.exit_code != 0
        assert "path separators" in result.output

    def test_rename_unknown_path_type(self, runner, tmp_path):
        txt = tmp_path / "not-h5.txt"
        txt.write_text("hello")
        result = runner.invoke(cli, ["rename", str(txt), "whatever"])
        assert result.exit_code != 0
        assert "experiment folder or an HDF5" in result.output
