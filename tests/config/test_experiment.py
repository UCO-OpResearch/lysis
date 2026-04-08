"""Unit tests for lysis.config.experiment.

Tests cover:
  - from_csv: basic two-row CSV creates correct folder structure and HDF5 files
  - from_csv: dry_run=True validates but creates no files
  - from_csv: experiment folder already exists → FileExistsError
  - from_csv: inconsistent CSV row → ValueError (ParameterConflict wrapped)
  - from_csv: empty CSV → ValueError
  - HDF5 round-trip: loaded params match what was written
  - run_code from CSV overrides auto-generated code
  - run_description stored in experiment.json
  - name defaults to CSV filename stem
"""

import csv
import json
import os
import textwrap

import pytest

from lysis.config.experiment import Experiment
from lysis.config.parameters import MacroParameters, MicroParameters
from lysis.dataio.datastore import DataStore


# ─── CSV fixtures ─────────────────────────────────────────────────────────────


def _write_minimal_csv(path, rows):
    """Write a CSV with headers from the first row dict and data rows."""
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _default_micro_row():
    """Return a minimal dict of micro CSV values (all independent, default values)."""
    mp = MicroParameters()
    d = mp.to_basedict()
    # Remove dependent params — we only need independent ones in the CSV
    import inspect
    micro_ind = set(inspect.signature(MicroParameters).parameters.keys())
    return {k: v for k, v in d.items() if k in micro_ind}


def _default_macro_row(micro=None):
    """Return a minimal dict of macro CSV values (all independent, default values)."""
    mp = micro or MicroParameters()
    d = MacroParameters(micro_params=mp).to_basedict()
    import inspect
    macro_ind = set(inspect.signature(MacroParameters).parameters.keys())
    macro_ind.discard("micro_params")
    return {k: v for k, v in d.items() if k in macro_ind}


def _two_row_csv(tmp_path):
    """Write a valid two-row CSV with slightly different total_molecules."""
    micro = _default_micro_row()
    macro1 = _default_macro_row()
    macro2 = dict(macro1)
    macro2["total_molecules"] = 86148  # double the default

    row1 = {**micro, **macro1, "run_code": "test-run-00", "run_description": "default"}
    row2 = {**micro, **macro2, "run_code": "test-run-01", "run_description": "double tPA"}
    csv_path = tmp_path / "test_exp.csv"
    _write_minimal_csv(csv_path, [row1, row2])
    return csv_path


# ─── Basic creation tests ─────────────────────────────────────────────────────


class TestFromCsvBasic:
    def test_experiment_folder_created(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        exp = Experiment.from_csv(csv_path, data_root, name="test-exp")
        assert os.path.isdir(exp.path)

    def test_experiment_json_written(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        exp = Experiment.from_csv(csv_path, data_root, name="test-exp")
        json_path = os.path.join(exp.path, "experiment.json")
        assert os.path.isfile(json_path)
        with open(json_path) as fh:
            meta = json.load(fh)
        assert meta["name"] == "test-exp"
        assert len(meta["runs"]) == 2

    def test_hdf5_files_created(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        exp = Experiment.from_csv(csv_path, data_root, name="test-exp")
        for run in exp.runs:
            h5 = os.path.join(exp.path, f"{run.run_code}.h5")
            assert os.path.isfile(h5), f"Missing HDF5: {h5}"

    def test_two_runs_created(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        exp = Experiment.from_csv(csv_path, data_root, name="test-exp")
        assert len(exp.runs) == 2

    def test_name_defaults_to_csv_stem(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        exp = Experiment.from_csv(csv_path, data_root)
        assert exp.name == "test_exp"

    def test_description_stored_in_json(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        exp = Experiment.from_csv(csv_path, data_root, name="x", description="my desc")
        json_path = os.path.join(exp.path, "experiment.json")
        with open(json_path) as fh:
            meta = json.load(fh)
        assert meta["description"] == "my desc"

    def test_run_code_from_csv(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        exp = Experiment.from_csv(csv_path, data_root, name="test-exp")
        codes = {run.run_code for run in exp.runs}
        assert "test-run-00" in codes
        assert "test-run-01" in codes

    def test_run_description_in_json(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        exp = Experiment.from_csv(csv_path, data_root, name="test-exp")
        json_path = os.path.join(exp.path, "experiment.json")
        with open(json_path) as fh:
            meta = json.load(fh)
        descs = {r["description"] for r in meta["runs"]}
        assert "default" in descs
        assert "double tPA" in descs


# ─── HDF5 content round-trip ──────────────────────────────────────────────────


class TestHdf5RoundTrip:
    def test_micro_params_readable(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        exp = Experiment.from_csv(csv_path, data_root, name="test-exp")
        run = exp.runs[0]
        ds = DataStore(run.run_code, exp.path, mode="r")
        loaded_micro = ds.micro_params
        ds.close()
        assert isinstance(loaded_micro, MicroParameters)
        assert loaded_micro.nodes_in_micro_row == run.micro_params.nodes_in_micro_row

    def test_macro_params_readable(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        exp = Experiment.from_csv(csv_path, data_root, name="test-exp")
        run0 = exp.runs[0]
        run1 = exp.runs[1]
        ds0 = DataStore(run0.run_code, exp.path, mode="r")
        ds1 = DataStore(run1.run_code, exp.path, mode="r")
        macro0 = ds0.macro_params
        macro1 = ds1.macro_params
        ds0.close()
        ds1.close()

        assert macro0.total_molecules == 43074
        assert macro1.total_molecules == 86148

    def test_params_match_written_values(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        exp = Experiment.from_csv(csv_path, data_root, name="test-exp")
        for run in exp.runs:
            ds = DataStore(run.run_code, exp.path, mode="r")
            loaded_micro = ds.micro_params
            loaded_macro = ds.macro_params
            ds.close()
            import math
            assert math.isclose(
                loaded_micro.fiber_radius.to("microns").magnitude,
                run.micro_params.fiber_radius.to("microns").magnitude,
                rel_tol=1e-6,
            )
            assert loaded_macro.cols == run.macro_params.cols
            assert loaded_macro.rows == run.macro_params.rows


# ─── dry_run ─────────────────────────────────────────────────────────────────


class TestDryRun:
    def test_dry_run_creates_no_folder(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        exp = Experiment.from_csv(csv_path, data_root, name="test-exp", dry_run=True)
        assert not os.path.exists(exp.path)

    def test_dry_run_populates_runs(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        exp = Experiment.from_csv(csv_path, data_root, name="test-exp", dry_run=True)
        assert len(exp.runs) == 2
        for run in exp.runs:
            assert isinstance(run.micro_params, MicroParameters)
            assert isinstance(run.macro_params, MacroParameters)

    def test_dry_run_then_real_run(self, tmp_path):
        """dry_run followed by real run should succeed (no folder created by dry_run)."""
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        Experiment.from_csv(csv_path, data_root, name="test-exp", dry_run=True)
        exp = Experiment.from_csv(csv_path, data_root, name="test-exp", dry_run=False)
        assert os.path.isdir(exp.path)


# ─── Error cases ─────────────────────────────────────────────────────────────


class TestErrorCases:
    def test_folder_already_exists(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()
        # Create the experiment folder in advance
        (data_root / "test-exp").mkdir()

        with pytest.raises(FileExistsError):
            Experiment.from_csv(csv_path, data_root, name="test-exp")

    def test_unrecognised_column_raises(self, tmp_path):
        micro = _default_micro_row()
        macro = _default_macro_row()
        row = {**micro, **macro, "not_a_param": "42"}
        csv_path = tmp_path / "bad.csv"
        _write_minimal_csv(csv_path, [row])
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        with pytest.raises(ValueError, match="Unrecognised CSV column"):
            Experiment.from_csv(csv_path, data_root)

    def test_empty_csv_raises(self, tmp_path):
        csv_path = tmp_path / "empty.csv"
        csv_path.write_text("fiber_radius\n")  # headers only, no data rows
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        with pytest.raises(ValueError, match="no data rows"):
            Experiment.from_csv(csv_path, data_root)

    def test_inconsistent_row_raises_with_row_info(self, tmp_path):
        micro = _default_micro_row()
        macro = _default_macro_row()
        # Make an inconsistent row: bind_rate * diss_const ≠ unbind_rate
        bad_micro = dict(micro)
        bad_micro["unbind_rate_tPA_woPLG"] = "9.99 1 / second"  # inconsistent
        row = {**bad_micro, **macro}
        csv_path = tmp_path / "bad.csv"
        _write_minimal_csv(csv_path, [row])
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        with pytest.raises(ValueError) as exc_info:
            Experiment.from_csv(csv_path, data_root)
        assert "Row 1" in str(exc_info.value)

    def test_all_row_errors_collected(self, tmp_path):
        """Errors in multiple rows are all reported at once."""
        micro = _default_micro_row()
        macro = _default_macro_row()
        bad_micro = dict(micro)
        bad_micro["unbind_rate_tPA_woPLG"] = "9.99 1 / second"
        row_bad = {**bad_micro, **macro}

        csv_path = tmp_path / "multi_bad.csv"
        _write_minimal_csv(csv_path, [row_bad, row_bad])
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        with pytest.raises(ValueError) as exc_info:
            Experiment.from_csv(csv_path, data_root)
        msg = str(exc_info.value)
        assert "Row 1" in msg
        assert "Row 2" in msg


# ─── Parameter inversion via CSV ─────────────────────────────────────────────


class TestInversionViaCsv:
    """Providing dependent micro params in CSV → solver fills in independent."""

    def test_unbind_rate_instead_of_diss_const(self, tmp_path):
        """Specify unbind_rate + bind_rate → diss_const solved automatically."""
        import inspect

        micro_ind = set(inspect.signature(MicroParameters).parameters.keys())
        # Start from default independent params, but remove diss_const_tPA_woPLG
        # and add unbind_rate_tPA_woPLG instead
        micro = _default_micro_row()
        del micro["diss_const_tPA_woPLG"]
        micro["unbind_rate_tPA_woPLG"] = "0.072 1 / second"
        # bind_rate_tPA is still present (default 0.1)
        # → diss_const_tPA_woPLG should be solved as 0.72 µM

        macro = _default_macro_row()
        row = {**micro, **macro}
        csv_path = tmp_path / "inv.csv"
        _write_minimal_csv(csv_path, [row])
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        exp = Experiment.from_csv(csv_path, data_root, name="inv-exp", dry_run=True)
        loaded = exp.runs[0].micro_params
        import math

        assert math.isclose(
            loaded.unbind_rate_tPA_woPLG.to("1/second").magnitude,
            0.072,
            rel_tol=1e-5,
        )
        assert math.isclose(
            loaded.diss_const_tPA_woPLG.to("micromolar").magnitude,
            0.72,
            rel_tol=1e-5,
        )
