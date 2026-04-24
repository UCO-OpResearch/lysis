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

import h5py
import pytest

from lysis.config.constants import CONST
from lysis.config.experiment import Experiment
from lysis.config.parameters import MacroParameters, MicroParameters
from lysis.dataio.datastore import DataStore


# ─── CSV fixtures ─────────────────────────────────────────────────────────────


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
    def test_hdf5_has_only_microscale_collection(self, tmp_path):
        """from_csv() produces HDF5 with microscale_out only; macro added later.

        initialize_macroscale() is deferred until after microscale Simulations
        complete so that forced_unbind can be computed from the results.
        """
        import h5py
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        exp = Experiment.from_csv(csv_path, data_root, name="test-exp")
        run = exp.runs[0]
        h5_path = os.path.join(exp.path, f"{run.run_code}.h5")
        with h5py.File(h5_path, "r") as f:
            assert "micro_data" in f
            assert "macro_data" not in f

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

    def test_macro_params_in_run(self, tmp_path):
        """from_csv() parses macro params into Run objects (not yet in HDF5).

        initialize_macroscale() is not called by from_csv(), so macro_data
        is absent from the HDF5 at this stage. The macro parameters live in
        run.macro_params until initialize_macroscale() is called after
        microscale Simulations complete.
        """
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        exp = Experiment.from_csv(csv_path, data_root, name="test-exp")
        run0 = exp.runs[0]
        run1 = exp.runs[1]

        assert isinstance(run0.macro_params, MacroParameters)
        assert isinstance(run1.macro_params, MacroParameters)
        assert run0.macro_params.total_molecules == 43074
        assert run1.macro_params.total_molecules == 86148

    def test_micro_params_match_written_values(self, tmp_path):
        """HDF5 micro_params round-trip: loaded values match what was written.

        Macro params are not in the HDF5 after from_csv() (initialize_macroscale()
        has not been called). Check micro params from HDF5 and macro params
        from the Run object instead.
        """
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()

        exp = Experiment.from_csv(csv_path, data_root, name="test-exp")
        import math
        for run in exp.runs:
            ds = DataStore(run.run_code, exp.path, mode="r")
            loaded_micro = ds.micro_params
            ds.close()
            assert math.isclose(
                loaded_micro.fiber_radius.to("microns").magnitude,
                run.micro_params.fiber_radius.to("microns").magnitude,
                rel_tol=1e-6,
            )
            # Macro params live in the Run object until initialize_macroscale()
            # is called after microscale Simulations complete.
            assert run.macro_params.cols == MacroParameters(
                micro_params=run.micro_params
            ).cols
            assert run.macro_params.rows == MacroParameters(
                micro_params=run.micro_params
            ).rows


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
        # Transposed format: header row has run column, but no parameter rows follow
        csv_path.write_text("parameter,run-00\n")
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
        assert "Run 1" in str(exc_info.value)

    def test_all_row_errors_collected(self, tmp_path):
        """Errors in multiple runs are all reported at once."""
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
        assert "Run 1" in msg
        assert "Run 2" in msg


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


# ─── Rename tests ────────────────────────────────────────────────────────────


class TestRunRename:
    def test_rename_renames_file_and_records_history(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()
        exp = Experiment.from_csv(csv_path, data_root, name="test-exp")

        run = exp.runs[0]
        old_code = run.run_code
        old_h5 = os.path.join(exp.path, f"{old_code}.h5")
        new_h5 = os.path.join(exp.path, "new-code.h5")
        assert os.path.isfile(old_h5)

        returned = run.rename("new-code")

        assert returned == old_code
        assert run.run_code == "new-code"
        assert not os.path.exists(old_h5)
        assert os.path.isfile(new_h5)

        with h5py.File(new_h5, "r") as f:
            assert f.attrs[CONST.RENAMED_FROM_ATTR] == old_code

    def test_rename_appends_history(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()
        exp = Experiment.from_csv(csv_path, data_root, name="test-exp")

        run = exp.runs[0]
        original = run.run_code
        run.rename("second")
        run.rename("third")

        final_h5 = os.path.join(exp.path, "third.h5")
        with h5py.File(final_h5, "r") as f:
            assert f.attrs[CONST.RENAMED_FROM_ATTR] == f"{original} -> second"

    def test_rename_to_existing_filename_raises(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()
        exp = Experiment.from_csv(csv_path, data_root, name="test-exp")

        other_code = exp.runs[1].run_code
        with pytest.raises(FileExistsError):
            exp.runs[0].rename(other_code)

    def test_rename_same_name_raises(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()
        exp = Experiment.from_csv(csv_path, data_root, name="test-exp")

        run = exp.runs[0]
        with pytest.raises(ValueError):
            run.rename(run.run_code)


class TestExperimentRename:
    def test_rename_folder_and_json(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()
        exp = Experiment.from_csv(csv_path, data_root, name="old-exp")
        old_path = exp.path

        exp.rename("new-exp")

        assert not os.path.exists(old_path)
        assert os.path.isdir(exp.path)
        assert exp.name == "new-exp"

        with open(os.path.join(exp.path, "experiment.json")) as fh:
            meta = json.load(fh)
        assert meta["name"] == "new-exp"

        for run in exp.runs:
            assert run.os_path == exp.path
            assert os.path.isfile(os.path.join(exp.path, f"{run.run_code}.h5"))

    def test_rename_then_load(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()
        exp = Experiment.from_csv(csv_path, data_root, name="old-exp")
        exp.rename("new-exp")

        reloaded = Experiment.load(exp.path)
        assert reloaded.name == "new-exp"
        assert len(reloaded.runs) == len(exp.runs)

    def test_rename_existing_target_raises(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()
        exp = Experiment.from_csv(csv_path, data_root, name="old-exp")
        (data_root / "taken").mkdir()

        with pytest.raises(FileExistsError):
            exp.rename("taken")


class TestExperimentRenameRun:
    def test_rename_run_updates_json_and_file(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()
        exp = Experiment.from_csv(csv_path, data_root, name="test-exp")

        old_code = exp.runs[0].run_code
        exp.rename_run(old_code, "brand-new")

        assert exp.runs[0].run_code == "brand-new"
        assert os.path.isfile(os.path.join(exp.path, "brand-new.h5"))
        assert not os.path.exists(os.path.join(exp.path, f"{old_code}.h5"))

        with open(os.path.join(exp.path, "experiment.json")) as fh:
            meta = json.load(fh)
        run_codes = [r["run_code"] for r in meta["runs"]]
        assert "brand-new" in run_codes
        assert old_code not in run_codes

        with h5py.File(os.path.join(exp.path, "brand-new.h5"), "r") as f:
            assert f.attrs[CONST.RENAMED_FROM_ATTR] == old_code

    def test_rename_run_unknown_code_raises(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()
        exp = Experiment.from_csv(csv_path, data_root, name="test-exp")

        with pytest.raises(KeyError):
            exp.rename_run("not-a-real-run", "anything")

    def test_rename_run_collision_raises(self, tmp_path):
        csv_path = _two_row_csv(tmp_path)
        data_root = tmp_path / "experiments"
        data_root.mkdir()
        exp = Experiment.from_csv(csv_path, data_root, name="test-exp")

        existing = exp.runs[1].run_code
        with pytest.raises(ValueError):
            exp.rename_run(exp.runs[0].run_code, existing)
