"""Unit tests for :mod:`lysis.execution.codeutil` — FortranMicro class.

Tests cover:

* :class:`FortranMicro` command construction (``exec_command``)
* Construction from HDF5 (``from_hdf5``)
* Execution in a working directory (``exec_in_workdir``)
* Result import into HDF5 (``import_results``)
* Full end-to-end workflow wrapper (``run_full``)
* Integration test with real fixture data (``TestImportResultsWithFixture``)
"""

import json
import os
import shutil
from pathlib import Path
from unittest.mock import MagicMock, call, patch

import h5py
import numpy as np
import pytest

from lysis.config.constants import CONST
from lysis.config.parameters import MicroParameters
from lysis.config.run import Run
from lysis.execution.codeutil import FortranMicro, MICRO_FORTRAN_DATASPEC_VERSION


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------


def _write_micro_hdf5(path: Path, mp: MicroParameters = None) -> None:
    """Write a minimal v2.0.0 HDF5 file with MicroParameters attributes."""
    mp = mp or MicroParameters()
    with h5py.File(str(path), "w") as f:
        f.attrs[CONST.DATASPEC_VERSION_ATTR] = "v2.0.0"
        grp = f.require_group("micro_data")
        for k, v in mp.to_basedict().items():
            grp.attrs[k] = str(v) if not isinstance(v, (int, float, bool)) else v


@pytest.fixture
def tmp_run(tmp_path):
    """A Run with default MicroParameters in a tmp directory."""
    r = Run(str(tmp_path))
    r.initialize_micro_param()
    return r


@pytest.fixture
def micro_hdf5(tmp_path):
    """A minimal HDF5 file with default MicroParameters; returns Path."""
    h5_path = tmp_path / "test-run.h5"
    _write_micro_hdf5(h5_path)
    return h5_path


@pytest.fixture
def fortran_micro(tmp_run):
    """A FortranMicro with default run and dummy executable."""
    return FortranMicro(run=tmp_run, executable="/bin/micro.exe")


# ---------------------------------------------------------------------------
# TestFortranMicroExecCommand
# ---------------------------------------------------------------------------


class TestFortranMicroExecCommand:
    """Tests for :meth:`FortranMicro.exec_command`."""

    def test_default_params_contains_run_code(self, tmp_run):
        """Command must include --runCode and --outFileCode."""
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe")
        cmd = fm.exec_command()
        assert "--runCode" in cmd
        assert tmp_run.run_code in cmd

    def test_default_params_no_extra_flags(self, tmp_run):
        """With all-default params the command should only have --runCode and --outFileCode."""
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe")
        cmd = fm.exec_command()
        # Only executable + 4 elements: --runCode, run_code, --outFileCode, ""
        assert len(cmd) == 5

    def test_non_default_param_included(self, tmp_path):
        """A non-default constructor param must appear in the command."""
        r = Run(str(tmp_path))
        r.initialize_micro_param({"micro_simulations": 999})
        fm = FortranMicro(run=r, executable="/bin/micro.exe")
        cmd = fm.exec_command()
        assert "--simulations" in cmd
        assert "999" in cmd

    def test_default_param_excluded(self, tmp_run):
        """A param equal to its default must not appear in the command."""
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe")
        cmd = fm.exec_command()
        # micro_simulations == default so --simulations must be absent
        assert "--simulations" not in cmd

    def test_unit_bearing_param_converted(self, tmp_path):
        """fiber_radius (a Pint Quantity) must be output in the expected SI unit."""
        from lysis.config.constants import Q_

        r = Run(str(tmp_path))
        r.initialize_micro_param({"fiber_radius": Q_("37 nanometer")})
        fm = FortranMicro(run=r, executable="/bin/micro.exe")
        cmd = fm.exec_command()
        assert "--radius" in cmd
        idx = cmd.index("--radius")
        # The Fortran code expects radius in microns; 37 nm = 0.037 µm
        val = float(cmd[idx + 1])
        assert pytest.approx(val, rel=1e-6) == 0.037

    def test_index_none_no_seed_split(self, tmp_run):
        """With index=None the seed must not be split and micro_simulations unchanged."""
        default_sims = tmp_run.micro_params.micro_simulations
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe", index=None)
        cmd = fm.exec_command()
        # micro_simulations should still be absent (equal to default)
        assert "--simulations" not in cmd
        assert tmp_run.micro_params.micro_simulations == default_sims

    def test_index_none_out_file_code_unchanged(self, tmp_run):
        """With index=None out_file_code must remain empty."""
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe",
                          out_file_code="suffix", index=None)
        fm.exec_command()
        assert fm.out_file_code == "suffix"

    def test_index_zero_sets_simulations_to_one(self, tmp_run):
        """index=0 must force micro_simulations to 1 in the command."""
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe", index=0)
        cmd = fm.exec_command()
        assert "--simulations" in cmd
        assert cmd[cmd.index("--simulations") + 1] == "1"

    def test_index_zero_appends_suffix(self, tmp_run):
        """index=0 must append '__00' to out_file_code."""
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe", index=0)
        fm.exec_command()
        assert fm.out_file_code == "__00"

    def test_index_two_appends_correct_suffix(self, tmp_run):
        """index=2 must append '__02' to out_file_code."""
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe", index=2)
        fm.exec_command()
        assert fm.out_file_code == "__02"

    def test_index_two_uses_correct_seed(self, tmp_run):
        """index=2 must derive seed from SeedSequence at position 2."""
        seed = tmp_run.micro_params.micro_seed
        stream = np.random.SeedSequence(seed)
        expected_seed = int(np.int32(stream.generate_state(3)[2]))

        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe", index=2)
        cmd = fm.exec_command()
        assert "--seed" in cmd
        assert cmd[cmd.index("--seed") + 1] == str(expected_seed)

    def test_executable_is_first_element(self, tmp_run):
        """The executable must be the first element of the command list."""
        fm = FortranMicro(run=tmp_run, executable="/path/to/micro.exe")
        cmd = fm.exec_command()
        assert cmd[0] == "/path/to/micro.exe"


# ---------------------------------------------------------------------------
# TestFortranMicroFromHdf5
# ---------------------------------------------------------------------------


class TestFortranMicroFromHdf5:
    """Tests for :meth:`FortranMicro.from_hdf5`."""

    def test_constructs_run_with_correct_run_code(self, micro_hdf5):
        fm = FortranMicro.from_hdf5(micro_hdf5, "/bin/micro.exe")
        assert fm.run.run_code == micro_hdf5.stem

    def test_loads_micro_params(self, micro_hdf5):
        fm = FortranMicro.from_hdf5(micro_hdf5, "/bin/micro.exe")
        assert fm.run.micro_params is not None
        assert isinstance(fm.run.micro_params, MicroParameters)

    def test_sets_executable(self, micro_hdf5):
        fm = FortranMicro.from_hdf5(micro_hdf5, "/bin/micro.exe")
        assert fm.executable == "/bin/micro.exe"

    def test_passes_out_file_code(self, micro_hdf5):
        fm = FortranMicro.from_hdf5(micro_hdf5, "/bin/micro.exe",
                                     out_file_code="_code")
        assert fm.out_file_code == "_code"

    def test_passes_index(self, micro_hdf5):
        fm = FortranMicro.from_hdf5(micro_hdf5, "/bin/micro.exe", index=3)
        assert fm.index == 3

    def test_raises_for_missing_file(self, tmp_path):
        missing = tmp_path / "nonexistent.h5"
        with pytest.raises((RuntimeError, OSError)):
            FortranMicro.from_hdf5(missing, "/bin/micro.exe")

    def test_accepts_path_or_str(self, micro_hdf5):
        fm_path = FortranMicro.from_hdf5(micro_hdf5, "/bin/micro.exe")
        fm_str = FortranMicro.from_hdf5(str(micro_hdf5), "/bin/micro.exe")
        assert fm_path.run.run_code == fm_str.run.run_code


# ---------------------------------------------------------------------------
# TestFortranMicroExecInWorkdir
# ---------------------------------------------------------------------------


class TestFortranMicroExecInWorkdir:
    """Tests for :meth:`FortranMicro.exec_in_workdir`."""

    def test_creates_data_run_code_dir(self, fortran_micro, tmp_path):
        with patch("subprocess.run"):
            fortran_micro.exec_in_workdir(tmp_path)
        expected = tmp_path / "data" / fortran_micro.run.run_code
        assert expected.is_dir()

    def test_returns_data_dir_path(self, fortran_micro, tmp_path):
        with patch("subprocess.run"):
            result = fortran_micro.exec_in_workdir(tmp_path)
        expected = tmp_path / "data" / fortran_micro.run.run_code
        assert result == expected

    def test_log_file_in_data_dir(self, fortran_micro, tmp_path):
        """stdout is redirected to micro{out_code}.txt inside data_dir."""
        with patch("subprocess.run") as mock_run:
            data_dir = fortran_micro.exec_in_workdir(tmp_path)
        log_path = data_dir / "micro.txt"
        # Check that open was called with the log file path (it's used as context
        # manager, so subprocess.run itself is only called once)
        assert mock_run.call_count == 1
        # stdout kwarg should be an open file object; we verify by checking
        # that the expected log file location exists in the call kwargs
        call_kwargs = mock_run.call_args.kwargs
        assert "cwd" in call_kwargs
        assert call_kwargs["cwd"] == str(tmp_path)

    def test_subprocess_run_called_once(self, fortran_micro, tmp_path):
        with patch("subprocess.run") as mock_run:
            fortran_micro.exec_in_workdir(tmp_path)
        assert mock_run.call_count == 1

    def test_subprocess_cwd_is_work_dir(self, fortran_micro, tmp_path):
        with patch("subprocess.run") as mock_run:
            fortran_micro.exec_in_workdir(tmp_path)
        call_kwargs = mock_run.call_args.kwargs
        assert call_kwargs["cwd"] == str(tmp_path)

    def test_params_json_written(self, fortran_micro, tmp_path):
        """params.json must be present in data_dir for v1.99.0 import pipeline."""
        with patch("subprocess.run"):
            data_dir = fortran_micro.exec_in_workdir(tmp_path)
        assert (data_dir / "params.json").exists()

    def test_accepts_str_work_dir(self, fortran_micro, tmp_path):
        with patch("subprocess.run"):
            result = fortran_micro.exec_in_workdir(str(tmp_path))
        assert isinstance(result, Path)


# ---------------------------------------------------------------------------
# TestFortranMicroImportResults
# ---------------------------------------------------------------------------


class TestFortranMicroImportResults:
    """Tests for :meth:`FortranMicro.import_results`."""

    def _mock_converted(self):
        """Return a minimal converted data dict with one array."""
        return {"tpa_leaving_time": np.zeros(10, dtype=np.float32)}

    @patch("lysis.dataio.dataconvert.convert_data")
    @patch("lysis.dataio.fileops.read_data_collection")
    def test_uses_module_constant_not_literal(self, mock_read, mock_convert, tmp_path):
        """import_results must use MICRO_FORTRAN_DATASPEC_VERSION, not 'v1.99.0'."""
        mock_read.return_value = {}
        mock_convert.return_value = {}

        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        hdf5_path = tmp_path / "run.h5"

        with h5py.File(str(hdf5_path), "w") as _:
            pass  # empty file

        FortranMicro.import_results(data_dir, hdf5_path, keep_tmpdir=True)

        # read_data_collection should have been called with the spec for the
        # module constant version
        assert mock_read.call_count == 1
        from lysis.dataio.dataspec import dataspec
        expected_spec = dataspec[MICRO_FORTRAN_DATASPEC_VERSION]["microscale_out"]
        actual_specs = mock_read.call_args.args[1]
        assert actual_specs[0] == expected_spec

    @patch("lysis.dataio.dataconvert.convert_data")
    @patch("lysis.dataio.fileops.read_data_collection")
    def test_cleanup_on_success(self, mock_read, mock_convert, tmp_path):
        """data_dir must be removed after successful import (keep_tmpdir=False)."""
        mock_read.return_value = {}
        mock_convert.return_value = {}

        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        hdf5_path = tmp_path / "run.h5"
        with h5py.File(str(hdf5_path), "w") as _:
            pass

        FortranMicro.import_results(data_dir, hdf5_path, keep_tmpdir=False)
        assert not data_dir.exists()

    @patch("lysis.dataio.dataconvert.convert_data")
    @patch("lysis.dataio.fileops.read_data_collection")
    def test_keep_tmpdir_preserves_on_success(self, mock_read, mock_convert, tmp_path):
        """data_dir must be preserved when keep_tmpdir=True."""
        mock_read.return_value = {}
        mock_convert.return_value = {}

        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        hdf5_path = tmp_path / "run.h5"
        with h5py.File(str(hdf5_path), "w") as _:
            pass

        FortranMicro.import_results(data_dir, hdf5_path, keep_tmpdir=True)
        assert data_dir.exists()

    @patch("lysis.dataio.dataconvert.convert_data")
    @patch("lysis.dataio.fileops.read_data_collection")
    def test_cleanup_on_failure_keep_on_failure_false(self, mock_read, mock_convert, tmp_path):
        """data_dir must be removed on failure when keep_on_failure=False."""
        mock_read.side_effect = RuntimeError("read failed")

        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        hdf5_path = tmp_path / "run.h5"
        with h5py.File(str(hdf5_path), "w") as _:
            pass

        with pytest.raises(RuntimeError, match="read failed"):
            FortranMicro.import_results(
                data_dir, hdf5_path, keep_on_failure=False, keep_tmpdir=False
            )
        assert not data_dir.exists()

    @patch("lysis.dataio.dataconvert.convert_data")
    @patch("lysis.dataio.fileops.read_data_collection")
    def test_preserve_on_failure_keep_on_failure_true(self, mock_read, mock_convert, tmp_path):
        """data_dir must be preserved on failure when keep_on_failure=True."""
        mock_read.side_effect = RuntimeError("read failed")

        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        hdf5_path = tmp_path / "run.h5"
        with h5py.File(str(hdf5_path), "w") as _:
            pass

        with pytest.raises(RuntimeError):
            FortranMicro.import_results(
                data_dir, hdf5_path, keep_on_failure=True, keep_tmpdir=False
            )
        assert data_dir.exists()

    @patch("lysis.dataio.dataconvert.convert_data")
    @patch("lysis.dataio.fileops.read_data_collection")
    def test_keep_tmpdir_preserves_on_failure(self, mock_read, mock_convert, tmp_path):
        """data_dir must be preserved on failure when keep_tmpdir=True."""
        mock_read.side_effect = RuntimeError("read failed")

        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        hdf5_path = tmp_path / "run.h5"
        with h5py.File(str(hdf5_path), "w") as _:
            pass

        with pytest.raises(RuntimeError):
            FortranMicro.import_results(
                data_dir, hdf5_path, keep_on_failure=False, keep_tmpdir=True
            )
        assert data_dir.exists()

    @patch("lysis.dataio.dataconvert.convert_data")
    @patch("lysis.dataio.fileops.read_data_collection")
    def test_passes_file_code_to_read(self, mock_read, mock_convert, tmp_path):
        """file_code must be forwarded to read_data_collection."""
        mock_read.return_value = {}
        mock_convert.return_value = {}

        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        hdf5_path = tmp_path / "run.h5"
        with h5py.File(str(hdf5_path), "w") as _:
            pass

        FortranMicro.import_results(data_dir, hdf5_path, file_code="_code",
                                    keep_tmpdir=True)
        actual_codes = mock_read.call_args.args[2]
        assert actual_codes == ["_code"]


# ---------------------------------------------------------------------------
# TestFortranMicroRunFull
# ---------------------------------------------------------------------------


class TestFortranMicroRunFull:
    """Tests for :meth:`FortranMicro.run_full`."""

    def test_exec_in_workdir_called(self, fortran_micro, tmp_path):
        hdf5_path = tmp_path / "run.h5"
        with h5py.File(str(hdf5_path), "w") as _:
            pass
        with (
            patch.object(FortranMicro, "exec_in_workdir",
                         return_value=tmp_path / "data_dir") as mock_exec,
            patch.object(FortranMicro, "import_results") as mock_import,
        ):
            fortran_micro.run_full(hdf5_path)
        assert mock_exec.call_count == 1

    def test_import_results_called_with_data_dir(self, fortran_micro, tmp_path):
        hdf5_path = tmp_path / "run.h5"
        data_dir = tmp_path / "data_dir"
        with h5py.File(str(hdf5_path), "w") as _:
            pass
        with (
            patch.object(FortranMicro, "exec_in_workdir",
                         return_value=data_dir),
            patch.object(FortranMicro, "import_results") as mock_import,
        ):
            fortran_micro.run_full(hdf5_path)
        assert mock_import.call_args.args[0] == data_dir

    def test_tmpdir_cleaned_on_success(self, fortran_micro, tmp_path):
        """The outer tmpdir must be removed after a successful run."""
        hdf5_path = tmp_path / "run.h5"
        with h5py.File(str(hdf5_path), "w") as _:
            pass

        created_tmpdirs = []

        def fake_exec(work_dir):
            created_tmpdirs.append(work_dir)
            return work_dir / "data"

        with (
            patch.object(FortranMicro, "exec_in_workdir", side_effect=fake_exec),
            patch.object(FortranMicro, "import_results"),
        ):
            fortran_micro.run_full(hdf5_path, keep_tmpdir=False)

        assert created_tmpdirs, "exec_in_workdir was never called"
        assert not created_tmpdirs[0].exists()

    def test_tmpdir_preserved_with_keep_tmpdir(self, fortran_micro, tmp_path):
        """The outer tmpdir must be preserved when keep_tmpdir=True."""
        hdf5_path = tmp_path / "run.h5"
        with h5py.File(str(hdf5_path), "w") as _:
            pass

        created_tmpdirs = []

        def fake_exec(work_dir):
            created_tmpdirs.append(work_dir)
            work_dir.mkdir(parents=True, exist_ok=True)
            return work_dir / "data"

        with (
            patch.object(FortranMicro, "exec_in_workdir", side_effect=fake_exec),
            patch.object(FortranMicro, "import_results"),
        ):
            fortran_micro.run_full(hdf5_path, keep_tmpdir=True)

        assert created_tmpdirs
        # keep_tmpdir=True means the directory is preserved
        assert created_tmpdirs[0].exists()

    def test_tmpdir_not_in_system_tmp(self, fortran_micro, tmp_path):
        """run_full must create its tmpdir alongside the HDF5 file, not in /tmp."""
        hdf5_path = tmp_path / "run.h5"
        with h5py.File(str(hdf5_path), "w") as _:
            pass

        created_tmpdirs = []

        def fake_exec(work_dir):
            created_tmpdirs.append(work_dir)
            return work_dir / "data"

        with (
            patch.object(FortranMicro, "exec_in_workdir", side_effect=fake_exec),
            patch.object(FortranMicro, "import_results"),
        ):
            fortran_micro.run_full(hdf5_path, keep_tmpdir=False)

        assert created_tmpdirs
        # tmpdir parent should be same as hdf5_path.parent
        assert created_tmpdirs[0].parent == tmp_path


# ---------------------------------------------------------------------------
# Integration test with real fixture data
# ---------------------------------------------------------------------------

# These constants mirror the ones in tests/dataio/test_real_data.py.
_MICRO_FILE_CODE = "_PLG2_tPA01_TB-xiii"
_PARAM_OVERRIDES = {
    "fibrinogen_length": "45nm",
    "fibrinogen_radius": "1.2nm",
    "micro_log_lvl": 40,
    "micro_version": "micro_rates",
    "snap_proportion": 0.66666667,
}


def _build_reference_v200(fixture_path: str) -> dict:
    """Read v1.95.0 fixture data, apply overrides, convert to v2.0.0.

    Returns a dict with v2.0.0 dataset names as keys and numpy arrays as
    values, plus a ``"params"`` key.
    """
    from lysis.dataio.dataspec import dataspec
    from lysis.dataio.fileops import read_data_collection
    from lysis.dataio.dataconvert import convert_data

    raw = read_data_collection(
        fixture_path,
        collections=[dataspec["v1.95.0"]["microscale_out"]],
        file_codes=[_MICRO_FILE_CODE],
    )
    for key, value in _PARAM_OVERRIDES.items():
        for section in raw["params"].values():
            if isinstance(section, dict):
                section[key] = value

    return convert_data(raw, "v1.95.0", "v2.0.0")


class TestImportResultsWithFixture:
    """Integration tests: import real fixture data and compare to reference.

    Uses the truncated Fortran output in ``tests/fixtures/fortran_sample/``
    (the same dataset exercised in ``tests/dataio/test_real_data.py``).

    The test simulates what happens after the Fortran binary finishes:
    copies the fixture files into a staging directory, writes a v1.99.0
    ``params.json``, then calls :meth:`FortranMicro.import_results`.
    The resulting HDF5 datasets are compared element-by-element against a
    reference produced by the established v1.95.0 → v2.0.0 conversion
    pipeline.
    """

    @pytest.fixture
    def reference(self, fortran_sample_path):
        """Reference v2.0.0 data dict via the established conversion path."""
        return _build_reference_v200(fortran_sample_path)

    @pytest.fixture
    def imported_hdf5(self, fortran_sample_path, reference, tmp_path):
        """Run import_results on a copy of the fixture and return the HDF5 path."""
        # Copy fixture to staging directory (import_results may delete it)
        data_dir = tmp_path / "staging" / "data" / "fortran_sample"
        shutil.copytree(fortran_sample_path, data_dir)

        # Write params.json with micro_params so the v1.99.0 spec can parse
        # parameters from JSON rather than from the log file.
        micro_params = reference["params"]["micro_params"]
        with open(data_dir / "params.json", "w") as fh:
            json.dump({"micro_params": micro_params}, fh, indent=4, default=str)

        hdf5_path = tmp_path / "result.h5"
        with h5py.File(str(hdf5_path), "w") as _:
            pass  # empty target file

        FortranMicro.import_results(
            data_dir, hdf5_path,
            file_code=_MICRO_FILE_CODE,
            keep_tmpdir=True,
        )
        return hdf5_path

    # ── Dataset presence ────────────────────────────────────────────

    def test_all_microscale_datasets_present(self, imported_hdf5, reference):
        """Every v2.0.0 microscale_out dataset in the reference must be in the HDF5 file."""
        from lysis.dataio.dataspec import dataspec
        spec = dataspec["v2.0.0"]["microscale_out"]
        with h5py.File(str(imported_hdf5), "r") as f:
            for name, ds_spec in spec.data.items():
                if name not in reference:
                    continue
                loc = ds_spec.data_location
                assert loc in f, f"Missing dataset: {loc} (name={name})"

    # ── Per-dataset value comparisons ───────────────────────────────

    def test_pli_first_time_matches(self, imported_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["pli_first_time"].data_location
        with h5py.File(str(imported_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["pli_first_time"])

    def test_tpa_final_num_matches(self, imported_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["tpa_final_num"].data_location
        with h5py.File(str(imported_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["tpa_final_num"])

    def test_fiber_degraded_matches(self, imported_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["fiber_degraded"].data_location
        with h5py.File(str(imported_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["fiber_degraded"])

    def test_sim_final_time_matches(self, imported_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["sim_final_time"].data_location
        with h5py.File(str(imported_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["sim_final_time"])

    def test_pli_generated_num_matches(self, imported_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["pli_generated_num"].data_location
        with h5py.File(str(imported_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["pli_generated_num"])

    def test_tpa_leaving_time_matches(self, imported_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["tpa_leaving_time"].data_location
        with h5py.File(str(imported_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["tpa_leaving_time"])

    def test_tpa_unbound_by_pli_matches(self, imported_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["tpa_unbound_by_pli"].data_location
        with h5py.File(str(imported_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["tpa_unbound_by_pli"])

    def test_tpa_unbound_kinetic_matches(self, imported_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["tpa_unbound_kinetic"].data_location
        with h5py.File(str(imported_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["tpa_unbound_kinetic"])

    # ── Dtype checks ────────────────────────────────────────────────

    def test_dtypes_match_spec(self, imported_hdf5):
        """Each dataset's HDF5 dtype must match the v2.0.0 dataspec."""
        from lysis.dataio.dataspec import dataspec
        spec = dataspec["v2.0.0"]["microscale_out"]
        with h5py.File(str(imported_hdf5), "r") as f:
            for name, ds_spec in spec.data.items():
                loc = ds_spec.data_location
                if loc not in f:
                    continue
                expected_dtype = np.dtype(ds_spec.dtype)
                actual_dtype = f[loc].dtype
                assert actual_dtype == expected_dtype, (
                    f"{name}: expected {expected_dtype}, got {actual_dtype}"
                )

    # ── Shape checks ────────────────────────────────────────────────

    def test_shapes_match_reference(self, imported_hdf5, reference):
        """Each dataset's shape must match the reference conversion."""
        from lysis.dataio.dataspec import dataspec
        spec = dataspec["v2.0.0"]["microscale_out"]
        with h5py.File(str(imported_hdf5), "r") as f:
            for name, ds_spec in spec.data.items():
                if name not in reference:
                    continue
                loc = ds_spec.data_location
                if loc not in f:
                    continue
                expected_shape = np.asarray(reference[name]).shape
                actual_shape = f[loc].shape
                assert actual_shape == expected_shape, (
                    f"{name}: expected shape {expected_shape}, got {actual_shape}"
                )
