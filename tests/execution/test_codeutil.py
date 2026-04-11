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
# TestFortranMicroArrayExecCommand
# ---------------------------------------------------------------------------


class TestFortranMicroArrayExecCommand:
    """Tests for array-mode behaviour in :meth:`FortranMicro.exec_command`."""

    def test_array_mode_appends_index_suffix(self, tmp_run):
        """index=3 with n_array_jobs must append '__03' to out_file_code."""
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe",
                          index=3, n_array_jobs=10)
        fm.exec_command()
        assert fm.out_file_code == "__03"

    def test_array_mode_distributes_simulations_evenly(self, tmp_path):
        """With 1000 sims / 10 jobs each job should get exactly 100."""
        r = Run(str(tmp_path))
        r.initialize_micro_param({"micro_simulations": 1000})
        fm = FortranMicro(run=r, executable="/bin/micro.exe",
                          index=0, n_array_jobs=10)
        cmd = fm.exec_command()
        assert "--simulations" in cmd
        assert cmd[cmd.index("--simulations") + 1] == "100"

    def test_array_mode_last_job_no_extra_when_divisible(self, tmp_path):
        """When total sims is divisible by n_jobs, all jobs get equal count."""
        r = Run(str(tmp_path))
        r.initialize_micro_param({"micro_simulations": 100})
        fm = FortranMicro(run=r, executable="/bin/micro.exe",
                          index=9, n_array_jobs=10)
        cmd = fm.exec_command()
        assert cmd[cmd.index("--simulations") + 1] == "10"

    def test_array_mode_remainder_goes_to_first_jobs(self, tmp_path):
        """With 103 sims / 10 jobs, jobs 0-2 get 11 and jobs 3-9 get 10."""
        r = Run(str(tmp_path))
        r.initialize_micro_param({"micro_simulations": 103})

        fm0 = FortranMicro(run=r, executable="/bin/micro.exe",
                           index=0, n_array_jobs=10)
        fm2 = FortranMicro(run=r, executable="/bin/micro.exe",
                           index=2, n_array_jobs=10)
        fm3 = FortranMicro(run=r, executable="/bin/micro.exe",
                           index=3, n_array_jobs=10)

        cmd0 = fm0.exec_command()
        cmd2 = fm2.exec_command()
        cmd3 = fm3.exec_command()

        assert cmd0[cmd0.index("--simulations") + 1] == "11"
        assert cmd2[cmd2.index("--simulations") + 1] == "11"
        assert cmd3[cmd3.index("--simulations") + 1] == "10"

    def test_array_mode_seed_from_generate_state_n(self, tmp_run):
        """Array-mode seed for job i must come from generate_state(N)[i]."""
        seed = tmp_run.micro_params.micro_seed
        n = 10
        stream = np.random.SeedSequence(seed)
        expected_seed = int(np.int32(stream.generate_state(n)[3]))

        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe",
                          index=3, n_array_jobs=n)
        cmd = fm.exec_command()
        assert "--seed" in cmd
        assert cmd[cmd.index("--seed") + 1] == str(expected_seed)

    def test_array_mode_seed_uses_generate_state_n_formula(self, tmp_run):
        """Array-mode seed for job i must use generate_state(N)[i], not generate_state(i+1)[i]."""
        seed = tmp_run.micro_params.micro_seed
        n = 10
        idx = 3
        # The correct array formula
        expected = int(np.int32(
            np.random.SeedSequence(seed).generate_state(n)[idx]
        ))

        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe",
                          index=idx, n_array_jobs=n)
        cmd = fm.exec_command()
        assert cmd[cmd.index("--seed") + 1] == str(expected)

    def test_array_mode_seeds_are_independent_per_job(self, tmp_run):
        """Each job's seed must be unique across all n_jobs."""
        seed = tmp_run.micro_params.micro_seed
        n = 10
        seeds = set()
        for i in range(n):
            fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe",
                              index=i, n_array_jobs=n)
            cmd = fm.exec_command()
            s = cmd[cmd.index("--seed") + 1]
            seeds.add(s)
        assert len(seeds) == n, "Not all array job seeds are unique"

    def test_legacy_mode_still_forces_one_simulation(self, tmp_run):
        """Without n_array_jobs, index alone must still force micro_simulations=1."""
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe", index=0)
        cmd = fm.exec_command()
        assert "--simulations" in cmd
        assert cmd[cmd.index("--simulations") + 1] == "1"


# ---------------------------------------------------------------------------
# TestFortranMicroFromHdf5Array
# ---------------------------------------------------------------------------


class TestFortranMicroFromHdf5Array:
    """Tests for :meth:`FortranMicro.from_hdf5` with ``n_array_jobs``."""

    def test_passes_n_array_jobs(self, micro_hdf5):
        fm = FortranMicro.from_hdf5(micro_hdf5, "/bin/micro.exe",
                                     index=2, n_array_jobs=8)
        assert fm.n_array_jobs == 8

    def test_n_array_jobs_none_by_default(self, micro_hdf5):
        fm = FortranMicro.from_hdf5(micro_hdf5, "/bin/micro.exe")
        assert fm.n_array_jobs is None


# ---------------------------------------------------------------------------
# TestFortranMicroImportArrayResults
# ---------------------------------------------------------------------------


class TestFortranMicroImportArrayResults:
    """Tests for :meth:`FortranMicro.import_array_results`."""

    # Binary dtypes for v1.99.0 microscale_out files
    _DATASETS = {
        "firstPLi": np.float64,
        "lysis":    np.float64,
        "tPA_time": np.float64,
        "lasttPA":  np.int32,
        "lyscomplete": np.uint32,
        "PLi":      np.int32,
        "tPAPLiunbd": np.int32,
        "tPAunbind":  np.int32,
    }

    def _write_chunk(self, data_dir: Path, index: int, n_sims: int,
                     run_code: str, base_code: str = "") -> None:
        """Write fake binary chunk files for job ``index``."""
        chunk_code = f"{base_code}__{index:02}"
        for name, dtype in self._DATASETS.items():
            arr = np.arange(index * n_sims, (index + 1) * n_sims, dtype=dtype)
            fname = f"{name}{chunk_code}.dat"
            arr.tofile(data_dir / fname)
        log_file = data_dir / f"micro{chunk_code}.txt"
        log_file.write_text(f"chunk {index}\n")

    def _write_hdf5(self, path: Path, total_sims: int) -> None:
        """Write a minimal HDF5 file with micro_simulations attribute."""
        mp = MicroParameters(micro_simulations=total_sims)
        with h5py.File(str(path), "w") as f:
            f.attrs[CONST.DATASPEC_VERSION_ATTR] = "v2.0.0"
            grp = f.require_group("micro_data")
            for k, v in mp.to_basedict().items():
                grp.attrs[k] = str(v) if not isinstance(v, (int, float, bool)) else v

    @patch.object(FortranMicro, "import_results")
    def test_calls_import_results(self, mock_import, tmp_path):
        """import_array_results must delegate to import_results after merging."""
        n_jobs, n_sims = 3, 4
        run_code = "run-01"
        staging = tmp_path / "staging"
        data_dir = staging / "data" / run_code
        data_dir.mkdir(parents=True)
        hdf5 = tmp_path / f"{run_code}.h5"
        self._write_hdf5(hdf5, n_jobs * n_sims)

        # Write params.json
        import json
        (data_dir / "params.json").write_text(
            json.dumps({"micro_params": {"micro_simulations": n_sims}})
        )
        for i in range(n_jobs):
            self._write_chunk(data_dir, i, n_sims, run_code)

        FortranMicro.import_array_results(
            staging, hdf5, n_jobs=n_jobs, run_code=run_code,
            keep_tmpdir=True,
        )
        mock_import.assert_called_once()

    @patch.object(FortranMicro, "import_results")
    def test_merged_binary_has_correct_length(self, mock_import, tmp_path):
        """The merged binary for 'firstPLi' must have n_jobs*n_sims values."""
        n_jobs, n_sims = 3, 4
        run_code = "run-01"
        staging = tmp_path / "staging"
        data_dir = staging / "data" / run_code
        data_dir.mkdir(parents=True)
        hdf5 = tmp_path / f"{run_code}.h5"
        self._write_hdf5(hdf5, n_jobs * n_sims)

        import json
        (data_dir / "params.json").write_text(
            json.dumps({"micro_params": {"micro_simulations": n_sims}})
        )
        for i in range(n_jobs):
            self._write_chunk(data_dir, i, n_sims, run_code)

        FortranMicro.import_array_results(
            staging, hdf5, n_jobs=n_jobs, run_code=run_code,
            keep_tmpdir=True,
        )
        merged_file = staging / "merged" / run_code / f"firstPLi.dat"
        arr = np.fromfile(merged_file, dtype=np.float64)
        assert len(arr) == n_jobs * n_sims

    @patch.object(FortranMicro, "import_results")
    def test_merged_binary_is_in_order(self, mock_import, tmp_path):
        """Chunks must be concatenated in ascending index order."""
        n_jobs, n_sims = 3, 4
        run_code = "run-01"
        staging = tmp_path / "staging"
        data_dir = staging / "data" / run_code
        data_dir.mkdir(parents=True)
        hdf5 = tmp_path / f"{run_code}.h5"
        self._write_hdf5(hdf5, n_jobs * n_sims)

        import json
        (data_dir / "params.json").write_text(
            json.dumps({"micro_params": {"micro_simulations": n_sims}})
        )
        for i in range(n_jobs):
            self._write_chunk(data_dir, i, n_sims, run_code)

        FortranMicro.import_array_results(
            staging, hdf5, n_jobs=n_jobs, run_code=run_code,
            keep_tmpdir=True,
        )
        merged_file = staging / "merged" / run_code / "firstPLi.dat"
        arr = np.fromfile(merged_file, dtype=np.float64)
        expected = np.arange(n_jobs * n_sims, dtype=np.float64)
        np.testing.assert_array_equal(arr, expected)

    @patch.object(FortranMicro, "import_results")
    def test_merged_params_json_has_total_simulations(self, mock_import, tmp_path):
        """merged/params.json must reflect the total simulation count from HDF5."""
        n_jobs, n_sims = 3, 4
        total_sims = n_jobs * n_sims
        run_code = "run-01"
        staging = tmp_path / "staging"
        data_dir = staging / "data" / run_code
        data_dir.mkdir(parents=True)
        hdf5 = tmp_path / f"{run_code}.h5"
        self._write_hdf5(hdf5, total_sims)

        import json
        (data_dir / "params.json").write_text(
            json.dumps({"micro_params": {"micro_simulations": n_sims}})
        )
        for i in range(n_jobs):
            self._write_chunk(data_dir, i, n_sims, run_code)

        FortranMicro.import_array_results(
            staging, hdf5, n_jobs=n_jobs, run_code=run_code,
            keep_tmpdir=True,
        )
        merged_params = json.loads(
            (staging / "merged" / run_code / "params.json").read_text()
        )
        assert merged_params["micro_params"]["micro_simulations"] == total_sims

    @patch.object(FortranMicro, "import_results")
    def test_merged_dir_removed_on_success(self, mock_import, tmp_path):
        """Merged run_code subdirectory must be removed on success when keep_tmpdir=False."""
        n_jobs, n_sims = 2, 5
        run_code = "run-01"
        staging = tmp_path / "staging"
        data_dir = staging / "data" / run_code
        data_dir.mkdir(parents=True)
        hdf5 = tmp_path / f"{run_code}.h5"
        self._write_hdf5(hdf5, n_jobs * n_sims)

        import json
        (data_dir / "params.json").write_text(
            json.dumps({"micro_params": {"micro_simulations": n_sims}})
        )
        for i in range(n_jobs):
            self._write_chunk(data_dir, i, n_sims, run_code)

        FortranMicro.import_array_results(
            staging, hdf5, n_jobs=n_jobs, run_code=run_code,
            keep_tmpdir=False,
        )
        assert not (staging / "merged" / run_code).exists()

    @patch.object(FortranMicro, "import_results")
    def test_merged_dir_preserved_with_keep_tmpdir(self, mock_import, tmp_path):
        """Merged run_code subdirectory must survive when keep_tmpdir=True."""
        n_jobs, n_sims = 2, 5
        run_code = "run-01"
        staging = tmp_path / "staging"
        data_dir = staging / "data" / run_code
        data_dir.mkdir(parents=True)
        hdf5 = tmp_path / f"{run_code}.h5"
        self._write_hdf5(hdf5, n_jobs * n_sims)

        import json
        (data_dir / "params.json").write_text(
            json.dumps({"micro_params": {"micro_simulations": n_sims}})
        )
        for i in range(n_jobs):
            self._write_chunk(data_dir, i, n_sims, run_code)

        FortranMicro.import_array_results(
            staging, hdf5, n_jobs=n_jobs, run_code=run_code,
            keep_tmpdir=True,
        )
        assert (staging / "merged" / run_code).exists()


# ---------------------------------------------------------------------------
# TestFortranMicroRunArrayFull
# ---------------------------------------------------------------------------


class TestFortranMicroRunArrayFull:
    """Tests for :meth:`FortranMicro.run_array_full`."""

    def test_calls_exec_in_workdir_n_times(self, fortran_micro, tmp_path):
        """run_array_full must call exec_in_workdir once per job."""
        hdf5_path = tmp_path / "run.h5"
        with h5py.File(str(hdf5_path), "w") as _:
            pass
        n_jobs = 4
        exec_calls = []

        def fake_exec(work_dir):
            exec_calls.append(work_dir)
            return work_dir / "data"

        with (
            patch.object(FortranMicro, "exec_in_workdir", side_effect=fake_exec),
            patch.object(FortranMicro, "import_array_results"),
        ):
            fortran_micro.run_array_full(hdf5_path, n_jobs=n_jobs)

        assert len(exec_calls) == n_jobs

    def test_each_child_has_correct_index(self, tmp_run, tmp_path):
        """Each child FortranMicro must have index equal to its position."""
        hdf5_path = tmp_path / "run.h5"
        with h5py.File(str(hdf5_path), "w") as _:
            pass
        n_jobs = 3
        created_indices = []

        original_init = FortranMicro.__init__

        def patched_exec(self_fm, work_dir):
            created_indices.append(self_fm.index)
            return work_dir / "data"

        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe")
        with (
            patch.object(FortranMicro, "exec_in_workdir", patched_exec),
            patch.object(FortranMicro, "import_array_results"),
        ):
            fm.run_array_full(hdf5_path, n_jobs=n_jobs)

        assert created_indices == list(range(n_jobs))

    def test_calls_import_array_results_once(self, fortran_micro, tmp_path):
        """run_array_full must call import_array_results exactly once."""
        hdf5_path = tmp_path / "run.h5"
        with h5py.File(str(hdf5_path), "w") as _:
            pass

        with (
            patch.object(FortranMicro, "exec_in_workdir",
                         return_value=tmp_path / "data"),
            patch.object(FortranMicro, "import_array_results") as mock_import,
        ):
            fortran_micro.run_array_full(hdf5_path, n_jobs=3)

        mock_import.assert_called_once()

    def test_tmpdir_not_in_system_tmp(self, fortran_micro, tmp_path):
        """run_array_full must create tmpdir alongside HDF5, not in /tmp."""
        hdf5_path = tmp_path / "run.h5"
        with h5py.File(str(hdf5_path), "w") as _:
            pass
        created_tmpdirs = []

        def fake_exec(work_dir):
            created_tmpdirs.append(work_dir)
            return work_dir / "data"

        with (
            patch.object(FortranMicro, "exec_in_workdir", side_effect=fake_exec),
            patch.object(FortranMicro, "import_array_results"),
        ):
            fortran_micro.run_array_full(hdf5_path, n_jobs=2)

        assert created_tmpdirs
        assert created_tmpdirs[0].parent == tmp_path

    def test_tmpdir_removed_on_success(self, fortran_micro, tmp_path):
        """run_array_full must remove tmpdir on success when keep_tmpdir=False."""
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
            patch.object(FortranMicro, "import_array_results"),
        ):
            fortran_micro.run_array_full(hdf5_path, n_jobs=2, keep_tmpdir=False)

        assert created_tmpdirs
        assert not created_tmpdirs[0].exists()

    def test_tmpdir_preserved_with_keep_tmpdir(self, fortran_micro, tmp_path):
        """run_array_full must preserve tmpdir when keep_tmpdir=True."""
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
            patch.object(FortranMicro, "import_array_results"),
        ):
            fortran_micro.run_array_full(hdf5_path, n_jobs=2, keep_tmpdir=True)

        assert created_tmpdirs
        assert created_tmpdirs[0].exists()

    def test_n_array_jobs_passed_to_children(self, tmp_run, tmp_path):
        """Each child FortranMicro must receive n_array_jobs equal to n_jobs."""
        hdf5_path = tmp_path / "run.h5"
        with h5py.File(str(hdf5_path), "w") as _:
            pass
        n_jobs = 5
        child_n_array_jobs = []

        def patched_exec(self_fm, work_dir):
            child_n_array_jobs.append(self_fm.n_array_jobs)
            return work_dir / "data"

        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe")
        with (
            patch.object(FortranMicro, "exec_in_workdir", patched_exec),
            patch.object(FortranMicro, "import_array_results"),
        ):
            fm.run_array_full(hdf5_path, n_jobs=n_jobs)

        assert all(n == n_jobs for n in child_n_array_jobs)


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


# ---------------------------------------------------------------------------
# Fortran binary execution integration test
# ---------------------------------------------------------------------------

# Number of simulations to run (subset of the 50,000 in the fixture).
_BINARY_TEST_SIMULATIONS = 1000

# Repo root, used to locate the compiled Fortran binary.
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent

# Fixture seed (from micro_PLG2_tPA01_TB-xiii.txt log file).
_FIXTURE_SEED = 2133256963


@pytest.mark.fortran_binary
class TestFortranBinaryExecution:
    """Integration tests that execute the real compiled Fortran binary.

    These tests compile-and-run the microscale simulation with the same
    parameters as the ``tests/fixtures/fortran_sample/`` fixture but only
    1,000 simulations (instead of 50,000) and compare the output to the
    first 1,000 elements of the fixture reference data.

    Requires:

    * The compiled Fortran binary at ``bin/micro_rates`` (repo root).
    * The ``tests/fixtures/fortran_sample/`` fixture.
    * The ``lysis`` conda environment on ``$PATH``.

    Skip with ``-m "not fortran_binary"`` to exclude these tests.
    """

    @pytest.fixture(autouse=True)
    def _skip_if_no_binary(self):
        binary = _REPO_ROOT / "bin" / "micro_rates"
        if not binary.exists():
            pytest.skip("Compiled Fortran binary not found at bin/micro_rates")

    # Datasets to exclude from comparisons: ``micro_log`` is a text log
    # whose length depends on total simulations (not sliceable), and
    # ``params`` is metadata, not numeric output.
    _SKIP_DATASETS = {"params", "micro_log"}

    @pytest.fixture
    def reference(self, fortran_sample_path):
        """Reference v2.0.0 data sliced to the first 1,000 simulations."""
        ref = _build_reference_v200(fortran_sample_path)
        sliced = {}
        for key, value in ref.items():
            if key in self._SKIP_DATASETS:
                continue
            arr = np.asarray(value)
            if arr.ndim >= 1 and arr.shape[0] >= _BINARY_TEST_SIMULATIONS:
                sliced[key] = arr[:_BINARY_TEST_SIMULATIONS]
            else:
                sliced[key] = arr
        return sliced

    @pytest.fixture
    def executed_hdf5(self, tmp_path):
        """Run the Fortran binary and import results into a fresh HDF5 file.

        Returns the path to the HDF5 file containing the imported results.
        """
        from lysis.config.constants import Q_

        binary = str(_REPO_ROOT / "bin" / "micro_rates")
        hdf5_path = tmp_path / "fortran-run.h5"

        # Create a minimal HDF5 with the fixture's micro parameters
        mp = MicroParameters(
            nodes_in_micro_row=13,
            fiber_radius=Q_("61.5 nanometer"),
            micro_seed=_FIXTURE_SEED,
            micro_simulations=_BINARY_TEST_SIMULATIONS,
        )
        _write_micro_hdf5(hdf5_path, mp)

        fm = FortranMicro.from_hdf5(hdf5_path, binary)
        data_dir = fm.exec_in_workdir(tmp_path)
        FortranMicro.import_results(
            data_dir, hdf5_path, file_code=fm.out_file_code, keep_tmpdir=True
        )
        return hdf5_path

    # ── Dataset presence ────────────────────────────────────────────

    def test_all_microscale_datasets_present(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        spec = dataspec["v2.0.0"]["microscale_out"]
        with h5py.File(str(executed_hdf5), "r") as f:
            for name in reference:
                ds_spec = spec.data.get(name)
                if ds_spec is None:
                    continue
                assert ds_spec.data_location in f, (
                    f"Missing dataset: {ds_spec.data_location} (name={name})"
                )

    # ── Per-dataset value comparisons ───────────────────────────────

    def test_pli_first_time_matches(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["pli_first_time"].data_location
        with h5py.File(str(executed_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["pli_first_time"])

    def test_tpa_final_num_matches(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["tpa_final_num"].data_location
        with h5py.File(str(executed_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["tpa_final_num"])

    def test_fiber_degraded_matches(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["fiber_degraded"].data_location
        with h5py.File(str(executed_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["fiber_degraded"])

    def test_sim_final_time_matches(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["sim_final_time"].data_location
        with h5py.File(str(executed_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["sim_final_time"])

    def test_pli_generated_num_matches(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["pli_generated_num"].data_location
        with h5py.File(str(executed_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["pli_generated_num"])

    def test_tpa_leaving_time_matches(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["tpa_leaving_time"].data_location
        with h5py.File(str(executed_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["tpa_leaving_time"])

    def test_tpa_unbound_by_pli_matches(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["tpa_unbound_by_pli"].data_location
        with h5py.File(str(executed_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["tpa_unbound_by_pli"])

    def test_tpa_unbound_kinetic_matches(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["tpa_unbound_kinetic"].data_location
        with h5py.File(str(executed_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["tpa_unbound_kinetic"])

    # ── Dtype checks ────────────────────────────────────────────────

    def test_dtypes_match_spec(self, executed_hdf5):
        from lysis.dataio.dataspec import dataspec
        spec = dataspec["v2.0.0"]["microscale_out"]
        with h5py.File(str(executed_hdf5), "r") as f:
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

    def test_shapes_match_expected(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        spec = dataspec["v2.0.0"]["microscale_out"]
        with h5py.File(str(executed_hdf5), "r") as f:
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


# ---------------------------------------------------------------------------
# Fortran array execution integration test
# ---------------------------------------------------------------------------

# Number of array jobs and simulations per job for array integration tests.
_ARRAY_N_JOBS = 10
_ARRAY_SIMS_PER_JOB = 100
_ARRAY_TOTAL_SIMS = _ARRAY_N_JOBS * _ARRAY_SIMS_PER_JOB


@pytest.mark.fortran_binary
class TestFortranArrayExecution:
    """Integration tests for the ``run_array_full`` workflow.

    Executes the real compiled Fortran binary as :data:`_ARRAY_N_JOBS`
    independent chunk jobs (each running :data:`_ARRAY_SIMS_PER_JOB`
    simulations), concatenates the results, and verifies:

    1. All datasets are present and have the correct total length.
    2. Results are reproducible — running again with the same seed yields
       identical output.
    3. The aggregated array results are statistically consistent with a
       single-process run of the same total number of simulations.

    .. note::

        Array jobs use ``numpy.random.SeedSequence`` to derive independent
        per-chunk seeds from the base ``micro_seed``.  The Fortran KISS32 RNG
        is initialised with a *different* seed for each chunk, so array output
        cannot be **bit-for-bit identical** to a single-process run that uses
        the base seed directly.  Tests therefore compare statistical properties
        (means, standard deviations) rather than exact values.

    Requires:

    * The compiled Fortran binary at ``bin/micro_rates`` (repo root).
    * The ``lysis`` conda environment on ``$PATH``.

    Skip with ``-m "not fortran_binary"`` to exclude these tests.
    """

    @pytest.fixture(autouse=True)
    def _skip_if_no_binary(self):
        binary = _REPO_ROOT / "bin" / "micro_rates"
        if not binary.exists():
            pytest.skip("Compiled Fortran binary not found at bin/micro_rates")

    _SKIP_DATASETS = {"params", "micro_log"}

    def _make_hdf5(self, tmp_path: Path, seed: int = _FIXTURE_SEED) -> Path:
        """Write a minimal HDF5 file using the fixture parameters."""
        from lysis.config.constants import Q_

        hdf5_path = tmp_path / "array-run.h5"
        mp = MicroParameters(
            nodes_in_micro_row=13,
            fiber_radius=Q_("61.5 nanometer"),
            micro_seed=seed,
            micro_simulations=_ARRAY_TOTAL_SIMS,
        )
        _write_micro_hdf5(hdf5_path, mp)
        return hdf5_path

    @pytest.fixture
    def array_hdf5(self, tmp_path):
        """Execute ``run_array_full`` and return the HDF5 path.

        Runs :data:`_ARRAY_N_JOBS` chunks of :data:`_ARRAY_SIMS_PER_JOB`
        simulations each using the fixture seed.
        """
        binary = str(_REPO_ROOT / "bin" / "micro_rates")
        hdf5_path = self._make_hdf5(tmp_path)
        fm = FortranMicro.from_hdf5(hdf5_path, binary)
        fm.run_array_full(hdf5_path, n_jobs=_ARRAY_N_JOBS)
        return hdf5_path

    @pytest.fixture
    def array_hdf5_second_run(self, tmp_path):
        """Independent second execution with the same seed for reproducibility.

        Uses a separate ``tmp_path`` sub-directory so there is no chance of
        file-name collisions with :attr:`array_hdf5`.
        """
        binary = str(_REPO_ROOT / "bin" / "micro_rates")
        sub = tmp_path / "second"
        sub.mkdir()
        hdf5_path = self._make_hdf5(sub)
        fm = FortranMicro.from_hdf5(hdf5_path, binary)
        fm.run_array_full(hdf5_path, n_jobs=_ARRAY_N_JOBS)
        return hdf5_path

    @pytest.fixture
    def single_hdf5(self, tmp_path):
        """Single-process run of :data:`_ARRAY_TOTAL_SIMS` simulations.

        Used as a statistical reference.  Results will not match the array run
        bit-for-bit because the RNG seeds differ, but aggregate statistics
        (mean, std) should be equivalent.
        """
        from lysis.config.constants import Q_

        binary = str(_REPO_ROOT / "bin" / "micro_rates")
        sub = tmp_path / "single"
        sub.mkdir()
        hdf5_path = sub / "single-run.h5"
        mp = MicroParameters(
            nodes_in_micro_row=13,
            fiber_radius=Q_("61.5 nanometer"),
            micro_seed=_FIXTURE_SEED,
            micro_simulations=_ARRAY_TOTAL_SIMS,
        )
        _write_micro_hdf5(hdf5_path, mp)
        fm = FortranMicro.from_hdf5(hdf5_path, binary)
        data_dir = fm.exec_in_workdir(sub)
        FortranMicro.import_results(
            data_dir, hdf5_path, file_code=fm.out_file_code, keep_tmpdir=False
        )
        return hdf5_path

    # ── Dataset presence / length ─────────────────────────────────────

    def test_all_datasets_present(self, array_hdf5):
        """All microscale_out datasets (except skipped ones) must be present."""
        from lysis.dataio.dataspec import dataspec
        spec = dataspec["v2.0.0"]["microscale_out"]
        with h5py.File(str(array_hdf5), "r") as f:
            for name, ds_spec in spec.data.items():
                if name in self._SKIP_DATASETS:
                    continue
                assert ds_spec.data_location in f, (
                    f"Missing dataset: {ds_spec.data_location} (name={name})"
                )

    def test_total_simulation_count(self, array_hdf5):
        """Leading dimension of every 1-D output must equal the total sims."""
        from lysis.dataio.dataspec import dataspec
        spec = dataspec["v2.0.0"]["microscale_out"]
        with h5py.File(str(array_hdf5), "r") as f:
            for name, ds_spec in spec.data.items():
                if name in self._SKIP_DATASETS:
                    continue
                loc = ds_spec.data_location
                if loc not in f:
                    continue
                arr = f[loc][...]
                if arr.ndim >= 1 and arr.shape[0] > 1:
                    assert arr.shape[0] == _ARRAY_TOTAL_SIMS, (
                        f"{name}: expected {_ARRAY_TOTAL_SIMS} rows, "
                        f"got {arr.shape[0]}"
                    )

    # ── Reproducibility ──────────────────────────────────────────────

    def test_reproducible_pli_first_time(self, array_hdf5, array_hdf5_second_run):
        """Two array runs with the same seed must produce identical pli_first_time."""
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["pli_first_time"].data_location
        with h5py.File(str(array_hdf5), "r") as f1, \
             h5py.File(str(array_hdf5_second_run), "r") as f2:
            np.testing.assert_array_equal(f1[loc][...], f2[loc][...])

    def test_reproducible_sim_final_time(self, array_hdf5, array_hdf5_second_run):
        """Two array runs with the same seed must produce identical sim_final_time."""
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["sim_final_time"].data_location
        with h5py.File(str(array_hdf5), "r") as f1, \
             h5py.File(str(array_hdf5_second_run), "r") as f2:
            np.testing.assert_array_equal(f1[loc][...], f2[loc][...])

    def test_reproducible_fiber_degraded(self, array_hdf5, array_hdf5_second_run):
        """Two array runs with the same seed must produce identical fiber_degraded."""
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["fiber_degraded"].data_location
        with h5py.File(str(array_hdf5), "r") as f1, \
             h5py.File(str(array_hdf5_second_run), "r") as f2:
            np.testing.assert_array_equal(f1[loc][...], f2[loc][...])

    # ── Statistical consistency with single-process run ──────────────

    def test_pli_first_time_mean_consistent(self, array_hdf5, single_hdf5):
        """Array mean of pli_first_time must be within 5 % of single-run mean.

        The array and single-process runs use different RNG seeds (by design),
        so exact equality is not expected.  With 1,000 simulations the sampling
        error is small enough that a 5 % tolerance is a meaningful consistency
        check without being fragile.
        """
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["pli_first_time"].data_location
        with h5py.File(str(array_hdf5), "r") as fa, \
             h5py.File(str(single_hdf5), "r") as fs:
            arr_mean = float(np.mean(fa[loc][...]))
            single_mean = float(np.mean(fs[loc][...]))
        rel_diff = abs(arr_mean - single_mean) / (abs(single_mean) + 1e-12)
        assert rel_diff < 0.05, (
            f"pli_first_time mean differs by {rel_diff:.1%}: "
            f"array={arr_mean:.4g}, single={single_mean:.4g}"
        )

    def test_sim_final_time_mean_consistent(self, array_hdf5, single_hdf5):
        """Array mean of sim_final_time must be within 5 % of single-run mean."""
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["sim_final_time"].data_location
        with h5py.File(str(array_hdf5), "r") as fa, \
             h5py.File(str(single_hdf5), "r") as fs:
            arr_mean = float(np.mean(fa[loc][...]))
            single_mean = float(np.mean(fs[loc][...]))
        rel_diff = abs(arr_mean - single_mean) / (abs(single_mean) + 1e-12)
        assert rel_diff < 0.05, (
            f"sim_final_time mean differs by {rel_diff:.1%}: "
            f"array={arr_mean:.4g}, single={single_mean:.4g}"
        )

    def test_fiber_degraded_fraction_consistent(self, array_hdf5, single_hdf5):
        """Fraction of degraded fibers must be within 5 % of single-run value."""
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["fiber_degraded"].data_location
        with h5py.File(str(array_hdf5), "r") as fa, \
             h5py.File(str(single_hdf5), "r") as fs:
            arr_frac = float(np.mean(fa[loc][...]))
            single_frac = float(np.mean(fs[loc][...]))
        rel_diff = abs(arr_frac - single_frac) / (abs(single_frac) + 1e-12)
        assert rel_diff < 0.05, (
            f"fiber_degraded fraction differs by {rel_diff:.1%}: "
            f"array={arr_frac:.4g}, single={single_frac:.4g}"
        )

    # ── Dtype / shape checks ─────────────────────────────────────────

    def test_dtypes_match_spec(self, array_hdf5):
        """All datasets must have the dtype specified by the dataspec."""
        from lysis.dataio.dataspec import dataspec
        spec = dataspec["v2.0.0"]["microscale_out"]
        with h5py.File(str(array_hdf5), "r") as f:
            for name, ds_spec in spec.data.items():
                loc = ds_spec.data_location
                if loc not in f:
                    continue
                expected_dtype = np.dtype(ds_spec.dtype)
                actual_dtype = f[loc].dtype
                assert actual_dtype == expected_dtype, (
                    f"{name}: expected {expected_dtype}, got {actual_dtype}"
                )
