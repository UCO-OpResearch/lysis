"""Tests for ``lysis run-micro`` CLI command."""

import json
import os
from pathlib import Path
from unittest.mock import MagicMock, patch

import h5py
import pytest
from click.testing import CliRunner

from lysis.cli import cli
from lysis.config.constants import CONST
from lysis.config.parameters import MicroParameters


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------


def _write_micro_hdf5(path: Path) -> None:
    """Write a minimal v2.0.0 HDF5 file with default MicroParameters."""
    mp = MicroParameters()
    with h5py.File(str(path), "w") as f:
        f.attrs[CONST.DATASPEC_VERSION_ATTR] = "v2.0.0"
        grp = f.require_group("micro_data")
        for k, v in mp.to_basedict().items():
            grp.attrs[k] = str(v) if not isinstance(v, (int, float, bool)) else v


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def micro_hdf5(tmp_path):
    """Minimal HDF5 file; returns Path."""
    h5_path = tmp_path / "run-01.h5"
    _write_micro_hdf5(h5_path)
    return h5_path


# ---------------------------------------------------------------------------
# Help
# ---------------------------------------------------------------------------


class TestRunMicroHelp:
    def test_help_exits_0(self, runner):
        result = runner.invoke(cli, ["run-micro", "--help"])
        assert result.exit_code == 0

    def test_help_mentions_executable(self, runner):
        result = runner.invoke(cli, ["run-micro", "--help"])
        assert "--executable" in result.output

    def test_help_mentions_slurm(self, runner):
        result = runner.invoke(cli, ["run-micro", "--help"])
        assert "--slurm" in result.output

    def test_help_mentions_keep_tmpdir(self, runner):
        result = runner.invoke(cli, ["run-micro", "--help"])
        assert "--keep-tmpdir" in result.output


# ---------------------------------------------------------------------------
# Missing required argument
# ---------------------------------------------------------------------------


class TestRunMicroMissingArgs:
    def test_missing_executable_exits_nonzero(self, runner, micro_hdf5):
        result = runner.invoke(cli, ["run-micro", str(micro_hdf5)])
        assert result.exit_code != 0

    def test_missing_hdf5_path_exits_nonzero(self, runner, tmp_path):
        result = runner.invoke(
            cli, ["run-micro", "--executable", "/bin/micro.exe"]
        )
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# Non-Slurm (local) execution
# ---------------------------------------------------------------------------


class TestRunMicroLocal:
    @patch("lysis.execution.fortran_micro.FortranMicro")
    def test_calls_from_hdf5(self, mock_cls, runner, micro_hdf5):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            ["run-micro", str(micro_hdf5), "--executable", "/bin/micro.exe"],
        )
        assert result.exit_code == 0, result.output
        mock_cls.from_hdf5.assert_called_once()
        call_args = mock_cls.from_hdf5.call_args
        assert str(micro_hdf5) in str(call_args.args[0])
        assert call_args.args[1] == "/bin/micro.exe"

    @patch("lysis.execution.fortran_micro.FortranMicro")
    def test_calls_run_full(self, mock_cls, runner, micro_hdf5):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            ["run-micro", str(micro_hdf5), "--executable", "/bin/micro.exe"],
        )
        assert result.exit_code == 0, result.output
        mock_fm.run_full.assert_called_once()

    @patch("lysis.execution.fortran_micro.FortranMicro")
    def test_keep_tmpdir_flag_passed_to_run_full(self, mock_cls, runner, micro_hdf5):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            [
                "run-micro", str(micro_hdf5),
                "--executable", "/bin/micro.exe",
                "--keep-tmpdir",
            ],
        )
        assert result.exit_code == 0, result.output
        call_kwargs = mock_fm.run_full.call_args.kwargs
        assert call_kwargs.get("keep_tmpdir") is True

    @patch("lysis.execution.fortran_micro.FortranMicro")
    def test_file_code_passed_to_from_hdf5(self, mock_cls, runner, micro_hdf5):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            [
                "run-micro", str(micro_hdf5),
                "--executable", "/bin/micro.exe",
                "--file-code", "_code",
            ],
        )
        assert result.exit_code == 0, result.output
        call_kwargs = mock_cls.from_hdf5.call_args.kwargs
        assert call_kwargs.get("out_file_code") == "_code"

    @patch("lysis.execution.fortran_micro.FortranMicro")
    def test_output_mentions_hdf5_path(self, mock_cls, runner, micro_hdf5):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            ["run-micro", str(micro_hdf5), "--executable", "/bin/micro.exe"],
        )
        assert result.exit_code == 0, result.output
        # Either the filename or path should appear in the output
        assert "run-01" in result.output

    @patch("lysis.execution.fortran_micro.FortranMicro")
    def test_exits_nonzero_when_from_hdf5_raises_valueerror(
        self, mock_cls, runner, micro_hdf5
    ):
        """A ValueError from from_hdf5 (e.g. wrong HDF5 state) must exit non-zero."""
        mock_cls.from_hdf5.side_effect = ValueError("microscale simulation has already been run")
        result = runner.invoke(
            cli,
            ["run-micro", str(micro_hdf5), "--executable", "/bin/micro.exe"],
        )
        assert result.exit_code != 0

    @patch("lysis.execution.fortran_micro.FortranMicro")
    def test_error_message_shown_when_from_hdf5_raises(
        self, mock_cls, runner, micro_hdf5
    ):
        """Error message from ValueError must appear in CLI output."""
        mock_cls.from_hdf5.side_effect = ValueError("microscale simulation has already been run")
        result = runner.invoke(
            cli,
            ["run-micro", str(micro_hdf5), "--executable", "/bin/micro.exe"],
        )
        assert "microscale simulation has already been run" in result.output


# ---------------------------------------------------------------------------
# Slurm dispatch
# ---------------------------------------------------------------------------


class TestRunMicroSlurm:
    @patch("lysis.tools.slurm.submit_micro_slurm_job", return_value=12345)
    def test_slurm_calls_submit(self, mock_submit, runner, micro_hdf5):
        result = runner.invoke(
            cli,
            [
                "run-micro", str(micro_hdf5),
                "--executable", "/bin/micro.exe",
                "--slurm",
            ],
        )
        assert result.exit_code == 0, result.output
        mock_submit.assert_called_once()

    @patch("lysis.tools.slurm.submit_micro_slurm_job", return_value=12345)
    def test_slurm_prints_job_id(self, mock_submit, runner, micro_hdf5):
        result = runner.invoke(
            cli,
            [
                "run-micro", str(micro_hdf5),
                "--executable", "/bin/micro.exe",
                "--slurm",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "12345" in result.output

    @patch("lysis.tools.slurm.submit_micro_slurm_job", return_value=1)
    def test_partition_forwarded_to_submit(self, mock_submit, runner, micro_hdf5):
        result = runner.invoke(
            cli,
            [
                "run-micro", str(micro_hdf5),
                "--executable", "/bin/micro.exe",
                "--slurm",
                "--partition", "long",
            ],
        )
        assert result.exit_code == 0, result.output
        call_kwargs = mock_submit.call_args.kwargs
        assert call_kwargs.get("partition") == "long"

    @patch("lysis.tools.slurm.submit_micro_slurm_job", return_value=1)
    def test_staging_root_forwarded(self, mock_submit, runner, micro_hdf5, tmp_path):
        staging = tmp_path / "staging"
        staging.mkdir()
        result = runner.invoke(
            cli,
            [
                "run-micro", str(micro_hdf5),
                "--executable", "/bin/micro.exe",
                "--slurm",
                "--staging-root", str(staging),
            ],
        )
        assert result.exit_code == 0, result.output
        call_kwargs = mock_submit.call_args.kwargs
        assert call_kwargs.get("staging_root") == str(staging)

    @patch("lysis.tools.slurm.submit_micro_slurm_job", return_value=1)
    def test_fast_tmp_root_forwarded(self, mock_submit, runner, micro_hdf5, tmp_path):
        result = runner.invoke(
            cli,
            [
                "run-micro", str(micro_hdf5),
                "--executable", "/bin/micro.exe",
                "--slurm",
                "--fast-tmp-root", "/nvme/scratch",
            ],
        )
        assert result.exit_code == 0, result.output
        call_kwargs = mock_submit.call_args.kwargs
        assert call_kwargs.get("fast_tmp_root") == "/nvme/scratch"

    @patch("lysis.tools.slurm.submit_micro_slurm_job", return_value=1)
    def test_keep_tmpdir_forwarded_to_submit(self, mock_submit, runner, micro_hdf5):
        result = runner.invoke(
            cli,
            [
                "run-micro", str(micro_hdf5),
                "--executable", "/bin/micro.exe",
                "--slurm",
                "--keep-tmpdir",
            ],
        )
        assert result.exit_code == 0, result.output
        call_kwargs = mock_submit.call_args.kwargs
        assert call_kwargs.get("keep_tmpdir") is True

    @patch("lysis.execution.fortran_micro.FortranMicro")
    def test_partition_without_slurm_does_not_crash(self, mock_cls, runner, micro_hdf5):
        """--partition without --slurm should not cause an error (it's ignored)."""
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            [
                "run-micro", str(micro_hdf5),
                "--executable", "/bin/micro.exe",
                "--partition", "long",
            ],
        )
        # Should succeed (partition is silently ignored in non-slurm mode)
        assert result.exit_code == 0, result.output
        mock_fm.run_full.assert_called_once()


# ---------------------------------------------------------------------------
# Batch (experiment / directory) helpers
# ---------------------------------------------------------------------------


def _make_experiment_dir(parent: Path, run_codes: list) -> Path:
    """Create a minimal experiment folder with HDF5 files and experiment.json."""
    exp_dir = parent / "my-experiment"
    exp_dir.mkdir()
    for rc in run_codes:
        _write_micro_hdf5(exp_dir / f"{rc}.h5")
    experiment_json = {
        "name": "my-experiment",
        "description": "",
        "created": "2026-01-01T00:00:00",
        "lysis_version": "test",
        "runs": [
            {"run_code": rc, "row_index": i, "description": "", "macro_params": None}
            for i, rc in enumerate(run_codes)
        ],
    }
    (exp_dir / "experiment.json").write_text(json.dumps(experiment_json))
    return exp_dir


@pytest.fixture
def experiment_dir(tmp_path):
    """Experiment folder with experiment.json and two HDF5 files."""
    return _make_experiment_dir(tmp_path, ["run-01", "run-02"])


@pytest.fixture
def hdf5_dir(tmp_path):
    """Directory with two HDF5 files but no experiment.json."""
    d = tmp_path / "raw-runs"
    d.mkdir()
    _write_micro_hdf5(d / "run-aa.h5")
    _write_micro_hdf5(d / "run-bb.h5")
    return d


# ---------------------------------------------------------------------------
# Batch local execution (experiment folder)
# ---------------------------------------------------------------------------


class TestRunMicroExperimentLocal:
    @patch("lysis.execution.fortran_micro.FortranMicro")
    def test_batch_with_experiment_json_calls_from_hdf5_for_each_run(
        self, mock_cls, runner, experiment_dir
    ):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            ["run-micro", str(experiment_dir), "--executable", "/bin/micro.exe"],
        )
        assert result.exit_code == 0, result.output
        assert mock_cls.from_hdf5.call_count == 2

    @patch("lysis.execution.fortran_micro.FortranMicro")
    def test_batch_with_experiment_json_calls_run_full_for_each_run(
        self, mock_cls, runner, experiment_dir
    ):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            ["run-micro", str(experiment_dir), "--executable", "/bin/micro.exe"],
        )
        assert result.exit_code == 0, result.output
        assert mock_fm.run_full.call_count == 2

    @patch("lysis.execution.fortran_micro.FortranMicro")
    def test_batch_glob_calls_from_hdf5_for_each_run(
        self, mock_cls, runner, hdf5_dir
    ):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            ["run-micro", str(hdf5_dir), "--executable", "/bin/micro.exe"],
        )
        assert result.exit_code == 0, result.output
        assert mock_cls.from_hdf5.call_count == 2

    @patch("lysis.execution.fortran_micro.FortranMicro")
    def test_batch_output_mentions_run_codes(
        self, mock_cls, runner, experiment_dir
    ):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            ["run-micro", str(experiment_dir), "--executable", "/bin/micro.exe"],
        )
        assert result.exit_code == 0, result.output
        # Rich may wrap long paths across lines; join output to check content
        flat = result.output.replace("\n", "")
        assert "run-01" in flat
        assert "run-02" in flat

    def test_file_code_with_directory_raises_error(self, runner, experiment_dir):
        result = runner.invoke(
            cli,
            [
                "run-micro", str(experiment_dir),
                "--executable", "/bin/micro.exe",
                "--file-code", "_x",
            ],
        )
        assert result.exit_code != 0
        assert "file-code" in result.output.lower() or "directory" in result.output.lower()

    def test_empty_directory_raises_error(self, runner, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        result = runner.invoke(
            cli,
            ["run-micro", str(empty), "--executable", "/bin/micro.exe"],
        )
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# Batch Slurm dispatch (experiment folder)
# ---------------------------------------------------------------------------


class TestRunMicroExperimentSlurm:
    @patch("lysis.tools.slurm.submit_micro_slurm_job", side_effect=[10001, 10002])
    def test_batch_slurm_calls_submit_for_each_run(
        self, mock_submit, runner, experiment_dir
    ):
        result = runner.invoke(
            cli,
            [
                "run-micro", str(experiment_dir),
                "--executable", "/bin/micro.exe",
                "--slurm",
            ],
        )
        assert result.exit_code == 0, result.output
        assert mock_submit.call_count == 2

    @patch("lysis.tools.slurm.submit_micro_slurm_job", side_effect=[10001, 10002])
    def test_batch_slurm_prints_all_job_ids(
        self, mock_submit, runner, experiment_dir
    ):
        result = runner.invoke(
            cli,
            [
                "run-micro", str(experiment_dir),
                "--executable", "/bin/micro.exe",
                "--slurm",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "10001" in result.output
        assert "10002" in result.output

    @patch("lysis.tools.slurm.submit_micro_slurm_job", side_effect=[10001, 10002])
    def test_batch_slurm_partition_forwarded_to_all_calls(
        self, mock_submit, runner, experiment_dir
    ):
        result = runner.invoke(
            cli,
            [
                "run-micro", str(experiment_dir),
                "--executable", "/bin/micro.exe",
                "--slurm",
                "--partition", "long",
            ],
        )
        assert result.exit_code == 0, result.output
        for c in mock_submit.call_args_list:
            assert c.kwargs.get("partition") == "long"

    @patch("lysis.tools.slurm.submit_micro_slurm_job", side_effect=[20001, 20002])
    def test_batch_glob_slurm_submits_for_each_h5(
        self, mock_submit, runner, hdf5_dir
    ):
        result = runner.invoke(
            cli,
            [
                "run-micro", str(hdf5_dir),
                "--executable", "/bin/micro.exe",
                "--slurm",
            ],
        )
        assert result.exit_code == 0, result.output
        assert mock_submit.call_count == 2
