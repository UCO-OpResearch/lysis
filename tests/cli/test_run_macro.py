"""Tests for ``lysis run-macro`` CLI command."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import h5py
import pytest
from click.testing import CliRunner

from lysis.cli import cli
from lysis.config.constants import CONST
from lysis.config.parameters import MacroParameters, MicroParameters


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------


def _write_macro_hdf5(path: Path) -> None:
    """Write a minimal v2.0.0 HDF5 file with both Micro- and MacroParameters."""
    mp = MicroParameters()
    macp = MacroParameters(micro_params=mp)
    with h5py.File(str(path), "w") as f:
        f.attrs[CONST.DATASPEC_VERSION_ATTR] = "v2.0.0"
        grp = f.require_group("micro_data")
        for k, v in mp.to_basedict().items():
            grp.attrs[k] = str(v) if not isinstance(v, (int, float, bool)) else v
        grp2 = f.require_group("macro_data")
        for k, v in macp.to_basedict().items():
            grp2.attrs[k] = str(v) if not isinstance(v, (int, float, bool)) else v


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def macro_hdf5(tmp_path):
    """Minimal HDF5 file with both params; returns Path."""
    h5_path = tmp_path / "run-01.h5"
    _write_macro_hdf5(h5_path)
    return h5_path


# ---------------------------------------------------------------------------
# Help
# ---------------------------------------------------------------------------


class TestRunMacroHelp:
    def test_help_exits_0(self, runner):
        result = runner.invoke(cli, ["run-macro", "--help"])
        assert result.exit_code == 0

    def test_help_mentions_executable(self, runner):
        result = runner.invoke(cli, ["run-macro", "--help"])
        assert "--executable" in result.output

    def test_help_mentions_slurm(self, runner):
        result = runner.invoke(cli, ["run-macro", "--help"])
        assert "--slurm" in result.output

    def test_help_mentions_keep_tmpdir(self, runner):
        result = runner.invoke(cli, ["run-macro", "--help"])
        assert "--keep-tmpdir" in result.output


# ---------------------------------------------------------------------------
# Missing required argument
# ---------------------------------------------------------------------------


class TestRunMacroMissingArgs:
    def test_missing_executable_exits_nonzero(self, runner, macro_hdf5):
        result = runner.invoke(cli, ["run-macro", str(macro_hdf5)])
        assert result.exit_code != 0

    def test_missing_hdf5_path_exits_nonzero(self, runner):
        result = runner.invoke(
            cli, ["run-macro", "--executable", "/bin/macro.exe"]
        )
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# Non-Slurm (local) execution
# ---------------------------------------------------------------------------


class TestRunMacroLocal:
    @patch("lysis.execution.codeutil.FortranMacro")
    def test_calls_from_hdf5(self, mock_cls, runner, macro_hdf5):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            ["run-macro", str(macro_hdf5), "--executable", "/bin/macro.exe"],
        )
        assert result.exit_code == 0, result.output
        mock_cls.from_hdf5.assert_called_once()
        call_args = mock_cls.from_hdf5.call_args
        assert str(macro_hdf5) in str(call_args.args[0])
        assert call_args.args[1] == "/bin/macro.exe"

    @patch("lysis.execution.codeutil.FortranMacro")
    def test_calls_run_full(self, mock_cls, runner, macro_hdf5):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            ["run-macro", str(macro_hdf5), "--executable", "/bin/macro.exe"],
        )
        assert result.exit_code == 0, result.output
        mock_fm.run_full.assert_called_once()

    @patch("lysis.execution.codeutil.FortranMacro")
    def test_keep_tmpdir_flag_passed_to_run_full(self, mock_cls, runner, macro_hdf5):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "/bin/macro.exe",
                "--keep-tmpdir",
            ],
        )
        assert result.exit_code == 0, result.output
        call_kwargs = mock_fm.run_full.call_args.kwargs
        assert call_kwargs.get("keep_tmpdir") is True

    @patch("lysis.execution.codeutil.FortranMacro")
    def test_file_code_passed_to_from_hdf5(self, mock_cls, runner, macro_hdf5):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "/bin/macro.exe",
                "--file-code", "_code",
            ],
        )
        assert result.exit_code == 0, result.output
        call_kwargs = mock_cls.from_hdf5.call_args.kwargs
        assert call_kwargs.get("out_file_code") == "_code"

    @patch("lysis.execution.codeutil.FortranMacro")
    def test_output_mentions_hdf5_path(self, mock_cls, runner, macro_hdf5):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            ["run-macro", str(macro_hdf5), "--executable", "/bin/macro.exe"],
        )
        assert result.exit_code == 0, result.output
        assert "run-01" in result.output


# ---------------------------------------------------------------------------
# Slurm dispatch
# ---------------------------------------------------------------------------


class TestRunMacroSlurm:
    @patch("lysis.tools.slurm.submit_macro_slurm_job", return_value=12345)
    def test_slurm_calls_submit(self, mock_submit, runner, macro_hdf5):
        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "/bin/macro.exe",
                "--slurm",
            ],
        )
        assert result.exit_code == 0, result.output
        mock_submit.assert_called_once()

    @patch("lysis.tools.slurm.submit_macro_slurm_job", return_value=12345)
    def test_slurm_prints_job_id(self, mock_submit, runner, macro_hdf5):
        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "/bin/macro.exe",
                "--slurm",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "12345" in result.output

    @patch("lysis.tools.slurm.submit_macro_slurm_job", return_value=1)
    def test_partition_forwarded_to_submit(self, mock_submit, runner, macro_hdf5):
        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "/bin/macro.exe",
                "--slurm",
                "--partition", "long",
            ],
        )
        assert result.exit_code == 0, result.output
        call_kwargs = mock_submit.call_args.kwargs
        assert call_kwargs.get("partition") == "long"

    @patch("lysis.tools.slurm.submit_macro_slurm_job", return_value=1)
    def test_staging_root_forwarded(self, mock_submit, runner, macro_hdf5, tmp_path):
        staging = tmp_path / "staging"
        staging.mkdir()
        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "/bin/macro.exe",
                "--slurm",
                "--staging-root", str(staging),
            ],
        )
        assert result.exit_code == 0, result.output
        call_kwargs = mock_submit.call_args.kwargs
        assert call_kwargs.get("staging_root") == str(staging)

    @patch("lysis.tools.slurm.submit_macro_slurm_job", return_value=1)
    def test_fast_tmp_root_forwarded(self, mock_submit, runner, macro_hdf5):
        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "/bin/macro.exe",
                "--slurm",
                "--fast-tmp-root", "/nvme/scratch",
            ],
        )
        assert result.exit_code == 0, result.output
        call_kwargs = mock_submit.call_args.kwargs
        assert call_kwargs.get("fast_tmp_root") == "/nvme/scratch"

    @patch("lysis.tools.slurm.submit_macro_slurm_job", return_value=1)
    def test_keep_tmpdir_forwarded_to_submit(self, mock_submit, runner, macro_hdf5):
        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "/bin/macro.exe",
                "--slurm",
                "--keep-tmpdir",
            ],
        )
        assert result.exit_code == 0, result.output
        call_kwargs = mock_submit.call_args.kwargs
        assert call_kwargs.get("keep_tmpdir") is True

    @patch("lysis.execution.codeutil.FortranMacro")
    def test_partition_without_slurm_does_not_crash(self, mock_cls, runner,
                                                     macro_hdf5):
        """--partition without --slurm should be silently ignored."""
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "/bin/macro.exe",
                "--partition", "long",
            ],
        )
        assert result.exit_code == 0, result.output
        mock_fm.run_full.assert_called_once()
