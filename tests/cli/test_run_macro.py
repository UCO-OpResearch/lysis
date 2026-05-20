"""Tests for ``lysis run-macro`` CLI command."""

import json
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
    """Write a minimal v2.0.0 HDF5 with default Micro and MacroParameters."""
    mp = MicroParameters()
    mcp = MacroParameters(micro_params=mp)
    with h5py.File(str(path), "w") as f:
        f.attrs[CONST.DATASPEC_VERSION_ATTR] = "v2.0.0"
        micro_grp = f.require_group("micro_data")
        for k, v in mp.to_basedict().items():
            micro_grp.attrs[k] = str(v) if not isinstance(v, (int, float, bool)) else v
        macro_grp = f.require_group("macro_data")
        for k, v in mcp.to_basedict().items():
            macro_grp.attrs[k] = str(v) if not isinstance(v, (int, float, bool)) else v


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def macro_hdf5(tmp_path):
    """Minimal HDF5 file with both Micro and MacroParameters."""
    h5_path = tmp_path / "run-01.h5"
    _write_macro_hdf5(h5_path)
    return h5_path


def _make_macro_experiment_dir(parent: Path, run_codes: list) -> Path:
    """Create a minimal experiment folder with macro HDF5 files and experiment.json."""
    exp_dir = parent / "my-experiment"
    exp_dir.mkdir()
    for rc in run_codes:
        _write_macro_hdf5(exp_dir / f"{rc}.h5")
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
    """Experiment folder with experiment.json and two macro HDF5 files."""
    return _make_macro_experiment_dir(tmp_path, ["run-01", "run-02"])


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

    def test_help_mentions_in_file_code(self, runner):
        result = runner.invoke(cli, ["run-macro", "--help"])
        assert "--in-file-code" in result.output

    def test_help_mentions_out_file_code(self, runner):
        result = runner.invoke(cli, ["run-macro", "--help"])
        assert "--out-file-code" in result.output


# ---------------------------------------------------------------------------
# Missing required argument
# ---------------------------------------------------------------------------


class TestRunMacroMissingArgs:
    def test_missing_executable_exits_nonzero(self, runner, macro_hdf5):
        result = runner.invoke(cli, ["run-macro", str(macro_hdf5)])
        assert result.exit_code != 0

    def test_missing_hdf5_path_exits_nonzero(self, runner, tmp_path):
        result = runner.invoke(
            cli, ["run-macro", "--executable", "/bin/macro.exe"]
        )
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# Non-Slurm (local) execution
# ---------------------------------------------------------------------------


class TestRunMacroLocal:
    @patch("lysis.execution.fortran_macro.FortranMacro")
    def test_calls_from_hdf5(self, mock_cls, runner, macro_hdf5):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            ["run-macro", str(macro_hdf5), "--executable", "/bin/macro.exe"],
        )
        assert result.exit_code == 0, result.output
        mock_cls.from_hdf5.assert_called_once()

    @patch("lysis.execution.fortran_macro.FortranMacro")
    def test_calls_run_full(self, mock_cls, runner, macro_hdf5):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            ["run-macro", str(macro_hdf5), "--executable", "/bin/macro.exe"],
        )
        assert result.exit_code == 0, result.output
        mock_fm.run_full.assert_called_once()

    @patch("lysis.execution.fortran_macro.FortranMacro")
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

    @patch("lysis.execution.fortran_macro.FortranMacro")
    def test_out_file_code_passed_to_from_hdf5(self, mock_cls, runner, macro_hdf5):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "/bin/macro.exe",
                "--out-file-code", "_out",
            ],
        )
        assert result.exit_code == 0, result.output
        call_kwargs = mock_cls.from_hdf5.call_args.kwargs
        assert call_kwargs.get("out_file_code") == "_out"

    @patch("lysis.execution.fortran_macro.FortranMacro")
    def test_in_file_code_passed_to_from_hdf5(self, mock_cls, runner, macro_hdf5):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "/bin/macro.exe",
                "--in-file-code", "_in",
            ],
        )
        assert result.exit_code == 0, result.output
        call_kwargs = mock_cls.from_hdf5.call_args.kwargs
        assert call_kwargs.get("in_file_code") == "_in"

    @patch("lysis.execution.fortran_macro.FortranMacro")
    def test_output_mentions_hdf5_path(self, mock_cls, runner, macro_hdf5):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            ["run-macro", str(macro_hdf5), "--executable", "/bin/macro.exe"],
        )
        assert result.exit_code == 0, result.output
        assert "run-01" in result.output

    @patch("lysis.execution.fortran_macro.FortranMacro")
    def test_exits_nonzero_when_from_hdf5_raises_valueerror(
        self, mock_cls, runner, macro_hdf5
    ):
        """A ValueError from from_hdf5 (e.g. wrong HDF5 state) must exit non-zero."""
        mock_cls.from_hdf5.side_effect = ValueError("macroscale simulation has already been run")
        result = runner.invoke(
            cli,
            ["run-macro", str(macro_hdf5), "--executable", "/bin/macro.exe"],
        )
        assert result.exit_code != 0

    @patch("lysis.execution.fortran_macro.FortranMacro")
    def test_error_message_shown_when_from_hdf5_raises(
        self, mock_cls, runner, macro_hdf5
    ):
        """Error message from ValueError must appear in CLI output."""
        mock_cls.from_hdf5.side_effect = ValueError("macroscale simulation has already been run")
        result = runner.invoke(
            cli,
            ["run-macro", str(macro_hdf5), "--executable", "/bin/macro.exe"],
        )
        assert "macroscale simulation has already been run" in result.output


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
    def test_in_code_forwarded_to_submit(self, mock_submit, runner, macro_hdf5):
        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "/bin/macro.exe",
                "--slurm",
                "--in-file-code", "_in",
            ],
        )
        assert result.exit_code == 0, result.output
        call_kwargs = mock_submit.call_args.kwargs
        assert call_kwargs.get("in_code") == "_in"

    @patch("lysis.tools.slurm.submit_macro_slurm_job", return_value=1)
    def test_out_code_forwarded_to_submit(self, mock_submit, runner, macro_hdf5):
        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "/bin/macro.exe",
                "--slurm",
                "--out-file-code", "_out",
            ],
        )
        assert result.exit_code == 0, result.output
        call_kwargs = mock_submit.call_args.kwargs
        assert call_kwargs.get("out_code") == "_out"

    @patch("lysis.tools.slurm.submit_macro_slurm_job", return_value=1)
    def test_compiler_default_is_forwarded(self, mock_submit, runner, macro_hdf5):
        """Without --compiler, the default module spec must be forwarded."""
        from lysis.tools.slurm import DEFAULT_COMPILER_MODULE
        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "/bin/macro.exe",
                "--slurm",
            ],
        )
        assert result.exit_code == 0, result.output
        assert (
            mock_submit.call_args.kwargs.get("compiler_module")
            == DEFAULT_COMPILER_MODULE
        )

    @patch("lysis.tools.slurm.submit_macro_slurm_job", return_value=1)
    def test_compiler_override_is_forwarded(self, mock_submit, runner, macro_hdf5):
        """--compiler=foo/2024 must reach submit_macro_slurm_job."""
        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "/bin/macro.exe",
                "--slurm",
                "--compiler", "intel-compilers/2024",
            ],
        )
        assert result.exit_code == 0, result.output
        assert (
            mock_submit.call_args.kwargs.get("compiler_module")
            == "intel-compilers/2024"
        )

    @patch("lysis.tools.slurm.submit_macro_slurm_job", return_value=1)
    def test_sbatch_default_is_empty_overrides(
        self, mock_submit, runner, macro_hdf5
    ):
        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "/bin/macro.exe",
                "--slurm",
            ],
        )
        assert result.exit_code == 0, result.output
        assert mock_submit.call_args.kwargs.get("sbatch_overrides") == {}

    @patch("lysis.tools.slurm.submit_macro_slurm_job", return_value=1)
    def test_sbatch_tokens_parsed_and_forwarded(
        self, mock_submit, runner, macro_hdf5
    ):
        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "/bin/macro.exe",
                "--slurm",
                "--sbatch", "mem=4G",
                "--sbatch", "hold",
                "--sbatch", "^exclusive=user",
            ],
        )
        assert result.exit_code == 0, result.output
        assert mock_submit.call_args.kwargs.get("sbatch_overrides") == {
            "mem": "4G",
            "hold": "",
            "exclusive=user": None,
        }

    def test_sbatch_malformed_token_rejected(self, runner, macro_hdf5):
        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "/bin/macro.exe",
                "--slurm",
                "--sbatch", "^",
            ],
        )
        assert result.exit_code != 0
        assert "--sbatch" in result.output

    def test_help_mentions_sbatch(self, runner):
        result = runner.invoke(cli, ["run-macro", "--help"])
        assert result.exit_code == 0
        assert "--sbatch" in result.output


# ---------------------------------------------------------------------------
# Batch (experiment / directory)
# ---------------------------------------------------------------------------


class TestRunMacroExperiment:
    @patch("lysis.execution.fortran_macro.FortranMacro")
    def test_batch_calls_from_hdf5_for_each_run(
        self, mock_cls, runner, experiment_dir
    ):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            ["run-macro", str(experiment_dir), "--executable", "/bin/macro.exe"],
        )
        assert result.exit_code == 0, result.output
        assert mock_cls.from_hdf5.call_count == 2

    @patch("lysis.execution.fortran_macro.FortranMacro")
    def test_batch_calls_run_full_for_each_run(
        self, mock_cls, runner, experiment_dir
    ):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            ["run-macro", str(experiment_dir), "--executable", "/bin/macro.exe"],
        )
        assert result.exit_code == 0, result.output
        assert mock_fm.run_full.call_count == 2

    @patch("lysis.execution.fortran_macro.FortranMacro")
    def test_batch_output_mentions_run_codes(
        self, mock_cls, runner, experiment_dir
    ):
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            ["run-macro", str(experiment_dir), "--executable", "/bin/macro.exe"],
        )
        assert result.exit_code == 0, result.output
        flat = result.output.replace("\n", "")
        assert "run-01" in flat
        assert "run-02" in flat

    def test_out_file_code_with_directory_raises_error(
        self, runner, experiment_dir
    ):
        result = runner.invoke(
            cli,
            [
                "run-macro", str(experiment_dir),
                "--executable", "/bin/macro.exe",
                "--out-file-code", "_x",
            ],
        )
        assert result.exit_code != 0
        assert (
            "out-file-code" in result.output.lower()
            or "directory" in result.output.lower()
        )

    @patch("lysis.tools.slurm.submit_macro_slurm_job", return_value=1)
    def test_slurm_batch_submits_for_each_run(
        self, mock_submit, runner, experiment_dir
    ):
        result = runner.invoke(
            cli,
            [
                "run-macro", str(experiment_dir),
                "--executable", "/bin/macro.exe",
                "--slurm",
            ],
        )
        assert result.exit_code == 0, result.output
        assert mock_submit.call_count == 2


# ---------------------------------------------------------------------------
# --fortran-commit (historical-build workflow)
# ---------------------------------------------------------------------------


class TestRunMacroFortranCommit:
    def test_bad_ref_exits_nonzero(self, runner, macro_hdf5):
        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "macro_diffuse_into_and_along__internal",
                "--fortran-commit", "this-ref-does-not-exist-deadbeef",
            ],
        )
        assert result.exit_code != 0
        assert "does not resolve" in result.output

    @patch("lysis.execution.fortran_macro.FortranMacro")
    @patch("lysis.execution.historical_build.build_historical_binary")
    def test_historical_attrs_threaded_to_runner(
        self, mock_build, mock_cls, runner, macro_hdf5, tmp_path
    ):
        from contextlib import contextmanager

        fake_binary = tmp_path / "bin" / "macro_diffuse_into_and_along__internal"
        fake_binary.parent.mkdir(parents=True)
        fake_binary.write_text("not a real binary")
        prov = {
            "binary_commit": "c" * 40,
            "binary_dirty": "clean",
            "binary_compiler": "GCC 11.4.0",
            "binary_source": "historical:" + "c" * 40,
        }

        @contextmanager
        def fake_cm(*args, **kwargs):
            # Local mode: compiler_module should be None.
            assert kwargs.get("compiler_module") is None
            yield fake_binary, prov

        mock_build.side_effect = fake_cm
        mock_fm = MagicMock()
        mock_cls.from_hdf5.return_value = mock_fm

        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "macro_diffuse_into_and_along__internal",
                "--fortran-commit", "HEAD",
            ],
        )
        assert result.exit_code == 0, result.output
        call_kwargs = mock_cls.from_hdf5.call_args.kwargs
        assert call_kwargs.get("historical_binary_attrs") == prov
        assert call_kwargs.get("skip_binary_verification") is True

    @patch("lysis.tools.slurm.submit_macro_slurm_job", return_value=99)
    @patch("lysis.execution.historical_build.build_historical_binary")
    def test_historical_attrs_threaded_to_slurm(
        self, mock_build, mock_submit, runner, macro_hdf5, tmp_path
    ):
        from contextlib import contextmanager

        fake_binary = tmp_path / "bin" / "macro_diffuse_into_and_along__internal"
        fake_binary.parent.mkdir(parents=True)
        fake_binary.write_text("not a real binary")
        prov = {
            "binary_commit": "d" * 40,
            "binary_dirty": "clean",
            "binary_compiler": "Intel(R) Fortran",
            "binary_source": "historical:" + "d" * 40,
        }

        @contextmanager
        def fake_cm(*args, **kwargs):
            assert kwargs.get("compiler_module") == "intel-compilers/2024"
            yield fake_binary, prov

        mock_build.side_effect = fake_cm
        result = runner.invoke(
            cli,
            [
                "run-macro", str(macro_hdf5),
                "--executable", "macro_diffuse_into_and_along__internal",
                "--fortran-commit", "HEAD",
                "--slurm",
                "--compiler", "intel-compilers/2024",
            ],
        )
        assert result.exit_code == 0, result.output
        call_kwargs = mock_submit.call_args.kwargs
        assert call_kwargs.get("historical_binary_attrs") == prov
