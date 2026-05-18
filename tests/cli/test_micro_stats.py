"""Tests for ``lysis micro-stats`` CLI command."""

from unittest.mock import patch

import pytest
from click.testing import CliRunner

from lysis.cli import cli

# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------


def _make_stats(
    fibers_degraded=1234,
    lysis_time_mean=42.0,
    lysis_time_std=3.5,
    lysis_time_median=41.5,
    tpa_leaving_mean=15.0,
    tpa_leaving_std=2.0,
    tpa_leaving_median=14.5,
):
    """Return a mock stats dict as returned by _load_run_micro_stats."""
    return {
        "fibers_degraded": fibers_degraded,
        "lysis_time_mean": lysis_time_mean,
        "lysis_time_std": lysis_time_std,
        "lysis_time_median": lysis_time_median,
        "tpa_leaving_mean": tpa_leaving_mean,
        "tpa_leaving_std": tpa_leaving_std,
        "tpa_leaving_median": tpa_leaving_median,
    }


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def mock_stats():
    return _make_stats()


# ---------------------------------------------------------------------------
# Help
# ---------------------------------------------------------------------------


class TestMicroStatsHelp:
    def test_help_exits_0(self, runner):
        result = runner.invoke(cli, ["micro-stats", "--help"])
        assert result.exit_code == 0

    def test_help_mentions_markdown(self, runner):
        result = runner.invoke(cli, ["micro-stats", "--help"])
        assert "--markdown" in result.output

    def test_help_mentions_no_progress(self, runner):
        result = runner.invoke(cli, ["micro-stats", "--help"])
        assert "--no-progress" in result.output

    def test_help_mentions_sort(self, runner):
        result = runner.invoke(cli, ["micro-stats", "--help"])
        assert "--sort" in result.output


# ---------------------------------------------------------------------------
# CLI integration — single-file mode
# ---------------------------------------------------------------------------


class TestMicroStatsMarkdownSingleFile:
    @patch("lysis.cli.micro_stats._load_run_micro_stats")
    def test_markdown_to_stdout(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(cli, ["micro-stats", str(h5), "--markdown", "-"])

        assert result.exit_code == 0
        assert "## run_A" in result.output

    @patch("lysis.cli.micro_stats._load_run_micro_stats")
    def test_markdown_to_stdout_has_table(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_B.h5"
        h5.touch()

        result = runner.invoke(cli, ["micro-stats", str(h5), "--markdown", "-"])

        assert result.exit_code == 0
        assert "| Metric | Value |" in result.output

    @patch("lysis.cli.micro_stats._load_run_micro_stats")
    def test_markdown_to_stdout_contains_fibers_degraded(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_C.h5"
        h5.touch()

        result = runner.invoke(cli, ["micro-stats", str(h5), "--markdown", "-"])

        assert result.exit_code == 0
        assert "Fibers Degraded" in result.output
        assert "1,234" in result.output

    @patch("lysis.cli.micro_stats._load_run_micro_stats")
    def test_markdown_to_file(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_D.h5"
        h5.touch()
        out_file = tmp_path / "stats.md"

        result = runner.invoke(cli, ["micro-stats", str(h5), "--markdown", str(out_file)])

        assert result.exit_code == 0
        assert out_file.exists()

    @patch("lysis.cli.micro_stats._load_run_micro_stats")
    def test_markdown_file_content(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_E.h5"
        h5.touch()
        out_file = tmp_path / "stats.md"

        runner.invoke(cli, ["micro-stats", str(h5), "--markdown", str(out_file)])

        content = out_file.read_text()
        assert "## run_E" in content
        assert "| Metric | Value |" in content

    @patch("lysis.cli.micro_stats._load_run_micro_stats")
    def test_markdown_file_confirmation_message(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_F.h5"
        h5.touch()
        out_file = tmp_path / "stats.md"

        result = runner.invoke(cli, ["micro-stats", str(h5), "--markdown", str(out_file)])

        assert "Markdown written to" in result.output

    @patch("lysis.cli.micro_stats._load_run_micro_stats")
    def test_no_markdown_flag_gives_normal_output(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_G.h5"
        h5.touch()

        result = runner.invoke(cli, ["micro-stats", str(h5), "--no-progress"])

        assert result.exit_code == 0
        assert "## run_G" not in result.output
        assert "run_G" in result.output
        assert "Fibers Degraded" in result.output


# ---------------------------------------------------------------------------
# CLI integration — directory mode
# ---------------------------------------------------------------------------


class TestMicroStatsMarkdownDirectory:
    @patch("lysis.cli.micro_stats._load_run_micro_stats")
    def test_markdown_to_stdout(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        (tmp_path / "run_1.h5").touch()
        (tmp_path / "run_2.h5").touch()

        result = runner.invoke(cli, ["micro-stats", str(tmp_path), "--markdown", "-"])

        assert result.exit_code == 0
        assert "|" in result.output

    @patch("lysis.cli.micro_stats._load_run_micro_stats")
    def test_markdown_to_stdout_run_codes_as_rows(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        (tmp_path / "run_A.h5").touch()
        (tmp_path / "run_B.h5").touch()

        result = runner.invoke(cli, ["micro-stats", str(tmp_path), "--markdown", "-"])

        assert result.exit_code == 0
        assert "run_A" in result.stdout
        assert "run_B" in result.stdout
        assert result.stdout.strip().startswith("| Run |")

    @patch("lysis.cli.micro_stats._load_run_micro_stats")
    def test_markdown_to_stdout_metrics_in_header(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        (tmp_path / "run_A.h5").touch()

        result = runner.invoke(cli, ["micro-stats", str(tmp_path), "--markdown", "-"])

        assert "Fibers Degraded" in result.output
        assert "Mean Lysis Time (min)" in result.output

    @patch("lysis.cli.micro_stats._load_run_micro_stats")
    def test_markdown_to_file(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        (tmp_path / "run_1.h5").touch()
        out_file = tmp_path / "stats.md"

        result = runner.invoke(
            cli, ["micro-stats", str(tmp_path), "--markdown", str(out_file)]
        )

        assert result.exit_code == 0
        assert out_file.exists()
        content = out_file.read_text()
        assert "run_1" in content

    @patch("lysis.cli.micro_stats._load_run_micro_stats")
    def test_no_markdown_flag_gives_normal_output(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        (tmp_path / "run_1.h5").touch()

        result = runner.invoke(cli, ["micro-stats", str(tmp_path), "--no-progress"])

        assert result.exit_code == 0
        assert not result.output.strip().startswith("#")
        assert "run_1" in result.output
