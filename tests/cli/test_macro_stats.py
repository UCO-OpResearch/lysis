"""Tests for ``lysis macro-stats`` CLI command."""

from unittest.mock import patch

import pandas as pd
import pytest
from click.testing import CliRunner

from lysis.cli import cli

# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------

_METRICS = [
    "Degradation rate (%/min)",
    "Lysis lag time (min)",
    "Time to full clot degradation (min)",
    "Percent of molecules that reached the back row",
    "First passage time (min)",
    "Front Velocity (microns/min)",
]


def _make_stats(mean=1.234, std=0.056):
    """Return a mock stats Series as returned by compute_run_statistics."""
    index = pd.MultiIndex.from_tuples(
        [(m, s) for m in _METRICS for s in ["Mean", "Standard Deviation"]]
    )
    values = [mean if s == "Mean" else std for m in _METRICS for s in ["Mean", "Standard Deviation"]]
    return pd.Series(values, index=index)


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def mock_stats():
    return _make_stats()


# ---------------------------------------------------------------------------
# Help
# ---------------------------------------------------------------------------


class TestSummarizeHelp:
    def test_help_exits_0(self, runner):
        result = runner.invoke(cli, ["macro-stats", "--help"])
        assert result.exit_code == 0

    def test_help_mentions_markdown(self, runner):
        result = runner.invoke(cli, ["macro-stats", "--help"])
        assert "--markdown" in result.output

    def test_help_mentions_no_progress(self, runner):
        result = runner.invoke(cli, ["macro-stats", "--help"])
        assert "--no-progress" in result.output


# ---------------------------------------------------------------------------
# Markdown helpers (unit tests — no CLI invocation needed)
# ---------------------------------------------------------------------------


class TestMdTable:
    def test_header_row(self):
        from lysis.cli.macro_stats import _md_table

        out = _md_table(["A", "B"], [["1", "2"]])
        lines = out.splitlines()
        assert lines[0] == "| A | B |"

    def test_separator_row(self):
        from lysis.cli.macro_stats import _md_table

        out = _md_table(["A", "B"], [["1", "2"]])
        lines = out.splitlines()
        assert lines[1] == "| --- | --- |"

    def test_data_row(self):
        from lysis.cli.macro_stats import _md_table

        out = _md_table(["A", "B"], [["hello", "world"]])
        assert "| hello | world |" in out

    def test_multiple_rows(self):
        from lysis.cli.macro_stats import _md_table

        out = _md_table(["X"], [["r1"], ["r2"], ["r3"]])
        assert out.count("| r") == 3


class TestSummarizeToMarkdown:
    def test_single_file_heading(self, mock_stats):
        from lysis.cli.macro_stats import _summarize_to_markdown

        md = _summarize_to_markdown({"run_A": mock_stats}, _METRICS, "run_A")
        assert md.startswith("## run_A")

    def test_single_file_metric_in_table(self, mock_stats):
        from lysis.cli.macro_stats import _summarize_to_markdown

        md = _summarize_to_markdown({"run_A": mock_stats}, _METRICS, "run_A")
        assert "Degradation rate (%/min)" in md

    def test_single_file_values_formatted(self, mock_stats):
        from lysis.cli.macro_stats import _summarize_to_markdown

        md = _summarize_to_markdown({"run_A": mock_stats}, _METRICS, "run_A")
        assert "1.234" in md
        assert "0.056" in md

    def test_single_file_three_columns(self, mock_stats):
        from lysis.cli.macro_stats import _summarize_to_markdown

        md = _summarize_to_markdown({"run_A": mock_stats}, _METRICS, "run_A")
        header_line = [l for l in md.splitlines() if "Metric" in l][0]
        assert header_line.count("|") == 4  # | Metric | Mean | Std Dev |

    def test_directory_run_codes_in_rows(self):
        """In directory mode, run codes appear as row values (first column)."""
        from lysis.cli.macro_stats import _summarize_to_markdown

        rows = {"run_X": _make_stats(1.0, 0.1), "run_Y": _make_stats(2.0, 0.2)}
        md = _summarize_to_markdown(rows, _METRICS)
        assert "run_X" in md
        assert "run_Y" in md
        # First column header is "Run"
        assert md.splitlines()[0].startswith("| Run |")

    def test_directory_mean_pm_std_format(self):
        from lysis.cli.macro_stats import _summarize_to_markdown

        rows = {"run_X": _make_stats(1.234, 0.056)}
        md = _summarize_to_markdown(rows, _METRICS)
        # Should contain "1.234 ± 0.056" in a cell
        assert "1.234" in md
        assert "\u00b1" in md  # ±

    def test_directory_all_metrics_present(self):
        from lysis.cli.macro_stats import _summarize_to_markdown

        rows = {"run_X": _make_stats()}
        md = _summarize_to_markdown(rows, _METRICS)
        for metric in _METRICS:
            assert metric in md


# ---------------------------------------------------------------------------
# CLI integration — single-file mode
# ---------------------------------------------------------------------------


class TestSummarizeMarkdownSingleFile:
    @patch("lysis.cli.macro_stats._load_run_stats")
    def test_markdown_to_stdout(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(cli, ["macro-stats", str(h5), "--markdown", "-"])

        assert result.exit_code == 0
        assert "## run_A" in result.output

    @patch("lysis.cli.macro_stats._load_run_stats")
    def test_markdown_to_stdout_has_table(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_B.h5"
        h5.touch()

        result = runner.invoke(cli, ["macro-stats", str(h5), "--markdown", "-"])

        assert result.exit_code == 0
        assert "| Metric | Mean | Std Dev |" in result.output

    @patch("lysis.cli.macro_stats._load_run_stats")
    def test_markdown_to_stdout_contains_metrics(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_C.h5"
        h5.touch()

        result = runner.invoke(cli, ["macro-stats", str(h5), "--markdown", "-"])

        assert result.exit_code == 0
        assert "Degradation rate (%/min)" in result.output

    @patch("lysis.cli.macro_stats._load_run_stats")
    def test_markdown_to_file(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_D.h5"
        h5.touch()
        out_file = tmp_path / "report.md"

        result = runner.invoke(
            cli, ["macro-stats", str(h5), "--markdown", str(out_file)]
        )

        assert result.exit_code == 0
        assert out_file.exists()

    @patch("lysis.cli.macro_stats._load_run_stats")
    def test_markdown_file_content(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_E.h5"
        h5.touch()
        out_file = tmp_path / "report.md"

        runner.invoke(cli, ["macro-stats", str(h5), "--markdown", str(out_file)])

        content = out_file.read_text()
        assert "## run_E" in content
        assert "| Metric | Mean | Std Dev |" in content

    @patch("lysis.cli.macro_stats._load_run_stats")
    def test_markdown_file_confirmation_message(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_F.h5"
        h5.touch()
        out_file = tmp_path / "report.md"

        result = runner.invoke(
            cli, ["macro-stats", str(h5), "--markdown", str(out_file)]
        )

        assert "Markdown written to" in result.output

    @patch("lysis.cli.macro_stats._load_run_stats")
    def test_no_markdown_flag_gives_normal_output(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_G.h5"
        h5.touch()

        result = runner.invoke(cli, ["macro-stats", str(h5), "--no-progress"])

        assert result.exit_code == 0
        # Normal output should NOT start with a markdown heading
        assert "## run_G" not in result.output
        # But the run code should appear
        assert "run_G" in result.output


# ---------------------------------------------------------------------------
# CLI integration — directory mode
# ---------------------------------------------------------------------------


class TestSummarizeMarkdownDirectory:
    @patch("lysis.cli.macro_stats._load_run_stats")
    def test_markdown_to_stdout(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        (tmp_path / "run_1.h5").touch()
        (tmp_path / "run_2.h5").touch()

        result = runner.invoke(
            cli, ["macro-stats", str(tmp_path), "--markdown", "-"]
        )

        assert result.exit_code == 0
        assert "|" in result.output  # markdown table

    @patch("lysis.cli.macro_stats._load_run_stats")
    def test_markdown_to_stdout_run_codes_in_rows(self, mock_load, runner, mock_stats, tmp_path):
        """Run codes appear as row values (first column) in directory markdown."""
        mock_load.return_value = mock_stats
        (tmp_path / "run_A.h5").touch()
        (tmp_path / "run_B.h5").touch()

        result = runner.invoke(
            cli, ["macro-stats", str(tmp_path), "--markdown", "-"]
        )

        assert result.exit_code == 0
        assert "run_A" in result.output
        assert "run_B" in result.output
        # First line is the header row beginning with "| Run |"
        assert result.output.strip().startswith("| Run |")

    @patch("lysis.cli.macro_stats._load_run_stats")
    def test_markdown_to_stdout_full_metric_names(self, mock_load, runner, mock_stats, tmp_path):
        """Directory markdown uses full metric names, not abbreviated ones."""
        mock_load.return_value = mock_stats
        (tmp_path / "run_A.h5").touch()

        result = runner.invoke(
            cli, ["macro-stats", str(tmp_path), "--markdown", "-"]
        )

        assert "Degradation rate (%/min)" in result.output

    @patch("lysis.cli.macro_stats._load_run_stats")
    def test_markdown_to_file(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        (tmp_path / "run_1.h5").touch()
        out_file = tmp_path / "summary.md"

        result = runner.invoke(
            cli, ["macro-stats", str(tmp_path), "--markdown", str(out_file)]
        )

        assert result.exit_code == 0
        assert out_file.exists()
        content = out_file.read_text()
        assert "run_1" in content

    @patch("lysis.cli.macro_stats._load_run_stats")
    def test_no_markdown_gives_no_pipe_chars_as_headings(
        self, mock_load, runner, mock_stats, tmp_path
    ):
        """Without --markdown the output is a Rich table, not raw markdown."""
        mock_load.return_value = mock_stats
        (tmp_path / "run_1.h5").touch()

        result = runner.invoke(
            cli, ["macro-stats", str(tmp_path), "--no-progress"]
        )

        assert result.exit_code == 0
        # The full metric name should appear in the normal output too
        # (Rich table columns), but the first line should not be a markdown header
        assert not result.output.strip().startswith("#")
