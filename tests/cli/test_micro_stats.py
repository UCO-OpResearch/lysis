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
# Wrappers around the new public API (replaces removed private helpers)
# ---------------------------------------------------------------------------


def _micro_stats_to_markdown(rows_dict, single_file_code=None):
    from lysis.analysis.summary import micro_stats_table
    from lysis.cli.display import stats_df_to_markdown
    df = micro_stats_table(rows_dict)
    return stats_df_to_markdown(df, "Run", single_code=single_file_code)


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
# Formatting helpers (unit tests)
# ---------------------------------------------------------------------------


class TestFmtMetric:
    def test_fibers_degraded_formatted_with_comma(self):
        from lysis.analysis.summary import _fmt_micro as _fmt_metric

        stats = _make_stats(fibers_degraded=1234)
        assert _fmt_metric("Fibers Degraded", stats) == "1,234"

    def test_mean_lysis_time_shows_pm(self):
        from lysis.analysis.summary import _fmt_micro as _fmt_metric

        stats = _make_stats(lysis_time_mean=42.0, lysis_time_std=3.5)
        result = _fmt_metric("Mean Lysis Time (min)", stats)
        assert "\u00b1" in result  # ±
        assert "42.000" in result
        assert "3.500" in result

    def test_median_lysis_time_no_pm(self):
        from lysis.analysis.summary import _fmt_micro as _fmt_metric

        stats = _make_stats(lysis_time_median=41.5)
        result = _fmt_metric("Median Lysis Time (min)", stats)
        assert "\u00b1" not in result
        assert "41.500" in result

    def test_mean_tpa_leaving_shows_pm(self):
        from lysis.analysis.summary import _fmt_micro as _fmt_metric

        stats = _make_stats(tpa_leaving_mean=15.0, tpa_leaving_std=2.0)
        result = _fmt_metric("Mean tPA Leaving Time (sec)", stats)
        assert "\u00b1" in result
        assert "15.000" in result
        assert "2.000" in result

    def test_median_tpa_leaving_no_pm(self):
        from lysis.analysis.summary import _fmt_micro as _fmt_metric

        stats = _make_stats(tpa_leaving_median=14.5)
        result = _fmt_metric("Median tPA Leaving Time (sec)", stats)
        assert "\u00b1" not in result
        assert "14.500" in result

    def test_unknown_key_raises(self):
        from lysis.analysis.summary import _fmt_micro as _fmt_metric

        with pytest.raises(KeyError):
            _fmt_metric("Unknown Metric", _make_stats())


# ---------------------------------------------------------------------------
# Markdown helpers (unit tests)
# ---------------------------------------------------------------------------


class TestMdTableMicroStats:
    def test_header_row(self):
        from lysis.cli.display import md_table as _md_table

        out = _md_table(["A", "B"], [["1", "2"]])
        assert out.splitlines()[0] == "| A | B |"

    def test_separator_row(self):
        from lysis.cli.display import md_table as _md_table

        out = _md_table(["A", "B"], [["1", "2"]])
        assert out.splitlines()[1] == "| --- | --- |"

    def test_data_row(self):
        from lysis.cli.display import md_table as _md_table

        out = _md_table(["A", "B"], [["hello", "world"]])
        assert "| hello | world |" in out


class TestMicroStatsToMarkdown:
    def test_single_file_heading(self, mock_stats):
        md = _micro_stats_to_markdown({"run_A": mock_stats}, single_file_code="run_A")
        assert md.startswith("## run_A")

    def test_single_file_two_columns(self, mock_stats):
        md = _micro_stats_to_markdown({"run_A": mock_stats}, single_file_code="run_A")
        header_line = [l for l in md.splitlines() if "Metric" in l][0]
        # | Metric | Value |
        assert header_line.count("|") == 3

    def test_single_file_all_metrics_present(self, mock_stats):
        from lysis.analysis.summary import MICRO_STATS_COLUMNS

        md = _micro_stats_to_markdown({"run_A": mock_stats}, single_file_code="run_A")
        for m in MICRO_STATS_COLUMNS:
            assert m in md

    def test_single_file_fibers_degraded_value(self, mock_stats):
        md = _micro_stats_to_markdown({"run_A": mock_stats}, single_file_code="run_A")
        assert "1,234" in md

    def test_single_file_mean_lysis_time_formatted(self, mock_stats):
        md = _micro_stats_to_markdown({"run_A": mock_stats}, single_file_code="run_A")
        assert "42.000" in md
        assert "\u00b1" in md

    def test_directory_run_codes_as_rows(self):
        rows = {"run_X": _make_stats(), "run_Y": _make_stats(fibers_degraded=500)}
        md = _micro_stats_to_markdown(rows)
        assert "run_X" in md
        assert "run_Y" in md

    def test_directory_header_starts_with_run(self):
        rows = {"run_X": _make_stats()}
        md = _micro_stats_to_markdown(rows)
        assert md.splitlines()[0].startswith("| Run |")

    def test_directory_all_metrics_in_header(self):
        from lysis.analysis.summary import MICRO_STATS_COLUMNS

        rows = {"run_X": _make_stats()}
        md = _micro_stats_to_markdown(rows)
        header = md.splitlines()[0]
        for m in MICRO_STATS_COLUMNS:
            assert m in header


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
        assert "run_A" in result.output
        assert "run_B" in result.output
        assert result.output.strip().startswith("| Run |")

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
