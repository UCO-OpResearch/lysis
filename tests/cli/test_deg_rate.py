"""Tests for ``lysis deg-rate`` CLI command."""

from unittest.mock import patch

import pytest
from click.testing import CliRunner

from lysis.cli import cli

# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------

_DEFAULT_INTERVALS = [(20, 80), (20, 50), (50, 80)]


def _make_stats(mean=0.1234, std=0.0056):
    """Return a mock stats dict as returned by _load_run_deg_rates."""
    return {ivl: (mean, std) for ivl in _DEFAULT_INTERVALS}


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def mock_stats():
    return _make_stats()


# ---------------------------------------------------------------------------
# Wrappers around the new public API (replaces removed private helpers)
# ---------------------------------------------------------------------------


def _deg_rate_to_markdown(rows_dict, intervals, single_file_code=None):
    from lysis.analysis.summary import deg_rate_table
    from lysis.tools.display import stats_df_to_markdown
    md_col_headers = {f"{s}% to {e}%": f"{s}% to {e}% (%/min)" for s, e in intervals}
    df = deg_rate_table(rows_dict, intervals).rename(columns=md_col_headers)
    return stats_df_to_markdown(df, "Run", single_code=single_file_code)


# ---------------------------------------------------------------------------
# Help
# ---------------------------------------------------------------------------


class TestDegRateHelp:
    def test_help_exits_0(self, runner):
        result = runner.invoke(cli, ["deg-rate", "--help"])
        assert result.exit_code == 0

    def test_help_mentions_markdown(self, runner):
        result = runner.invoke(cli, ["deg-rate", "--help"])
        assert "--markdown" in result.output

    def test_help_mentions_no_progress(self, runner):
        result = runner.invoke(cli, ["deg-rate", "--help"])
        assert "--no-progress" in result.output

    def test_help_mentions_add_drop(self, runner):
        result = runner.invoke(cli, ["deg-rate", "--help"])
        assert "--add" in result.output
        assert "--drop" in result.output


# ---------------------------------------------------------------------------
# Interval helpers (unit tests)
# ---------------------------------------------------------------------------


class TestParseInterval:
    def test_valid_interval(self):
        from lysis.cli.deg_rate import _parse_interval

        assert _parse_interval("20-80") == (20, 80)

    def test_valid_interval_0_100(self):
        from lysis.cli.deg_rate import _parse_interval

        assert _parse_interval("0-100") == (0, 100)

    def test_invalid_format_raises(self):
        import click
        from lysis.cli.deg_rate import _parse_interval

        with pytest.raises(click.BadParameter):
            _parse_interval("20")

    def test_non_integer_raises(self):
        import click
        from lysis.cli.deg_rate import _parse_interval

        with pytest.raises(click.BadParameter):
            _parse_interval("a-b")

    def test_reversed_raises(self):
        import click
        from lysis.cli.deg_rate import _parse_interval

        with pytest.raises(click.BadParameter):
            _parse_interval("80-20")

    def test_equal_raises(self):
        import click
        from lysis.cli.deg_rate import _parse_interval

        with pytest.raises(click.BadParameter):
            _parse_interval("50-50")

    def test_out_of_range_raises(self):
        import click
        from lysis.cli.deg_rate import _parse_interval

        with pytest.raises(click.BadParameter):
            _parse_interval("0-101")


# ---------------------------------------------------------------------------
# Markdown helpers (unit tests)
# ---------------------------------------------------------------------------


class TestMdTableDegRate:
    def test_header_row(self):
        from lysis.tools.display import md_table as _md_table

        out = _md_table(["A", "B"], [["1", "2"]])
        assert out.splitlines()[0] == "| A | B |"

    def test_separator_row(self):
        from lysis.tools.display import md_table as _md_table

        out = _md_table(["A", "B"], [["1", "2"]])
        assert out.splitlines()[1] == "| --- | --- |"

    def test_data_row(self):
        from lysis.tools.display import md_table as _md_table

        out = _md_table(["A", "B"], [["hello", "world"]])
        assert "| hello | world |" in out


class TestDegRateToMarkdown:
    def test_single_file_heading(self, mock_stats):
        md = _deg_rate_to_markdown(
            {"run_A": mock_stats}, _DEFAULT_INTERVALS, single_file_code="run_A"
        )
        assert md.startswith("## run_A")

    def test_single_file_two_columns(self, mock_stats):
        # Single-file markdown now shows combined "mean ± std" in one Value column
        md = _deg_rate_to_markdown(
            {"run_A": mock_stats}, _DEFAULT_INTERVALS, single_file_code="run_A"
        )
        header_line = [l for l in md.splitlines() if "Metric" in l][0]
        assert header_line.count("|") == 3  # | Metric | Value |

    def test_single_file_interval_label(self, mock_stats):
        md = _deg_rate_to_markdown(
            {"run_A": mock_stats}, _DEFAULT_INTERVALS, single_file_code="run_A"
        )
        assert "20% to 80%" in md

    def test_single_file_values_formatted(self, mock_stats):
        md = _deg_rate_to_markdown(
            {"run_A": mock_stats}, _DEFAULT_INTERVALS, single_file_code="run_A"
        )
        assert "0.1234" in md
        assert "0.0056" in md

    def test_directory_run_codes_as_rows(self):
        rows = {"run_X": _make_stats(0.1, 0.01), "run_Y": _make_stats(0.2, 0.02)}
        md = _deg_rate_to_markdown(rows, _DEFAULT_INTERVALS)
        assert "run_X" in md
        assert "run_Y" in md

    def test_directory_header_starts_with_run(self):
        rows = {"run_X": _make_stats()}
        md = _deg_rate_to_markdown(rows, _DEFAULT_INTERVALS)
        assert md.splitlines()[0].startswith("| Run |")

    def test_directory_mean_pm_std_format(self):
        rows = {"run_X": _make_stats(0.1234, 0.0056)}
        md = _deg_rate_to_markdown(rows, _DEFAULT_INTERVALS)
        assert "0.1234" in md
        assert "\u00b1" in md  # ±

    def test_directory_all_intervals_present(self):
        rows = {"run_X": _make_stats()}
        md = _deg_rate_to_markdown(rows, _DEFAULT_INTERVALS)
        for s, e in _DEFAULT_INTERVALS:
            assert f"{s}% to {e}%" in md


# ---------------------------------------------------------------------------
# CLI integration — single-file mode
# ---------------------------------------------------------------------------


class TestDegRateMarkdownSingleFile:
    @patch("lysis.cli.deg_rate._load_run_deg_rates")
    def test_markdown_to_stdout(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(cli, ["deg-rate", str(h5), "--markdown", "-"])

        assert result.exit_code == 0
        assert "## run_A" in result.output

    @patch("lysis.cli.deg_rate._load_run_deg_rates")
    def test_markdown_to_stdout_has_table(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_B.h5"
        h5.touch()

        result = runner.invoke(cli, ["deg-rate", str(h5), "--markdown", "-"])

        assert result.exit_code == 0
        assert "| Metric | Value |" in result.output

    @patch("lysis.cli.deg_rate._load_run_deg_rates")
    def test_markdown_to_stdout_contains_interval(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_C.h5"
        h5.touch()

        result = runner.invoke(cli, ["deg-rate", str(h5), "--markdown", "-"])

        assert result.exit_code == 0
        assert "20% to 80%" in result.output

    @patch("lysis.cli.deg_rate._load_run_deg_rates")
    def test_markdown_to_file(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_D.h5"
        h5.touch()
        out_file = tmp_path / "rates.md"

        result = runner.invoke(cli, ["deg-rate", str(h5), "--markdown", str(out_file)])

        assert result.exit_code == 0
        assert out_file.exists()

    @patch("lysis.cli.deg_rate._load_run_deg_rates")
    def test_markdown_file_content(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_E.h5"
        h5.touch()
        out_file = tmp_path / "rates.md"

        runner.invoke(cli, ["deg-rate", str(h5), "--markdown", str(out_file)])

        content = out_file.read_text()
        assert "## run_E" in content
        assert "| Metric | Value |" in content

    @patch("lysis.cli.deg_rate._load_run_deg_rates")
    def test_markdown_file_confirmation_message(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_F.h5"
        h5.touch()
        out_file = tmp_path / "rates.md"

        result = runner.invoke(cli, ["deg-rate", str(h5), "--markdown", str(out_file)])

        assert "Markdown written to" in result.output

    @patch("lysis.cli.deg_rate._load_run_deg_rates")
    def test_no_markdown_flag_gives_normal_output(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_G.h5"
        h5.touch()

        result = runner.invoke(cli, ["deg-rate", str(h5), "--no-progress"])

        assert result.exit_code == 0
        assert "## run_G" not in result.output
        assert "run_G" in result.output


# ---------------------------------------------------------------------------
# CLI integration — directory mode
# ---------------------------------------------------------------------------


class TestDegRateMarkdownDirectory:
    @patch("lysis.cli.deg_rate._load_run_deg_rates")
    def test_markdown_to_stdout(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        (tmp_path / "run_1.h5").touch()
        (tmp_path / "run_2.h5").touch()

        result = runner.invoke(cli, ["deg-rate", str(tmp_path), "--markdown", "-"])

        assert result.exit_code == 0
        assert "|" in result.output

    @patch("lysis.cli.deg_rate._load_run_deg_rates")
    def test_markdown_to_stdout_run_codes_as_rows(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        (tmp_path / "run_A.h5").touch()
        (tmp_path / "run_B.h5").touch()

        result = runner.invoke(cli, ["deg-rate", str(tmp_path), "--markdown", "-"])

        assert result.exit_code == 0
        assert "run_A" in result.output
        assert "run_B" in result.output
        assert result.output.strip().startswith("| Run |")

    @patch("lysis.cli.deg_rate._load_run_deg_rates")
    def test_markdown_to_stdout_intervals_in_header(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        (tmp_path / "run_A.h5").touch()

        result = runner.invoke(cli, ["deg-rate", str(tmp_path), "--markdown", "-"])

        assert "20% to 80%" in result.output

    @patch("lysis.cli.deg_rate._load_run_deg_rates")
    def test_markdown_to_file(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        (tmp_path / "run_1.h5").touch()
        out_file = tmp_path / "rates.md"

        result = runner.invoke(cli, ["deg-rate", str(tmp_path), "--markdown", str(out_file)])

        assert result.exit_code == 0
        assert out_file.exists()
        content = out_file.read_text()
        assert "run_1" in content

    @patch("lysis.cli.deg_rate._load_run_deg_rates")
    def test_no_markdown_flag_gives_normal_output(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        (tmp_path / "run_1.h5").touch()

        result = runner.invoke(cli, ["deg-rate", str(tmp_path), "--no-progress"])

        assert result.exit_code == 0
        assert not result.output.strip().startswith("#")


# ---------------------------------------------------------------------------
# --add / --drop interval options
# ---------------------------------------------------------------------------


class TestDegRateAddDrop:
    @patch("lysis.cli.deg_rate._load_run_deg_rates")
    def test_add_interval(self, mock_load, runner, tmp_path):
        stats = {**_make_stats(), (0, 100): (0.05, 0.005)}
        mock_load.return_value = stats
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(
            cli, ["deg-rate", str(h5), "--add", "0-100", "--markdown", "-"]
        )

        assert result.exit_code == 0
        assert "0% to 100%" in result.output

    @patch("lysis.cli.deg_rate._load_run_deg_rates")
    def test_drop_interval(self, mock_load, runner, tmp_path):
        stats = {(20, 50): (0.1, 0.01), (50, 80): (0.1, 0.01)}
        mock_load.return_value = stats
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(
            cli, ["deg-rate", str(h5), "--drop", "20-80", "--markdown", "-"]
        )

        assert result.exit_code == 0
        assert "20% to 80%" not in result.output

    def test_drop_all_intervals_exits_nonzero(self, runner, tmp_path):
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(
            cli,
            [
                "deg-rate",
                str(h5),
                "--drop", "20-80",
                "--drop", "20-50",
                "--drop", "50-80",
            ],
        )

        assert result.exit_code != 0

    def test_invalid_add_interval_exits_nonzero(self, runner, tmp_path):
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(cli, ["deg-rate", str(h5), "--add", "invalid"])

        assert result.exit_code != 0

    def test_invalid_drop_interval_exits_nonzero(self, runner, tmp_path):
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(cli, ["deg-rate", str(h5), "--drop", "80-20"])

        assert result.exit_code != 0
