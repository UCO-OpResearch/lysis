"""Tests for ``lysis deg-time`` CLI command."""

from unittest.mock import patch

import pytest
from click.testing import CliRunner

from lysis.cli import cli

# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------

_DEFAULT_MARKERS = [5, 20, 50, 80, 100]


def _make_stats(mean=42.0, std=3.5):
    """Return a mock stats dict as returned by _load_run_deg_times."""
    return {m: (mean, std) for m in _DEFAULT_MARKERS}


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def mock_stats():
    return _make_stats()


# ---------------------------------------------------------------------------
# Wrappers around the new public API (replaces removed private helpers)
# ---------------------------------------------------------------------------


def _deg_time_to_markdown(rows_dict, markers, single_file_code=None):
    from lysis.analysis.summary import deg_time_table
    from lysis.cli.display import stats_df_to_markdown
    md_col_headers = {f"{m}%": f"{m}% (min)" for m in markers}
    df = deg_time_table(rows_dict, markers).rename(columns=md_col_headers)
    return stats_df_to_markdown(df, "Run", single_code=single_file_code)


# ---------------------------------------------------------------------------
# Help
# ---------------------------------------------------------------------------


class TestDegTimeHelp:
    def test_help_exits_0(self, runner):
        result = runner.invoke(cli, ["deg-time", "--help"])
        assert result.exit_code == 0

    def test_help_mentions_markdown(self, runner):
        result = runner.invoke(cli, ["deg-time", "--help"])
        assert "--markdown" in result.output

    def test_help_mentions_no_progress(self, runner):
        result = runner.invoke(cli, ["deg-time", "--help"])
        assert "--no-progress" in result.output

    def test_help_mentions_add_drop(self, runner):
        result = runner.invoke(cli, ["deg-time", "--help"])
        assert "--add" in result.output
        assert "--drop" in result.output


# ---------------------------------------------------------------------------
# Marker helpers (unit tests)
# ---------------------------------------------------------------------------


class TestParseMarker:
    def test_plain_integer(self):
        from lysis.cli.deg_time import _parse_marker

        assert _parse_marker("50") == 50

    def test_percent_suffix(self):
        from lysis.cli.deg_time import _parse_marker

        assert _parse_marker("80%") == 80

    def test_zero(self):
        from lysis.cli.deg_time import _parse_marker

        assert _parse_marker("0") == 0

    def test_hundred(self):
        from lysis.cli.deg_time import _parse_marker

        assert _parse_marker("100") == 100

    def test_non_integer_raises(self):
        import click
        from lysis.cli.deg_time import _parse_marker

        with pytest.raises(click.BadParameter):
            _parse_marker("abc")

    def test_negative_raises(self):
        import click
        from lysis.cli.deg_time import _parse_marker

        with pytest.raises(click.BadParameter):
            _parse_marker("-10")

    def test_over_100_raises(self):
        import click
        from lysis.cli.deg_time import _parse_marker

        with pytest.raises(click.BadParameter):
            _parse_marker("101")


# ---------------------------------------------------------------------------
# Markdown helpers (unit tests)
# ---------------------------------------------------------------------------


class TestMdTableDegTime:
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


class TestDegTimeToMarkdown:
    def test_single_file_heading(self, mock_stats):
        md = _deg_time_to_markdown(
            {"run_A": mock_stats}, _DEFAULT_MARKERS, single_file_code="run_A"
        )
        assert md.startswith("## run_A")

    def test_single_file_two_columns(self, mock_stats):
        # Single-file markdown now shows combined "mean ± std" in one Value column
        md = _deg_time_to_markdown(
            {"run_A": mock_stats}, _DEFAULT_MARKERS, single_file_code="run_A"
        )
        header_line = [l for l in md.splitlines() if "Metric" in l][0]
        assert header_line.count("|") == 3  # | Metric | Value |

    def test_single_file_milestone_label(self, mock_stats):
        md = _deg_time_to_markdown(
            {"run_A": mock_stats}, _DEFAULT_MARKERS, single_file_code="run_A"
        )
        assert "50%" in md

    def test_single_file_values_formatted(self, mock_stats):
        md = _deg_time_to_markdown(
            {"run_A": mock_stats}, _DEFAULT_MARKERS, single_file_code="run_A"
        )
        assert "42.00" in md
        assert "3.50" in md

    def test_directory_run_codes_as_rows(self):
        rows = {"run_X": _make_stats(40.0, 2.0), "run_Y": _make_stats(50.0, 4.0)}
        md = _deg_time_to_markdown(rows, _DEFAULT_MARKERS)
        assert "run_X" in md
        assert "run_Y" in md

    def test_directory_header_starts_with_run(self):
        rows = {"run_X": _make_stats()}
        md = _deg_time_to_markdown(rows, _DEFAULT_MARKERS)
        assert md.splitlines()[0].startswith("| Run |")

    def test_directory_mean_pm_std_format(self):
        rows = {"run_X": _make_stats(42.0, 3.5)}
        md = _deg_time_to_markdown(rows, _DEFAULT_MARKERS)
        assert "42.00" in md
        assert "\u00b1" in md  # ±

    def test_directory_all_markers_present(self):
        rows = {"run_X": _make_stats()}
        md = _deg_time_to_markdown(rows, _DEFAULT_MARKERS)
        for m in _DEFAULT_MARKERS:
            assert f"{m}%" in md


# ---------------------------------------------------------------------------
# CLI integration — single-file mode
# ---------------------------------------------------------------------------


class TestDegTimeMarkdownSingleFile:
    @patch("lysis.cli.deg_time._load_run_deg_times")
    def test_markdown_to_stdout(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(cli, ["deg-time", str(h5), "--markdown", "-"])

        assert result.exit_code == 0
        assert "## run_A" in result.output

    @patch("lysis.cli.deg_time._load_run_deg_times")
    def test_markdown_to_stdout_has_table(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_B.h5"
        h5.touch()

        result = runner.invoke(cli, ["deg-time", str(h5), "--markdown", "-"])

        assert result.exit_code == 0
        assert "| Metric | Value |" in result.output

    @patch("lysis.cli.deg_time._load_run_deg_times")
    def test_markdown_to_stdout_contains_milestone(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_C.h5"
        h5.touch()

        result = runner.invoke(cli, ["deg-time", str(h5), "--markdown", "-"])

        assert result.exit_code == 0
        assert "50%" in result.output

    @patch("lysis.cli.deg_time._load_run_deg_times")
    def test_markdown_to_file(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_D.h5"
        h5.touch()
        out_file = tmp_path / "times.md"

        result = runner.invoke(cli, ["deg-time", str(h5), "--markdown", str(out_file)])

        assert result.exit_code == 0
        assert out_file.exists()

    @patch("lysis.cli.deg_time._load_run_deg_times")
    def test_markdown_file_content(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_E.h5"
        h5.touch()
        out_file = tmp_path / "times.md"

        runner.invoke(cli, ["deg-time", str(h5), "--markdown", str(out_file)])

        content = out_file.read_text()
        assert "## run_E" in content
        assert "| Metric | Value |" in content

    @patch("lysis.cli.deg_time._load_run_deg_times")
    def test_markdown_file_confirmation_message(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_F.h5"
        h5.touch()
        out_file = tmp_path / "times.md"

        result = runner.invoke(cli, ["deg-time", str(h5), "--markdown", str(out_file)])

        assert "Markdown written to" in result.output

    @patch("lysis.cli.deg_time._load_run_deg_times")
    def test_no_markdown_flag_gives_normal_output(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        h5 = tmp_path / "run_G.h5"
        h5.touch()

        result = runner.invoke(cli, ["deg-time", str(h5), "--no-progress"])

        assert result.exit_code == 0
        assert "## run_G" not in result.output
        assert "run_G" in result.output


# ---------------------------------------------------------------------------
# CLI integration — directory mode
# ---------------------------------------------------------------------------


class TestDegTimeMarkdownDirectory:
    @patch("lysis.cli.deg_time._load_run_deg_times")
    def test_markdown_to_stdout(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        (tmp_path / "run_1.h5").touch()
        (tmp_path / "run_2.h5").touch()

        result = runner.invoke(cli, ["deg-time", str(tmp_path), "--markdown", "-"])

        assert result.exit_code == 0
        assert "|" in result.output

    @patch("lysis.cli.deg_time._load_run_deg_times")
    def test_markdown_to_stdout_run_codes_as_rows(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        (tmp_path / "run_A.h5").touch()
        (tmp_path / "run_B.h5").touch()

        result = runner.invoke(cli, ["deg-time", str(tmp_path), "--markdown", "-"])

        assert result.exit_code == 0
        assert "run_A" in result.output
        assert "run_B" in result.output
        assert result.output.strip().startswith("| Run |")

    @patch("lysis.cli.deg_time._load_run_deg_times")
    def test_markdown_to_stdout_milestones_in_header(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        (tmp_path / "run_A.h5").touch()

        result = runner.invoke(cli, ["deg-time", str(tmp_path), "--markdown", "-"])

        assert "50%" in result.output

    @patch("lysis.cli.deg_time._load_run_deg_times")
    def test_markdown_to_file(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        (tmp_path / "run_1.h5").touch()
        out_file = tmp_path / "times.md"

        result = runner.invoke(cli, ["deg-time", str(tmp_path), "--markdown", str(out_file)])

        assert result.exit_code == 0
        assert out_file.exists()
        content = out_file.read_text()
        assert "run_1" in content

    @patch("lysis.cli.deg_time._load_run_deg_times")
    def test_no_markdown_flag_gives_normal_output(self, mock_load, runner, mock_stats, tmp_path):
        mock_load.return_value = mock_stats
        (tmp_path / "run_1.h5").touch()

        result = runner.invoke(cli, ["deg-time", str(tmp_path), "--no-progress"])

        assert result.exit_code == 0
        assert not result.output.strip().startswith("#")


# ---------------------------------------------------------------------------
# --add / --drop marker options
# ---------------------------------------------------------------------------


class TestDegTimeAddDrop:
    @patch("lysis.cli.deg_time._load_run_deg_times")
    def test_add_marker(self, mock_load, runner, tmp_path):
        stats = {**_make_stats(), 10: (35.0, 2.0)}
        mock_load.return_value = stats
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(
            cli, ["deg-time", str(h5), "--add", "10", "--markdown", "-"]
        )

        assert result.exit_code == 0
        assert "10%" in result.output

    @patch("lysis.cli.deg_time._load_run_deg_times")
    def test_add_marker_with_percent_suffix(self, mock_load, runner, tmp_path):
        stats = {**_make_stats(), 90: (60.0, 5.0)}
        mock_load.return_value = stats
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(
            cli, ["deg-time", str(h5), "--add", "90%", "--markdown", "-"]
        )

        assert result.exit_code == 0
        assert "90%" in result.output

    @patch("lysis.cli.deg_time._load_run_deg_times")
    def test_drop_marker(self, mock_load, runner, tmp_path):
        # Drop 5% so only 20, 50, 80, 100 remain
        stats = {20: (45.0, 2.0), 50: (55.0, 3.0), 80: (65.0, 4.0), 100: (75.0, 5.0)}
        mock_load.return_value = stats
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(
            cli, ["deg-time", str(h5), "--drop", "5", "--markdown", "-"]
        )

        assert result.exit_code == 0
        # 5% should not appear as a standalone milestone label
        lines = result.output.splitlines()
        milestone_cells = [l for l in lines if "| 5% |" in l or l.strip() == "| 5% |"]
        assert not milestone_cells

    def test_drop_all_markers_exits_nonzero(self, runner, tmp_path):
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(
            cli,
            [
                "deg-time", str(h5),
                "--drop", "5",
                "--drop", "20",
                "--drop", "50",
                "--drop", "80",
                "--drop", "100",
            ],
        )

        assert result.exit_code != 0

    def test_invalid_add_marker_exits_nonzero(self, runner, tmp_path):
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(cli, ["deg-time", str(h5), "--add", "invalid"])

        assert result.exit_code != 0

    def test_add_marker_over_100_exits_nonzero(self, runner, tmp_path):
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(cli, ["deg-time", str(h5), "--add", "110"])

        assert result.exit_code != 0

    @patch("lysis.cli.deg_time._load_run_deg_times")
    def test_markers_always_sorted(self, mock_load, runner, tmp_path):
        # Add 10 which should appear before the defaults after sorting.
        # Single-file markdown uses rows (not columns) per milestone,
        # so check that the "10%" row appears before the "20%" row.
        stats = {**_make_stats(), 10: (35.0, 2.0)}
        mock_load.return_value = stats
        h5 = tmp_path / "run_A.h5"
        h5.touch()

        result = runner.invoke(
            cli, ["deg-time", str(h5), "--add", "10", "--markdown", "-"]
        )

        assert result.exit_code == 0
        lines = result.output.splitlines()
        row_10 = next((i for i, l in enumerate(lines) if "| 10%" in l), None)
        row_20 = next((i for i, l in enumerate(lines) if "| 20%" in l), None)
        assert row_10 is not None and row_20 is not None
        assert row_10 < row_20
