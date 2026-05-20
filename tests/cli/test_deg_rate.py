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

        # ``result.stdout`` (not ``result.output``) excludes the top-level
        # src/lysis/ dirty warning that is emitted to stderr.
        assert result.exit_code == 0
        assert "run_A" in result.stdout
        assert "run_B" in result.stdout
        assert result.stdout.strip().startswith("| Run |")

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
