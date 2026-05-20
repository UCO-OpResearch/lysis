"""Tests for ``lysis compare`` CLI command.

``lysis compare`` now covers only the two statistical subcommands
(``micro-stats`` / ``macro-stats``) — element-wise data-table diffing
was moved to the ``lysis diff`` command (see ``tests/cli/test_diff.py``).
These tests check path-type validation and command registration; the
statistical comparisons themselves are exercised at the analysis layer
(``tests/analysis/test_compare.py``).
"""

import pytest
from click.testing import CliRunner

from lysis.cli import cli


@pytest.fixture
def runner():
    return CliRunner()


class TestCompareRegistration:
    def test_command_is_registered(self, runner):
        result = runner.invoke(cli, ["compare", "--help"])
        assert result.exit_code == 0
        assert "compare" in result.output

    def test_only_stats_subcommands_accepted(self, runner):
        result = runner.invoke(cli, ["compare", "--help"])
        assert result.exit_code == 0
        # WHICH argument only offers the two stats keys now.
        assert "micro-stats" in result.output
        assert "macro-stats" in result.output
        # Data-table modes moved to lysis diff.
        assert "micro-data" not in result.output
        assert "[data|" not in result.output and "|data]" not in result.output

    def test_data_choice_rejected(self, runner, tmp_path):
        (tmp_path / "a").mkdir()
        (tmp_path / "b").mkdir()
        result = runner.invoke(
            cli,
            [
                "compare",
                "data",
                str(tmp_path / "a"),
                str(tmp_path / "b"),
                "--no-progress",
            ],
        )
        # Click rejects at parse time with exit_code=2.
        assert result.exit_code != 0
        assert "Invalid value" in result.output or "invalid choice" in result.output.lower()


class TestComparePathTypes:
    def test_mixed_file_and_dir_rejected(self, runner, tmp_path):
        # Create a minimal h5 file and a directory.
        h5_path = tmp_path / "a.h5"
        h5_path.write_bytes(b"")
        (tmp_path / "b").mkdir()
        result = runner.invoke(
            cli,
            [
                "compare",
                "micro-stats",
                str(h5_path),
                str(tmp_path / "b"),
                "--no-progress",
            ],
        )
        assert result.exit_code != 0
        assert "both be directories or both be .h5 files" in result.output
