"""Tests for ``lysis convert`` CLI command."""

from unittest.mock import patch, MagicMock

import numpy as np
import pytest
from click.testing import CliRunner

from lysis.cli import cli


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def mock_data():
    """Minimal data dict returned by read_data_collection."""
    return {
        "params": {"micro_params": {"micro_simulations": 2}},
        "lysis": np.array([100.0, 200.0]),
    }


class TestConvertArgParsing:
    """Test argument parsing and validation."""

    def test_help_exits_0(self, runner):
        result = runner.invoke(cli, ["convert", "--help"])
        assert result.exit_code == 0
        assert "Convert simulation data" in result.output

    def test_missing_required_args(self, runner):
        result = runner.invoke(cli, ["convert"])
        assert result.exit_code != 0
        assert "Missing" in result.output or "Usage" in result.output

    def test_missing_from_flag(self, runner):
        result = runner.invoke(cli, ["convert", "/in", "/out", "-t", "hdf5"])
        assert result.exit_code != 0

    def test_missing_to_flag(self, runner):
        result = runner.invoke(cli, ["convert", "/in", "/out", "-f", "fortran"])
        assert result.exit_code != 0


class TestResolveSpec:
    """Test spec resolution and error messages."""

    def test_unknown_spec_prints_error(self, runner):
        """Unknown spec version produces a helpful error."""
        result = runner.invoke(
            cli,
            ["convert", ".", ".", "-f", "nonexistent", "-t", "hdf5"],
        )
        assert result.exit_code != 0
        assert "Unknown spec" in result.output

    def test_unknown_spec_lists_available(self, runner):
        """Error message lists available versions and tags."""
        result = runner.invoke(
            cli,
            ["convert", ".", ".", "-f", "nonexistent", "-t", "hdf5"],
        )
        assert "v1.99.0" in result.output or "v2.0.0" in result.output

    def test_unknown_collection_prints_error(self, runner):
        """Unknown collection name produces a helpful error."""
        with patch("lysis.data.fileops.read_data_collection"):
            result = runner.invoke(
                cli,
                ["convert", ".", ".", "-f", "fortran", "-t", "hdf5",
                 "-c", "nonexistent_collection"],
            )
        assert result.exit_code != 0
        assert "Unknown collection" in result.output


class TestConvertExecution:
    """Test the convert command execution paths."""

    @patch("lysis.data.fileops.write_data_collection")
    @patch("lysis.data.dataconvert.convert_data")
    @patch("lysis.data.fileops.read_data_collection")
    def test_successful_conversion(self, mock_read, mock_convert, mock_write,
                                   runner, mock_data):
        """Successful conversion reads, converts, and writes."""
        mock_read.return_value = mock_data
        mock_convert.return_value = mock_data

        result = runner.invoke(
            cli,
            ["convert", ".", ".", "-f", "fortran", "-t", "hdf5"],
        )
        assert result.exit_code == 0
        mock_read.assert_called_once()
        mock_convert.assert_called_once()
        mock_write.assert_called_once()

    @patch("lysis.data.fileops.write_data_collection")
    @patch("lysis.data.dataconvert.convert_data")
    @patch("lysis.data.fileops.read_data_collection")
    def test_dry_run_skips_write(self, mock_read, mock_convert, mock_write,
                                 runner, mock_data):
        """--dry-run reads and converts but does not write."""
        mock_read.return_value = mock_data
        mock_convert.return_value = mock_data

        result = runner.invoke(
            cli,
            ["convert", ".", ".", "-f", "fortran", "-t", "hdf5", "--dry-run"],
        )
        assert result.exit_code == 0
        mock_read.assert_called_once()
        mock_convert.assert_called_once()
        mock_write.assert_not_called()
        assert "Dry run" in result.output

    @patch("lysis.data.fileops.read_data_collection",
           side_effect=FileNotFoundError("sim_00.dat"))
    def test_file_not_found_returns_1(self, mock_read, runner):
        """Missing input file returns exit code 1."""
        result = runner.invoke(
            cli,
            ["convert", ".", ".", "-f", "fortran", "-t", "hdf5"],
        )
        assert result.exit_code == 1
        assert "not found" in result.output

    @patch("lysis.data.fileops.read_data_collection")
    @patch("lysis.data.dataconvert.convert_data",
           side_effect=NotImplementedError("converter stub"))
    def test_not_implemented_returns_1(self, mock_convert, mock_read,
                                       runner, mock_data):
        """Unimplemented converter returns exit code 1."""
        mock_read.return_value = mock_data

        result = runner.invoke(
            cli,
            ["convert", ".", ".", "-f", "fortran", "-t", "hdf5"],
        )
        assert result.exit_code == 1
        assert "not implemented" in result.output.lower()

    @patch("lysis.data.fileops.read_data_collection")
    @patch("lysis.data.dataconvert.convert_data",
           side_effect=OverflowError("int32 overflow"))
    def test_overflow_returns_1(self, mock_convert, mock_read,
                                runner, mock_data):
        """Type overflow during conversion returns exit code 1."""
        mock_read.return_value = mock_data

        result = runner.invoke(
            cli,
            ["convert", ".", ".", "-f", "fortran", "-t", "hdf5"],
        )
        assert result.exit_code == 1
        assert "conversion failed" in result.output.lower()

    @patch("lysis.data.fileops.write_data_collection")
    @patch("lysis.data.dataconvert.convert_data")
    @patch("lysis.data.fileops.read_data_collection")
    def test_verbose_shows_details(self, mock_read, mock_convert, mock_write,
                                   runner, mock_data):
        """Verbose flag shows conversion details."""
        mock_read.return_value = mock_data
        mock_convert.return_value = mock_data

        result = runner.invoke(
            cli,
            ["-v", "convert", ".", ".", "-f", "fortran", "-t", "hdf5"],
        )
        assert result.exit_code == 0
        assert "v1.99.0" in result.output
        assert "v2.0.0" in result.output

    @patch("lysis.data.fileops.write_data_collection")
    @patch("lysis.data.dataconvert.convert_data")
    @patch("lysis.data.fileops.read_data_collection")
    def test_collections_filter(self, mock_read, mock_convert, mock_write,
                                runner, mock_data):
        """--collections filters which collections are processed."""
        mock_read.return_value = mock_data
        mock_convert.return_value = mock_data

        result = runner.invoke(
            cli,
            ["convert", ".", ".", "-f", "fortran", "-t", "hdf5",
             "-c", "microscale_out"],
        )
        assert result.exit_code == 0
        # read_data_collection should receive only the microscale_out collection
        call_args = mock_read.call_args
        collections_arg = call_args[0][1]  # second positional arg
        assert len(collections_arg) == 1
