"""Tests for ``lysis validate`` CLI command."""

from unittest.mock import patch

import numpy as np
import pytest
from click.testing import CliRunner

from lysis.cli import cli
from lysis.dataio.dataspec import dataspec


@pytest.fixture
def runner():
    return CliRunner()


def _make_valid_hdf5_micro_data():
    """Build a minimal data dict that passes v2.0.0 microscale_out validation.

    Creates arrays with correct dtypes and plausible shapes for each dataset
    in the v2.0.0 microscale_out spec.
    """
    micro_spec = dataspec["v2.0.0"]["microscale_out"]
    data = {
        "params": {
            "micro_params": {"micro_simulations": 5},
            "macro_params": {
                "rows": 3,
                "cols": 2,
                "total_edges": 17,
                "total_molecules": 5,
                "empty_rows": 0,
                "full_row": 5,
                "pore_size": 1.0,
            },
        },
    }
    # Build a valid array for each dataset in the spec
    for name, spec in micro_spec.data.items():
        arr = np.zeros(5, dtype=spec.dtype)
        data[name] = arr
    return data


class TestValidateArgParsing:
    """Test argument parsing and validation."""

    def test_help_exits_0(self, runner):
        result = runner.invoke(cli, ["validate", "--help"])
        assert result.exit_code == 0
        assert "Validate simulation data" in result.output

    def test_missing_spec(self, runner):
        result = runner.invoke(cli, ["validate", "."])
        assert result.exit_code != 0

    def test_unknown_spec_prints_error(self, runner):
        result = runner.invoke(
            cli, ["validate", ".", "-s", "nonexistent"]
        )
        assert result.exit_code != 0
        assert "Unknown spec" in result.output

    def test_unknown_collection_prints_error(self, runner):
        with patch("lysis.dataio.fileops.read_data_collection"):
            result = runner.invoke(
                cli,
                ["validate", ".", "-s", "hdf5", "-c", "bogus_collection"],
            )
        assert result.exit_code != 0
        assert "Unknown collection" in result.output


class TestValidateExecution:
    """Test the validate command execution paths."""

    @patch("lysis.dataio.fileops.read_data_collection")
    def test_all_pass_returns_0(self, mock_read, runner):
        """All datasets passing returns exit code 0."""
        mock_read.return_value = _make_valid_hdf5_micro_data()

        result = runner.invoke(
            cli,
            ["validate", ".", "-s", "hdf5", "-c", "microscale_out"],
        )
        assert result.exit_code == 0
        assert "passed" in result.output
        # No FAIL or MISSING in output
        assert "FAIL" not in result.output
        assert "MISSING" not in result.output

    @patch("lysis.dataio.fileops.read_data_collection")
    def test_missing_dataset_returns_1(self, mock_read, runner):
        """Dataset present in spec but missing from data counts as FAIL."""
        data = _make_valid_hdf5_micro_data()
        # Remove one dataset
        first_dataset = next(
            k for k in data if k != "params"
        )
        del data[first_dataset]
        mock_read.return_value = data

        result = runner.invoke(
            cli,
            ["validate", ".", "-s", "hdf5", "-c", "microscale_out"],
        )
        assert result.exit_code == 1
        assert "MISSING" in result.output

    @patch("lysis.dataio.fileops.read_data_collection")
    def test_wrong_dtype_returns_1(self, mock_read, runner):
        """Dataset with wrong dtype counts as FAIL."""
        data = _make_valid_hdf5_micro_data()
        # Find a dataset with a numeric dtype (not object) and replace it
        micro_spec = dataspec["v2.0.0"]["microscale_out"]
        numeric_dataset = next(
            name for name, spec in micro_spec.data.items()
            if spec.dtype != object and name in data
        )
        data[numeric_dataset] = np.array(["a", "b"], dtype=object)
        mock_read.return_value = data

        result = runner.invoke(
            cli,
            ["validate", ".", "-s", "hdf5", "-c", "microscale_out"],
        )
        assert result.exit_code == 1
        assert "FAIL" in result.output

    @patch("lysis.dataio.fileops.read_data_collection")
    def test_output_shows_ok_labels(self, mock_read, runner):
        """Output contains OK labels for passing datasets."""
        mock_read.return_value = _make_valid_hdf5_micro_data()

        result = runner.invoke(
            cli,
            ["validate", ".", "-s", "hdf5", "-c", "microscale_out"],
        )
        assert "OK" in result.output

    @patch("lysis.dataio.fileops.read_data_collection")
    def test_output_shows_collection_header(self, mock_read, runner):
        """Output includes collection name as a header."""
        mock_read.return_value = _make_valid_hdf5_micro_data()

        result = runner.invoke(
            cli,
            ["validate", ".", "-s", "hdf5", "-c", "microscale_out"],
        )
        assert "microscale_out" in result.output

    @patch("lysis.dataio.fileops.read_data_collection")
    def test_per_sim_list_validated(self, mock_read, runner):
        """Per-simulation data (list of arrays) is validated correctly."""
        data = _make_valid_hdf5_micro_data()
        # Make one dataset a list (per-simulation)
        first_dataset = next(
            k for k in data if k != "params"
        )
        spec = dataspec["v2.0.0"]["microscale_out"].data[first_dataset]
        arr = np.zeros(5, dtype=spec.dtype)
        data[first_dataset] = [arr, arr, arr]
        mock_read.return_value = data

        result = runner.invoke(
            cli,
            ["validate", ".", "-s", "hdf5", "-c", "microscale_out"],
        )
        # Should show simulation count
        assert "3 sims" in result.output

    @patch("lysis.dataio.fileops.read_data_collection",
           side_effect=FileNotFoundError("output.h5"))
    def test_file_not_found_returns_1(self, mock_read, runner):
        """Missing data file returns exit code 1."""
        result = runner.invoke(
            cli,
            ["validate", ".", "-s", "hdf5"],
        )
        assert result.exit_code == 1
        assert "not found" in result.output

    @patch("lysis.dataio.fileops.read_data_collection")
    def test_result_summary_counts(self, mock_read, runner):
        """Result summary shows correct pass/fail counts."""
        data = _make_valid_hdf5_micro_data()
        # Remove one to create a failure
        first_dataset = next(
            k for k in data if k != "params"
        )
        del data[first_dataset]
        mock_read.return_value = data

        result = runner.invoke(
            cli,
            ["validate", ".", "-s", "hdf5", "-c", "microscale_out"],
        )
        assert "passed" in result.output
        assert "failed" in result.output
