"""Tests for ``lysis compare`` CLI command (file-pair mode)."""

import os
import shutil

import h5py
import numpy as np
import pytest
from click.testing import CliRunner

from lysis.cli import cli
from lysis.config.parameters import MicroParameters
from lysis.dataio.dataspec import dataspec
from lysis.dataio.datastore import DataStore


@pytest.fixture
def runner():
    return CliRunner()


def _fill_microscale_data(filepath, n_sims=5):
    """Populate microscale datasets in an empty micro-only HDF5 file."""
    spec = dataspec["v2.0.0"]["microscale_out"]
    with h5py.File(filepath, "a") as f:
        for name, ds_spec in spec.data.items():
            if ds_spec.data_location is None:
                continue
            if ds_spec.dtype == h5py.string_dtype():
                continue
            path = ds_spec.data_location
            if path not in f:
                continue
            dataset = f[path]
            dataset.resize((n_sims,) + dataset.shape[1:])
            if ds_spec.dtype == np.bool_:
                dataset[:] = np.array([True, False, True, False, True])
            elif np.issubdtype(ds_spec.dtype, np.integer):
                dataset[:] = np.arange(n_sims, dtype=ds_spec.dtype)
            else:
                dataset[:] = np.arange(n_sims, dtype=ds_spec.dtype) * 1.5


@pytest.fixture
def two_identical_files(tmp_path):
    """Build two HDF5 files with identical microscale data."""
    micro = MicroParameters()
    dir1 = tmp_path / "a"
    dir2 = tmp_path / "b"
    dir1.mkdir()
    dir2.mkdir()

    ds = DataStore.create("run01", str(dir1), micro)
    ds.close()
    _fill_microscale_data(str(dir1 / "run01.h5"))

    shutil.copy(str(dir1 / "run01.h5"), str(dir2 / "run01.h5"))
    return dir1 / "run01.h5", dir2 / "run01.h5"


@pytest.fixture
def two_divergent_files(two_identical_files):
    """Start from two identical files, then mutate one dataset in the second."""
    _, path2 = two_identical_files
    with h5py.File(str(path2), "a") as f:
        arr = f["micro_data/tpa_leaving_time"][:]
        arr[2] = arr[2] + 1.0  # one element differs
        f["micro_data/tpa_leaving_time"][:] = arr
    return two_identical_files


class TestComparePathTypes:
    def test_mixed_file_and_dir_rejected(self, runner, tmp_path, two_identical_files):
        path1, _ = two_identical_files
        result = runner.invoke(
            cli,
            [
                "compare",
                "micro-data",
                str(path1),
                str(tmp_path),  # a directory
                "--no-progress",
            ],
        )
        assert result.exit_code != 0
        assert "both be directories or both be .h5 files" in result.output


class TestCompareFilePairIdentical:
    def test_micro_data_reports_ok_for_all_tables(self, runner, two_identical_files):
        path1, path2 = two_identical_files
        result = runner.invoke(
            cli,
            ["compare", "micro-data", str(path1), str(path2), "--no-progress"],
        )
        assert result.exit_code == 0, result.output
        # Every table should render as OK; no DIFF rows.
        assert "OK" in result.output
        assert "DIFF" not in result.output

    def test_detail_table_lists_tables_not_runs(self, runner, two_identical_files):
        path1, path2 = two_identical_files
        result = runner.invoke(
            cli,
            ["compare", "micro-data", str(path1), str(path2), "--no-progress"],
        )
        assert result.exit_code == 0, result.output
        # First column header is "Table", not "Run".
        assert "Table" in result.output


class TestCompareFilePairDivergent:
    def test_diff_table_shows_diff_row(self, runner, two_divergent_files):
        path1, path2 = two_divergent_files
        result = runner.invoke(
            cli,
            ["compare", "micro-data", str(path1), str(path2), "--no-progress"],
        )
        assert result.exit_code == 0, result.output
        assert "DIFF" in result.output
        assert "microscale_out/tpa_leaving_time" in result.output

    def test_location_is_rendered(self, runner, two_divergent_files):
        path1, path2 = two_divergent_files
        result = runner.invoke(
            cli,
            ["compare", "micro-data", str(path1), str(path2), "--no-progress"],
        )
        assert result.exit_code == 0, result.output
        assert "(2,)" in result.output


class TestCompareDiffFlag:
    def test_diff_rejects_folder_mode(self, runner, tmp_path):
        (tmp_path / "a").mkdir()
        (tmp_path / "b").mkdir()
        result = runner.invoke(
            cli,
            [
                "compare",
                "micro-data",
                str(tmp_path / "a"),
                str(tmp_path / "b"),
                "--diff",
                "microscale_out/tpa_leaving_time",
                "--no-progress",
            ],
        )
        assert result.exit_code != 0
        assert "file-pair mode" in result.output

    def test_diff_unknown_table_errors_with_available(
        self, runner, two_identical_files
    ):
        path1, path2 = two_identical_files
        result = runner.invoke(
            cli,
            [
                "compare",
                "data",
                str(path1),
                str(path2),
                "--diff",
                "microscale_out/does_not_exist",
                "--no-progress",
            ],
        )
        assert result.exit_code != 0
        assert "not found" in result.output
        assert "Available tables" in result.output
        assert "microscale_out/tpa_leaving_time" in result.output

    def test_diff_renders_side_by_side(self, runner, two_divergent_files):
        path1, path2 = two_divergent_files
        result = runner.invoke(
            cli,
            [
                "compare",
                "data",
                str(path1),
                str(path2),
                "--diff",
                "microscale_out/tpa_leaving_time",
                "--no-progress",
            ],
        )
        assert result.exit_code == 0, result.output
        # Title is the dataset label.
        assert "microscale_out/tpa_leaving_time" in result.output
        # Columns: Index, file names, % Diff.
        assert "Index" in result.output
        assert "% Diff" in result.output
        # Element [2] was mutated; some pct-diff should appear.
        assert "%" in result.output

    def test_diff_matches_leave_no_pct(self, runner, two_identical_files):
        path1, path2 = two_identical_files
        result = runner.invoke(
            cli,
            [
                "compare",
                "data",
                str(path1),
                str(path2),
                "--diff",
                "microscale_out/tpa_leaving_time",
                "--no-progress",
            ],
        )
        assert result.exit_code == 0, result.output
        # Matching rows render an em-dash in the % Diff column; no
        # signed percent value (e.g. "+1.23%") should appear anywhere.
        assert "—" in result.output
        import re

        assert re.search(r"[+-]\d+\.\d+%", result.output) is None


class TestCompareFilePairMarkdown:
    def test_markdown_output_contains_all_columns(
        self, runner, two_divergent_files, tmp_path
    ):
        path1, path2 = two_divergent_files
        out = tmp_path / "report.md"
        result = runner.invoke(
            cli,
            [
                "compare",
                "micro-data",
                str(path1),
                str(path2),
                "--markdown",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        text = out.read_text()
        for header in ("Table", "Status", "Mismatches", "Max % Diff", "Location"):
            assert header in text
        assert "DIFF" in text
