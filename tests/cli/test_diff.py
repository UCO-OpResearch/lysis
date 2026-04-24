"""Tests for the ``lysis diff`` CLI command.

File-pair mode fixtures are ported from the old ``lysis compare
data`` / ``micro-data`` / ``--diff`` tests; folder-pair and
scale-mismatch coverage is new.
"""

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


# ---------------------------------------------------------------------------
# Fixtures — micro-only file-pair builders
# ---------------------------------------------------------------------------


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


@pytest.fixture
def ulp_divergent_files(two_identical_files):
    """Start from two identical files, then shift one float dataset by 2 ULPs."""
    _, path2 = two_identical_files
    with h5py.File(str(path2), "a") as f:
        arr = f["micro_data/tpa_leaving_time"][:]
        shifted = np.nextafter(np.nextafter(arr, np.inf), np.inf)
        f["micro_data/tpa_leaving_time"][:] = shifted
    return two_identical_files


@pytest.fixture
def high_precision_files(two_identical_files):
    """Mutate with a value whose full-precision float64 string differs from :.6g."""
    _, path2 = two_identical_files
    with h5py.File(str(path2), "a") as f:
        arr = f["micro_data/tpa_leaving_time"][:]
        arr[0] = 0.1 + 0.2
        f["micro_data/tpa_leaving_time"][:] = arr
    return two_identical_files


# ---------------------------------------------------------------------------
# Fixtures — folder-pair builders
# ---------------------------------------------------------------------------


def _build_run_in_folder(folder, run_code):
    """Create a micro-only Run inside *folder* and populate it."""
    micro = MicroParameters()
    ds = DataStore.create(run_code, str(folder), micro)
    ds.close()
    _fill_microscale_data(str(folder / f"{run_code}.h5"))


@pytest.fixture
def two_folders_identical(tmp_path):
    """Two folders with the same three micro-only Runs, all identical."""
    dir1 = tmp_path / "foldA"
    dir2 = tmp_path / "foldB"
    dir1.mkdir()
    dir2.mkdir()
    for rc in ("run01", "run02", "run03"):
        _build_run_in_folder(dir1, rc)
        shutil.copy(str(dir1 / f"{rc}.h5"), str(dir2 / f"{rc}.h5"))
    return dir1, dir2


@pytest.fixture
def two_folders_divergent(two_folders_identical):
    """Perturb one run in the second folder so the diff has something to report."""
    _, dir2 = two_folders_identical
    path = dir2 / "run02.h5"
    with h5py.File(str(path), "a") as f:
        arr = f["micro_data/tpa_leaving_time"][:]
        arr[1] = arr[1] + 1.0
        f["micro_data/tpa_leaving_time"][:] = arr
    return two_folders_identical


# ---------------------------------------------------------------------------
# Path-type validation
# ---------------------------------------------------------------------------


class TestDiffPathTypes:
    def test_mixed_file_and_dir_rejected(
        self, runner, tmp_path, two_identical_files
    ):
        path1, _ = two_identical_files
        result = runner.invoke(
            cli,
            ["diff", str(path1), str(tmp_path), "--no-progress"],
        )
        assert result.exit_code != 0
        assert "both be directories or both be .h5 files" in result.output

    def test_file_pair_rejects_non_h5_extension(self, runner, tmp_path):
        bad1 = tmp_path / "a.dat"
        bad2 = tmp_path / "b.dat"
        bad1.write_bytes(b"")
        bad2.write_bytes(b"")
        result = runner.invoke(
            cli,
            ["diff", str(bad1), str(bad2), "--no-progress"],
        )
        assert result.exit_code != 0
        assert ".h5 extensions" in result.output

    def test_table_flag_requires_file_pair(self, runner, tmp_path):
        (tmp_path / "a").mkdir()
        (tmp_path / "b").mkdir()
        result = runner.invoke(
            cli,
            [
                "diff",
                str(tmp_path / "a"),
                str(tmp_path / "b"),
                "--table",
                "microscale_out/tpa_leaving_time",
                "--no-progress",
            ],
        )
        assert result.exit_code != 0
        assert "file-pair mode" in result.output


# ---------------------------------------------------------------------------
# File-pair mode
# ---------------------------------------------------------------------------


class TestDiffFilePairIdentical:
    def test_reports_ok_for_all_tables(self, runner, two_identical_files):
        path1, path2 = two_identical_files
        result = runner.invoke(
            cli, ["diff", str(path1), str(path2), "--no-progress"]
        )
        assert result.exit_code == 0, result.output
        assert "OK" in result.output
        assert "DIFF" not in result.output

    def test_detail_table_lists_tables_not_runs(
        self, runner, two_identical_files
    ):
        path1, path2 = two_identical_files
        result = runner.invoke(
            cli, ["diff", str(path1), str(path2), "--no-progress"]
        )
        assert result.exit_code == 0, result.output
        assert "Table" in result.output


class TestDiffFilePairUlpTolerance:
    def test_detail_table_reports_ok(self, runner, ulp_divergent_files):
        path1, path2 = ulp_divergent_files
        result = runner.invoke(
            cli, ["diff", str(path1), str(path2), "--no-progress"]
        )
        assert result.exit_code == 0, result.output
        assert "DIFF" not in result.output
        assert "OK" in result.output

    def test_table_side_by_side_shows_no_differences(
        self, runner, ulp_divergent_files
    ):
        path1, path2 = ulp_divergent_files
        result = runner.invoke(
            cli,
            [
                "diff",
                str(path1),
                str(path2),
                "--table",
                "microscale_out/tpa_leaving_time",
                "--no-progress",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "No differences found" in result.output


class TestDiffFilePairDivergent:
    def test_detail_table_shows_diff_row(self, runner, two_divergent_files):
        path1, path2 = two_divergent_files
        result = runner.invoke(
            cli, ["diff", str(path1), str(path2), "--no-progress"]
        )
        assert result.exit_code == 0, result.output
        assert "DIFF" in result.output
        assert "microscale_out/tpa_leaving_time" in result.output

    def test_location_is_rendered(self, runner, two_divergent_files):
        path1, path2 = two_divergent_files
        result = runner.invoke(
            cli, ["diff", str(path1), str(path2), "--no-progress"]
        )
        assert result.exit_code == 0, result.output
        assert "(2,)" in result.output


class TestDiffTableFlag:
    def test_unknown_table_errors_with_available(
        self, runner, two_identical_files
    ):
        path1, path2 = two_identical_files
        result = runner.invoke(
            cli,
            [
                "diff",
                str(path1),
                str(path2),
                "--table",
                "microscale_out/does_not_exist",
                "--no-progress",
            ],
        )
        assert result.exit_code != 0
        assert "not found" in result.output
        assert "Available tables" in result.output
        assert "microscale_out/tpa_leaving_time" in result.output

    def test_renders_side_by_side(self, runner, two_divergent_files):
        path1, path2 = two_divergent_files
        result = runner.invoke(
            cli,
            [
                "diff",
                str(path1),
                str(path2),
                "--table",
                "microscale_out/tpa_leaving_time",
                "--no-progress",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "microscale_out/tpa_leaving_time" in result.output
        assert "Index" in result.output
        assert "% Diff" in result.output
        assert "%" in result.output

    def test_identical_reports_no_differences(
        self, runner, two_identical_files
    ):
        path1, path2 = two_identical_files
        result = runner.invoke(
            cli,
            [
                "diff",
                str(path1),
                str(path2),
                "--table",
                "microscale_out/tpa_leaving_time",
                "--no-progress",
            ],
        )
        assert result.exit_code == 0, result.output
        assert (
            "No differences found in microscale_out/tpa_leaving_time"
            in result.output
        )
        import re

        assert re.search(r"[+-]\d+\.\d+%", result.output) is None

    def test_renders_full_precision(self, runner, high_precision_files):
        path1, path2 = high_precision_files
        result = runner.invoke(
            cli,
            [
                "diff",
                str(path1),
                str(path2),
                "--table",
                "microscale_out/tpa_leaving_time",
                "--no-progress",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "0.30000000000000004" in result.output

    def test_skips_matching_rows(self, runner, two_divergent_files):
        path1, path2 = two_divergent_files
        result = runner.invoke(
            cli,
            [
                "diff",
                str(path1),
                str(path2),
                "--table",
                "microscale_out/tpa_leaving_time",
                "--no-progress",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "28.57" in result.output
        assert "1.5" not in result.output
        assert "4.5" not in result.output


class TestDiffFilePairMarkdown:
    def test_markdown_output_contains_all_columns(
        self, runner, two_divergent_files, tmp_path
    ):
        path1, path2 = two_divergent_files
        out = tmp_path / "report.md"
        result = runner.invoke(
            cli,
            [
                "diff",
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


# ---------------------------------------------------------------------------
# Folder-pair mode
# ---------------------------------------------------------------------------


class TestDiffFolderPairSummary:
    def test_compact_summary_has_one_row_per_run(
        self, runner, two_folders_divergent
    ):
        dir1, dir2 = two_folders_divergent
        result = runner.invoke(
            cli, ["diff", str(dir1), str(dir2), "--no-progress"]
        )
        assert result.exit_code == 0, result.output
        # Compact summary columns.
        assert "Result" in result.output
        assert "Max % Diff" in result.output
        assert "Worst Table" in result.output
        # All three run codes appear.
        for rc in ("run01", "run02", "run03"):
            assert rc in result.output
        # run02 is the perturbed one — exactly one of the rows should
        # show a mismatch result.
        assert "Exact Match" in result.output
        assert "differ" in result.output


class TestDiffFolderPairVerbose:
    def test_verbose_emits_full_matrix(
        self, runner, two_folders_divergent, tmp_path
    ):
        dir1, dir2 = two_folders_divergent
        # Use markdown to bypass Rich's column truncation, so we can see
        # the full per-table headers.
        out = tmp_path / "report.md"
        result = runner.invoke(
            cli,
            [
                "diff",
                str(dir1),
                str(dir2),
                "--verbose",
                "--markdown",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        text = out.read_text()
        # Full matrix uses dataset labels as column headers.
        assert "microscale_out/" in text
        # All three run codes present as row labels.
        for rc in ("run01", "run02", "run03"):
            assert rc in text
        # The compact-summary columns should NOT appear.
        assert "Worst Table" not in text


class TestDiffFolderPairNoCommonCodes:
    def test_empty_intersection_errors(self, runner, tmp_path):
        dir1 = tmp_path / "x"
        dir2 = tmp_path / "y"
        dir1.mkdir()
        dir2.mkdir()
        _build_run_in_folder(dir1, "run_a")
        _build_run_in_folder(dir2, "run_b")
        result = runner.invoke(
            cli, ["diff", str(dir1), str(dir2), "--no-progress"]
        )
        assert result.exit_code != 0
        assert "No run codes in common" in result.output
