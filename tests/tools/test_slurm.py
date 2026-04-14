"""Unit tests for :mod:`lysis.tools.slurm`.

Tests cover:

* :func:`generate_micro_child_script` — single-tier and two-tier bash script generation
* :func:`submit_micro_child_job` — writing and submitting a script via ``gs.sbatch``
* :func:`wait_for_jobs` — polling squeue until jobs complete or fail
* :func:`submit_micro_slurm_job` — end-to-end staging, script generation, and submission
"""

import os
from pathlib import Path
from unittest.mock import MagicMock, call, patch

import h5py
import pytest

from lysis.config.constants import CONST
from lysis.config.parameters import MicroParameters
from lysis.tools.slurm import (
    generate_micro_child_script,
    submit_micro_child_job,
    submit_micro_slurm_job,
    wait_for_jobs,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_micro_hdf5(path: Path) -> None:
    """Write a minimal v2.0.0 HDF5 file with default MicroParameters."""
    mp = MicroParameters()
    with h5py.File(str(path), "w") as f:
        f.attrs[CONST.DATASPEC_VERSION_ATTR] = "v2.0.0"
        grp = f.require_group("micro_data")
        for k, v in mp.to_basedict().items():
            grp.attrs[k] = str(v) if not isinstance(v, (int, float, bool)) else v


@pytest.fixture
def micro_hdf5(tmp_path):
    """Minimal HDF5 file; returns (hdf5_path, run_code)."""
    h5_path = tmp_path / "run-01.h5"
    _write_micro_hdf5(h5_path)
    return h5_path


# ---------------------------------------------------------------------------
# TestGenerateMicroChildScript
# ---------------------------------------------------------------------------


class TestGenerateMicroChildScript:
    """Tests for :func:`generate_micro_child_script`."""

    def test_single_tier_contains_staging_dir(self, tmp_path):
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", ""
        )
        assert str(staging) in script

    def test_single_tier_no_local_workdir(self, tmp_path):
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", ""
        )
        assert "local_work_dir" not in script

    def test_single_tier_contains_binary_copy(self, tmp_path):
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", ""
        )
        assert "cp" in script
        assert "micro.exe" in script

    def test_single_tier_contains_hdf5_path(self, tmp_path):
        staging = tmp_path / "staging"
        h5_path = tmp_path / "run-01.h5"
        script = generate_micro_child_script(
            staging, "run-01", h5_path, "/bin/micro.exe", ""
        )
        assert str(h5_path) in script

    def test_single_tier_contains_out_code(self, tmp_path):
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "_mycode"
        )
        assert "_mycode" in script

    def test_single_tier_partition_in_sbatch_header(self, tmp_path):
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "", partition="long"
        )
        assert "partition" in script
        assert "long" in script

    def test_single_tier_no_partition_when_not_set(self, tmp_path):
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", ""
        )
        assert "--partition" not in script

    def test_two_tier_contains_fast_tmp_root(self, tmp_path):
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "",
            fast_tmp_root="/nvme/scratch",
        )
        assert "/nvme/scratch" in script

    def test_two_tier_contains_mktemp(self, tmp_path):
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "",
            fast_tmp_root="/nvme/scratch",
        )
        assert "mktemp" in script

    def test_two_tier_contains_mv_step(self, tmp_path):
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "",
            fast_tmp_root="/nvme/scratch",
        )
        assert "mv" in script

    def test_two_tier_contains_cleanup(self, tmp_path):
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "",
            fast_tmp_root="/nvme/scratch",
        )
        assert "rm -rf" in script

    def test_two_tier_contains_binary_copy(self, tmp_path):
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "",
            fast_tmp_root="/nvme/scratch",
        )
        assert "cp" in script
        assert "micro.exe" in script

    def test_two_tier_staging_dir_also_present(self, tmp_path):
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "",
            fast_tmp_root="/nvme/scratch",
        )
        assert str(staging) in script

    def test_returns_string(self, tmp_path):
        staging = tmp_path / "staging"
        result = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", ""
        )
        assert isinstance(result, str)

    def test_has_sbatch_shebang(self, tmp_path):
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", ""
        )
        assert script.startswith("#!/bin/bash")

    def test_single_tier_uses_fortran_micro_import(self, tmp_path):
        """Child script must import from fortran_micro, not the old codeutil."""
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", ""
        )
        assert "lysis.execution.fortran_micro" in script
        assert "codeutil" not in script

    def test_two_tier_uses_fortran_micro_import(self, tmp_path):
        """Two-tier child script must import from fortran_micro, not codeutil."""
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "",
            fast_tmp_root="/nvme/scratch",
        )
        assert "lysis.execution.fortran_micro" in script
        assert "codeutil" not in script


# ---------------------------------------------------------------------------
# TestSubmitMicroChildJob
# ---------------------------------------------------------------------------


class TestSubmitMicroChildJob:
    """Tests for :func:`submit_micro_child_job`."""

    @patch("lysis.tools.slurm.gs.sbatch", return_value=42)
    def test_returns_job_id(self, mock_sbatch, tmp_path):
        script_path = tmp_path / "child.sh"
        job_id = submit_micro_child_job("#!/bin/bash\necho hi", script_path)
        assert job_id == 42

    @patch("lysis.tools.slurm.gs.sbatch", return_value=7)
    def test_writes_script_to_path(self, mock_sbatch, tmp_path):
        script_path = tmp_path / "child.sh"
        content = "#!/bin/bash\necho hello"
        submit_micro_child_job(content, script_path)
        assert script_path.read_text() == content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=7)
    def test_sbatch_called_with_script_path(self, mock_sbatch, tmp_path):
        script_path = tmp_path / "child.sh"
        submit_micro_child_job("#!/bin/bash\necho hi", script_path)
        mock_sbatch.assert_called_once_with([str(script_path)])

    @patch("lysis.tools.slurm.gs.sbatch", return_value=7)
    def test_accepts_str_path(self, mock_sbatch, tmp_path):
        script_path = tmp_path / "child.sh"
        job_id = submit_micro_child_job("#!/bin/bash\necho hi", str(script_path))
        assert job_id == 7


# ---------------------------------------------------------------------------
# TestWaitForJobs
# ---------------------------------------------------------------------------


class TestWaitForJobs:
    """Tests for :func:`wait_for_jobs`."""

    @patch("lysis.tools.slurm.time.sleep")
    @patch("lysis.tools.slurm.gs.squeue")
    def test_returns_when_jobs_gone(self, mock_squeue, mock_sleep):
        """Should return normally when all job IDs leave the queue."""
        mock_squeue.read.side_effect = [
            [{"JOBID": "101", "STATE": "RUNNING"}],
            [],
        ]
        wait_for_jobs([101], poll_interval=1)
        assert mock_sleep.call_count >= 1

    @patch("lysis.tools.slurm.time.sleep")
    @patch("lysis.tools.slurm.gs.squeue")
    def test_no_error_when_jobs_complete_normally(self, mock_squeue, mock_sleep):
        """Should not raise when job exits queue with no FAILED/CANCELLED."""
        mock_squeue.read.side_effect = [
            [],  # job already gone on first poll
        ]
        # Should not raise
        wait_for_jobs([101], poll_interval=1)

    @patch("lysis.tools.slurm.time.sleep")
    @patch("lysis.tools.slurm.gs.squeue")
    def test_raises_on_failed_state(self, mock_squeue, mock_sleep):
        mock_squeue.read.return_value = [
            {"JOBID": "101", "STATE": "FAILED"}
        ]
        with pytest.raises(RuntimeError, match="FAILED"):
            wait_for_jobs([101], poll_interval=1)

    @patch("lysis.tools.slurm.time.sleep")
    @patch("lysis.tools.slurm.gs.squeue")
    def test_raises_on_cancelled_state(self, mock_squeue, mock_sleep):
        mock_squeue.read.return_value = [
            {"JOBID": "101", "STATE": "CANCELLED"}
        ]
        with pytest.raises(RuntimeError, match="CANCELLED"):
            wait_for_jobs([101], poll_interval=1)

    @patch("lysis.tools.slurm.time.sleep")
    @patch("lysis.tools.slurm.gs.squeue")
    def test_waits_for_multiple_jobs(self, mock_squeue, mock_sleep):
        """All job IDs must disappear before returning."""
        mock_squeue.read.side_effect = [
            # First poll: both jobs still running
            [
                {"JOBID": "101", "STATE": "RUNNING"},
                {"JOBID": "102", "STATE": "RUNNING"},
            ],
            # Second poll: only job 102 remains
            [{"JOBID": "102", "STATE": "RUNNING"}],
            # Third poll: all done
            [],
        ]
        wait_for_jobs([101, 102], poll_interval=1)
        assert mock_sleep.call_count == 3

    @patch("lysis.tools.slurm.time.sleep")
    @patch("lysis.tools.slurm.gs.squeue")
    def test_empty_job_list_returns_immediately(self, mock_squeue, mock_sleep):
        """An empty job list should return without calling sleep or squeue."""
        wait_for_jobs([], poll_interval=1)
        mock_sleep.assert_not_called()
        mock_squeue.read.assert_not_called()


# ---------------------------------------------------------------------------
# TestSubmitMicroSlurmJob
# ---------------------------------------------------------------------------


class TestSubmitMicroSlurmJob:
    """Tests for :func:`submit_micro_slurm_job`."""

    @patch("lysis.tools.slurm.gs.sbatch", return_value=999)
    def test_returns_master_job_id(self, mock_sbatch, micro_hdf5, tmp_path):
        job_id = submit_micro_slurm_job(
            micro_hdf5, "/bin/micro.exe", staging_root=tmp_path
        )
        assert job_id == 999

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_staging_dir_created_in_staging_root(self, mock_sbatch, micro_hdf5, tmp_path):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(micro_hdf5, "/bin/micro.exe",
                               staging_root=staging_root)
        # A unique subdirectory should have been created in staging_root
        children = list(staging_root.iterdir())
        assert len(children) == 1
        assert children[0].is_dir()
        assert "lysis-micro" in children[0].name

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_staging_dir_defaults_to_hdf5_parent(self, mock_sbatch, micro_hdf5, tmp_path):
        """When staging_root is None, staging dir goes in hdf5_path.parent."""
        submit_micro_slurm_job(micro_hdf5, "/bin/micro.exe")
        staging_dirs = [
            p for p in tmp_path.iterdir()
            if p.is_dir() and "lysis-micro" in p.name
        ]
        assert len(staging_dirs) >= 1

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_child_script_written_to_staging_dir(self, mock_sbatch, micro_hdf5, tmp_path):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(micro_hdf5, "/bin/micro.exe",
                               staging_root=staging_root)
        staging_dir = list(staging_root.iterdir())[0]
        assert (staging_dir / "child_000.sh").exists()

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_script_written_to_staging_dir(self, mock_sbatch, micro_hdf5, tmp_path):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(micro_hdf5, "/bin/micro.exe",
                               staging_root=staging_root)
        staging_dir = list(staging_root.iterdir())[0]
        assert (staging_dir / "master.sh").exists()
        assert (staging_dir / "master.py").exists()

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_sbatch_called_with_master_sh(self, mock_sbatch, micro_hdf5, tmp_path):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(micro_hdf5, "/bin/micro.exe",
                               staging_root=staging_root)
        staging_dir = list(staging_root.iterdir())[0]
        expected_path = str(staging_dir / "master.sh")
        mock_sbatch.assert_called_once_with([expected_path])

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_partition_in_master_script(self, mock_sbatch, micro_hdf5, tmp_path):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(micro_hdf5, "/bin/micro.exe",
                               staging_root=staging_root, partition="long")
        staging_dir = list(staging_root.iterdir())[0]
        master_content = (staging_dir / "master.sh").read_text()
        assert "partition" in master_content
        assert "long" in master_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_partition_forwarded_to_child_script(self, mock_sbatch, micro_hdf5, tmp_path):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(micro_hdf5, "/bin/micro.exe",
                               staging_root=staging_root, partition="long")
        staging_dir = list(staging_root.iterdir())[0]
        child_content = (staging_dir / "child_000.sh").read_text()
        assert "partition" in child_content
        assert "long" in child_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_fast_tmp_root_forwarded_to_child_script(self, mock_sbatch, micro_hdf5, tmp_path):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(micro_hdf5, "/bin/micro.exe",
                               staging_root=staging_root,
                               fast_tmp_root="/nvme/scratch")
        staging_dir = list(staging_root.iterdir())[0]
        child_content = (staging_dir / "child_000.sh").read_text()
        assert "/nvme/scratch" in child_content
        assert "mktemp" in child_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_staging_dir_name_contains_run_code(self, mock_sbatch, micro_hdf5, tmp_path):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(micro_hdf5, "/bin/micro.exe",
                               staging_root=staging_root)
        staging_dir = list(staging_root.iterdir())[0]
        assert "run-01" in staging_dir.name

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_py_imports_fortran_micro_module(self, mock_sbatch, micro_hdf5, tmp_path):
        """master.py must import from fortran_micro, not the old codeutil."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(micro_hdf5, "/bin/micro.exe",
                               staging_root=staging_root)
        staging_dir = list(staging_root.iterdir())[0]
        master_content = (staging_dir / "master.py").read_text()
        assert "lysis.execution.fortran_micro" in master_content
        assert "codeutil" not in master_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_py_calls_import_results_as_instance_method(self, mock_sbatch, micro_hdf5, tmp_path):
        """master.py must call import_results on an instance, not as a static method."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(micro_hdf5, "/bin/micro.exe",
                               staging_root=staging_root)
        staging_dir = list(staging_root.iterdir())[0]
        master_content = (staging_dir / "master.py").read_text()
        assert "FortranMicro.import_results(" not in master_content
        assert "fm.import_results(" in master_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_relative_executable_resolved_to_absolute(self, mock_sbatch, micro_hdf5, tmp_path):
        """Relative executable paths must be resolved to absolute in generated scripts."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        # Pass a relative path by making it relative to cwd
        import os
        rel_exe = os.path.relpath("/bin/micro.exe")
        submit_micro_slurm_job(micro_hdf5, rel_exe, staging_root=staging_root)
        staging_dir = list(staging_root.iterdir())[0]
        child_content = (staging_dir / "child_000.sh").read_text()
        assert "/bin/micro.exe" in child_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_relative_hdf5_path_resolved_to_absolute(self, mock_sbatch, micro_hdf5, tmp_path):
        """Relative HDF5 paths must be resolved to absolute in generated scripts."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        import os
        rel_h5 = os.path.relpath(str(micro_hdf5))
        submit_micro_slurm_job(rel_h5, "/bin/micro.exe", staging_root=staging_root)
        staging_dir = list(staging_root.iterdir())[0]
        child_content = (staging_dir / "child_000.sh").read_text()
        assert str(micro_hdf5) in child_content
