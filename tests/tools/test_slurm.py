"""Unit tests for :mod:`lysis.tools.slurm`.

Tests cover:

* :func:`generate_micro_child_script` — single-tier and two-tier bash script generation
* :func:`submit_micro_child_job` — writing and submitting a script via ``gs.sbatch``
* :func:`wait_for_jobs` — polling squeue until jobs complete or fail
* :func:`submit_micro_slurm_job` — end-to-end staging, script generation, and submission
"""

import os
import re
from pathlib import Path
from unittest.mock import MagicMock, call, patch

import h5py
import pytest

from lysis.config.constants import CONST
from lysis.config.parameters import MacroParameters, MicroParameters
from lysis.tools.slurm import (
    _python_macro_spec,
    _snapshot_lysis_src,
    apply_sbatch_overrides,
    generate_macro_array_script,
    generate_micro_array_script,
    generate_micro_child_script,
    parse_sbatch_tokens,
    submit_macro_slurm_job,
    submit_micro_child_job,
    submit_micro_slurm_job,
    submit_python_macro_slurm_job,
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

    def test_single_tier_no_runtime_binary_copy(self, tmp_path):
        """Single-tier must NOT cp the binary at runtime.

        The binary is pre-staged by submit_micro_slurm_job before sbatch
        is called; a runtime cp is redundant and, under --fortran-commit,
        actively broken because the build dir is already gone.
        """
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", ""
        )
        assert "cp " not in script
        # The binary name still appears — inside the python heredoc, as
        # the resolved staging-dir path passed to FortranMicro.from_hdf5.
        assert "micro.exe" in script
        assert f"{staging}/micro.exe" in script

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

    def test_two_tier_copies_binary_from_staging_not_original(self, tmp_path):
        """Two-tier cp must source the binary from staging_dir, not the
        caller's --executable path.

        Regression test for the --fortran-commit bug: the caller's path
        points into a temp build dir that is removed as soon as sbatch
        returns, so a later runtime cp from that path fails with
        ``cp: cannot stat …``.  The pre-staged copy in staging_dir is
        always present.
        """
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "",
            fast_tmp_root="/nvme/scratch",
        )
        assert f'cp "{staging}/micro.exe" "${{local_work_dir}}/"' in script
        assert "/bin/micro.exe" not in script

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

    def test_default_modules_in_preamble(self, tmp_path):
        """Without modules, the default LMod spec is loaded."""
        from lysis.tools.slurm import DEFAULT_MODULES
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "",
        )
        assert f"module load {DEFAULT_MODULES}" in script

    def test_custom_modules_in_preamble(self, tmp_path):
        """modules override must appear in the module-load line."""
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "",
            modules="intel-compilers/2024.2.0",
        )
        assert "module load intel-compilers/2024.2.0" in script
        assert "module load intel-compilers/2023" not in script

    def test_space_separated_modules_in_preamble(self, tmp_path):
        """A space-separated module list is forwarded verbatim to module load."""
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "",
            modules="intel-compilers/2024 SciPy-bundle/2023.07",
        )
        assert (
            "module load intel-compilers/2024 SciPy-bundle/2023.07" in script
        )


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

    @pytest.fixture(autouse=True)
    def _patch_copy2(self):
        """Patch shutil.copy2 so the executable pre-stage step doesn't need a real binary."""
        with patch("lysis.tools.slurm.shutil.copy2"):
            yield

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
        assert (staging_dir / "lysis-micro-child__run-01.sh").exists()

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_script_written_to_staging_dir(self, mock_sbatch, micro_hdf5, tmp_path):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(micro_hdf5, "/bin/micro.exe",
                               staging_root=staging_root)
        staging_dir = list(staging_root.iterdir())[0]
        assert (staging_dir / "lysis-micro-master__run-01.sh").exists()
        assert (staging_dir / "lysis-micro-master__run-01.py").exists()

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_sbatch_called_with_master_sh(self, mock_sbatch, micro_hdf5, tmp_path):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(micro_hdf5, "/bin/micro.exe",
                               staging_root=staging_root)
        staging_dir = list(staging_root.iterdir())[0]
        expected_path = str(staging_dir / "lysis-micro-master__run-01.sh")
        mock_sbatch.assert_called_once_with([expected_path])

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_partition_in_master_script(self, mock_sbatch, micro_hdf5, tmp_path):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(micro_hdf5, "/bin/micro.exe",
                               staging_root=staging_root, partition="long")
        staging_dir = list(staging_root.iterdir())[0]
        master_content = (staging_dir / "lysis-micro-master__run-01.sh").read_text()
        assert "partition" in master_content
        assert "long" in master_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_partition_forwarded_to_child_script(self, mock_sbatch, micro_hdf5, tmp_path):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(micro_hdf5, "/bin/micro.exe",
                               staging_root=staging_root, partition="long")
        staging_dir = list(staging_root.iterdir())[0]
        child_content = (staging_dir / "lysis-micro-child__run-01.sh").read_text()
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
        child_content = (staging_dir / "lysis-micro-child__run-01.sh").read_text()
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
        master_content = (staging_dir / "lysis-micro-master__run-01.py").read_text()
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
        master_content = (staging_dir / "lysis-micro-master__run-01.py").read_text()
        assert "FortranMicro.import_results(" not in master_content
        assert "fm.import_results(" in master_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_py_passes_dispatcher_log_dir(
        self, mock_sbatch, micro_hdf5, tmp_path
    ):
        """master.py passes dispatcher_log_dir=STAGING_DIR so worker .out logs
        are ingested into log_files/dispatcher/."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(
            micro_hdf5, "/bin/micro.exe", staging_root=staging_root
        )
        staging_dir = list(staging_root.iterdir())[0]
        master_content = (staging_dir / "lysis-micro-master__run-01.py").read_text()
        assert "dispatcher_log_dir=STAGING_DIR" in master_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_py_passes_executable_to_runner(
        self, mock_sbatch, micro_hdf5, tmp_path
    ):
        """master.py must construct FortranMicro with executable= so import_results can stamp binary_* attrs."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(micro_hdf5, "/bin/micro.exe",
                               staging_root=staging_root)
        staging_dir = list(staging_root.iterdir())[0]
        master_content = (staging_dir / "lysis-micro-master__run-01.py").read_text()
        assert "BINARY_NAME = 'micro.exe'" in master_content
        assert "executable=str(STAGING_DIR / BINARY_NAME)" in master_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_py_historical_attrs_absent_by_default(
        self, mock_sbatch, micro_hdf5, tmp_path
    ):
        """Without --fortran-commit, the master script keeps HISTORICAL_BACKEND_ATTRS = None."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(micro_hdf5, "/bin/micro.exe",
                               staging_root=staging_root)
        staging_dir = list(staging_root.iterdir())[0]
        master_content = (staging_dir / "lysis-micro-master__run-01.py").read_text()
        assert "HISTORICAL_BACKEND_ATTRS = None" in master_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_py_bakes_historical_attrs(
        self, mock_sbatch, micro_hdf5, tmp_path
    ):
        """When historical_backend_attrs is given, the dict is repr'd into the master script."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        prov = {
            "backend_commit": "a" * 40,
            "backend_dirty": "clean",
            "backend_compiler": "GCC 11.4.0",
            "backend_historical": True,
        }
        submit_micro_slurm_job(
            micro_hdf5, "/bin/micro.exe",
            staging_root=staging_root,
            historical_backend_attrs=prov,
        )
        staging_dir = list(staging_root.iterdir())[0]
        master_content = (staging_dir / "lysis-micro-master__run-01.py").read_text()
        # The literal dict must round-trip through repr() into the script.
        assert "HISTORICAL_BACKEND_ATTRS = {" in master_content
        assert "'backend_historical': True" in master_content
        # And the runner construction must opt in to skip_binary_verification.
        assert (
            "skip_binary_verification=(HISTORICAL_BACKEND_ATTRS is not None)"
            in master_content
        )

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_relative_executable_resolved_to_absolute(self, mock_sbatch, micro_hdf5, tmp_path):
        """Relative executable paths must be resolved before pre-staging.

        Generated single-tier scripts reference the pre-staged copy at
        ``{staging_dir}/{binary_basename}`` (not the caller's --executable
        path).  Resolution still matters because ``submit_micro_slurm_job``
        uses the resolved path for ``shutil.copy2(executable, staging_dir)``;
        the test checks the binary basename appears at an absolute path
        inside the script.
        """
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        # Pass a relative path by making it relative to cwd
        import os
        rel_exe = os.path.relpath("/bin/micro.exe")
        submit_micro_slurm_job(micro_hdf5, rel_exe, staging_root=staging_root)
        staging_dir = list(staging_root.iterdir())[0]
        child_content = (staging_dir / "lysis-micro-child__run-01.sh").read_text()
        # The staged binary is referenced by absolute path inside the heredoc.
        assert f"{staging_dir}/micro.exe" in child_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_relative_hdf5_path_resolved_to_absolute(self, mock_sbatch, micro_hdf5, tmp_path):
        """Relative HDF5 paths must be resolved to absolute in generated scripts."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        import os
        rel_h5 = os.path.relpath(str(micro_hdf5))
        submit_micro_slurm_job(rel_h5, "/bin/micro.exe", staging_root=staging_root)
        staging_dir = list(staging_root.iterdir())[0]
        child_content = (staging_dir / "lysis-micro-child__run-01.sh").read_text()
        assert str(micro_hdf5) in child_content


# ---------------------------------------------------------------------------
# Macro Slurm helpers and fixtures
# ---------------------------------------------------------------------------


def _write_macro_hdf5(path: Path) -> None:
    """Write a minimal v2.0.0 HDF5 with default Micro and MacroParameters.

    Includes a non-empty ``micro_data/tpa_leaving_time`` to satisfy the
    post-run-micro state required by :meth:`FortranMacro.from_hdf5`.
    """
    import numpy as np

    mp = MicroParameters()
    mcp = MacroParameters(micro_params=mp)
    with h5py.File(str(path), "w") as f:
        f.attrs[CONST.DATASPEC_VERSION_ATTR] = "v2.0.0"
        micro_grp = f.require_group("micro_data")
        for k, v in mp.to_basedict().items():
            micro_grp.attrs[k] = str(v) if not isinstance(v, (int, float, bool)) else v
        f.create_dataset(
            "micro_data/tpa_leaving_time",
            data=np.array([1.0]),
            dtype=np.float64,
        )
        macro_grp = f.require_group("macro_data")
        for k, v in mcp.to_basedict().items():
            macro_grp.attrs[k] = str(v) if not isinstance(v, (int, float, bool)) else v


@pytest.fixture
def macro_hdf5(tmp_path):
    """Minimal HDF5 file with both Micro and MacroParameters."""
    h5_path = tmp_path / "run-01.h5"
    _write_macro_hdf5(h5_path)
    return h5_path


# ---------------------------------------------------------------------------
# TestGenerateMacroArrayScript
# ---------------------------------------------------------------------------


class TestGenerateMacroArrayScript:
    """Tests for :func:`generate_macro_array_script`."""

    def test_returns_string(self, tmp_path):
        script = generate_macro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/macro.exe", n_sims=5,
        )
        assert isinstance(script, str)

    def test_has_sbatch_shebang(self, tmp_path):
        script = generate_macro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/macro.exe", n_sims=5,
        )
        assert script.startswith("#!/bin/bash")

    def test_array_directive_present(self, tmp_path):
        """Script must contain #SBATCH --array 0-{n_sims-1}."""
        script = generate_macro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/macro.exe", n_sims=5,
        )
        assert "--array" in script
        assert "0-4" in script

    def test_single_tier_no_local_workdir(self, tmp_path):
        """Single-tier script must not use a local working directory."""
        script = generate_macro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/macro.exe", n_sims=3,
        )
        assert "local_work_dir" not in script

    def test_single_tier_uses_fortran_macro_import(self, tmp_path):
        script = generate_macro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/macro.exe", n_sims=3,
        )
        assert "lysis.execution.fortran_macro" in script

    def test_single_tier_uses_slurm_array_task_id(self, tmp_path):
        script = generate_macro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/macro.exe", n_sims=3,
        )
        assert "SLURM_ARRAY_TASK_ID" in script

    def test_single_tier_contains_in_code(self, tmp_path):
        script = generate_macro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/macro.exe", n_sims=3,
            in_code="_incode",
        )
        assert "_incode" in script

    def test_single_tier_contains_out_code(self, tmp_path):
        script = generate_macro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/macro.exe", n_sims=3,
            out_code="_outcode",
        )
        assert "_outcode" in script

    def test_single_tier_partition_in_header(self, tmp_path):
        script = generate_macro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/macro.exe", n_sims=3,
            partition="long",
        )
        assert "partition" in script
        assert "long" in script

    def test_single_tier_no_partition_when_not_set(self, tmp_path):
        script = generate_macro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/macro.exe", n_sims=3,
        )
        assert "--partition" not in script

    def test_two_tier_contains_mktemp(self, tmp_path):
        script = generate_macro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/macro.exe", n_sims=3,
            fast_tmp_root="/nvme/scratch",
        )
        assert "mktemp" in script

    def test_two_tier_contains_mv_step(self, tmp_path):
        script = generate_macro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/macro.exe", n_sims=3,
            fast_tmp_root="/nvme/scratch",
        )
        assert "mv" in script

    def test_two_tier_contains_cleanup(self, tmp_path):
        script = generate_macro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/macro.exe", n_sims=3,
            fast_tmp_root="/nvme/scratch",
        )
        assert "rm -rf" in script

    def test_two_tier_contains_fast_tmp_root(self, tmp_path):
        script = generate_macro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/macro.exe", n_sims=3,
            fast_tmp_root="/nvme/scratch",
        )
        assert "/nvme/scratch" in script

    def test_single_tier_no_cp_of_executable(self, tmp_path):
        """Single-tier script must not copy the executable (pre-staged by submit_macro_slurm_job).

        Regression test: all array tasks share the same staging dir, so any
        ``cp`` in the array script causes a race where tasks 1..N fail with
        "File exists" after task 0 wins the copy.
        """
        script = generate_macro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/macro.exe", n_sims=3,
        )
        assert "cp " not in script


# ---------------------------------------------------------------------------
# TestSubmitMacroSlurmJob
# ---------------------------------------------------------------------------


class TestSubmitMacroSlurmJob:
    """Tests for :func:`submit_macro_slurm_job`."""

    @pytest.fixture
    def mock_write_setup(self):
        """Patch FortranMacro._write_setup_files and shutil.copy2 to avoid needing real files."""
        with patch(
            "lysis.execution.fortran_macro.FortranMacro._write_setup_files"
        ) as mock, patch("lysis.tools.slurm.shutil.copy2"):
            yield mock

    @patch("lysis.tools.slurm.gs.sbatch", return_value=999)
    def test_returns_master_job_id(self, mock_sbatch, macro_hdf5, tmp_path, mock_write_setup):
        job_id = submit_macro_slurm_job(
            macro_hdf5, "/bin/macro.exe", staging_root=tmp_path
        )
        assert job_id == 999

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_staging_dir_created_in_staging_root(
        self, mock_sbatch, macro_hdf5, tmp_path, mock_write_setup
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_macro_slurm_job(
            macro_hdf5, "/bin/macro.exe", staging_root=staging_root
        )
        children = list(staging_root.iterdir())
        assert len(children) == 1
        assert children[0].is_dir()
        assert "lysis-macro" in children[0].name

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_staging_dir_name_contains_run_code(
        self, mock_sbatch, macro_hdf5, tmp_path, mock_write_setup
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_macro_slurm_job(
            macro_hdf5, "/bin/macro.exe", staging_root=staging_root
        )
        staging_dir = list(staging_root.iterdir())[0]
        assert "run-01" in staging_dir.name

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_array_script_written(
        self, mock_sbatch, macro_hdf5, tmp_path, mock_write_setup
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_macro_slurm_job(
            macro_hdf5, "/bin/macro.exe", staging_root=staging_root
        )
        staging_dir = list(staging_root.iterdir())[0]
        assert (staging_dir / "lysis-macro-array__run-01.sh").exists()

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_scripts_written(
        self, mock_sbatch, macro_hdf5, tmp_path, mock_write_setup
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_macro_slurm_job(
            macro_hdf5, "/bin/macro.exe", staging_root=staging_root
        )
        staging_dir = list(staging_root.iterdir())[0]
        assert (staging_dir / "lysis-macro-master__run-01.sh").exists()
        assert (staging_dir / "lysis-macro-master__run-01.py").exists()

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_sbatch_called_with_master_sh(
        self, mock_sbatch, macro_hdf5, tmp_path, mock_write_setup
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_macro_slurm_job(
            macro_hdf5, "/bin/macro.exe", staging_root=staging_root
        )
        staging_dir = list(staging_root.iterdir())[0]
        expected_path = str(staging_dir / "lysis-macro-master__run-01.sh")
        mock_sbatch.assert_called_once_with([expected_path])

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_array_sh_has_array_directive(
        self, mock_sbatch, macro_hdf5, tmp_path, mock_write_setup
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_macro_slurm_job(
            macro_hdf5, "/bin/macro.exe", staging_root=staging_root
        )
        staging_dir = list(staging_root.iterdir())[0]
        array_content = (staging_dir / "lysis-macro-array__run-01.sh").read_text()
        assert "--array" in array_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_py_imports_fortran_macro(
        self, mock_sbatch, macro_hdf5, tmp_path, mock_write_setup
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_macro_slurm_job(
            macro_hdf5, "/bin/macro.exe", staging_root=staging_root
        )
        staging_dir = list(staging_root.iterdir())[0]
        master_content = (staging_dir / "lysis-macro-master__run-01.py").read_text()
        assert "lysis.execution.fortran_macro" in master_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_py_passes_dispatcher_log_dir(
        self, mock_sbatch, macro_hdf5, tmp_path, mock_write_setup
    ):
        """master.py passes dispatcher_log_dir=STAGING_DIR so worker .out logs
        are ingested into log_files/dispatcher/."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_macro_slurm_job(
            macro_hdf5, "/bin/macro.exe", staging_root=staging_root
        )
        staging_dir = list(staging_root.iterdir())[0]
        master_content = (staging_dir / "lysis-macro-master__run-01.py").read_text()
        assert "dispatcher_log_dir=STAGING_DIR" in master_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_py_passes_executable_to_runner(
        self, mock_sbatch, macro_hdf5, tmp_path, mock_write_setup
    ):
        """master.py must construct FortranMacro with executable= so import_results can stamp binary_* attrs."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_macro_slurm_job(
            macro_hdf5, "/bin/macro.exe", staging_root=staging_root
        )
        staging_dir = list(staging_root.iterdir())[0]
        master_content = (staging_dir / "lysis-macro-master__run-01.py").read_text()
        assert "BINARY_NAME = 'macro.exe'" in master_content
        assert "executable=str(STAGING_DIR / BINARY_NAME)" in master_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_py_bakes_historical_attrs(
        self, mock_sbatch, macro_hdf5, tmp_path, mock_write_setup
    ):
        """Macro slurm path: historical_backend_attrs dict gets repr'd into the master script."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        prov = {
            "backend_commit": "b" * 40,
            "backend_dirty": "clean",
            "backend_compiler": "Intel(R) Fortran 2023.0",
            "backend_historical": True,
        }
        submit_macro_slurm_job(
            macro_hdf5, "/bin/macro.exe",
            staging_root=staging_root,
            historical_backend_attrs=prov,
        )
        staging_dir = list(staging_root.iterdir())[0]
        master_content = (staging_dir / "lysis-macro-master__run-01.py").read_text()
        assert "HISTORICAL_BACKEND_ATTRS = {" in master_content
        assert "'backend_historical': True" in master_content
        assert (
            "skip_binary_verification=(HISTORICAL_BACKEND_ATTRS is not None)"
            in master_content
        )

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_partition_in_master_script(
        self, mock_sbatch, macro_hdf5, tmp_path, mock_write_setup
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_macro_slurm_job(
            macro_hdf5, "/bin/macro.exe",
            staging_root=staging_root, partition="long"
        )
        staging_dir = list(staging_root.iterdir())[0]
        master_content = (staging_dir / "lysis-macro-master__run-01.sh").read_text()
        assert "partition" in master_content
        assert "long" in master_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_write_setup_files_called(
        self, mock_sbatch, macro_hdf5, tmp_path, mock_write_setup
    ):
        """setup files must be pre-staged in the staging data dir."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_macro_slurm_job(
            macro_hdf5, "/bin/macro.exe", staging_root=staging_root
        )
        assert mock_write_setup.call_count == 1

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_executable_pre_staged_in_staging_dir(
        self, mock_sbatch, macro_hdf5, tmp_path
    ):
        """Executable must be copied to staging dir before the array job is submitted.

        Regression test: without pre-staging, each array task races to cp the
        binary to the shared staging dir and all but the first fail with
        "File exists".
        """
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        fake_exe = tmp_path / "macro.exe"
        fake_exe.write_bytes(b"fake binary")
        with patch("lysis.execution.fortran_macro.FortranMacro._write_setup_files"):
            submit_macro_slurm_job(macro_hdf5, str(fake_exe), staging_root=staging_root)
        staging_dir = list(staging_root.iterdir())[0]
        assert (staging_dir / "macro.exe").exists()

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_py_filters_squeue_by_array_job_id_equality(
        self, mock_sbatch, macro_hdf5, tmp_path, mock_write_setup
    ):
        """Master script must use ``==`` not ``startswith`` for ARRAY_JOB_ID.

        Regression test: ``str(123).startswith(str(123))`` is True, but it
        also matches array_job_id 1234 / 12345 / etc. — a prefix collision
        that would block the master from ever exiting the polling loop on
        large clusters where neighbouring jobs share a prefix.
        """
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_macro_slurm_job(
            macro_hdf5, "/bin/macro.exe", staging_root=staging_root
        )
        staging_dir = list(staging_root.iterdir())[0]
        master_py = (staging_dir / "lysis-macro-master__run-01.py").read_text()
        assert "ARRAY_JOB_ID" in master_py
        assert ".startswith(" not in master_py
        assert 'row["ARRAY_JOB_ID"] == str(array_job_id)' in master_py


# ---------------------------------------------------------------------------
# TestGenerateMicroArrayScript
# ---------------------------------------------------------------------------


class TestGenerateMicroArrayScript:
    """Tests for :func:`generate_micro_array_script`."""

    def test_returns_string(self, tmp_path):
        script = generate_micro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/micro.exe", num_children=10,
        )
        assert isinstance(script, str)

    def test_has_sbatch_shebang(self, tmp_path):
        script = generate_micro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/micro.exe", num_children=10,
        )
        assert script.startswith("#!/bin/bash")

    def test_array_directive_uses_num_children(self, tmp_path):
        script = generate_micro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/micro.exe", num_children=10,
        )
        assert "--array" in script
        assert "0-9" in script

    def test_single_tier_no_local_workdir(self, tmp_path):
        script = generate_micro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/micro.exe", num_children=10,
        )
        assert "local_work_dir" not in script

    def test_single_tier_uses_fortran_micro_import(self, tmp_path):
        script = generate_micro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/micro.exe", num_children=10,
        )
        assert "lysis.execution.fortran_micro" in script
        assert "FortranMicro" in script

    def test_single_tier_uses_slurm_array_task_id(self, tmp_path):
        script = generate_micro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/micro.exe", num_children=10,
        )
        assert "SLURM_ARRAY_TASK_ID" in script

    def test_single_tier_threads_num_children_into_from_hdf5(self, tmp_path):
        """The python -c body must pass num_children=K to from_hdf5."""
        script = generate_micro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/micro.exe", num_children=10,
        )
        assert "num_children=10" in script

    def test_single_tier_contains_out_code(self, tmp_path):
        script = generate_micro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/micro.exe", num_children=4,
            out_code="_outcode",
        )
        assert "_outcode" in script

    def test_single_tier_partition_in_header(self, tmp_path):
        script = generate_micro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/micro.exe", num_children=4,
            partition="long",
        )
        assert "partition" in script
        assert "long" in script

    def test_single_tier_no_partition_when_not_set(self, tmp_path):
        script = generate_micro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/micro.exe", num_children=4,
        )
        assert "--partition" not in script

    def test_single_tier_no_cp_of_executable(self, tmp_path):
        """Single-tier script must not copy the executable (pre-staged by submit_micro_slurm_job).

        Regression test: all array tasks share the same staging dir, so any
        ``cp`` in the array script causes a race where tasks 1..N fail with
        "File exists" after task 0 wins the copy.
        """
        script = generate_micro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/micro.exe", num_children=4,
        )
        assert "cp " not in script

    def test_two_tier_contains_mktemp(self, tmp_path):
        script = generate_micro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/micro.exe", num_children=4,
            fast_tmp_root="/nvme/scratch",
        )
        assert "mktemp" in script

    def test_two_tier_contains_mv_step(self, tmp_path):
        script = generate_micro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/micro.exe", num_children=4,
            fast_tmp_root="/nvme/scratch",
        )
        assert "mv " in script

    def test_two_tier_mv_uses_glob_not_per_sim_subdir(self, tmp_path):
        """Micro array tasks write flat files; the mv must glob, not target {SIM}."""
        script = generate_micro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/micro.exe", num_children=4,
            fast_tmp_root="/nvme/scratch",
        )
        assert '"${local_datadir}"/*' in script
        assert "${SIM}" not in script.split("# Move")[1]

    def test_two_tier_contains_cleanup(self, tmp_path):
        script = generate_micro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/micro.exe", num_children=4,
            fast_tmp_root="/nvme/scratch",
        )
        assert "rm -rf" in script

    def test_two_tier_contains_fast_tmp_root(self, tmp_path):
        script = generate_micro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/micro.exe", num_children=4,
            fast_tmp_root="/nvme/scratch",
        )
        assert "/nvme/scratch" in script


# ---------------------------------------------------------------------------
# TestSubmitMicroSlurmJob — array path (num_children >= 1)
# ---------------------------------------------------------------------------


class TestSubmitMicroSlurmJobArrayPath:
    """Tests for the array branch of :func:`submit_micro_slurm_job`.

    Activated by passing ``num_children >= 1``; these complement the
    legacy single-child tests in :class:`TestSubmitMicroSlurmJob`.
    """

    @pytest.fixture
    def mock_write_setup(self):
        """Patch FortranMicro._write_setup_files and shutil.copy2 to avoid real files."""
        with patch(
            "lysis.execution.fortran_micro.FortranMicro._write_setup_files"
        ) as mock, patch("lysis.tools.slurm.shutil.copy2"):
            yield mock

    @patch("lysis.tools.slurm.gs.sbatch", return_value=999)
    def test_returns_master_job_id(
        self, mock_sbatch, micro_hdf5, tmp_path, mock_write_setup
    ):
        job_id = submit_micro_slurm_job(
            micro_hdf5, "/bin/micro.exe",
            staging_root=tmp_path, num_children=10,
        )
        assert job_id == 999

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_array_script_written(
        self, mock_sbatch, micro_hdf5, tmp_path, mock_write_setup
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(
            micro_hdf5, "/bin/micro.exe",
            staging_root=staging_root, num_children=10,
        )
        staging_dir = list(staging_root.iterdir())[0]
        assert (staging_dir / "lysis-micro-array__run-01.sh").exists()
        # Legacy child script must NOT be written in array mode
        assert not (staging_dir / "lysis-micro-child__run-01.sh").exists()

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_scripts_written(
        self, mock_sbatch, micro_hdf5, tmp_path, mock_write_setup
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(
            micro_hdf5, "/bin/micro.exe",
            staging_root=staging_root, num_children=10,
        )
        staging_dir = list(staging_root.iterdir())[0]
        assert (staging_dir / "lysis-micro-master__run-01.sh").exists()
        assert (staging_dir / "lysis-micro-master__run-01.py").exists()

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_py_calls_concatenate(
        self, mock_sbatch, micro_hdf5, tmp_path, mock_write_setup
    ):
        """Master must concat per-task outputs before importing."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(
            micro_hdf5, "/bin/micro.exe",
            staging_root=staging_root, num_children=10,
        )
        staging_dir = list(staging_root.iterdir())[0]
        master_py = (staging_dir / "lysis-micro-master__run-01.py").read_text()
        assert "concatenate_child_outputs" in master_py
        assert "num_children=10" in master_py

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_array_sh_has_array_directive(
        self, mock_sbatch, micro_hdf5, tmp_path, mock_write_setup
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(
            micro_hdf5, "/bin/micro.exe",
            staging_root=staging_root, num_children=10,
        )
        staging_dir = list(staging_root.iterdir())[0]
        array_content = (staging_dir / "lysis-micro-array__run-01.sh").read_text()
        assert "--array" in array_content
        assert "0-9" in array_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_py_imports_fortran_micro(
        self, mock_sbatch, micro_hdf5, tmp_path, mock_write_setup
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(
            micro_hdf5, "/bin/micro.exe",
            staging_root=staging_root, num_children=10,
        )
        staging_dir = list(staging_root.iterdir())[0]
        master_content = (staging_dir / "lysis-micro-master__run-01.py").read_text()
        assert "lysis.execution.fortran_micro" in master_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_py_passes_executable_to_runner(
        self, mock_sbatch, micro_hdf5, tmp_path, mock_write_setup
    ):
        """master.py must construct FortranMicro with executable= so import_results can stamp binary_* attrs."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(
            micro_hdf5, "/bin/micro.exe",
            staging_root=staging_root, num_children=10,
        )
        staging_dir = list(staging_root.iterdir())[0]
        master_content = (staging_dir / "lysis-micro-master__run-01.py").read_text()
        assert "BINARY_NAME = 'micro.exe'" in master_content
        assert "executable=str(STAGING_DIR / BINARY_NAME)" in master_content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_executable_pre_staged_in_staging_dir(
        self, mock_sbatch, micro_hdf5, tmp_path
    ):
        """Executable must be copied to staging dir before the array tasks run.

        Regression test: identical to the macro version — without
        pre-staging, each array task would race to ``cp`` the binary
        and all but the first would fail with "File exists".
        """
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        fake_exe = tmp_path / "micro.exe"
        fake_exe.write_bytes(b"fake binary")
        with patch("lysis.execution.fortran_micro.FortranMicro._write_setup_files"):
            submit_micro_slurm_job(
                micro_hdf5, str(fake_exe),
                staging_root=staging_root, num_children=10,
            )
        staging_dir = list(staging_root.iterdir())[0]
        assert (staging_dir / "micro.exe").exists()

    def test_num_children_zero_raises(self, micro_hdf5, tmp_path):
        with pytest.raises(ValueError):
            submit_micro_slurm_job(
                micro_hdf5, "/bin/micro.exe",
                staging_root=tmp_path, num_children=0,
            )

    def test_num_children_negative_raises(self, micro_hdf5, tmp_path):
        with pytest.raises(ValueError):
            submit_micro_slurm_job(
                micro_hdf5, "/bin/micro.exe",
                staging_root=tmp_path, num_children=-1,
            )

    def test_num_children_exceeds_micro_simulations_raises(
        self, micro_hdf5, tmp_path
    ):
        """Cannot split N simulations across more than N tasks."""
        # Default MicroParameters().micro_simulations == 50_000.
        with pytest.raises(ValueError, match="exceeds micro_simulations"):
            submit_micro_slurm_job(
                micro_hdf5, "/bin/micro.exe",
                staging_root=tmp_path, num_children=50_001,
            )

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_num_children_none_uses_legacy_single_child_path(
        self, mock_sbatch, micro_hdf5, tmp_path, mock_write_setup
    ):
        """num_children=None (the default) keeps the legacy single-child path."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(
            micro_hdf5, "/bin/micro.exe",
            staging_root=staging_root,  # no num_children
        )
        staging_dir = list(staging_root.iterdir())[0]
        assert (staging_dir / "lysis-micro-child__run-01.sh").exists()
        assert not (staging_dir / "lysis-micro-array__run-01.sh").exists()


# ---------------------------------------------------------------------------
# TestParseSbatchTokens
# ---------------------------------------------------------------------------


class TestParseSbatchTokens:
    """Tests for :func:`parse_sbatch_tokens`."""

    def test_key_value_pair(self):
        assert parse_sbatch_tokens(["mem=4G"]) == {"mem": "4G"}

    def test_bare_key_is_flag_style(self):
        assert parse_sbatch_tokens(["hold"]) == {"hold": ""}

    def test_caret_marks_removal(self):
        assert parse_sbatch_tokens(["^exclusive=user"]) == {"exclusive=user": None}

    def test_value_can_contain_equals(self):
        """Split is on the FIRST '=' only; value keeps any remaining '='."""
        assert parse_sbatch_tokens(["time=01:00:00"]) == {"time": "01:00:00"}
        assert parse_sbatch_tokens(["foo=a=b=c"]) == {"foo": "a=b=c"}

    def test_users_full_example(self):
        """End-to-end: mirror the user's original dict-literal example."""
        tokens = ["mem=4G", "hold", "^exclusive=user"]
        assert parse_sbatch_tokens(tokens) == {
            "mem": "4G",
            "hold": "",
            "exclusive=user": None,
        }

    def test_empty_token_raises(self):
        with pytest.raises(ValueError):
            parse_sbatch_tokens([""])

    def test_empty_key_in_set_raises(self):
        with pytest.raises(ValueError):
            parse_sbatch_tokens(["=4G"])

    def test_caret_alone_raises(self):
        with pytest.raises(ValueError):
            parse_sbatch_tokens(["^"])

    def test_empty_input_returns_empty(self):
        assert parse_sbatch_tokens([]) == {}


# ---------------------------------------------------------------------------
# TestApplySbatchOverrides
# ---------------------------------------------------------------------------


class TestApplySbatchOverrides:
    """Tests for :func:`apply_sbatch_overrides`."""

    def test_none_overrides_returns_base_copy(self):
        base = {"mem": 3096, "exclusive=user": ""}
        result = apply_sbatch_overrides(base, None)
        assert result == base
        assert result is not base

    def test_empty_overrides_returns_base_copy(self):
        base = {"mem": 3096}
        assert apply_sbatch_overrides(base, {}) == base

    def test_override_sets_value(self):
        result = apply_sbatch_overrides({"mem": 3096}, {"mem": "4G"})
        assert result["mem"] == "4G"

    def test_override_adds_new_key(self):
        result = apply_sbatch_overrides({"mem": 3096}, {"hold": ""})
        assert result == {"mem": 3096, "hold": ""}

    def test_none_value_removes_key(self):
        result = apply_sbatch_overrides(
            {"mem": 3096, "exclusive=user": ""},
            {"exclusive=user": None},
        )
        assert "exclusive=user" not in result
        assert result == {"mem": 3096}

    def test_none_value_on_missing_key_is_noop(self):
        """Removing a key that isn't in base must not raise."""
        result = apply_sbatch_overrides({"mem": 3096}, {"nonexistent": None})
        assert result == {"mem": 3096}

    def test_does_not_mutate_base(self):
        base = {"mem": 3096}
        apply_sbatch_overrides(base, {"mem": "4G", "hold": ""})
        assert base == {"mem": 3096}

    def test_users_full_example_end_to_end(self):
        """parse + apply on the user's example produces the documented result."""
        base = {
            "array": "0-9",
            "job-name": "foo",
            "mem": 3096,
            "exclusive=user": "",
        }
        overrides = parse_sbatch_tokens(["mem=4G", "hold", "^exclusive=user"])
        result = apply_sbatch_overrides(base, overrides)
        assert result == {
            "array": "0-9",
            "job-name": "foo",
            "mem": "4G",
            "hold": "",
        }


# ---------------------------------------------------------------------------
# TestSbatchOverridesInGeneratedScripts
# ---------------------------------------------------------------------------


class TestSbatchOverridesInGeneratedScripts:
    """End-to-end: sbatch_overrides flows into the rendered #SBATCH header."""

    def test_macro_array_override_adds_flag(self, tmp_path):
        script = generate_macro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/macro.exe", n_sims=3,
            sbatch_overrides={"hold": ""},
        )
        assert "#SBATCH --hold" in script

    def test_macro_array_override_removes_default(self, tmp_path):
        """Removing the built-in exclusive=user default drops it from the header."""
        script = generate_macro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/macro.exe", n_sims=3,
            sbatch_overrides={"exclusive=user": None},
        )
        assert "--exclusive=user" not in script

    def test_micro_array_override_propagates(self, tmp_path):
        script = generate_micro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/micro.exe", num_children=4,
            sbatch_overrides={"hold": "", "exclusive=user": None},
        )
        assert "#SBATCH --hold" in script
        assert "--exclusive=user" not in script

    def test_micro_child_override_propagates(self, tmp_path):
        script = generate_micro_child_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/micro.exe",
            sbatch_overrides={"hold": ""},
        )
        assert "#SBATCH --hold" in script

    def test_no_overrides_leaves_defaults_intact(self, tmp_path):
        """Default behaviour is preserved when sbatch_overrides is None."""
        script = generate_macro_array_script(
            tmp_path / "staging", "run-01",
            tmp_path / "run-01.h5", "/bin/macro.exe", n_sims=3,
        )
        assert "--exclusive=user" in script
        assert "#SBATCH --hold" not in script


# ---------------------------------------------------------------------------
# TestSbatchOverridesInSubmitMicroSlurmJob
# ---------------------------------------------------------------------------


class TestSbatchOverridesInSubmitMicroSlurmJob:
    """Tests that sbatch_overrides reaches BOTH master and child/array scripts."""

    @pytest.fixture(autouse=True)
    def _patch_copy2(self):
        with patch("lysis.tools.slurm.shutil.copy2"):
            yield

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_legacy_single_child_overrides_master_and_child(
        self, mock_sbatch, micro_hdf5, tmp_path,
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(
            micro_hdf5, "/bin/micro.exe",
            staging_root=staging_root,
            sbatch_overrides={"hold": "", "exclusive=user": None},
        )
        staging_dir = list(staging_root.iterdir())[0]
        master = (staging_dir / "lysis-micro-master__run-01.sh").read_text()
        child = (staging_dir / "lysis-micro-child__run-01.sh").read_text()
        for content in (master, child):
            assert "#SBATCH --hold" in content
            assert "--exclusive=user" not in content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_array_path_overrides_master_and_array(
        self, mock_sbatch, micro_hdf5, tmp_path,
    ):
        # _write_setup_files needs to run; the default fixture's HDF5 file
        # has micro_simulations=400 so num_children up to 400 is valid.
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        with patch.object(
            __import__("lysis.execution.fortran_micro", fromlist=["FortranMicro"]).FortranMicro,
            "_write_setup_files",
        ):
            submit_micro_slurm_job(
                micro_hdf5, "/bin/micro.exe",
                staging_root=staging_root,
                num_children=4,
                sbatch_overrides={"hold": "", "exclusive=user": None},
            )
        staging_dir = list(staging_root.iterdir())[0]
        master = (staging_dir / "lysis-micro-master__run-01.sh").read_text()
        array = (staging_dir / "lysis-micro-array__run-01.sh").read_text()
        for content in (master, array):
            assert "#SBATCH --hold" in content
            assert "--exclusive=user" not in content


# ---------------------------------------------------------------------------
# TestHistoricalBinaryAttrsInChildAndArrayScripts
# ---------------------------------------------------------------------------


class TestHistoricalBinaryAttrsInChildAndArrayScripts:
    """Regression tests for the --fortran-commit child/array script bugs.

    Two bugs the cherry-picked commit introduced:

    * Generated bash scripts did a runtime ``cp "{executable}"`` that
      pointed at the historical-build temp dir, which the CLI removes
      as soon as ``sbatch`` returns.  Child/array tasks ran later and
      hit ``cp: cannot stat …``.
    * Generated ``FortranMicro.from_hdf5(...)`` / ``FortranMacro.from_hdf5(...)``
      calls didn't pass ``skip_binary_verification=True``, so the
      runner's ``_verify_binary_version`` raised
      ``StaleBinaryError`` on the intentional binary↔source mismatch.
    """

    _HIST = {
        "backend_commit": "a" * 40,
        "backend_dirty": False,
        "backend_historical": True,
    }

    # -- generate_micro_child_script (legacy single-child path) ------------

    def test_legacy_single_tier_no_runtime_cp_under_historical(self, tmp_path):
        """Single-tier legacy child must not cp from the (deleted) build dir."""
        script = generate_micro_child_script(
            tmp_path / "stage", "run-01", tmp_path / "run-01.h5",
            "/build/bin/micro_rates", "",
            historical_backend_attrs=self._HIST,
        )
        assert "/build/bin/micro_rates" not in script
        assert "cp " not in script

    def test_legacy_single_tier_emits_skip_verify_under_historical(self, tmp_path):
        script = generate_micro_child_script(
            tmp_path / "stage", "run-01", tmp_path / "run-01.h5",
            "/build/bin/micro_rates", "",
            historical_backend_attrs=self._HIST,
        )
        assert "skip_binary_verification=True" in script

    def test_legacy_two_tier_cp_uses_staging_under_historical(self, tmp_path):
        """Two-tier legacy child cp must source from staging_dir under --fortran-commit."""
        staging = tmp_path / "stage"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/build/bin/micro_rates", "",
            fast_tmp_root="/nvme/scratch",
            historical_backend_attrs=self._HIST,
        )
        assert "/build/bin/micro_rates" not in script
        assert f'cp "{staging}/micro_rates" "${{local_work_dir}}/"' in script

    def test_legacy_two_tier_emits_skip_verify_under_historical(self, tmp_path):
        script = generate_micro_child_script(
            tmp_path / "stage", "run-01", tmp_path / "run-01.h5",
            "/build/bin/micro_rates", "",
            fast_tmp_root="/nvme/scratch",
            historical_backend_attrs=self._HIST,
        )
        assert "skip_binary_verification=True" in script

    def test_no_skip_verify_when_historical_absent(self, tmp_path):
        """Non-historical case must not emit skip_binary_verification."""
        script = generate_micro_child_script(
            tmp_path / "stage", "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "",
        )
        assert "skip_binary_verification" not in script

    # -- generate_micro_array_script (array micro path) --------------------

    def test_micro_array_single_tier_emits_skip_verify_under_historical(self, tmp_path):
        script = generate_micro_array_script(
            tmp_path / "stage", "run-01", tmp_path / "run-01.h5",
            "/build/bin/micro_rates", num_children=4,
            historical_backend_attrs=self._HIST,
        )
        assert "skip_binary_verification=True" in script

    def test_micro_array_two_tier_cp_uses_staging_under_historical(self, tmp_path):
        staging = tmp_path / "stage"
        script = generate_micro_array_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/build/bin/micro_rates", num_children=4,
            fast_tmp_root="/nvme/scratch",
            historical_backend_attrs=self._HIST,
        )
        assert "/build/bin/micro_rates" not in script
        assert f'cp "{staging}/micro_rates" "${{local_work_dir}}/"' in script
        assert "skip_binary_verification=True" in script

    # -- generate_macro_array_script (array macro path) --------------------

    def test_macro_array_single_tier_emits_skip_verify_under_historical(self, tmp_path):
        script = generate_macro_array_script(
            tmp_path / "stage", "run-01", tmp_path / "run-01.h5",
            "/build/bin/macro.exe", n_sims=3,
            historical_backend_attrs=self._HIST,
        )
        assert "skip_binary_verification=True" in script

    def test_macro_array_two_tier_cp_uses_staging_under_historical(self, tmp_path):
        staging = tmp_path / "stage"
        script = generate_macro_array_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/build/bin/macro.exe", n_sims=3,
            fast_tmp_root="/nvme/scratch",
            historical_backend_attrs=self._HIST,
        )
        assert "/build/bin/macro.exe" not in script
        assert f'cp "{staging}/macro.exe" "${{local_work_dir}}/"' in script
        assert "skip_binary_verification=True" in script

    # -- two-tier cp source is correct even WITHOUT --fortran-commit -------

    def test_two_tier_cp_uses_staging_in_normal_mode_too(self, tmp_path):
        """Two-tier cp must always source from staging_dir (not just under historical).

        Regression: the original two-tier path cp'd from the caller's
        --executable, which is fine in normal mode (the binary's still
        there) but pointlessly leaks the caller's path into the script.
        """
        staging = tmp_path / "stage"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "",
            fast_tmp_root="/nvme/scratch",
        )
        assert f'cp "{staging}/micro.exe" "${{local_work_dir}}/"' in script
        assert "/bin/micro.exe" not in script


# ---------------------------------------------------------------------------
# TestSourceStampInGeneratedScripts
# ---------------------------------------------------------------------------


class TestSourceStampInGeneratedScripts:
    """The master Slurm job resolves the ``src/fortran/`` stamp on the
    submit host (where the repo is) and bakes it into the generated child
    and array task scripts as a ``source_stamp=(commit, dirty)`` kwarg to
    ``FortranMicro.from_hdf5`` / ``FortranMacro.from_hdf5``.  The compute
    node then compares the binary against that stamp directly instead of
    running ``git`` from the staging directory (which is outside the
    checkout and would fail).
    """

    _STAMP = ("abc123def456", "clean")

    # -- generate_micro_child_script (legacy single-child path) ------------

    def test_legacy_single_tier_bakes_source_stamp(self, tmp_path):
        script = generate_micro_child_script(
            tmp_path / "stage", "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "",
            source_stamp=self._STAMP,
        )
        assert "source_stamp=('abc123def456', 'clean')" in script

    def test_legacy_two_tier_bakes_source_stamp(self, tmp_path):
        script = generate_micro_child_script(
            tmp_path / "stage", "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "",
            fast_tmp_root="/nvme/scratch",
            source_stamp=self._STAMP,
        )
        assert "source_stamp=('abc123def456', 'clean')" in script

    def test_legacy_no_source_stamp_when_unset(self, tmp_path):
        script = generate_micro_child_script(
            tmp_path / "stage", "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "",
        )
        assert "source_stamp=" not in script

    def test_legacy_source_stamp_suppressed_under_historical(self, tmp_path):
        """Under --fortran-commit the staleness check is bypassed, so any
        passed source_stamp must not be baked into the script (it would be
        dead code at best, confusing at worst)."""
        script = generate_micro_child_script(
            tmp_path / "stage", "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "",
            historical_backend_attrs={"backend_commit": "h" * 40},
            source_stamp=self._STAMP,
        )
        assert "source_stamp=" not in script
        assert "skip_binary_verification=True" in script

    # -- generate_micro_array_script (array micro path) --------------------

    def test_micro_array_single_tier_bakes_source_stamp(self, tmp_path):
        from lysis.tools.slurm import generate_array_script, _micro_spec
        spec = _micro_spec(4, "")
        script = generate_array_script(
            spec, tmp_path / "stage", "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe",
            source_stamp=self._STAMP,
        )
        assert "source_stamp=('abc123def456', 'clean')" in script

    def test_micro_array_two_tier_bakes_source_stamp(self, tmp_path):
        from lysis.tools.slurm import generate_array_script, _micro_spec
        spec = _micro_spec(4, "")
        script = generate_array_script(
            spec, tmp_path / "stage", "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe",
            fast_tmp_root="/nvme/scratch",
            source_stamp=self._STAMP,
        )
        assert "source_stamp=('abc123def456', 'clean')" in script

    def test_array_source_stamp_suppressed_under_historical(self, tmp_path):
        from lysis.tools.slurm import generate_array_script, _micro_spec
        spec = _micro_spec(4, "")
        script = generate_array_script(
            spec, tmp_path / "stage", "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe",
            historical_backend_attrs={"backend_commit": "h" * 40},
            source_stamp=self._STAMP,
        )
        assert "source_stamp=" not in script
        assert "skip_binary_verification=True" in script


# ---------------------------------------------------------------------------
# TestSubmitWithHistoricalBinaryAttrsEndToEnd
# ---------------------------------------------------------------------------


class TestSubmitWithHistoricalBinaryAttrsEndToEnd:
    """End-to-end: submit_*_slurm_job threads historical_backend_attrs into
    every generated child/array script.
    """

    _HIST = {
        "backend_commit": "b" * 40,
        "backend_dirty": False,
        "backend_historical": True,
    }

    @pytest.fixture(autouse=True)
    def _patch_copy2(self):
        with patch("lysis.tools.slurm.shutil.copy2"):
            yield

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_legacy_single_child_propagates_to_child_script(
        self, mock_sbatch, micro_hdf5, tmp_path,
    ):
        """Bug reproducer: user's command was --num-children=0 --fortran-commit ..."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(
            micro_hdf5, "/build/bin/micro_rates",
            staging_root=staging_root,
            historical_backend_attrs=self._HIST,
        )
        staging_dir = list(staging_root.iterdir())[0]
        child = (staging_dir / "lysis-micro-child__run-01.sh").read_text()
        assert "skip_binary_verification=True" in child
        # The original --executable path must not appear (would imply
        # the script tries to cp from the now-deleted build dir).
        assert "/build/bin/micro_rates" not in child

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_array_path_propagates_to_array_script(
        self, mock_sbatch, micro_hdf5, tmp_path,
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        with patch.object(
            __import__("lysis.execution.fortran_micro", fromlist=["FortranMicro"]).FortranMicro,
            "_write_setup_files",
        ):
            submit_micro_slurm_job(
                micro_hdf5, "/build/bin/micro_rates",
                staging_root=staging_root,
                num_children=4,
                historical_backend_attrs=self._HIST,
            )
        staging_dir = list(staging_root.iterdir())[0]
        array = (staging_dir / "lysis-micro-array__run-01.sh").read_text()
        assert "skip_binary_verification=True" in array


# ---------------------------------------------------------------------------
# TestSourcePinning
# ---------------------------------------------------------------------------


class TestSnapshotLysisSrc:
    """Tests for :func:`_snapshot_lysis_src`."""

    def test_returns_python_src_directory(self, tmp_path):
        from lysis.tools.slurm import _repo_root
        staging = tmp_path / "staging"
        staging.mkdir()
        result = _snapshot_lysis_src(_repo_root(), staging)
        assert result == staging / "python_src"

    def test_creates_lysis_package_under_python_src(self, tmp_path):
        from lysis.tools.slurm import _repo_root
        staging = tmp_path / "staging"
        staging.mkdir()
        _snapshot_lysis_src(_repo_root(), staging)
        assert (staging / "python_src" / "lysis" / "__init__.py").exists()

    def test_copies_subpackages(self, tmp_path):
        """Subpackages such as tools/ and execution/ must be in the snapshot."""
        from lysis.tools.slurm import _repo_root
        staging = tmp_path / "staging"
        staging.mkdir()
        _snapshot_lysis_src(_repo_root(), staging)
        snap = staging / "python_src" / "lysis"
        assert (snap / "tools" / "slurm.py").exists()
        assert (snap / "execution" / "fortran_micro.py").exists()

    def test_excludes_pycache(self, tmp_path):
        """__pycache__ directories must be filtered out of the snapshot."""
        from lysis.tools.slurm import _repo_root
        # Ensure at least one __pycache__ exists in the source tree by importing.
        import lysis.tools.slurm  # noqa: F401
        staging = tmp_path / "staging"
        staging.mkdir()
        _snapshot_lysis_src(_repo_root(), staging)
        snap = staging / "python_src" / "lysis"
        pycache_dirs = list(snap.rglob("__pycache__"))
        assert pycache_dirs == []


class TestSourcePinningInGeneratedScripts:
    """PYTHONPATH source pinning must reach every generated script."""

    @pytest.fixture(autouse=True)
    def _patch_copy2(self):
        """The binary pre-stage copy doesn't need a real file."""
        with patch("lysis.tools.slurm.shutil.copy2"):
            yield

    def _staging_dir(self, staging_root):
        return list(staging_root.iterdir())[0]

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_legacy_micro_master_exports_pythonpath(
        self, mock_sbatch, micro_hdf5, tmp_path,
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(
            micro_hdf5, "/bin/micro.exe", staging_root=staging_root,
        )
        staging_dir = self._staging_dir(staging_root)
        master = (staging_dir / "lysis-micro-master__run-01.sh").read_text()
        assert f'export PYTHONPATH="{staging_dir}/python_src"' in master

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_legacy_micro_child_exports_pythonpath(
        self, mock_sbatch, micro_hdf5, tmp_path,
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(
            micro_hdf5, "/bin/micro.exe", staging_root=staging_root,
        )
        staging_dir = self._staging_dir(staging_root)
        child = (staging_dir / "lysis-micro-child__run-01.sh").read_text()
        assert f'export PYTHONPATH="{staging_dir}/python_src"' in child

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_legacy_micro_creates_snapshot_directory(
        self, mock_sbatch, micro_hdf5, tmp_path,
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(
            micro_hdf5, "/bin/micro.exe", staging_root=staging_root,
        )
        staging_dir = self._staging_dir(staging_root)
        assert (staging_dir / "python_src" / "lysis" / "__init__.py").exists()

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_array_micro_master_and_array_export_pythonpath(
        self, mock_sbatch, micro_hdf5, tmp_path,
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        submit_micro_slurm_job(
            micro_hdf5, "/bin/micro.exe",
            staging_root=staging_root,
            num_children=4,
        )
        staging_dir = self._staging_dir(staging_root)
        master = (staging_dir / "lysis-micro-master__run-01.sh").read_text()
        array = (staging_dir / "lysis-micro-array__run-01.sh").read_text()
        expected = f'export PYTHONPATH="{staging_dir}/python_src"'
        assert expected in master
        assert expected in array

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_macro_master_and_array_export_pythonpath(
        self, mock_sbatch, macro_hdf5, tmp_path,
    ):
        with patch(
            "lysis.execution.fortran_macro.FortranMacro._write_setup_files"
        ):
            staging_root = tmp_path / "staging_root"
            staging_root.mkdir()
            submit_macro_slurm_job(
                macro_hdf5, "/bin/macro.exe", staging_root=staging_root,
            )
        staging_dir = self._staging_dir(staging_root)
        master = (staging_dir / "lysis-macro-master__run-01.sh").read_text()
        array = (staging_dir / "lysis-macro-array__run-01.sh").read_text()
        expected = f'export PYTHONPATH="{staging_dir}/python_src"'
        assert expected in master
        assert expected in array

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_macro_creates_snapshot_directory(
        self, mock_sbatch, macro_hdf5, tmp_path,
    ):
        with patch(
            "lysis.execution.fortran_macro.FortranMacro._write_setup_files"
        ):
            staging_root = tmp_path / "staging_root"
            staging_root.mkdir()
            submit_macro_slurm_job(
                macro_hdf5, "/bin/macro.exe", staging_root=staging_root,
            )
        staging_dir = self._staging_dir(staging_root)
        assert (staging_dir / "python_src" / "lysis" / "__init__.py").exists()

    def test_standalone_generate_micro_child_omits_pythonpath_by_default(
        self, tmp_path,
    ):
        """The public generate_* helpers leave PYTHONPATH untouched unless asked."""
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "",
        )
        assert "PYTHONPATH" not in script

    def test_standalone_generate_micro_child_honours_pythonpath_arg(
        self, tmp_path,
    ):
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5",
            "/bin/micro.exe", "",
            pythonpath="/snap/python_src",
        )
        assert 'export PYTHONPATH="/snap/python_src"' in script


# ---------------------------------------------------------------------------
# Python backend (--backend python --slurm)
# ---------------------------------------------------------------------------


class TestPythonMacroSpec:
    """:func:`_python_macro_spec` is the contract between Python CLI dispatch
    and the shared array-mode internals."""

    def test_backend_is_python(self):
        assert _python_macro_spec().backend == "python"

    def test_single_task_array(self):
        """HDF5 single-writer constraint forces all sims in series → 1 task."""
        assert _python_macro_spec().num_array_tasks == 1

    def test_runner_class_is_python_macro(self):
        spec = _python_macro_spec()
        assert spec.runner_module == "lysis.execution.python_macro"
        assert spec.runner_class == "PythonMacro"

    def test_no_binary_or_codes(self):
        """No Fortran executable to pre-stage, no in/out file codes."""
        spec = _python_macro_spec()
        assert spec.prestage_executable is False
        assert spec.in_file_code is None
        assert spec.out_file_code == ""
        assert spec.num_children is None

    def test_no_concat_no_nfs_wait(self):
        """Task writes directly to HDF5 — nothing to concat, nothing to flush."""
        spec = _python_macro_spec()
        assert spec.needs_concat_step is False
        assert spec.nfs_wait_seconds == 0


class TestSubmitPythonMacroSlurmJob:
    """End-to-end script generation for the Python --slurm path.

    All tests monkeypatch ``gs.sbatch`` and the venv-sync side effect so
    nothing leaves the test tmpdir.  We then inspect the generated array
    and master scripts to assert the contract.
    """

    @pytest.fixture(autouse=True)
    def _patch_env_sync(self):
        """The venv sync is expensive and irrelevant here."""
        with patch("lysis.tools.slurm._ensure_env_synced"):
            yield

    def _staging_dir(self, staging_root):
        return list(staging_root.iterdir())[0]

    @patch("lysis.tools.slurm.gs.sbatch", return_value=8765)
    def test_returns_master_job_id(self, mock_sbatch, tmp_path):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        h5 = tmp_path / "py-run.h5"
        h5.write_text("")  # placeholder; PythonMacro is never instantiated here
        job_id = submit_python_macro_slurm_job(
            h5, staging_root=staging_root,
        )
        assert job_id == 8765

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_no_setup_data_dir_pre_staged(self, mock_sbatch, tmp_path):
        """Python backend has no setup files — staging dir must not contain
        a data/{run_code}/ skeleton that the Fortran path pre-stages."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        h5 = tmp_path / "py-run.h5"
        h5.write_text("")
        submit_python_macro_slurm_job(h5, staging_root=staging_root)
        staging_dir = self._staging_dir(staging_root)
        assert not (staging_dir / "data").exists()

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_array_script_calls_python_macro_execute(
        self, mock_sbatch, tmp_path,
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        h5 = tmp_path / "py-run.h5"
        h5.write_text("")
        submit_python_macro_slurm_job(h5, staging_root=staging_root)
        staging_dir = self._staging_dir(staging_root)
        array = (staging_dir / "lysis-macro-array__py-run.sh").read_text()
        assert "from lysis.execution.python_macro import PythonMacro" in array
        assert f"PythonMacro.from_hdf5('{h5.resolve()}').execute()" in array

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_array_script_has_array_zero_zero_directive(
        self, mock_sbatch, tmp_path,
    ):
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        h5 = tmp_path / "py-run.h5"
        h5.write_text("")
        submit_python_macro_slurm_job(h5, staging_root=staging_root)
        staging_dir = self._staging_dir(staging_root)
        array = (staging_dir / "lysis-macro-array__py-run.sh").read_text()
        assert "0-0" in array

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_array_script_does_not_reference_binary_or_index(
        self, mock_sbatch, tmp_path,
    ):
        """PythonMacro.from_hdf5 takes no binary, no index, no file codes."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        h5 = tmp_path / "py-run.h5"
        h5.write_text("")
        submit_python_macro_slurm_job(h5, staging_root=staging_root)
        staging_dir = self._staging_dir(staging_root)
        array = (staging_dir / "lysis-macro-array__py-run.sh").read_text()
        # The runner call must not index by task id (single 0-0 array runs
        # every sim in series).  The identifying header (#81) legitimately
        # prints ${SLURM_ARRAY_TASK_ID}, so target the Fortran runner-call
        # form specifically rather than the bare env-var name.
        assert "os.environ['SLURM_ARRAY_TASK_ID']" not in array
        assert "exec_in_workdir" not in array
        assert "out_file_code" not in array

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_single_tier_no_local_workdir_or_h5_copy(
        self, mock_sbatch, tmp_path,
    ):
        """Without --fast-tmp-root the task writes straight to the original h5."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        h5 = tmp_path / "py-run.h5"
        h5.write_text("")
        submit_python_macro_slurm_job(h5, staging_root=staging_root)
        staging_dir = self._staging_dir(staging_root)
        array = (staging_dir / "lysis-macro-array__py-run.sh").read_text()
        assert "local_work_dir" not in array
        assert "local_h5_path" not in array
        assert "mktemp" not in array

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_two_tier_copies_h5_to_fast_tmp_and_back(
        self, mock_sbatch, tmp_path,
    ):
        """--fast-tmp-root must cp the .h5 in before running and back after."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        h5 = tmp_path / "py-run.h5"
        h5.write_text("")
        submit_python_macro_slurm_job(
            h5, staging_root=staging_root, fast_tmp_root="/nvme/scratch",
        )
        staging_dir = self._staging_dir(staging_root)
        array = (staging_dir / "lysis-macro-array__py-run.sh").read_text()
        h5_resolved = str(h5.resolve())
        assert 'mktemp -d -p "/nvme/scratch"' in array
        assert 'local_h5_path="${local_work_dir}/py-run.h5"' in array
        # cp out (h5 → fast tmp)
        assert f'cp "{h5_resolved}" "${{local_h5_path}}"' in array
        # cp back (fast tmp → h5)
        assert f'cp "${{local_h5_path}}" "{h5_resolved}"' in array
        # final cleanup of the local work dir
        assert 'rm -rf "${local_work_dir}"' in array

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_two_tier_array_script_uses_local_h5_in_python_call(
        self, mock_sbatch, tmp_path,
    ):
        """The python -c must run PythonMacro against the local h5 copy."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        h5 = tmp_path / "py-run.h5"
        h5.write_text("")
        submit_python_macro_slurm_job(
            h5, staging_root=staging_root, fast_tmp_root="/nvme/scratch",
        )
        staging_dir = self._staging_dir(staging_root)
        array = (staging_dir / "lysis-macro-array__py-run.sh").read_text()
        assert "PythonMacro.from_hdf5('${local_h5_path}').execute()" in array

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_script_does_not_import_results(
        self, mock_sbatch, tmp_path,
    ):
        """Task wrote to HDF5 itself — master must not try to import per-sim
        outputs (no Fortran-style data_dir / runner construction)."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        h5 = tmp_path / "py-run.h5"
        h5.write_text("")
        submit_python_macro_slurm_job(h5, staging_root=staging_root)
        staging_dir = self._staging_dir(staging_root)
        master = (staging_dir / "lysis-macro-master__py-run.py").read_text()
        assert "import_results" not in master
        assert "concatenate_child_outputs" not in master
        assert "BINARY_NAME" not in master
        assert "Python task wrote results directly into HDF5" in master

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_master_script_still_polls_squeue(
        self, mock_sbatch, tmp_path,
    ):
        """Single-task array still needs the squeue poll — the master is
        the only handle on task completion before staging cleanup."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        h5 = tmp_path / "py-run.h5"
        h5.write_text("")
        submit_python_macro_slurm_job(h5, staging_root=staging_root)
        staging_dir = self._staging_dir(staging_root)
        master = (staging_dir / "lysis-macro-master__py-run.py").read_text()
        assert "gs.squeue.read()" in master
        assert "ARRAY_JOB_ID" in master

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_python_backend_exports_pythonpath_to_snapshot(
        self, mock_sbatch, tmp_path,
    ):
        """Source pinning applies to the Python backend too — same snapshot
        directory, same PYTHONPATH export, immune to source edits while queued."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        h5 = tmp_path / "py-run.h5"
        h5.write_text("")
        submit_python_macro_slurm_job(h5, staging_root=staging_root)
        staging_dir = self._staging_dir(staging_root)
        array = (staging_dir / "lysis-macro-array__py-run.sh").read_text()
        master = (staging_dir / "lysis-macro-master__py-run.sh").read_text()
        expected = f'export PYTHONPATH="{staging_dir}/python_src"'
        assert expected in array
        assert expected in master
        assert (staging_dir / "python_src" / "lysis" / "__init__.py").exists()

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_partition_and_sbatch_overrides_apply(
        self, mock_sbatch, tmp_path,
    ):
        """--partition and --sbatch must reach both array and master headers."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        h5 = tmp_path / "py-run.h5"
        h5.write_text("")
        submit_python_macro_slurm_job(
            h5,
            staging_root=staging_root,
            partition="long",
            sbatch_overrides={"hold": "", "exclusive=user": None},
        )
        staging_dir = self._staging_dir(staging_root)
        array = (staging_dir / "lysis-macro-array__py-run.sh").read_text()
        master = (staging_dir / "lysis-macro-master__py-run.sh").read_text()
        for content in (array, master):
            assert "long" in content
            assert "#SBATCH --hold" in content
            assert "--exclusive=user" not in content

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_python_uses_uv_run_prefix(self, mock_sbatch, tmp_path):
        """Generated python invocations must go through ``uv run --frozen
        python`` so the project's .venv is on sys.path."""
        staging_root = tmp_path / "staging_root"
        staging_root.mkdir()
        h5 = tmp_path / "py-run.h5"
        h5.write_text("")
        submit_python_macro_slurm_job(h5, staging_root=staging_root)
        staging_dir = self._staging_dir(staging_root)
        array = (staging_dir / "lysis-macro-array__py-run.sh").read_text()
        assert "uv run" in array
        assert "--frozen" in array


# ---------------------------------------------------------------------------
# Worker log relocation (#84) + identifying header & shell tracing (#81)
# ---------------------------------------------------------------------------


def _sbatch_value_line(script: str, flag: str) -> str:
    """Return the ``#SBATCH {flag} ...`` line from a generated script."""
    for line in script.splitlines():
        s = line.strip()
        if s.startswith(f"#SBATCH {flag} "):
            return s
    raise AssertionError(f"no '#SBATCH {flag}' directive in script")


class TestWorkerLogOutputAndHeader:
    """#84: worker ``--output`` → staging dir. #81: header + ``set -x`` tracing."""

    # ---- #84: worker --output lands in the staging dir, NOT in .slurm -----

    def test_child_output_in_staging_not_slurm(self, tmp_path):
        staging = tmp_path / "staging"
        script = generate_micro_child_script(
            staging, "run-01", tmp_path / "run-01.h5", "/bin/micro.exe", "",
        )
        out_line = _sbatch_value_line(script, "--out")
        assert str(staging / "lysis-micro-child__run-01.out") in out_line
        assert ".slurm" not in out_line

    def test_micro_array_output_in_staging_not_slurm(self, tmp_path):
        staging = tmp_path / "staging"
        script = generate_micro_array_script(
            staging, "run-01", tmp_path / "run-01.h5", "/bin/micro.exe", 5,
        )
        out_line = _sbatch_value_line(script, "--out")
        assert str(staging / "lysis-micro-array__run-01__%a.out") in out_line
        assert ".slurm" not in out_line

    def test_macro_array_output_in_staging_not_slurm(self, tmp_path):
        staging = tmp_path / "staging"
        script = generate_macro_array_script(
            staging, "run-01", tmp_path / "run-01.h5", "/bin/macro.exe", n_sims=5,
        )
        out_line = _sbatch_value_line(script, "--out")
        assert str(staging / "lysis-macro-array__run-01__%a.out") in out_line
        assert ".slurm" not in out_line

    # ---- #81 comment: array --job-name carries no literal %a --------------

    def test_micro_array_job_name_has_no_literal_percent_a(self, tmp_path):
        script = generate_micro_array_script(
            tmp_path / "s", "run-01", tmp_path / "run-01.h5", "/bin/micro.exe", 5,
        )
        name_line = _sbatch_value_line(script, "--job-name")
        assert name_line.endswith("lysis-micro-array__run-01")
        assert "%a" not in name_line

    def test_macro_array_job_name_has_no_literal_percent_a(self, tmp_path):
        script = generate_macro_array_script(
            tmp_path / "s", "run-01", tmp_path / "run-01.h5", "/bin/macro.exe",
            n_sims=5,
        )
        name_line = _sbatch_value_line(script, "--job-name")
        assert name_line.endswith("lysis-macro-array__run-01")
        assert "%a" not in name_line

    # ---- #81: identifying header content ----------------------------------

    def test_header_present_with_baked_and_runtime_fields(self, tmp_path):
        script = generate_micro_array_script(
            tmp_path / "s", "run-01", tmp_path / "run-01.h5", "/bin/micro.exe", 5,
        )
        assert "lysis slurm task" in script
        # Baked submit-time fields (emitted via printf, hence not line-anchored).
        assert re.search(r"run\s+: run-01", script)
        assert re.search(r"scale\s+: micro", script)
        assert re.search(r"role\s+: array task", script)
        # Runtime ${SLURM_*} expansions (evaluated on the compute node).
        assert "${SLURM_JOB_ID" in script
        assert "${SLURM_ARRAY_TASK_ID" in script

    def test_child_header_role_is_micro_child(self, tmp_path):
        script = generate_micro_child_script(
            tmp_path / "s", "run-01", tmp_path / "run-01.h5", "/bin/micro.exe", "",
        )
        assert re.search(r"role\s+: micro child", script)

    def test_experiment_field_default_unknown(self, tmp_path):
        script = generate_micro_child_script(
            tmp_path / "s", "run-01", tmp_path / "run-01.h5", "/bin/micro.exe", "",
        )
        assert re.search(r"experiment\s+: \(unknown\)", script)

    def test_experiment_field_threaded(self, tmp_path):
        script = generate_micro_child_script(
            tmp_path / "s", "run-01", tmp_path / "run-01.h5", "/bin/micro.exe", "",
            experiment="EXP-42",
        )
        assert re.search(r"experiment\s+: EXP-42", script)

    def test_experiment_value_is_shell_quoted_against_injection(self, tmp_path):
        """A caller-controlled experiment with $(...) must be baked verbatim.

        The header's submit-time fields are emitted via ``printf`` with each
        line shell-quoted, so bash never performs command/parameter
        substitution on them.  The malicious sequence must appear inside a
        single-quoted printf argument, where ``$(...)`` is inert.
        """
        script = generate_micro_child_script(
            tmp_path / "s", "run-01", tmp_path / "run-01.h5", "/bin/micro.exe", "",
            experiment="EXP-$(touch /tmp/pwned)",
        )
        # The whole header line is a single-quoted printf literal — $() inert.
        assert "'experiment  : EXP-$(touch /tmp/pwned)'" in script
        # And it is NOT emitted via the unquoted runtime heredoc.
        assert "cat <<HEADER" not in script

    # ---- #81: strict mode + tracing, enabled AFTER the header -------------

    def test_strict_mode_and_trace_enabled_after_header(self, tmp_path):
        script = generate_micro_child_script(
            tmp_path / "s", "run-01", tmp_path / "run-01.h5", "/bin/micro.exe", "",
        )
        assert "set -euo pipefail" in script
        assert "\nset -x" in script
        header_pos = script.index("lysis slurm task")
        # Strict mode (incl. nounset) is enabled only after the header prints,
        # so the header's optional ${SLURM_*} reads can't trip ``set -u``.
        assert header_pos < script.index("set -euo pipefail")
        assert header_pos < script.index("\nset -x")

    def test_two_tier_setup_runs_under_trace(self, tmp_path):
        """SIM/cp/mv setup must come AFTER ``set -x`` so it is traced."""
        script = generate_macro_array_script(
            tmp_path / "s", "run-01", tmp_path / "run-01.h5", "/bin/macro.exe",
            n_sims=5, fast_tmp_root="/scratch",
        )
        assert script.index("\nset -x") < script.index("SIM=")

    # ---- B: two-tier fast-scratch is reclaimed on any exit (set -e safe) ---

    @pytest.mark.parametrize(
        "make_script",
        [
            lambda p: generate_micro_array_script(
                p / "s", "run-01", p / "run-01.h5", "/bin/micro.exe", 5,
                fast_tmp_root="/scratch",
            ),
            lambda p: generate_macro_array_script(
                p / "s", "run-01", p / "run-01.h5", "/bin/macro.exe",
                n_sims=5, fast_tmp_root="/scratch",
            ),
            lambda p: generate_micro_child_script(
                p / "s", "run-01", p / "run-01.h5", "/bin/micro.exe", "",
                fast_tmp_root="/scratch",
            ),
        ],
        ids=["micro-array", "macro-array", "micro-child"],
    )
    def test_two_tier_registers_cleanup_trap_before_mv(self, tmp_path, make_script):
        script = make_script(tmp_path)
        assert 'trap \'rm -rf "${local_work_dir}"\' EXIT' in script
        # Trap must be registered before the risky mv so a set -e abort there
        # still reclaims the fast-tmp dir (no permanent node-local leak).
        assert script.index("trap ") < script.index("\nmv ")
        # The old unconditional post-mv ``rm -rf`` is gone (trap owns cleanup);
        # only the trap references rm -rf of the work dir.
        assert script.count('rm -rf "${local_work_dir}"') == 1

    def test_two_tier_setup_cp_glob_is_guarded(self, tmp_path):
        """The staged-files copy must not use a bare glob under ``set -e``.

        A bare ``cp "${staging_datadir}/"*`` aborts on an unexpanded ``*`` when
        the staging dir is empty.  The guarded form collects matches via
        nullglob and fails loudly with a clear message instead.
        """
        script = generate_macro_array_script(
            tmp_path / "s", "run-01", tmp_path / "run-01.h5", "/bin/macro.exe",
            n_sims=5, fast_tmp_root="/scratch",
        )
        # No bare glob copy remains.
        assert 'cp "${staging_datadir}/"*' not in script
        # Guarded form: nullglob collection + emptiness check + array copy.
        assert "shopt -s nullglob" in script
        assert 'staged_files=("${staging_datadir}"/*)' in script
        assert 'cp "${staged_files[@]}" "${local_datadir}/"' in script
        # The emptiness guard precedes the cp (fail loud, don't skip setup).
        assert script.index('${#staged_files[@]}') < script.index(
            'cp "${staged_files[@]}"'
        )

    @pytest.mark.parametrize(
        "make_script",
        [
            lambda p: generate_micro_array_script(
                p / "s", "run-01", p / "run-01.h5", "/bin/micro.exe", 5,
                fast_tmp_root="/scratch",
            ),
            lambda p: generate_micro_child_script(
                p / "s", "run-01", p / "run-01.h5", "/bin/micro.exe", "",
                fast_tmp_root="/scratch",
            ),
        ],
        ids=["micro-array", "micro-child"],
    )
    def test_micro_two_tier_output_mv_glob_is_guarded(self, tmp_path, make_script):
        """Micro paths move flat ``__NN`` files via a glob — it must be guarded.

        Unlike the macro per-sim ``mv`` (a single named subdir), the micro
        output move is a ``${local_datadir}/*`` glob with the same unexpanded-
        ``*`` hazard under ``set -e`` as the staged-files cp.
        """
        script = make_script(tmp_path)
        # No bare glob mv remains.
        assert 'mv "${local_datadir}"/*' not in script
        # Guarded form: nullglob collection + emptiness check + array mv.
        assert 'output_files=("${local_datadir}"/*)' in script
        assert 'mv "${output_files[@]}" "${staging_datadir}/"' in script
        assert script.index('${#output_files[@]}') < script.index(
            'mv "${output_files[@]}"'
        )

    def test_macro_two_tier_mv_is_named_path_not_glob(self, tmp_path):
        """The macro per-sim move is a named subdir — no glob, no guard needed."""
        script = generate_macro_array_script(
            tmp_path / "s", "run-01", tmp_path / "run-01.h5", "/bin/macro.exe",
            n_sims=5, fast_tmp_root="/scratch",
        )
        assert 'mv "${local_datadir}/${SIM}" "${staging_datadir}/"' in script
        assert 'mv "${local_datadir}"/*' not in script
        assert "output_files=" not in script


class TestMasterLogStaysInSlurm:
    """#84 scope: master ``--output`` stays in ``.slurm`` (worker-only move)."""

    @pytest.fixture(autouse=True)
    def _patch_copy2(self):
        with patch("lysis.tools.slurm.shutil.copy2"):
            yield

    def _staging_dir(self, root):
        return next(
            p for p in root.iterdir() if p.is_dir() and "lysis-micro" in p.name
        )

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_legacy_master_output_in_slurm_with_header(
        self, mock_sbatch, micro_hdf5, tmp_path
    ):
        submit_micro_slurm_job(micro_hdf5, "/bin/micro.exe", staging_root=tmp_path)
        master = (
            self._staging_dir(tmp_path) / "lysis-micro-master__run-01.sh"
        ).read_text()
        out_line = _sbatch_value_line(master, "--out")
        assert ".slurm" in out_line
        # The master also gets the identifying header + tracing.
        assert "lysis slurm task" in master
        assert re.search(r"role\s+: master", master)
        assert "set -x" in master

    @patch("lysis.tools.slurm.gs.sbatch", return_value=1)
    def test_array_master_output_in_slurm_with_header(
        self, mock_sbatch, micro_hdf5, tmp_path
    ):
        with patch(
            "lysis.execution.fortran_micro.FortranMicro._write_setup_files"
        ):
            submit_micro_slurm_job(
                micro_hdf5, "/bin/micro.exe",
                staging_root=tmp_path, num_children=5,
            )
        master = (
            self._staging_dir(tmp_path) / "lysis-micro-master__run-01.sh"
        ).read_text()
        out_line = _sbatch_value_line(master, "--out")
        assert ".slurm" in out_line
        assert re.search(r"role\s+: master", master)
