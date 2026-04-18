"""Slurm dispatch helpers for the Lysis microscale workflow.

This module provides utilities for submitting and monitoring the Fortran
microscale simulation via the Slurm workload manager, using
`GooseSLURM <https://gooseslurm.readthedocs.io/>`_ for script generation
and job submission.

Workflow overview
-----------------
:func:`submit_micro_slurm_job` is the high-level entry point:

1. Creates a unique shared *staging directory* alongside the HDF5 file
   (never in ``/tmp``, which is a ramdisk on many HPC nodes).
2. Generates one or more *child* bash scripts (one per seed-split index,
   currently always 1) and writes them to the staging directory.
3. Generates a *master* Python script and a thin bash wrapper, both written
   to the staging directory.
4. Submits the bash wrapper as a Slurm job.  The master job runs on a
   compute node and orchestrates the children:

   a. Submits child script(s) via ``sbatch``.
   b. Polls ``squeue`` until all children finish.
   c. Calls :meth:`~lysis.execution.fortran_micro.FortranMicro.import_results`
      to convert Fortran output into the HDF5 file.
   d. Optionally removes the staging directory.

Two-tier storage (opt-in)
~~~~~~~~~~~~~~~~~~~~~~~~~
When *fast_tmp_root* is provided each child writes Fortran output to fast
node-local scratch (``mktemp -d -p <fast_tmp_root>``) rather than the
shared staging directory, then moves the results to staging before exiting.
The master job then reads from staging as usual.  This optimises I/O on
clusters with NVMe scratch on compute nodes.

.. important::

   All temporary directories are created with an explicit *dir* argument
   (never defaulting to ``/tmp``) because most HPC nodes use a ramdisk
   for ``/tmp`` that is too small for simulation data.
"""

import shutil
import tempfile
import time
from pathlib import Path
from string import Template
from typing import List, Optional

import GooseSLURM as gs


__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2026, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#: Template for the micro master Python script.  Substitution keys use ``$name``
#: syntax (not ``{name}``) so Python braces in the template body are left
#: untouched.
_MICRO_MASTER_PY_TEMPLATE = Template('''\
#!/usr/bin/env python3
"""Master Slurm job: submit child Fortran jobs, poll for completion, import results."""
import time
import shutil
from pathlib import Path

import GooseSLURM as gs
from lysis.config.run import Run
from lysis.execution.fortran_micro import FortranMicro

STAGING_DIR = Path($staging_dir)
HDF5_PATH = Path($hdf5_path)
RUN_CODE = $run_code
FILE_CODE = $file_code
KEEP_TMPDIR = $keep_tmpdir

# ---------------------------------------------------------------------------
# Submit child jobs
# ---------------------------------------------------------------------------
child_ids = []
for script in sorted(STAGING_DIR.glob("child_*.sh")):
    job_id = gs.sbatch([str(script)])
    child_ids.append(job_id)
    print(f"Submitted child job {job_id}: {script.name}", flush=True)

if not child_ids:
    raise RuntimeError(f"No child scripts found in {STAGING_DIR}")

# ---------------------------------------------------------------------------
# Poll until all children complete
# ---------------------------------------------------------------------------
remaining = list(child_ids)
while remaining:
    time.sleep(30)
    squeue_rows = gs.squeue.read()
    running_ids = {row["JOBID"] for row in squeue_rows}
    for row in squeue_rows:
        if row["JOBID"] in {str(j) for j in remaining}:
            state = row.get("STATE", "").upper()
            if state in ("FAILED", "CANCELLED"):
                raise RuntimeError(
                    f"Child job {row['JOBID']} entered state {state}"
                )
    remaining = [j for j in remaining if str(j) in running_ids]

print("All child jobs complete.", flush=True)

# ---------------------------------------------------------------------------
# Import results into HDF5
# ---------------------------------------------------------------------------
data_dir = STAGING_DIR / "data" / RUN_CODE
run = Run(str(HDF5_PATH.parent), run_code=RUN_CODE)
fm = FortranMicro(run=run, out_file_code=FILE_CODE)
fm.import_results(
    data_dir,
    keep_on_failure=True,
    keep_tmpdir=KEEP_TMPDIR,
)
print("Results imported successfully.", flush=True)

# ---------------------------------------------------------------------------
# Clean up staging directory
# ---------------------------------------------------------------------------
if not KEEP_TMPDIR:
    shutil.rmtree(STAGING_DIR, ignore_errors=True)
    print("Staging directory removed.", flush=True)

print("Master job complete.", flush=True)
''')


def _generate_micro_master_py(
    staging_dir: Path,
    hdf5_path: Path,
    run_code: str,
    file_code: str,
    keep_tmpdir: bool,
) -> str:
    """Return the content of the master Python script with values substituted.

    :param staging_dir: Absolute path to the shared staging directory.
    :param hdf5_path: Absolute path to the target ``.h5`` file.
    :param run_code: Run identifier (HDF5 file stem).
    :param file_code: Output file code suffix for the Fortran binary.
    :param keep_tmpdir: Whether to preserve the staging directory on success.
    :return: Python script text.
    :rtype: str
    """
    return _MICRO_MASTER_PY_TEMPLATE.substitute(
        staging_dir=repr(str(staging_dir)),
        hdf5_path=repr(str(hdf5_path)),
        run_code=repr(run_code),
        file_code=repr(file_code),
        keep_tmpdir=repr(keep_tmpdir),
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_micro_child_script(
    staging_dir: "Path | str",
    run_code: str,
    hdf5_path: "Path | str",
    executable: "Path | str",
    out_code: str = "",
    *,
    partition: Optional[str] = None,
    fast_tmp_root: Optional[str] = None,
    slurm_log_dir: Optional["Path | str"] = None,
) -> str:
    """Generate a Slurm bash script for a single child microscale Fortran job.

    The child script sources the user's shell environment, copies the
    executable into the working directory, and calls
    :meth:`~lysis.execution.fortran_micro.FortranMicro.exec_in_workdir` via
    ``python -c``.

    **Single-tier** (default, ``fast_tmp_root=None``): Fortran writes its
    output directly to ``{staging_dir}/data/{run_code}/``.

    **Two-tier** (``fast_tmp_root`` provided): Fortran writes to a private
    ``mktemp -d -p {fast_tmp_root}`` directory on the compute node, then the
    script moves the results to ``{staging_dir}/data/{run_code}/`` before
    exiting.

    :param staging_dir: Path to the shared staging directory (must already
        exist when the script runs).
    :type staging_dir: Path or str
    :param run_code: Run identifier (used to name the output subdirectory).
    :type run_code: str
    :param hdf5_path: Full path to the run's ``.h5`` file (passed to
        :meth:`~lysis.execution.fortran_micro.FortranMicro.from_hdf5`).
    :type hdf5_path: Path or str
    :param executable: Path to the compiled Fortran microscale binary.
    :type executable: Path or str
    :param out_code: Output file code suffix, defaults to ``""``.
    :type out_code: str, optional
    :param partition: Slurm partition for ``#SBATCH --partition``, defaults to
        ``None`` (no partition directive).
    :type partition: str, optional
    :param fast_tmp_root: Root directory for fast node-local scratch storage.
        When ``None`` (default) single-tier mode is used.
    :type fast_tmp_root: str, optional
    :param slurm_log_dir: Directory for Slurm ``.out`` logs.  Defaults to
        ``hdf5_path.parent / ".slurm"``.
    :type slurm_log_dir: Path or str, optional
    :return: Slurm bash script text.
    :rtype: str
    """
    staging_dir = Path(staging_dir)
    hdf5_path = Path(hdf5_path)
    executable = Path(executable)
    binary_name = executable.name
    if slurm_log_dir is None:
        slurm_log_dir = hdf5_path.parent / ".slurm"
    slurm_log_dir = Path(slurm_log_dir)

    sbatch_opts = {
        "out": str(slurm_log_dir / "child.slurm-%j.out"),
        "nodes": 1,
        "mem": 3096,
        "ntasks": 1,
        "cpus-per-task": 1,
        "exclusive=user": "",
    }
    if partition:
        sbatch_opts["partition"] = partition

    if fast_tmp_root is None:
        # ------------------------------------------------------------------
        # Single-tier: Fortran writes directly to the shared staging dir
        # ------------------------------------------------------------------
        setup = f"""\
# Setup — create output directory and copy binary to staging dir
staging_datadir="{staging_dir}/data/{run_code}"
mkdir -p "${{staging_datadir}}"
cp "{executable}" "{staging_dir}/" """

        execute = f"""\
# Execute Fortran binary via lysis
source ~/.bashrc && source ~/lysis.sh
python -c "
from pathlib import Path
from lysis.execution.fortran_micro import FortranMicro
fm = FortranMicro.from_hdf5('{hdf5_path}', '{staging_dir}/{binary_name}', out_file_code='{out_code}')
fm.exec_in_workdir(Path('{staging_dir}'))
" """

        sections = [setup, execute]

    else:
        # ------------------------------------------------------------------
        # Two-tier: Fortran → fast local dir, then mv to shared staging dir
        # ------------------------------------------------------------------
        setup = f"""\
# Setup — create local fast dir and shared staging dir
local_work_dir=$(mktemp -d -p "{fast_tmp_root}")
local_datadir="${{local_work_dir}}/data/{run_code}"
staging_datadir="{staging_dir}/data/{run_code}"
mkdir -p "${{local_datadir}}"
mkdir -p "${{staging_datadir}}"
cp "{executable}" "${{local_work_dir}}/" """

        execute = f"""\
# Execute Fortran binary into local fast storage
source ~/.bashrc && source ~/lysis.sh
python -c "
from pathlib import Path
from lysis.execution.fortran_micro import FortranMicro
fm = FortranMicro.from_hdf5('{hdf5_path}', '${{local_work_dir}}/{binary_name}', out_file_code='{out_code}')
fm.exec_in_workdir(Path('${{local_work_dir}}'))
" """

        move = """\
# Move output to shared staging, then remove local fast dir
mv "${local_datadir}"/* "${staging_datadir}/"
rm -rf "${local_work_dir}" """

        sections = [setup, execute, move]

    return gs.scripts.plain(sections, **sbatch_opts)


def submit_micro_child_job(script_text: str, script_path: "Path | str") -> int:
    """Write *script_text* to *script_path* and submit it to Slurm.

    :param script_text: Slurm bash script content (e.g. from
        :func:`generate_micro_child_script`).
    :type script_text: str
    :param script_path: File path to write the script to.
    :type script_path: Path or str
    :return: Slurm job ID.
    :rtype: int
    :raises subprocess.CalledProcessError: If ``sbatch`` fails.
    :raises ValueError: If ``sbatch`` output cannot be parsed for a job ID.
    """
    script_path = Path(script_path)
    script_path.write_text(script_text)
    script_path.chmod(0o755)
    return gs.sbatch([str(script_path)])


def wait_for_jobs(job_ids: List[int], poll_interval: int = 30) -> None:
    """Block until all *job_ids* have left the Slurm queue.

    Polls ``squeue`` every *poll_interval* seconds.  A job is considered
    done when its ID is absent from the squeue output.

    :param job_ids: Slurm job IDs to wait for.
    :type job_ids: list[int]
    :param poll_interval: Seconds between squeue polls, defaults to ``30``.
    :type poll_interval: int, optional
    :raises RuntimeError: If any job appears in squeue with state
        ``FAILED`` or ``CANCELLED``.
    """
    remaining = list(job_ids)
    while remaining:
        time.sleep(poll_interval)
        squeue_rows = gs.squeue.read()
        running_ids = {row["JOBID"] for row in squeue_rows}
        for row in squeue_rows:
            if row["JOBID"] in {str(j) for j in remaining}:
                state = row.get("STATE", "").upper()
                if state in ("FAILED", "CANCELLED"):
                    raise RuntimeError(
                        f"Job {row['JOBID']} entered state {state}"
                    )
        remaining = [j for j in remaining if str(j) in running_ids]


def submit_micro_slurm_job(
    hdf5_path: "Path | str",
    executable: "Path | str",
    *,
    staging_root: Optional["Path | str"] = None,
    partition: Optional[str] = None,
    fast_tmp_root: Optional[str] = None,
    keep_tmpdir: bool = False,
    out_code: str = "",
) -> int:
    """Submit the full HDF5-integrated microscale workflow as a Slurm master job.

    Creates a unique staging directory, writes child and master scripts into
    it, and submits a master Slurm job.  The master job (running on a compute
    node) orchestrates the child Fortran jobs, imports results into the HDF5
    file, and optionally cleans up.

    Directory naming
    ~~~~~~~~~~~~~~~~
    The staging directory is created with::

        tempfile.mkdtemp(
            prefix=f"lysis-micro-{run_code}-",
            dir=staging_root or hdf5_path.parent,
        )

    The ``dir`` argument is **always explicit** — ``/tmp`` is never used.

    :param hdf5_path: Full path to the run's ``.h5`` file.  Must contain
        ``micro_params``.
    :type hdf5_path: Path or str
    :param executable: Path to the compiled Fortran microscale binary.
    :type executable: Path or str
    :param staging_root: Root directory under which the staging temp dir is
        created.  Defaults to ``hdf5_path.parent``.
    :type staging_root: Path or str, optional
    :param partition: Slurm partition for both master and child jobs.
    :type partition: str, optional
    :param fast_tmp_root: Root directory for fast node-local scratch (opt-in
        two-tier storage).  When ``None`` (default) the child writes directly
        to the shared staging directory.
    :type fast_tmp_root: str, optional
    :param keep_tmpdir: If ``True``, preserve the staging directory after the
        master job completes (useful for debugging).
    :type keep_tmpdir: bool, optional
    :param out_code: Output file code suffix for the Fortran binary,
        defaults to ``""``.
    :type out_code: str, optional
    :return: Master Slurm job ID.
    :rtype: int
    :raises subprocess.CalledProcessError: If ``sbatch`` fails.
    """
    from lysis.execution.fortran_micro import FortranMicro  # noqa: PLC0415

    hdf5_path = Path(hdf5_path).resolve()
    executable = Path(executable).resolve()
    run_code = hdf5_path.stem

    # Create unique staging directory — never use /tmp
    staging_root_dir = Path(staging_root) if staging_root is not None else hdf5_path.parent
    staging_dir = Path(tempfile.mkdtemp(
        prefix=f"lysis-micro-{run_code}-",
        dir=str(staging_root_dir),
    ))

    # Slurm log directory (sibling to the HDF5 file)
    slurm_log_dir = hdf5_path.parent / ".slurm"
    slurm_log_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Pre-stage setup files (params.json)
    # ------------------------------------------------------------------
    staging_data_dir = staging_dir / "data" / run_code
    staging_data_dir.mkdir(parents=True, exist_ok=True)
    fm_setup = FortranMicro.from_hdf5(
        hdf5_path, str(executable), out_file_code=out_code,
    )
    fm_setup._write_setup_files(staging_data_dir)

    # ------------------------------------------------------------------
    # Write child script
    # ------------------------------------------------------------------
    child_script = generate_micro_child_script(
        staging_dir, run_code, hdf5_path, executable, out_code,
        partition=partition, fast_tmp_root=fast_tmp_root,
        slurm_log_dir=slurm_log_dir,
    )
    child_path = staging_dir / "child_000.sh"
    child_path.write_text(child_script)
    child_path.chmod(0o755)

    # ------------------------------------------------------------------
    # Write master.py
    # ------------------------------------------------------------------
    master_py_content = _generate_micro_master_py(
        staging_dir, hdf5_path, run_code, out_code, keep_tmpdir,
    )
    master_py_path = staging_dir / "master.py"
    master_py_path.write_text(master_py_content)

    # ------------------------------------------------------------------
    # Write master bash wrapper and submit
    # ------------------------------------------------------------------
    sbatch_opts = {
        "job-name": f"lysis-micro-{run_code}",
        "out": str(slurm_log_dir / "master.slurm-%j.out"),
        "nodes": 1,
        "mem": 3096,
        "ntasks": 1,
        "cpus-per-task": 1,
        "exclusive=user": "",
    }
    if partition:
        sbatch_opts["partition"] = partition

    master_sh_content = gs.scripts.plain(
        [
            "source ~/.bashrc && source ~/lysis.sh",
            f'python "{master_py_path}"',
        ],
        **sbatch_opts,
    )
    master_sh_path = staging_dir / "master.sh"
    master_sh_path.write_text(master_sh_content)
    master_sh_path.chmod(0o755)

    return gs.sbatch([str(master_sh_path)])


# ---------------------------------------------------------------------------
# Macroscale Slurm support
# ---------------------------------------------------------------------------

#: Template for the macro master Python script.  Substitution keys use
#: ``$name`` syntax so Python braces in the template body are untouched.
_MACRO_MASTER_PY_TEMPLATE = Template('''\
#!/usr/bin/env python3
"""Master Slurm job: submit macro array job, poll for completion, import results."""
import time
import shutil
from pathlib import Path

import GooseSLURM as gs
from lysis.config.run import Run
from lysis.execution.fortran_macro import FortranMacro

STAGING_DIR = Path($staging_dir)
HDF5_PATH = Path($hdf5_path)
RUN_CODE = $run_code
OUT_FILE_CODE = $out_file_code
KEEP_TMPDIR = $keep_tmpdir

# ---------------------------------------------------------------------------
# Submit array job
# ---------------------------------------------------------------------------
array_script = STAGING_DIR / "array.sh"
array_job_id = gs.sbatch([str(array_script)])
print(f"Submitted array job {array_job_id}: {array_script.name}", flush=True)

# ---------------------------------------------------------------------------
# Poll until all array tasks complete
# ---------------------------------------------------------------------------
# Initial sleep to allow tasks to appear in the queue
time.sleep(30)
while True:
    squeue_rows = gs.squeue.read()
    child_job_rows = [
        row for row in squeue_rows if row["ARRAY_JOB_ID"].startswith(str(array_job_id))
    ]
    for row in child_job_rows:
        state = row.get("STATE", "").upper()
        if state in ("FAILED", "CANCELLED"):
            raise RuntimeError(
                f"Array task {row[\'JOBID\']} entered state {state}"
            )
    still_running = len(child_job_rows) > 0
    if not still_running:
        break
    time.sleep(30)

print("All array tasks complete.", flush=True)

# Allow NFS to flush writes from the compute nodes before reading.
time.sleep(60)

# ---------------------------------------------------------------------------
# Import results into HDF5
# ---------------------------------------------------------------------------
data_dir = STAGING_DIR / "data" / RUN_CODE
run = Run(str(HDF5_PATH.parent), run_code=RUN_CODE)
fm = FortranMacro(run=run, out_file_code=OUT_FILE_CODE)
fm.import_results(
    data_dir,
    keep_on_failure=True,
    keep_tmpdir=KEEP_TMPDIR,
)
print("Results imported successfully.", flush=True)

# ---------------------------------------------------------------------------
# Clean up staging directory
# ---------------------------------------------------------------------------
if not KEEP_TMPDIR:
    shutil.rmtree(STAGING_DIR, ignore_errors=True)
    print("Staging directory removed.", flush=True)

print("Master job complete.", flush=True)
''')


def _generate_macro_master_py(
    staging_dir: Path,
    hdf5_path: Path,
    run_code: str,
    out_file_code: str,
    keep_tmpdir: bool,
) -> str:
    """Return the content of the macro master Python script with values substituted.

    :param staging_dir: Absolute path to the shared staging directory.
    :param hdf5_path: Absolute path to the target ``.h5`` file.
    :param run_code: Run identifier (HDF5 file stem).
    :param out_file_code: Output file code suffix for the Fortran binary.
    :param keep_tmpdir: Whether to preserve the staging directory on success.
    :return: Python script text.
    :rtype: str
    """
    return _MACRO_MASTER_PY_TEMPLATE.substitute(
        staging_dir=repr(str(staging_dir)),
        hdf5_path=repr(str(hdf5_path)),
        run_code=repr(run_code),
        out_file_code=repr(out_file_code),
        keep_tmpdir=repr(keep_tmpdir),
    )


def generate_macro_array_script(
    staging_dir: "Path | str",
    run_code: str,
    hdf5_path: "Path | str",
    executable: "Path | str",
    n_sims: int,
    in_code: str = "",
    out_code: str = "",
    *,
    partition: Optional[str] = None,
    fast_tmp_root: Optional[str] = None,
    slurm_log_dir: Optional["Path | str"] = None,
) -> str:
    """Generate a Slurm array job script for the macroscale Fortran simulation.

    Each array task handles one simulation (indexed by
    ``SLURM_ARRAY_TASK_ID``), calling
    :meth:`~lysis.execution.fortran_macro.FortranMacro.exec_in_workdir` with
    the task index.

    Each task passes ``--runCode {run_code}/{sim:02}`` to the Fortran binary so
    it writes directly to ``data/{run_code}/{sim:02}/`` relative to its CWD,
    isolating each task in its own subdirectory.  Setup files pre-staged in
    ``{staging_dir}/data/{run_code}/`` are symlinked into each per-sim subdir
    before execution and removed afterwards.

    **Single-tier** (default, ``fast_tmp_root=None``): Fortran writes its
    output directly to ``{staging_dir}/data/{run_code}/{sim:02}/``.

    **Two-tier** (``fast_tmp_root`` provided): Fortran writes to a private
    ``mktemp -d -p {fast_tmp_root}`` directory on the compute node.  Setup
    files are copied from staging first, the simulation runs locally, and the
    per-simulation output directory is then moved to staging.

    :param staging_dir: Path to the shared staging directory (must already
        exist when the script runs).
    :type staging_dir: Path or str
    :param run_code: Run identifier (used to name the output subdirectory).
    :type run_code: str
    :param hdf5_path: Full path to the run's ``.h5`` file.
    :type hdf5_path: Path or str
    :param executable: Path to the compiled Fortran macroscale binary.
    :type executable: Path or str
    :param n_sims: Total number of simulations (sets ``--array=0-{n_sims-1}``).
    :type n_sims: int
    :param in_code: Input file code suffix, defaults to ``""``.
    :type in_code: str, optional
    :param out_code: Output file code suffix, defaults to ``""``.
    :type out_code: str, optional
    :param partition: Slurm partition for ``#SBATCH --partition``, defaults
        to ``None`` (no partition directive).
    :type partition: str, optional
    :param fast_tmp_root: Root directory for fast node-local scratch.  When
        ``None`` (default) single-tier mode is used.
    :type fast_tmp_root: str, optional
    :param slurm_log_dir: Directory for Slurm ``.out`` logs.  Defaults to
        ``hdf5_path.parent / ".slurm"``.
    :type slurm_log_dir: Path or str, optional
    :return: Slurm array job bash script text.
    :rtype: str
    """
    staging_dir = Path(staging_dir)
    hdf5_path = Path(hdf5_path)
    executable = Path(executable)
    binary_name = executable.name
    if slurm_log_dir is None:
        slurm_log_dir = hdf5_path.parent / ".slurm"
    slurm_log_dir = Path(slurm_log_dir)

    sbatch_opts = {
        "array": f"0-{n_sims - 1}",
        "out": str(slurm_log_dir / "array.slurm-%A_%a.out"),
        "nodes": 1,
        "mem": 3096,
        "ntasks": 1,
        "cpus-per-task": 1,
        "exclusive=user": "",
    }
    if partition:
        sbatch_opts["partition"] = partition

    if fast_tmp_root is None:
        # ------------------------------------------------------------------
        # Single-tier: Fortran writes directly to the shared staging dir.
        # The executable is pre-staged by submit_macro_slurm_job() before
        # any array tasks start, so no cp is needed here.
        # ------------------------------------------------------------------
        execute = f"""\
# Execute macroscale Fortran binary via lysis
source ~/.bashrc && source ~/lysis.sh
python -c "
import os
from pathlib import Path
from lysis.execution.fortran_macro import FortranMacro
fm = FortranMacro.from_hdf5(
    '{hdf5_path}', '{staging_dir}/{binary_name}',
    in_file_code='{in_code}', out_file_code='{out_code}',
    index=int(os.environ['SLURM_ARRAY_TASK_ID'])
)
fm.exec_in_workdir(Path('{staging_dir}'))
" """

        sections = [execute]

    else:
        # ------------------------------------------------------------------
        # Two-tier: copy setup files → local fast dir, run, then mv to staging
        # ------------------------------------------------------------------
        setup = f"""\
# Setup — create local fast dir, copy setup files and binary
SIM=$(printf "%02d" ${{SLURM_ARRAY_TASK_ID}})
local_work_dir=$(mktemp -d -p "{fast_tmp_root}")
local_datadir="${{local_work_dir}}/data/{run_code}"
staging_datadir="{staging_dir}/data/{run_code}"
mkdir -p "${{local_datadir}}"
cp "${{staging_datadir}}/"* "${{local_datadir}}/"
cp "{executable}" "${{local_work_dir}}/" """

        execute = f"""\
# Execute macroscale Fortran binary into local fast storage
source ~/.bashrc && source ~/lysis.sh
python -c "
import os
from pathlib import Path
from lysis.execution.fortran_macro import FortranMacro
fm = FortranMacro.from_hdf5(
    '{hdf5_path}', '${{local_work_dir}}/{binary_name}',
    in_file_code='{in_code}', out_file_code='{out_code}',
    index=int(os.environ['SLURM_ARRAY_TASK_ID'])
)
fm.exec_in_workdir(Path('${{local_work_dir}}'))
" """

        move = """\
# Move per-simulation output to shared staging, then clean up local dir
mv "${local_datadir}/${SIM}" "${staging_datadir}/"
rm -rf "${local_work_dir}" """

        sections = [setup, execute, move]

    return gs.scripts.plain(sections, **sbatch_opts)


def submit_macro_slurm_job(
    hdf5_path: "Path | str",
    executable: "Path | str",
    *,
    staging_root: Optional["Path | str"] = None,
    partition: Optional[str] = None,
    fast_tmp_root: Optional[str] = None,
    keep_tmpdir: bool = False,
    in_code: str = "",
    out_code: str = "",
) -> int:
    """Submit the full HDF5-integrated macroscale workflow as a Slurm master job.

    Creates a unique staging directory, pre-stages setup files (macroscale_in
    text files and ``params.json``), writes an array job script and a master
    orchestration script, and submits the master.  The master job orchestrates
    the array tasks, imports results into the HDF5 file, and optionally cleans
    up.

    Directory naming
    ~~~~~~~~~~~~~~~~
    The staging directory is created with::

        tempfile.mkdtemp(
            prefix=f"lysis-macro-{run_code}-",
            dir=staging_root or hdf5_path.parent,
        )

    The ``dir`` argument is **always explicit** — ``/tmp`` is never used.

    :param hdf5_path: Full path to the run's ``.h5`` file.  Must contain both
        ``micro_params`` and ``macro_params``.
    :type hdf5_path: Path or str
    :param executable: Path to the compiled Fortran macroscale binary.
    :type executable: Path or str
    :param staging_root: Root directory under which the staging temp dir is
        created.  Defaults to ``hdf5_path.parent``.
    :type staging_root: Path or str, optional
    :param partition: Slurm partition for both master and array jobs.
    :type partition: str, optional
    :param fast_tmp_root: Root directory for fast node-local scratch (opt-in
        two-tier storage).  When ``None`` (default) array tasks write directly
        to the shared staging directory.
    :type fast_tmp_root: str, optional
    :param keep_tmpdir: If ``True``, preserve the staging directory after the
        master job completes (useful for debugging).
    :type keep_tmpdir: bool, optional
    :param in_code: Input file code suffix for the Fortran binary,
        defaults to ``""``.
    :type in_code: str, optional
    :param out_code: Output file code suffix for the Fortran binary,
        defaults to ``""``.
    :type out_code: str, optional
    :return: Master Slurm job ID.
    :rtype: int
    :raises subprocess.CalledProcessError: If ``sbatch`` fails.
    :raises ValueError: If the HDF5 file does not contain ``macro_params``.
    """
    from lysis.execution.fortran_macro import FortranMacro  # noqa: PLC0415

    hdf5_path = Path(hdf5_path).resolve()
    executable = Path(executable).resolve()
    run_code = hdf5_path.stem

    # Create unique staging directory — never use /tmp
    staging_root_dir = (
        Path(staging_root) if staging_root is not None else hdf5_path.parent
    )
    staging_dir = Path(tempfile.mkdtemp(
        prefix=f"lysis-macro-{run_code}-",
        dir=str(staging_root_dir),
    ))

    # Slurm log directory (sibling to the HDF5 file)
    slurm_log_dir = hdf5_path.parent / ".slurm"
    slurm_log_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Pre-stage setup files (macroscale_in + params.json)
    # ------------------------------------------------------------------
    staging_data_dir = staging_dir / "data" / run_code
    staging_data_dir.mkdir(parents=True, exist_ok=True)
    fm_setup = FortranMacro.from_hdf5(
        hdf5_path, str(executable), in_file_code=in_code, out_file_code=out_code,
    )
    fm_setup._write_setup_files(staging_data_dir)
    # Pre-stage the executable so no array task needs to cp it — avoids a
    # race where multiple tasks simultaneously try to create the same file.
    shutil.copy2(executable, staging_dir)

    # Read n_sims from the run's macro_params
    n_sims = fm_setup.run.macro_params.macro_simulations

    # ------------------------------------------------------------------
    # Write array job script
    # ------------------------------------------------------------------
    array_script = generate_macro_array_script(
        staging_dir, run_code, hdf5_path, executable, n_sims,
        in_code=in_code, out_code=out_code,
        partition=partition, fast_tmp_root=fast_tmp_root,
        slurm_log_dir=slurm_log_dir,
    )
    array_path = staging_dir / "array.sh"
    array_path.write_text(array_script)
    array_path.chmod(0o755)

    # ------------------------------------------------------------------
    # Write master.py
    # ------------------------------------------------------------------
    master_py_content = _generate_macro_master_py(
        staging_dir, hdf5_path, run_code, out_code, keep_tmpdir,
    )
    master_py_path = staging_dir / "master.py"
    master_py_path.write_text(master_py_content)

    # ------------------------------------------------------------------
    # Write master bash wrapper and submit
    # ------------------------------------------------------------------
    sbatch_opts = {
        "job-name": f"lysis-macro-{run_code}",
        "out": str(slurm_log_dir / "master.slurm-%j.out"),
        "nodes": 1,
        "mem": 3096,
        "ntasks": 1,
        "cpus-per-task": 1,
        "exclusive=user": "",
    }
    if partition:
        sbatch_opts["partition"] = partition

    master_sh_content = gs.scripts.plain(
        [
            "source ~/.bashrc && source ~/lysis.sh",
            f'python "{master_py_path}"',
        ],
        **sbatch_opts,
    )
    master_sh_path = staging_dir / "master.sh"
    master_sh_path.write_text(master_sh_content)
    master_sh_path.chmod(0o755)

    return gs.sbatch([str(master_sh_path)])
