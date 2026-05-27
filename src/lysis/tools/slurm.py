"""Slurm dispatch helpers for the Lysis microscale and macroscale workflows.

This module provides utilities for submitting and monitoring the Fortran
simulations via the Slurm workload manager, using
`GooseSLURM <https://gooseslurm.readthedocs.io/>`_ for script generation
and job submission.

Workflow overview
-----------------
Two scales (``micro`` / ``macro``) and two execution patterns:

* **Single-child micro** (legacy): :func:`submit_micro_slurm_job` with
  ``num_children=None`` (its default) launches one Fortran child that
  runs all microscale simulations sequentially, then imports.
* **Array (both scales)**: :func:`submit_macro_slurm_job` and
  :func:`submit_micro_slurm_job` (``num_children >= 1``) submit a Slurm
  array of tasks via :func:`generate_array_script`; the master job polls
  ``squeue``, optionally concatenates per-task outputs (microscale only),
  and imports.

The array-mode internals — :class:`_SlurmJobSpec`,
:data:`_ARRAY_MASTER_PY_TEMPLATE`, :func:`generate_array_script`, and
:func:`submit_slurm_job` — are shared between the two scales.

Two-tier storage (opt-in)
~~~~~~~~~~~~~~~~~~~~~~~~~
When *fast_tmp_root* is provided each task writes Fortran output to fast
node-local scratch (``mktemp -d -p <fast_tmp_root>``) rather than the
shared staging directory, then moves the results to staging before exiting.
This optimises I/O on clusters with NVMe scratch on compute nodes.

Source pinning
~~~~~~~~~~~~~~
At submit time both ``submit_slurm_job`` and ``submit_micro_slurm_job``
snapshot ``src/lysis/`` into ``{staging_dir}/python_src/`` (see
:func:`_snapshot_lysis_src`) and bake ``PYTHONPATH=<snapshot>`` into the
generated master and per-task scripts.  ``PYTHONPATH`` precedes
site-packages on ``sys.path``, so the snapshot shadows the live editable
install — every Python invocation in the job (master and array tasks)
imports the frozen package, immune to source edits made while the job
is queued or running.  Only ``src/lysis/`` is pinned; third-party
dependencies still come from ``.venv`` and are pinned separately by
``uv sync --frozen`` (see :func:`_ensure_env_synced`).

.. important::

   All temporary directories are created with an explicit *dir* argument
   (never defaulting to ``/tmp``) because most HPC nodes use a ramdisk
   for ``/tmp`` that is too small for simulation data.
"""

import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from string import Template
from typing import Dict, List, Mapping, Optional

import GooseSLURM as gs


__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


# ---------------------------------------------------------------------------
# Environment setup for generated Slurm scripts
# ---------------------------------------------------------------------------
#
# Generated scripts no longer depend on per-user dotfiles (``~/.bashrc`` /
# ``~/lysis.sh``).  Instead, every script emits a self-contained preamble
# that loads the Intel compiler runtime via LMod and invokes Python through
# the project's ``uv``-managed environment.  The submitting user's repo
# root is baked in at script-generation time from ``lysis.__file__``.

#: Default LMod module providing the Fortran binary's runtime libraries.
#: Assumed to be available on the cluster's LMod stack for all users.
#: Callers can override per-job via the ``compiler_module`` keyword argument
#: on the public submit/generate functions (exposed as ``--compiler`` on the
#: ``lysis run-micro`` and ``lysis run-macro`` CLI commands).
DEFAULT_COMPILER_MODULE: str = "intel-compilers/2023"


def _repo_root() -> Path:
    """Return the absolute path to the lysis repository root.

    Computed from the installed ``lysis`` package location
    (``src/lysis/__init__.py`` → repo root is two parents up).  This is the
    submitting user's *own* checkout — the path baked into generated scripts
    is correct for that user without depending on any shared convention.
    """
    import lysis  # noqa: PLC0415 — avoids pulling lysis at module import

    return Path(lysis.__file__).resolve().parents[2]


def _env_preamble(
    repo_root: Path,
    compiler_module: str,
    *,
    pythonpath: "Path | str | None" = None,
) -> str:
    """Return the bash preamble that prepares a generated script's environment.

    Loads the requested Fortran compiler runtime via LMod (required by the
    Fortran binary) and puts ``~/.local/bin`` on PATH so ``uv`` is
    discoverable on compute nodes.  Defensive against non-login shells by
    sourcing ``lmod.sh`` when the ``module`` function isn't already defined.

    When *pythonpath* is supplied, an ``export PYTHONPATH=...`` line is
    appended.  ``uv run`` inherits environment variables from the calling
    shell, and PYTHONPATH entries come before site-packages on ``sys.path``,
    so a snapshot directory containing ``lysis/`` will shadow the editable
    install — pinning the Python source for the duration of the job.

    :param repo_root: Absolute path to the lysis repo root (currently
        unused inside the preamble itself; included for future extensions
        like ``cd`` into the repo).
    :type repo_root: Path
    :param compiler_module: LMod module spec to load (e.g.
        ``"intel-compilers/2023"`` or ``"intel-compilers/2024"``).  See
        :data:`DEFAULT_COMPILER_MODULE` for the project default.
    :type compiler_module: str
    :param pythonpath: When set, exported as ``PYTHONPATH`` so generated
        scripts import ``lysis`` from that path rather than the live
        editable install (see :func:`_snapshot_lysis_src`).  ``None``
        (default) leaves ``PYTHONPATH`` untouched.
    :type pythonpath: Path or str, optional
    :return: Multi-line bash snippet, no trailing newline.
    :rtype: str
    """
    lines = [
        "# Lysis environment setup (user-independent)",
        "[ -z \"${LMOD_CMD:-}\" ] && [ -f /etc/profile.d/lmod.sh ] "
        "&& source /etc/profile.d/lmod.sh",
        f"module purge && module load {compiler_module}",
        'export PATH="$HOME/.local/bin:$PATH"',
    ]
    if pythonpath is not None:
        lines.append(f'export PYTHONPATH="{pythonpath}"')
    return "\n".join(lines)


def _snapshot_lysis_src(repo_root: Path, staging_dir: Path) -> Path:
    """Copy ``src/lysis/`` into the staging dir and return the importable root.

    Produces ``{staging_dir}/python_src/lysis/`` from
    ``{repo_root}/src/lysis/``, excluding ``__pycache__`` and ``.pyc``
    files so the snapshot starts with no stale bytecode.  The returned
    path is the directory to be placed on ``PYTHONPATH``
    (``{staging_dir}/python_src``), so generated Slurm scripts import
    ``lysis`` from this frozen copy instead of the live editable install.

    Pins only the Python package — third-party dependencies still come
    from the project's ``.venv`` (and from ``uv.lock``, frozen at submit
    time by :func:`_ensure_env_synced`).

    :param repo_root: Absolute path to the lysis repo root.
    :type repo_root: Path
    :param staging_dir: Existing per-job staging directory.
    :type staging_dir: Path
    :return: The directory to assign to ``PYTHONPATH`` (the parent of the
        snapshot ``lysis/`` package).
    :rtype: Path
    """
    src = repo_root / "src" / "lysis"
    pythonpath = staging_dir / "python_src"
    dst = pythonpath / "lysis"
    pythonpath.mkdir(parents=True, exist_ok=True)
    shutil.copytree(
        src, dst, ignore=shutil.ignore_patterns("__pycache__", "*.pyc")
    )
    return pythonpath


def _uv_python_prefix(repo_root: Path) -> str:
    """Return the command prefix used in place of bare ``python``.

    Invokes Python through ``uv run`` pinned to the repo's ``uv.lock``
    (``--frozen`` refuses any resolver activity), so generated scripts
    use exactly the dependency set committed alongside the source.

    :param repo_root: Absolute path to the lysis repo root.
    :type repo_root: Path
    :return: Shell command prefix, e.g.
        ``uv run --project /home/.../lysis --frozen python``.
    :rtype: str
    """
    return f'uv run --project "{repo_root}" --frozen python'


def _ensure_env_synced(repo_root: Path) -> None:
    """Run ``uv sync --frozen`` once on the submit host before sbatch.

    Guarantees the project's ``.venv`` matches ``uv.lock`` before any
    Slurm tasks fire — array tasks then start in a known-good state
    without racing each other to create or populate the venv.

    :param repo_root: Absolute path to the lysis repo root.
    :type repo_root: Path
    :raises subprocess.CalledProcessError: If ``uv sync --frozen`` fails
        (e.g. lockfile drift, uv not on PATH).
    """
    subprocess.run(
        ["uv", "sync", "--frozen"],
        cwd=str(repo_root),
        check=True,
    )


# ---------------------------------------------------------------------------
# User-facing #SBATCH override helpers
# ---------------------------------------------------------------------------


def parse_sbatch_tokens(tokens) -> Dict[str, Optional[str]]:
    """Parse a sequence of ``--sbatch`` CLI tokens into an overrides mapping.

    Token syntax (mirrors what users type after ``--sbatch``):

    * ``KEY=VALUE`` → ``{KEY: VALUE}``.  Splits on the **first** ``=``;
      anything after the first ``=`` is the value.  Suitable for the common
      ``--mem=4MB``, ``--time=01:00:00`` style.
    * ``KEY`` (no ``=``) → ``{KEY: ""}``.  Renders as a flag-only ``#SBATCH
      --KEY`` (e.g. ``--sbatch hold``).
    * ``^KEY`` (leading caret) → ``{KEY: None}``.  Marks the default whose
      dict key matches ``KEY`` exactly for removal.  Use this to drop a
      default whose key itself contains ``=`` — for example the built-in
      ``exclusive=user`` is dropped via ``--sbatch '^exclusive=user'``.

    The returned mapping is in the shape consumed by
    :func:`apply_sbatch_overrides` (``None`` = remove, anything else = set).

    :param tokens: Iterable of raw ``--sbatch`` token strings.
    :type tokens: Iterable[str]
    :return: Overrides mapping suitable for :func:`apply_sbatch_overrides`.
    :rtype: dict[str, str or None]
    :raises ValueError: If any token is empty or has an empty key after
        the caret/``=`` split.
    """
    result: Dict[str, Optional[str]] = {}
    for tok in tokens:
        if not tok:
            raise ValueError("empty --sbatch token")
        if tok.startswith("^"):
            key = tok[1:]
            if not key:
                raise ValueError(f"empty key in --sbatch override: {tok!r}")
            result[key] = None
        elif "=" in tok:
            key, _, value = tok.partition("=")
            if not key:
                raise ValueError(f"empty key in --sbatch override: {tok!r}")
            result[key] = value
        else:
            result[tok] = ""
    return result


def apply_sbatch_overrides(
    base: Mapping[str, object],
    overrides: Optional[Mapping[str, Optional[str]]],
) -> Dict[str, object]:
    """Return a copy of *base* with *overrides* merged in.

    A value of ``None`` removes that key from the result (silent no-op if
    the key isn't present in *base*).  Any other value sets the key.

    :param base: The default ``sbatch_opts`` dict to start from.
    :type base: Mapping[str, object]
    :param overrides: Mapping of user overrides; ``None`` value = remove.
    :type overrides: Mapping[str, str or None] or None
    :return: New dict with overrides applied.  *base* is not mutated.
    :rtype: dict[str, object]
    """
    result: Dict[str, object] = dict(base)
    if not overrides:
        return result
    for key, value in overrides.items():
        if value is None:
            result.pop(key, None)
        else:
            result[key] = value
    return result


# ---------------------------------------------------------------------------
# Single-child microscale internals (legacy path)
# ---------------------------------------------------------------------------

#: Template for the single-child micro master Python script.  Substitution
#: keys use ``$name`` syntax so Python braces in the template body are
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
BINARY_NAME = $binary_name
KEEP_TMPDIR = $keep_tmpdir
HISTORICAL_BINARY_ATTRS = $historical_binary_attrs

# ---------------------------------------------------------------------------
# Submit child jobs
# ---------------------------------------------------------------------------
child_ids = []
for script in sorted(STAGING_DIR.glob(f"lysis-micro-child__{RUN_CODE}*.sh")):
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
fm = FortranMicro(
    run=run,
    out_file_code=FILE_CODE,
    executable=str(STAGING_DIR / BINARY_NAME),
    skip_binary_verification=(HISTORICAL_BINARY_ATTRS is not None),
    historical_binary_attrs=HISTORICAL_BINARY_ATTRS,
)
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
    binary_name: str,
    keep_tmpdir: bool,
    historical_binary_attrs: Optional[dict] = None,
) -> str:
    """Return the content of the single-child micro master Python script.

    :param historical_binary_attrs: Pre-computed binary-provenance dict
        from
        :func:`~lysis.tools.provenance.gather_historical_binary_provenance`
        for the historical-build workflow.  Baked into the master script
        as a literal so the master can stamp it in place of the default
        ``<binary> --version`` query.  ``None`` (default) preserves the
        normal-mode behaviour.
    :type historical_binary_attrs: dict or None
    """
    return _MICRO_MASTER_PY_TEMPLATE.substitute(
        staging_dir=repr(str(staging_dir)),
        hdf5_path=repr(str(hdf5_path)),
        run_code=repr(run_code),
        file_code=repr(file_code),
        binary_name=repr(binary_name),
        keep_tmpdir=repr(keep_tmpdir),
        historical_binary_attrs=repr(historical_binary_attrs),
    )


# ---------------------------------------------------------------------------
# Unified array-mode internals (used by both micro and macro array paths)
# ---------------------------------------------------------------------------


@dataclass
class _SlurmJobSpec:
    """Description of one array-mode Slurm dispatch.

    Holds the scale-specific knobs that distinguish a microscale array
    job from a macroscale array job.  Consumed by
    :func:`generate_array_script` (bash) and :func:`submit_slurm_job`
    (master Python script + sbatch).

    :ivar scale: ``"micro"`` or ``"macro"`` — used in script and job names.
    :vartype scale: str
    :ivar runner_module: Dotted import path of the Fortran runner class
        (e.g. ``"lysis.execution.fortran_micro"``).
    :vartype runner_module: str
    :ivar runner_class: Class name (e.g. ``"FortranMicro"``).
    :vartype runner_class: str
    :ivar num_array_tasks: Total number of array tasks (sets
        ``--array=0-{num_array_tasks - 1}``).
    :vartype num_array_tasks: int
    :ivar needs_concat_step: Whether the master script must call
        ``fm.concatenate_child_outputs(data_dir, num_array_tasks)`` before
        ``import_results`` (microscale only).
    :vartype needs_concat_step: bool
    :ivar nfs_wait_seconds: Seconds the master sleeps after the last task
        leaves the queue, to let NFS flush writes from compute nodes.
    :vartype nfs_wait_seconds: int
    :ivar prestage_executable: Whether :func:`submit_slurm_job` should copy
        the executable into the staging dir before submission (avoids a
        race where multiple array tasks try to ``cp`` the same file).
    :vartype prestage_executable: bool
    :ivar out_file_code: Output-file suffix passed to the runner.
    :vartype out_file_code: str
    :ivar in_file_code: Input-file suffix (macro only; ``None`` for micro).
    :vartype in_file_code: str or None
    :ivar num_children: Total number of children in the partition (micro
        array only).  Threaded into ``from_hdf5`` so the runner sees the
        partition contract; ``None`` for macro.
    :vartype num_children: int or None
    """

    scale: str
    runner_module: str
    runner_class: str
    num_array_tasks: int
    needs_concat_step: bool
    nfs_wait_seconds: int
    prestage_executable: bool
    out_file_code: str
    in_file_code: Optional[str] = None
    num_children: Optional[int] = None


#: Template for the master Python script in array mode (both scales).
#: Substitution keys use ``$name`` syntax so Python braces in the
#: template body are untouched.
_ARRAY_MASTER_PY_TEMPLATE = Template('''\
#!/usr/bin/env python3
"""Master Slurm job: submit array tasks, poll for completion, $concat_blurb import."""
import time
import shutil
from pathlib import Path

import GooseSLURM as gs
from lysis.config.run import Run
from $runner_module import $runner_class

STAGING_DIR = Path($staging_dir)
HDF5_PATH = Path($hdf5_path)
RUN_CODE = $run_code
OUT_FILE_CODE = $out_file_code
BINARY_NAME = $binary_name
KEEP_TMPDIR = $keep_tmpdir
NUM_ARRAY_TASKS = $num_array_tasks
NFS_WAIT_SECONDS = $nfs_wait_seconds
HISTORICAL_BINARY_ATTRS = $historical_binary_attrs

# ---------------------------------------------------------------------------
# Submit array job
# ---------------------------------------------------------------------------
array_script = STAGING_DIR / f"lysis-$scale-array__{RUN_CODE}.sh"
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
        row for row in squeue_rows if row["ARRAY_JOB_ID"] == str(array_job_id)
    ]
    for row in child_job_rows:
        state = row.get("STATE", "").upper()
        if state in ("FAILED", "CANCELLED"):
            raise RuntimeError(
                f"Array task {row['JOBID']} entered state {state}"
            )
    still_running = len(child_job_rows) > 0
    if not still_running:
        break
    time.sleep(30)

print("All array tasks complete.", flush=True)

# Allow NFS to flush writes from the compute nodes before reading.
time.sleep(NFS_WAIT_SECONDS)

# ---------------------------------------------------------------------------
# Import results into HDF5
# ---------------------------------------------------------------------------
data_dir = STAGING_DIR / "data" / RUN_CODE
run = Run(str(HDF5_PATH.parent), run_code=RUN_CODE)
fm = $runner_class(
    run=run,
    out_file_code=OUT_FILE_CODE,
    executable=str(STAGING_DIR / BINARY_NAME),
    skip_binary_verification=(HISTORICAL_BINARY_ATTRS is not None),
    historical_binary_attrs=HISTORICAL_BINARY_ATTRS,
)
$concat_block
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


def _generate_array_master_py(
    spec: _SlurmJobSpec,
    staging_dir: Path,
    hdf5_path: Path,
    run_code: str,
    binary_name: str,
    keep_tmpdir: bool,
    historical_binary_attrs: Optional[dict] = None,
) -> str:
    """Return the master Python script for an array-mode dispatch.

    :param historical_binary_attrs: Pre-computed binary-provenance dict
        from
        :func:`~lysis.tools.provenance.gather_historical_binary_provenance`
        for the historical-build workflow.  Baked into the master script
        as a literal so the master can stamp it in place of the default
        ``<binary> --version`` query.  ``None`` (default) preserves the
        normal-mode behaviour.
    :type historical_binary_attrs: dict or None
    """
    if spec.needs_concat_step:
        concat_block = (
            f"fm.concatenate_child_outputs(data_dir, num_children={spec.num_array_tasks})"
        )
        concat_blurb = "concatenate per-task outputs, then"
    else:
        concat_block = "# (no concatenation step — runner imports per-task subdirs directly)"
        concat_blurb = "then"

    return _ARRAY_MASTER_PY_TEMPLATE.substitute(
        scale=spec.scale,
        runner_module=spec.runner_module,
        runner_class=spec.runner_class,
        staging_dir=repr(str(staging_dir)),
        hdf5_path=repr(str(hdf5_path)),
        run_code=repr(run_code),
        out_file_code=repr(spec.out_file_code),
        binary_name=repr(binary_name),
        keep_tmpdir=repr(keep_tmpdir),
        num_array_tasks=spec.num_array_tasks,
        nfs_wait_seconds=spec.nfs_wait_seconds,
        concat_block=concat_block,
        concat_blurb=concat_blurb,
        historical_binary_attrs=repr(historical_binary_attrs),
    )


def _runner_init_line(
    spec: _SlurmJobSpec,
    hdf5_path: Path,
    binary_path: str,
    *,
    skip_binary_verification: bool = False,
) -> str:
    """Build the ``from_hdf5(...)`` call used inside the array task ``python -c``.

    Single-line so it composes cleanly inside the bash heredoc.  When
    *skip_binary_verification* is True the call emits
    ``skip_binary_verification=True`` so the historical-build workflow's
    binary↔source mismatch does not crash ``_verify_binary_version`` at
    task start.
    """
    args = [f"'{hdf5_path}'", f"'{binary_path}'"]
    if spec.in_file_code is not None:
        args.append(f"in_file_code='{spec.in_file_code}'")
    args.append(f"out_file_code='{spec.out_file_code}'")
    args.append("index=int(os.environ['SLURM_ARRAY_TASK_ID'])")
    if spec.num_children is not None:
        args.append(f"num_children={spec.num_children}")
    if skip_binary_verification:
        args.append("skip_binary_verification=True")
    return f"{spec.runner_class}.from_hdf5({', '.join(args)})"


def generate_array_script(
    spec: _SlurmJobSpec,
    staging_dir: "Path | str",
    run_code: str,
    hdf5_path: "Path | str",
    executable: "Path | str",
    *,
    partition: Optional[str] = None,
    fast_tmp_root: Optional[str] = None,
    slurm_log_dir: Optional["Path | str"] = None,
    compiler_module: str = DEFAULT_COMPILER_MODULE,
    sbatch_overrides: Optional[Mapping[str, Optional[str]]] = None,
    historical_binary_attrs: Optional[dict] = None,
    pythonpath: "Path | str | None" = None,
) -> str:
    """Generate a Slurm array job bash script for an array-mode dispatch.

    Each array task indexed by ``SLURM_ARRAY_TASK_ID`` calls
    ``{spec.runner_class}.from_hdf5(...).exec_in_workdir(...)``.

    For the *macroscale* path each task isolates its output in a
    per-simulation subdirectory ``data/{run_code}/{sim:02}/`` (managed by
    :class:`~lysis.execution.fortran_macro.FortranMacro` itself).  For the
    *microscale* array path each task writes flat files with a ``__NN``
    suffix into ``data/{run_code}/``; the master concatenates them after
    all tasks complete.

    **Single-tier** (default, ``fast_tmp_root=None``): Fortran writes its
    output directly into ``{staging_dir}/data/{run_code}/``.  The
    executable is *not* copied here — :func:`submit_slurm_job` pre-stages
    it before any task starts (avoids a cp race).

    **Two-tier** (``fast_tmp_root`` provided): Fortran writes to a private
    ``mktemp -d -p {fast_tmp_root}`` directory on the compute node; the
    output is moved to the shared staging dir before exit.

    :param spec: Scale-specific dispatch description.
    :type spec: _SlurmJobSpec
    :param staging_dir: Path to the shared staging directory.
    :type staging_dir: Path or str
    :param run_code: Run identifier.
    :type run_code: str
    :param hdf5_path: Full path to the run's ``.h5`` file.
    :type hdf5_path: Path or str
    :param executable: Path to the compiled Fortran binary.
    :type executable: Path or str
    :param partition: Slurm partition for ``#SBATCH --partition``.
    :type partition: str, optional
    :param fast_tmp_root: Root directory for fast node-local scratch.
    :type fast_tmp_root: str, optional
    :param slurm_log_dir: Directory for Slurm ``.out`` logs.  Defaults to
        ``hdf5_path.parent / ".slurm"``.
    :type slurm_log_dir: Path or str, optional
    :param compiler_module: LMod module spec for the Fortran compiler
        runtime, baked into each generated task script.  Defaults to
        :data:`DEFAULT_COMPILER_MODULE`.
    :type compiler_module: str, optional
    :param sbatch_overrides: User overrides for the per-task ``#SBATCH``
        header, in the shape returned by :func:`parse_sbatch_tokens`.
        A ``None`` value removes the matching default; any other value
        sets/overrides the key.  Defaults to ``None`` (no overrides).
    :type sbatch_overrides: Mapping[str, str or None], optional
    :param historical_binary_attrs: When non-``None``, indicates the
        binary was built from a historical commit (see
        :func:`~lysis.execution.historical_build.build_historical_binary`).
        The generated ``from_hdf5(...)`` call receives
        ``skip_binary_verification=True`` so ``_verify_binary_version``
        does not crash on the intentional binary↔source mismatch.  The
        dict itself is not baked into the task — the master script
        owns the actual provenance stamping.
    :type historical_binary_attrs: dict or None
    :param pythonpath: When set, baked into the generated script as
        ``export PYTHONPATH=...`` so each task imports ``lysis`` from
        that path rather than the live editable install.  ``None``
        (default) leaves ``PYTHONPATH`` untouched.
    :type pythonpath: Path or str, optional
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
    skip_binary_verification = historical_binary_attrs is not None

    repo_root = _repo_root()
    env_preamble = _env_preamble(
        repo_root, compiler_module, pythonpath=pythonpath
    )
    py = _uv_python_prefix(repo_root)

    sbatch_opts = {
        "array": f"0-{spec.num_array_tasks - 1}",
        "job-name": f"lysis-{spec.scale}-array__{run_code}__%a",
        "out": str(slurm_log_dir / f"lysis-{spec.scale}-array__{run_code}__%a.out"),
        "nodes": 1,
        "mem": 3096,
        "ntasks": 1,
        "cpus-per-task": 1,
        "exclusive=user": "",
    }
    if partition:
        sbatch_opts["partition"] = partition
    sbatch_opts = apply_sbatch_overrides(sbatch_opts, sbatch_overrides)

    if fast_tmp_root is None:
        # ------------------------------------------------------------------
        # Single-tier: Fortran writes directly to the shared staging dir.
        # The executable is pre-staged by submit_slurm_job() before any
        # array tasks start, so no cp is needed here.
        # ------------------------------------------------------------------
        init = _runner_init_line(
            spec, hdf5_path, f"{staging_dir}/{binary_name}",
            skip_binary_verification=skip_binary_verification,
        )
        execute = f"""\
# Execute {spec.scale}scale Fortran binary via lysis
{env_preamble}
{py} -c "
import os
from pathlib import Path
from {spec.runner_module} import {spec.runner_class}
fm = {init}
fm.exec_in_workdir(Path('{staging_dir}'))
" """

        sections = [execute]

    else:
        # ------------------------------------------------------------------
        # Two-tier: copy setup files → local fast dir, run, then mv to staging
        # ------------------------------------------------------------------
        if spec.scale == "macro":
            # Macro tasks write per-sim subdirs; mv just the one this task owns.
            mv_section = 'mv "${local_datadir}/${SIM}" "${staging_datadir}/"'
        else:
            # Micro array tasks write flat __NN-suffixed files; mv all of them.
            mv_section = 'mv "${local_datadir}"/* "${staging_datadir}/"'

        # Source the binary from the pre-staged staging-dir copy, not the
        # caller's --executable path: under --fortran-commit the build dir
        # is already gone by the time this script runs.
        setup = f"""\
# Setup — create local fast dir, copy setup files and binary
SIM=$(printf "%02d" ${{SLURM_ARRAY_TASK_ID}})
local_work_dir=$(mktemp -d -p "{fast_tmp_root}")
local_datadir="${{local_work_dir}}/data/{run_code}"
staging_datadir="{staging_dir}/data/{run_code}"
mkdir -p "${{local_datadir}}"
cp "${{staging_datadir}}/"* "${{local_datadir}}/"
cp "{staging_dir}/{binary_name}" "${{local_work_dir}}/" """

        init = _runner_init_line(
            spec, hdf5_path, f"${{local_work_dir}}/{binary_name}",
            skip_binary_verification=skip_binary_verification,
        )
        execute = f"""\
# Execute {spec.scale}scale Fortran binary into local fast storage
{env_preamble}
{py} -c "
import os
from pathlib import Path
from {spec.runner_module} import {spec.runner_class}
fm = {init}
fm.exec_in_workdir(Path('${{local_work_dir}}'))
" """

        move = f"""\
# Move per-task output to shared staging, then clean up local dir
{mv_section}
rm -rf "${{local_work_dir}}" """

        sections = [setup, execute, move]

    return gs.scripts.plain(sections, **sbatch_opts)


def submit_slurm_job(
    spec: _SlurmJobSpec,
    hdf5_path: Path,
    executable: Path,
    staging_root_dir: Path,
    *,
    partition: Optional[str] = None,
    fast_tmp_root: Optional[str] = None,
    keep_tmpdir: bool = False,
    compiler_module: str = DEFAULT_COMPILER_MODULE,
    sbatch_overrides: Optional[Mapping[str, Optional[str]]] = None,
    historical_binary_attrs: Optional[dict] = None,
) -> int:
    """Stage scripts and submit a master Slurm job for an array-mode dispatch.

    Pre-stages setup files (and optionally the executable, per ``spec``),
    writes the array bash script and master Python+bash scripts into a
    unique staging directory, and ``sbatch``-submits the master.

    :param spec: Scale-specific dispatch description.
    :type spec: _SlurmJobSpec
    :param hdf5_path: Resolved absolute path to the ``.h5`` file.
    :type hdf5_path: Path
    :param executable: Resolved absolute path to the Fortran binary.
    :type executable: Path
    :param staging_root_dir: Root directory under which the staging temp
        dir is created (must already exist; the temp dir itself is
        created here).
    :type staging_root_dir: Path
    :param partition: Slurm partition for both master and array jobs.
    :type partition: str, optional
    :param fast_tmp_root: Root directory for fast node-local scratch.
    :type fast_tmp_root: str, optional
    :param keep_tmpdir: Preserve staging dir after the master job.
    :type keep_tmpdir: bool, optional
    :param compiler_module: LMod module spec for the Fortran compiler
        runtime, baked into both the master and array task scripts.
        Defaults to :data:`DEFAULT_COMPILER_MODULE`.
    :type compiler_module: str, optional
    :param sbatch_overrides: User overrides for the ``#SBATCH`` header,
        applied to BOTH the master job and each array task (shape returned
        by :func:`parse_sbatch_tokens`).  Defaults to ``None``.
    :type sbatch_overrides: Mapping[str, str or None], optional
    :param historical_binary_attrs: Pre-computed binary-provenance dict
        from
        :func:`~lysis.tools.provenance.gather_historical_binary_provenance`
        for the historical-build workflow.  Baked into the master Python
        script so the master stamps it (with binary verification skipped)
        in place of the default ``<binary> --version`` query.  ``None``
        (default) preserves normal-mode provenance.
    :type historical_binary_attrs: dict, optional
    :return: Master Slurm job ID.
    :rtype: int
    """
    run_code = hdf5_path.stem

    # Sync the project venv against uv.lock before any Slurm tasks fire,
    # so array tasks never race each other on first-use venv creation.
    repo_root = _repo_root()
    _ensure_env_synced(repo_root)

    # Create unique staging directory — never use /tmp
    staging_dir = Path(tempfile.mkdtemp(
        prefix=f"lysis-{spec.scale}-{run_code}-",
        dir=str(staging_root_dir),
    ))

    # Slurm log directory (sibling to the HDF5 file)
    slurm_log_dir = hdf5_path.parent / ".slurm"
    slurm_log_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Pre-stage setup files (params.json, plus macroscale_in for macro)
    # ------------------------------------------------------------------
    staging_data_dir = staging_dir / "data" / run_code
    staging_data_dir.mkdir(parents=True, exist_ok=True)

    # Build a setup runner via from_hdf5 to call _write_setup_files.
    # We import lazily to avoid pulling in numpy etc. at module load.
    import importlib
    runner_module = importlib.import_module(spec.runner_module)
    runner_class = getattr(runner_module, spec.runner_class)
    setup_kwargs = {"out_file_code": spec.out_file_code}
    if spec.in_file_code is not None:
        setup_kwargs["in_file_code"] = spec.in_file_code
    fm_setup = runner_class.from_hdf5(
        hdf5_path, str(executable), **setup_kwargs,
    )
    fm_setup._write_setup_files(staging_data_dir)

    # Pre-stage the executable so no array task needs to cp it — avoids a
    # race where multiple tasks simultaneously try to create the same file.
    if spec.prestage_executable:
        shutil.copy2(executable, staging_dir)

    # Snapshot src/lysis/ so every task (master + array) imports the
    # frozen copy via PYTHONPATH, immune to source edits made while the
    # job is queued or running.
    pythonpath = _snapshot_lysis_src(repo_root, staging_dir)

    # ------------------------------------------------------------------
    # Write array job script
    # ------------------------------------------------------------------
    array_script = generate_array_script(
        spec, staging_dir, run_code, hdf5_path, executable,
        partition=partition, fast_tmp_root=fast_tmp_root,
        slurm_log_dir=slurm_log_dir,
        compiler_module=compiler_module,
        sbatch_overrides=sbatch_overrides,
        historical_binary_attrs=historical_binary_attrs,
        pythonpath=pythonpath,
    )
    array_path = staging_dir / f"lysis-{spec.scale}-array__{run_code}.sh"
    array_path.write_text(array_script)
    array_path.chmod(0o755)

    # ------------------------------------------------------------------
    # Write master.py
    # ------------------------------------------------------------------
    master_py_content = _generate_array_master_py(
        spec, staging_dir, hdf5_path, run_code, executable.name, keep_tmpdir,
        historical_binary_attrs=historical_binary_attrs,
    )
    master_py_path = staging_dir / f"lysis-{spec.scale}-master__{run_code}.py"
    master_py_path.write_text(master_py_content)

    # ------------------------------------------------------------------
    # Write master bash wrapper and submit
    # ------------------------------------------------------------------
    sbatch_opts = {
        "job-name": f"lysis-{spec.scale}-master__{run_code}",
        "out": str(slurm_log_dir / f"lysis-{spec.scale}-master__{run_code}.out"),
        "nodes": 1,
        "mem": 3096,
        "ntasks": 1,
        "cpus-per-task": 1,
        "exclusive=user": "",
    }
    if partition:
        sbatch_opts["partition"] = partition
    sbatch_opts = apply_sbatch_overrides(sbatch_opts, sbatch_overrides)

    master_sh_content = gs.scripts.plain(
        [
            _env_preamble(repo_root, compiler_module, pythonpath=pythonpath),
            f'{_uv_python_prefix(repo_root)} "{master_py_path}"',
        ],
        **sbatch_opts,
    )
    master_sh_path = staging_dir / f"lysis-{spec.scale}-master__{run_code}.sh"
    master_sh_path.write_text(master_sh_content)
    master_sh_path.chmod(0o755)

    return gs.sbatch([str(master_sh_path)])


# ---------------------------------------------------------------------------
# Microscale public API
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
    compiler_module: str = DEFAULT_COMPILER_MODULE,
    sbatch_overrides: Optional[Mapping[str, Optional[str]]] = None,
    historical_binary_attrs: Optional[dict] = None,
    pythonpath: "Path | str | None" = None,
) -> str:
    """Generate a Slurm bash script for a single child microscale Fortran job.

    Used by the legacy single-child path (``submit_micro_slurm_job``
    with ``num_children=None``).

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

    :param staging_dir: Path to the shared staging directory.
    :type staging_dir: Path or str
    :param run_code: Run identifier.
    :type run_code: str
    :param hdf5_path: Full path to the run's ``.h5`` file.
    :type hdf5_path: Path or str
    :param executable: Path to the compiled Fortran microscale binary.
    :type executable: Path or str
    :param out_code: Output file code suffix, defaults to ``""``.
    :type out_code: str, optional
    :param partition: Slurm partition for ``#SBATCH --partition``.
    :type partition: str, optional
    :param fast_tmp_root: Root directory for fast node-local scratch.
    :type fast_tmp_root: str, optional
    :param slurm_log_dir: Directory for Slurm ``.out`` logs.  Defaults to
        ``hdf5_path.parent / ".slurm"``.
    :type slurm_log_dir: Path or str, optional
    :param compiler_module: LMod module spec for the Fortran compiler
        runtime, baked into the generated script.  Defaults to
        :data:`DEFAULT_COMPILER_MODULE`.
    :type compiler_module: str, optional
    :param sbatch_overrides: User overrides for the child's ``#SBATCH``
        header (shape returned by :func:`parse_sbatch_tokens`).  Defaults
        to ``None``.
    :type sbatch_overrides: Mapping[str, str or None], optional
    :param historical_binary_attrs: When non-``None``, indicates the
        binary was built from a historical commit (see
        :func:`~lysis.execution.historical_build.build_historical_binary`).
        The generated ``from_hdf5(...)`` call receives
        ``skip_binary_verification=True`` so ``_verify_binary_version``
        does not crash on the intentional binary↔source mismatch.  The
        dict itself is not baked into the child — the master script
        owns the actual provenance stamping.
    :type historical_binary_attrs: dict or None
    :param pythonpath: When set, baked into the generated script as
        ``export PYTHONPATH=...`` so the child imports ``lysis`` from
        that path rather than the live editable install.  ``None``
        (default) leaves ``PYTHONPATH`` untouched.
    :type pythonpath: Path or str, optional
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

    repo_root = _repo_root()
    env_preamble = _env_preamble(
        repo_root, compiler_module, pythonpath=pythonpath
    )
    py = _uv_python_prefix(repo_root)

    sbatch_opts = {
        "job-name": f"lysis-micro-child__{run_code}",
        "out": str(slurm_log_dir / f"lysis-micro-child__{run_code}.out"),
        "nodes": 1,
        "mem": 3096,
        "ntasks": 1,
        "cpus-per-task": 1,
        "exclusive=user": "",
    }
    if partition:
        sbatch_opts["partition"] = partition
    sbatch_opts = apply_sbatch_overrides(sbatch_opts, sbatch_overrides)

    skip_verify_kwarg = (
        ", skip_binary_verification=True"
        if historical_binary_attrs is not None
        else ""
    )

    if fast_tmp_root is None:
        # ------------------------------------------------------------------
        # Single-tier: Fortran writes directly to the shared staging dir.
        # The binary is pre-staged by submit_micro_slurm_job before this
        # script runs, so no runtime cp is needed — and would in fact be
        # broken under --fortran-commit, which removes the build dir as
        # soon as sbatch returns.
        # ------------------------------------------------------------------
        setup = f"""\
# Setup — create output directory
staging_datadir="{staging_dir}/data/{run_code}"
mkdir -p "${{staging_datadir}}" """

        execute = f"""\
# Execute Fortran binary via lysis
{env_preamble}
{py} -c "
from pathlib import Path
from lysis.execution.fortran_micro import FortranMicro
fm = FortranMicro.from_hdf5('{hdf5_path}', '{staging_dir}/{binary_name}', out_file_code='{out_code}'{skip_verify_kwarg})
fm.exec_in_workdir(Path('{staging_dir}'))
" """

        sections = [setup, execute]

    else:
        # ------------------------------------------------------------------
        # Two-tier: Fortran → fast local dir, then mv to shared staging dir.
        # Source the binary from the pre-staged copy in staging_dir, not
        # the caller's --executable path: under --fortran-commit the build
        # dir is already gone by the time this script runs.
        # ------------------------------------------------------------------
        setup = f"""\
# Setup — create local fast dir and shared staging dir
local_work_dir=$(mktemp -d -p "{fast_tmp_root}")
local_datadir="${{local_work_dir}}/data/{run_code}"
staging_datadir="{staging_dir}/data/{run_code}"
mkdir -p "${{local_datadir}}"
mkdir -p "${{staging_datadir}}"
cp "{staging_dir}/{binary_name}" "${{local_work_dir}}/" """

        execute = f"""\
# Execute Fortran binary into local fast storage
{env_preamble}
{py} -c "
from pathlib import Path
from lysis.execution.fortran_micro import FortranMicro
fm = FortranMicro.from_hdf5('{hdf5_path}', '${{local_work_dir}}/{binary_name}', out_file_code='{out_code}'{skip_verify_kwarg})
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


def _micro_spec(num_children: int, out_code: str) -> _SlurmJobSpec:
    """Build the microscale array-mode dispatch spec."""
    return _SlurmJobSpec(
        scale="micro",
        runner_module="lysis.execution.fortran_micro",
        runner_class="FortranMicro",
        num_array_tasks=num_children,
        needs_concat_step=True,
        nfs_wait_seconds=60,
        prestage_executable=True,
        out_file_code=out_code,
        in_file_code=None,
        num_children=num_children,
    )


def generate_micro_array_script(
    staging_dir: "Path | str",
    run_code: str,
    hdf5_path: "Path | str",
    executable: "Path | str",
    num_children: int,
    out_code: str = "",
    *,
    partition: Optional[str] = None,
    fast_tmp_root: Optional[str] = None,
    slurm_log_dir: Optional["Path | str"] = None,
    compiler_module: str = DEFAULT_COMPILER_MODULE,
    sbatch_overrides: Optional[Mapping[str, Optional[str]]] = None,
    historical_binary_attrs: Optional[dict] = None,
) -> str:
    """Generate a Slurm array job script for the microscale Fortran simulation.

    Thin wrapper over :func:`generate_array_script` with a micro
    :class:`_SlurmJobSpec`.  See that function for the bash script
    contract; micro-specific notes:

    * Each array task handles ``micro_simulations / num_children`` of the
      total simulations (with the lowest-indexed tasks absorbing the
      remainder, partitioned by
      :meth:`~lysis.execution.fortran.FortranRunner.exec_command`).
    * Each task writes flat ``__NN``-suffixed files into
      ``data/{run_code}/``; the master concatenates them after all
      tasks complete.

    :param num_children: Number of array tasks (sets
        ``--array=0-{num_children-1}`` and is threaded into
        :meth:`~lysis.execution.fortran_micro.FortranMicro.from_hdf5` so
        the partition contract is honoured).
    :type num_children: int
    :param out_code: Output file code suffix, defaults to ``""``.
    :type out_code: str, optional
    """
    spec = _micro_spec(num_children, out_code)
    return generate_array_script(
        spec, staging_dir, run_code, hdf5_path, executable,
        partition=partition, fast_tmp_root=fast_tmp_root,
        slurm_log_dir=slurm_log_dir,
        compiler_module=compiler_module,
        sbatch_overrides=sbatch_overrides,
        historical_binary_attrs=historical_binary_attrs,
    )


def submit_micro_slurm_job(
    hdf5_path: "Path | str",
    executable: "Path | str",
    *,
    staging_root: Optional["Path | str"] = None,
    partition: Optional[str] = None,
    fast_tmp_root: Optional[str] = None,
    keep_tmpdir: bool = False,
    out_code: str = "",
    num_children: Optional[int] = None,
    compiler_module: str = DEFAULT_COMPILER_MODULE,
    sbatch_overrides: Optional[Mapping[str, Optional[str]]] = None,
    historical_binary_attrs: Optional[dict] = None,
) -> int:
    """Submit the full HDF5-integrated microscale workflow as a Slurm master job.

    Two execution modes:

    * **Single-child (legacy)**: ``num_children=None`` (default).  Submits
      one child Fortran job that runs all microscale simulations
      sequentially, then imports.  Bit-for-bit reproducible with
      pre-existing single-task artifacts.
    * **Array (opt-in)**: ``num_children >= 1``.  Submits a Slurm array
      of ``num_children`` tasks; the master concatenates their per-task
      outputs and imports.  Each task runs
      ``micro_simulations // num_children`` simulations (with the
      lowest-indexed tasks absorbing the remainder).

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
    :param partition: Slurm partition for both master and child/array jobs.
    :type partition: str, optional
    :param fast_tmp_root: Root directory for fast node-local scratch (opt-in
        two-tier storage).  When ``None`` (default) tasks write directly
        to the shared staging directory.
    :type fast_tmp_root: str, optional
    :param keep_tmpdir: If ``True``, preserve the staging directory after the
        master job completes (useful for debugging).
    :type keep_tmpdir: bool, optional
    :param out_code: Output file code suffix for the Fortran binary,
        defaults to ``""``.
    :type out_code: str, optional
    :param num_children: When set to a positive integer, switches to the
        Slurm array execution mode and partitions ``micro_simulations``
        across that many tasks.  ``None`` (default) selects the legacy
        single-child path.
    :type num_children: int, optional
    :param compiler_module: LMod module spec for the Fortran compiler
        runtime, baked into every generated script (master, array tasks,
        and the legacy single child).  Defaults to
        :data:`DEFAULT_COMPILER_MODULE`.
    :type compiler_module: str, optional
    :param sbatch_overrides: User overrides for the ``#SBATCH`` header,
        applied to BOTH the master and child/array task scripts (shape
        returned by :func:`parse_sbatch_tokens`).  Defaults to ``None``.
    :type sbatch_overrides: Mapping[str, str or None], optional
    :param historical_binary_attrs: Pre-computed binary-provenance dict
        from
        :func:`~lysis.tools.provenance.gather_historical_binary_provenance`
        for the historical-build workflow.  Threaded through to the
        master Python script (in both array and legacy single-child
        paths) so the master stamps it in place of the default
        ``<binary> --version`` query.  ``None`` (default) preserves
        normal-mode provenance.
    :type historical_binary_attrs: dict, optional
    :return: Master Slurm job ID.
    :rtype: int
    :raises subprocess.CalledProcessError: If ``sbatch`` fails.
    :raises ValueError: If ``num_children`` is set and is either less than
        1 or greater than ``micro_simulations``.
    """
    from lysis.execution.fortran_micro import FortranMicro  # noqa: PLC0415

    hdf5_path = Path(hdf5_path).resolve()
    executable = Path(executable).resolve()
    run_code = hdf5_path.stem

    # Resolve staging root once; both branches use it.
    staging_root_dir = (
        Path(staging_root) if staging_root is not None else hdf5_path.parent
    )

    # ------------------------------------------------------------------
    # Array path: delegate to the unified array-mode dispatch.
    # ------------------------------------------------------------------
    if num_children is not None:
        if num_children < 1:
            raise ValueError(
                f"num_children must be >= 1 when set, got {num_children}"
            )
        from lysis.config.run import Run  # noqa: PLC0415
        run_check = Run(str(hdf5_path.parent), run_code=run_code)
        run_check.load_params_from_hdf5()
        total_sims = int(run_check.micro_params.micro_simulations)
        if num_children > total_sims:
            raise ValueError(
                f"num_children ({num_children}) exceeds micro_simulations "
                f"({total_sims}) for run {run_code!r}; cannot split "
                f"{total_sims} simulations across {num_children} tasks."
            )
        spec = _micro_spec(num_children, out_code)
        return submit_slurm_job(
            spec, hdf5_path, executable, staging_root_dir,
            partition=partition, fast_tmp_root=fast_tmp_root,
            keep_tmpdir=keep_tmpdir,
            compiler_module=compiler_module,
            sbatch_overrides=sbatch_overrides,
            historical_binary_attrs=historical_binary_attrs,
        )

    # ------------------------------------------------------------------
    # Legacy single-child path.
    # ------------------------------------------------------------------

    # Sync the project venv against uv.lock before any Slurm tasks fire.
    repo_root = _repo_root()
    _ensure_env_synced(repo_root)

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

    # Pre-stage the binary into staging_dir so the master can resolve it
    # at import time (its `--version` output supplies binary provenance).
    # Mirrors the pre-stage step in submit_slurm_job.  Always present so
    # two-tier (fast_tmp_root) doesn't leave the master without a binary.
    shutil.copy2(executable, staging_dir)

    # Snapshot src/lysis/ so the master + child import the frozen copy
    # via PYTHONPATH, immune to source edits made while the job runs.
    pythonpath = _snapshot_lysis_src(repo_root, staging_dir)

    # ------------------------------------------------------------------
    # Write child script
    # ------------------------------------------------------------------
    child_script = generate_micro_child_script(
        staging_dir, run_code, hdf5_path, executable, out_code,
        partition=partition, fast_tmp_root=fast_tmp_root,
        slurm_log_dir=slurm_log_dir,
        compiler_module=compiler_module,
        sbatch_overrides=sbatch_overrides,
        historical_binary_attrs=historical_binary_attrs,
        pythonpath=pythonpath,
    )
    child_path = staging_dir / f"lysis-micro-child__{run_code}.sh"
    child_path.write_text(child_script)
    child_path.chmod(0o755)

    # ------------------------------------------------------------------
    # Write master.py
    # ------------------------------------------------------------------
    master_py_content = _generate_micro_master_py(
        staging_dir, hdf5_path, run_code, out_code, executable.name, keep_tmpdir,
        historical_binary_attrs=historical_binary_attrs,
    )
    master_py_path = staging_dir / f"lysis-micro-master__{run_code}.py"
    master_py_path.write_text(master_py_content)

    # ------------------------------------------------------------------
    # Write master bash wrapper and submit
    # ------------------------------------------------------------------
    sbatch_opts = {
        "job-name": f"lysis-micro-master__{run_code}",
        "out": str(slurm_log_dir / f"lysis-micro-master__{run_code}.out"),
        "nodes": 1,
        "mem": 3096,
        "ntasks": 1,
        "cpus-per-task": 1,
        "exclusive=user": "",
    }
    if partition:
        sbatch_opts["partition"] = partition
    sbatch_opts = apply_sbatch_overrides(sbatch_opts, sbatch_overrides)

    master_sh_content = gs.scripts.plain(
        [
            _env_preamble(repo_root, compiler_module, pythonpath=pythonpath),
            f'{_uv_python_prefix(repo_root)} "{master_py_path}"',
        ],
        **sbatch_opts,
    )
    master_sh_path = staging_dir / f"lysis-micro-master__{run_code}.sh"
    master_sh_path.write_text(master_sh_content)
    master_sh_path.chmod(0o755)

    return gs.sbatch([str(master_sh_path)])


# ---------------------------------------------------------------------------
# Macroscale public API (thin wrappers over the unified array internals)
# ---------------------------------------------------------------------------


def _macro_spec(n_sims: int, in_code: str, out_code: str) -> _SlurmJobSpec:
    """Build the macroscale array-mode dispatch spec."""
    return _SlurmJobSpec(
        scale="macro",
        runner_module="lysis.execution.fortran_macro",
        runner_class="FortranMacro",
        num_array_tasks=n_sims,
        needs_concat_step=False,
        nfs_wait_seconds=60,
        prestage_executable=True,
        out_file_code=out_code,
        in_file_code=in_code,
        num_children=None,
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
    compiler_module: str = DEFAULT_COMPILER_MODULE,
    sbatch_overrides: Optional[Mapping[str, Optional[str]]] = None,
    historical_binary_attrs: Optional[dict] = None,
) -> str:
    """Generate a Slurm array job script for the macroscale Fortran simulation.

    Thin wrapper over :func:`generate_array_script` with a macro
    :class:`_SlurmJobSpec`.  See that function for the bash script
    contract; macro-specific notes:

    * Each array task handles one simulation indexed by
      ``SLURM_ARRAY_TASK_ID`` and writes to its own
      ``data/{run_code}/{sim:02}/`` subdirectory.
    * Setup files pre-staged in ``{staging_dir}/data/{run_code}/`` are
      symlinked into each per-sim subdir by
      :class:`~lysis.execution.fortran_macro.FortranMacro`.

    :param n_sims: Total number of simulations (sets ``--array=0-{n_sims-1}``).
    :type n_sims: int
    :param in_code: Input file code suffix, defaults to ``""``.
    :type in_code: str, optional
    :param out_code: Output file code suffix, defaults to ``""``.
    :type out_code: str, optional
    """
    spec = _macro_spec(n_sims, in_code, out_code)
    return generate_array_script(
        spec, staging_dir, run_code, hdf5_path, executable,
        partition=partition, fast_tmp_root=fast_tmp_root,
        slurm_log_dir=slurm_log_dir,
        compiler_module=compiler_module,
        sbatch_overrides=sbatch_overrides,
        historical_binary_attrs=historical_binary_attrs,
    )


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
    compiler_module: str = DEFAULT_COMPILER_MODULE,
    sbatch_overrides: Optional[Mapping[str, Optional[str]]] = None,
    historical_binary_attrs: Optional[dict] = None,
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
    :param compiler_module: LMod module spec for the Fortran compiler
        runtime, baked into the generated master and array task scripts.
        Defaults to :data:`DEFAULT_COMPILER_MODULE`.
    :type compiler_module: str, optional
    :param sbatch_overrides: User overrides for the ``#SBATCH`` header,
        applied to BOTH the master and array task scripts (shape returned
        by :func:`parse_sbatch_tokens`).  Defaults to ``None``.
    :type sbatch_overrides: Mapping[str, str or None], optional
    :param historical_binary_attrs: Pre-computed binary-provenance dict
        from
        :func:`~lysis.tools.provenance.gather_historical_binary_provenance`
        for the historical-build workflow.  Threaded through to the
        master Python script so the master stamps it in place of the
        default ``<binary> --version`` query.  ``None`` (default)
        preserves normal-mode provenance.
    :type historical_binary_attrs: dict, optional
    :return: Master Slurm job ID.
    :rtype: int
    :raises subprocess.CalledProcessError: If ``sbatch`` fails.
    :raises ValueError: If the HDF5 file does not contain ``macro_params``.
    """
    from lysis.execution.fortran_macro import FortranMacro  # noqa: PLC0415

    hdf5_path = Path(hdf5_path).resolve()
    executable = Path(executable).resolve()

    staging_root_dir = (
        Path(staging_root) if staging_root is not None else hdf5_path.parent
    )

    # Read n_sims from the run's macro_params (peek before submit_slurm_job
    # builds its own setup runner — submit_slurm_job will rebuild and call
    # _write_setup_files itself).
    fm_peek = FortranMacro.from_hdf5(
        hdf5_path, str(executable), in_file_code=in_code, out_file_code=out_code,
    )
    n_sims = fm_peek.run.macro_params.macro_simulations

    spec = _macro_spec(n_sims, in_code, out_code)
    return submit_slurm_job(
        spec, hdf5_path, executable, staging_root_dir,
        partition=partition, fast_tmp_root=fast_tmp_root,
        keep_tmpdir=keep_tmpdir,
        compiler_module=compiler_module,
        sbatch_overrides=sbatch_overrides,
        historical_binary_attrs=historical_binary_attrs,
    )
