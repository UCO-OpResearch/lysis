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

import shlex
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from string import Template
from typing import Dict, List, Mapping, Optional

# ``GooseSLURM`` is an optional HPC dependency, only needed when actually
# submitting or polling Slurm jobs.  Guard the import the same way the rest of
# the package guards optional deps (cf. ``cupy`` in ``lysis/__init__.py`` and
# ``matplotlib`` in ``lysis/tools/__init__.py``) so that ``import lysis``
# succeeds in environments without it — e.g. the docs build.  The submission
# helpers below check :func:`_require_gs` before using ``gs``.
try:
    import GooseSLURM as gs
except ImportError:  # pragma: no cover - exercised only off-cluster
    gs = None


def _require_gs():
    """Return the :mod:`GooseSLURM` module, or raise a clear error if absent.

    :return: The imported ``GooseSLURM`` module.
    :rtype: module
    :raises ImportError: If ``GooseSLURM`` is not installed.
    """
    if gs is None:  # pragma: no cover - exercised only off-cluster
        raise ImportError(
            "GooseSLURM is required for Slurm job submission/polling but is "
            "not installed. Install it (HPC environments only) to use the "
            "lysis.tools.slurm submission helpers."
        )
    return gs


# ---------------------------------------------------------------------------
# Environment setup for generated Slurm scripts
# ---------------------------------------------------------------------------
#
# Generated scripts no longer depend on per-user dotfiles (``~/.bashrc`` /
# ``~/lysis.sh``).  Instead, every script emits a self-contained preamble
# that loads the requested LMod modules and invokes Python through the
# project's ``uv``-managed environment.  The submitting user's repo root is
# baked in at script-generation time from ``lysis.__file__``.

#: Default space-separated list of LMod modules providing the Fortran
#: binary's runtime libraries.  Assumed to be available on the cluster's
#: LMod stack for all users.  Callers can override per-job via the
#: ``modules`` keyword argument on the public submit/generate functions
#: (exposed as ``--modules`` on the ``lysis run-micro`` and
#: ``lysis run-macro`` CLI commands).
DEFAULT_MODULES: str = "intel-compilers/2023"


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
    modules: str,
    *,
    pythonpath: "Path | str | None" = None,
) -> str:
    """Return the bash preamble that prepares a generated script's environment.

    Loads the requested LMod modules (which should include the Fortran
    compiler runtime required by the Fortran binary) and puts
    ``~/.local/bin`` on PATH so ``uv`` is discoverable on compute nodes.
    Defensive against non-login shells by sourcing ``lmod.sh`` when the
    ``module`` function isn't already defined.

    When *pythonpath* is supplied, an ``export PYTHONPATH=...`` line is
    appended.  ``uv run`` inherits environment variables from the calling
    shell, and PYTHONPATH entries come before site-packages on ``sys.path``,
    so a snapshot directory containing ``lysis/`` will shadow the editable
    install — pinning the Python source for the duration of the job.

    :param repo_root: Absolute path to the lysis repo root (currently
        unused inside the preamble itself; included for future extensions
        like ``cd`` into the repo).
    :type repo_root: Path
    :param modules: Space-separated list of LMod module specs to load (e.g.
        ``"intel-compilers/2023"`` or
        ``"intel-compilers/2024 SciPy-bundle/2023.07"``).  Forwarded
        verbatim to ``module load``.  See :data:`DEFAULT_MODULES` for the
        project default.
    :type modules: str
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
        f"module purge && module load {modules}",
        'export PATH="$HOME/.local/bin:$PATH"',
    ]
    if pythonpath is not None:
        lines.append(f'export PYTHONPATH="{pythonpath}"')
    return "\n".join(lines)


def _log_header(
    *,
    scale: str,
    role: str,
    run_code: str,
    experiment: Optional[str] = None,
    pythonpath: "Path | str | None" = None,
) -> str:
    """Return the bash block that prints a job header, then enables strict tracing.

    Emitted in every generated Slurm script (master, array task, legacy
    single child) immediately *after* :func:`_env_preamble`.  The block has
    two parts:

    1. A fenced **identifying header** written to stdout (so it lands in the
       job's ``.out`` log).  *Submit-time* fields baked in here (*scale*,
       *role*, *run_code*, *experiment*, the generator's ``lysis.__version__``,
       and the pinned *pythonpath* snapshot) are emitted via ``printf`` with
       each line passed as a ``shlex.quote``-escaped shell literal, so bash
       performs **no** parameter/command substitution on them — a
       caller-controlled value such as ``experiment="EXP-$(cmd)"`` is printed
       verbatim, never executed.  *Runtime* fields are emitted via a separate
       heredoc whose body contains only our own ``${SLURM_*}``/``$(...)``
       expansions (cluster, host, partition, node list, job/array ids,
       allocated cpus/mem, submit dir, start time) — no caller input.
    2. ``set -euo pipefail`` followed by ``set -x``.

    The header is printed **before** strict mode is enabled so its optional
    ``${SLURM_*}`` reads cannot trip ``set -u`` (each is additionally guarded
    with a ``${VAR:-…}`` default for readable output); ``set -x`` then traces
    every subsequent command into the same log.

    :param scale: ``"micro"`` or ``"macro"``.
    :type scale: str
    :param role: Job role for the header, e.g. ``"master"``,
        ``"array task"``, ``"micro child"``.
    :type role: str
    :param run_code: Run identifier (HDF5 file stem).
    :type run_code: str
    :param experiment: Best-effort experiment identifier; ``None`` (default)
        renders as ``(unknown)`` (no first-class experiment concept yet).
    :type experiment: str, optional
    :param pythonpath: The pinned ``python_src`` snapshot directory baked
        onto ``PYTHONPATH``; ``None`` (default) renders as
        ``(editable install)``.
    :type pythonpath: Path or str, optional
    :return: Bash snippet (header + strict-mode lines), no trailing newline.
    :rtype: str
    """
    import lysis  # noqa: PLC0415 — avoids pulling lysis at module import

    version = getattr(lysis, "__version__", "(unknown)")
    experiment_field = experiment if experiment else "(unknown)"
    pythonpath_field = (
        str(pythonpath) if pythonpath is not None else "(editable install)"
    )
    # Submit-time fields: emitted via ``printf`` with each line shell-quoted,
    # so caller-controlled values (e.g. *experiment*) cannot inject command or
    # parameter substitution into the unprivileged-looking header.
    baked_lines = [
        "============================================================",
        "lysis slurm task",
        "------------------------------------------------------------",
        f"scale       : {scale}",
        f"role        : {role}",
        f"run         : {run_code}",
        f"experiment  : {experiment_field}",
        f"generator   : lysis {version}",
        f"pythonpath  : {pythonpath_field}",
    ]
    baked_block = "printf '%s\\n' \\\n" + " \\\n".join(
        "    " + shlex.quote(line) for line in baked_lines
    )
    # Runtime fields: intentionally expanded by bash on the compute node.  The
    # heredoc body contains only our own ${SLURM_*}/$(...) — never caller input.
    runtime_block = (
        "cat <<RUNTIME\n"
        "job name    : ${SLURM_JOB_NAME:-(unknown)}\n"
        "cluster     : ${SLURM_CLUSTER_NAME:-(unknown)}\n"
        "partition   : ${SLURM_JOB_PARTITION:-(unknown)}\n"
        "host        : $(hostname -f 2>/dev/null || hostname)\n"
        "node list   : ${SLURM_JOB_NODELIST:-(unknown)}\n"
        "user        : ${USER:-$(id -un)}\n"
        "submit dir  : ${SLURM_SUBMIT_DIR:-(unknown)}\n"
        "job id      : ${SLURM_JOB_ID:-(unknown)}\n"
        "array       : ${SLURM_ARRAY_JOB_ID:-n/a}_${SLURM_ARRAY_TASK_ID:-n/a}\n"
        "cpus        : ${SLURM_CPUS_ON_NODE:-(unknown)}\n"
        "mem (MB)    : ${SLURM_MEM_PER_NODE:-(unknown)}\n"
        "started     : $(date -Iseconds)\n"
        "============================================================\n"
        "RUNTIME"
    )
    return (
        "# Identifying header (to stdout) printed BEFORE strict mode so its\n"
        "# optional ${SLURM_*} reads can't abort the job under ``set -u``;\n"
        "# ``set -x`` then traces the rest of the script into the same .out\n"
        "# log.  Submit-time fields are shell-quoted (printf) to prevent\n"
        "# command injection; runtime ${SLURM_*} fields are expanded on the\n"
        "# compute node (cat heredoc).\n"
        f"{baked_block}\n"
        f"{runtime_block}\n"
        "\n"
        "# Strict mode + command tracing for the remainder of the script.\n"
        "set -euo pipefail\n"
        "set -x"
    )


def _guarded_glob_xfer(
    verb: str,
    glob: str,
    dest: str,
    *,
    array_name: str,
    empty_msg: str,
) -> str:
    """Return a bash block that ``cp``/``mv``-es a glob's matches, or fails loud.

    Generated Slurm scripts run under ``set -euo pipefail``.  A bare
    ``cp``/``mv`` of ``${dir}/*`` is unsafe there: with ``nullglob`` off an
    empty *dir* leaves the ``*`` unexpanded, so the command aborts the task
    with a cryptic ``cannot stat '…/*'``.  This helper collects the matches
    into an array under ``nullglob``, then either prints a clear,
    lysis-branded error and exits non-zero when there are none (so an
    unexpectedly empty source is never *silently* skipped), or runs *verb* on
    the collected files.

    Used for every glob-based transfer in the two-tier (fast-tmp) paths — the
    staged-setup-file copy and the micro per-task output move.  Named-path
    transfers (the macro per-sim ``mv`` and the single-file binary ``cp``) do
    not need it: they already fail clearly on a missing operand.

    :param verb: ``"cp"`` or ``"mv"``.
    :type verb: str
    :param glob: Source glob, already quoted for bash (e.g.
        ``'"${local_datadir}"/*'``).
    :type glob: str
    :param dest: Destination, already quoted for bash (e.g.
        ``'"${staging_datadir}/"'``).
    :type dest: str
    :param array_name: Bash array variable to collect matches into; must be
        unique within the script (e.g. ``"staged_files"``, ``"output_files"``).
    :type array_name: str
    :param empty_msg: Reason shown in the "nothing matched" error; may
        reference bash vars (e.g. ``"no output in ${local_datadir}"``).
    :type empty_msg: str
    :return: Multi-line bash snippet, no trailing newline.
    :rtype: str
    """
    return (
        "shopt -s nullglob\n"
        f"{array_name}=({glob})\n"
        "shopt -u nullglob\n"
        f'if [ "${{#{array_name}[@]}}" -eq 0 ]; then\n'
        f'    echo "lysis: {empty_msg}" >&2\n'
        "    exit 1\n"
        "fi\n"
        f'{verb} "${{{array_name}[@]}}" {dest}'
    )


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
HISTORICAL_BACKEND_ATTRS = $historical_backend_attrs

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
    skip_binary_verification=(HISTORICAL_BACKEND_ATTRS is not None),
    historical_backend_attrs=HISTORICAL_BACKEND_ATTRS,
)
fm.import_results(
    data_dir,
    keep_on_failure=True,
    keep_tmpdir=KEEP_TMPDIR,
    dispatcher_log_dir=STAGING_DIR,
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
    historical_backend_attrs: Optional[dict] = None,
) -> str:
    """Return the content of the single-child micro master Python script.

    :param historical_backend_attrs: Pre-computed binary-provenance dict
        from
        :func:`~lysis.tools.provenance.gather_historical_backend_provenance`
        for the historical-build workflow.  Baked into the master script
        as a literal so the master can stamp it in place of the default
        ``<binary> --version`` query.  ``None`` (default) preserves the
        normal-mode behaviour.
    :type historical_backend_attrs: dict or None
    """
    return _MICRO_MASTER_PY_TEMPLATE.substitute(
        staging_dir=repr(str(staging_dir)),
        hdf5_path=repr(str(hdf5_path)),
        run_code=repr(run_code),
        file_code=repr(file_code),
        binary_name=repr(binary_name),
        keep_tmpdir=repr(keep_tmpdir),
        historical_backend_attrs=repr(historical_backend_attrs),
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
    :ivar backend: ``"fortran"`` (default) or ``"python"``.  Selects the
        body of the array-task bash script (Fortran binary vs. in-process
        :class:`~lysis.execution.python_macro.PythonMacro`) and of the
        master Python script's post-poll block (Fortran imports per-task
        outputs into HDF5; Python tasks write to HDF5 themselves and the
        master has nothing left to do).  Python backends ignore
        ``out_file_code``, ``in_file_code``, ``num_children``,
        ``needs_concat_step`` and ``prestage_executable``.
    :vartype backend: str
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
    backend: str = "fortran"
    #: Forward the legacy ``seed_scheme="direct"`` path into each array
    #: task's ``from_hdf5(...)`` call (micro only).  Defaults to False.
    direct: bool = False


#: Template for the master Python script in array mode (both backends).
#: Substitution keys use ``$name`` syntax so Python braces in the
#: template body are untouched.  The post-poll body is supplied per-backend
#: via the ``$import_block`` substitution (Fortran: instantiate runner and
#: call ``import_results``; Python: print a status line, since the array
#: task wrote to HDF5 itself).
_ARRAY_MASTER_PY_TEMPLATE = Template('''\
#!/usr/bin/env python3
"""Master Slurm job: submit array tasks, poll for completion, $post_poll_blurb."""
import time
import shutil
from pathlib import Path

import GooseSLURM as gs

STAGING_DIR = Path($staging_dir)
HDF5_PATH = Path($hdf5_path)
RUN_CODE = $run_code
KEEP_TMPDIR = $keep_tmpdir
NUM_ARRAY_TASKS = $num_array_tasks
NFS_WAIT_SECONDS = $nfs_wait_seconds

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
# Per-backend post-poll body (Fortran: import per-task outputs into HDF5;
# Python: nothing to do — array task already wrote directly to HDF5).
# ---------------------------------------------------------------------------
$import_block

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
    binary_name: Optional[str],
    keep_tmpdir: bool,
    historical_backend_attrs: Optional[dict] = None,
) -> str:
    """Return the master Python script for an array-mode dispatch.

    The post-poll body is assembled per-backend:

    * ``"fortran"``: imports the runner class, instantiates it, optionally
      concatenates per-task outputs, and calls ``import_results`` to merge
      the per-task Fortran output into the HDF5 file.
    * ``"python"``: prints a status line.  The array task ran
      :meth:`~lysis.execution.python_macro.PythonRunner.execute` itself
      and wrote directly to HDF5, so the master has nothing left to do
      before the staging-dir cleanup.

    :param binary_name: Basename of the pre-staged Fortran binary inside
        the staging dir.  Required for ``backend="fortran"``; ignored
        (and may be ``None``) for ``backend="python"``.
    :type binary_name: str or None
    :param historical_backend_attrs: Pre-computed backend-provenance dict
        from
        :func:`~lysis.tools.provenance.gather_historical_backend_provenance`
        for the historical-build workflow.  Baked into the master script
        as a literal so the master can stamp it in place of the default
        ``<binary> --version`` query.  ``None`` (default) preserves the
        normal-mode behaviour.  Only meaningful for ``backend="fortran"``.
    :type historical_backend_attrs: dict or None
    """
    if spec.backend == "python":
        import_block = (
            'print("Python task wrote results directly into HDF5.",'
            " flush=True)"
        )
        post_poll_blurb = "exit (HDF5 already written by the array task)"
    else:
        if spec.needs_concat_step:
            concat_block = (
                "fm.concatenate_child_outputs("
                f"data_dir, num_children={spec.num_array_tasks})"
            )
            post_poll_blurb = "concatenate per-task outputs, then import"
        else:
            concat_block = (
                "# (no concatenation step — runner imports per-task subdirs"
                " directly)"
            )
            post_poll_blurb = "import"

        import_block = (
            "from lysis.config.run import Run\n"
            f"from {spec.runner_module} import {spec.runner_class}\n"
            "\n"
            f"OUT_FILE_CODE = {repr(spec.out_file_code)}\n"
            f"BINARY_NAME = {repr(binary_name)}\n"
            f"HISTORICAL_BACKEND_ATTRS = {repr(historical_backend_attrs)}\n"
            "\n"
            "data_dir = STAGING_DIR / \"data\" / RUN_CODE\n"
            "run = Run(str(HDF5_PATH.parent), run_code=RUN_CODE)\n"
            f"fm = {spec.runner_class}(\n"
            "    run=run,\n"
            "    out_file_code=OUT_FILE_CODE,\n"
            "    executable=str(STAGING_DIR / BINARY_NAME),\n"
            "    skip_binary_verification=(HISTORICAL_BACKEND_ATTRS is not None),\n"
            "    historical_backend_attrs=HISTORICAL_BACKEND_ATTRS,\n"
            ")\n"
            f"{concat_block}\n"
            "fm.import_results(\n"
            "    data_dir,\n"
            "    keep_on_failure=True,\n"
            "    keep_tmpdir=KEEP_TMPDIR,\n"
            "    dispatcher_log_dir=STAGING_DIR,\n"
            ")\n"
            'print("Results imported successfully.", flush=True)'
        )

    return _ARRAY_MASTER_PY_TEMPLATE.substitute(
        scale=spec.scale,
        staging_dir=repr(str(staging_dir)),
        hdf5_path=repr(str(hdf5_path)),
        run_code=repr(run_code),
        keep_tmpdir=repr(keep_tmpdir),
        num_array_tasks=spec.num_array_tasks,
        nfs_wait_seconds=spec.nfs_wait_seconds,
        import_block=import_block,
        post_poll_blurb=post_poll_blurb,
    )


def _runner_init_line(
    spec: _SlurmJobSpec,
    hdf5_path: str,
    binary_path: Optional[str],
    *,
    skip_binary_verification: bool = False,
    source_stamp: Optional[tuple] = None,
) -> str:
    """Build the ``from_hdf5(...)`` call used inside the array task ``python -c``.

    Single-line so it composes cleanly inside the bash heredoc.

    For ``spec.backend == "fortran"`` the call carries the binary path,
    in/out file codes, the per-task ``SLURM_ARRAY_TASK_ID``, optionally
    ``num_children`` for the microscale array partition, and (when
    *skip_binary_verification* is True under the historical-build
    workflow) ``skip_binary_verification=True`` so the intentional
    binary↔source mismatch does not crash ``_verify_binary_version`` at
    task start.  When *source_stamp* is provided it is baked into the
    call as ``source_stamp=(commit, dirty)`` so the array task does not
    have to invoke git from the staging directory (which is outside the
    source repository on most HPC layouts).

    For ``spec.backend == "python"`` only ``hdf5_path`` is emitted —
    :meth:`~lysis.execution.python_macro.PythonRunner.from_hdf5` takes
    no other arguments, the single array task runs every simulation in
    series, and there is no binary to verify.

    :param hdf5_path: HDF5 path as it should appear inside the ``python
        -c`` string.  May be a literal absolute path or a ``${...}`` bash
        variable reference (the caller is responsible for picking which).
    :type hdf5_path: str
    :param binary_path: Fortran binary path inside the bash script.
        Required for ``backend="fortran"``; ignored (and may be ``None``)
        for ``backend="python"``.
    :type binary_path: str or None
    """
    if spec.backend == "python":
        return f"{spec.runner_class}.from_hdf5('{hdf5_path}')"
    args = [f"'{hdf5_path}'", f"'{binary_path}'"]
    if spec.in_file_code is not None:
        args.append(f"in_file_code='{spec.in_file_code}'")
    args.append(f"out_file_code='{spec.out_file_code}'")
    args.append("index=int(os.environ['SLURM_ARRAY_TASK_ID'])")
    if spec.num_children is not None:
        args.append(f"num_children={spec.num_children}")
    if spec.direct:
        args.append("direct=True")
    if skip_binary_verification:
        args.append("skip_binary_verification=True")
    if source_stamp is not None and not skip_binary_verification:
        commit, dirty = source_stamp
        args.append(f"source_stamp=('{commit}', '{dirty}')")
    return f"{spec.runner_class}.from_hdf5({', '.join(args)})"


def generate_array_script(
    spec: _SlurmJobSpec,
    staging_dir: "Path | str",
    run_code: str,
    hdf5_path: "Path | str",
    executable: "Path | str | None",
    *,
    partition: Optional[str] = None,
    fast_tmp_root: Optional[str] = None,
    slurm_log_dir: Optional["Path | str"] = None,
    modules: str = DEFAULT_MODULES,
    sbatch_overrides: Optional[Mapping[str, Optional[str]]] = None,
    historical_backend_attrs: Optional[dict] = None,
    pythonpath: "Path | str | None" = None,
    source_stamp: Optional[tuple] = None,
    experiment: Optional[str] = None,
) -> str:
    """Generate a Slurm array job bash script for an array-mode dispatch.

    For ``spec.backend == "fortran"`` each array task indexed by
    ``SLURM_ARRAY_TASK_ID`` calls
    ``{spec.runner_class}.from_hdf5(...).exec_in_workdir(...)``.

    Every generated task begins with an identifying header (see
    :func:`_log_header`) and then enables ``set -euo pipefail`` + ``set -x``,
    so the per-task ``.out`` log records which job/host/allocation ran and
    traces each shell step.  Because Slurm does not substitute ``%a`` in
    ``--job-name`` (only in ``--output``), the array shares a single job
    name and the per-task index is surfaced via that header.

    For the *macroscale* path each task isolates its output in a
    per-simulation subdirectory ``data/{run_code}/{sim:02}/`` (managed by
    :class:`~lysis.execution.fortran_macro.FortranMacro` itself).  For the
    *microscale* array path each task writes flat files with a ``__NN``
    suffix into ``data/{run_code}/``; the master concatenates them after
    all tasks complete.

    For ``spec.backend == "python"`` the single array task calls
    ``{spec.runner_class}.from_hdf5(hdf5).execute()`` directly — the
    Python runner writes every simulation, in series, straight into the
    HDF5 file via its open DataStore (HDF5 has a single-writer
    constraint, so an n-task array is not possible at this scale).
    ``executable`` is ignored and may be ``None``.

    **Single-tier** (default, ``fast_tmp_root=None``):

    * Fortran: writes its output directly into
      ``{staging_dir}/data/{run_code}/``.  The executable is *not* copied
      here — :func:`submit_slurm_job` pre-stages it before any task
      starts (avoids a cp race).
    * Python: writes directly into the original ``hdf5_path`` (typically
      on a shared filesystem).

    **Two-tier** (``fast_tmp_root`` provided):

    * Fortran: writes to a private ``mktemp -d -p {fast_tmp_root}``
      directory on the compute node; the output is moved to the shared
      staging dir before exit.
    * Python: ``cp``\\ s ``hdf5_path`` into the local fast dir, runs the
      Python task against that copy, then ``cp``\\ s the result back
      over the original on success.  A mid-job crash leaves the
      original ``.h5`` in its ``MACRO_EMPTY`` state for clean
      resubmission.

    :param spec: Scale-specific dispatch description.
    :type spec: _SlurmJobSpec
    :param staging_dir: Path to the shared staging directory.
    :type staging_dir: Path or str
    :param run_code: Run identifier.
    :type run_code: str
    :param hdf5_path: Full path to the run's ``.h5`` file.
    :type hdf5_path: Path or str
    :param executable: Path to the compiled Fortran binary.  May be
        ``None`` when ``spec.backend == "python"``.
    :type executable: Path or str or None
    :param partition: Slurm partition for ``#SBATCH --partition``.
    :type partition: str, optional
    :param fast_tmp_root: Root directory for fast node-local scratch.
    :type fast_tmp_root: str, optional
    :param slurm_log_dir: Retained for API compatibility; no longer used
        by this generator.  The array task's ``--output`` now writes into
        *staging_dir* (alongside the generated scripts), so worker
        stdout/stderr — including the identifying header and ``set -x``
        trace — is co-located with the run's scripts and data rather than
        in ``hdf5_path.parent / ".slurm"``.
    :type slurm_log_dir: Path or str, optional
    :param modules: Space-separated list of LMod module specs to load in
        each generated task script (include a Fortran compiler module so
        the binary finds its runtime).  Defaults to
        :data:`DEFAULT_MODULES`.
    :type modules: str, optional
    :param sbatch_overrides: User overrides for the per-task ``#SBATCH``
        header, in the shape returned by :func:`parse_sbatch_tokens`.
        A ``None`` value removes the matching default; any other value
        sets/overrides the key.  Defaults to ``None`` (no overrides).
    :type sbatch_overrides: Mapping[str, str or None], optional
    :param historical_backend_attrs: When non-``None``, indicates the
        binary was built from a historical commit (see
        :func:`~lysis.execution.historical_build.build_historical_binary`).
        The generated ``from_hdf5(...)`` call receives
        ``skip_binary_verification=True`` so ``_verify_binary_version``
        does not crash on the intentional binary↔source mismatch.  The
        dict itself is not baked into the task — the master script
        owns the actual provenance stamping.
    :type historical_backend_attrs: dict or None
    :param pythonpath: When set, baked into the generated script as
        ``export PYTHONPATH=...`` so each task imports ``lysis`` from
        that path rather than the live editable install.  ``None``
        (default) leaves ``PYTHONPATH`` untouched.
    :type pythonpath: Path or str, optional
    :param source_stamp: Pre-resolved ``(commit, dirty)`` for
        ``src/fortran/``.  Baked into the array task's ``from_hdf5``
        call so the binary↔source check on the compute node compares
        against the master's stamp instead of running ``git`` from the
        staging dir (which is outside the repo and would fail).
        Ignored when ``historical_backend_attrs`` is set (the staleness
        check is bypassed in that workflow).  ``None`` (default)
        preserves the in-process git lookup.
    :type source_stamp: tuple[str, str] or None
    :param experiment: Best-effort experiment identifier shown in the log
        header; ``None`` (default) renders as ``(unknown)``.
    :type experiment: str, optional
    :return: Slurm array job bash script text.
    :rtype: str
    :raises ImportError: If ``GooseSLURM`` is not installed.
    """
    _require_gs()
    staging_dir = Path(staging_dir)
    hdf5_path = Path(hdf5_path)
    if executable is not None:
        executable = Path(executable)
        binary_name = executable.name
    else:
        binary_name = None
    skip_binary_verification = historical_backend_attrs is not None

    repo_root = _repo_root()
    # Leading section: env setup (module/PATH/PYTHONPATH) runs first, then the
    # header prints and strict mode + tracing engage for the rest of the task.
    lead = "{}\n{}".format(
        _env_preamble(repo_root, modules, pythonpath=pythonpath),
        _log_header(
            scale=spec.scale, role="array task", run_code=run_code,
            experiment=experiment, pythonpath=pythonpath,
        ),
    )
    py = _uv_python_prefix(repo_root)

    sbatch_opts = {
        "array": f"0-{spec.num_array_tasks - 1}",
        # Slurm does not substitute ``%a`` in ``--job-name`` (only in
        # ``--output``/``--error``), so the whole array shares one name; the
        # per-task index is surfaced via the log header instead.
        "job-name": f"lysis-{spec.scale}-array__{run_code}",
        # Worker stdout/stderr (incl. the header + ``set -x`` trace) lands in
        # the staging dir alongside the generated scripts — not in ``.slurm``.
        "out": str(staging_dir / f"lysis-{spec.scale}-array__{run_code}__%a.out"),
        "nodes": 1,
        "mem": 3096,
        "ntasks": 1,
        "cpus-per-task": 1,
        "exclusive=user": "",
    }
    if partition:
        sbatch_opts["partition"] = partition
    sbatch_opts = apply_sbatch_overrides(sbatch_opts, sbatch_overrides)

    if spec.backend == "python":
        # ------------------------------------------------------------------
        # Python backend: single array task runs every simulation in series
        # via PythonRunner.execute(), writing directly to HDF5.  No binary,
        # no per-task subdirs, no setup files.
        # ------------------------------------------------------------------
        if fast_tmp_root is None:
            # Single-tier: write directly to the original .h5 path.
            init = _runner_init_line(
                spec, str(hdf5_path), None,
            )
            execute = f"""\
# Execute {spec.scale}scale Python runner via lysis
{py} -c "
from {spec.runner_module} import {spec.runner_class}
{init}.execute()
" """
            sections = [lead, execute]
        else:
            # Two-tier: cp the .h5 to node-local scratch, run there, cp back.
            # Use cp (not mv) so a mid-job crash leaves the original .h5 in
            # its MACRO_EMPTY state for clean resubmission.
            setup = f"""\
# Setup — create local fast dir, copy HDF5 file in
local_work_dir=$(mktemp -d -p "{fast_tmp_root}")
local_h5_path="${{local_work_dir}}/{hdf5_path.name}"
cp "{hdf5_path}" "${{local_h5_path}}" """

            init = _runner_init_line(
                spec, "${local_h5_path}", None,
            )
            execute = f"""\
# Execute {spec.scale}scale Python runner against the local HDF5 copy
{py} -c "
from {spec.runner_module} import {spec.runner_class}
{init}.execute()
" """

            move = f"""\
# Copy populated HDF5 back over the original, then clean up local dir
cp "${{local_h5_path}}" "{hdf5_path}"
rm -rf "${{local_work_dir}}" """

            sections = [lead, setup, execute, move]

    elif fast_tmp_root is None:
        # ------------------------------------------------------------------
        # Single-tier: Fortran writes directly to the shared staging dir.
        # The executable is pre-staged by submit_slurm_job() before any
        # array tasks start, so no cp is needed here.
        # ------------------------------------------------------------------
        init = _runner_init_line(
            spec, str(hdf5_path), f"{staging_dir}/{binary_name}",
            skip_binary_verification=skip_binary_verification,
            source_stamp=source_stamp,
        )
        execute = f"""\
# Execute {spec.scale}scale Fortran binary via lysis
{py} -c "
import os
from pathlib import Path
from {spec.runner_module} import {spec.runner_class}
fm = {init}
fm.exec_in_workdir(Path('{staging_dir}'))
" """

        sections = [lead, execute]

    else:
        # ------------------------------------------------------------------
        # Two-tier: copy setup files → local fast dir, run, then mv to staging
        # ------------------------------------------------------------------
        if spec.scale == "macro":
            # Macro tasks write a single named per-sim subdir — mv by name
            # (no glob, so it fails clearly on a missing operand).
            mv_section = 'mv "${local_datadir}/${SIM}" "${staging_datadir}/"'
        else:
            # Micro array tasks write flat __NN-suffixed files; guard the glob
            # mv so an empty local dir fails loud instead of choking on ``*``.
            mv_section = _guarded_glob_xfer(
                "mv", '"${local_datadir}"/*', '"${staging_datadir}/"',
                array_name="output_files",
                empty_msg="no per-task output in ${local_datadir} to move",
            )

        # Guarded copy of the staged setup files into the local fast dir (the
        # ``*`` would otherwise choke under ``set -e`` on an empty staging dir).
        guarded_setup_cp = _guarded_glob_xfer(
            "cp", '"${staging_datadir}"/*', '"${local_datadir}/"',
            array_name="staged_files",
            empty_msg="no staged setup files in ${staging_datadir}",
        )

        # Source the binary from the pre-staged staging-dir copy, not the
        # caller's --executable path: under --fortran-commit the build dir
        # is already gone by the time this script runs.
        setup = f"""\
# Setup — create local fast dir, copy setup files and binary
SIM=$(printf "%02d" ${{SLURM_ARRAY_TASK_ID}})
local_work_dir=$(mktemp -d -p "{fast_tmp_root}")
# Always reclaim node-local scratch on exit — including a ``set -e`` abort in
# the mv below — so a failed move never leaks the fast-tmp work dir.
trap 'rm -rf "${{local_work_dir}}"' EXIT
local_datadir="${{local_work_dir}}/data/{run_code}"
staging_datadir="{staging_dir}/data/{run_code}"
mkdir -p "${{local_datadir}}"
# Copy staged setup files into the local fast dir (glob guarded).
{guarded_setup_cp}
cp "{staging_dir}/{binary_name}" "${{local_work_dir}}/" """

        init = _runner_init_line(
            spec, str(hdf5_path), f"${{local_work_dir}}/{binary_name}",
            skip_binary_verification=skip_binary_verification,
            source_stamp=source_stamp,
        )
        execute = f"""\
# Execute {spec.scale}scale Fortran binary into local fast storage
{py} -c "
import os
from pathlib import Path
from {spec.runner_module} import {spec.runner_class}
fm = {init}
fm.exec_in_workdir(Path('${{local_work_dir}}'))
" """

        move = f"""\
# Move per-task output to shared staging (local dir reclaimed by EXIT trap)
{mv_section}"""

        sections = [lead, setup, execute, move]

    return gs.scripts.plain(sections, **sbatch_opts)


def submit_slurm_job(
    spec: _SlurmJobSpec,
    hdf5_path: Path,
    executable: Optional[Path],
    staging_root_dir: Path,
    *,
    partition: Optional[str] = None,
    fast_tmp_root: Optional[str] = None,
    keep_tmpdir: bool = False,
    modules: str = DEFAULT_MODULES,
    sbatch_overrides: Optional[Mapping[str, Optional[str]]] = None,
    historical_backend_attrs: Optional[dict] = None,
    experiment: Optional[str] = None,
) -> int:
    """Stage scripts and submit a master Slurm job for an array-mode dispatch.

    For ``spec.backend == "fortran"`` pre-stages setup files (and
    optionally the executable, per ``spec``), then writes the array bash
    script and master Python+bash scripts into a unique staging directory
    and ``sbatch``-submits the master.

    For ``spec.backend == "python"`` skips the setup-file and executable
    pre-staging — the Python runner writes directly to the HDF5 file and
    has no per-task input files — and otherwise follows the same flow:
    snapshot ``src/lysis/`` for PYTHONPATH pinning, write the array bash
    script and master Python+bash scripts, ``sbatch``-submit the master.

    :param spec: Scale-specific dispatch description.
    :type spec: _SlurmJobSpec
    :param hdf5_path: Resolved absolute path to the ``.h5`` file.
    :type hdf5_path: Path
    :param executable: Resolved absolute path to the Fortran binary.
        May be ``None`` when ``spec.backend == "python"``.
    :type executable: Path or None
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
    :param modules: Space-separated list of LMod module specs (include a
        Fortran compiler module so the binary finds its runtime), baked into both the master and array task scripts.
        Defaults to :data:`DEFAULT_MODULES`.
    :type modules: str, optional
    :param sbatch_overrides: User overrides for the ``#SBATCH`` header,
        applied to BOTH the master job and each array task (shape returned
        by :func:`parse_sbatch_tokens`).  Defaults to ``None``.
    :type sbatch_overrides: Mapping[str, str or None], optional
    :param historical_backend_attrs: Pre-computed binary-provenance dict
        from
        :func:`~lysis.tools.provenance.gather_historical_backend_provenance`
        for the historical-build workflow.  Baked into the master Python
        script so the master stamps it (with binary verification skipped)
        in place of the default ``<binary> --version`` query.  ``None``
        (default) preserves normal-mode provenance.
    :type historical_backend_attrs: dict, optional
    :param experiment: Best-effort experiment identifier shown in the log
        header of the master and array-task scripts; ``None`` (default)
        renders as ``(unknown)``.
    :type experiment: str, optional
    :return: Master Slurm job ID.
    :rtype: int
    :raises ImportError: If ``GooseSLURM`` is not installed.
    """
    _require_gs()
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
    # Pre-stage setup files (Fortran only — params.json, plus
    # macroscale_in for macro).  Python runners write directly to HDF5
    # and have no per-task setup files.
    # ------------------------------------------------------------------
    if spec.backend != "python":
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

    # Resolve src/fortran/ stamp once on the submit host (still inside the
    # repo) and bake it into the array script so compute-node tasks compare
    # the binary against this stamp instead of running git in the staging
    # dir (which is outside the checkout — the lookup would fail there and
    # raise StaleBinaryError).  Skipped under --fortran-commit, where the
    # binary↔source mismatch is intentional and the check is bypassed.
    source_stamp: Optional[tuple] = None
    if historical_backend_attrs is None:
        from lysis.tools.provenance import (  # noqa: PLC0415
            gather_fortran_source_provenance,
        )
        source_stamp = gather_fortran_source_provenance()

    # ------------------------------------------------------------------
    # Write array job script
    # ------------------------------------------------------------------
    array_script = generate_array_script(
        spec, staging_dir, run_code, hdf5_path, executable,
        partition=partition, fast_tmp_root=fast_tmp_root,
        slurm_log_dir=slurm_log_dir,
        modules=modules,
        sbatch_overrides=sbatch_overrides,
        historical_backend_attrs=historical_backend_attrs,
        pythonpath=pythonpath,
        source_stamp=source_stamp,
        experiment=experiment,
    )
    array_path = staging_dir / f"lysis-{spec.scale}-array__{run_code}.sh"
    array_path.write_text(array_script)
    array_path.chmod(0o755)

    # ------------------------------------------------------------------
    # Write master.py
    # ------------------------------------------------------------------
    binary_name = executable.name if executable is not None else None
    master_py_content = _generate_array_master_py(
        spec, staging_dir, hdf5_path, run_code, binary_name, keep_tmpdir,
        historical_backend_attrs=historical_backend_attrs,
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
            _env_preamble(repo_root, modules, pythonpath=pythonpath),
            _log_header(
                scale=spec.scale, role="master", run_code=run_code,
                experiment=experiment, pythonpath=pythonpath,
            ),
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
    modules: str = DEFAULT_MODULES,
    sbatch_overrides: Optional[Mapping[str, Optional[str]]] = None,
    historical_backend_attrs: Optional[dict] = None,
    pythonpath: "Path | str | None" = None,
    source_stamp: Optional[tuple] = None,
    experiment: Optional[str] = None,
) -> str:
    """Generate a Slurm bash script for a single child microscale Fortran job.

    Used by the legacy single-child path (``submit_micro_slurm_job``
    with ``num_children=None``).

    The child script sources the user's shell environment, copies the
    executable into the working directory, and calls
    :meth:`~lysis.execution.fortran_micro.FortranMicro.exec_in_workdir` via
    ``python -c``.  It begins with an identifying header (see
    :func:`_log_header`) and then enables ``set -euo pipefail`` + ``set -x``
    so the ``.out`` log records the job context and traces each shell step.

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
    :param slurm_log_dir: Retained for API compatibility; no longer used
        by this generator.  The child's ``--output`` now writes into
        *staging_dir* (alongside the generated scripts), so worker
        stdout/stderr — including the identifying header and ``set -x``
        trace — is co-located with the run's scripts and data rather than
        in ``hdf5_path.parent / ".slurm"``.
    :type slurm_log_dir: Path or str, optional
    :param modules: Space-separated list of LMod module specs (include a
        Fortran compiler module so the binary finds its runtime), baked into the generated script.  Defaults to
        :data:`DEFAULT_MODULES`.
    :type modules: str, optional
    :param sbatch_overrides: User overrides for the child's ``#SBATCH``
        header (shape returned by :func:`parse_sbatch_tokens`).  Defaults
        to ``None``.
    :type sbatch_overrides: Mapping[str, str or None], optional
    :param historical_backend_attrs: When non-``None``, indicates the
        binary was built from a historical commit (see
        :func:`~lysis.execution.historical_build.build_historical_binary`).
        The generated ``from_hdf5(...)`` call receives
        ``skip_binary_verification=True`` so ``_verify_binary_version``
        does not crash on the intentional binary↔source mismatch.  The
        dict itself is not baked into the child — the master script
        owns the actual provenance stamping.
    :type historical_backend_attrs: dict or None
    :param pythonpath: When set, baked into the generated script as
        ``export PYTHONPATH=...`` so the child imports ``lysis`` from
        that path rather than the live editable install.  ``None``
        (default) leaves ``PYTHONPATH`` untouched.
    :type pythonpath: Path or str, optional
    :param source_stamp: Pre-resolved ``(commit, dirty)`` for
        ``src/fortran/``.  Baked into the child's ``from_hdf5`` call so
        the binary↔source check on the compute node compares against the
        master's stamp instead of running ``git`` from the staging dir
        (which is outside the repo and would fail).  Ignored when
        ``historical_backend_attrs`` is set.  ``None`` (default)
        preserves the in-process git lookup.
    :type source_stamp: tuple[str, str] or None
    :param experiment: Best-effort experiment identifier shown in the log
        header; ``None`` (default) renders as ``(unknown)``.
    :type experiment: str, optional
    :return: Slurm bash script text.
    :rtype: str
    :raises ImportError: If ``GooseSLURM`` is not installed.
    """
    _require_gs()
    staging_dir = Path(staging_dir)
    hdf5_path = Path(hdf5_path)
    executable = Path(executable)
    binary_name = executable.name

    repo_root = _repo_root()
    # Leading section: env setup runs first, then the header prints and strict
    # mode + tracing engage for the rest of the child script.
    lead = "{}\n{}".format(
        _env_preamble(repo_root, modules, pythonpath=pythonpath),
        _log_header(
            scale="micro", role="micro child", run_code=run_code,
            experiment=experiment, pythonpath=pythonpath,
        ),
    )
    py = _uv_python_prefix(repo_root)

    sbatch_opts = {
        "job-name": f"lysis-micro-child__{run_code}",
        # Worker stdout/stderr (incl. the header + ``set -x`` trace) lands in
        # the staging dir alongside the generated scripts — not in ``.slurm``.
        "out": str(staging_dir / f"lysis-micro-child__{run_code}.out"),
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
        if historical_backend_attrs is not None
        else ""
    )
    if source_stamp is not None and historical_backend_attrs is None:
        commit, dirty = source_stamp
        source_stamp_kwarg = f", source_stamp=('{commit}', '{dirty}')"
    else:
        source_stamp_kwarg = ""

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
{py} -c "
from pathlib import Path
from lysis.execution.fortran_micro import FortranMicro
fm = FortranMicro.from_hdf5('{hdf5_path}', '{staging_dir}/{binary_name}', out_file_code='{out_code}'{skip_verify_kwarg}{source_stamp_kwarg})
fm.exec_in_workdir(Path('{staging_dir}'))
" """

        sections = [lead, setup, execute]

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
# Always reclaim node-local scratch on exit — including a ``set -e`` abort in
# the mv below — so a failed move never leaks the fast-tmp work dir.
trap 'rm -rf "${{local_work_dir}}"' EXIT
local_datadir="${{local_work_dir}}/data/{run_code}"
staging_datadir="{staging_dir}/data/{run_code}"
mkdir -p "${{local_datadir}}"
mkdir -p "${{staging_datadir}}"
cp "{staging_dir}/{binary_name}" "${{local_work_dir}}/" """

        execute = f"""\
# Execute Fortran binary into local fast storage
{py} -c "
from pathlib import Path
from lysis.execution.fortran_micro import FortranMicro
fm = FortranMicro.from_hdf5('{hdf5_path}', '${{local_work_dir}}/{binary_name}', out_file_code='{out_code}'{skip_verify_kwarg}{source_stamp_kwarg})
fm.exec_in_workdir(Path('${{local_work_dir}}'))
" """

        # Guard the glob mv so an empty local dir fails loud (clear lysis
        # message + non-zero exit) instead of choking on an unexpanded ``*``
        # under ``set -e``; the EXIT trap still reclaims the fast-tmp dir.
        guarded_move = _guarded_glob_xfer(
            "mv", '"${local_datadir}"/*', '"${staging_datadir}/"',
            array_name="output_files",
            empty_msg="no output in ${local_datadir} to move",
        )
        move = f"""\
# Move output to shared staging (local fast dir reclaimed by EXIT trap)
{guarded_move}"""

        sections = [lead, setup, execute, move]

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
    _require_gs()
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
    _require_gs()
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


def _micro_spec(num_children: int, out_code: str, direct: bool = False) -> _SlurmJobSpec:
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
        direct=direct,
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
    modules: str = DEFAULT_MODULES,
    sbatch_overrides: Optional[Mapping[str, Optional[str]]] = None,
    historical_backend_attrs: Optional[dict] = None,
    experiment: Optional[str] = None,
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
        modules=modules,
        sbatch_overrides=sbatch_overrides,
        historical_backend_attrs=historical_backend_attrs,
        experiment=experiment,
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
    direct: bool = False,
    modules: str = DEFAULT_MODULES,
    sbatch_overrides: Optional[Mapping[str, Optional[str]]] = None,
    historical_backend_attrs: Optional[dict] = None,
    experiment: Optional[str] = None,
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
    :param modules: Space-separated list of LMod module specs (include a
        Fortran compiler module so the binary finds its runtime), baked into every generated script (master, array tasks,
        and the legacy single child).  Defaults to
        :data:`DEFAULT_MODULES`.
    :type modules: str, optional
    :param sbatch_overrides: User overrides for the ``#SBATCH`` header,
        applied to BOTH the master and child/array task scripts (shape
        returned by :func:`parse_sbatch_tokens`).  Defaults to ``None``.
    :type sbatch_overrides: Mapping[str, str or None], optional
    :param historical_backend_attrs: Pre-computed binary-provenance dict
        from
        :func:`~lysis.tools.provenance.gather_historical_backend_provenance`
        for the historical-build workflow.  Threaded through to the
        master Python script (in both array and legacy single-child
        paths) so the master stamps it in place of the default
        ``<binary> --version`` query.  ``None`` (default) preserves
        normal-mode provenance.
    :type historical_backend_attrs: dict, optional
    :param experiment: Best-effort experiment identifier shown in the log
        header of every generated script; ``None`` (default) renders as
        ``(unknown)``.
    :type experiment: str, optional
    :return: Master Slurm job ID.
    :rtype: int
    :raises subprocess.CalledProcessError: If ``sbatch`` fails.
    :raises ValueError: If ``num_children`` is set and is either less than
        1 or greater than ``micro_simulations``.
    """
    _require_gs()
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
        if num_children < 2:
            raise ValueError(
                f"num_children must be >= 2 when set, got {num_children}"
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
        spec = _micro_spec(num_children, out_code, direct=direct)
        return submit_slurm_job(
            spec, hdf5_path, executable, staging_root_dir,
            partition=partition, fast_tmp_root=fast_tmp_root,
            keep_tmpdir=keep_tmpdir,
            modules=modules,
            sbatch_overrides=sbatch_overrides,
            historical_backend_attrs=historical_backend_attrs,
            experiment=experiment,
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

    # Resolve the src/fortran/ stamp once before the child runs from
    # outside the repo (see submit_slurm_job for the analogous comment in
    # the array path).
    source_stamp: Optional[tuple] = None
    if historical_backend_attrs is None:
        from lysis.tools.provenance import (  # noqa: PLC0415
            gather_fortran_source_provenance,
        )
        source_stamp = gather_fortran_source_provenance()

    # ------------------------------------------------------------------
    # Write child script
    # ------------------------------------------------------------------
    child_script = generate_micro_child_script(
        staging_dir, run_code, hdf5_path, executable, out_code,
        partition=partition, fast_tmp_root=fast_tmp_root,
        slurm_log_dir=slurm_log_dir,
        modules=modules,
        sbatch_overrides=sbatch_overrides,
        historical_backend_attrs=historical_backend_attrs,
        pythonpath=pythonpath,
        source_stamp=source_stamp,
        experiment=experiment,
    )
    child_path = staging_dir / f"lysis-micro-child__{run_code}.sh"
    child_path.write_text(child_script)
    child_path.chmod(0o755)

    # ------------------------------------------------------------------
    # Write master.py
    # ------------------------------------------------------------------
    master_py_content = _generate_micro_master_py(
        staging_dir, hdf5_path, run_code, out_code, executable.name, keep_tmpdir,
        historical_backend_attrs=historical_backend_attrs,
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
            _env_preamble(repo_root, modules, pythonpath=pythonpath),
            _log_header(
                scale="micro", role="master", run_code=run_code,
                experiment=experiment, pythonpath=pythonpath,
            ),
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
    modules: str = DEFAULT_MODULES,
    sbatch_overrides: Optional[Mapping[str, Optional[str]]] = None,
    historical_backend_attrs: Optional[dict] = None,
    experiment: Optional[str] = None,
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
        modules=modules,
        sbatch_overrides=sbatch_overrides,
        historical_backend_attrs=historical_backend_attrs,
        experiment=experiment,
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
    modules: str = DEFAULT_MODULES,
    sbatch_overrides: Optional[Mapping[str, Optional[str]]] = None,
    historical_backend_attrs: Optional[dict] = None,
    experiment: Optional[str] = None,
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
    :param modules: Space-separated list of LMod module specs (include a
        Fortran compiler module so the binary finds its runtime), baked into the generated master and array task scripts.
        Defaults to :data:`DEFAULT_MODULES`.
    :type modules: str, optional
    :param sbatch_overrides: User overrides for the ``#SBATCH`` header,
        applied to BOTH the master and array task scripts (shape returned
        by :func:`parse_sbatch_tokens`).  Defaults to ``None``.
    :type sbatch_overrides: Mapping[str, str or None], optional
    :param historical_backend_attrs: Pre-computed binary-provenance dict
        from
        :func:`~lysis.tools.provenance.gather_historical_backend_provenance`
        for the historical-build workflow.  Threaded through to the
        master Python script so the master stamps it in place of the
        default ``<binary> --version`` query.  ``None`` (default)
        preserves normal-mode provenance.
    :type historical_backend_attrs: dict, optional
    :param experiment: Best-effort experiment identifier shown in the log
        header of the master and array-task scripts; ``None`` (default)
        renders as ``(unknown)``.
    :type experiment: str, optional
    :return: Master Slurm job ID.
    :rtype: int
    :raises subprocess.CalledProcessError: If ``sbatch`` fails.
    :raises ValueError: If the HDF5 file does not contain ``macro_params``.
    """
    _require_gs()
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
        modules=modules,
        sbatch_overrides=sbatch_overrides,
        historical_backend_attrs=historical_backend_attrs,
        experiment=experiment,
    )


def _python_macro_spec() -> _SlurmJobSpec:
    """Build the macroscale Python-backend dispatch spec.

    Single-task array (``num_array_tasks=1``) — HDF5's single-writer
    constraint forces every simulation in one Run to run in series within
    one job.  Each Run is a separate ``.h5`` file, so the CLI batch loop
    will fan out one independent master job per Run when invoked against
    a directory.

    No executable to pre-stage, no per-task setup files, no per-task
    output to concatenate, no NFS flush wait — the task writes directly
    to HDF5 via :class:`~lysis.execution.python_macro.PythonMacro` and
    the master has nothing left to do once the task completes.
    """
    return _SlurmJobSpec(
        scale="macro",
        runner_module="lysis.execution.python_macro",
        runner_class="PythonMacro",
        num_array_tasks=1,
        needs_concat_step=False,
        nfs_wait_seconds=0,
        prestage_executable=False,
        out_file_code="",
        in_file_code=None,
        num_children=None,
        backend="python",
    )


def submit_python_macro_slurm_job(
    hdf5_path: "Path | str",
    *,
    staging_root: Optional["Path | str"] = None,
    partition: Optional[str] = None,
    fast_tmp_root: Optional[str] = None,
    keep_tmpdir: bool = False,
    modules: str = DEFAULT_MODULES,
    sbatch_overrides: Optional[Mapping[str, Optional[str]]] = None,
) -> int:
    """Submit the Python-backend macroscale workflow as a Slurm master job.

    The master job submits a single array task (``--array=0-0``) that
    runs :meth:`~lysis.execution.python_macro.PythonRunner.execute` —
    every simulation for this Run, in series, writing straight into the
    HDF5 file's open DataStore.  The master then exits without touching
    the file.  ``src/lysis/`` is snapshotted into the staging directory
    and bound to ``PYTHONPATH`` so the job imports the frozen package
    (immune to source edits while queued).

    Directory naming
    ~~~~~~~~~~~~~~~~
    Same as :func:`submit_macro_slurm_job` — staging dir is::

        tempfile.mkdtemp(
            prefix=f"lysis-macro-{run_code}-",
            dir=staging_root or hdf5_path.parent,
        )

    HDF5 staging
    ~~~~~~~~~~~~
    * **Without** ``fast_tmp_root``: the array task writes directly to
      ``hdf5_path`` (typically a shared filesystem).
    * **With** ``fast_tmp_root``: the array task ``cp``\\ s ``hdf5_path``
      into ``mktemp -d -p {fast_tmp_root}`` on the compute node, runs
      against that copy, then ``cp``\\ s the populated file back over
      ``hdf5_path`` on success.  A mid-job crash leaves the original
      ``.h5`` in its ``MACRO_EMPTY`` state for clean resubmission.

    :param hdf5_path: Full path to the run's ``.h5`` file.  Must contain
        both ``micro_params`` and ``macro_params``.
    :type hdf5_path: Path or str
    :param staging_root: Root directory under which the staging temp dir
        is created.  Defaults to ``hdf5_path.parent``.
    :type staging_root: Path or str, optional
    :param partition: Slurm partition for both master and array jobs.
    :type partition: str, optional
    :param fast_tmp_root: Root directory for node-local fast scratch
        (opt-in two-tier HDF5 staging).  When ``None`` (default) the
        array task writes directly to ``hdf5_path``.
    :type fast_tmp_root: str, optional
    :param keep_tmpdir: If ``True``, preserve the staging directory
        after the master job completes (useful for debugging).
    :type keep_tmpdir: bool, optional
    :param modules: Space-separated list of LMod module specs loaded by the
        generated master and array task scripts — exported by the
        Python-only path too so the preamble matches the Fortran path
        verbatim and node module state remains predictable.  Defaults to
        :data:`DEFAULT_MODULES`.
    :type modules: str, optional
    :param sbatch_overrides: User overrides for the ``#SBATCH`` header,
        applied to BOTH the master and array task scripts (shape returned
        by :func:`parse_sbatch_tokens`).  Defaults to ``None``.
    :type sbatch_overrides: Mapping[str, str or None], optional
    :return: Master Slurm job ID.
    :rtype: int
    :raises subprocess.CalledProcessError: If ``sbatch`` fails.
    """
    hdf5_path = Path(hdf5_path).resolve()

    staging_root_dir = (
        Path(staging_root) if staging_root is not None else hdf5_path.parent
    )

    spec = _python_macro_spec()
    return submit_slurm_job(
        spec, hdf5_path, None, staging_root_dir,
        partition=partition, fast_tmp_root=fast_tmp_root,
        keep_tmpdir=keep_tmpdir,
        modules=modules,
        sbatch_overrides=sbatch_overrides,
    )
