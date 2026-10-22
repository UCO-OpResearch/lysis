"""Build Fortran binaries from an arbitrary historical commit.

Backs the ``--fortran-commit <ref>`` flag on ``lysis run-micro`` and
``lysis run-macro`` (see :mod:`lysis.cli.run_micro`,
:mod:`lysis.cli.run_macro`).  The workflow:

1. Resolve *ref* to a full SHA via ``git rev-parse <ref>^{commit}``.
2. ``git archive`` ``src/fortran/``, ``src/c/``, and ``Makefile`` from
   that SHA into a fresh temp directory (no worktree state, no
   collision with parallel CLI invocations).
3. Run ``make`` inside the temp directory, capturing stdout/stderr to
   ``.build_log``.  When *modules* is supplied, the build is
   wrapped in ``bash -c 'module purge && module load <modules> && …'``
   so the binary links against the same toolchain that Slurm jobs will
   load at runtime.  When *compiler* is supplied, it is passed as
   ``make LYSIS_COMPILER=<compiler>``; otherwise the commit's own
   Makefile default applies (ifx since #137).  Commits from before #137
   have a Makefile that hard-codes ``ifort`` and ignores
   ``LYSIS_COMPILER``: they build only with ifort, which needs a module
   that still ships it (``intel-compilers/2024.2.0`` or older).
4. Compute a synthesised provenance dict via
   :func:`~lysis.tools.provenance.gather_historical_backend_provenance`.
5. Yield ``(binary_path, provenance_dict)`` for the caller to thread
   through to :class:`~lysis.execution.fortran_macro.FortranMacro` /
   :class:`~lysis.execution.fortran_micro.FortranMicro`.
6. On context exit, remove the build dir (unless ``keep_dir=True``).

The build dir lives in ``$TMPDIR`` (falling back to ``/tmp``) — never
``run.os_path`` — because a single CLI invocation produces one binary
shared by every run it processes, and ``run.os_path`` is per-run.  HPC
users whose ``/tmp`` is small/ramdisk are already expected to point
``$TMPDIR`` at scratch.
"""

import contextlib
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Iterator, Optional

from ..tools.provenance import gather_historical_backend_provenance
from ..tools.provenance._git import _git, _package_repo_root

__all__ = [
    "DEFAULT_FORTRAN_COMPILER",
    "FORTRAN_COMPILERS",
    "HistoricalBuildError",
    "resolve_ref",
    "build_historical_binary",
]

#: Fortran compilers the Makefile's ``LYSIS_COMPILER`` accepts.
FORTRAN_COMPILERS: tuple[str, ...] = ("ifx", "ifort", "gfortran")

#: The Makefile's default ``LYSIS_COMPILER``.  Informational: when no
#: compiler is requested, nothing is passed to ``make`` and the built
#: commit's own Makefile default applies.  A test keeps this in step with
#: ``LYSIS_COMPILER_DEFAULT`` in the Makefile.
DEFAULT_FORTRAN_COMPILER: str = "ifx"

#: Module suggested when a pre-#137 commit needs ifort.  The newest
#: ``intel-compilers`` module on the cluster that still ships ifort.
_IFORT_MODULE_HINT = "intel-compilers/2024.2.0"


class HistoricalBuildError(RuntimeError):
    """Raised when the historical-build pipeline cannot produce a binary."""


def resolve_ref(ref: str, repo_root: Optional[Path] = None) -> str:
    """Resolve a git ref to its full commit SHA.

    :param ref: Any git ref understood by ``git rev-parse`` — a full or
        abbreviated SHA, a branch, or a tag.  ``"HEAD"`` and similar
        symbolic refs are accepted.
    :type ref: str
    :param repo_root: Repository root containing ``.git/``.  Defaults to
        the lysis package's own repo (resolved from ``lysis.__file__``).
    :type repo_root: pathlib.Path, optional
    :return: Full 40-character SHA of the commit *ref* points to.
    :rtype: str
    :raises HistoricalBuildError: If *ref* does not resolve to a commit,
        or if the lysis package repo cannot be located.
    """
    if repo_root is None:
        repo_root = _package_repo_root()
    if repo_root is None:
        raise HistoricalBuildError(
            "Cannot locate the lysis package's git repository root; "
            "--fortran-commit needs the source tree to extract from."
        )
    sha = _git(["rev-parse", "--verify", f"{ref}^{{commit}}"], repo_root)
    if not sha:
        raise HistoricalBuildError(
            f"Git ref {ref!r} does not resolve to a commit in {repo_root}."
        )
    return sha


def _resolve_binary_name(executable_name: str) -> str:
    """Strip a leading ``bin/`` and reject names with path separators.

    When ``--fortran-commit`` is in play, ``--executable`` is interpreted
    as a basename within the build's ``bin/`` directory.  Callers may
    spell it with or without a leading ``bin/``; anything more elaborate
    is a misconfiguration.
    """
    name = executable_name.lstrip("/")
    if name.startswith("bin/"):
        name = name[len("bin/") :]
    if "/" in name or name in ("", ".", ".."):
        raise HistoricalBuildError(
            f"--executable must be a binary name (with or without a leading "
            f"'bin/'), not a path; got {executable_name!r}."
        )
    return name


def _archive_into(sha: str, repo_root: Path, build_dir: Path) -> None:
    """Extract ``src/fortran``, ``src/c``, and ``Makefile`` at *sha* into *build_dir*.

    Uses ``git archive | tar -x`` so no extra worktree state is created.
    """
    archive = subprocess.Popen(
        [
            "git", "-C", str(repo_root),
            "archive", sha, "--",
            "src/fortran", "src/c", "Makefile",
        ],
        stdout=subprocess.PIPE,
    )
    try:
        result = subprocess.run(
            ["tar", "-x", "-C", str(build_dir)],
            stdin=archive.stdout,
            capture_output=True,
            text=True,
            check=False,
        )
    finally:
        if archive.stdout is not None:
            archive.stdout.close()
        archive.wait()
    if archive.returncode != 0 or result.returncode != 0:
        raise HistoricalBuildError(
            f"git archive | tar -x failed for {sha}: "
            f"git rc={archive.returncode}, tar rc={result.returncode}; "
            f"tar stderr: {result.stderr.strip()}"
        )


def _makefile_supports_compiler_choice(build_dir: Path) -> bool:
    """Return ``True`` if the extracted Makefile understands ``LYSIS_COMPILER``.

    Makefiles from before #137 hard-code ``ifort`` and silently ignore the
    variable.  A missing or unreadable Makefile counts as not supporting it.
    """
    try:
        text = (build_dir / "Makefile").read_text(errors="replace")
    except OSError:
        return False
    return "LYSIS_COMPILER" in text


def _make_argv(
    build_dir: Path,
    modules: Optional[str],
    compiler: Optional[str] = None,
) -> list[str]:
    """Build the argv that runs ``make`` in *build_dir*, optionally module-wrapped.

    When *modules* is given, the build is invoked inside a
    short bash script that ``module load``s the requested
    (space-separated) LMod modules first.  When *compiler* is given, it is
    passed on the ``make`` command line as ``LYSIS_COMPILER=<compiler>``,
    which the Makefile ranks above any other source of the choice.
    """
    make_args = ["make", "-C", str(build_dir)]
    if compiler is not None:
        make_args.append(f"LYSIS_COMPILER={compiler}")
    if modules is None:
        return make_args
    script = (
        "[ -z \"${LMOD_CMD:-}\" ] && [ -f /etc/profile.d/lmod.sh ] "
        "&& source /etc/profile.d/lmod.sh; "
        f"module purge && module load {modules} && "
        + " ".join(f'"{a}"' for a in make_args)
    )
    return ["bash", "-c", script]


def _make_env() -> dict:
    """Return the environment for ``make``, without ``LYSIS_COMPILER``.

    The Makefile ignores ``LYSIS_COMPILER`` from the environment but warns
    about it, and that warning would name a compiler in the build log ahead
    of the one that ran, misleading
    :func:`~lysis.tools.provenance.binary._identify_compiler_binary`.
    """
    env = dict(os.environ)
    env.pop("LYSIS_COMPILER", None)
    return env


def _check_compiler_supported(
    sha: str, compiler: Optional[str], supports_choice: bool
) -> None:
    """Reject a compiler choice the commit's Makefile cannot honour.

    :raises HistoricalBuildError: If *compiler* is not a known value, or if
        a compiler other than ifort is requested for a pre-#137 commit.
    """
    if compiler is None:
        return
    if compiler not in FORTRAN_COMPILERS:
        raise HistoricalBuildError(
            f"Unknown Fortran compiler {compiler!r}; valid values: "
            f"{', '.join(FORTRAN_COMPILERS)}."
        )
    if not supports_choice and compiler != "ifort":
        raise HistoricalBuildError(
            f"Commit {sha[:7]} predates Fortran compiler selection: its "
            f"Makefile hard-codes ifort and cannot build with {compiler}.  "
            f"Omit the compiler choice (or choose ifort) and load a module "
            f"that ships ifort, e.g. {_IFORT_MODULE_HINT}."
        )


def _legacy_ifort_hint(build_log: Path, modules: Optional[str]) -> str:
    """Explain a failed pre-#137 build whose log shows ifort was missing."""
    try:
        text = build_log.read_text(errors="replace")
    except OSError:
        return ""
    missing = (
        "ifort: command not found" in text
        or "ifort: No such file or directory" in text
        or ("ifort" in text and "not found" in text)
    )
    if not missing:
        return ""
    where = f"in modules {modules!r}" if modules else "on PATH"
    return (
        f"\nThis commit predates Fortran compiler selection (#137): its "
        f"Makefile hard-codes ifort, which was not found {where}.  Load a "
        f"module that still ships ifort, e.g. {_IFORT_MODULE_HINT}; "
        f"intel-compilers 2025.x and 2026.x have only ifx."
    )


def _format_build_failure(build_log: Path, returncode: int) -> str:
    """Format a ClickException-friendly message including the build log tail."""
    try:
        tail = "\n".join(build_log.read_text(errors="replace").splitlines()[-40:])
    except OSError:
        tail = "(build log unreadable)"
    return (
        f"Historical Fortran build failed (make exit {returncode}).\n"
        f"Build log: {build_log}\n"
        f"--- last 40 lines ---\n{tail}"
    )


@contextlib.contextmanager
def build_historical_binary(
    ref: str,
    executable_name: str,
    *,
    modules: Optional[str] = None,
    compiler: Optional[str] = None,
    keep_dir: bool = False,
    repo_root: Optional[Path] = None,
) -> Iterator[tuple[Path, dict]]:
    """Build Fortran binaries at *ref* and yield ``(binary_path, provenance_dict)``.

    On context entry: resolves *ref*, extracts source via ``git archive``,
    runs ``make`` (optionally inside ``module load <modules>``, and with
    ``LYSIS_COMPILER=<compiler>`` when *compiler* is given),
    and computes a synthesised provenance dict.

    On context exit: removes the build directory unless ``keep_dir=True``.

    :param ref: Git ref to build from — any form accepted by
        :func:`resolve_ref`.
    :type ref: str
    :param executable_name: Basename (or ``bin/<basename>``) of the
        binary the caller wants to run.  Verified to exist after the
        build; the full path returned in the yielded tuple is
        ``<build_dir>/bin/<basename>``.
    :type executable_name: str
    :param modules: Space-separated list of LMod module specs (e.g.
        ``"intel-compilers/2024.2.0"`` or
        ``"intel-compilers/2024.2.0 GCC/12.3.0"``) to load around
        ``make`` and the compiler-version probe.  Include a module that
        provides the chosen Fortran compiler so the build succeeds.
        ``None`` (default) uses the current environment; appropriate for
        local-mode runs where the same environment will execute the binary.
    :type modules: str or None
    :param compiler: Fortran compiler to build with: one of
        :data:`FORTRAN_COMPILERS`.  ``None`` (default) passes nothing to
        ``make``, so the commit's own Makefile default applies: ifx
        (:data:`DEFAULT_FORTRAN_COMPILER`) for commits with compiler
        selection, and ifort for commits from before it.  Requesting
        anything but ifort for such an older commit raises
        :class:`HistoricalBuildError`.
    :type compiler: str or None
    :param keep_dir: Preserve the build directory after the context
        exits.  Useful for debugging.  Defaults to ``False``.
    :type keep_dir: bool, optional
    :param repo_root: Repository to extract from.  Defaults to the
        lysis package's own repo.
    :type repo_root: pathlib.Path, optional
    :yields: ``(binary_path, provenance_dict)`` where *binary_path* is
        the absolute path to the requested binary inside the build dir,
        and *provenance_dict* is the result of
        :func:`~lysis.tools.provenance.gather_historical_backend_provenance`.
    :raises HistoricalBuildError: If ref resolution, source extraction,
        ``make``, or binary lookup fails, or if *compiler* is unknown or
        cannot be used for the commit.
    """
    if repo_root is None:
        repo_root = _package_repo_root()
    if repo_root is None:
        raise HistoricalBuildError(
            "Cannot locate the lysis package's git repository root; "
            "--fortran-commit needs the source tree to extract from."
        )

    binary_basename = _resolve_binary_name(executable_name)
    resolved_sha = resolve_ref(ref, repo_root=repo_root)
    short_sha = resolved_sha[:7]
    build_dir = Path(tempfile.mkdtemp(prefix=f"lysis-fortran-{short_sha}-"))
    try:
        _archive_into(resolved_sha, repo_root, build_dir)
        # The Makefile's .build-pre target is empty (live builds rely on
        # ``bin/`` and ``lib/`` already existing in the developer's
        # checkout); ``git archive`` does not carry empty dirs, so we
        # create them here to keep the historical build self-contained.
        (build_dir / "bin").mkdir(exist_ok=True)
        (build_dir / "lib").mkdir(exist_ok=True)
        supports_choice = _makefile_supports_compiler_choice(build_dir)
        _check_compiler_supported(resolved_sha, compiler, supports_choice)
        build_log = build_dir / ".build_log"
        with build_log.open("w") as fh:
            result = subprocess.run(
                _make_argv(
                    build_dir,
                    modules,
                    compiler if supports_choice else None,
                ),
                stdout=fh,
                stderr=subprocess.STDOUT,
                check=False,
                env=_make_env(),
            )
        if result.returncode != 0:
            message = _format_build_failure(build_log, result.returncode)
            if not supports_choice:
                message += _legacy_ifort_hint(build_log, modules)
            raise HistoricalBuildError(message)

        binary_path = build_dir / "bin" / binary_basename
        if not binary_path.is_file():
            available = sorted(
                p.name for p in (build_dir / "bin").iterdir()
            ) if (build_dir / "bin").is_dir() else []
            raise HistoricalBuildError(
                f"Build produced no binary named {binary_basename!r} in "
                f"{build_dir / 'bin'}; available: {available}."
            )

        provenance = gather_historical_backend_provenance(
            binary_path,
            resolved_sha=resolved_sha,
            build_log=build_log,
            modules=modules,
        )
        yield binary_path, provenance
    except BaseException:
        if not keep_dir:
            shutil.rmtree(build_dir, ignore_errors=True)
        raise
    else:
        if not keep_dir:
            shutil.rmtree(build_dir, ignore_errors=True)
