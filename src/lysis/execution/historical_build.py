"""Build Fortran binaries from an arbitrary historical commit.

Backs the ``--fortran-commit <ref>`` flag on ``lysis run-micro`` and
``lysis run-macro`` (see :mod:`lysis.cli.run_micro`,
:mod:`lysis.cli.run_macro`).  The workflow:

1. Resolve *ref* to a full SHA via ``git rev-parse <ref>^{commit}``.
2. ``git archive`` ``src/fortran/``, ``src/c/``, and ``Makefile`` from
   that SHA into a fresh temp directory (no worktree state, no
   collision with parallel CLI invocations).
3. Run ``make`` inside the temp directory, capturing stdout/stderr to
   ``.build_log``.  When a *compiler_module* is supplied, the build is
   wrapped in ``bash -c 'module purge && module load <module> && …'``
   so the binary links against the same toolchain that Slurm jobs will
   load at runtime.
4. Compute a synthesised provenance dict via
   :func:`~lysis.tools.provenance.gather_historical_binary_provenance`.
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
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Iterator, Optional

from ..tools.provenance import gather_historical_binary_provenance
from ..tools.provenance._git import _git, _package_repo_root

__all__ = [
    "HistoricalBuildError",
    "resolve_ref",
    "build_historical_binary",
]


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


def _make_argv(build_dir: Path, compiler_module: Optional[str]) -> list[str]:
    """Build the argv that runs ``make`` in *build_dir*, optionally module-wrapped.

    When *compiler_module* is given, the build is invoked inside a
    short bash script that loads the requested LMod module first.
    """
    if compiler_module is None:
        return ["make", "-C", str(build_dir)]
    script = (
        "[ -z \"${LMOD_CMD:-}\" ] && [ -f /etc/profile.d/lmod.sh ] "
        "&& source /etc/profile.d/lmod.sh; "
        f"module purge && module load {compiler_module} && "
        f'make -C "{build_dir}"'
    )
    return ["bash", "-c", script]


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
    compiler_module: Optional[str] = None,
    keep_dir: bool = False,
    repo_root: Optional[Path] = None,
) -> Iterator[tuple[Path, dict]]:
    """Build Fortran binaries at *ref* and yield ``(binary_path, provenance_dict)``.

    On context entry: resolves *ref*, extracts source via ``git archive``,
    runs ``make`` (optionally inside ``module load <compiler_module>``),
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
    :param compiler_module: LMod module spec (e.g.
        ``"intel-compilers/2023"``) to load around ``make`` and the
        compiler-version probe.  ``None`` (default) uses the current
        environment; appropriate for local-mode runs where the same
        environment will execute the binary.
    :type compiler_module: str or None
    :param keep_dir: Preserve the build directory after the context
        exits.  Useful for debugging.  Defaults to ``False``.
    :type keep_dir: bool, optional
    :param repo_root: Repository to extract from.  Defaults to the
        lysis package's own repo.
    :type repo_root: pathlib.Path, optional
    :yields: ``(binary_path, provenance_dict)`` where *binary_path* is
        the absolute path to the requested binary inside the build dir,
        and *provenance_dict* is the result of
        :func:`~lysis.tools.provenance.gather_historical_binary_provenance`.
    :raises HistoricalBuildError: If ref resolution, source extraction,
        ``make``, or binary lookup fails.
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
        build_log = build_dir / ".build_log"
        with build_log.open("w") as fh:
            result = subprocess.run(
                _make_argv(build_dir, compiler_module),
                stdout=fh,
                stderr=subprocess.STDOUT,
                check=False,
            )
        if result.returncode != 0:
            raise HistoricalBuildError(
                _format_build_failure(build_log, result.returncode)
            )

        binary_path = build_dir / "bin" / binary_basename
        if not binary_path.is_file():
            available = sorted(
                p.name for p in (build_dir / "bin").iterdir()
            ) if (build_dir / "bin").is_dir() else []
            raise HistoricalBuildError(
                f"Build produced no binary named {binary_basename!r} in "
                f"{build_dir / 'bin'}; available: {available}."
            )

        provenance = gather_historical_binary_provenance(
            binary_path,
            resolved_sha=resolved_sha,
            build_log=build_log,
            compiler_module=compiler_module,
        )
        yield binary_path, provenance
    except BaseException:
        if not keep_dir:
            shutil.rmtree(build_dir, ignore_errors=True)
        raise
    else:
        if not keep_dir:
            shutil.rmtree(build_dir, ignore_errors=True)
