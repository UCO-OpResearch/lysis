"""Pipeline-time and init-time provenance stamps for HDF5 outputs.

Captures *when*, *where*, and *with what code* the ``src/lysis/`` Python
pipeline initialised or orchestrated a microscale or macroscale simulation,
so the resulting HDF5 file is self-describing for future reproducibility
checks.  (The simulation *engine* that produced the output is recorded
separately as the backend; see
:mod:`lysis.tools.provenance.binary`.)

Two public gather functions are exposed:

- :func:`gather_init_provenance` — written by ``init-experiment`` /
  ``init-macroscale``, keyed by ``CONST.INIT_*_ATTR``.
- :func:`gather_pipeline_provenance` — written by ``run-micro`` /
  ``run-macro`` at HDF5 import time, keyed by ``CONST.PIPELINE_*_ATTR``.

Both functions scope their version/dirty check to the ``src/lysis/``
subtree only (mirroring how :func:`~lysis.tools.provenance.binary.
gather_fortran_source_provenance` scopes to ``src/fortran/``).  Changes
to ``src/fortran/``, ``docs/``, or ``tests/`` do not flag a clean
``src/lysis/`` as dirty.

Version resolution order (path-scoped to ``src/lysis/``):

1. Exact-match git tag at ``HEAD`` (``git describe --exact-match --tags``),
   but only when HEAD is the most recent commit touching ``src/lysis/``.
2. Most recent commit hash that touched ``src/lysis/``
   (``git log -1 --format=%H HEAD -- src/lysis``).
3. Installed package version from ``importlib.metadata.version("lysis")``.
4. ``"unknown"`` if all of the above fail.

Dirty bit: ``git diff --quiet HEAD -- src/lysis`` (exit 0 → ``"clean"``;
non-zero → ``"dirty"``; subprocess failure → ``"unknown"``).

A dirty working tree emits a :class:`UserWarning` once per process the
first time either gather function is called.  The lysis CLI top-level
callback prints a coloured warning of its own and then calls
:func:`mark_dirty_warning_emitted` to suppress the redundant inline
warning; non-CLI callers (notebooks, ad-hoc scripts) still see the
warning naturally.
"""

import socket
import subprocess
import warnings
from datetime import datetime
from importlib.metadata import PackageNotFoundError, version as _pkg_version

from ...config.constants import CONST
from ._git import _git, _package_repo_root

__all__ = [
    "gather_pipeline_provenance",
    "gather_init_provenance",
    "allow_dirty_from_env",
    "allow_commit_mismatch_from_env",
    "mark_dirty_warning_emitted",
]


# Module-level dedup flag.  Process-local; reset across pytest runs by
# tests that exercise the warning behaviour.
_dirty_warning_emitted: bool = False


def _resolve_version() -> str:
    """Resolve the version stamp for the current ``src/lysis/`` source tree.

    Resolution order (path-scoped to ``src/lysis/``):

    1. Exact-match tag at HEAD, *only* when HEAD == the most recent
       commit that touched ``src/lysis/``.
    2. Most recent commit hash touching ``src/lysis/``.
    3. Installed package version.
    4. ``"unknown"``.
    """
    repo_root = _package_repo_root()
    if repo_root is not None:
        lysis_path = repo_root / "src" / "lysis"
        path_commit = _git(
            ["log", "-1", "--format=%H", "HEAD", "--", str(lysis_path)],
            repo_root,
        )
        if path_commit:
            head_commit = _git(["rev-parse", "HEAD"], repo_root)
            if head_commit == path_commit:
                tag = _git(
                    ["describe", "--exact-match", "--tags", "HEAD"], repo_root
                )
                if tag:
                    return tag
            return path_commit
    try:
        return _pkg_version("lysis")
    except PackageNotFoundError:
        return "unknown"


def _resolve_dirty() -> str:
    """Return ``"clean"``, ``"dirty"``, or ``"unknown"`` for ``src/lysis/``.

    Uses ``git diff --quiet HEAD -- src/lysis`` — mirrors the pattern in
    :func:`~lysis.tools.provenance.binary.gather_fortran_source_provenance`.
    """
    repo_root = _package_repo_root()
    if repo_root is None:
        return "unknown"
    lysis_path = repo_root / "src" / "lysis"
    try:
        result = subprocess.run(
            [
                "git", "-C", str(repo_root),
                "diff", "--quiet", "HEAD", "--", str(lysis_path),
            ],
            check=False,
            capture_output=True,
        )
    except (FileNotFoundError, OSError):
        return "unknown"
    return "dirty" if result.returncode != 0 else "clean"


def _resolve_timestamp() -> str:
    # %Z yields the local tz abbreviation on Linux systems with tzdata
    # (CDT/CST/EST/…); falls back to a numeric offset on minimal containers.
    return datetime.now().astimezone().strftime("%Y-%m-%dT%H:%M:%S %Z")


def _resolve_hostname() -> str:
    return socket.getfqdn()


def _emit_dirty_warning_if_needed(dirty: str) -> None:
    """Emit a :class:`UserWarning` the first time ``src/lysis/`` is dirty.

    Subsequent calls in the same process are silent (mirrors Python's
    default warning-filter dedup behaviour).
    """
    global _dirty_warning_emitted
    if dirty == "dirty" and not _dirty_warning_emitted:
        warnings.warn(
            "src/lysis/ has uncommitted changes; recorded version does not "
            "uniquely identify the code that ran.",
            UserWarning,
            stacklevel=3,
        )
        _dirty_warning_emitted = True


def mark_dirty_warning_emitted() -> None:
    """Mark the dirty warning as already emitted for the current process.

    Called by the lysis CLI top-level callback after it prints its own
    coloured warning, so the inline :func:`warnings.warn` inside
    :func:`gather_pipeline_provenance` / :func:`gather_init_provenance`
    does not re-fire for the same invocation.  Non-CLI callers
    (notebooks, scripts) never call this and therefore receive the
    inline warning normally.
    """
    global _dirty_warning_emitted
    _dirty_warning_emitted = True


def gather_pipeline_provenance() -> dict:
    """Collect pipeline-time provenance stamps for ``src/lysis/``.

    :return: Dict keyed by the ``PIPELINE_*_ATTR`` names in
        :data:`lysis.config.constants.CONST`, ready to write as HDF5 attrs.
        ``PIPELINE_DIRTY_ATTR`` is the 3-state string
        ``"clean"|"dirty"|"unknown"`` (matching ``INIT_DIRTY_ATTR`` and
        ``BACKEND_DIRTY_ATTR``), so the ``"unknown"`` case (git unavailable
        / repo root not found) is never collapsed into ``"clean"``.
    :rtype: dict

    Emits a :class:`UserWarning` the first time a dirty ``src/lysis/``
    is observed in the current process (see
    :func:`mark_dirty_warning_emitted`).
    """
    dirty_str = _resolve_dirty()
    _emit_dirty_warning_if_needed(dirty_str)
    return {
        CONST.PIPELINE_VERSION_ATTR: _resolve_version(),
        CONST.PIPELINE_DIRTY_ATTR: dirty_str,
        CONST.PIPELINE_TIMESTAMP_ATTR: _resolve_timestamp(),
        CONST.PIPELINE_HOSTNAME_ATTR: _resolve_hostname(),
    }


def gather_init_provenance() -> dict:
    """Collect init-time provenance stamps for ``src/lysis/``.

    Same shape as :func:`gather_pipeline_provenance` but keyed by
    ``CONST.INIT_*_ATTR`` and stores the dirty value as the 3-state
    string ``"clean"|"dirty"|"unknown"`` (matches
    :func:`~lysis.tools.provenance.binary.gather_fortran_source_provenance`).

    :return: Dict keyed by ``CONST.INIT_VERSION_ATTR``,
        ``CONST.INIT_DIRTY_ATTR``, ``CONST.INIT_TIMESTAMP_ATTR``,
        ``CONST.INIT_HOSTNAME_ATTR``.
    :rtype: dict

    Emits a :class:`UserWarning` once per process if ``src/lysis/`` is
    dirty, just like :func:`gather_pipeline_provenance`.
    """
    dirty_str = _resolve_dirty()
    _emit_dirty_warning_if_needed(dirty_str)
    return {
        CONST.INIT_VERSION_ATTR: _resolve_version(),
        CONST.INIT_DIRTY_ATTR: dirty_str,
        CONST.INIT_TIMESTAMP_ATTR: _resolve_timestamp(),
        CONST.INIT_HOSTNAME_ATTR: _resolve_hostname(),
    }


def _env_truthy(env_var: str) -> bool:
    import os

    raw = os.environ.get(env_var, "").strip().lower()
    return raw in {"1", "true", "yes", "on"}


def allow_dirty_from_env() -> bool:
    """Return ``True`` if ``LYSIS_ALLOW_DIRTY`` is set to a truthy value.

    Accepted truthy strings (case-insensitive): ``1``, ``true``, ``yes``,
    ``on``.  Anything else (including unset) is ``False``.
    """
    return _env_truthy(CONST.LYSIS_ALLOW_DIRTY_ENV)


def allow_commit_mismatch_from_env() -> bool:
    """Return ``True`` if ``LYSIS_ALLOW_COMMIT_MISMATCH`` is set truthy.

    Same truthy-string rules as :func:`allow_dirty_from_env`.
    """
    return _env_truthy(CONST.LYSIS_ALLOW_COMMIT_MISMATCH_ENV)
