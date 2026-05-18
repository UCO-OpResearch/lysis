"""Execution provenance stamps written to HDF5 scale-output groups.

Captures *when*, *where*, and *with what code* a microscale or macroscale
simulation was executed so that the resulting HDF5 file is self-describing
for future reproducibility checks.

The single public entry point is :func:`gather_execution_provenance`, which
returns a dict keyed by the attribute names defined in
:mod:`lysis.config.constants`. :class:`~lysis.dataio.datastore.DataStore`
writes these attributes on the ``microscale_out`` / ``macroscale_out``
group at the end of ``import_collection``.

Version resolution order:

1. Exact-match git tag at ``HEAD`` (``git describe --exact-match --tags``);
   branch pointers do not count.
2. Full commit hash at ``HEAD`` (``git rev-parse HEAD``).
3. Installed package version from ``importlib.metadata.version("lysis")``.
4. ``"unknown"`` if all of the above fail.

Dirty trees emit a :class:`UserWarning` but do not block execution.
"""

import socket
import warnings
from datetime import datetime
from importlib.metadata import PackageNotFoundError, version as _pkg_version

from ...config.constants import CONST
from ._git import _git, _package_repo_root

__all__ = ["gather_execution_provenance"]


def _resolve_version() -> str:
    repo_root = _package_repo_root()
    if repo_root is not None:
        tag = _git(["describe", "--exact-match", "--tags", "HEAD"], repo_root)
        if tag:
            return tag
        commit = _git(["rev-parse", "HEAD"], repo_root)
        if commit:
            return commit
    try:
        return _pkg_version("lysis")
    except PackageNotFoundError:
        return "unknown"


def _resolve_dirty() -> bool:
    repo_root = _package_repo_root()
    if repo_root is None:
        return False
    status = _git(["status", "--porcelain"], repo_root)
    return bool(status)


def _resolve_timestamp() -> str:
    # %Z yields the local tz abbreviation on Linux systems with tzdata
    # (CDT/CST/EST/…); falls back to a numeric offset on minimal containers.
    return datetime.now().astimezone().strftime("%Y-%m-%dT%H:%M:%S %Z")


def _resolve_hostname() -> str:
    return socket.getfqdn()


def gather_execution_provenance() -> dict:
    """Collect provenance stamps for the current execution.

    :return: Dict keyed by the ``EXECUTION_*_ATTR`` names in
        :data:`lysis.config.constants.CONST`, ready to write as HDF5 attrs.
    :rtype: dict

    Emits a :class:`UserWarning` if the lysis working tree is dirty, so the
    user knows the recorded version does not uniquely identify the code that
    actually ran.
    """
    dirty = _resolve_dirty()
    if dirty:
        warnings.warn(
            "lysis working tree is dirty at execution time; "
            "recorded execution_version does not uniquely identify the code that ran.",
            UserWarning,
            stacklevel=2,
        )
    return {
        CONST.EXECUTION_VERSION_ATTR: _resolve_version(),
        CONST.EXECUTION_DIRTY_ATTR: dirty,
        CONST.EXECUTION_TIMESTAMP_ATTR: _resolve_timestamp(),
        CONST.EXECUTION_HOSTNAME_ATTR: _resolve_hostname(),
    }
