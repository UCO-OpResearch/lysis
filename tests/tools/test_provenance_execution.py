"""Unit tests for :mod:`lysis.tools.provenance.execution`.

Covers :func:`gather_pipeline_provenance`, :func:`gather_init_provenance`,
the once-per-process dirty-warning machinery, and the path-scoped
``src/lysis/`` resolvers:

* ``_resolve_version``   — tag (only when HEAD == path commit) > path-
  scoped commit hash > pip metadata > ``"unknown"``.
* ``_resolve_dirty``     — 3-state string from ``git diff --quiet
  HEAD -- src/lysis``.
* ``_resolve_timestamp`` — format and tz suffix.
* ``_resolve_hostname``  — FQDN.
"""

import re
import warnings

import pytest

from lysis.config.constants import CONST
from lysis.tools.provenance import (
    execution,
    gather_pipeline_provenance,
    gather_init_provenance,
    mark_dirty_warning_emitted,
)
from lysis.tools.provenance.execution import (
    _resolve_dirty,
    _resolve_hostname,
    _resolve_timestamp,
    _resolve_version,
)


@pytest.fixture(autouse=True)
def reset_dirty_warning_flag():
    """Reset the module-level dedup flag between tests."""
    execution._dirty_warning_emitted = False
    yield
    execution._dirty_warning_emitted = False


def test_gather_returns_expected_keys():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        result = gather_pipeline_provenance()
    assert set(result.keys()) == {
        CONST.PIPELINE_VERSION_ATTR,
        CONST.PIPELINE_DIRTY_ATTR,
        CONST.PIPELINE_TIMESTAMP_ATTR,
        CONST.PIPELINE_HOSTNAME_ATTR,
    }


def test_gather_init_returns_expected_keys():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        result = gather_init_provenance()
    assert set(result.keys()) == {
        CONST.INIT_VERSION_ATTR,
        CONST.INIT_DIRTY_ATTR,
        CONST.INIT_TIMESTAMP_ATTR,
        CONST.INIT_HOSTNAME_ATTR,
    }


def test_version_is_nonempty_string():
    assert isinstance(_resolve_version(), str)
    assert _resolve_version()


def test_version_matches_path_scoped_commit_when_no_tag(monkeypatch, tmp_path):
    repo_root = tmp_path
    lysis_path = repo_root / "src" / "lysis"
    path_commit = "a" * 40

    monkeypatch.setattr(execution, "_package_repo_root", lambda: repo_root)

    def fake_git(args, root):
        assert root == repo_root
        if args == ["log", "-1", "--format=%H", "HEAD", "--", str(lysis_path)]:
            return path_commit
        if args == ["rev-parse", "HEAD"]:
            return path_commit
        if args == ["describe", "--exact-match", "--tags", "HEAD"]:
            return None
        raise AssertionError(f"unexpected git call: {args}")

    monkeypatch.setattr(execution, "_git", fake_git)
    assert _resolve_version() == path_commit


def test_version_falls_back_to_package_metadata(monkeypatch):
    monkeypatch.setattr(execution, "_package_repo_root", lambda: None)
    version = _resolve_version()
    # Package is installed editable; importlib.metadata should resolve it.
    assert version
    assert version != "unknown"


def test_dirty_is_string_state():
    assert _resolve_dirty() in {"clean", "dirty", "unknown"}


def test_dirty_unknown_when_not_a_git_checkout(monkeypatch):
    monkeypatch.setattr(execution, "_package_repo_root", lambda: None)
    assert _resolve_dirty() == "unknown"


def test_timestamp_has_iso_prefix_and_tz_suffix():
    ts = _resolve_timestamp()
    assert re.match(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\s+\S+$", ts), ts


def test_hostname_is_nonempty():
    assert _resolve_hostname()


def test_dirty_emits_warning(monkeypatch):
    monkeypatch.setattr(execution, "_resolve_dirty", lambda: "dirty")
    with pytest.warns(UserWarning, match="uncommitted changes"):
        result = gather_pipeline_provenance()
    # Stored as the 3-state string (matches init_* and backend_*).
    assert result[CONST.PIPELINE_DIRTY_ATTR] == "dirty"


def test_clean_tree_does_not_warn(monkeypatch):
    monkeypatch.setattr(execution, "_resolve_dirty", lambda: "clean")
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        result = gather_pipeline_provenance()
    assert result[CONST.PIPELINE_DIRTY_ATTR] == "clean"


def test_dirty_warning_fires_once_per_process(monkeypatch):
    """Second consecutive dirty gather is silent."""
    monkeypatch.setattr(execution, "_resolve_dirty", lambda: "dirty")
    with pytest.warns(UserWarning, match="uncommitted changes"):
        gather_pipeline_provenance()
    # No warning on the second call.
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        gather_pipeline_provenance()


def test_mark_dirty_warning_emitted_suppresses_subsequent(monkeypatch):
    """The CLI's explicit mark prevents the inline warning entirely."""
    monkeypatch.setattr(execution, "_resolve_dirty", lambda: "dirty")
    mark_dirty_warning_emitted()
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        gather_pipeline_provenance()
        gather_init_provenance()


def test_init_dirty_stored_as_string(monkeypatch):
    """Init attrs use the 3-state string, not bool."""
    monkeypatch.setattr(execution, "_resolve_dirty", lambda: "dirty")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        result = gather_init_provenance()
    assert result[CONST.INIT_DIRTY_ATTR] == "dirty"


def test_init_and_execution_share_dedup_flag(monkeypatch):
    """First gather (either kind) emits; second (either kind) is silent."""
    monkeypatch.setattr(execution, "_resolve_dirty", lambda: "dirty")
    with pytest.warns(UserWarning, match="uncommitted changes"):
        gather_init_provenance()
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        gather_pipeline_provenance()


def test_allow_dirty_env_truthy(monkeypatch):
    from lysis.tools.provenance import allow_dirty_from_env
    monkeypatch.setenv(CONST.LYSIS_ALLOW_DIRTY_ENV, "1")
    assert allow_dirty_from_env() is True
    monkeypatch.setenv(CONST.LYSIS_ALLOW_DIRTY_ENV, "yes")
    assert allow_dirty_from_env() is True
    monkeypatch.setenv(CONST.LYSIS_ALLOW_DIRTY_ENV, "FALSE")
    assert allow_dirty_from_env() is False
    monkeypatch.delenv(CONST.LYSIS_ALLOW_DIRTY_ENV, raising=False)
    assert allow_dirty_from_env() is False


def test_allow_commit_mismatch_env_truthy(monkeypatch):
    from lysis.tools.provenance import allow_commit_mismatch_from_env
    monkeypatch.setenv(CONST.LYSIS_ALLOW_COMMIT_MISMATCH_ENV, "on")
    assert allow_commit_mismatch_from_env() is True
    monkeypatch.delenv(CONST.LYSIS_ALLOW_COMMIT_MISMATCH_ENV, raising=False)
    assert allow_commit_mismatch_from_env() is False
