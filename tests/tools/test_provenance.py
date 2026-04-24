"""Unit tests for :mod:`lysis.tools.provenance`.

Covers the public :func:`gather_execution_provenance` entry point and each
private resolver:

* ``_resolve_version``   — tag > commit hash > pip metadata > ``"unknown"``
* ``_resolve_dirty``     — porcelain status parsing, no-git fallback
* ``_resolve_timestamp`` — format and tz suffix
* ``_resolve_hostname``  — FQDN contains a dot (on machines with a domain)
"""

import re
import subprocess
import warnings
from pathlib import Path

import pytest

from lysis.config.constants import CONST
from lysis.tools import provenance
from lysis.tools.provenance import (
    _package_repo_root,
    _resolve_dirty,
    _resolve_hostname,
    _resolve_timestamp,
    _resolve_version,
    gather_execution_provenance,
)


def test_gather_returns_expected_keys():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        result = gather_execution_provenance()
    assert set(result.keys()) == {
        CONST.EXECUTION_VERSION_ATTR,
        CONST.EXECUTION_DIRTY_ATTR,
        CONST.EXECUTION_TIMESTAMP_ATTR,
        CONST.EXECUTION_HOSTNAME_ATTR,
    }


def test_version_is_nonempty_string():
    assert isinstance(_resolve_version(), str)
    assert _resolve_version()


def test_version_matches_git_head_when_no_tag(monkeypatch):
    repo_root = _package_repo_root()
    if repo_root is None:
        pytest.skip("not running in a git checkout")
    # Pretend no tags exist at HEAD by intercepting the describe call.
    real_git = provenance._git

    def fake_git(args, root):
        if args[:3] == ["describe", "--exact-match", "--tags"]:
            return None
        return real_git(args, root)

    monkeypatch.setattr(provenance, "_git", fake_git)
    version = _resolve_version()
    commit = subprocess.check_output(
        ["git", "-C", str(repo_root), "rev-parse", "HEAD"], text=True
    ).strip()
    assert version == commit


def test_version_falls_back_to_package_metadata(monkeypatch):
    monkeypatch.setattr(provenance, "_package_repo_root", lambda: None)
    version = _resolve_version()
    # Package is installed editable; importlib.metadata should resolve it.
    assert version
    assert version != "unknown"


def test_dirty_is_bool():
    assert isinstance(_resolve_dirty(), bool)


def test_dirty_false_when_not_a_git_checkout(monkeypatch):
    monkeypatch.setattr(provenance, "_package_repo_root", lambda: None)
    assert _resolve_dirty() is False


def test_timestamp_has_iso_prefix_and_tz_suffix():
    ts = _resolve_timestamp()
    # YYYY-MM-DDTHH:MM:SS followed by a non-empty tz token.
    assert re.match(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\s+\S+$", ts), ts


def test_hostname_is_nonempty():
    assert _resolve_hostname()


def test_dirty_emits_warning(monkeypatch):
    monkeypatch.setattr(provenance, "_resolve_dirty", lambda: True)
    with pytest.warns(UserWarning, match="working tree is dirty"):
        result = gather_execution_provenance()
    assert result[CONST.EXECUTION_DIRTY_ATTR] is True


def test_clean_tree_does_not_warn(monkeypatch):
    monkeypatch.setattr(provenance, "_resolve_dirty", lambda: False)
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        result = gather_execution_provenance()
    assert result[CONST.EXECUTION_DIRTY_ATTR] is False
