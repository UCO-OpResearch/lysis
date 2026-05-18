"""Unit tests for :mod:`lysis.tools.provenance.binary`.

Exercises each layer of the preflight check:

* :func:`query_binary_version` — parses ``--version`` output and treats
  every failure mode as ``("unknown", "unknown")``.
* :func:`gather_fortran_source_provenance` — wraps git rev-parse and
  ``git diff --quiet HEAD -- src/fortran/``.
* :func:`allow_stale_from_env` — env-var truthiness rules.
* :func:`verify_binary_matches_source` — match/mismatch/override matrix
  and the banner + HDF5 attrs returned on overridden mismatch.
"""

import subprocess
import warnings

import pytest

from lysis.config.constants import CONST
from lysis.tools.provenance import (
    StaleBinaryError,
    allow_stale_from_env,
    gather_fortran_source_provenance,
    query_binary_version,
    verify_binary_matches_source,
)
from lysis.tools.provenance import binary as binary_mod


# ----------------------------------------------------------------------
# query_binary_version
# ----------------------------------------------------------------------

class _FakeCompleted:
    def __init__(self, returncode=0, stdout=""):
        self.returncode = returncode
        self.stdout = stdout


def _patch_run(monkeypatch, *, returncode=0, stdout="", raises=None):
    def fake_run(*args, **kwargs):
        if raises is not None:
            raise raises
        return _FakeCompleted(returncode=returncode, stdout=stdout)
    monkeypatch.setattr(binary_mod.subprocess, "run", fake_run)


def test_query_binary_version_happy_path(monkeypatch):
    _patch_run(monkeypatch, stdout="abc123 clean\n")
    assert query_binary_version("/fake/bin") == ("abc123", "clean")


def test_query_binary_version_dirty_state(monkeypatch):
    _patch_run(monkeypatch, stdout="abc123 dirty\n")
    assert query_binary_version("/fake/bin") == ("abc123", "dirty")


def test_query_binary_version_nonzero_exit_is_unknown(monkeypatch):
    _patch_run(monkeypatch, returncode=2, stdout="abc123 clean\n")
    assert query_binary_version("/fake/bin") == ("unknown", "unknown")


def test_query_binary_version_empty_output_is_unknown(monkeypatch):
    _patch_run(monkeypatch, stdout="")
    assert query_binary_version("/fake/bin") == ("unknown", "unknown")


def test_query_binary_version_too_few_fields_is_unknown(monkeypatch):
    _patch_run(monkeypatch, stdout="abc123\n")
    assert query_binary_version("/fake/bin") == ("unknown", "unknown")


def test_query_binary_version_too_many_fields_is_unknown(monkeypatch):
    _patch_run(monkeypatch, stdout="abc123 clean extra\n")
    assert query_binary_version("/fake/bin") == ("unknown", "unknown")


def test_query_binary_version_file_not_found_is_unknown(monkeypatch):
    _patch_run(monkeypatch, raises=FileNotFoundError())
    assert query_binary_version("/fake/bin") == ("unknown", "unknown")


def test_query_binary_version_timeout_is_unknown(monkeypatch):
    _patch_run(monkeypatch, raises=subprocess.TimeoutExpired(cmd="x", timeout=1))
    assert query_binary_version("/fake/bin") == ("unknown", "unknown")


# ----------------------------------------------------------------------
# gather_fortran_source_provenance
# ----------------------------------------------------------------------

def test_source_provenance_no_repo_is_unknown(monkeypatch):
    monkeypatch.setattr(binary_mod, "_package_repo_root", lambda: None)
    assert gather_fortran_source_provenance() == ("unknown", "unknown")


def test_source_provenance_clean(monkeypatch, tmp_path):
    monkeypatch.setattr(binary_mod, "_package_repo_root", lambda: tmp_path)
    monkeypatch.setattr(binary_mod, "_git", lambda args, root: "abc123")
    # `git diff --quiet` exit 0 means clean.
    _patch_run(monkeypatch, returncode=0)
    assert gather_fortran_source_provenance() == ("abc123", "clean")


def test_source_provenance_dirty(monkeypatch, tmp_path):
    monkeypatch.setattr(binary_mod, "_package_repo_root", lambda: tmp_path)
    monkeypatch.setattr(binary_mod, "_git", lambda args, root: "abc123")
    # `git diff --quiet` exit 1 means dirty.
    _patch_run(monkeypatch, returncode=1)
    assert gather_fortran_source_provenance() == ("abc123", "dirty")


def test_source_provenance_git_unavailable_is_unknown(monkeypatch, tmp_path):
    monkeypatch.setattr(binary_mod, "_package_repo_root", lambda: tmp_path)
    monkeypatch.setattr(binary_mod, "_git", lambda args, root: None)
    assert gather_fortran_source_provenance() == ("unknown", "unknown")


# ----------------------------------------------------------------------
# allow_stale_from_env
# ----------------------------------------------------------------------

@pytest.mark.parametrize("val", ["1", "true", "TRUE", "yes", "YES", "on", "On"])
def test_allow_stale_truthy(monkeypatch, val):
    monkeypatch.setenv(CONST.LYSIS_ALLOW_STALE_BINARY_ENV, val)
    assert allow_stale_from_env() is True


@pytest.mark.parametrize("val", ["0", "false", "no", "off", "", "maybe"])
def test_allow_stale_falsy(monkeypatch, val):
    monkeypatch.setenv(CONST.LYSIS_ALLOW_STALE_BINARY_ENV, val)
    assert allow_stale_from_env() is False


def test_allow_stale_unset_is_false(monkeypatch):
    monkeypatch.delenv(CONST.LYSIS_ALLOW_STALE_BINARY_ENV, raising=False)
    assert allow_stale_from_env() is False


# ----------------------------------------------------------------------
# verify_binary_matches_source
# ----------------------------------------------------------------------

def _patch_check(monkeypatch, *, binary, source):
    """Stub out both the binary query and source-tree query."""
    monkeypatch.setattr(
        binary_mod, "query_binary_version", lambda exe, **kw: binary,
    )
    monkeypatch.setattr(
        binary_mod, "gather_fortran_source_provenance", lambda: source,
    )


def test_verify_match_returns_empty(monkeypatch):
    _patch_check(monkeypatch, binary=("abc123", "clean"), source=("abc123", "clean"))
    assert verify_binary_matches_source("/fake/bin", allow_stale=False) == {}


def test_verify_match_dirty_both_sides_returns_empty(monkeypatch):
    # If both source and binary are dirty AND on the same commit, that's
    # still a match — we trust the user knows what they're doing locally.
    _patch_check(monkeypatch, binary=("abc123", "dirty"), source=("abc123", "dirty"))
    assert verify_binary_matches_source("/fake/bin", allow_stale=False) == {}


def test_verify_commit_mismatch_raises(monkeypatch):
    _patch_check(monkeypatch, binary=("abc123", "clean"), source=("def456", "clean"))
    with pytest.raises(StaleBinaryError) as excinfo:
        verify_binary_matches_source("/fake/bin", allow_stale=False)
    assert "abc123" in str(excinfo.value)
    assert "def456" in str(excinfo.value)


def test_verify_dirty_mismatch_raises(monkeypatch):
    # Binary built from a clean commit, but source has uncommitted edits.
    _patch_check(monkeypatch, binary=("abc123", "clean"), source=("abc123", "dirty"))
    with pytest.raises(StaleBinaryError):
        verify_binary_matches_source("/fake/bin", allow_stale=False)


def test_verify_unknown_binary_raises(monkeypatch):
    _patch_check(
        monkeypatch, binary=("unknown", "unknown"), source=("abc123", "clean"),
    )
    with pytest.raises(StaleBinaryError):
        verify_binary_matches_source("/fake/bin", allow_stale=False)


def test_verify_unknown_source_raises(monkeypatch):
    # Even matched unknowns count as mismatch — we never claim a clean
    # match without a real commit to anchor on.
    _patch_check(
        monkeypatch, binary=("unknown", "unknown"), source=("unknown", "unknown"),
    )
    with pytest.raises(StaleBinaryError):
        verify_binary_matches_source("/fake/bin", allow_stale=False)


def test_verify_override_warns_and_returns_dict(monkeypatch):
    _patch_check(monkeypatch, binary=("abc123", "clean"), source=("def456", "clean"))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", UserWarning)
        info = verify_binary_matches_source("/fake/bin", allow_stale=True)
    # Loud warning emitted.
    assert any("abc123" in str(w.message) and "def456" in str(w.message)
               for w in caught)
    # Banner present and contains both stamps.
    assert "banner" in info
    assert "abc123" in info["banner"]
    assert "def456" in info["banner"]
    assert "STALE FORTRAN BINARY" in info["banner"]
    # HDF5 attrs present.
    assert info[CONST.STALE_BINARY_OVERRIDE_ATTR] is True
    assert info[CONST.BINARY_COMMIT_ATTR] == "abc123"
    assert info[CONST.BINARY_DIRTY_ATTR] == "clean"


def test_verify_env_var_acts_as_override(monkeypatch):
    _patch_check(monkeypatch, binary=("abc123", "clean"), source=("def456", "clean"))
    monkeypatch.setenv(CONST.LYSIS_ALLOW_STALE_BINARY_ENV, "1")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        info = verify_binary_matches_source("/fake/bin", allow_stale=None)
    assert info[CONST.STALE_BINARY_OVERRIDE_ATTR] is True


def test_verify_explicit_false_overrides_env(monkeypatch):
    # allow_stale=False must take precedence over the env var.
    _patch_check(monkeypatch, binary=("abc123", "clean"), source=("def456", "clean"))
    monkeypatch.setenv(CONST.LYSIS_ALLOW_STALE_BINARY_ENV, "1")
    with pytest.raises(StaleBinaryError):
        verify_binary_matches_source("/fake/bin", allow_stale=False)
