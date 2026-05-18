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
    # No compiler token → compiler is "unknown".
    assert query_binary_version("/fake/bin") == ("abc123", "clean", "unknown")


def test_query_binary_version_dirty_state(monkeypatch):
    _patch_run(monkeypatch, stdout="abc123 dirty\n")
    assert query_binary_version("/fake/bin") == ("abc123", "dirty", "unknown")


def test_query_binary_version_nonzero_exit_is_unknown(monkeypatch):
    _patch_run(monkeypatch, returncode=2, stdout="abc123 clean\n")
    assert query_binary_version("/fake/bin") == ("unknown", "unknown", "unknown")


def test_query_binary_version_empty_output_is_unknown(monkeypatch):
    _patch_run(monkeypatch, stdout="")
    assert query_binary_version("/fake/bin") == ("unknown", "unknown", "unknown")


def test_query_binary_version_too_few_fields_is_unknown(monkeypatch):
    _patch_run(monkeypatch, stdout="abc123\n")
    assert query_binary_version("/fake/bin") == ("unknown", "unknown", "unknown")


def test_query_binary_version_returns_compiler_suffix(monkeypatch):
    # Real binaries emit "<commit> <state> <compiler-info>" where the
    # compiler string itself contains spaces (e.g. "Intel(R) Fortran ...").
    # The parser must accept the line and return the compiler suffix
    # verbatim as the third tuple element.
    _patch_run(
        monkeypatch,
        stdout="abc123 clean Intel(R) Fortran Classic 2021.9.0\n",
    )
    assert query_binary_version("/fake/bin") == (
        "abc123", "clean", "Intel(R) Fortran Classic 2021.9.0",
    )


def test_query_binary_version_file_not_found_is_unknown(monkeypatch):
    _patch_run(monkeypatch, raises=FileNotFoundError())
    assert query_binary_version("/fake/bin") == ("unknown", "unknown", "unknown")


def test_query_binary_version_timeout_is_unknown(monkeypatch):
    _patch_run(monkeypatch, raises=subprocess.TimeoutExpired(cmd="x", timeout=1))
    assert query_binary_version("/fake/bin") == ("unknown", "unknown", "unknown")


# ----------------------------------------------------------------------
# gather_binary_provenance
# ----------------------------------------------------------------------

def test_gather_binary_provenance_keys(monkeypatch):
    from lysis.tools.provenance import gather_binary_provenance

    _patch_run(monkeypatch, stdout="abc123 clean Intel\n")
    result = gather_binary_provenance("/fake/bin")
    assert result == {
        CONST.BINARY_COMMIT_ATTR: "abc123",
        CONST.BINARY_DIRTY_ATTR: "clean",
        CONST.BINARY_COMPILER_ATTR: "Intel",
    }


def test_gather_binary_provenance_failure_yields_unknown(monkeypatch):
    from lysis.tools.provenance import gather_binary_provenance

    _patch_run(monkeypatch, raises=FileNotFoundError())
    result = gather_binary_provenance("/fake/bin")
    assert result == {
        CONST.BINARY_COMMIT_ATTR: "unknown",
        CONST.BINARY_DIRTY_ATTR: "unknown",
        CONST.BINARY_COMPILER_ATTR: "unknown",
    }


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
    """Stub out both the binary query and source-tree query.

    *binary* is a 3-tuple ``(commit, dirty, compiler)``.  *source* is the
    2-tuple ``(commit, dirty)`` returned by
    :func:`gather_fortran_source_provenance`.
    """
    monkeypatch.setattr(
        binary_mod, "query_binary_version", lambda exe, **kw: binary,
    )
    monkeypatch.setattr(
        binary_mod, "gather_fortran_source_provenance", lambda: source,
    )


def test_verify_match_returns_empty(monkeypatch):
    _patch_check(
        monkeypatch,
        binary=("abc123", "clean", "Intel"),
        source=("abc123", "clean"),
    )
    assert verify_binary_matches_source("/fake/bin", allow_stale=False) == {}


def test_verify_match_dirty_both_sides_returns_empty(monkeypatch):
    _patch_check(
        monkeypatch,
        binary=("abc123", "dirty", "Intel"),
        source=("abc123", "dirty"),
    )
    assert verify_binary_matches_source("/fake/bin", allow_stale=False) == {}


def test_verify_commit_mismatch_raises(monkeypatch):
    _patch_check(
        monkeypatch,
        binary=("abc123", "clean", "Intel"),
        source=("def456", "clean"),
    )
    with pytest.raises(StaleBinaryError) as excinfo:
        verify_binary_matches_source("/fake/bin", allow_stale=False)
    assert "abc123" in str(excinfo.value)
    assert "def456" in str(excinfo.value)


def test_verify_dirty_mismatch_raises(monkeypatch):
    _patch_check(
        monkeypatch,
        binary=("abc123", "clean", "Intel"),
        source=("abc123", "dirty"),
    )
    with pytest.raises(StaleBinaryError):
        verify_binary_matches_source("/fake/bin", allow_stale=False)


def test_verify_unknown_binary_raises(monkeypatch):
    _patch_check(
        monkeypatch,
        binary=("unknown", "unknown", "unknown"),
        source=("abc123", "clean"),
    )
    with pytest.raises(StaleBinaryError):
        verify_binary_matches_source("/fake/bin", allow_stale=False)


def test_verify_unknown_source_raises(monkeypatch):
    _patch_check(
        monkeypatch,
        binary=("unknown", "unknown", "unknown"),
        source=("unknown", "unknown"),
    )
    with pytest.raises(StaleBinaryError):
        verify_binary_matches_source("/fake/bin", allow_stale=False)


def test_verify_override_warns_and_returns_dict(monkeypatch):
    _patch_check(
        monkeypatch,
        binary=("abc123", "clean", "Intel"),
        source=("def456", "clean"),
    )
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
    # Only the override flag is forwarded as an HDF5 attr — the binary's
    # commit/dirty/compiler are now stamped unconditionally by
    # gather_binary_provenance(), not via this return value.
    assert info[CONST.STALE_BINARY_OVERRIDE_ATTR] is True
    assert CONST.BINARY_COMMIT_ATTR not in info
    assert CONST.BINARY_DIRTY_ATTR not in info


def test_verify_env_var_acts_as_override(monkeypatch):
    _patch_check(
        monkeypatch,
        binary=("abc123", "clean", "Intel"),
        source=("def456", "clean"),
    )
    monkeypatch.setenv(CONST.LYSIS_ALLOW_STALE_BINARY_ENV, "1")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        info = verify_binary_matches_source("/fake/bin", allow_stale=None)
    assert info[CONST.STALE_BINARY_OVERRIDE_ATTR] is True


def test_verify_explicit_false_overrides_env(monkeypatch):
    # allow_stale=False must take precedence over the env var.
    _patch_check(
        monkeypatch,
        binary=("abc123", "clean", "Intel"),
        source=("def456", "clean"),
    )
    monkeypatch.setenv(CONST.LYSIS_ALLOW_STALE_BINARY_ENV, "1")
    with pytest.raises(StaleBinaryError):
        verify_binary_matches_source("/fake/bin", allow_stale=False)
