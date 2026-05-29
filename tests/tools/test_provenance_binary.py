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
    # Output is "<tag> <commit> <state> [compiler]"; the tag is dropped and
    # with no compiler token the compiler is "unknown".
    _patch_run(monkeypatch, stdout="v0.3.0 abc123 clean\n")
    assert query_binary_version("/fake/bin") == ("abc123", "clean", "unknown")


def test_query_binary_version_dirty_state(monkeypatch):
    _patch_run(monkeypatch, stdout="v0.3.0 abc123 dirty\n")
    assert query_binary_version("/fake/bin") == ("abc123", "dirty", "unknown")


def test_query_binary_version_nonzero_exit_is_unknown(monkeypatch):
    _patch_run(monkeypatch, returncode=2, stdout="v0.3.0 abc123 clean\n")
    assert query_binary_version("/fake/bin") == ("unknown", "unknown", "unknown")


def test_query_binary_version_empty_output_is_unknown(monkeypatch):
    _patch_run(monkeypatch, stdout="")
    assert query_binary_version("/fake/bin") == ("unknown", "unknown", "unknown")


def test_query_binary_version_too_few_fields_is_unknown(monkeypatch):
    # Only tag + commit (no dirty state) → fewer than three tokens → unknown.
    _patch_run(monkeypatch, stdout="v0.3.0 abc123\n")
    assert query_binary_version("/fake/bin") == ("unknown", "unknown", "unknown")


def test_query_binary_version_returns_compiler_suffix(monkeypatch):
    # Real binaries emit "<tag> <commit> <state> <compiler-info>" where the
    # compiler string itself contains spaces (e.g. "Intel(R) Fortran ...").
    # The parser must drop the tag and return the compiler suffix verbatim
    # as the third tuple element.
    _patch_run(
        monkeypatch,
        stdout="v0.3.0 abc123 clean Intel(R) Fortran Classic 2021.9.0\n",
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

    _patch_run(monkeypatch, stdout="v0.3.0 abc123 clean Intel\n")
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


# ----------------------------------------------------------------------
# verify_binary_matches_source — source_stamp override
# ----------------------------------------------------------------------

def _sentinel_source_provenance():
    raise AssertionError(
        "gather_fortran_source_provenance must not be called when "
        "source_stamp is provided"
    )


def test_verify_source_stamp_match_skips_git_lookup(monkeypatch):
    # When source_stamp matches the binary, we accept without calling git.
    monkeypatch.setattr(
        binary_mod, "query_binary_version",
        lambda exe, **kw: ("abc123", "clean", "Intel"),
    )
    monkeypatch.setattr(
        binary_mod, "gather_fortran_source_provenance",
        _sentinel_source_provenance,
    )
    result = verify_binary_matches_source(
        "/fake/bin",
        allow_stale=False,
        source_stamp=("abc123", "clean"),
    )
    assert result == {}


def test_verify_source_stamp_mismatch_raises(monkeypatch):
    monkeypatch.setattr(
        binary_mod, "query_binary_version",
        lambda exe, **kw: ("abc123", "clean", "Intel"),
    )
    monkeypatch.setattr(
        binary_mod, "gather_fortran_source_provenance",
        _sentinel_source_provenance,
    )
    with pytest.raises(StaleBinaryError) as excinfo:
        verify_binary_matches_source(
            "/fake/bin",
            allow_stale=False,
            source_stamp=("def456", "clean"),
        )
    assert "abc123" in str(excinfo.value)
    assert "def456" in str(excinfo.value)


def test_verify_source_stamp_none_falls_back_to_git_lookup(monkeypatch):
    # The fallback path (source_stamp=None) must still call
    # gather_fortran_source_provenance.
    calls = {"n": 0}

    def fake_source():
        calls["n"] += 1
        return ("abc123", "clean")

    monkeypatch.setattr(
        binary_mod, "query_binary_version",
        lambda exe, **kw: ("abc123", "clean", "Intel"),
    )
    monkeypatch.setattr(
        binary_mod, "gather_fortran_source_provenance", fake_source,
    )
    assert verify_binary_matches_source(
        "/fake/bin", allow_stale=False, source_stamp=None,
    ) == {}
    assert calls["n"] == 1


def test_verify_source_stamp_dirty_match(monkeypatch):
    monkeypatch.setattr(
        binary_mod, "query_binary_version",
        lambda exe, **kw: ("abc123", "dirty", "Intel"),
    )
    monkeypatch.setattr(
        binary_mod, "gather_fortran_source_provenance",
        _sentinel_source_provenance,
    )
    assert verify_binary_matches_source(
        "/fake/bin",
        allow_stale=False,
        source_stamp=("abc123", "dirty"),
    ) == {}


# ----------------------------------------------------------------------
# gather_historical_binary_provenance
# ----------------------------------------------------------------------

from lysis.tools.provenance import gather_historical_binary_provenance
from lysis.tools.provenance.binary import (
    _identify_compiler_binary,
    _probe_iso_fortran_env_compiler_version,
)


def test_identify_compiler_gfortran(tmp_path):
    log = tmp_path / "build.log"
    log.write_text("gfortran -mcmodel=medium -fbacktrace src/fortran/foo.f90\n")
    assert _identify_compiler_binary(log) == "gfortran"


def test_identify_compiler_ifort(tmp_path):
    log = tmp_path / "build.log"
    log.write_text(
        "ifort -r8 -mcmodel medium -traceback src/fortran/foo.f90\n"
        "gcc -std=c99 src/c/kiss.c\n"
    )
    assert _identify_compiler_binary(log) == "ifort"


def test_identify_compiler_ifx_prefers_first_known(tmp_path):
    log = tmp_path / "build.log"
    log.write_text("ifx src/fortran/foo.f90\nifort other.f90\n")
    assert _identify_compiler_binary(log) == "ifx"


def test_identify_compiler_no_log_returns_none():
    assert _identify_compiler_binary(None) is None


def test_identify_compiler_no_match_returns_none(tmp_path):
    log = tmp_path / "build.log"
    log.write_text("make: nothing to do\n")
    assert _identify_compiler_binary(log) is None


def test_identify_compiler_missing_file_returns_none(tmp_path):
    assert _identify_compiler_binary(tmp_path / "nope.log") is None


def test_historical_provenance_uses_binary_version_on_sha_match(monkeypatch):
    sha = "a" * 40
    _patch_run(monkeypatch, stdout=f"v0.3.0 {sha} clean Intel(R) 2021.9\n")
    result = gather_historical_binary_provenance(
        "/fake/bin", resolved_sha=sha
    )
    assert result == {
        CONST.BINARY_COMMIT_ATTR: sha,
        CONST.BINARY_DIRTY_ATTR: "clean",
        CONST.BINARY_COMPILER_ATTR: "Intel(R) 2021.9",
        CONST.BINARY_SOURCE_ATTR: f"historical:{sha}",
    }


def test_historical_provenance_synthesises_when_version_missing(
    monkeypatch, tmp_path
):
    sha = "b" * 40
    _patch_run(monkeypatch, raises=FileNotFoundError())
    log = tmp_path / "build.log"
    log.write_text("gfortran -c src/fortran/foo.f90\n")
    # Stub the iso_fortran_env probe so the test doesn't need a real compiler.
    monkeypatch.setattr(
        binary_mod,
        "_probe_iso_fortran_env_compiler_version",
        lambda compiler, **kw: "GCC version 11.4.0",
    )
    result = gather_historical_binary_provenance(
        "/fake/bin", resolved_sha=sha, build_log=log
    )
    assert result[CONST.BINARY_COMMIT_ATTR] == sha
    assert result[CONST.BINARY_DIRTY_ATTR] == "clean"
    assert result[CONST.BINARY_COMPILER_ATTR] == "GCC version 11.4.0"
    assert result[CONST.BINARY_SOURCE_ATTR] == f"historical:{sha}"


def test_historical_provenance_synthesises_when_sha_mismatches(
    monkeypatch, tmp_path
):
    requested = "c" * 40
    embedded = "d" * 40
    _patch_run(monkeypatch, stdout=f"v0.3.0 {embedded} clean ifort\n")
    log = tmp_path / "build.log"
    log.write_text("ifort src/fortran/foo.f90\n")
    monkeypatch.setattr(
        binary_mod,
        "_probe_iso_fortran_env_compiler_version",
        lambda compiler, **kw: "Intel(R) Fortran 2023.0",
    )
    result = gather_historical_binary_provenance(
        "/fake/bin", resolved_sha=requested, build_log=log
    )
    # Synthesised values, not the embedded ones.
    assert result[CONST.BINARY_COMMIT_ATTR] == requested
    assert result[CONST.BINARY_COMPILER_ATTR] == "Intel(R) Fortran 2023.0"
    assert result[CONST.BINARY_SOURCE_ATTR] == f"historical:{requested}"


def test_historical_provenance_compiler_unknown_when_log_absent(monkeypatch):
    sha = "e" * 40
    _patch_run(monkeypatch, raises=FileNotFoundError())
    result = gather_historical_binary_provenance(
        "/fake/bin", resolved_sha=sha, build_log=None
    )
    assert result[CONST.BINARY_COMPILER_ATTR] == "unknown"
    assert result[CONST.BINARY_SOURCE_ATTR] == f"historical:{sha}"


def test_historical_provenance_compiler_unknown_when_probe_fails(
    monkeypatch, tmp_path
):
    sha = "f" * 40
    _patch_run(monkeypatch, raises=FileNotFoundError())
    log = tmp_path / "build.log"
    log.write_text("gfortran src/fortran/foo.f90\n")
    monkeypatch.setattr(
        binary_mod,
        "_probe_iso_fortran_env_compiler_version",
        lambda compiler, **kw: None,
    )
    result = gather_historical_binary_provenance(
        "/fake/bin", resolved_sha=sha, build_log=log
    )
    assert result[CONST.BINARY_COMPILER_ATTR] == "unknown"


@pytest.mark.fortran_binary
def test_probe_iso_fortran_env_with_real_gfortran():
    # Smoke test that the probe actually compiles and runs against
    # whichever gfortran happens to be on PATH.  Skipped on machines
    # without it (the fortran_binary marker is documented as
    # "requires a compiled Fortran executable").
    import shutil as _sh
    if _sh.which("gfortran") is None:
        pytest.skip("gfortran not on PATH")
    out = _probe_iso_fortran_env_compiler_version("gfortran")
    assert out is not None
    assert "GCC" in out or "gfortran" in out.lower()
