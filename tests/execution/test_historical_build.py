"""Tests for :mod:`lysis.execution.historical_build`.

Most behaviour is exercised with stubs — the real ``git archive``/``make``
roundtrip needs a Fortran compiler and several seconds, so it lives behind
the ``fortran_binary`` marker and is skipped when no compiler is on PATH.
"""

import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from lysis.execution import historical_build
from lysis.execution.historical_build import (
    HistoricalBuildError,
    _resolve_binary_name,
    build_historical_binary,
    resolve_ref,
)


# ----------------------------------------------------------------------
# _resolve_binary_name
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("micro_rates", "micro_rates"),
        ("bin/micro_rates", "micro_rates"),
        ("/bin/macro_diffuse_into_and_along__internal",
         "macro_diffuse_into_and_along__internal"),
    ],
)
def test_resolve_binary_name_strips_prefix(raw, expected):
    assert _resolve_binary_name(raw) == expected


@pytest.mark.parametrize("bad", ["a/b/c", "", ".", "..", "sub/bin/micro"])
def test_resolve_binary_name_rejects_paths(bad):
    with pytest.raises(HistoricalBuildError, match="binary name"):
        _resolve_binary_name(bad)


# ----------------------------------------------------------------------
# resolve_ref
# ----------------------------------------------------------------------


def test_resolve_ref_head_returns_sha():
    sha = resolve_ref("HEAD")
    assert isinstance(sha, str)
    assert len(sha) == 40
    assert all(c in "0123456789abcdef" for c in sha)


def test_resolve_ref_bad_raises():
    with pytest.raises(HistoricalBuildError, match="does not resolve"):
        resolve_ref("this-ref-does-not-exist-deadbeef")


def test_resolve_ref_no_repo_raises(monkeypatch):
    monkeypatch.setattr(
        historical_build, "_package_repo_root", lambda: None
    )
    with pytest.raises(HistoricalBuildError, match="repository root"):
        resolve_ref("HEAD")


# ----------------------------------------------------------------------
# build_historical_binary — error paths (don't need a real toolchain)
# ----------------------------------------------------------------------


def test_build_failure_surfaces_log_tail(monkeypatch, tmp_path):
    # Pretend ref resolves and archive succeeds; make `make` return
    # non-zero so the failure path runs.
    fake_sha = "c" * 40
    monkeypatch.setattr(
        historical_build, "resolve_ref", lambda ref, repo_root=None: fake_sha
    )
    monkeypatch.setattr(
        historical_build, "_archive_into", lambda sha, root, build_dir: None
    )

    def fake_make_subprocess(*args, **kwargs):
        # Write something to the log via the file handle the caller passes
        # in (kwargs["stdout"]) so the tail can be read back.
        kwargs["stdout"].write("FATAL: compiler not found\n")
        kwargs["stdout"].flush()

        class _R:
            returncode = 2
        return _R()

    monkeypatch.setattr(historical_build.subprocess, "run", fake_make_subprocess)

    with pytest.raises(HistoricalBuildError) as exc_info:
        with build_historical_binary("HEAD", "micro_rates"):
            pytest.fail("Should not enter the context on build failure")
    assert "make exit 2" in str(exc_info.value)
    assert "FATAL: compiler not found" in str(exc_info.value)


def test_build_missing_binary_raises(monkeypatch, tmp_path):
    fake_sha = "d" * 40
    monkeypatch.setattr(
        historical_build, "resolve_ref", lambda ref, repo_root=None: fake_sha
    )
    monkeypatch.setattr(
        historical_build, "_archive_into", lambda sha, root, build_dir: None
    )

    # `make` "succeeds" but produces nothing in bin/.
    def fake_make(*args, **kwargs):
        kwargs["stdout"].write("ok\n")

        class _R:
            returncode = 0
        return _R()

    monkeypatch.setattr(historical_build.subprocess, "run", fake_make)

    with pytest.raises(HistoricalBuildError, match="produced no binary"):
        with build_historical_binary("HEAD", "micro_rates"):
            pytest.fail("Should not enter the context when binary missing")


def test_build_dir_cleaned_up_on_failure(monkeypatch):
    # Capture the build dir that gets created so we can confirm it's gone
    # after the failure path runs.
    captured = []
    real_mkdtemp = historical_build.tempfile.mkdtemp

    def spy_mkdtemp(*args, **kwargs):
        d = real_mkdtemp(*args, **kwargs)
        captured.append(Path(d))
        return d

    monkeypatch.setattr(historical_build.tempfile, "mkdtemp", spy_mkdtemp)
    monkeypatch.setattr(
        historical_build, "resolve_ref", lambda ref, repo_root=None: "e" * 40
    )
    monkeypatch.setattr(
        historical_build, "_archive_into", lambda sha, root, bd: None
    )

    def fake_make(*args, **kwargs):
        kwargs["stdout"].write("nope\n")

        class _R:
            returncode = 1
        return _R()

    monkeypatch.setattr(historical_build.subprocess, "run", fake_make)

    with pytest.raises(HistoricalBuildError):
        with build_historical_binary("HEAD", "micro_rates"):
            pass

    assert len(captured) == 1
    assert not captured[0].exists()


def test_build_dir_preserved_when_keep_dir(monkeypatch):
    captured = []
    real_mkdtemp = historical_build.tempfile.mkdtemp

    def spy_mkdtemp(*args, **kwargs):
        d = real_mkdtemp(*args, **kwargs)
        captured.append(Path(d))
        return d

    monkeypatch.setattr(historical_build.tempfile, "mkdtemp", spy_mkdtemp)
    monkeypatch.setattr(
        historical_build, "resolve_ref", lambda ref, repo_root=None: "f" * 40
    )
    monkeypatch.setattr(
        historical_build, "_archive_into", lambda sha, root, bd: None
    )

    def fake_make(*args, **kwargs):
        kwargs["stdout"].write("nope\n")

        class _R:
            returncode = 1
        return _R()

    monkeypatch.setattr(historical_build.subprocess, "run", fake_make)

    with pytest.raises(HistoricalBuildError):
        with build_historical_binary("HEAD", "micro_rates", keep_dir=True):
            pass

    assert len(captured) == 1
    assert captured[0].exists()
    # Clean up so we don't leave test debris.
    import shutil
    shutil.rmtree(captured[0], ignore_errors=True)


# ----------------------------------------------------------------------
# End-to-end smoke: build current HEAD as "historical" (needs a compiler)
# ----------------------------------------------------------------------


@pytest.mark.fortran_binary
def test_build_head_as_historical_end_to_end():
    """Build whichever micro binary the current Makefile produces.

    Marker note: requires a Fortran compiler on PATH.  Skipped when none
    is found, so it stays runnable on dev machines without compilers
    installed.
    """
    import shutil as _sh
    if not any(_sh.which(c) for c in ("ifx", "ifort", "gfortran")):
        pytest.skip("no Fortran compiler on PATH")

    with build_historical_binary("HEAD", "micro_rates") as (
        binary_path, provenance
    ):
        assert binary_path.is_file()
        assert provenance["binary_source"].startswith("historical:")
        # binary_compiler may come from --version (if HEAD's binary supports
        # it) OR from the iso_fortran_env probe — either is acceptable.
        assert provenance["binary_compiler"] != "unknown"
