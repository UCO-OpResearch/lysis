"""Tests for the Makefile's Fortran compiler selection (#137).

Each test runs the real ``Makefile`` with GNU make in a throwaway git
repository holding a copy of ``src/fortran``, ``src/c`` and ``.gitignore``.
The compilers are stub scripts placed first on ``PATH``: each answers
``--version`` with a fake version line and "compiles" by writing
``built by <name>: <args>`` to its ``-o`` file.  So the tests need no real
compiler, and every artifact says which compiler made it and with which flags.
"""

import os
import re
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

from lysis.execution.historical_build import (
    DEFAULT_FORTRAN_COMPILER,
    FORTRAN_COMPILERS,
)

_REPO_ROOT = Path(__file__).resolve().parent.parent
_BINARIES = (
    "micro_rates",
    "macro_diffuse_into_and_along__internal",
    "macro_diffuse_into_and_along__external",
)

pytestmark = pytest.mark.skipif(
    shutil.which("make") is None or shutil.which("git") is None,
    reason="needs GNU make and git",
)

_STUB = """#!/bin/sh
name=$(basename "$0")
if [ "$1" = "--version" ]; then
    [ -n "$STUB_MISSING" ] && [ "$STUB_MISSING" = "$name" ] && exit 127
    echo "$name (STUB) ${STUB_VERSION:-1.0}"
    exit 0
fi
out=""; prev=""
for a in "$@"; do
    [ "$prev" = "-o" ] && out="$a"
    prev="$a"
done
[ -n "$out" ] && echo "built by $name: $*" > "$out"
exit 0
"""


@pytest.fixture
def build_tree(tmp_path):
    """A git repo with the Makefile and sources, plus a dir of stub compilers."""
    repo = tmp_path / "repo"
    (repo / "src").mkdir(parents=True)
    shutil.copy(_REPO_ROOT / "Makefile", repo / "Makefile")
    shutil.copy(_REPO_ROOT / ".gitignore", repo / ".gitignore")
    shutil.copytree(
        _REPO_ROOT / "src" / "fortran", repo / "src" / "fortran",
        ignore=shutil.ignore_patterns("version_stamp.f90", "*.mod"),
    )
    shutil.copytree(_REPO_ROOT / "src" / "c", repo / "src" / "c")
    git = ["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t"]
    subprocess.run([*git, "init", "-q"], check=True)
    subprocess.run([*git, "add", "."], check=True)
    subprocess.run([*git, "commit", "-q", "-m", "fixture"], check=True)

    stubs = tmp_path / "stubs"
    stubs.mkdir()
    for name in (*FORTRAN_COMPILERS, "gcc"):
        path = stubs / name
        path.write_text(_STUB)
        path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return repo, stubs


def _make(build_tree, *args, env_extra=None):
    repo, stubs = build_tree
    env = {k: v for k, v in os.environ.items()
           if k not in ("LYSIS_COMPILER", "MAKEFLAGS", "MFLAGS", "MAKELEVEL")}
    env["PATH"] = f"{stubs}{os.pathsep}{env['PATH']}"
    env.update(env_extra or {})
    return subprocess.run(
        ["make", "-C", str(repo), *args],
        capture_output=True, text=True, env=env, check=False,
    )


def _built_by(build_tree):
    """Map each Fortran binary to the first word of its stub contents' compiler."""
    repo, _ = build_tree
    out = {}
    for name in _BINARIES:
        text = (repo / "bin" / name).read_text()
        out[name] = re.match(r"built by (\S+):", text).group(1)
    return out


def _mtimes(build_tree):
    repo, _ = build_tree
    return {name: (repo / "bin" / name).stat().st_mtime_ns for name in _BINARIES}


# ---------------------------------------------------------------------------
# Choosing the compiler
# ---------------------------------------------------------------------------


def test_default_is_ifx(build_tree):
    result = _make(build_tree)
    assert result.returncode == 0, result.stderr
    assert set(_built_by(build_tree).values()) == {"ifx"}
    assert "Fortran compiler: ifx (default)" in result.stdout


def test_python_default_matches_makefile():
    text = (_REPO_ROOT / "Makefile").read_text()
    default = re.search(r"^LYSIS_COMPILER_DEFAULT\s*:=\s*(\S+)", text, re.M).group(1)
    valid = re.search(r"^LYSIS_COMPILERS\s*:=\s*(.+)$", text, re.M).group(1).split()
    assert default == DEFAULT_FORTRAN_COMPILER
    assert tuple(valid) == FORTRAN_COMPILERS


@pytest.mark.parametrize("compiler", ["ifort", "gfortran"])
def test_command_line_overrides_default(build_tree, compiler):
    result = _make(build_tree, f"LYSIS_COMPILER={compiler}")
    assert result.returncode == 0, result.stderr
    assert set(_built_by(build_tree).values()) == {compiler}
    assert f"Fortran compiler: {compiler} (command line)" in result.stdout


def test_gfortran_flags_match_intel_semantics(build_tree):
    """8-byte reals and doubles, no line limit, FP traps, -O2."""
    assert _make(build_tree, "LYSIS_COMPILER=gfortran").returncode == 0
    repo, _ = build_tree
    for name in _BINARIES:
        args = (repo / "bin" / name).read_text().split()
        for flag in ("-fdefault-real-8", "-fdefault-double-8",
                     "-ffree-line-length-none",
                     "-ffpe-trap=invalid,zero,overflow", "-O2"):
            assert flag in args, f"{name}: {flag} missing"


def test_configure_saves_choice_for_plain_make(build_tree):
    repo, _ = build_tree
    result = _make(build_tree, "configure", "LYSIS_COMPILER=gfortran")
    assert result.returncode == 0, result.stderr
    assert (repo / "config.mk").is_file()
    result = _make(build_tree)
    assert result.returncode == 0, result.stderr
    assert set(_built_by(build_tree).values()) == {"gfortran"}
    assert "gfortran (saved in config.mk" in _make(build_tree, "help").stdout


def test_command_line_beats_saved_choice(build_tree):
    assert _make(build_tree, "configure", "LYSIS_COMPILER=gfortran").returncode == 0
    assert _make(build_tree, "LYSIS_COMPILER=ifort").returncode == 0
    assert set(_built_by(build_tree).values()) == {"ifort"}
    # The one-off choice is not saved: plain make goes back to gfortran.
    assert _make(build_tree).returncode == 0
    assert set(_built_by(build_tree).values()) == {"gfortran"}


def test_unconfigure_returns_to_default(build_tree):
    repo, _ = build_tree
    assert _make(build_tree, "configure", "LYSIS_COMPILER=gfortran").returncode == 0
    assert _make(build_tree, "unconfigure").returncode == 0
    assert not (repo / "config.mk").exists()
    assert _make(build_tree).returncode == 0
    assert set(_built_by(build_tree).values()) == {"ifx"}


def test_configure_requires_a_value(build_tree):
    repo, _ = build_tree
    result = _make(build_tree, "configure")
    assert result.returncode != 0
    assert "usage: make configure LYSIS_COMPILER=" in result.stderr
    assert not (repo / "config.mk").exists()


def test_environment_is_ignored(build_tree):
    result = _make(build_tree, env_extra={"LYSIS_COMPILER": "gfortran"})
    assert result.returncode == 0, result.stderr
    assert set(_built_by(build_tree).values()) == {"ifx"}
    assert "in the environment is ignored" in result.stderr


@pytest.mark.parametrize("args", [
    ["LYSIS_COMPILER=nvfortran"],
    ["configure", "LYSIS_COMPILER=nvfortran"],
    ["fort", "LYSIS_COMPILER=ifx gfortran"],
])
def test_unknown_value_fails_listing_valid_ones(build_tree, args):
    repo, _ = build_tree
    result = _make(build_tree, *args)
    assert result.returncode != 0
    assert "Unknown LYSIS_COMPILER" in result.stderr
    assert "valid values: ifx ifort gfortran" in result.stderr
    assert not (repo / "bin" / "micro_rates").exists()
    assert not (repo / "config.mk").exists()


def test_bad_saved_value_fails_but_can_be_inspected_and_removed(build_tree):
    repo, _ = build_tree
    (repo / "config.mk").write_text("LYSIS_SAVED_COMPILER := nvfortran\n")
    result = _make(build_tree)
    assert result.returncode != 0
    assert "Unknown LYSIS_COMPILER 'nvfortran'" in result.stderr
    assert _make(build_tree, "help").returncode == 0
    assert _make(build_tree, "unconfigure").returncode == 0
    assert not (repo / "config.mk").exists()


def test_missing_compiler_says_how_to_load_it(build_tree):
    result = _make(build_tree, env_extra={"STUB_MISSING": "ifx"})
    assert result.returncode != 0
    assert "'ifx' not found on PATH" in result.stderr
    assert "intel-compilers/2024.2.0" in result.stderr


# ---------------------------------------------------------------------------
# Rebuilding on a compiler change
# ---------------------------------------------------------------------------


def test_unchanged_compiler_rebuilds_nothing(build_tree):
    assert _make(build_tree).returncode == 0
    before = _mtimes(build_tree)
    assert _make(build_tree).returncode == 0
    assert _mtimes(build_tree) == before


def test_switching_compiler_rebuilds_every_binary(build_tree):
    assert _make(build_tree).returncode == 0
    before = _mtimes(build_tree)
    assert _make(build_tree, "LYSIS_COMPILER=gfortran").returncode == 0
    after = _mtimes(build_tree)
    assert all(after[n] != before[n] for n in _BINARIES)
    assert set(_built_by(build_tree).values()) == {"gfortran"}


def test_compiler_version_change_rebuilds(build_tree):
    assert _make(build_tree, env_extra={"STUB_VERSION": "1.0"}).returncode == 0
    before = _mtimes(build_tree)
    assert _make(build_tree, env_extra={"STUB_VERSION": "2.0"}).returncode == 0
    after = _mtimes(build_tree)
    assert all(after[n] != before[n] for n in _BINARIES)


def test_switching_keeps_tree_and_version_stamp_clean(build_tree):
    """config.mk and the compiler stamp stay untracked; the stamp reads clean."""
    repo, _ = build_tree
    for args in ([], ["LYSIS_COMPILER=gfortran"],
                 ["configure", "LYSIS_COMPILER=ifort"], []):
        assert _make(build_tree, *args).returncode == 0
    assert (repo / "config.mk").is_file()
    status = subprocess.run(
        ["git", "-C", str(repo), "status", "--porcelain"],
        capture_output=True, text=True, check=True,
    ).stdout
    assert status == ""
    stamp = (repo / "src" / "fortran" / "version_stamp.f90").read_text()
    assert 'BUILD_DIRTY  = "clean"' in stamp


@pytest.mark.parametrize("target", ["clean", "clobber"])
def test_clean_and_clobber_keep_saved_choice(build_tree, target):
    repo, _ = build_tree
    assert _make(build_tree, "configure", "LYSIS_COMPILER=gfortran").returncode == 0
    assert _make(build_tree).returncode == 0
    assert _make(build_tree, target).returncode == 0
    assert not (repo / "bin").exists()
    assert (repo / "config.mk").is_file()
