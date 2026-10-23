"""Compiler choice in :mod:`lysis.execution.historical_build` (#137).

Stubbed tests cover how the choice reaches ``make`` and how commits from
before compiler selection, whose Makefile hard-codes ifort, are handled.
The ``fortran_binary`` tests build for real: the current ``HEAD`` with each
compiler, and a commit from before the change with the default choice.
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

from lysis.execution import historical_build
from lysis.execution.historical_build import (
    DEFAULT_FORTRAN_COMPILER,
    FORTRAN_COMPILERS,
    HistoricalBuildError,
    _make_argv,
    _make_env,
    build_historical_binary,
)
from lysis.tools.provenance.binary import _identify_compiler_binary

#: Last ``main`` commit before compiler selection: its Makefile hard-codes
#: ifort.  Kept as a fixed SHA because "before #137" is a fixed point.
_PRE_SELECTION_COMMIT = "a6a626d9cbfed74afab27fe1ad62a8793d26ab86"

_NEW_MAKEFILE = "LYSIS_COMPILER_DEFAULT := ifx\n# uses LYSIS_COMPILER\n"
_LEGACY_MAKEFILE = "IFORT = ifort -r8\nFORT = $(IFORT)\n"


# ----------------------------------------------------------------------
# _make_argv / _make_env
# ----------------------------------------------------------------------


def test_make_argv_without_compiler_passes_nothing():
    assert _make_argv(Path("/b"), None) == ["make", "-C", "/b"]
    assert _make_argv(Path("/b"), None, None) == ["make", "-C", "/b"]


@pytest.mark.parametrize("compiler", FORTRAN_COMPILERS)
def test_make_argv_passes_compiler_on_command_line(compiler):
    assert _make_argv(Path("/b"), None, compiler) == [
        "make", "-C", "/b", f"LYSIS_COMPILER={compiler}",
    ]


def test_make_argv_module_wrapped_keeps_compiler():
    argv = _make_argv(Path("/b"), "intel-compilers/2024.2.0", "gfortran")
    assert argv[:2] == ["bash", "-c"]
    assert "module load intel-compilers/2024.2.0 &&" in argv[2]
    assert argv[2].endswith('"make" "-C" "/b" "LYSIS_COMPILER=gfortran"')


def test_make_env_drops_lysis_compiler(monkeypatch):
    monkeypatch.setenv("LYSIS_COMPILER", "gfortran")
    env = _make_env()
    assert "LYSIS_COMPILER" not in env
    assert env["PATH"] == os.environ["PATH"]


# ----------------------------------------------------------------------
# build_historical_binary with a stubbed archive and make
# ----------------------------------------------------------------------


@pytest.fixture
def stub_build(monkeypatch):
    """Stub ref resolution, archive and make; record the make call.

    ``state["makefile"]`` is written as the archived Makefile;
    ``state["log"]`` and ``state["rc"]`` are what the fake make produces.
    """
    state = {"makefile": _NEW_MAKEFILE, "log": "ok\n", "rc": 0, "calls": []}
    monkeypatch.setattr(
        historical_build, "resolve_ref", lambda ref, repo_root=None: "9" * 40
    )

    def fake_archive(sha, root, build_dir):
        if state["makefile"] is not None:
            (build_dir / "Makefile").write_text(state["makefile"])

    monkeypatch.setattr(historical_build, "_archive_into", fake_archive)

    def fake_run(argv, **kwargs):
        state["calls"].append((argv, kwargs))
        kwargs["stdout"].write(state["log"])
        kwargs["stdout"].flush()

        class _R:
            returncode = state["rc"]
        return _R()

    monkeypatch.setattr(historical_build.subprocess, "run", fake_run)
    return state


def _enter(**kwargs):
    with build_historical_binary("HEAD", "micro_rates", **kwargs):
        pytest.fail("the stub produces no binary")


def test_new_makefile_gets_requested_compiler(stub_build):
    with pytest.raises(HistoricalBuildError, match="produced no binary"):
        _enter(compiler="gfortran")
    argv, kwargs = stub_build["calls"][0]
    assert argv[-1] == "LYSIS_COMPILER=gfortran"
    assert "LYSIS_COMPILER" not in kwargs["env"]


def test_default_leaves_choice_to_makefile(stub_build):
    with pytest.raises(HistoricalBuildError, match="produced no binary"):
        _enter()
    argv, _ = stub_build["calls"][0]
    assert not any(a.startswith("LYSIS_COMPILER=") for a in argv)


def test_unknown_compiler_rejected_before_make(stub_build):
    with pytest.raises(HistoricalBuildError, match="valid values: ifx, ifort, gfortran"):
        _enter(compiler="nvfortran")
    assert stub_build["calls"] == []


@pytest.mark.parametrize("compiler", ["ifx", "gfortran"])
def test_legacy_commit_refuses_non_ifort(stub_build, compiler):
    stub_build["makefile"] = _LEGACY_MAKEFILE
    with pytest.raises(HistoricalBuildError) as exc:
        _enter(compiler=compiler)
    msg = str(exc.value)
    assert "predates Fortran compiler selection" in msg
    assert "hard-codes ifort" in msg
    assert "intel-compilers/2024.2.0" in msg
    assert stub_build["calls"] == []


@pytest.mark.parametrize("compiler", [None, "ifort"])
def test_legacy_commit_builds_with_ifort_choice(stub_build, compiler):
    """Default or explicit ifort: make runs, with nothing it would ignore."""
    stub_build["makefile"] = _LEGACY_MAKEFILE
    with pytest.raises(HistoricalBuildError, match="produced no binary"):
        _enter(compiler=compiler)
    argv, _ = stub_build["calls"][0]
    assert not any(a.startswith("LYSIS_COMPILER=") for a in argv)


@pytest.mark.parametrize("log", [
    "/bin/sh: ifort: command not found\n",
    "make: ifort: No such file or directory\n",
])
def test_legacy_commit_missing_ifort_explains(stub_build, log):
    stub_build["makefile"] = _LEGACY_MAKEFILE
    stub_build["log"] = log
    stub_build["rc"] = 2
    with pytest.raises(HistoricalBuildError) as exc:
        _enter(modules="intel-compilers/2026.1.0")
    msg = str(exc.value)
    assert "make exit 2" in msg
    assert "hard-codes ifort, which was not found in modules" in msg
    assert "'intel-compilers/2026.1.0'" in msg
    assert "intel-compilers/2024.2.0" in msg


def test_new_commit_failure_has_no_legacy_hint(stub_build):
    stub_build["log"] = "error: Fortran compiler 'ifx' not found on PATH\n"
    stub_build["rc"] = 2
    with pytest.raises(HistoricalBuildError) as exc:
        _enter()
    assert "predates" not in str(exc.value)


# ----------------------------------------------------------------------
# Real builds (need the compilers)
# ----------------------------------------------------------------------

#: Substring of ``compiler_version()`` that identifies each compiler.
_COMPILER_SIGNATURE = {
    "ifx": "Intel(R) Fortran Compiler for",
    "ifort": "Compiler Classic",
    "gfortran": "GCC version",
}


@pytest.mark.fortran_binary
@pytest.mark.parametrize("compiler", FORTRAN_COMPILERS)
def test_build_head_with_each_compiler(compiler):
    """HEAD builds with each compiler; provenance names the one that ran."""
    if shutil.which(compiler) is None:
        pytest.skip(f"{compiler} not on PATH")
    with build_historical_binary(
        "HEAD", "micro_rates", compiler=compiler, keep_dir=True
    ) as (binary, prov):
        build_dir = binary.parent.parent
        try:
            assert _identify_compiler_binary(build_dir / ".build_log") == compiler
            assert _COMPILER_SIGNATURE[compiler] in prov["backend_compiler"]
            if compiler != "ifort":
                assert "Classic" not in prov["backend_compiler"]
            assert prov["backend_historical"] is True
        finally:
            shutil.rmtree(build_dir, ignore_errors=True)


@pytest.mark.fortran_binary
def test_build_head_default_is_ifx():
    if shutil.which(DEFAULT_FORTRAN_COMPILER) is None:
        pytest.skip(f"{DEFAULT_FORTRAN_COMPILER} not on PATH")
    with build_historical_binary("HEAD", "micro_rates") as (_, prov):
        assert _COMPILER_SIGNATURE["ifx"] in prov["backend_compiler"]


def _have_pre_selection_commit():
    root = Path(historical_build.__file__).resolve().parents[3]
    return subprocess.run(
        ["git", "-C", str(root), "cat-file", "-e",
         f"{_PRE_SELECTION_COMMIT}^{{commit}}"],
        capture_output=True, check=False,
    ).returncode == 0


@pytest.mark.fortran_binary
def test_pre_selection_commit_builds_with_ifort_by_default():
    if not _have_pre_selection_commit():
        pytest.skip("pre-#137 commit not in this clone")
    if shutil.which("ifort") is None:
        pytest.skip("ifort not on PATH")
    with build_historical_binary(_PRE_SELECTION_COMMIT, "micro_rates") as (_, prov):
        assert prov["backend_commit"] == _PRE_SELECTION_COMMIT
        assert _COMPILER_SIGNATURE["ifort"] in prov["backend_compiler"]


@pytest.mark.fortran_binary
def test_pre_selection_commit_without_ifort_says_what_to_load(monkeypatch, tmp_path):
    """With no ifort on PATH the build fails and says which module ships it."""
    if not _have_pre_selection_commit():
        pytest.skip("pre-#137 commit not in this clone")
    # The system PATH alone: make, git and gcc, but (on the cluster) no
    # Intel compiler, which comes only from modules.
    system_path = os.pathsep.join(["/usr/bin", "/bin"])
    if shutil.which("ifort", path=system_path) is not None:
        pytest.skip("ifort is on the system PATH")
    for tool in ("make", "git", "gcc"):
        if shutil.which(tool, path=system_path) is None:
            pytest.skip(f"{tool} not on the system PATH")
    monkeypatch.setenv("PATH", system_path)
    with pytest.raises(HistoricalBuildError) as exc:
        with build_historical_binary(_PRE_SELECTION_COMMIT, "micro_rates"):
            pass
    msg = str(exc.value)
    assert "hard-codes ifort, which was not found on PATH" in msg
    assert "intel-compilers/2024.2.0" in msg
