"""Tests for ``--fortran-compiler`` on ``lysis run-micro`` and ``lysis run-macro``.

The option chooses the compiler for a ``--fortran-commit`` build (#137).
The build itself is stubbed; these tests check only that the choice reaches
:func:`~lysis.execution.historical_build.build_historical_binary`, that the
default leaves the choice to the commit's Makefile, and that the option is
refused without ``--fortran-commit``.
"""

from contextlib import contextmanager
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from lysis.cli import cli
from lysis.execution.historical_build import FORTRAN_COMPILERS
from tests.cli.test_run_macro import _write_macro_hdf5
from tests.cli.test_run_micro import _write_micro_hdf5

_PROV = {
    "backend_commit": "a" * 40,
    "backend_dirty": "clean",
    "backend_compiler": "stub",
    "backend_historical": True,
}

# (command, FortranX class to stub, HDF5 writer)
_COMMANDS = [
    ("run-micro", "lysis.execution.fortran_micro.FortranMicro", _write_micro_hdf5),
    ("run-macro", "lysis.execution.fortran_macro.FortranMacro", _write_macro_hdf5),
]


@pytest.fixture
def runner():
    return CliRunner()


def _invoke_with_stubbed_build(runner, tmp_path, command, cls_path, writer, extra):
    """Run *command* with a stubbed historical build; return (result, build kwargs)."""
    h5 = tmp_path / "run-01.h5"
    writer(h5)
    fake_binary = tmp_path / "bin" / "micro_rates"
    fake_binary.parent.mkdir(parents=True)
    fake_binary.write_text("not a real binary")
    seen = {}

    @contextmanager
    def fake_cm(*args, **kwargs):
        seen.update(kwargs)
        yield fake_binary, _PROV

    with patch(cls_path) as mock_cls, patch(
        "lysis.execution.historical_build.build_historical_binary",
        side_effect=fake_cm,
    ):
        mock_cls.from_hdf5.return_value = MagicMock()
        result = runner.invoke(
            cli,
            [command, str(h5), "--executable", "micro_rates",
             "--fortran-commit", "HEAD", *extra],
        )
    return result, seen


@pytest.mark.parametrize("command, cls_path, writer", _COMMANDS)
@pytest.mark.parametrize("compiler", FORTRAN_COMPILERS)
def test_compiler_reaches_historical_build(
    runner, tmp_path, command, cls_path, writer, compiler
):
    result, seen = _invoke_with_stubbed_build(
        runner, tmp_path, command, cls_path, writer,
        ["--fortran-compiler", compiler],
    )
    assert result.exit_code == 0, result.output
    assert seen.get("compiler") == compiler


@pytest.mark.parametrize("command, cls_path, writer", _COMMANDS)
def test_default_leaves_choice_to_makefile(runner, tmp_path, command, cls_path, writer):
    """No option -> compiler=None, so the commit's own Makefile default applies."""
    result, seen = _invoke_with_stubbed_build(
        runner, tmp_path, command, cls_path, writer, [],
    )
    assert result.exit_code == 0, result.output
    assert "compiler" in seen
    assert seen["compiler"] is None


@pytest.mark.parametrize("command, cls_path, writer", _COMMANDS)
def test_unknown_compiler_rejected(runner, tmp_path, command, cls_path, writer):
    result, _ = _invoke_with_stubbed_build(
        runner, tmp_path, command, cls_path, writer,
        ["--fortran-compiler", "nvfortran"],
    )
    assert result.exit_code != 0
    assert "nvfortran" in result.output


@pytest.mark.parametrize("command, cls_path, writer", _COMMANDS)
def test_requires_fortran_commit(runner, tmp_path, command, cls_path, writer):
    h5 = tmp_path / "run-01.h5"
    writer(h5)
    with patch(cls_path):
        result = runner.invoke(
            cli,
            [command, str(h5), "--executable", "bin/micro_rates",
             "--fortran-compiler", "gfortran"],
        )
    assert result.exit_code != 0
    assert "--fortran-commit" in result.output
    assert "LYSIS_COMPILER" in result.output


@pytest.mark.parametrize("command", ["run-micro", "run-macro"])
def test_modules_default_is_new_default_module(runner, command):
    """--modules defaults to intel-compilers/2024.2.0, which ships ifx and ifort."""
    from lysis.tools.slurm import DEFAULT_MODULES

    assert DEFAULT_MODULES == "intel-compilers/2024.2.0"
    params = {p.name: p for p in cli.commands[command].params}
    assert params["modules"].default == DEFAULT_MODULES
    assert params["fortran_compiler"].default is None
    result = runner.invoke(cli, [command, "--help"])
    assert result.exit_code == 0
    assert "--fortran-compiler" in result.output
