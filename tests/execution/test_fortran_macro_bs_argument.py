"""``--bs`` is passed to the macroscale binary as an integer (#137).

``bs`` is an ``INTEGER`` in the Fortran programs, but ``binding_sites`` is
computed as a real (426.69... µM by default).  ifx and ifort truncate a real
literal read into an integer; gfortran rejects it, and the program then
stops without writing output.  Python therefore passes the truncated value,
which is what the Intel builds have always used.
"""

import math

import pytest

from lysis.config.parameters import MicroParameters
from lysis.config.run import Run
from lysis.execution.fortran_macro import FortranMacro


@pytest.fixture(autouse=True)
def _stub_preflight_binary_check(monkeypatch):
    monkeypatch.setattr(FortranMacro, "_verify_binary_version", lambda self: {})


@pytest.fixture
def default_macro(tmp_path):
    run = Run(str(tmp_path))
    run.initialize_micro_param()
    run.initialize_macro_param()
    return FortranMacro(run=run, executable="/bin/macro.exe")


def _bs_value(cmd):
    return cmd[cmd.index("--bs") + 1]


def test_bs_is_an_integer_literal(default_macro):
    value = _bs_value(default_macro.exec_command())
    assert value.lstrip("-").isdigit(), value


def test_bs_is_truncated_not_rounded(default_macro):
    """The default 426.69... becomes 426, as ifx and ifort read it."""
    units = MicroParameters.units()
    exact = default_macro.run.micro_params.binding_sites.m_as(units["binding_sites"])
    assert exact != int(exact)  # the default really is non-integral
    assert _bs_value(default_macro.exec_command()) == str(math.trunc(exact))
    assert _bs_value(default_macro.exec_command()) == "426"
