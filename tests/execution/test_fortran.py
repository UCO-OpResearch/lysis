"""Unit tests for :mod:`lysis.execution.fortran` — FortranRunner.

Tests exercise :meth:`FortranRunner.exec_command` via a minimal concrete
stub subclass that supplies the parameter-class hooks.  This lets us test
the shared template logic (unit conversion, seed splitting, Fortran-name
suffixes) in isolation from the macro/micro specialisations.
"""

import pytest
import numpy as np

from lysis.config.constants import Q_
from lysis.config.parameters import MicroParameters
from lysis.config.run import Run
from lysis.execution.fortran import FortranRunner


# ---------------------------------------------------------------------------
# Minimal concrete stub
# ---------------------------------------------------------------------------


def _make_stub_class(params_override=None):
    """Return a FortranRunner subclass using MicroParameters as the param class."""

    class _StubRunner(FortranRunner):
        """Concrete FortranRunner stub backed by MicroParameters."""

        def _get_params(self):
            return self.run.micro_params

        def _get_params_class(self):
            return MicroParameters

        def _seed_field(self):
            return "micro_seed"

        def _simulations_field(self):
            return "micro_simulations"

        def _seed_split_count(self, params):
            return self.index + 1

        def _base_arguments(self):
            return ["--runCode", self.run.run_code, "--outFileCode", self.out_file_code]

        def _log_prefix(self):
            return "micro"

    return _StubRunner


@pytest.fixture
def tmp_run(tmp_path):
    """Run with default MicroParameters."""
    r = Run(str(tmp_path))
    r.initialize_micro_param()
    return r


@pytest.fixture
def stub_runner(tmp_run):
    """Default-params stub runner."""
    cls = _make_stub_class()
    return cls(run=tmp_run, executable="/bin/stub.exe")


# ---------------------------------------------------------------------------
# TestFortranRunnerAbstract
# ---------------------------------------------------------------------------


class TestFortranRunnerAbstract:
    """FortranRunner itself must be abstract."""

    def test_cannot_instantiate_directly(self, tmp_run):
        with pytest.raises(TypeError):
            FortranRunner(run=tmp_run, executable="/bin/exe")

    def test_stub_instantiates(self, stub_runner):
        assert isinstance(stub_runner, FortranRunner)


# ---------------------------------------------------------------------------
# TestExecCommandTemplate
# ---------------------------------------------------------------------------


class TestExecCommandTemplate:
    """Tests for the shared exec_command() template method."""

    def test_executable_is_first_element(self, stub_runner):
        cmd = stub_runner.exec_command()
        assert cmd[0] == "/bin/stub.exe"

    def test_run_code_present(self, stub_runner):
        cmd = stub_runner.exec_command()
        assert "--runCode" in cmd
        assert stub_runner.run.run_code in cmd

    def test_out_file_code_present(self, stub_runner):
        cmd = stub_runner.exec_command()
        assert "--outFileCode" in cmd

    def test_default_params_only_base_args(self, stub_runner):
        """With all-default params only base args should appear."""
        cmd = stub_runner.exec_command()
        # base: executable + --runCode + code + --outFileCode + ""  = 5 elements
        assert len(cmd) == 5

    def test_non_default_param_included(self, tmp_path):
        """A non-default param must appear as a CLI flag."""
        r = Run(str(tmp_path))
        r.initialize_micro_param({"micro_simulations": 999})
        cls = _make_stub_class()
        runner = cls(run=r, executable="/bin/stub.exe")
        cmd = runner.exec_command()
        assert "--simulations" in cmd
        assert "999" in cmd

    def test_default_param_excluded(self, stub_runner):
        """A param equal to its default must not appear."""
        cmd = stub_runner.exec_command()
        assert "--simulations" not in cmd

    def test_quantity_unit_conversion(self, tmp_path):
        """Pint Quantity params must be converted to the expected SI unit."""
        r = Run(str(tmp_path))
        r.initialize_micro_param({"fiber_radius": Q_("37 nanometer")})
        cls = _make_stub_class()
        runner = cls(run=r, executable="/bin/stub.exe")
        cmd = runner.exec_command()
        assert "--radius" in cmd
        val = float(cmd[cmd.index("--radius") + 1])
        # MicroParameters units["fiber_radius"] is microns; 37 nm = 0.037 µm
        assert pytest.approx(val, rel=1e-6) == 0.037

    def test_minus_one_suffix_adjusts_index(self, tmp_path):
        """Params with a '-1' Fortran-name suffix must have 1 added to their value."""
        import inspect as _inspect

        sig = _inspect.signature(MicroParameters)
        fname_map = MicroParameters.fortran_names()
        # Only consider keys that are actually MicroParameters fields
        minus_one_keys = [
            k for k, v in fname_map.items()
            if v.endswith("-1") and k in sig.parameters
        ]
        if not minus_one_keys:
            pytest.skip("No '-1' suffix params in MicroParameters")

        key = minus_one_keys[0]
        default_val = sig.parameters[key].default
        override_val = default_val + 1

        r = Run(str(tmp_path))
        r.initialize_micro_param({key: override_val})
        cls = _make_stub_class()
        runner = cls(run=r, executable="/bin/stub.exe")
        cmd = runner.exec_command()

        fortran_flag = "--" + fname_map[key][:-2]
        assert fortran_flag in cmd
        idx = cmd.index(fortran_flag)
        assert cmd[idx + 1] == str(override_val + 1)

    # ------------------------------------------------------------------
    # Seed splitting
    # ------------------------------------------------------------------

    def test_index_none_no_split(self, stub_runner):
        stub_runner.index = None
        cmd = stub_runner.exec_command()
        assert "--simulations" not in cmd

    def test_index_zero_forces_simulations_to_one(self, tmp_run):
        cls = _make_stub_class()
        runner = cls(run=tmp_run, executable="/bin/stub.exe", index=0)
        cmd = runner.exec_command()
        assert "--simulations" in cmd
        assert cmd[cmd.index("--simulations") + 1] == "1"

    def test_index_zero_appends_suffix(self, tmp_run):
        cls = _make_stub_class()
        runner = cls(run=tmp_run, executable="/bin/stub.exe", index=0)
        runner.exec_command()
        assert runner.out_file_code == "__00"

    def test_index_two_appends_correct_suffix(self, tmp_run):
        cls = _make_stub_class()
        runner = cls(run=tmp_run, executable="/bin/stub.exe", index=2)
        runner.exec_command()
        assert runner.out_file_code == "__02"

    def test_index_two_correct_seed(self, tmp_run):
        seed = tmp_run.micro_params.micro_seed
        stream = np.random.SeedSequence(seed)
        expected_seed = int(np.int32(stream.generate_state(3)[2]))

        cls = _make_stub_class()
        runner = cls(run=tmp_run, executable="/bin/stub.exe", index=2)
        cmd = runner.exec_command()
        assert "--seed" in cmd
        assert cmd[cmd.index("--seed") + 1] == str(expected_seed)

    def test_no_index_out_file_code_unchanged(self, tmp_run):
        cls = _make_stub_class()
        runner = cls(run=tmp_run, executable="/bin/stub.exe",
                     out_file_code="orig", index=None)
        runner.exec_command()
        assert runner.out_file_code == "orig"
