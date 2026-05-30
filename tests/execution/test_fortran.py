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
            return self.num_children if self.num_children is not None else self.index + 1

        def _base_arguments(self):
            return ["--runCode", self.run.run_code, "--outFileCode", self.out_file_code]

        def _log_prefix(self):
            return "micro"

        def _collection_name(self):
            return "microscale_out"

        def _fortran_dataspec_version(self):
            from lysis.execution.fortran import MICRO_FORTRAN_DATASPEC_VERSION
            return MICRO_FORTRAN_DATASPEC_VERSION

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

    # ------------------------------------------------------------------
    # Array partition (num_children set)
    # ------------------------------------------------------------------

    def test_partition_clean_division(self, tmp_path):
        """N=50000, k=10 → every task gets exactly 5000."""
        r = Run(str(tmp_path))
        r.initialize_micro_param({"micro_simulations": 50000})
        cls = _make_stub_class()
        for i in range(10):
            runner = cls(run=r, executable="/bin/stub.exe", index=i, num_children=10)
            cmd = runner.exec_command()
            assert "--simulations" in cmd, f"task {i} missing --simulations"
            assert cmd[cmd.index("--simulations") + 1] == "5000"

    def test_partition_with_remainder(self, tmp_path):
        """N=53, k=10 → tasks 0..2 get 6, tasks 3..9 get 5; total = N."""
        r = Run(str(tmp_path))
        r.initialize_micro_param({"micro_simulations": 53})
        cls = _make_stub_class()
        per_task = []
        for i in range(10):
            runner = cls(run=r, executable="/bin/stub.exe", index=i, num_children=10)
            cmd = runner.exec_command()
            assert "--simulations" in cmd
            per_task.append(int(cmd[cmd.index("--simulations") + 1]))

        assert per_task[:3] == [6, 6, 6]
        assert per_task[3:] == [5] * 7
        assert sum(per_task) == 53

    def test_num_children_none_preserves_legacy_one_sim(self, tmp_run):
        """index set, num_children=None → simulations=1 (legacy path)."""
        cls = _make_stub_class()
        runner = cls(
            run=tmp_run, executable="/bin/stub.exe", index=2, num_children=None
        )
        cmd = runner.exec_command()
        assert "--simulations" in cmd
        assert cmd[cmd.index("--simulations") + 1] == "1"

    def test_seed_split_count_returns_num_children_when_set(self, tmp_run):
        """When num_children is set, every sibling draws from the same-sized stream."""
        seed = tmp_run.micro_params.micro_seed
        stream = np.random.SeedSequence(seed)
        seeds = stream.generate_state(10)

        cls = _make_stub_class()
        for i in range(10):
            runner = cls(run=tmp_run, executable="/bin/stub.exe",
                         index=i, num_children=10)
            cmd = runner.exec_command()
            assert "--seed" in cmd
            expected = int(np.int32(seeds[i]))
            assert cmd[cmd.index("--seed") + 1] == str(expected)

    def test_partition_with_num_children_one(self, tmp_path):
        """num_children=1 → the single task gets all N simulations."""
        r = Run(str(tmp_path))
        r.initialize_micro_param({"micro_simulations": 42})
        cls = _make_stub_class()
        runner = cls(run=r, executable="/bin/stub.exe", index=0, num_children=1)
        cmd = runner.exec_command()
        assert cmd[cmd.index("--simulations") + 1] == "42"

    # ------------------------------------------------------------------
    # |uint32 seed tag (Python uint32 → signed INTEGER*4 CLI token)
    # ------------------------------------------------------------------

    def test_high_bit_seed_cli_arg_is_signed_int32(self, tmp_path):
        """A uint32 seed with the high bit set is emitted as the signed int32 decimal."""
        high_bit_seed = np.uint32(0x80000001)
        r = Run(str(tmp_path))
        r.initialize_micro_param({"micro_seed": high_bit_seed})
        cls = _make_stub_class()
        # index=None → no SeedSequence split, seed flows through as-is.
        runner = cls(run=r, executable="/bin/stub.exe", index=None)
        cmd = runner.exec_command()

        expected = str(
            int(np.array(high_bit_seed, dtype=np.uint32).astype(np.int32))
        )
        assert expected == "-2147483647"
        assert "--seed" in cmd
        assert cmd[cmd.index("--seed") + 1] == expected

    def test_deadbeef_seed_cli_arg_matches_fortran_signed(self, tmp_path):
        """0xDEADBEEF uint32 emits -559038737 (Fortran's signed int32 view)."""
        r = Run(str(tmp_path))
        r.initialize_micro_param({"micro_seed": np.uint32(0xDEADBEEF)})
        cls = _make_stub_class()
        runner = cls(run=r, executable="/bin/stub.exe", index=None)
        cmd = runner.exec_command()
        assert "--seed" in cmd
        assert cmd[cmd.index("--seed") + 1] == "-559038737"


# ---------------------------------------------------------------------------
# skip_binary_verification (historical-build path)
# ---------------------------------------------------------------------------


class TestSkipBinaryVerification:
    """Behaviour of the FortranRunner.skip_binary_verification field."""

    def test_default_runs_the_check(self, monkeypatch, tmp_path):
        """Without the flag, _verify_binary_version calls into the provenance helper."""
        import lysis.tools.provenance as prov_pkg
        called = []
        monkeypatch.setattr(
            prov_pkg,
            "verify_binary_matches_source",
            lambda exe, **kw: called.append(exe) or {},
        )
        r = Run(str(tmp_path))
        r.initialize_micro_param()
        cls = _make_stub_class()
        runner = cls(run=r, executable="/bin/stub.exe")
        runner._verify_binary_version()
        assert called == ["/bin/stub.exe"]

    def test_skip_returns_empty_dict_no_subprocess(
        self, monkeypatch, tmp_path
    ):
        """skip_binary_verification=True must NOT invoke the provenance helper."""
        import lysis.tools.provenance as prov_pkg
        called = []
        monkeypatch.setattr(
            prov_pkg,
            "verify_binary_matches_source",
            lambda exe, **kw: called.append(exe) or {},
        )
        r = Run(str(tmp_path))
        r.initialize_micro_param()
        cls = _make_stub_class()
        runner = cls(
            run=r,
            executable="/bin/does-not-exist",
            skip_binary_verification=True,
        )
        info = runner._verify_binary_version()
        assert info == {}
        assert called == []
        # No banner attrs either.
        assert runner._backend_hdf5_attrs == {}


# ---------------------------------------------------------------------------
# source_stamp (Slurm-master copy-out path)
# ---------------------------------------------------------------------------


class TestSourceStampPlumbing:
    """A runner constructed with a ``source_stamp`` must hand it to
    ``verify_binary_matches_source`` so the binary check on a compute
    node does not need to invoke git from outside the source repo.
    """

    def test_runner_forwards_source_stamp(self, monkeypatch, tmp_path):
        import lysis.tools.provenance as prov_pkg
        captured = {}

        def fake_verify(exe, **kw):
            captured["exe"] = exe
            captured["kw"] = kw
            return {}

        monkeypatch.setattr(
            prov_pkg, "verify_binary_matches_source", fake_verify,
        )
        r = Run(str(tmp_path))
        r.initialize_micro_param()
        cls = _make_stub_class()
        runner = cls(
            run=r,
            executable="/bin/stub.exe",
            source_stamp=("abc123", "clean"),
        )
        runner._verify_binary_version()
        assert captured["exe"] == "/bin/stub.exe"
        assert captured["kw"]["source_stamp"] == ("abc123", "clean")

    def test_runner_default_source_stamp_is_none(
        self, monkeypatch, tmp_path,
    ):
        import lysis.tools.provenance as prov_pkg
        captured = {}

        def fake_verify(exe, **kw):
            captured["kw"] = kw
            return {}

        monkeypatch.setattr(
            prov_pkg, "verify_binary_matches_source", fake_verify,
        )
        r = Run(str(tmp_path))
        r.initialize_micro_param()
        cls = _make_stub_class()
        runner = cls(run=r, executable="/bin/stub.exe")
        runner._verify_binary_version()
        assert captured["kw"]["source_stamp"] is None
