"""Unit tests for :mod:`lysis.execution.fortran_macro` — FortranMacro."""

from pathlib import Path
from unittest.mock import call, patch

import numpy as np
import pytest

from lysis.config.constants import Q_
from lysis.config.parameters import MacroParameters, MicroParameters
from lysis.config.run import Run
from lysis.execution.fortran_macro import FortranMacro


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def tmp_run(tmp_path):
    """A Run with default Micro and Macro parameters."""
    r = Run(str(tmp_path))
    r.initialize_micro_param()
    r.initialize_macro_param()
    return r


@pytest.fixture
def fortran_macro(tmp_run):
    """A FortranMacro with default run and dummy executable."""
    return FortranMacro(run=tmp_run, executable="/bin/macro.exe")


# ---------------------------------------------------------------------------
# TestFortranMacroExecCommand
# ---------------------------------------------------------------------------


class TestFortranMacroExecCommand:
    """Tests for :meth:`FortranMacro.exec_command`."""

    def test_executable_first(self, fortran_macro):
        cmd = fortran_macro.exec_command()
        assert cmd[0] == "/bin/macro.exe"

    def test_includes_run_code(self, fortran_macro):
        cmd = fortran_macro.exec_command()
        assert "--runCode" in cmd
        assert fortran_macro.run.run_code in cmd

    def test_includes_in_file_code(self, tmp_run):
        fm = FortranMacro(run=tmp_run, executable="/bin/macro.exe",
                          in_file_code="_in")
        cmd = fm.exec_command()
        assert "--inFileCode" in cmd
        idx = cmd.index("--inFileCode")
        assert cmd[idx + 1] == "_in"

    def test_includes_out_file_code(self, tmp_run):
        fm = FortranMacro(run=tmp_run, executable="/bin/macro.exe",
                          out_file_code="_out")
        cmd = fm.exec_command()
        assert "--outFileCode" in cmd
        idx = cmd.index("--outFileCode")
        assert cmd[idx + 1] == "_out"

    def test_radius_appended(self, fortran_macro):
        """--radius from micro_params must appear in the command."""
        cmd = fortran_macro.exec_command()
        assert "--radius" in cmd

    def test_bs_appended(self, fortran_macro):
        """--bs (binding sites) from micro_params must appear in the command."""
        cmd = fortran_macro.exec_command()
        assert "--bs" in cmd

    def test_radius_value_correct(self, fortran_macro):
        """--radius must carry the fiber_radius in the expected SI unit."""
        cmd = fortran_macro.exec_command()
        idx = cmd.index("--radius")
        expected = fortran_macro.run.micro_params.fiber_radius.m_as(
            MicroParameters.units()["fiber_radius"]
        )
        assert pytest.approx(float(cmd[idx + 1]), rel=1e-6) == expected

    def test_seed_split_uses_macro_simulations(self, tmp_run):
        """Seed split count for FortranMacro must equal macro_simulations."""
        n_sims = tmp_run.macro_params.macro_simulations
        seed = tmp_run.macro_params.macro_seed

        stream = np.random.SeedSequence(seed)
        expected_seed = int(np.int32(stream.generate_state(n_sims)[0]))

        fm = FortranMacro(run=tmp_run, executable="/bin/macro.exe", index=0)
        cmd = fm.exec_command()
        # After split, --simulations should be 1
        assert "--simulations" in cmd
        assert cmd[cmd.index("--simulations") + 1] == "1"
        # Seed value must match
        assert "--seed" in cmd
        assert cmd[cmd.index("--seed") + 1] == str(expected_seed)

    def test_index_appends_suffix(self, tmp_run):
        fm = FortranMacro(run=tmp_run, executable="/bin/macro.exe", index=3)
        fm.exec_command()
        assert fm.out_file_code == "__03"


# ---------------------------------------------------------------------------
# TestFortranMacroGenerateNeighborhoods
# ---------------------------------------------------------------------------


class TestFortranMacroGenerateNeighborhoods:
    """Tests for :meth:`FortranMacro.generate_neighborhoods`."""

    def test_creates_neighbors_file(self, fortran_macro, tmp_path):
        """generate_neighborhoods() must create neighbors.dat in run.os_path."""
        fortran_macro.generate_neighborhoods()
        assert (tmp_path / "neighbors.dat").exists()

    def test_file_is_non_empty(self, fortran_macro, tmp_path):
        fortran_macro.generate_neighborhoods()
        assert (tmp_path / "neighbors.dat").stat().st_size > 0

    def test_file_contains_integers(self, fortran_macro, tmp_path):
        """neighbors.dat must be newline-delimited integers (1-based)."""
        fortran_macro.generate_neighborhoods()
        content = (tmp_path / "neighbors.dat").read_text()
        lines = [ln for ln in content.splitlines() if ln.strip()]
        # All lines should parse as integers
        values = [int(ln) for ln in lines]
        # All values should be >= 1 (Fortran 1-based indexing)
        assert all(v >= 1 for v in values)


# ---------------------------------------------------------------------------
# TestFortranMacroExecute
# ---------------------------------------------------------------------------


class TestFortranMacroExecute:
    """Tests for :meth:`FortranMacro.execute`."""

    def test_generate_neighborhoods_called(self, fortran_macro, tmp_path):
        """execute() must call generate_neighborhoods() before the subprocess."""
        with (
            patch.object(FortranMacro, "generate_neighborhoods") as mock_gen,
            patch("subprocess.run"),
        ):
            fortran_macro.execute()
        assert mock_gen.call_count == 1

    def test_subprocess_called_once(self, fortran_macro, tmp_path):
        with (
            patch.object(FortranMacro, "generate_neighborhoods"),
            patch("subprocess.run") as mock_run,
        ):
            fortran_macro.execute()
        assert mock_run.call_count == 1

    def test_log_file_created_in_os_path(self, fortran_macro, tmp_path):
        """The macro log file must be created in run.os_path."""
        with (
            patch.object(FortranMacro, "generate_neighborhoods"),
            patch("subprocess.run"),
        ):
            fortran_macro.execute()
        log_files = list(tmp_path.glob("macro*.txt"))
        assert len(log_files) == 1
