"""Unit tests for :mod:`lysis.execution.fortran_macro` — FortranMacro."""

from pathlib import Path
from unittest.mock import MagicMock, call, patch

import h5py
import numpy as np
import pytest

from lysis.config.constants import CONST, Q_
from lysis.config.parameters import MacroParameters, MicroParameters
from lysis.config.run import Run
from lysis.execution.fortran import MACRO_FORTRAN_DATASPEC_VERSION
from lysis.execution.fortran_macro import FortranMacro


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_macro_hdf5(path: Path) -> None:
    """Write a minimal v2.0.0 HDF5 file with default Micro and MacroParameters."""
    mp = MicroParameters()
    mcp = MacroParameters(micro_params=mp)
    with h5py.File(str(path), "w") as f:
        f.attrs[CONST.DATASPEC_VERSION_ATTR] = "v2.0.0"
        micro_grp = f.require_group("micro_data")
        for k, v in mp.to_basedict().items():
            micro_grp.attrs[k] = str(v) if not isinstance(v, (int, float, bool)) else v
        macro_grp = f.require_group("macro_data")
        for k, v in mcp.to_basedict().items():
            macro_grp.attrs[k] = str(v) if not isinstance(v, (int, float, bool)) else v


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

    def test_output_matches_legacy_tofile_content(self, tmp_run, tmp_path):
        """write_dataset output must contain the same integers as the old tofile path.

        The new implementation uses np.savetxt (one integer per line) instead
        of ndarray.tofile(sep=os.linesep).  Both write integer values in the
        same row-major order; the only allowed difference is a trailing newline.
        """
        import os
        import numpy as np
        from lysis.geometry.edge_grid import generate_fortran_neighborhood_structure

        fm = FortranMacro(run=tmp_run, executable="/bin/macro.exe")

        # --- new path ---
        fm.generate_neighborhoods()
        new_content = (tmp_path / "neighbors.dat").read_text()
        new_values = [int(ln) for ln in new_content.splitlines() if ln.strip()]

        # --- legacy path ---
        fort_neighbors = (
            generate_fortran_neighborhood_structure(
                tmp_run.macro_params.rows, tmp_run.macro_params.cols
            )
            + 1
        )
        legacy_str = os.linesep.join(str(v) for v in fort_neighbors.flatten())
        legacy_values = [int(v) for v in legacy_str.split() if v.strip()]

        assert new_values == legacy_values


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


# ---------------------------------------------------------------------------
# TestFortranMacroFromHdf5
# ---------------------------------------------------------------------------


@pytest.fixture
def macro_hdf5(tmp_path):
    """Minimal HDF5 file with both Micro and MacroParameters."""
    h5_path = tmp_path / "run-01.h5"
    _write_macro_hdf5(h5_path)
    return h5_path


class TestFortranMacroFromHdf5:
    """Tests for :meth:`FortranMacro.from_hdf5`."""

    def test_returns_fortran_macro_instance(self, macro_hdf5):
        fm = FortranMacro.from_hdf5(macro_hdf5, "/bin/macro.exe")
        assert isinstance(fm, FortranMacro)

    def test_sets_run_run_code(self, macro_hdf5):
        fm = FortranMacro.from_hdf5(macro_hdf5, "/bin/macro.exe")
        assert fm.run.run_code == "run-01"

    def test_sets_micro_params(self, macro_hdf5):
        fm = FortranMacro.from_hdf5(macro_hdf5, "/bin/macro.exe")
        assert fm.run.micro_params is not None

    def test_sets_macro_params(self, macro_hdf5):
        fm = FortranMacro.from_hdf5(macro_hdf5, "/bin/macro.exe")
        assert fm.run.macro_params is not None

    def test_sets_in_file_code(self, macro_hdf5):
        fm = FortranMacro.from_hdf5(macro_hdf5, "/bin/macro.exe",
                                     in_file_code="_in")
        assert fm.in_file_code == "_in"

    def test_sets_out_file_code(self, macro_hdf5):
        fm = FortranMacro.from_hdf5(macro_hdf5, "/bin/macro.exe",
                                     out_file_code="_out")
        assert fm.out_file_code == "_out"

    def test_sets_index(self, macro_hdf5):
        fm = FortranMacro.from_hdf5(macro_hdf5, "/bin/macro.exe", index=3)
        assert fm.index == 3

    def test_index_defaults_to_none(self, macro_hdf5):
        fm = FortranMacro.from_hdf5(macro_hdf5, "/bin/macro.exe")
        assert fm.index is None

    def test_raises_without_macro_params(self, tmp_path):
        """Should raise ValueError if the HDF5 file has no macro_params."""
        h5_path = tmp_path / "micro-only.h5"
        mp = MicroParameters()
        with h5py.File(str(h5_path), "w") as f:
            f.attrs[CONST.DATASPEC_VERSION_ATTR] = "v2.0.0"
            micro_grp = f.require_group("micro_data")
            for k, v in mp.to_basedict().items():
                micro_grp.attrs[k] = str(v) if not isinstance(v, (int, float, bool)) else v
        with pytest.raises(ValueError, match="macro_params"):
            FortranMacro.from_hdf5(h5_path, "/bin/macro.exe")

    def test_accepts_str_path(self, macro_hdf5):
        fm = FortranMacro.from_hdf5(str(macro_hdf5), "/bin/macro.exe")
        assert isinstance(fm, FortranMacro)


# ---------------------------------------------------------------------------
# TestFortranMacroExecInWorkdir
# ---------------------------------------------------------------------------


class TestFortranMacroExecInWorkdir:
    """Tests for :meth:`FortranMacro.exec_in_workdir`."""

    @pytest.fixture
    def fortran_macro_mock_setup(self, tmp_run):
        """FortranMacro with _write_setup_files mocked out."""
        fm = FortranMacro(run=tmp_run, executable="/bin/macro.exe")
        return fm

    def test_creates_data_run_code_dir(self, fortran_macro_mock_setup, tmp_path):
        fm = fortran_macro_mock_setup
        with (
            patch.object(FortranMacro, "_write_setup_files"),
            patch("subprocess.run"),
        ):
            fm.exec_in_workdir(tmp_path)
        assert (tmp_path / "data" / fm.run.run_code).is_dir()

    def test_returns_data_dir_path(self, fortran_macro_mock_setup, tmp_path):
        fm = fortran_macro_mock_setup
        with (
            patch.object(FortranMacro, "_write_setup_files"),
            patch("subprocess.run"),
        ):
            result = fm.exec_in_workdir(tmp_path)
        expected = tmp_path / "data" / fm.run.run_code
        assert result == expected

    def test_calls_write_setup_files_when_no_index(self, fortran_macro_mock_setup, tmp_path):
        fm = fortran_macro_mock_setup
        with (
            patch.object(FortranMacro, "_write_setup_files") as mock_setup,
            patch("subprocess.run"),
        ):
            fm.exec_in_workdir(tmp_path)
        assert mock_setup.call_count == 1

    def test_skips_write_setup_files_when_index_set(self, tmp_run, tmp_path):
        fm = FortranMacro(run=tmp_run, executable="/bin/macro.exe", index=0)
        with (
            patch.object(FortranMacro, "_write_setup_files") as mock_setup,
            patch("subprocess.run"),
        ):
            fm.exec_in_workdir(tmp_path)
        assert mock_setup.call_count == 0

    def test_subprocess_called_n_times_for_all_sims(self, fortran_macro_mock_setup, tmp_path):
        fm = fortran_macro_mock_setup
        n_sims = fm.run.macro_params.macro_simulations
        with (
            patch.object(FortranMacro, "_write_setup_files"),
            patch("subprocess.run") as mock_run,
        ):
            fm.exec_in_workdir(tmp_path)
        assert mock_run.call_count == n_sims

    def test_subprocess_called_once_when_index_set(self, tmp_run, tmp_path):
        fm = FortranMacro(run=tmp_run, executable="/bin/macro.exe", index=2)
        with (
            patch.object(FortranMacro, "_write_setup_files"),
            patch("subprocess.run") as mock_run,
        ):
            fm.exec_in_workdir(tmp_path)
        assert mock_run.call_count == 1

    def test_per_sim_subdirs_created(self, fortran_macro_mock_setup, tmp_path):
        fm = fortran_macro_mock_setup
        n_sims = fm.run.macro_params.macro_simulations
        with (
            patch.object(FortranMacro, "_write_setup_files"),
            patch("subprocess.run"),
        ):
            data_dir = fm.exec_in_workdir(tmp_path)
        for sim in range(n_sims):
            assert (data_dir / f"{sim:02}").is_dir()

    def test_log_file_moved_to_sim_subdir(self, fortran_macro_mock_setup, tmp_path):
        """Log file for sim 0 must land in data_dir/00/."""
        fm = fortran_macro_mock_setup
        with (
            patch.object(FortranMacro, "_write_setup_files"),
            patch("subprocess.run"),
        ):
            data_dir = fm.exec_in_workdir(tmp_path)
        # Log file for sim 0
        log_name = f"macro{fm.out_file_code}_00.txt"
        assert (data_dir / "00" / log_name).exists()

    def test_out_file_code_per_sim_suffix(self, tmp_run, tmp_path):
        """Binary must be called with --outFileCode {out_code}_{sim:02}."""
        fm = FortranMacro(run=tmp_run, executable="/bin/macro.exe",
                          out_file_code="_test", index=0)
        with (
            patch.object(FortranMacro, "_write_setup_files"),
            patch("subprocess.run") as mock_run,
        ):
            fm.exec_in_workdir(tmp_path)
        command = mock_run.call_args.args[0]
        idx = command.index("--outFileCode")
        assert command[idx + 1] == "_test_00"

    def test_accepts_str_work_dir(self, fortran_macro_mock_setup, tmp_path):
        fm = fortran_macro_mock_setup
        with (
            patch.object(FortranMacro, "_write_setup_files"),
            patch("subprocess.run"),
        ):
            result = fm.exec_in_workdir(str(tmp_path))
        assert isinstance(result, Path)


# ---------------------------------------------------------------------------
# TestFortranMacroImportResults
# ---------------------------------------------------------------------------


class TestFortranMacroImportResults:
    """Tests for :meth:`FortranMacro.import_results` (inherited from FortranRunner)."""

    @pytest.fixture
    def mock_datastore(self):
        """Return (mock_cls, mock_ds): a mocked DataStore class + open instance."""
        mock_ds = MagicMock()
        mock_cm = MagicMock()
        mock_cm.__enter__ = MagicMock(return_value=mock_ds)
        mock_cm.__exit__ = MagicMock(return_value=False)
        mock_cls = MagicMock(return_value=mock_cm)
        return mock_cls, mock_ds

    def test_delegates_to_import_collection_macroscale_out(
        self, fortran_macro, tmp_path, mock_datastore
    ):
        """import_results must call import_collection with 'macroscale_out'."""
        mock_cls, mock_ds = mock_datastore
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        with patch("lysis.execution.fortran.DataStore", mock_cls):
            fortran_macro.import_results(data_dir, keep_tmpdir=True)
        args = mock_ds.import_collection.call_args.args
        assert args[0] == "macroscale_out"

    def test_uses_macro_fortran_dataspec_version(
        self, fortran_macro, tmp_path, mock_datastore
    ):
        """import_results must pass MACRO_FORTRAN_DATASPEC_VERSION."""
        mock_cls, mock_ds = mock_datastore
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        with patch("lysis.execution.fortran.DataStore", mock_cls):
            fortran_macro.import_results(data_dir, keep_tmpdir=True)
        args = mock_ds.import_collection.call_args.args
        assert args[1] == MACRO_FORTRAN_DATASPEC_VERSION

    def test_passes_out_file_code(self, tmp_run, tmp_path, mock_datastore):
        """self.out_file_code must be forwarded to import_collection."""
        mock_cls, mock_ds = mock_datastore
        fm = FortranMacro(run=tmp_run, executable="/bin/macro.exe",
                          out_file_code="_code")
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        with patch("lysis.execution.fortran.DataStore", mock_cls):
            fm.import_results(data_dir, keep_tmpdir=True)
        assert mock_ds.import_collection.call_args.args[3] == ["_code"]

    def test_cleanup_on_success(self, fortran_macro, tmp_path, mock_datastore):
        """data_dir must be removed after successful import."""
        mock_cls, _ = mock_datastore
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        with patch("lysis.execution.fortran.DataStore", mock_cls):
            fortran_macro.import_results(data_dir, keep_tmpdir=False)
        assert not data_dir.exists()

    def test_keep_tmpdir_preserves_on_success(self, fortran_macro, tmp_path, mock_datastore):
        """data_dir must be preserved when keep_tmpdir=True."""
        mock_cls, _ = mock_datastore
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        with patch("lysis.execution.fortran.DataStore", mock_cls):
            fortran_macro.import_results(data_dir, keep_tmpdir=True)
        assert data_dir.exists()


# ---------------------------------------------------------------------------
# TestFortranMacroRunFull
# ---------------------------------------------------------------------------


class TestFortranMacroRunFull:
    """Tests for :meth:`FortranMacro.run_full` (inherited from FortranRunner)."""

    def test_exec_in_workdir_called(self, fortran_macro, tmp_path):
        with (
            patch.object(FortranMacro, "exec_in_workdir",
                         return_value=tmp_path / "data_dir") as mock_exec,
            patch.object(FortranMacro, "import_results"),
        ):
            fortran_macro.run_full()
        assert mock_exec.call_count == 1

    def test_import_results_called_with_data_dir(self, fortran_macro, tmp_path):
        data_dir = tmp_path / "data_dir"
        with (
            patch.object(FortranMacro, "exec_in_workdir", return_value=data_dir),
            patch.object(FortranMacro, "import_results") as mock_import,
        ):
            fortran_macro.run_full()
        assert mock_import.call_args.args[0] == data_dir

    def test_tmpdir_in_run_os_path(self, fortran_macro, tmp_path):
        """Temporary directory must be created inside run.os_path."""
        created_dirs = []

        def capture_exec(work_dir):
            created_dirs.append(work_dir)
            return work_dir / "data"

        with (
            patch.object(FortranMacro, "exec_in_workdir",
                         side_effect=capture_exec),
            patch.object(FortranMacro, "import_results"),
        ):
            fortran_macro.run_full()

        assert len(created_dirs) == 1
        tmpdir = created_dirs[0]
        assert str(tmpdir).startswith(fortran_macro.run.os_path)
        assert "lysis-macro" in str(tmpdir)

    def test_tmpdir_cleaned_on_success(self, fortran_macro, tmp_path):
        captured = []

        def capture_exec(work_dir):
            captured.append(work_dir)
            return work_dir / "data"

        with (
            patch.object(FortranMacro, "exec_in_workdir",
                         side_effect=capture_exec),
            patch.object(FortranMacro, "import_results"),
        ):
            fortran_macro.run_full()

        assert not captured[0].exists()

    def test_tmpdir_preserved_with_keep_tmpdir(self, fortran_macro, tmp_path):
        captured = []

        def capture_exec(work_dir):
            captured.append(work_dir)
            return work_dir / "data"

        with (
            patch.object(FortranMacro, "exec_in_workdir",
                         side_effect=capture_exec),
            patch.object(FortranMacro, "import_results"),
        ):
            fortran_macro.run_full(keep_tmpdir=True)

        assert captured[0].exists()
        captured[0].rmdir()  # cleanup
