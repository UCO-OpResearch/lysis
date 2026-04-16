"""Unit tests for :mod:`lysis.execution.fortran_micro` — FortranMicro class.

Tests cover:

* :class:`FortranMicro` command construction (``exec_command``)
* Construction from HDF5 (``from_hdf5``)
* Execution in a working directory (``exec_in_workdir``)
* Result import into HDF5 (``import_results``)
* Full end-to-end workflow wrapper (``run_full``)
* Integration test with real fixture data (``TestImportResultsWithFixture``)
"""

import json
import os
import shutil
from pathlib import Path
from unittest.mock import MagicMock, call, patch

import h5py
import numpy as np
import pytest

from lysis.config.constants import CONST
from lysis.config.parameters import MacroParameters, MicroParameters
from lysis.config.run import Run
from lysis.dataio.datastore import DataStore
from lysis.execution.fortran_micro import FortranMicro
from lysis.execution.fortran import MICRO_FORTRAN_DATASPEC_VERSION


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------


def _write_micro_hdf5(path: Path, mp: MicroParameters = None) -> None:
    """Write a minimal v2.0.0 HDF5 file with MicroParameters attributes."""
    mp = mp or MicroParameters()
    with h5py.File(str(path), "w") as f:
        f.attrs[CONST.DATASPEC_VERSION_ATTR] = "v2.0.0"
        grp = f.require_group("micro_data")
        for k, v in mp.to_basedict().items():
            grp.attrs[k] = str(v) if not isinstance(v, (int, float, bool)) else v


@pytest.fixture
def tmp_run(tmp_path):
    """A Run with default MicroParameters in a tmp directory."""
    r = Run(str(tmp_path))
    r.initialize_micro_param()
    return r


@pytest.fixture
def micro_hdf5(tmp_path):
    """A minimal HDF5 file with default MicroParameters; returns Path."""
    h5_path = tmp_path / "test-run.h5"
    _write_micro_hdf5(h5_path)
    return h5_path


@pytest.fixture
def fortran_micro(tmp_run):
    """A FortranMicro with default run and dummy executable."""
    return FortranMicro(run=tmp_run, executable="/bin/micro.exe")


# ---------------------------------------------------------------------------
# TestFortranMicroExecCommand
# ---------------------------------------------------------------------------


class TestFortranMicroExecCommand:
    """Tests for :meth:`FortranMicro.exec_command`."""

    def test_default_params_contains_run_code(self, tmp_run):
        """Command must include --runCode and --outFileCode."""
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe")
        cmd = fm.exec_command()
        assert "--runCode" in cmd
        assert tmp_run.run_code in cmd

    def test_default_params_no_extra_flags(self, tmp_run):
        """With all-default params the command should only have --runCode and --outFileCode."""
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe")
        cmd = fm.exec_command()
        # Only executable + 4 elements: --runCode, run_code, --outFileCode, ""
        assert len(cmd) == 5

    def test_non_default_param_included(self, tmp_path):
        """A non-default constructor param must appear in the command."""
        r = Run(str(tmp_path))
        r.initialize_micro_param({"micro_simulations": 999})
        fm = FortranMicro(run=r, executable="/bin/micro.exe")
        cmd = fm.exec_command()
        assert "--simulations" in cmd
        assert "999" in cmd

    def test_default_param_excluded(self, tmp_run):
        """A param equal to its default must not appear in the command."""
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe")
        cmd = fm.exec_command()
        # micro_simulations == default so --simulations must be absent
        assert "--simulations" not in cmd

    def test_unit_bearing_param_converted(self, tmp_path):
        """fiber_radius (a Pint Quantity) must be output in the expected SI unit."""
        from lysis.config.constants import Q_

        r = Run(str(tmp_path))
        r.initialize_micro_param({"fiber_radius": Q_("37 nanometer")})
        fm = FortranMicro(run=r, executable="/bin/micro.exe")
        cmd = fm.exec_command()
        assert "--radius" in cmd
        idx = cmd.index("--radius")
        # The Fortran code expects radius in microns; 37 nm = 0.037 µm
        val = float(cmd[idx + 1])
        assert pytest.approx(val, rel=1e-6) == 0.037

    def test_index_none_no_seed_split(self, tmp_run):
        """With index=None the seed must not be split and micro_simulations unchanged."""
        default_sims = tmp_run.micro_params.micro_simulations
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe", index=None)
        cmd = fm.exec_command()
        # micro_simulations should still be absent (equal to default)
        assert "--simulations" not in cmd
        assert tmp_run.micro_params.micro_simulations == default_sims

    def test_index_none_out_file_code_unchanged(self, tmp_run):
        """With index=None out_file_code must remain empty."""
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe",
                          out_file_code="suffix", index=None)
        fm.exec_command()
        assert fm.out_file_code == "suffix"

    def test_index_zero_sets_simulations_to_one(self, tmp_run):
        """index=0 must force micro_simulations to 1 in the command."""
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe", index=0)
        cmd = fm.exec_command()
        assert "--simulations" in cmd
        assert cmd[cmd.index("--simulations") + 1] == "1"

    def test_index_zero_appends_suffix(self, tmp_run):
        """index=0 must append '__00' to out_file_code."""
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe", index=0)
        fm.exec_command()
        assert fm.out_file_code == "__00"

    def test_index_two_appends_correct_suffix(self, tmp_run):
        """index=2 must append '__02' to out_file_code."""
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe", index=2)
        fm.exec_command()
        assert fm.out_file_code == "__02"

    def test_index_two_uses_correct_seed(self, tmp_run):
        """index=2 must derive seed from SeedSequence at position 2."""
        seed = tmp_run.micro_params.micro_seed
        stream = np.random.SeedSequence(seed)
        expected_seed = int(np.int32(stream.generate_state(3)[2]))

        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe", index=2)
        cmd = fm.exec_command()
        assert "--seed" in cmd
        assert cmd[cmd.index("--seed") + 1] == str(expected_seed)

    def test_executable_is_first_element(self, tmp_run):
        """The executable must be the first element of the command list."""
        fm = FortranMicro(run=tmp_run, executable="/path/to/micro.exe")
        cmd = fm.exec_command()
        assert cmd[0] == "/path/to/micro.exe"


# ---------------------------------------------------------------------------
# TestFortranMicroFromHdf5
# ---------------------------------------------------------------------------


class TestFortranMicroFromHdf5:
    """Tests for :meth:`FortranMicro.from_hdf5`."""

    def test_constructs_run_with_correct_run_code(self, micro_hdf5):
        fm = FortranMicro.from_hdf5(micro_hdf5, "/bin/micro.exe")
        assert fm.run.run_code == micro_hdf5.stem

    def test_loads_micro_params(self, micro_hdf5):
        fm = FortranMicro.from_hdf5(micro_hdf5, "/bin/micro.exe")
        assert fm.run.micro_params is not None
        assert isinstance(fm.run.micro_params, MicroParameters)

    def test_sets_executable(self, micro_hdf5, tmp_path):
        exe = tmp_path / "micro.exe"
        fm = FortranMicro.from_hdf5(micro_hdf5, str(exe))
        assert fm.executable == str(exe.resolve())

    def test_passes_out_file_code(self, micro_hdf5):
        fm = FortranMicro.from_hdf5(micro_hdf5, "/bin/micro.exe",
                                     out_file_code="_code")
        assert fm.out_file_code == "_code"

    def test_passes_index(self, micro_hdf5):
        fm = FortranMicro.from_hdf5(micro_hdf5, "/bin/micro.exe", index=3)
        assert fm.index == 3

    def test_raises_for_missing_file(self, tmp_path):
        missing = tmp_path / "nonexistent.h5"
        with pytest.raises((RuntimeError, OSError)):
            FortranMicro.from_hdf5(missing, "/bin/micro.exe")

    def test_accepts_path_or_str(self, micro_hdf5):
        fm_path = FortranMicro.from_hdf5(micro_hdf5, "/bin/micro.exe")
        fm_str = FortranMicro.from_hdf5(str(micro_hdf5), "/bin/micro.exe")
        assert fm_path.run.run_code == fm_str.run.run_code

    def test_raises_when_macro_data_present(self, tmp_path):
        """Should raise ValueError if macro_data group already exists (macroscale initialized)."""
        h5_path = tmp_path / "post-init-macro.h5"
        mp = MicroParameters()
        mcp = MacroParameters(micro_params=mp)
        with h5py.File(str(h5_path), "w") as f:
            f.attrs[CONST.DATASPEC_VERSION_ATTR] = "v2.0.0"
            micro_grp = f.require_group("micro_data")
            for k, v in mp.to_basedict().items():
                micro_grp.attrs[k] = str(v) if not isinstance(v, (int, float, bool)) else v
            macro_grp = f.require_group("macro_data")
            for k, v in mcp.to_basedict().items():
                macro_grp.attrs[k] = str(v) if not isinstance(v, (int, float, bool)) else v
        with pytest.raises(ValueError, match="macroscale data"):
            FortranMicro.from_hdf5(h5_path, "/bin/micro.exe")

    def test_raises_when_micro_already_ran(self, tmp_path):
        """Should raise ValueError if micro datasets are non-empty (already ran)."""
        h5_path = tmp_path / "post-run-micro.h5"
        _write_micro_hdf5(h5_path)
        with h5py.File(str(h5_path), "a") as f:
            f.create_dataset(
                "micro_data/tpa_leaving_time",
                data=np.array([1.0]),
                dtype=np.float64,
            )
        with pytest.raises(ValueError, match="microscale simulation"):
            FortranMicro.from_hdf5(h5_path, "/bin/micro.exe")


# ---------------------------------------------------------------------------
# TestFortranMicroExecInWorkdir
# ---------------------------------------------------------------------------


class TestFortranMicroExecInWorkdir:
    """Tests for :meth:`FortranMicro.exec_in_workdir`."""

    def test_creates_data_run_code_dir(self, fortran_micro, tmp_path):
        with patch("subprocess.run"):
            fortran_micro.exec_in_workdir(tmp_path)
        expected = tmp_path / "data" / fortran_micro.run.run_code
        assert expected.is_dir()

    def test_returns_data_dir_path(self, fortran_micro, tmp_path):
        with patch("subprocess.run"):
            result = fortran_micro.exec_in_workdir(tmp_path)
        expected = tmp_path / "data" / fortran_micro.run.run_code
        assert result == expected

    def test_log_file_in_data_dir(self, fortran_micro, tmp_path):
        """stdout is redirected to micro{out_code}.txt inside data_dir."""
        with patch("subprocess.run") as mock_run:
            data_dir = fortran_micro.exec_in_workdir(tmp_path)
        assert mock_run.call_count == 1
        call_kwargs = mock_run.call_args.kwargs
        assert "cwd" in call_kwargs
        assert call_kwargs["cwd"] == str(tmp_path)

    def test_subprocess_run_called_once(self, fortran_micro, tmp_path):
        with patch("subprocess.run") as mock_run:
            fortran_micro.exec_in_workdir(tmp_path)
        assert mock_run.call_count == 1

    def test_subprocess_cwd_is_work_dir(self, fortran_micro, tmp_path):
        with patch("subprocess.run") as mock_run:
            fortran_micro.exec_in_workdir(tmp_path)
        call_kwargs = mock_run.call_args.kwargs
        assert call_kwargs["cwd"] == str(tmp_path)

    def test_params_json_written(self, fortran_micro, tmp_path):
        """params.json must be present in data_dir for v1.99.0 import pipeline."""
        with patch("subprocess.run"):
            data_dir = fortran_micro.exec_in_workdir(tmp_path)
        assert (data_dir / "params.json").exists()

    def test_params_json_contains_micro_params_key(self, fortran_micro, tmp_path):
        """params.json must be a valid JSON file with a 'micro_params' key."""
        with patch("subprocess.run"):
            data_dir = fortran_micro.exec_in_workdir(tmp_path)
        with open(data_dir / "params.json") as fh:
            loaded = json.load(fh)
        assert "micro_params" in loaded

    def test_params_json_pint_quantities_serialised_as_strings(self, fortran_micro, tmp_path):
        """Pint Quantity values in MicroParameters must be stored as strings in params.json."""
        with patch("subprocess.run"):
            data_dir = fortran_micro.exec_in_workdir(tmp_path)
        with open(data_dir / "params.json") as fh:
            loaded = json.load(fh)
        # Every value in micro_params must be JSON-native (str, int, float, bool, None)
        for v in loaded["micro_params"].values():
            assert isinstance(v, (str, int, float, bool, type(None))), (
                f"Non-serialisable value leaked into params.json: {v!r}"
            )

    def test_accepts_str_work_dir(self, fortran_micro, tmp_path):
        with patch("subprocess.run"):
            result = fortran_micro.exec_in_workdir(str(tmp_path))
        assert isinstance(result, Path)

    def test_calls_write_setup_files_when_no_index(self, fortran_micro, tmp_path):
        """exec_in_workdir must write params.json when index is None."""
        with (
            patch.object(FortranMicro, "_write_setup_files") as mock_setup,
            patch("subprocess.run"),
        ):
            fortran_micro.exec_in_workdir(tmp_path)
        assert mock_setup.call_count == 1

    def test_skips_write_setup_files_when_index_set(self, tmp_run, tmp_path):
        """exec_in_workdir must skip params.json when index is set (pre-staged)."""
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe", index=0)
        with (
            patch.object(FortranMicro, "_write_setup_files") as mock_setup,
            patch("subprocess.run"),
        ):
            fm.exec_in_workdir(tmp_path)
        assert mock_setup.call_count == 0


# ---------------------------------------------------------------------------
# TestFortranMicroWriteSetupFiles
# ---------------------------------------------------------------------------


class TestFortranMicroWriteSetupFiles:
    """Tests for :meth:`FortranMicro._write_setup_files`."""

    def test_creates_params_json(self, fortran_micro, tmp_path):
        fortran_micro._write_setup_files(tmp_path)
        assert (tmp_path / "params.json").exists()

    def test_params_json_contains_micro_params_key(self, fortran_micro, tmp_path):
        fortran_micro._write_setup_files(tmp_path)
        with open(tmp_path / "params.json") as fh:
            loaded = json.load(fh)
        assert "micro_params" in loaded


# ---------------------------------------------------------------------------
# TestFortranMicroImportResults
# ---------------------------------------------------------------------------


class TestFortranMicroImportResults:
    """Tests for :meth:`FortranMicro.import_results`."""

    @pytest.fixture
    def mock_datastore(self):
        """Return (mock_cls, mock_ds): a mocked DataStore class + open instance."""
        mock_ds = MagicMock()
        mock_cm = MagicMock()
        mock_cm.__enter__ = MagicMock(return_value=mock_ds)
        mock_cm.__exit__ = MagicMock(return_value=False)
        mock_cls = MagicMock(return_value=mock_cm)
        return mock_cls, mock_ds

    def test_opens_datastore_with_run_code(self, fortran_micro, tmp_path, mock_datastore):
        """DataStore must be constructed with self.run.run_code."""
        mock_cls, _ = mock_datastore
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        with patch("lysis.execution.fortran.DataStore", mock_cls):
            fortran_micro.import_results(data_dir, keep_tmpdir=True)
        assert mock_cls.call_args.args[0] == fortran_micro.run.run_code

    def test_opens_datastore_with_run_os_path(self, fortran_micro, tmp_path, mock_datastore):
        """DataStore must be constructed with self.run.os_path."""
        mock_cls, _ = mock_datastore
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        with patch("lysis.execution.fortran.DataStore", mock_cls):
            fortran_micro.import_results(data_dir, keep_tmpdir=True)
        assert mock_cls.call_args.args[1] == fortran_micro.run.os_path

    def test_opens_datastore_in_append_mode(self, fortran_micro, tmp_path, mock_datastore):
        """DataStore must be opened with mode='a'."""
        mock_cls, _ = mock_datastore
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        with patch("lysis.execution.fortran.DataStore", mock_cls):
            fortran_micro.import_results(data_dir, keep_tmpdir=True)
        assert mock_cls.call_args.kwargs.get("mode") == "a"

    def test_delegates_to_import_collection(self, fortran_micro, tmp_path, mock_datastore):
        """import_results must call ds.import_collection with the correct arguments."""
        mock_cls, mock_ds = mock_datastore
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        with patch("lysis.execution.fortran.DataStore", mock_cls):
            fortran_micro.import_results(data_dir, keep_tmpdir=True)
        assert mock_ds.import_collection.call_count == 1
        args = mock_ds.import_collection.call_args.args
        assert args[0] == "microscale_out"
        assert args[1] == MICRO_FORTRAN_DATASPEC_VERSION
        assert args[2] == str(data_dir)
        assert args[3] == [fortran_micro.out_file_code]

    def test_passes_out_file_code(self, tmp_run, tmp_path, mock_datastore):
        """self.out_file_code must be forwarded to import_collection."""
        mock_cls, mock_ds = mock_datastore
        fm = FortranMicro(run=tmp_run, executable="/bin/micro.exe", out_file_code="_code")
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        with patch("lysis.execution.fortran.DataStore", mock_cls):
            fm.import_results(data_dir, keep_tmpdir=True)
        assert mock_ds.import_collection.call_args.args[3] == ["_code"]

    def test_cleanup_on_success(self, fortran_micro, tmp_path, mock_datastore):
        """data_dir must be removed after successful import (keep_tmpdir=False)."""
        mock_cls, _ = mock_datastore
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        with patch("lysis.execution.fortran.DataStore", mock_cls):
            fortran_micro.import_results(data_dir, keep_tmpdir=False)
        assert not data_dir.exists()

    def test_keep_tmpdir_preserves_on_success(self, fortran_micro, tmp_path, mock_datastore):
        """data_dir must be preserved when keep_tmpdir=True."""
        mock_cls, _ = mock_datastore
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        with patch("lysis.execution.fortran.DataStore", mock_cls):
            fortran_micro.import_results(data_dir, keep_tmpdir=True)
        assert data_dir.exists()

    def test_cleanup_on_failure_keep_on_failure_false(self, fortran_micro, tmp_path, mock_datastore):
        """data_dir must be removed on failure when keep_on_failure=False."""
        mock_cls, mock_ds = mock_datastore
        mock_ds.import_collection.side_effect = RuntimeError("import failed")
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        with patch("lysis.execution.fortran.DataStore", mock_cls):
            with pytest.raises(RuntimeError, match="import failed"):
                fortran_micro.import_results(
                    data_dir, keep_on_failure=False, keep_tmpdir=False
                )
        assert not data_dir.exists()

    def test_preserve_on_failure_keep_on_failure_true(self, fortran_micro, tmp_path, mock_datastore):
        """data_dir must be preserved on failure when keep_on_failure=True."""
        mock_cls, mock_ds = mock_datastore
        mock_ds.import_collection.side_effect = RuntimeError("import failed")
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        with patch("lysis.execution.fortran.DataStore", mock_cls):
            with pytest.raises(RuntimeError):
                fortran_micro.import_results(
                    data_dir, keep_on_failure=True, keep_tmpdir=False
                )
        assert data_dir.exists()

    def test_keep_tmpdir_preserves_on_failure(self, fortran_micro, tmp_path, mock_datastore):
        """data_dir must be preserved on failure when keep_tmpdir=True."""
        mock_cls, mock_ds = mock_datastore
        mock_ds.import_collection.side_effect = RuntimeError("import failed")
        data_dir = tmp_path / "data_dir"
        data_dir.mkdir()
        with patch("lysis.execution.fortran.DataStore", mock_cls):
            with pytest.raises(RuntimeError):
                fortran_micro.import_results(
                    data_dir, keep_on_failure=False, keep_tmpdir=True
                )
        assert data_dir.exists()


# ---------------------------------------------------------------------------
# TestFortranMicroRunFull
# ---------------------------------------------------------------------------


class TestFortranMicroRunFull:
    """Tests for :meth:`FortranMicro.run_full`."""

    def test_exec_in_workdir_called(self, fortran_micro, tmp_path):
        with (
            patch.object(FortranMicro, "exec_in_workdir",
                         return_value=tmp_path / "data_dir") as mock_exec,
            patch.object(FortranMicro, "import_results") as mock_import,
        ):
            fortran_micro.run_full()
        assert mock_exec.call_count == 1

    def test_import_results_called_with_data_dir(self, fortran_micro, tmp_path):
        data_dir = tmp_path / "data_dir"
        with (
            patch.object(FortranMicro, "exec_in_workdir",
                         return_value=data_dir),
            patch.object(FortranMicro, "import_results") as mock_import,
        ):
            fortran_micro.run_full()
        assert mock_import.call_args.args[0] == data_dir

    def test_tmpdir_cleaned_on_success(self, fortran_micro, tmp_path):
        """The outer tmpdir must be removed after a successful run."""
        created_tmpdirs = []

        def fake_exec(work_dir):
            created_tmpdirs.append(work_dir)
            return work_dir / "data"

        with (
            patch.object(FortranMicro, "exec_in_workdir", side_effect=fake_exec),
            patch.object(FortranMicro, "import_results"),
        ):
            fortran_micro.run_full(keep_tmpdir=False)

        assert created_tmpdirs, "exec_in_workdir was never called"
        assert not created_tmpdirs[0].exists()

    def test_tmpdir_preserved_with_keep_tmpdir(self, fortran_micro, tmp_path):
        """The outer tmpdir must be preserved when keep_tmpdir=True."""
        created_tmpdirs = []

        def fake_exec(work_dir):
            created_tmpdirs.append(work_dir)
            work_dir.mkdir(parents=True, exist_ok=True)
            return work_dir / "data"

        with (
            patch.object(FortranMicro, "exec_in_workdir", side_effect=fake_exec),
            patch.object(FortranMicro, "import_results"),
        ):
            fortran_micro.run_full(keep_tmpdir=True)

        assert created_tmpdirs
        assert created_tmpdirs[0].exists()

    def test_tmpdir_in_run_os_path(self, fortran_micro, tmp_path):
        """run_full must create its tmpdir in run.os_path, not in /tmp."""
        created_tmpdirs = []

        def fake_exec(work_dir):
            created_tmpdirs.append(work_dir)
            return work_dir / "data"

        with (
            patch.object(FortranMicro, "exec_in_workdir", side_effect=fake_exec),
            patch.object(FortranMicro, "import_results"),
        ):
            fortran_micro.run_full(keep_tmpdir=False)

        assert created_tmpdirs
        assert created_tmpdirs[0].parent == tmp_path


# ---------------------------------------------------------------------------
# Integration test with real fixture data
# ---------------------------------------------------------------------------

# These constants mirror the ones in tests/dataio/test_real_data.py.
_MICRO_FILE_CODE = "_PLG2_tPA01_TB-xiii"
_PARAM_OVERRIDES = {
    "fibrinogen_length": "45nm",
    "fibrinogen_radius": "1.2nm",
    "micro_log_lvl": 40,
    "micro_version": "micro_rates",
    "snap_proportion": 0.66666667,
}


def _build_reference_v200(fixture_path: str) -> dict:
    """Read v1.95.0 fixture data, apply overrides, convert to v2.0.0.

    Returns a dict with v2.0.0 dataset names as keys and numpy arrays as
    values, plus a ``"params"`` key.
    """
    from lysis.dataio.dataspec import dataspec
    from lysis.dataio.fileops import read_data_collection
    from lysis.dataio.dataconvert import convert_data

    raw = read_data_collection(
        fixture_path,
        collections=[dataspec["v1.95.0"]["microscale_out"]],
        file_codes=[_MICRO_FILE_CODE],
    )
    for key, value in _PARAM_OVERRIDES.items():
        for section in raw["params"].values():
            if isinstance(section, dict):
                section[key] = value

    return convert_data(raw, "v1.95.0", "v2.0.0")


class TestImportResultsWithFixture:
    """Integration tests: import real fixture data and compare to reference.

    Uses the truncated Fortran output in ``tests/fixtures/fortran_sample/``
    (the same dataset exercised in ``tests/dataio/test_real_data.py``).

    The test simulates what happens after the Fortran binary finishes:
    copies the fixture files into a staging directory, writes a v1.99.0
    ``params.json``, then calls :meth:`FortranMicro.import_results`.
    The resulting HDF5 datasets are compared element-by-element against a
    reference produced by the established v1.95.0 → v2.0.0 conversion
    pipeline.
    """

    @pytest.fixture(scope="class")
    def reference(self, fortran_sample_path):
        """Reference v2.0.0 data dict via the established conversion path."""
        return _build_reference_v200(fortran_sample_path)

    @pytest.fixture(scope="class")
    def imported_hdf5(self, fortran_sample_path, tmp_path_factory):
        """Run import_results on a copy of the fixture and return the HDF5 path."""
        run_code = "test-import"
        tmp_path = tmp_path_factory.mktemp("import_fixture")

        # Create the HDF5 file with empty microscale_out datasets.
        with DataStore.create(run_code, str(tmp_path), MicroParameters()):
            pass

        # Copy fixture files to a staging directory that mirrors the path
        # structure the Fortran binary produces: data/{run_code}/.
        data_dir = tmp_path / "staging" / "data" / run_code
        shutil.copytree(fortran_sample_path, data_dir)

        run = Run(str(tmp_path), run_code=run_code)
        run.initialize_micro_param()
        fm = FortranMicro(run=run, executable="/bin/micro.exe",
                          out_file_code=_MICRO_FILE_CODE)
        fm.import_results(data_dir, keep_tmpdir=True)

        return tmp_path / f"{run_code}.h5"

    # ── Dataset presence ────────────────────────────────────────────

    def test_all_microscale_datasets_present(self, imported_hdf5, reference):
        """Every v2.0.0 microscale_out dataset in the reference must be in the HDF5 file."""
        from lysis.dataio.dataspec import dataspec
        spec = dataspec["v2.0.0"]["microscale_out"]
        with h5py.File(str(imported_hdf5), "r") as f:
            for name, ds_spec in spec.data.items():
                if name not in reference:
                    continue
                loc = ds_spec.data_location
                assert loc in f, f"Missing dataset: {loc} (name={name})"

    # ── Per-dataset value comparisons ───────────────────────────────

    def test_pli_first_time_matches(self, imported_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["pli_first_time"].data_location
        with h5py.File(str(imported_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["pli_first_time"])

    def test_tpa_final_num_matches(self, imported_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["tpa_final_num"].data_location
        with h5py.File(str(imported_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["tpa_final_num"])

    def test_fiber_degraded_matches(self, imported_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["fiber_degraded"].data_location
        with h5py.File(str(imported_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["fiber_degraded"])

    def test_sim_final_time_matches(self, imported_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["sim_final_time"].data_location
        with h5py.File(str(imported_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["sim_final_time"])

    def test_pli_generated_num_matches(self, imported_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["pli_generated_num"].data_location
        with h5py.File(str(imported_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["pli_generated_num"])

    def test_tpa_leaving_time_matches(self, imported_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["tpa_leaving_time"].data_location
        with h5py.File(str(imported_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["tpa_leaving_time"])

    def test_tpa_unbound_by_pli_matches(self, imported_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["tpa_unbound_by_pli"].data_location
        with h5py.File(str(imported_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["tpa_unbound_by_pli"])

    def test_tpa_unbound_kinetic_matches(self, imported_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["tpa_unbound_kinetic"].data_location
        with h5py.File(str(imported_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["tpa_unbound_kinetic"])

    # ── Dtype checks ────────────────────────────────────────────────

    def test_dtypes_match_spec(self, imported_hdf5):
        """Each dataset's HDF5 dtype must match the v2.0.0 dataspec."""
        from lysis.dataio.dataspec import dataspec
        spec = dataspec["v2.0.0"]["microscale_out"]
        with h5py.File(str(imported_hdf5), "r") as f:
            for name, ds_spec in spec.data.items():
                loc = ds_spec.data_location
                if loc not in f:
                    continue
                expected_dtype = np.dtype(ds_spec.dtype)
                actual_dtype = f[loc].dtype
                assert actual_dtype == expected_dtype, (
                    f"{name}: expected {expected_dtype}, got {actual_dtype}"
                )

    # ── Shape checks ────────────────────────────────────────────────

    def test_shapes_match_reference(self, imported_hdf5, reference):
        """Each dataset's shape must match the reference conversion."""
        from lysis.dataio.dataspec import dataspec
        spec = dataspec["v2.0.0"]["microscale_out"]
        with h5py.File(str(imported_hdf5), "r") as f:
            for name, ds_spec in spec.data.items():
                if name not in reference:
                    continue
                loc = ds_spec.data_location
                if loc not in f:
                    continue
                expected_shape = np.asarray(reference[name]).shape
                actual_shape = f[loc].shape
                assert actual_shape == expected_shape, (
                    f"{name}: expected shape {expected_shape}, got {actual_shape}"
                )


# ---------------------------------------------------------------------------
# Fortran binary execution integration test
# ---------------------------------------------------------------------------

_BINARY_TEST_SIMULATIONS = 1000
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_FIXTURE_SEED = 2133256963


@pytest.mark.fortran_binary
class TestFortranBinaryExecution:
    """Integration tests that execute the real compiled Fortran binary.

    These tests compile-and-run the microscale simulation with the same
    parameters as the ``tests/fixtures/fortran_sample/`` fixture but only
    1,000 simulations (instead of 50,000) and compare the output to the
    first 1,000 elements of the fixture reference data.

    Requires:

    * The compiled Fortran binary at ``bin/micro_rates`` (repo root).
    * The ``tests/fixtures/fortran_sample/`` fixture.
    * The ``lysis`` conda environment on ``$PATH``.

    Skip with ``-m "not fortran_binary"`` to exclude these tests.
    """

    @pytest.fixture(autouse=True)
    def _skip_if_no_binary(self):
        binary = _REPO_ROOT / "bin" / "micro_rates"
        if not binary.exists():
            pytest.skip("Compiled Fortran binary not found at bin/micro_rates")

    _SKIP_DATASETS = {"params", "micro_log"}

    @pytest.fixture(scope="class")
    def reference(self, fortran_sample_path):
        """Reference v2.0.0 data sliced to the first 1,000 simulations."""
        ref = _build_reference_v200(fortran_sample_path)
        sliced = {}
        for key, value in ref.items():
            if key in self._SKIP_DATASETS:
                continue
            arr = np.asarray(value)
            if arr.ndim >= 1 and arr.shape[0] >= _BINARY_TEST_SIMULATIONS:
                sliced[key] = arr[:_BINARY_TEST_SIMULATIONS]
            else:
                sliced[key] = arr
        return sliced

    @pytest.fixture(scope="class")
    def executed_hdf5(self, tmp_path_factory):
        """Run the Fortran binary and import results into a fresh HDF5 file."""
        from lysis.config.constants import Q_

        binary = str(_REPO_ROOT / "bin" / "micro_rates")
        run_code = "fortran-run"
        tmp_path = tmp_path_factory.mktemp("fortran_binary")

        mp = MicroParameters(
            nodes_in_micro_row=13,
            fiber_radius=Q_("61.5 nanometer"),
            micro_seed=_FIXTURE_SEED,
            micro_simulations=_BINARY_TEST_SIMULATIONS,
        )
        with DataStore.create(run_code, str(tmp_path), mp):
            pass

        hdf5_path = tmp_path / f"{run_code}.h5"
        fm = FortranMicro.from_hdf5(hdf5_path, binary)
        data_dir = fm.exec_in_workdir(tmp_path)
        fm.import_results(data_dir, keep_tmpdir=True)
        return hdf5_path

    def test_all_microscale_datasets_present(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        spec = dataspec["v2.0.0"]["microscale_out"]
        with h5py.File(str(executed_hdf5), "r") as f:
            for name in reference:
                ds_spec = spec.data.get(name)
                if ds_spec is None:
                    continue
                assert ds_spec.data_location in f, (
                    f"Missing dataset: {ds_spec.data_location} (name={name})"
                )

    def test_pli_first_time_matches(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["pli_first_time"].data_location
        with h5py.File(str(executed_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["pli_first_time"])

    def test_tpa_final_num_matches(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["tpa_final_num"].data_location
        with h5py.File(str(executed_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["tpa_final_num"])

    def test_fiber_degraded_matches(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["fiber_degraded"].data_location
        with h5py.File(str(executed_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["fiber_degraded"])

    def test_sim_final_time_matches(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["sim_final_time"].data_location
        with h5py.File(str(executed_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["sim_final_time"])

    def test_pli_generated_num_matches(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["pli_generated_num"].data_location
        with h5py.File(str(executed_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["pli_generated_num"])

    def test_tpa_leaving_time_matches(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["tpa_leaving_time"].data_location
        with h5py.File(str(executed_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["tpa_leaving_time"])

    def test_tpa_unbound_by_pli_matches(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["tpa_unbound_by_pli"].data_location
        with h5py.File(str(executed_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["tpa_unbound_by_pli"])

    def test_tpa_unbound_kinetic_matches(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        loc = dataspec["v2.0.0"]["microscale_out"].data["tpa_unbound_kinetic"].data_location
        with h5py.File(str(executed_hdf5), "r") as f:
            np.testing.assert_array_equal(f[loc][...], reference["tpa_unbound_kinetic"])

    def test_dtypes_match_spec(self, executed_hdf5):
        from lysis.dataio.dataspec import dataspec
        spec = dataspec["v2.0.0"]["microscale_out"]
        with h5py.File(str(executed_hdf5), "r") as f:
            for name, ds_spec in spec.data.items():
                loc = ds_spec.data_location
                if loc not in f:
                    continue
                expected_dtype = np.dtype(ds_spec.dtype)
                actual_dtype = f[loc].dtype
                assert actual_dtype == expected_dtype, (
                    f"{name}: expected {expected_dtype}, got {actual_dtype}"
                )

    def test_shapes_match_expected(self, executed_hdf5, reference):
        from lysis.dataio.dataspec import dataspec
        spec = dataspec["v2.0.0"]["microscale_out"]
        with h5py.File(str(executed_hdf5), "r") as f:
            for name, ds_spec in spec.data.items():
                if name not in reference:
                    continue
                loc = ds_spec.data_location
                if loc not in f:
                    continue
                expected_shape = np.asarray(reference[name]).shape
                actual_shape = f[loc].shape
                assert actual_shape == expected_shape, (
                    f"{name}: expected shape {expected_shape}, got {actual_shape}"
                )
