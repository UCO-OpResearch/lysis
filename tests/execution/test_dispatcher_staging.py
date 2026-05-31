"""Unit tests for ``stage_dispatcher_logs`` on FortranMicro / FortranMacro (#55).

Synthesises Slurm worker ``.out`` files (no Slurm / Fortran involvement) and
checks that they are staged into ``data_dir`` under the dispatcher-log spec
filenames that :meth:`import_results` then ingests.
"""

from pathlib import Path

import pytest

from lysis.config.run import Run
from lysis.dataio.dataspec import dataspec
from lysis.execution.fortran_micro import FortranMicro
from lysis.execution.fortran_macro import FortranMacro


@pytest.fixture
def micro(tmp_path):
    r = Run(str(tmp_path))
    r.initialize_micro_param()
    return FortranMicro(run=r, executable="/bin/micro.exe", out_file_code="_TST")


@pytest.fixture
def macro(tmp_path):
    r = Run(str(tmp_path))
    r.initialize_micro_param()
    r.initialize_macro_param()
    return FortranMacro(
        run=r,
        executable="/bin/macro.exe",
        out_file_code="_TST",
        skip_binary_verification=True,
    )


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


class TestMicroStageDispatcherLogs:
    def test_concatenates_array_and_child_in_task_order(self, micro, tmp_path):
        src = tmp_path / "staging"
        data = tmp_path / "data"
        src.mkdir()
        data.mkdir()
        rc = micro.run.run_code
        # Out-of-order, with a two-digit index to catch lexical-vs-int sorting.
        _write(src / f"lysis-micro-array__{rc}__0.out", "A0\n")
        _write(src / f"lysis-micro-array__{rc}__10.out", "A10\n")
        _write(src / f"lysis-micro-array__{rc}__2.out", "A2\n")
        _write(src / f"lysis-micro-child__{rc}.out", "CHILD\n")
        micro.stage_dispatcher_logs(src, data)
        combined = (data / f"micro_dispatcher{micro.out_file_code}.out").read_text()
        assert combined == "A0\nA2\nA10\nCHILD\n"

    def test_no_worker_logs_stages_nothing(self, micro, tmp_path):
        src = tmp_path / "staging"
        data = tmp_path / "data"
        src.mkdir()
        data.mkdir()
        micro.stage_dispatcher_logs(src, data)
        assert list(data.iterdir()) == []

    def test_staged_filename_matches_spec(self, micro, tmp_path):
        # The bridge output must match the v1.99.0 source spec location so
        # import_results' read_data_collection picks it up.
        src = tmp_path / "staging"
        data = tmp_path / "data"
        src.mkdir()
        data.mkdir()
        _write(src / f"lysis-micro-child__{micro.run.run_code}.out", "X\n")
        micro.stage_dispatcher_logs(src, data)
        spec = dataspec["v1.99.0"]["microscale_out"].data["micro_dispatcher_log"]
        expected = spec.data_location.format(file_code=micro.out_file_code)
        assert (data / expected).exists()

    def test_master_log_is_not_captured(self, micro, tmp_path):
        src = tmp_path / "staging"
        data = tmp_path / "data"
        src.mkdir()
        data.mkdir()
        rc = micro.run.run_code
        _write(src / f"lysis-micro-master__{rc}.out", "MASTER\n")
        micro.stage_dispatcher_logs(src, data)
        assert list(data.iterdir()) == []


class TestMacroStageDispatcherLogs:
    def test_dense_per_sim_with_placeholder_for_missing_task(self, macro, tmp_path):
        src = tmp_path / "staging"
        data = tmp_path / "data"
        src.mkdir()
        data.mkdir()
        rc = macro.run.run_code
        for s in (0, 1, 2):
            (data / f"{s:02}").mkdir()
        _write(src / f"lysis-macro-array__{rc}__0.out", "S0\n")
        _write(src / f"lysis-macro-array__{rc}__1.out", "S1\n")
        # sim 2's task produced no .out file
        macro.stage_dispatcher_logs(src, data)
        ofc = macro.out_file_code
        assert (data / f"macro_dispatcher{ofc}_00.out").read_text() == "S0\n"
        assert (data / f"macro_dispatcher{ofc}_01.out").read_text() == "S1\n"
        # placeholder keeps the per-sim count dense (== macro_simulations)
        assert (data / f"macro_dispatcher{ofc}_02.out").read_text() == "\n"

    def test_no_worker_logs_stages_nothing(self, macro, tmp_path):
        src = tmp_path / "staging"
        data = tmp_path / "data"
        src.mkdir()
        data.mkdir()
        (data / "00").mkdir()
        macro.stage_dispatcher_logs(src, data)
        assert not list(data.glob("macro_dispatcher*"))

    def test_staged_filename_matches_spec(self, macro, tmp_path):
        src = tmp_path / "staging"
        data = tmp_path / "data"
        src.mkdir()
        data.mkdir()
        (data / "00").mkdir()
        _write(src / f"lysis-macro-array__{macro.run.run_code}__0.out", "X\n")
        macro.stage_dispatcher_logs(src, data)
        spec = dataspec["v1.99.0"]["macroscale_out"].data["macro_dispatcher_log"]
        expected = spec.data_location.format(file_code=macro.out_file_code, sim=0)
        assert (data / expected).exists()
