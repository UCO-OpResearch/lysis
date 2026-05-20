"""Pure-unit tests for :meth:`FortranMicro.concatenate_child_outputs`.

Synthesizes per-task ``.dat`` files with known byte contents (no Fortran
binary involvement) and verifies that the concatenation is byte-correct,
that the pre-staged ``params.json`` is not disturbed, and that error
paths raise informative exceptions.
"""

import json
from pathlib import Path

import pytest

from lysis.config.parameters import MicroParameters
from lysis.config.run import Run
from lysis.execution.fortran_micro import (
    FortranMicro,
    MICROSCALE_BINARY_DATASETS,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def tmp_run(tmp_path):
    r = Run(str(tmp_path))
    r.initialize_micro_param()
    return r


@pytest.fixture
def fortran_micro(tmp_run):
    return FortranMicro(
        run=tmp_run, executable="/bin/micro.exe", out_file_code="_TST"
    )


def _write_per_task_files(
    data_dir: Path,
    out_file_code: str,
    num_children: int,
    payload_size: int = 16,
):
    """Write deterministic ``__NN.dat`` files for every dataset.

    Each task's payload is ``bytes([i] * payload_size)`` so the
    concatenated output is trivially predictable.
    """
    data_dir.mkdir(parents=True, exist_ok=True)
    expected_concat = {}
    for dataset in MICROSCALE_BINARY_DATASETS:
        chunks = []
        for i in range(num_children):
            payload = bytes([i % 256] * payload_size)
            (data_dir / f"{dataset}{out_file_code}__{i:02}.dat").write_bytes(
                payload
            )
            chunks.append(payload)
        expected_concat[dataset] = b"".join(chunks)
    return expected_concat


def _write_per_task_logs(
    data_dir: Path, out_file_code: str, num_children: int
):
    """Write per-task ``micro{code}__NN.txt`` log files; return concat bytes."""
    chunks = []
    for i in range(num_children):
        line = f"task {i} log\n".encode()
        (data_dir / f"micro{out_file_code}__{i:02}.txt").write_bytes(line)
        chunks.append(line)
    return b"".join(chunks)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestConcatenateChildOutputs:
    """:meth:`FortranMicro.concatenate_child_outputs` byte-correctness tests."""

    def test_concat_binary_equal_sizes(self, tmp_path, fortran_micro):
        """Three equal-sized chunks per dataset must concatenate in task order."""
        expected = _write_per_task_files(
            tmp_path, fortran_micro.out_file_code, num_children=3
        )
        fortran_micro.concatenate_child_outputs(tmp_path, num_children=3)

        for dataset in MICROSCALE_BINARY_DATASETS:
            combined = (
                tmp_path / f"{dataset}{fortran_micro.out_file_code}.dat"
            ).read_bytes()
            assert combined == expected[dataset], dataset

    def test_concat_all_eight_datasets(self, tmp_path, fortran_micro):
        """Every dataset in MICROSCALE_BINARY_DATASETS must yield an output file."""
        _write_per_task_files(
            tmp_path, fortran_micro.out_file_code, num_children=4, payload_size=8
        )
        fortran_micro.concatenate_child_outputs(tmp_path, num_children=4)

        for dataset in MICROSCALE_BINARY_DATASETS:
            combined = tmp_path / f"{dataset}{fortran_micro.out_file_code}.dat"
            assert combined.exists(), f"missing combined output for {dataset}"
            assert combined.stat().st_size == 4 * 8

    def test_concat_log_files(self, tmp_path, fortran_micro):
        """The per-task log files must concatenate into the unsuffixed log."""
        _write_per_task_files(
            tmp_path, fortran_micro.out_file_code, num_children=3
        )
        expected_log = _write_per_task_logs(
            tmp_path, fortran_micro.out_file_code, num_children=3
        )

        fortran_micro.concatenate_child_outputs(tmp_path, num_children=3)

        log = tmp_path / f"micro{fortran_micro.out_file_code}.txt"
        assert log.exists()
        assert log.read_bytes() == expected_log

    def test_concat_missing_file_raises(self, tmp_path, fortran_micro):
        """A missing per-task .dat file must raise FileNotFoundError naming the index."""
        _write_per_task_files(
            tmp_path, fortran_micro.out_file_code, num_children=4
        )
        # Drop task 02 of the first dataset.
        missing = (
            tmp_path
            / f"{MICROSCALE_BINARY_DATASETS[0]}{fortran_micro.out_file_code}__02.dat"
        )
        missing.unlink()

        with pytest.raises(FileNotFoundError) as excinfo:
            fortran_micro.concatenate_child_outputs(tmp_path, num_children=4)
        msg = str(excinfo.value)
        assert "02" in msg
        assert MICROSCALE_BINARY_DATASETS[0] in msg

    def test_concat_preserves_pre_staged_params_json(
        self, tmp_path, fortran_micro
    ):
        """Pre-staged params.json (with aggregate count) must not be touched."""
        _write_per_task_files(
            tmp_path, fortran_micro.out_file_code, num_children=2
        )
        params_path = tmp_path / "params.json"
        sentinel = {
            "micro_params": {"micro_simulations": 12345, "micro_seed": 42}
        }
        params_path.write_text(json.dumps(sentinel))

        fortran_micro.concatenate_child_outputs(tmp_path, num_children=2)

        assert json.loads(params_path.read_text()) == sentinel

    def test_concat_single_child_renames_only(self, tmp_path, fortran_micro):
        """num_children=1 must rename __00 → unsuffixed, content preserved."""
        expected = _write_per_task_files(
            tmp_path, fortran_micro.out_file_code, num_children=1
        )
        _write_per_task_logs(tmp_path, fortran_micro.out_file_code, num_children=1)

        fortran_micro.concatenate_child_outputs(tmp_path, num_children=1)

        for dataset in MICROSCALE_BINARY_DATASETS:
            unsuffixed = (
                tmp_path / f"{dataset}{fortran_micro.out_file_code}.dat"
            )
            suffixed = (
                tmp_path / f"{dataset}{fortran_micro.out_file_code}__00.dat"
            )
            assert unsuffixed.exists()
            assert not suffixed.exists()
            assert unsuffixed.read_bytes() == expected[dataset]

        log = tmp_path / f"micro{fortran_micro.out_file_code}.txt"
        log_src = tmp_path / f"micro{fortran_micro.out_file_code}__00.txt"
        assert log.exists()
        assert not log_src.exists()

    def test_concat_zero_children_raises_value_error(
        self, tmp_path, fortran_micro
    ):
        with pytest.raises(ValueError):
            fortran_micro.concatenate_child_outputs(tmp_path, num_children=0)
