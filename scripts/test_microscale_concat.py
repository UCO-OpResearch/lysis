#!/usr/bin/env python3
"""Test script for amalgamating microscale Fortran binary outputs.

Motivation: we'd like to execute the microscale Fortran simulation as a
slurm array (one task per seed) and then amalgamate the per-task binary
output files before handing them to :func:`read_data_collection` for
HDF5 import.  It would be wasteful to give each batch its own ``Run``/
HDF5 file (50k per simulation x many seeds), and we don't want to invent
a new dataspec just for this case.  This script verifies that the naive
"concatenate the raw ``.dat`` files end-to-end" approach is
byte-equivalent to what dataio would produce from a single larger run.

What this script does:

1. Executes ``N_RUNS`` microscale Fortran simulations, each in its own
   working directory, each with ``micro_simulations=SIMS_PER_RUN`` and
   a distinct ``micro_seed``.  All other parameters are defaults.
2. Reads each run's output through :func:`read_data_collection` and
   stores the per-dataset arrays.
3. Byte-concatenates the raw ``.dat`` files (``firstPLi``, ``lasttPA``,
   ``lyscomplete``, ``lysis``, ``PLi``, ``tPA_time``, ``tPAPLiunbd``,
   ``tPAunbind``) from the N runs into a single amalgamated directory,
   and writes a fresh ``params.json`` there with
   ``micro_simulations = N_RUNS * SIMS_PER_RUN``.
4. Reads the amalgamated directory with :func:`read_data_collection`.
5. Compares the amalgamated arrays to ``numpy.concatenate`` of the
   individual arrays, element-by-element, per dataset.

Scope: only the eight per-simulation binary datasets in
``microscale_out`` are checked.  The ``macroscale_in`` files
(``tPAleave``, ``tsectPA``, ``lysismat``, ``lenlysisvect``) are
per-run *aggregates* (histograms / percentile tables) and CANNOT be
byte-concatenated; a slurm amalgamation pipeline must recompute them
from the combined microscale output.

Usage:

    python scripts/test_microscale_concat.py [--binary BIN] [--sims N] [--runs N]

Exits 0 on success, 1 on any mismatch.
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

from lysis.config.parameters import MicroParameters
from lysis.config.run import Run
from lysis.dataio.dataspec import dataspec
from lysis.dataio.fileops import read_data_collection, write_dataset
from lysis.execution.fortran import MICRO_FORTRAN_DATASPEC_VERSION
from lysis.execution.fortran_micro import FortranMicro


__author__ = "Bradley Paynter"
__copyright__ = "Copyright 2026, Brittany Bannish"
__license__ = "GPLv3"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


# Binary per-simulation datasets in microscale_out — these are the only
# datasets whose raw files concatenate meaningfully.  Every one is shape
# (-1,) so numpy.fromfile auto-infers length from file size.
BINARY_DATASETS = [
    "firstPLi",
    "lasttPA",
    "lyscomplete",
    "lysis",
    "PLi",
    "tPA_time",
    "tPAPLiunbd",
    "tPAunbind",
]


def _execute_one_run(
    workdir: Path, binary: str, run_code: str, seed: int, sims: int, file_code: str
) -> Path:
    """Run the Fortran micro binary once and return the output data dir."""
    workdir.mkdir(parents=True, exist_ok=True)
    run = Run(str(workdir), run_code=run_code)
    run.initialize_micro_param({"micro_simulations": sims, "micro_seed": seed})

    fm = FortranMicro(
        run=run,
        executable=binary,
        out_file_code=file_code,
    )
    data_dir = fm.exec_in_workdir(workdir)
    return data_dir


def _concatenate_binaries(
    run_dirs: list, combined_dir: Path, file_code: str, total_sims: int
) -> None:
    """Byte-concatenate each binary dataset from ``run_dirs`` into
    ``combined_dir`` and write a combined ``params.json``.
    """
    combined_dir.mkdir(parents=True, exist_ok=True)

    for dataset in BINARY_DATASETS:
        filename = f"{dataset}{file_code}.dat"
        dst = combined_dir / filename
        with open(dst, "wb") as out_fh:
            for src_dir in run_dirs:
                src = src_dir / filename
                if not src.exists():
                    raise FileNotFoundError(f"Expected {src} to exist after run")
                with open(src, "rb") as in_fh:
                    shutil.copyfileobj(in_fh, out_fh)

    # read_data_collection also reads the micro_log text file; concatenate it
    # so the combined dir is a complete microscale_out collection.
    log_name = f"micro{file_code}.txt"
    with open(combined_dir / log_name, "wb") as out_fh:
        for src_dir in run_dirs:
            src = src_dir / log_name
            if not src.exists():
                continue
            with open(src, "rb") as in_fh:
                shutil.copyfileobj(in_fh, out_fh)

    # Write a fresh params.json with the amalgamated simulation count so
    # read_data_collection resolves shapes against the correct total.
    combined_params = MicroParameters(
        micro_simulations=total_sims,
        # Seed here is informational for the combined file — the actual
        # seeds that produced the bytes are preserved in each input file.
        micro_seed=0,
    )
    write_dataset(
        {"micro_params": combined_params.to_basedict()},
        str(combined_dir),
        dataspec[MICRO_FORTRAN_DATASPEC_VERSION]["microscale_out"].params,
    )


def _read_binary_datasets(data_dir: Path, file_code: str) -> dict:
    """Read microscale_out from ``data_dir`` and return binary datasets only."""
    data = read_data_collection(
        path=str(data_dir),
        collections=[dataspec[MICRO_FORTRAN_DATASPEC_VERSION]["microscale_out"]],
        file_codes=[file_code],
    )
    return {name: np.asarray(data[name]) for name in BINARY_DATASETS}


def run_test(binary: str, n_runs: int, sims_per_run: int, base_seed: int) -> bool:
    print("=" * 80)
    print("MICROSCALE BINARY CONCATENATION TEST")
    print("=" * 80)
    print(f"  binary       : {binary}")
    print(f"  N runs       : {n_runs}")
    print(f"  sims per run : {sims_per_run}")
    print(f"  base seed    : {base_seed}")
    print()

    if not Path(binary).exists():
        print(f"  ✗ Binary not found: {binary}")
        return False

    file_code = "_CAT"

    with tempfile.TemporaryDirectory(prefix="lysis_cat_test_") as tmp:
        tmp_path = Path(tmp)

        # Step 1: Run N microscale simulations with different seeds.
        print(f"Step 1: Executing {n_runs} microscale simulations...")
        run_data_dirs = []
        try:
            for i in range(n_runs):
                run_code = f"sim-{i:02d}"
                seed = base_seed + i
                workdir = tmp_path / run_code
                data_dir = _execute_one_run(
                    workdir, binary, run_code, seed, sims_per_run, file_code
                )
                run_data_dirs.append(data_dir)
                print(f"  ✓ {run_code}: seed={seed}, output={data_dir}")
        except subprocess.CalledProcessError as exc:
            print(f"  ✗ Fortran execution failed: {exc}")
            return False
        except Exception as exc:
            print(f"  ✗ Failed to execute runs: {exc}")
            import traceback
            traceback.print_exc()
            return False

        # Step 2: Read each run individually with dataio.
        print(f"\nStep 2: Reading each run with dataio...")
        try:
            individual = [_read_binary_datasets(d, file_code) for d in run_data_dirs]
            for i, arrs in enumerate(individual):
                sample = arrs[BINARY_DATASETS[0]]
                print(
                    f"  ✓ run {i:02d}: "
                    f"{BINARY_DATASETS[0]}.shape={sample.shape}, dtype={sample.dtype}"
                )
        except Exception as exc:
            print(f"  ✗ Failed to read individual runs: {exc}")
            import traceback
            traceback.print_exc()
            return False

        # Step 3: Byte-concatenate the raw binary files.
        print(f"\nStep 3: Byte-concatenating raw binary files...")
        combined_dir = tmp_path / "combined"
        total_sims = n_runs * sims_per_run
        try:
            _concatenate_binaries(run_data_dirs, combined_dir, file_code, total_sims)
            print(f"  ✓ Combined dir: {combined_dir}")
            for dataset in BINARY_DATASETS:
                size = (combined_dir / f"{dataset}{file_code}.dat").stat().st_size
                print(f"    - {dataset}{file_code}.dat: {size} bytes")
        except Exception as exc:
            print(f"  ✗ Failed to concatenate: {exc}")
            import traceback
            traceback.print_exc()
            return False

        # Step 4: Read the combined directory with dataio.
        print(f"\nStep 4: Reading concatenated directory with dataio...")
        try:
            combined = _read_binary_datasets(combined_dir, file_code)
            for dataset in BINARY_DATASETS:
                print(
                    f"  ✓ {dataset}: shape={combined[dataset].shape}, "
                    f"dtype={combined[dataset].dtype}"
                )
        except Exception as exc:
            print(f"  ✗ Failed to read combined dir: {exc}")
            import traceback
            traceback.print_exc()
            return False

        # Step 5: Compare combined vs np.concatenate(individuals).
        print(f"\nStep 5: Comparing combined read vs numpy-concatenated reads...")
        all_match = True
        for dataset in BINARY_DATASETS:
            expected = np.concatenate([ind[dataset] for ind in individual])
            actual = combined[dataset]

            if actual.shape != expected.shape:
                print(
                    f"  ✗ {dataset}: shape mismatch — "
                    f"combined {actual.shape} vs expected {expected.shape}"
                )
                all_match = False
                continue

            if actual.dtype != expected.dtype:
                print(
                    f"  ✗ {dataset}: dtype mismatch — "
                    f"combined {actual.dtype} vs expected {expected.dtype}"
                )
                all_match = False
                continue

            if np.array_equal(actual, expected):
                print(
                    f"  ✓ {dataset}: MATCH  "
                    f"(n={actual.shape[0]}, dtype={actual.dtype})"
                )
            else:
                n_diff = int(np.sum(actual != expected))
                print(
                    f"  ✗ {dataset}: {n_diff} of {actual.size} elements differ"
                )
                all_match = False

    print("\n" + "=" * 80)
    if all_match:
        print("TEST PASSED — raw byte-concatenation is dataio-equivalent")
        print("=" * 80)
        return True
    else:
        print("TEST FAILED — see mismatches above")
        print("=" * 80)
        return False


def parse_arguments():
    repo_root = Path(__file__).resolve().parent.parent
    default_binary = repo_root / "bin" / "micro_rates"

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--binary",
        type=str,
        default=str(default_binary),
        help="Path to the compiled microscale Fortran binary.",
    )
    parser.add_argument(
        "--runs", type=int, default=5, help="Number of microscale runs to execute."
    )
    parser.add_argument(
        "--sims",
        type=int,
        default=10,
        help="Number of internal simulations per run (keep small — real Fortran).",
    )
    parser.add_argument(
        "--seed", type=int, default=9999, help="Base seed (run i uses seed+i)."
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    ok = run_test(
        binary=args.binary,
        n_runs=args.runs,
        sims_per_run=args.sims,
        base_seed=args.seed,
    )
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
