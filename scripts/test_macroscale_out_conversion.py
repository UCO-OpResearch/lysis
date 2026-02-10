#!/usr/bin/env python3
"""Test script for Macroscale Output data conversion.

This script tests the macroscale output conversion pipeline by:
1. Reading Fortran v1.99.0 macroscale output files
2. Converting to v2.0.0 format
3. Converting back to v1.99.0 format (round-trip)
4. Comparing round-trip results with originals in memory
5. Writing round-trip results to a test directory
6. Reading written files back and comparing with originals (file-level)

The round-trip test validates that:
- All macroscale output converters work correctly in both directions
- Data can survive v1.99.0 -> v2.0.0 -> v1.99.0 without loss
- m_bound reconstruction from tpa_bind_events matches original m_bound

Test data directory should contain:
- Simulation folders (00, 01, ...) with macroscale output files
- params.json with simulation parameters
"""

import sys
import shutil
import numpy as np
from pathlib import Path

# Add src directory to path so we can import lysis
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from lysis.util.fileops import read_data_collection, write_data_collection
from lysis.util.dataconvert import convert_data
from lysis.util.dataspec import dataspec, DataSetSpec, parse_shape
from lysis.util.constants import CONST


def compare_arrays(
    original, converted, dataset_name: str, tolerance: float = 1e-9
) -> bool:
    """Compare two arrays or lists of arrays for equality.

    For list data (simulations_combined=False), compares each simulation
    element-by-element and reports per-simulation results only on failure.

    :param original: Original data (array or list of arrays)
    :param converted: Converted data (array or list of arrays)
    :param dataset_name: Name of dataset for reporting
    :type dataset_name: str
    :param tolerance: Tolerance for floating point comparisons
    :type tolerance: float
    :return: True if all data matches, False otherwise
    :rtype: bool
    """
    if isinstance(original, list) and isinstance(converted, list):
        if len(original) != len(converted):
            print(
                f"  \u2717 {dataset_name}: SIMULATION COUNT MISMATCH "
                f"({len(original)} vs {len(converted)})"
            )
            return False
        all_match = True
        for i, (orig_arr, conv_arr) in enumerate(zip(original, converted)):
            if not _compare_single_array(
                orig_arr, conv_arr, f"{dataset_name}[sim {i:02}]", tolerance
            ):
                all_match = False
        if all_match:
            print(
                f"  \u2713 {dataset_name}: ALL SIMULATIONS MATCH "
                f"({len(original)} sims)"
            )
        return all_match
    else:
        return _compare_single_array(original, converted, dataset_name, tolerance)


def _compare_single_array(
    arr1, arr2, label: str, tolerance: float = 1e-9
) -> bool:
    """Compare two single arrays for equality.

    Handles structured arrays by comparing field by field, string arrays by
    exact comparison, and numeric arrays by tolerance-based comparison.

    :param arr1: First array (reference)
    :param arr2: Second array (test)
    :param label: Label for reporting
    :type label: str
    :param tolerance: Tolerance for floating point comparisons
    :type tolerance: float
    :return: True if arrays match, False otherwise
    :rtype: bool
    """
    arr1 = np.asarray(arr1)
    arr2 = np.asarray(arr2)

    if arr1.shape != arr2.shape:
        print(f"  \u2717 {label}: SHAPE MISMATCH ({arr1.shape} vs {arr2.shape})")
        return False

    # Structured arrays: compare field by field
    if arr1.dtype.names:
        all_match = True
        for field in arr1.dtype.names:
            if not _compare_single_array(
                arr1[field], arr2[field], f"{label}.{field}", tolerance
            ):
                all_match = False
        return all_match

    # String arrays: exact comparison
    if arr1.dtype.kind in ("U", "S", "O"):
        if np.array_equal(arr1, arr2):
            return True
        mismatches = np.sum(arr1 != arr2)
        print(f"  \u2717 {label}: {mismatches} STRING MISMATCHES")
        return False

    # Numeric arrays
    if np.array_equal(arr1, arr2):
        return True

    if np.allclose(arr1, arr2, rtol=tolerance, atol=tolerance):
        max_diff = np.max(np.abs(arr1.astype(float) - arr2.astype(float)))
        print(f"  ~ {label}: CLOSE MATCH (max diff: {max_diff:.2e})")
        return True

    # Report mismatch details
    diff = np.abs(arr1.astype(float) - arr2.astype(float))
    max_diff = np.max(diff)
    max_idx = np.unravel_index(np.argmax(diff), diff.shape)
    print(f"  \u2717 {label}: VALUES DIFFER")
    print(f"    Max difference: {max_diff:.2e} at index {max_idx}")
    print(f"    Original: {arr1[max_idx]}")
    print(f"    Converted: {arr2[max_idx]}")
    return False


def compare_files(
    file1: Path,
    file2: Path,
    spec: DataSetSpec,
    dataset_name: str,
    params=None,
    tolerance: float = 1e-9,
) -> bool:
    """Compare two data files for equality using the dataset specification.

    Reads both files using the appropriate numpy reader for the storage type
    (binary via ``np.fromfile``, text via ``np.loadtxt``), then compares the
    resulting arrays. This avoids false negatives from whitespace or
    formatting differences in text files.

    :param file1: Path to first file (reference)
    :type file1: Path
    :param file2: Path to second file (generated)
    :type file2: Path
    :param spec: Dataset specification (provides dtype, shape, storage type,
        and delimiter)
    :type spec: DataSetSpec
    :param dataset_name: Name of dataset for reporting
    :type dataset_name: str
    :param params: Simulation parameters for resolving dynamic shapes
    :param tolerance: Tolerance for floating point comparisons
    :type tolerance: float
    :return: True if files match, False otherwise
    :rtype: bool
    """
    storage = spec.dataset_storage_type
    try:
        if storage == CONST.DATASET_STORAGE_TYPE.FILE_BINARY:
            arr1 = np.fromfile(str(file1), dtype=spec.dtype)
            arr2 = np.fromfile(str(file2), dtype=spec.dtype)
            if spec.shape:
                shape = parse_shape(spec.shape, params=params)
                arr1 = arr1.reshape(shape)
                arr2 = arr2.reshape(shape)
        elif storage == CONST.DATASET_STORAGE_TYPE.FILE_TEXT:
            delimiter = spec.delimiter if spec.delimiter is not None else " "
            if spec.dtype == str:
                # String data: compare raw bytes
                with open(file1, "rb") as f:
                    raw1 = f.read()
                with open(file2, "rb") as f:
                    raw2 = f.read()
                if raw1 == raw2:
                    print(f"  \u2713 {dataset_name}: EXACT FILE MATCH")
                    return True
                print(f"  \u2717 {dataset_name}: TEXT FILES DIFFER")
                return False
            elif hasattr(spec.dtype, "names") and spec.dtype.names:
                # Structured array: read as text with delimiter
                arr1 = np.loadtxt(
                    str(file1), dtype=spec.dtype, delimiter=delimiter
                )
                arr2 = np.loadtxt(
                    str(file2), dtype=spec.dtype, delimiter=delimiter
                )
            else:
                arr1 = np.loadtxt(
                    str(file1), dtype=spec.dtype, delimiter=delimiter
                )
                arr2 = np.loadtxt(
                    str(file2), dtype=spec.dtype, delimiter=delimiter
                )
        else:
            # Fallback: raw byte comparison
            with open(file1, "rb") as f:
                raw1 = f.read()
            with open(file2, "rb") as f:
                raw2 = f.read()
            if raw1 == raw2:
                print(f"  \u2713 {dataset_name}: EXACT FILE MATCH")
                return True
            print(f"  \u2717 {dataset_name}: FILES DIFFER")
            return False
    except Exception as e:
        print(f"  \u2717 {dataset_name}: FAILED TO READ ({e})")
        return False

    result = _compare_single_array(arr1, arr2, dataset_name, tolerance)
    if result:
        print(f"  \u2713 {dataset_name}: DATA MATCH")
    return result


def test_macroscale_out_conversion(data_path: str, file_code: str) -> bool:
    """Test macroscale output conversion pipeline.

    :param data_path: Path to directory containing test data
    :type data_path: str
    :param file_code: File code for the data files
    :type file_code: str
    :return: True if all tests passed, False otherwise
    :rtype: bool
    """
    print("=" * 80)
    print("MACROSCALE OUTPUT CONVERSION TEST")
    print("=" * 80)
    print(f"Data directory: {data_path}")
    print(f"File code: {file_code}\n")

    data_path = Path(data_path)
    test_output_dir = data_path / "test_macro_out_roundtrip"
    all_passed = True

    # Create test output directory (clean it if it exists)
    if test_output_dir.exists():
        print(f"Cleaning existing test output directory: {test_output_dir}")
        shutil.rmtree(test_output_dir)

    # ---- Step 1: Read v1.99.0 macroscale output ----
    print("Step 1: Reading Fortran v1.99.0 macroscale output files...")
    try:
        macro_out_spec_v199 = dataspec["v1.99.0"]["macroscale_out"]
        original_data = read_data_collection(
            path=str(data_path),
            collections=[macro_out_spec_v199],
            file_codes=[file_code],
        )
        # Report what was read
        n_datasets = sum(1 for k in original_data if k != "params")
        print(f"  \u2713 Successfully read {n_datasets} datasets")
        for key in sorted(original_data.keys()):
            if key == "params":
                continue
            val = original_data[key]
            if isinstance(val, list):
                if len(val) > 0:
                    first = val[0]
                    print(
                        f"    {key}: {len(val)} sims, "
                        f"shape={first.shape}, dtype={first.dtype}"
                    )
                else:
                    print(f"    {key}: 0 sims")
            elif isinstance(val, np.ndarray):
                print(f"    {key}: shape={val.shape}, dtype={val.dtype}")
    except Exception as e:
        print(f"  \u2717 Failed to read macroscale output: {e}")
        import traceback

        traceback.print_exc()
        return False

    # ---- Step 2: Convert v1.99.0 -> v2.0.0 ----
    print("\nStep 2: Converting macroscale output v1.99.0 \u2192 v2.0.0...")
    try:
        v2_data = convert_data(
            input_data=original_data,
            input_set_spec="v1.99.0",
            output_set_spec="v2.0.0",
        )
        n_datasets = sum(1 for k in v2_data if k != "params")
        print(f"  \u2713 Successfully converted to v2.0.0 ({n_datasets} datasets)")
        for key in sorted(v2_data.keys()):
            if key == "params":
                continue
            val = v2_data[key]
            if isinstance(val, list):
                if len(val) > 0:
                    first = val[0]
                    print(
                        f"    {key}: {len(val)} sims, "
                        f"shape={first.shape}, dtype={first.dtype}"
                    )
            elif isinstance(val, np.ndarray):
                print(f"    {key}: shape={val.shape}, dtype={val.dtype}")
    except Exception as e:
        print(f"  \u2717 Failed to convert to v2.0.0: {e}")
        import traceback

        traceback.print_exc()
        return False

    # ---- Step 3: Convert v2.0.0 -> v1.99.0 ----
    print("\nStep 3: Converting macroscale output v2.0.0 \u2192 v1.99.0...")
    try:
        roundtrip_data = convert_data(
            input_data=v2_data,
            input_set_spec="v2.0.0",
            output_set_spec="v1.99.0",
        )
        n_datasets = sum(1 for k in roundtrip_data if k != "params")
        print(
            f"  \u2713 Successfully converted back to v1.99.0 ({n_datasets} datasets)"
        )
        for key in sorted(roundtrip_data.keys()):
            if key == "params":
                continue
            val = roundtrip_data[key]
            if isinstance(val, list):
                if len(val) > 0:
                    first = val[0]
                    print(
                        f"    {key}: {len(val)} sims, "
                        f"shape={first.shape}, dtype={first.dtype}"
                    )
            elif isinstance(val, np.ndarray):
                print(f"    {key}: shape={val.shape}, dtype={val.dtype}")
    except Exception as e:
        print(f"  \u2717 Failed to convert back to v1.99.0: {e}")
        import traceback

        traceback.print_exc()
        return False

    # ---- Step 4: In-memory comparison ----
    print("\nStep 4: Comparing round-trip data with originals (in-memory)...")

    # Datasets that should survive the round-trip exactly
    datasets_to_compare = [
        "macro_log",
        "Nsave",
        "tsave",
        "f_deg_list",
        "m_bind_t",
        "m_loc",
        "m_bound",
        "mfpt",
    ]

    for dataset_name in datasets_to_compare:
        if dataset_name not in original_data:
            print(f"  ! {dataset_name}: Not in original data (skipping)")
            continue
        if dataset_name not in roundtrip_data:
            print(f"  \u2717 {dataset_name}: Not in round-trip data")
            all_passed = False
            continue
        if not compare_arrays(
            original_data[dataset_name],
            roundtrip_data[dataset_name],
            dataset_name,
        ):
            all_passed = False

    # ---- Step 5: Write round-trip files ----
    print(f"\nStep 5: Writing round-trip files to {test_output_dir}...")
    write_succeeded = False
    try:
        # Create simulation subdirectories
        n_sims = len(original_data["tsave"])
        for sim in range(n_sims):
            (test_output_dir / f"{sim:02}").mkdir(parents=True, exist_ok=True)

        write_data_collection(
            data=roundtrip_data,
            path=str(test_output_dir),
            collections=[macro_out_spec_v199],
            file_codes=[file_code],
        )
        print(f"  \u2713 Successfully wrote round-trip files")
        write_succeeded = True
    except NotImplementedError as e:
        print(f"  ! Write skipped (JSON writer not yet implemented)")
        print(f"    File-level comparison will be skipped.")
    except Exception as e:
        print(f"  \u2717 Failed to write files: {e}")
        import traceback

        traceback.print_exc()
        all_passed = False

    # ---- Step 6: File-level comparison ----
    if write_succeeded and test_output_dir.exists():
        print("\nStep 6: Comparing written files with originals...")

        params = original_data["params"]
        for dataset_name, spec in macro_out_spec_v199.data.items():
            for sim in range(n_sims):
                loc = spec.data_location.format(
                    sim=sim, file_code=file_code
                )
                orig_file = data_path / loc
                test_file = test_output_dir / loc

                if not orig_file.exists():
                    continue
                if not test_file.exists():
                    print(
                        f"  \u2717 {dataset_name}[sim {sim:02}]: "
                        f"Generated file not found"
                    )
                    all_passed = False
                    continue
                if not compare_files(
                    orig_file,
                    test_file,
                    spec,
                    f"{dataset_name}[sim {sim:02}]",
                    params=params,
                ):
                    all_passed = False
    else:
        print(
            "\nStep 6: Skipped (file write not available; "
            "JSON writer needed for params.json)"
        )

    # ---- Summary ----
    print("\n" + "=" * 80)
    if all_passed:
        print("ALL TESTS PASSED")
        print(
            "Round-trip conversion v1.99.0 \u2192 v2.0.0 \u2192 v1.99.0 "
            "preserves all data!"
        )
    else:
        print("SOME TESTS FAILED")
        print("Check the output above for details on failures.")
    print("=" * 80)

    return all_passed


def main():
    """Main entry point."""
    # Default test data path and file code
    data_path = "/home/bpaynter/git/UCO-OpResearch/lysis/data/2024-09-02-1412"
    file_code = "_TB-xiii__21_105"

    # Allow override from command line
    if len(sys.argv) > 1:
        data_path = sys.argv[1]
    if len(sys.argv) > 2:
        file_code = sys.argv[2]

    success = test_macroscale_out_conversion(data_path, file_code)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
