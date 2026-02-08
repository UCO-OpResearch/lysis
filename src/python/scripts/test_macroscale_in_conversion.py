#!/usr/bin/env python3
"""Test script for Macroscale Input data generation and conversion.

This script tests the macroscale input generation pipeline by:
1. Reading Fortran v1.99.0 microscale output files
2. Converting to v2.0.0 HDF5 format
3. Generating macroscale input data using generate_macroscale_in()
4. Converting macroscale input back to v1.99.0 Fortran format
5. Writing generated files to a test directory
6. Round-trip test: reading files back and converting (v1.99.0 → v2.0.0 → v1.99.0)
7. Reading original Fortran macroscale input files
8. Comparing generated files with original Fortran macroscale input files

The script validates that:
- The macroscale input generation pipeline works correctly
- Data can be written and read back without loss (round-trip integrity)
- Generated output matches the original Fortran-generated files

Test data directory should contain:
- Microscale output: {dataset_name}_PLG2_tPA01_TB-xiii.dat
- Macroscale input: tPAleave_PLG2_tPA01_TB-xiii.dat, etc.
"""

import sys
import os
import shutil
import numpy as np
from pathlib import Path

# Add parent directory (src/python) to path so we can import lysis
sys.path.insert(0, str(Path(__file__).parent.parent))

from lysis.util.fileops import read_data_collection, write_data_collection
from lysis.util.dataconvert import convert_data, generate_macroscale_in
from lysis.util.dataspec import dataspec


def compare_files(file1: Path, file2: Path, dataset_name: str, tolerance: float = 1e-9):
    """Compare two data files for equality.

    :param file1: Path to first file (reference)
    :type file1: Path
    :param file2: Path to second file (generated)
    :type file2: Path
    :param dataset_name: Name of dataset for reporting
    :type dataset_name: str
    :param tolerance: Tolerance for floating point comparisons
    :type tolerance: float
    :return: True if files match, False otherwise
    :rtype: bool
    """
    # Read both files
    with open(file1, "r") as f:
        data1 = np.fromfile(f, dtype=np.float64)
    with open(file2, "r") as f:
        data2 = np.fromfile(f, dtype=np.float64)

    # Check shapes match
    if data1.shape != data2.shape:
        print(f"  ✗ {dataset_name}: SHAPE MISMATCH")
        print(f"    Reference: {data1.shape}, Generated: {data2.shape}")
        return False

    # Check for exact match first
    if np.array_equal(data1, data2):
        print(f"  ✓ {dataset_name}: EXACT MATCH")
        return True

    # Check with tolerance
    if np.allclose(data1, data2, rtol=tolerance, atol=tolerance):
        max_diff = np.max(np.abs(data1 - data2))
        print(f"  ✓ {dataset_name}: MATCH (max diff: {max_diff:.2e})")
        return True

    # Report mismatch
    diff = np.abs(data1 - data2)
    max_diff = np.max(diff)
    max_idx = np.argmax(diff)
    print(f"  ✗ {dataset_name}: VALUES DIFFER")
    print(f"    Max difference: {max_diff:.2e} at index {max_idx}")
    print(f"    Reference[{max_idx}]: {data1.flat[max_idx]}")
    print(f"    Generated[{max_idx}]: {data2.flat[max_idx]}")
    return False


def test_macroscale_in_generation(data_path: str):
    """Test macroscale input generation and conversion pipeline.

    :param data_path: Path to directory containing test data
    :type data_path: str
    :return: True if all tests passed, False otherwise
    :rtype: bool
    """
    print("=" * 80)
    print("MACROSCALE INPUT GENERATION TEST")
    print("=" * 80)
    print(f"Data directory: {data_path}\n")

    data_path = Path(data_path)
    test_output_dir = data_path / "test_macro_in_generated"

    # Create test output directory (clean it if it exists)
    if test_output_dir.exists():
        print(f"Cleaning existing test output directory: {test_output_dir}")
        shutil.rmtree(test_output_dir)
    test_output_dir.mkdir()
    print(f"Created test output directory: {test_output_dir}\n")

    # Step 1: Read Fortran v1.99.0 microscale output
    print("Step 1: Reading Fortran v1.99.0 microscale output files...")
    try:
        microscale_spec = dataspec["v1.99.0"]["microscale_out"]
        microscale_data = read_data_collection(
            path=str(data_path),
            collections=[microscale_spec],
            file_codes=["__TF-v__9_951"],
        )
        print(f"  ✓ Successfully read {len(microscale_data) - 1} microscale datasets")
    except Exception as e:
        print(f"  ✗ Failed to read microscale data: {e}")
        import traceback

        traceback.print_exc()
        return False

    # Step 2: Convert v1.99.0 → v2.0.0
    print("\nStep 2: Converting microscale output v1.99.0 → v2.0.0...")
    try:
        microscale_v2 = convert_data(
            input_data=microscale_data,
            input_set_spec="v1.99.0",
            output_set_spec="v2.0.0",
        )
        print(f"  ✓ Successfully converted to v2.0.0 format")
    except Exception as e:
        print(f"  ✗ Failed to convert data: {e}")
        import traceback

        traceback.print_exc()
        return False

    # Step 3: Generate macroscale input from microscale output
    print("\nStep 3: Generating macroscale input data...")
    try:
        macroscale_in_v2 = generate_macroscale_in(microscale_v2)
        print(f"  ✓ Successfully generated macroscale input")
        print(f"    Generated datasets:")
        for key in macroscale_in_v2.keys():
            if key != "params":
                if isinstance(macroscale_in_v2[key], np.ndarray):
                    print(
                        f"      - {key}: shape={macroscale_in_v2[key].shape}, dtype={macroscale_in_v2[key].dtype}"
                    )
    except Exception as e:
        print(f"  ✗ Failed to generate macroscale input: {e}")
        import traceback

        traceback.print_exc()
        return False

    # Step 4: Convert macroscale input v2.0.0 → v1.99.0
    print("\nStep 4: Converting macroscale input v2.0.0 → v1.99.0...")
    try:
        macroscale_in_v1 = convert_data(
            input_data=macroscale_in_v2,
            input_set_spec="v2.0.0",
            output_set_spec="v1.99.0",
        )
        print(f"  ✓ Successfully converted to v1.99.0 Fortran format")
        print(f"    Generated datasets:")
        for key in macroscale_in_v1.keys():
            if key != "params":
                if isinstance(macroscale_in_v1[key], np.ndarray):
                    print(
                        f"      - {key}: shape={macroscale_in_v1[key].shape}, dtype={macroscale_in_v1[key].dtype}"
                    )
    except Exception as e:
        print(f"  ✗ Failed to convert to Fortran format: {e}")
        import traceback

        traceback.print_exc()
        return False

    # Step 5: Write generated files to test directory
    print(f"\nStep 5: Writing generated files to {test_output_dir}...")
    try:
        macroscale_in_spec = dataspec["v1.99.0"]["macroscale_in"]
        write_data_collection(
            path=str(test_output_dir),
            data=macroscale_in_v1,
            collections=[macroscale_in_spec],
            file_codes=["_TF-v__9_951"],
        )
        print(f"  ✓ Successfully wrote generated files")

        # List generated files
        generated_files = sorted(test_output_dir.glob("*.dat"))
        print(f"    Generated {len(generated_files)} files:")
        for f in generated_files:
            print(f"      - {f.name}")
    except Exception as e:
        print(f"  ✗ Failed to write files: {e}")
        import traceback

        traceback.print_exc()
        return False

    # Step 6: Round-trip test - read generated files and convert back
    print(f"\nStep 6: Testing round-trip conversion (v1.99.0 → v2.0.0 → v1.99.0)...")
    try:
        # Read the files we just wrote
        macroscale_in_read = read_data_collection(
            path=str(test_output_dir),
            collections=[macroscale_in_spec],
            file_codes=["_TF-v__9_951"],
        )
        print(f"  ✓ Read back generated files")
        print(f"    Generated datasets:")
        for key in macroscale_in_read.keys():
            if key != "params":
                if isinstance(macroscale_in_read[key], np.ndarray):
                    print(
                        f"      - {key}: shape={macroscale_in_read[key].shape}, dtype={macroscale_in_read[key].dtype}"
                    )
    except Exception as e:
        print(f"  ✗ Round-trip test failed: {e}")
        import traceback

        traceback.print_exc()
        return False

    try:
        # Convert to v2.0.0
        macroscale_in_v2_roundtrip = convert_data(
            input_data=macroscale_in_read,
            input_set_spec="v1.99.0",
            output_set_spec="v2.0.0",
        )
        print(f"  ✓ Converted to v2.0.0")
        print(f"    Generated datasets:")
        for key in macroscale_in_v2_roundtrip.keys():
            if key != "params":
                if isinstance(macroscale_in_v2_roundtrip[key], np.ndarray):
                    print(
                        f"      - {key}: shape={macroscale_in_v2_roundtrip[key].shape}, dtype={macroscale_in_v2_roundtrip[key].dtype}"
                    )
    except Exception as e:
        print(f"  ✗ Round-trip test failed: {e}")
        import traceback

        traceback.print_exc()
        return False

    try:
        # Convert back to v1.99.0
        macroscale_in_v1_roundtrip = convert_data(
            input_data=macroscale_in_v2_roundtrip,
            input_set_spec="v2.0.0",
            output_set_spec="v1.99.0",
        )
        print(f"  ✓ Converted back to v1.99.0")
        print(f"    Generated datasets:")
        for key in macroscale_in_v1_roundtrip.keys():
            if key != "params":
                if isinstance(macroscale_in_v1_roundtrip[key], np.ndarray):
                    print(
                        f"      - {key}: shape={macroscale_in_v1_roundtrip[key].shape}, dtype={macroscale_in_v1_roundtrip[key].dtype}"
                    )
    except Exception as e:
        print(f"  ✗ Round-trip test failed: {e}")
        import traceback

        traceback.print_exc()
        return False

    try:
        # Compare with original in-memory data
        datasets_for_roundtrip = [
            "tPAleave",
            "tsectPA",
            "lysismat",
            "lenlysisvect",
            "neighbors",
        ]

        roundtrip_passed = True
        for dataset_name in datasets_for_roundtrip:
            if (
                dataset_name in macroscale_in_v1
                and dataset_name in macroscale_in_v1_roundtrip
            ):
                if np.array_equal(
                    macroscale_in_v1[dataset_name],
                    macroscale_in_v1_roundtrip[dataset_name],
                ):
                    print(f"  ✓ {dataset_name}: round-trip MATCH")
                elif np.allclose(
                    macroscale_in_v1[dataset_name],
                    macroscale_in_v1_roundtrip[dataset_name],
                ):
                    print(f"  ✓ {dataset_name}: round-trip MATCH (within tolerance)")
                else:
                    print(f"  ✗ {dataset_name}: round-trip MISMATCH")
                    max_diff = np.max(
                        np.abs(
                            macroscale_in_v1[dataset_name]
                            - macroscale_in_v1_roundtrip[dataset_name]
                        )
                    )
                    print(f"    Max difference: {max_diff}")
                    roundtrip_passed = False

        if not roundtrip_passed:
            print(f"\n  ! Round-trip test failed - data integrity issue")
            return False

        print(f"\n  ✓ Round-trip conversion successful!")

    except Exception as e:
        print(f"  ✗ Round-trip test failed: {e}")
        import traceback

        traceback.print_exc()
        return False

    # Step 7: Read original Fortran macroscale input files
    print(f"\nStep 7: Reading original Fortran macroscale input files...")
    try:
        macroscale_in_original = read_data_collection(
            path=str(data_path),
            collections=[macroscale_in_spec],
            file_codes=["_TF-v__9_951"],
        )
        print(f"  ✓ Successfully read original macroscale input files")
    except Exception as e:
        print(f"  ✗ Failed to read original files: {e}")
        import traceback

        traceback.print_exc()
        return False

    # Step 8: Compare generated files with original files
    print("\nStep 8: Comparing generated files with original files...")

    # Get list of datasets to compare (exclude params and neighbors)
    # neighbors is generated from params, not from microscale output
    datasets_to_compare = [
        "tPAleave",
        "tsectPA",
        "lysismat",
        "lenlysisvect",
        "neighbors",
    ]

    all_match = True
    for dataset_name in datasets_to_compare:
        original_file = data_path / f"{dataset_name}_TF-v__9_951.dat"
        generated_file = test_output_dir / f"{dataset_name}_TF-v__9_951.dat"

        if not original_file.exists():
            print(f"  ! {dataset_name}: Original file not found (skipping)")
            continue

        if not generated_file.exists():
            print(f"  ✗ {dataset_name}: Generated file not found")
            all_match = False
            continue

        if not compare_files(original_file, generated_file, dataset_name):
            all_match = False

    # Summary
    print("\n" + "=" * 80)
    if all_match:
        print("TEST COMPLETED SUCCESSFULLY")
        print("All generated files match the original Fortran files!")
    else:
        print("TEST COMPLETED WITH ERRORS")
        print("Some generated files differ from the original Fortran files.")
    print("=" * 80)

    return all_match


def main():
    """Main entry point."""
    # Default test data path
    data_path = "/home/bpaynter/git/UCO-OpResearch/lysis/data/_TF-v__9_951"

    # Allow override from command line
    if len(sys.argv) > 1:
        data_path = sys.argv[1]

    # Run the test
    success = test_macroscale_in_generation(data_path)

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
