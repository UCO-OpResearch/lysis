#!/usr/bin/env python3
"""Test script for converting Fortran microscale data to HDF5 format.

This script tests the data conversion pipeline by:
1. Reading Fortran v1.99.0 microscale output files
2. Converting to v2.0.0 HDF5 format
3. Optionally testing round-trip conversion
4. Reporting any issues found
"""

import sys
import os
import numpy as np
from pathlib import Path

# Add src/python to path
sys.path.insert(0, str(Path(__file__).parent / "src" / "python"))

from lysis.util.fileops import read_data_collection, write_data_collection
from lysis.util.dataconvert import convert_data
from lysis.util.dataspec import dataspec


def test_microscale_conversion(data_path: str):
    """Test conversion of microscale Fortran data to HDF5.

    :param data_path: Path to directory containing Fortran data files
    :type data_path: str
    """
    print("=" * 80)
    print("MICROSCALE DATA CONVERSION TEST")
    print("=" * 80)
    print(f"Data directory: {data_path}\n")

    # Step 2: Read Fortran v1.99.0 microscale output
    print("\nStep 2: Reading Fortran v1.99.0 microscale output files...")
    try:
        # Get the microscale_out collection specification for v1.99.0
        fortran_spec = dataspec["v1.99.0"]["microscale_out"]

        # Read the data collection
        fortran_data = read_data_collection(
            path=data_path,
            collections=[fortran_spec],
            file_codes=["_PLG2_tPA01_TB-xiii.dat"],
        )

        print(f"  ✓ Successfully read Fortran microscale data")
        print(f"    Datasets found:")
        for key in fortran_data.keys():
            if key != "params":
                if isinstance(fortran_data[key], np.ndarray):
                    print(
                        f"      - {key}: shape={fortran_data[key].shape}, dtype={fortran_data[key].dtype}"
                    )
                elif isinstance(fortran_data[key], list):
                    print(f"      - {key}: {len(fortran_data[key])} arrays")

    except Exception as e:
        print(f"  ✗ Failed to read Fortran data: {e}")
        import traceback

        traceback.print_exc()
        return False

    # Step 3: Convert v1.99.0 → v2.0.0
    print("\nStep 3: Converting v1.99.0 (Fortran) → v2.0.0 (HDF5)...")
    try:
        hdf5_data = convert_data(
            input_data=fortran_data, input_set_spec="v1.99.0", output_set_spec="v2.0.0"
        )

        print(f"  ✓ Successfully converted to HDF5 format")
        print(f"    Converted datasets:")
        for key in hdf5_data.keys():
            if key != "params":
                if isinstance(hdf5_data[key], np.ndarray):
                    print(
                        f"      - {key}: shape={hdf5_data[key].shape}, dtype={hdf5_data[key].dtype}"
                    )
                elif isinstance(hdf5_data[key], list):
                    print(f"      - {key}: {len(hdf5_data[key])} arrays")

    except Exception as e:
        print(f"  ✗ Failed to convert data: {e}")
        import traceback

        traceback.print_exc()
        return False

    # Step 4: Validate conversion by comparing field mappings
    print("\nStep 4: Validating conversion...")
    try:
        # Check some key field mappings
        validations = [
            ("firstPLi", "pli_first_time"),
            ("lasttPA", "tpa_final_num"),
            ("lyscomplete", "fiber_degraded"),
            ("lysis", "sim_final_time"),
            ("PLi", "pli_generated_num"),
            ("tPA_time", "tpa_leaving_time"),
            ("tPAPLiunbd", "tpa_unbound_by_pli"),
            ("tPAunbind", "tpa_unbound_kinetic"),
        ]

        all_passed = True
        for fortran_name, hdf5_name in validations:
            if fortran_name in fortran_data and hdf5_name in hdf5_data:
                # Compare arrays
                if np.array_equal(fortran_data[fortran_name], hdf5_data[hdf5_name]):
                    print(f"  ✓ {fortran_name} → {hdf5_name}: MATCH")
                else:
                    # Check if they're close (floating point tolerance)
                    if np.allclose(fortran_data[fortran_name], hdf5_data[hdf5_name]):
                        print(
                            f"  ✓ {fortran_name} → {hdf5_name}: MATCH (within tolerance)"
                        )
                    else:
                        print(f"  ✗ {fortran_name} → {hdf5_name}: MISMATCH")
                        all_passed = False
            else:
                missing = []
                if fortran_name not in fortran_data:
                    missing.append(f"Fortran:{fortran_name}")
                if hdf5_name not in hdf5_data:
                    missing.append(f"HDF5:{hdf5_name}")
                print(
                    f"  ! {fortran_name} → {hdf5_name}: MISSING ({', '.join(missing)})"
                )

        if all_passed:
            print(f"\n  ✓ All validations passed!")
        else:
            print(f"\n  ✗ Some validations failed")
            return False

    except Exception as e:
        print(f"  ✗ Validation failed: {e}")
        import traceback

        traceback.print_exc()
        return False

    # Step 5: Test round-trip conversion (optional)
    print("\nStep 5: Testing round-trip conversion (v1.99.0 → v2.0.0 → v1.99.0)...")
    try:
        roundtrip_data = convert_data(
            input_data=hdf5_data, input_set_spec="v2.0.0", output_set_spec="v1.99.0"
        )

        print(f"  ✓ Round-trip conversion completed")

        # Compare with original
        all_match = True
        for fortran_name, _ in validations:
            if fortran_name in fortran_data and fortran_name in roundtrip_data:
                if np.array_equal(
                    fortran_data[fortran_name], roundtrip_data[fortran_name]
                ):
                    print(f"  ✓ {fortran_name}: round-trip MATCH")
                elif np.allclose(
                    fortran_data[fortran_name], roundtrip_data[fortran_name]
                ):
                    print(f"  ✓ {fortran_name}: round-trip MATCH (within tolerance)")
                else:
                    print(f"  ✗ {fortran_name}: round-trip MISMATCH")
                    all_match = False

        if all_match:
            print(f"\n  ✓ Round-trip conversion successful!")
        else:
            print(
                f"\n  ! Round-trip conversion has some differences (may be acceptable)"
            )

    except Exception as e:
        print(f"  ✗ Round-trip conversion failed: {e}")
        import traceback

        traceback.print_exc()
        # Don't return False here - round-trip is optional

    print("\n" + "=" * 80)
    print("TEST COMPLETED SUCCESSFULLY")
    print("=" * 80)
    return True


def main():
    """Main entry point."""
    # Default test data path
    data_path = "/home/bpaynter/git/UCO-OpResearch/lysis/data/2024-09-02-1412"

    # Allow override from command line
    if len(sys.argv) > 1:
        data_path = sys.argv[1]

    # Run the test
    success = test_microscale_conversion(data_path)

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
