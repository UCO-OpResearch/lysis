"""Data conversion utilities for lysis simulation data.

This module provides tools for converting simulation data between different
specification versions and formats. The primary use cases are:

1. **Version conversion**: Converting data between v1.99.0 (Fortran file-based)
   and v2.0.0 (HDF5 unified) specifications

2. **Format transformation**: Transforming data structures to meet the requirements
   of different simulation stages (microscale → macroscale conversion)

3. **Safe type conversion**: Ensuring numerical data types are converted without
   overflow or precision loss

Conversion Architecture
-----------------------

The module uses a dictionary-based converter registry (``data_converters``) that
maps (input_spec, output_spec) tuples to dictionaries of dataset conversion
functions. Each converter function takes a ``DataCollectionType`` and returns a
``DataSetType`` (or list thereof).

Conversion functions can be:
- Simple lambdas for direct field mapping: ``lambda data: data["field_name"]``
- Partial functions with fixed parameters: ``functools.partial(converter, input_spec="v1.99.0")``
- Full converter functions for complex transformations: ``convert_fiber_degrade_time()``
- Not-implemented stubs: ``_not_implemented()`` for pending converters

Key Functions
-------------

- :func:`convert_data`: Main entry point for converting complete data collections
- :func:`generate_macroscale_in`: Generates macroscale input from microscale output
- :func:`convert_fiber_degrade_time`: Converts fiber degradation timing data
- :func:`safe_np_int_conversion`: Safely converts integer arrays with bounds checking
- :func:`safe_np_bool_conversion`: Safely converts boolean arrays with validation

Example Usage
-------------

Converting from v1.99.0 to v2.0.0 format::

    from lysis.util.dataconvert import convert_data

    # Load v1.99.0 data
    old_data = load_fortran_data()

    # Convert to v2.0.0 HDF5 format
    new_data = convert_data(old_data, "v1.99.0", "v2.0.0")

    # Write to HDF5 file
    write_hdf5(new_data, "output.h5")

Generating macroscale input from microscale output::

    from lysis.util.dataconvert import generate_macroscale_in

    # microscale_output contains simulation results
    macroscale_input = generate_macroscale_in(microscale_output)

    # Convert to Fortran format and write
    fortran_data = convert_data(macroscale_input, "v2.0.0", "v1.99.0")
    write_fortran_files(fortran_data)

Notes
-----

- Tag aliases (like "current", "fortran") are automatically resolved to version numbers
- Type conversions use safe methods that check bounds to prevent overflow
- Missing converters raise ``NotImplementedError`` with clear messages
- The module assumes data has already been validated against its input specification

See Also
--------

:mod:`lysis.util.dataspec` : Data specification definitions and validation
:mod:`lysis.util.edge_grid` : Grid neighborhood structure utilities
"""

import functools
import warnings

from dataclasses import asdict
from enum import Flag, auto, unique
from typing import Any, AnyStr, List, Mapping, Union, Callable

import numpy as np
import h5py

from .constants import CONST
from .dataspec import (
    DataCollectionSpec,
    DataCollectionType,
    DataSetSpec,
    DataSetType,
    dataspec,
    parse_shape,
    check_dataset_spec,
    tags,
)
from .edge_grid import (
    generate_fortran_neighborhood_structure,
    from_fortran_edge_index_array,
)

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2025, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


def safe_np_int_conversion(int_array, dtype=np.uint8, copy=True):
    """Safely convert integer arrays to a target dtype with bounds checking.

    NumPy's default casting can cause silent wraparound errors when converting
    between integer types (e.g., 300 → 44 when converting to uint8). This function
    checks that all values fit within the target dtype's range before converting.

    This is particularly important when converting from signed to unsigned types
    or when downcasting from larger to smaller integer types.

    :param int_array: Input array or array-like of integers to convert
    :type int_array: array_like
    :param dtype: Target NumPy integer dtype (e.g., np.int32, np.uint8)
    :type dtype: numpy.dtype, optional
    :param copy: Whether to copy the data (True) or reuse memory if possible (False)
    :type copy: bool, optional
    :return: Array converted to target dtype
    :rtype: numpy.ndarray
    :raises TypeError: If input cannot be safely converted to target dtype
    :raises OverflowError: If any values exceed target dtype's min/max bounds

    Examples
    --------
    Safe conversion within bounds::

        >>> arr = np.array([0, 100, 255])
        >>> safe_np_int_conversion(arr, dtype=np.uint8)
        array([  0, 100, 255], dtype=uint8)

    Overflow detection::

        >>> arr = np.array([0, 100, 300])
        >>> safe_np_int_conversion(arr, dtype=np.uint8)
        OverflowError: Cannot convert safely to uint8 type

    Signed to unsigned conversion (when values are non-negative)::

        >>> arr = np.array([-1, 0, 1], dtype=np.int32)
        >>> safe_np_int_conversion(arr[1:], dtype=np.uint8)
        array([0, 1], dtype=uint8)

    Notes
    -----
    Based on: https://stackoverflow.com/questions/56684893/
    Modified to handle signed→unsigned conversions with bounds checking.
    """
    int_array = np.array(int_array)
    if int_array.size == 0:
        return int_array.astype(dtype, copy=copy)  # Allow empty arrays of any type
    try:
        return int_array.astype(dtype, casting="safe", copy=copy)
    except TypeError:
        bounds = np.iinfo(dtype)
        if np.all(int_array >= bounds.min) and np.all(int_array <= bounds.max):
            if int_array.dtype.kind == "i" and np.dtype(dtype).kind == "u":
                # Allow casting from int to unsigned int, since we have checked bounds
                casting = "unsafe"
            else:
                # Raise a TypeError when we try to convert from, e.g., a float.
                casting = "same_kind"
            return int_array.astype(dtype, casting=casting, copy=copy)
        else:
            raise OverflowError("Cannot convert safely to {} type".format(dtype))


def safe_np_bool_conversion(int_array, copy=True):
    """Safely convert integer arrays to boolean with validation.

    NumPy's default boolean conversion treats any non-zero value as True,
    which can mask data errors. This function enforces that values must be
    exactly 0 or 1 before converting to boolean, catching unexpected values
    that might indicate bugs or data corruption.

    :param int_array: Input array or array-like of integers (must contain only 0 or 1)
    :type int_array: array_like
    :param copy: Whether to copy the data (True) or reuse memory if possible (False)
    :type copy: bool, optional
    :return: Array converted to np.bool dtype
    :rtype: numpy.ndarray
    :raises TypeError: If input cannot be safely converted to boolean
    :raises OverflowError: If any values are not 0 or 1

    Examples
    --------
    Valid boolean conversion::

        >>> arr = np.array([0, 1, 1, 0])
        >>> safe_np_bool_conversion(arr)
        array([False,  True,  True, False])

    Invalid values detected::

        >>> arr = np.array([0, 1, 2])
        >>> safe_np_bool_conversion(arr)
        OverflowError: Cannot convert safely to np.bool type

    Catching data errors::

        >>> # Simulation returned unexpected flag value
        >>> completion_flags = np.array([0, 1, -1, 1])
        >>> safe_np_bool_conversion(completion_flags)
        OverflowError: Cannot convert safely to np.bool type

    Notes
    -----
    Based on: https://stackoverflow.com/questions/56684893/
    Modified to validate boolean values before conversion.
    """
    int_array = np.array(int_array)
    if int_array.size == 0:
        return int_array.astype(np.bool, copy=copy)  # Allow empty arrays of any type
    try:
        return int_array.astype(np.bool, casting="safe", copy=copy)
    except TypeError:
        if np.all(int_array >= 0) and np.all(int_array <= 1):
            if int_array.dtype.kind == "i":
                # Allow casting from int to bool, since we have checked bounds
                casting = "unsafe"
            else:
                # Raise a TypeError when we try to convert from, e.g., a float.
                casting = "same_kind"
            return int_array.astype(np.bool, casting=casting, copy=copy)
        else:
            raise OverflowError("Cannot convert safely to np.bool type")


def generate_macroscale_in(in_data: DataCollectionType) -> DataCollectionType:
    """Generate macroscale input data from microscale simulation output.

    This function transforms microscale simulation results into the input format
    required by the Fortran macroscale model. It performs several key operations:

    1. Bins microscale simulations (typically 50,000) into 100 groups based on
       tPA leaving time distribution
    2. Computes the cumulative distribution function (CDF) of tPA leaving times
    3. Organizes fiber degradation times by bin and simulation outcome
    4. Calculates forced vs kinetic unbinding proportions
    5. Generates grid neighborhood structure for spatial modeling

    **Important**: While this function uses the HDF5/v2.0.0 data specification
    internally, the output is intended for the Fortran macroscale model and
    **must be converted to v1.99.0 specification before writing to disk**.

    :param in_data: Microscale simulation output containing results from many
                    individual fiber simulations. Must include fields like
                    pli_first_time, tpa_leaving_time, fiber_degraded, etc.
    :type in_data: DataCollectionType
    :return: Macroscale input data with binned statistics and grid structure.
             Fields include bin_edge_proportions, bin_edge_tpa_leaving_time,
             binned_fiber_degrade_time, binned_fiber_degraded, edge_grid_neighbors.
    :rtype: DataCollectionType

    Examples
    --------
    Typical workflow::

        # Load microscale results
        micro_data = load_microscale_output("microscale_results.h5")

        # Generate macroscale input
        macro_input = generate_macroscale_in(micro_data)

        # Convert to Fortran format and write
        fortran_data = convert_data(macro_input, "v2.0.0", "v1.99.0")
        write_fortran_files(fortran_data)

    Notes
    -----
    - Microscale simulations will always be divided into 100 bins
    - Uses 6,000 seconds as "infinity" marker for incomplete fiber degradation
    - Preserves all parameters from input data in the output

    See Also
    --------
    :func:`convert_data` : Convert between specification versions
    :func:`generate_fortran_neighborhood_structure` : Create grid neighbor lists
    """

    # TODO: Add code to check that `in_data` meets the "current" data specification
    out_data = in_data.copy()

    # Determine bin size: divide microscale simulations into 100 bins
    # Typically 50,000 microscale runs → 500 simulations per bin
    bin_size = in_data["pli_first_time"].size // 100

    # Create the CDF of tPA leaving time distribution
    # bin_edge_proportions = [0.00, 0.01, 0.02, ..., 0.99, 1.00] (101 values)
    # These represent cumulative proportions: 0%, 1%, 2%, ..., 100% of tPA have left
    out_data["bin_edge_proportions"] = np.append(np.arange(0, 1, 0.01), [1.0])

    # Sort microscale simulations by tPA leaving time
    # This ordering will be used to assign simulations to bins
    indices = in_data["tpa_leaving_time"].argsort()

    # Find the actual tPA leaving times at bin boundaries
    # E.g., bin 0 ends at time T₁, bin 1 ends at time T₂, etc.
    # indices[bin_size-1::bin_size] selects the last simulation in each bin
    out_data["bin_edge_tpa_leaving_time"] = np.append(
        [0], in_data["tpa_leaving_time"][indices[bin_size - 1 :: bin_size]]
    )

    # Prepare fiber degradation time data for binning
    # Start with boolean array indicating which simulations completed
    lysis_complete = in_data["fiber_degraded"]
    # Get the final simulation time for each microscale run
    lysis_time = in_data["sim_final_time"]

    # For simulations where lysis did NOT complete, mark with sentinel value
    # 6,000 seconds represents "infinity" (much longer than typical simulation times)
    # This allows incomplete simulations to be identified and handled separately
    lysis_time[~lysis_complete] = 6_000

    # Organize fiber degradation times into 100 bins (columns)
    # 1. Split simulations into bins using the sorted ordering (indices)
    # 2. Sort lysis times within each bin (fastest to slowest degradation)
    # 3. Stack bins as rows, then transpose so bins become columns
    # Result: each column = one bin, sorted by degradation time
    out_data["binned_fiber_degrade_time"] = np.stack(
        [
            np.sort(lysis_time[indices[i * bin_size : (i + 1) * bin_size]])
            for i in range(100)
        ]
    ).T

    # Count how many simulations in each bin successfully degraded
    # argmax() finds first occurrence of 6000 (sentinel) in each column
    # If no 6000 exists (all degraded), argmax returns 0 → all succeeded
    # TODO: Handle case where all fibers degrade
    # This gives the count of successful degradations per bin
    out_data["binned_fiber_degraded"] = out_data["binned_fiber_degrade_time"].argmax(
        axis=0
    )

    # Generate spatial grid neighborhood structure for macroscale model
    # Creates list of neighbors for each grid location (edge) in the 2D grid
    out_data["edge_grid_neighbors"] = generate_fortran_neighborhood_structure(
        in_data["params"]["macro_params"]["rows"],
        in_data["params"]["macro_params"]["cols"],
    )

    # Calculate the proportion of tPA unbinding events that were forced (by PLi)
    # vs kinetic (spontaneous). This ratio affects macroscale binding dynamics.
    # Formula: forced_unbind_rate = n_forced / (n_forced + n_kinetic)
    out_data["params"]["macro_params"]["forced_unbind"] = np.count_nonzero(
        in_data["tPA_forced_unbind"]
    ) / (
        np.count_nonzero(in_data["tPA_forced_unbind"])
        + np.count_nonzero(in_data["tPA_kinetic_unbind"])
    )

    return out_data


def convert_fiber_degrade_time(
    input_data: DataCollectionType, input_spec: str, output_spec: str
) -> DataSetType:
    """Convert fiber degradation time data between specification formats.

    This function handles the transformation of fiber degradation event data
    from Fortran's format (using 1D grid indices) to HDF5 format (using 2D
    row/column coordinates). The data tracks when and where fibers degrade
    during macroscale simulation.

    The Fortran format uses a single "Grid Location Index" that must be
    converted to separate "Grid Location Row" and "Grid Location Rank" fields
    in the HDF5 format. The conversion accounts for Fortran's 1-based indexing.

    :param input_data: Complete data collection containing fiber degradation
                       events and grid parameters
    :type input_data: DataCollectionType
    :param input_spec: Input specification version (e.g., "v1.99.0")
    :type input_spec: str
    :param output_spec: Output specification version (e.g., "v2.0.0")
    :type output_spec: str
    :return: List of structured arrays containing converted degradation events,
             one array per macroscale simulation
    :rtype: list[numpy.ndarray]

    Examples
    --------
    Converting Fortran macroscale output to HDF5 format::

        # Load Fortran macroscale output
        fortran_data = load_fortran_macro_output()

        # Convert fiber degradation time data
        hdf5_degrade_times = convert_fiber_degrade_time(
            fortran_data, "v1.99.0", "v2.0.0"
        )

    Notes
    -----
    - Assumes that input data is in Fortran format while output is in HDF5 format
    - Converts from 1-based (Fortran) to 0-based (Python) indexing
    - Splits 1D grid indices into 2D row/column coordinates
    - Preserves timing information unchanged
    - TODO: Check if simulations_combined = True and handle accordingly
    - TODO: Implement for other data specifications.

    See Also
    --------
    :func:`from_fortran_edge_index_array` : Converts 1D grid indices to 2D coordinates
    """
    # TODO: Check if simulations_combined = True
    output_data = []

    # Process each macroscale simulation's degradation events
    for data in input_data["f_deg_list"]:
        # Create empty array with correct dtype for output specification
        output_data.append(
            np.empty(
                data.shape,
                dtype=dataspec["v2.0.0"]["macroscale_out"]
                .data["fiber_degrade_time"]
                .dtype,
            )
        )

        # Copy timing fields directly (no conversion needed)
        output_data[-1][["Simulation Time Elapsed", "Fiber New Degrade Time"]] = data[
            ["Simulation Time Elapsed", "Fiber New Degrade Time"]
        ]

        # Convert 1D Fortran grid index to 2D (row, column) coordinates
        # Subtract 1 to convert from Fortran's 1-based to Python's 0-based indexing
        locations = from_fortran_edge_index_array(
            data["Grid Location Index"] - 1,
            input_data["params"]["macro_params"]["rows"],
            input_data["params"]["macro_params"]["cols"],
        )

        # Assign row and column coordinates to output fields
        output_data[-1]["Grid Location Row"] = locations[:, 0]
        output_data[-1]["Grid Location Rank"] = locations[:, 1]

    return output_data


# TODO Do the same thing with the tpa_bind_events from cell 10 of H5-File-Builder.ipynb


def _not_implemented(dataset_name: str, input_spec: str, output_spec: str):
    """Create a converter stub that raises NotImplementedError with a clear message.

    :param dataset_name: Name of the dataset being converted
    :type dataset_name: str
    :param input_spec: Input specification version
    :type input_spec: str
    :param output_spec: Output specification version
    :type output_spec: str
    :return: Function that raises NotImplementedError when called
    :rtype: Callable
    """

    def _raise_error(data):
        raise NotImplementedError(
            f"Conversion for dataset '{dataset_name}' from {input_spec} to {output_spec} "
            f"is not yet implemented."
        )

    return _raise_error


# Dictionary mapping (input_spec, output_spec) pairs to dataset converter functions
# Structure: {(input_version, output_version): {dataset_name: converter_function}}
#
# Converter functions take a DataCollectionType and return a DataSetType (or list thereof)
# They can be:
#   - Lambda functions for direct field mapping: lambda data: data["field_name"]
#   - Partial functions with preset parameters: functools.partial(func, param="value")
#   - Full converter functions: convert_fiber_degrade_time
#   - Not-implemented stubs: _not_implemented("name", "in_spec", "out_spec")
data_converters: dict[
    tuple[str, str], dict[str, Callable[[DataCollectionType], DataSetType]]
] = {
    # Convert from v2.0.0 (HDF5 unified format) to v1.99.0 (Fortran file-based format)
    ("v2.0.0", "v1.99.0"): {
        # Microscale output datasets (fully implemented - direct field mapping)
        "micro_log": lambda data: data["micro_log"],
        "firstPLi": lambda data: data["pli_first_time"],
        "lasttPA": lambda data: data["tpa_final_num"],
        "lyscomplete": lambda data: data["fiber_degraded"],
        "lysis": lambda data: data["sim_final_time"],
        "PLi": lambda data: data["pli_generated_num"],
        "tPA_time": lambda data: data["tpa_leaving_time"],
        "tPAPLiunbd": lambda data: data["tpa_unbound_by_pli"],
        "tPAunbind": lambda data: data["tpa_unbound_kinetic"],
        # Macroscale input datasets (not yet implemented - complex binned data structures)
        "tPAleave": _not_implemented("tPAleave", "v2.0.0", "v1.99.0"),
        "tsectPA": _not_implemented("tsectPA", "v2.0.0", "v1.99.0"),
        "lysismat": _not_implemented("lysismat", "v2.0.0", "v1.99.0"),
        "lenlysisvect": _not_implemented("lenlysisvect", "v2.0.0", "v1.99.0"),
        "neighbors": _not_implemented("neighbors", "v2.0.0", "v1.99.0"),
        # Macroscale output datasets (partially implemented - requires grid index conversion)
        "macro_log": lambda data: data["macro_log"],  # Direct mapping
        "Nsave": _not_implemented("Nsave", "v2.0.0", "v1.99.0"),
        "tsave": _not_implemented("tsave", "v2.0.0", "v1.99.0"),
        "f_deg_list": _not_implemented("f_deg_list", "v2.0.0", "v1.99.0"),
        "m_bind_t": _not_implemented("m_bind_t", "v2.0.0", "v1.99.0"),
        "m_loc": _not_implemented("m_loc", "v2.0.0", "v1.99.0"),
        "m_bound": _not_implemented("m_bound", "v2.0.0", "v1.99.0"),
        "mfpt": _not_implemented("mfpt", "v2.0.0", "v1.99.0"),
    },
    # Convert from v1.99.0 (Fortran file-based format) to v2.0.0 (HDF5 unified format)
    ("v1.99.0", "v2.0.0"): {
        # Microscale output datasets (fully implemented - direct field mapping)
        "micro_log": lambda data: data["micro_log"],
        "pli_first_time": lambda data: data["firstPLi"],
        "tpa_final_num": lambda data: data["lasttPA"],
        "fiber_degraded": lambda data: data["lyscomplete"],
        "sim_final_time": lambda data: data["lysis"],
        "pli_generated_num": lambda data: data["PLi"],
        "tpa_leaving_time": lambda data: data["tPA_time"],
        "tpa_unbound_by_pli": lambda data: data["tPAPLiunbd"],
        "tpa_unbound_kinetic": lambda data: data["tPAunbind"],
        # Macroscale input datasets (not yet implemented - reverse binning is non-trivial)
        "bin_edge_proportions": _not_implemented(
            "bin_edge_proportions", "v1.99.0", "v2.0.0"
        ),
        "bin_edge_tpa_leaving_time": _not_implemented(
            "bin_edge_tpa_leaving_time", "v1.99.0", "v2.0.0"
        ),
        "binned_fiber_degrade_time": _not_implemented(
            "binned_fiber_degrade_time", "v1.99.0", "v2.0.0"
        ),
        "binned_fiber_degraded": _not_implemented(
            "binned_fiber_degraded", "v1.99.0", "v2.0.0"
        ),
        "edge_grid_neighbors": _not_implemented(
            "edge_grid_neighbors", "v1.99.0", "v2.0.0"
        ),
        # Macroscale output datasets (partially implemented - grid coordinate conversion)
        "macro_log": lambda data: data["macro_log"],  # Direct mapping
        "snapshot_time": _not_implemented("snapshot_time", "v1.99.0", "v2.0.0"),
        "fiber_degrade_time": functools.partial(
            convert_fiber_degrade_time, input_spec="v1.99.0", output_spec="v2.0.0"
        ),
        "tpa_bind_events": _not_implemented("tpa_bind_events", "v1.99.0", "v2.0.0"),
        "tpa_location_snapshot": _not_implemented(
            "tpa_location_snapshot", "v1.99.0", "v2.0.0"
        ),
        "tpa_transit_time": _not_implemented("tpa_transit_time", "v1.99.0", "v2.0.0"),
    },
}


def convert_data(
    input_data: DataCollectionType,
    input_set_spec: str,
    output_set_spec: str,
) -> DataCollectionType:
    """Convert a complete data collection between specification versions.

    This is the main entry point for data conversion. It converts an entire
    data collection from one specification format to another by applying
    appropriate converter functions to each dataset. The function:

    1. Resolves tag aliases (like "current", "fortran") to version numbers
    2. Iterates through all datasets required by the output specification
    3. Applies the appropriate converter function for each dataset
    4. Performs safe type conversions to match output specification dtypes
    5. Handles both combined and per-simulation data structures

    The function uses the ``data_converters`` registry to find appropriate
    conversion functions. Missing converters trigger warnings but don't stop
    the conversion process.

    :param input_data: Complete data collection to convert, including all
                       datasets and parameters
    :type input_data: DataCollectionType
    :param input_set_spec: Input specification version (e.g., "v1.99.0", "v2.0.0")
                           or tag alias (e.g., "fortran", "current")
    :type input_set_spec: str
    :param output_set_spec: Output specification version (e.g., "v1.99.0", "v2.0.0")
                            or tag alias (e.g., "fortran", "current")
    :type output_set_spec: str
    :return: Converted data collection meeting the output specification
    :rtype: DataCollectionType
    :raises NotImplementedError: If a required converter is not implemented
    :raises TypeError: If type conversion fails
    :raises OverflowError: If numerical values don't fit in target dtype

    Examples
    --------
    Convert Fortran format to HDF5 format::

        >>> fortran_data = load_fortran_files()
        >>> hdf5_data = convert_data(fortran_data, "v1.99.0", "v2.0.0")
        >>> write_hdf5(hdf5_data, "output.h5")

    Using tag aliases::

        >>> convert_data(data, "fortran", "current")  # fortran → latest version

    Convert microscale results for macroscale input::

        >>> micro_results = load_microscale_output()
        >>> macro_input = generate_macroscale_in(micro_results)
        >>> fortran_input = convert_data(macro_input, "v2.0.0", "v1.99.0")

    Notes
    -----
    - Tag aliases are automatically resolved before conversion
    - Parameters are copied directly without conversion
    - Missing datasets generate warnings but don't stop conversion
    - Type conversions use safe methods that check bounds
    - TODO: Validate input_data against input_set_spec before conversion
    - TODO: Validate that input_set_spec and output_set_spec exist

    See Also
    --------
    :func:`generate_macroscale_in` : Generate macroscale input from microscale output
    :func:`safe_np_int_conversion` : Safe integer type conversion
    :func:`safe_np_bool_conversion` : Safe boolean type conversion
    :data:`data_converters` : Registry of converter functions
    """
    # Resolve any tag aliases to actual version numbers
    # E.g., "fortran" → "v1.99.0", "current" → "v2.0.0"
    while input_set_spec in tags:
        input_set_spec = tags[input_set_spec]
    while output_set_spec in tags:
        output_set_spec = tags[output_set_spec]

    # TODO: Check that `input_set_spec` and `output_set_spec` exist
    # TODO: Add code to check that `input_data` meets the specifications of `input_set_spec`.

    # Initialize output dictionary and copy parameters (unchanged across formats)
    out_data = {}
    out_data["params"] = input_data["params"]

    # Iterate through all data collections in the output specification
    # (e.g., microscale_out, macroscale_in, macroscale_out)
    for collection in dataspec[output_set_spec].values():
        # Process each dataset within the collection
        for dataset_needed in collection.data.keys():
            # Apply the appropriate converter function for this dataset
            try:
                out = data_converters[input_set_spec, output_set_spec][dataset_needed](
                    input_data
                )
            except KeyError as e:
                # Converter not found - warn but continue with other datasets
                warnings.warn(f"Missing data for {dataset_needed}")
                continue

            # Normalize data structure: ensure we're working with a list
            # If simulations_combined=True, converter returns single array → wrap in list
            if collection.simulations_combined:
                out = [out]

            # Convert each data table to the correct dtype
            out_data[dataset_needed] = []
            for data_table in out:
                try:
                    # Attempt standard NumPy type conversion (same_kind casting)
                    output_table = data_table.astype(
                        collection.data[dataset_needed].dtype,
                        casting="same_kind",
                    )
                except TypeError as e:
                    # Standard conversion failed - use safe converters based on dtype kind
                    dt = np.dtype(collection.data[dataset_needed].dtype)
                    match dt.kind:
                        case "u":  # Unsigned integer - check bounds
                            output_table = safe_np_int_conversion(data_table, dtype=dt)
                        case "b":  # Boolean - validate values are 0 or 1
                            output_table = safe_np_bool_conversion(data_table)
                        case "f":  # Float - safe conversion not yet implemented
                            # TODO: Create a function that does the same thing as the safe_np_int_conversion, but for floats.
                            raise NotImplementedError("Not implemented yet")
                        case _:  # Other types - re-raise original error
                            raise e
                out_data[dataset_needed].append(output_table)

            # Unwrap list if simulations_combined=True (single array, not list of arrays)
            if collection.simulations_combined:
                out_data[dataset_needed] = out_data[dataset_needed][0]

    return out_data
