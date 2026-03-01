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

Converter functions can be:

- Simple lambdas for direct field mapping: ``lambda data: data["field_name"]``
- Partial functions with fixed parameters: ``functools.partial(converter, ...)``
- Generic converters for structured arrays: :func:`convert_structured_grid_fields`
- Generic converters for location snapshots: :func:`convert_location_snapshot`
- Not-implemented stubs: ``_not_implemented()`` for pending converters

Conversion Routing
------------------

Conversions are routed automatically via a spanning tree built from the
``data_converters`` registry.  At module load time the module:

1. Builds an undirected graph whose vertices are spec versions and whose
   edges are ``(a, b)`` pairs that have *both* directions in
   ``data_converters``.
2. Computes a BFS spanning tree rooted at ``min(versions)`` for
   determinism.
3. Precomputes every pairwise path and stores them in the public
   ``conversion_paths`` dict.

``convert_data()`` looks up the path and chains through intermediate
single-step conversions automatically.  For example, converting from
v1.95.0 to v2.0.0 chains through v1.99.0 without the caller needing to
know about intermediate versions.

To add a new spec version, register ``(new, existing)`` and
``(existing, new)`` converters in ``data_converters``; the spanning tree
and paths will be rebuilt on next import.

If the spanning-tree routing strategy proves insufficient (e.g., if a
future spec's grid indexing cannot be expressed as either 1D Fortran
indices or 2D row/rank coordinates), consider refactoring the grid
conversion helpers to use a normalizer/denormalizer registry, where each
spec version registers functions to convert its grid indices to/from a
canonical intermediate form.

Key Functions
-------------

- :func:`convert_data`: Main entry point for converting complete data collections
- :func:`generate_macroscale_in`: Generates macroscale input from microscale output
- :func:`convert_structured_grid_fields`: Generic converter for structured arrays
  with grid-location fields (e.g., fiber_degrade_time, tpa_bind_events)
- :func:`convert_location_snapshot`: Generic converter for molecule location
  snapshot arrays (m_loc / tpa_location_snapshot)
- :func:`safe_np_int_conversion`: Safely converts integer arrays with bounds checking
- :func:`safe_np_bool_conversion`: Safely converts boolean arrays with validation
- :func:`safe_np_string_conversion`: Safely converts string/object arrays between dtypes

Example Usage
-------------

Converting from v1.99.0 to v2.0.0 format::

    from lysis.dataio.dataconvert import convert_data

    # Load v1.99.0 data
    old_data = load_fortran_data()

    # Convert to v2.0.0 HDF5 format
    new_data = convert_data(old_data, "v1.99.0", "v2.0.0")

    # Write to HDF5 file
    write_hdf5(new_data, "output.h5")

Generating macroscale input from microscale output::

    from lysis.dataio.dataconvert import generate_macroscale_in

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

:mod:`lysis.dataio.dataspec` : Data specification definitions and validation
:mod:`lysis.geometry.edge_grid` : Grid neighborhood structure utilities
"""

import functools
import warnings
from collections import deque

from dataclasses import asdict
from enum import Flag, auto, unique
from typing import Any, AnyStr, List, Mapping, Union, Callable

import numpy as np
import h5py

from pint import Quantity

import dataclasses

from ..config.constants import CONST, Q_
from ..config.parameters import MacroParameters, MicroParameters, Parameters
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
from ..geometry.edge_grid import (
    generate_fortran_neighborhood_structure,
    from_fortran_edge_index_array,
    to_fortran_edge_index_array,
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

    This is particularly important when converting from signed to unsigned types,
    downcasting from larger to smaller integer types, or converting from float
    arrays that contain integer values (common when reading from text files).

    :param int_array: Input array or array-like of integers (or floats with integer values)
    :type int_array: array_like
    :param dtype: Target NumPy integer dtype (e.g., np.int32, np.uint8)
    :type dtype: numpy.dtype, optional
    :param copy: Whether to copy the data (True) or reuse memory if possible (False)
    :type copy: bool, optional
    :return: Array converted to target dtype
    :rtype: numpy.ndarray
    :raises TypeError: If input cannot be safely converted to target dtype (e.g., floats with non-integer values)
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

    Float array with integer values (e.g., from text file I/O)::

        >>> arr = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        >>> safe_np_int_conversion(arr, dtype=np.uint8)
        array([1, 2, 3], dtype=uint8)

    Notes
    -----
    Based on: https://stackoverflow.com/questions/56684893/
    Modified to handle:
    - Signed→unsigned conversions with bounds checking
    - Float→integer conversions when float values are whole numbers
    """
    int_array = np.array(int_array)
    if int_array.size == 0:
        return int_array.astype(dtype, copy=copy)  # Allow empty arrays of any type

    # Handle float arrays that contain integer values (e.g., from text file I/O)
    if int_array.dtype.kind == "f":
        # Check if all values are whole numbers
        if not np.all(np.equal(np.mod(int_array, 1), 0)):
            raise TypeError(
                f"Cannot convert float array with non-integer values to {dtype}. "
                f"Values must be whole numbers."
            )
        # Convert to int64 first to preserve values, then proceed with bounds checking
        int_array = int_array.astype(np.int64, copy=False)

    try:
        return int_array.astype(dtype, casting="safe", copy=copy)
    except TypeError:
        bounds = np.iinfo(dtype)
        if np.all(int_array >= bounds.min) and np.all(int_array <= bounds.max):
            if int_array.dtype.kind == "i" and np.dtype(dtype).kind == "u":
                # Allow casting from int to unsigned int, since we have checked bounds
                casting = "unsafe"
            else:
                # Raise a TypeError when we try to convert from other types
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
            if int_array.dtype.kind in ("i", "u"):
                # Allow casting from int/uint to bool, since we have checked bounds
                casting = "unsafe"
            else:
                # Raise a TypeError when we try to convert from, e.g., a float.
                casting = "same_kind"
            return int_array.astype(np.bool, casting=casting, copy=copy)
        else:
            raise OverflowError("Cannot convert safely to np.bool type")


def safe_np_string_conversion(str_array, dtype, copy=True):
    """Safely convert string/object arrays to a target string dtype.

    Handles conversions between object arrays containing strings and fixed-width
    Unicode string dtypes (e.g., '<U75'). This is needed when round-trip converting
    data between formats that use different string representations (e.g., v1.99.0
    uses fixed-width Unicode strings like '<U75', while v2.0.0 may use object dtype).

    The function allows 'unsafe' casting for string-to-string conversions because
    we're only changing the dtype representation, not the actual string data.

    :param str_array: Input array or array-like of strings
    :type str_array: array_like
    :param dtype: Target NumPy dtype (e.g., 'U75', np.dtype('<U75'), object)
    :type dtype: numpy.dtype or str
    :param copy: Whether to copy the data (True) or reuse memory if possible (False)
    :type copy: bool, optional
    :return: Array converted to target dtype
    :rtype: numpy.ndarray
    :raises TypeError: If input cannot be safely converted to target dtype

    Examples
    --------
    Object to Unicode string (common in round-trip conversions)::

        >>> arr = np.array(['hello', 'world'], dtype=object)
        >>> safe_np_string_conversion(arr, dtype='U10')
        array(['hello', 'world'], dtype='<U10')

    Unicode to object::

        >>> arr = np.array(['hello', 'world'], dtype='U5')
        >>> safe_np_string_conversion(arr, dtype=object)
        array(['hello', 'world'], dtype=object)

    Unicode size change::

        >>> arr = np.array(['test'], dtype='U4')
        >>> safe_np_string_conversion(arr, dtype='U10')
        array(['test'], dtype='<U10')

    Notes
    -----
    This function is particularly important for round-trip conversions where:
    - v1.99.0 → v2.0.0: Fixed-width Unicode ('<U75') may become object dtype
    - v2.0.0 → v1.99.0: Object dtype needs to convert back to fixed-width Unicode

    The 'unsafe' casting is safe for string types because NumPy will:
    - Truncate strings that are too long for the target dtype (by design)
    - Preserve all data when converting object → Unicode if strings fit
    """
    str_array = np.array(str_array)
    target_dtype = np.dtype(dtype)

    if str_array.size == 0:
        return str_array.astype(target_dtype, copy=copy)

    # For string/object conversions, we can use 'unsafe' casting
    # because we're not changing the actual data, just the dtype representation
    # Check if source and target are both string-like types
    if str_array.dtype.kind in ["U", "S", "O"] and target_dtype.kind in ["U", "S", "O"]:
        return str_array.astype(target_dtype, casting="unsafe", copy=copy)
    else:
        # Not a string-to-string conversion, raise TypeError
        raise TypeError(
            f"Cannot convert from dtype {str_array.dtype} to {target_dtype}"
        )


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
    - Uses infinity as marker for incomplete fiber degradation
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
    # This allows incomplete simulations to be identified and handled separately
    lysis_time[~lysis_complete] = float("inf")

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
    # argmax() finds first occurrence of infinity (sentinel) in each column
    # If no infinity exists (all degraded), argmax returns 0 → all succeeded
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
        in_data["tpa_unbound_by_pli"]
    ) / (
        np.count_nonzero(in_data["tpa_unbound_by_pli"])
        + np.count_nonzero(in_data["tpa_unbound_kinetic"])
    )

    return out_data


def convert_structured_grid_fields(
    input_data: DataCollectionType,
    input_dataset: str,
    output_dataset: str,
    input_spec: str,
    output_spec: str,
    field_offsets: dict[str, int] | None = None,
) -> list[np.ndarray]:
    """Convert structured arrays between grid indexing systems.

    Generic converter for structured arrays that differ primarily in their
    grid-location fields. Copies fields with matching names directly and
    converts grid-location fields between indexing systems.

    Currently handles conversion between:

    - v1.99.0: 1D "Grid Location Index" (1-based Fortran indexing)
    - v2.0.0: 2D "Grid Location Row" + "Grid Location Rank" (0-based)

    This single function replaces the need for separate converters for each
    structured-array dataset (e.g., fiber_degrade_time, tpa_bind_events).

    :param input_data: Complete data collection
    :type input_data: DataCollectionType
    :param input_dataset: Key for the source dataset in input_data
    :type input_dataset: str
    :param output_dataset: Key for the dataset in the output dataspec
        (used to look up the output dtype)
    :type output_dataset: str
    :param input_spec: Input specification version (e.g., "v1.99.0")
    :type input_spec: str
    :param output_spec: Output specification version (e.g., "v2.0.0")
    :type output_spec: str
    :param field_offsets: Integer offsets to apply to common fields during
        copy. For example, ``{"tPA Molecule Index": -1}`` converts a
        1-based Fortran index to 0-based. Applied before dtype casting.
    :type field_offsets: dict[str, int] or None
    :return: List of structured arrays (one per simulation) with converted
        grid-location fields and matching fields copied directly
    :rtype: list[numpy.ndarray]
    :raises NotImplementedError: If either dtype contains fields that are
        neither common to both dtypes nor recognized grid-location fields
        ("Grid Location Index", "Grid Location Row", "Grid Location Rank").
        This prevents silent data loss when dtypes diverge unexpectedly.
    :raises NotImplementedError: If the input dtype contains no recognized
        grid-location fields.

    Notes
    -----
    Grid-location fields are identified by name: "Grid Location Index"
    (v1.99.0) and "Grid Location Row"/"Grid Location Rank" (v2.0.0).

    Conversion routing: all spec versions should convert through the
    v1.99.0 <-> v2.0.0 pair. Specs older than v1.99.0 should first convert
    to v1.99.0; specs newer than v2.0.0 should first convert from v2.0.0.
    If this routing strategy proves insufficient, consider refactoring to
    use a normalizer/denormalizer registry where each spec version registers
    how to convert its grid indices to/from a canonical form.

    See Also
    --------
    :func:`from_fortran_edge_index_array` : 1D Fortran indices to 2D coordinates
    :func:`to_fortran_edge_index_array` : 2D coordinates to 1D Fortran indices
    """
    rows = input_data["params"]["macro_params"]["rows"]
    cols = input_data["params"]["macro_params"]["cols"]
    out_dtype = dataspec[output_spec]["macroscale_out"].data[output_dataset].dtype

    # Known grid-location field names that this function can convert between
    grid_location_fields = {
        "Grid Location Index",
        "Grid Location Row",
        "Grid Location Rank",
    }

    output = []
    for table in input_data[input_dataset]:
        out_table = np.empty(table.shape, dtype=out_dtype)
        in_dtype = table.dtype

        # Copy fields that exist in both dtypes, applying offsets if specified
        common_fields = set(in_dtype.names) & set(out_dtype.names)
        for field in common_fields:
            if field_offsets and field in field_offsets:
                out_table[field] = table[field] + field_offsets[field]
            else:
                out_table[field] = table[field]

        # Check for unhandled fields: any field that is not common to both
        # dtypes and not a recognized grid-location field
        handled_fields = common_fields | grid_location_fields
        unhandled_in = set(in_dtype.names) - handled_fields
        unhandled_out = set(out_dtype.names) - handled_fields
        if unhandled_in or unhandled_out:
            msg = (
                f"Cannot automatically convert all fields between "
                f"{input_dataset} and {output_dataset}."
            )
            if unhandled_in:
                msg += f"\n  Unhandled input fields: {unhandled_in}"
            if unhandled_out:
                msg += f"\n  Unhandled output fields: {unhandled_out}"
            raise NotImplementedError(msg)

        # Convert grid-location fields between indexing systems
        if "Grid Location Index" in in_dtype.names:
            # 1D 1-based Fortran index → 2D 0-based (row, rank)
            coords = from_fortran_edge_index_array(
                table["Grid Location Index"] - 1, rows, cols
            )
            out_table["Grid Location Row"] = coords[:, 0]
            out_table["Grid Location Rank"] = coords[:, 1]
        elif "Grid Location Row" in in_dtype.names and "Grid Location Rank" in in_dtype.names:
            # 2D 0-based (row, rank) → 1D 1-based Fortran index
            coords = np.column_stack(
                [
                    table["Grid Location Row"],
                    table["Grid Location Rank"],
                ]
            )
            out_table["Grid Location Index"] = (
                to_fortran_edge_index_array(coords, rows, cols) + 1
            )
        else:
            raise NotImplementedError(
                "Only route v1.99.0 and v2.0.0 specifications through this function. " \
                "Other specifications should be converted into one of those first."
            )

        output.append(out_table)
    return output


def convert_location_snapshot(
    input_data: DataCollectionType,
    input_dataset: str,
    input_spec: str,
    output_spec: str,
) -> list[np.ndarray]:
    """Convert tPA location snapshot arrays between grid indexing systems.

    Handles conversion of per-timestep molecule location arrays between
    v1.99.0's flat 1D index format and v2.0.0's 3D coordinate format.
    The two formats also differ in axis ordering.

    - v1.99.0 ``m_loc``: shape ``(n_snapshots, n_molecules)`` of 1-based
      1D Fortran edge indices
    - v2.0.0 ``tpa_location_snapshot``: shape ``(n_molecules, 2, n_snapshots)``
      of 0-based 2D (row, rank) coordinates

    :param input_data: Complete data collection
    :type input_data: DataCollectionType
    :param input_dataset: Key for the source dataset in input_data
    :type input_dataset: str
    :param input_spec: Input specification version — must be "v1.99.0" or "v2.0.0"
    :type input_spec: str
    :param output_spec: Output specification version (e.g., "v2.0.0")
    :type output_spec: str
    :return: List of converted location arrays (one per simulation)
    :rtype: list[numpy.ndarray]
    :raises NotImplementedError: If input_spec is not "v1.99.0" or "v2.0.0".
        Other specifications should be converted to one of these first.

    Notes
    -----
    Only v1.99.0 and v2.0.0 are supported directly. Other spec versions
    should be converted to one of these before calling this function.
    See :func:`convert_structured_grid_fields` for the conversion routing
    strategy and future extensibility notes.

    See Also
    --------
    :func:`convert_structured_grid_fields` : For structured arrays with
        grid-location fields (e.g., fiber_degrade_time, tpa_bind_events)
    """
    rows = input_data["params"]["macro_params"]["rows"]
    cols = input_data["params"]["macro_params"]["cols"]
    n_mol = input_data["params"]["macro_params"]["total_molecules"]

    output = []
    for table in input_data[input_dataset]:
        if input_spec == "v1.99.0":
            # (n_snapshots, n_molecules) of 1D 1-based indices
            # → (n_molecules, 2, n_snapshots) of 2D 0-based coordinates
            coords = from_fortran_edge_index_array(
                (table - 1).ravel(),
                rows,
                cols,
            ).reshape(-1, n_mol, 2)
            # Rearrange axes: (snapshot, molecule, coord) → (molecule, coord, snapshot)
            output.append(np.moveaxis(coords, [0, 1, 2], [2, 0, 1]))
        elif input_spec == "v2.0.0":
            # (n_molecules, 2, n_snapshots) of 2D 0-based coordinates
            # → (n_snapshots, n_molecules) of 1D 1-based indices
            # Rearrange axes: (molecule, coord, snapshot) → (snapshot, molecule, coord)
            coords = np.moveaxis(table, [0, 1, 2], [1, 2, 0])
            output.append(
                to_fortran_edge_index_array(
                    coords.reshape(-1, 2),
                    rows,
                    cols,
                ).reshape(-1, n_mol)
                + 1
            )
        else:
            raise NotImplementedError(
                "Only route v1.99.0 and v2.0.0 specifications through this function. " \
                "Other specifications should be converted into one of those first."
            )
    return output


def replay_event_log_to_snapshot(
    event_log: np.ndarray,
    snapshot_times: np.ndarray,
    n_entities: int,
    entity_field: str,
    time_field: str,
    state_mapper: Callable[[np.ndarray], np.ndarray],
    initial_value: int = 0,
    output_dtype: np.dtype = np.int32,
) -> np.ndarray:
    """Replay a chronologically-sorted event log into a snapshot array.

    Given an event log where each row records a state change for some entity,
    produces a 2D array of shape ``(n_snapshots, n_entities)`` where cell
    ``[i, j]`` contains the state of entity ``j`` at snapshot time ``i``.

    The state of each entity at a given snapshot is determined by applying
    ``state_mapper`` to the entity's most recent event before that snapshot.
    Entities with no events retain ``initial_value``.

    Events are processed incrementally: for each snapshot, only the events
    that occurred since the previous snapshot are applied. This gives
    O(n_events + n_snapshots * n_entities) complexity.

    :param event_log: Structured array of events, sorted chronologically by
        ``time_field``.
    :type event_log: np.ndarray
    :param snapshot_times: 1D array of times at which to capture state.
    :type snapshot_times: np.ndarray
    :param n_entities: Total number of entities (determines output width).
    :type n_entities: int
    :param entity_field: Name of the structured array field containing
        0-based entity IDs.
    :type entity_field: str
    :param time_field: Name of the field containing event timestamps.
    :type time_field: str
    :param state_mapper: Function that accepts a slice of ``event_log``
        (structured array rows) and returns an array of output values.
        Called once per batch of new events at each snapshot.
    :type state_mapper: Callable[[np.ndarray], np.ndarray]
    :param initial_value: Value to fill the output array with before any
        events are applied. Defaults to 0.
    :type initial_value: int
    :param output_dtype: NumPy dtype for the output array.
    :type output_dtype: np.dtype
    :return: Array of shape ``(len(snapshot_times), n_entities)``.
    :rtype: np.ndarray
    """
    n_snapshots = len(snapshot_times)
    output = np.full((n_snapshots, n_entities), initial_value, dtype=output_dtype)

    if event_log.size == 0:
        return output

    event_times = event_log[time_field]
    entity_ids = event_log[entity_field]

    # For each snapshot, find the index of the first event AFTER that time
    cutoffs = np.searchsorted(event_times, snapshot_times, side="right")

    # Incrementally replay: maintain current state and apply only new events
    current_state = np.full(n_entities, initial_value, dtype=output_dtype)
    prev_cutoff = 0

    for i in range(n_snapshots):
        if cutoffs[i] > prev_cutoff:
            # Apply events that occurred between previous and current snapshot
            event_slice = event_log[prev_cutoff : cutoffs[i]]
            new_ids = entity_ids[prev_cutoff : cutoffs[i]]
            new_values = state_mapper(event_slice)
            # Scatter assignment: last write wins for duplicate entity IDs
            current_state[new_ids] = new_values
            prev_cutoff = cutoffs[i]

        output[i] = current_state

    return output


def convert_bind_events_to_bound(
    input_data: DataCollectionType,
    input_dataset: str = "tpa_bind_events",
    snapshot_dataset: str = "snapshot_time",
) -> list[np.ndarray]:
    """Convert tPA binding event logs to binding status snapshot arrays.

    Reconstructs the ``m_bound`` dataset by replaying ``tpa_bind_events``
    against ``snapshot_time``. For each simulation, produces an array of
    shape ``(n_snapshots, total_molecules)`` where each cell is 1 if the
    molecule is bound at that snapshot time, 0 otherwise.

    A molecule is considered bound if its most recent event before the
    snapshot had ``Molecule New Status == CONST.MOL_STATUS.BOUND``.
    All molecules start unbound.

    :param input_data: Complete data collection containing the event log,
        snapshot times, and parameters.
    :type input_data: DataCollectionType
    :param input_dataset: Key for the binding events dataset in
        ``input_data``.
    :type input_dataset: str
    :param snapshot_dataset: Key for the snapshot times dataset in
        ``input_data``.
    :type snapshot_dataset: str
    :return: List of arrays (one per simulation), each of shape
        ``(n_snapshots, total_molecules)`` with dtype ``np.int32``.
    :rtype: list[np.ndarray]
    """
    n_mol = input_data["params"]["macro_params"]["total_molecules"]

    def _is_bound(events: np.ndarray) -> np.ndarray:
        return (
            events["Molecule New Status"] == CONST.MOL_STATUS.BOUND
        ).astype(np.int32)

    output = []
    for events, snap_times in zip(
        input_data[input_dataset], input_data[snapshot_dataset]
    ):
        m_bound = replay_event_log_to_snapshot(
            event_log=events,
            snapshot_times=snap_times,
            n_entities=n_mol,
            entity_field="tPA Molecule Index",
            time_field="Simulation Time Elapsed",
            state_mapper=_is_bound,
            initial_value=0,
            output_dtype=np.int32,
        )
        output.append(m_bound)

    return output


def convert_f_deg_list_to_f_deg_time(input_data: DataCollectionType) -> list[np.ndarray]:
    """Convert v1.95.0 ``f_deg_list`` event log to v1.90.0 ``f_deg_time`` snapshot array.

    Replays the fiber degradation event log against snapshot times using
    :func:`replay_event_log_to_snapshot`. Because ``f_deg_list`` uses 1-based
    Fortran grid indices, the ``Grid Location Index`` column is decremented
    before passing to the replay function.

    :param input_data: Data collection containing ``f_deg_list``, ``tsave``,
        and ``params``.
    :type input_data: DataCollectionType
    :return: List of arrays (one per simulation), each of shape
        ``(n_snapshots, total_edges)`` with dtype ``np.float64``.  Entries
        for edges not yet scheduled for degradation are filled with
        ``9.9e100`` (the Fortran sentinel value).
    :rtype: list[np.ndarray]
    """
    n_edges = input_data["params"]["macro_params"]["total_edges"]
    initial_degrade_time = 9.9e100  # Fortran sentinel: t_degrade = 9.9d+100

    output = []
    for events, tsave in zip(input_data["f_deg_list"], input_data["tsave"]):
        # Adjust Grid Location Index: 1-based Fortran → 0-based Python
        zero_based = events.copy()
        zero_based["Grid Location Index"] -= 1
        f_deg_time = replay_event_log_to_snapshot(
            event_log=zero_based,
            snapshot_times=tsave,
            n_entities=n_edges,
            entity_field="Grid Location Index",
            time_field="Simulation Time Elapsed",
            state_mapper=lambda e: e["Fiber New Degrade Time"],
            initial_value=initial_degrade_time,
            output_dtype=np.float64,
        )
        output.append(f_deg_time)
    return output


def convert_f_deg_time_to_f_deg_list(input_data: DataCollectionType) -> list[np.ndarray]:
    """Convert v1.90.0 ``f_deg_time`` snapshot array to v1.95.0 ``f_deg_list`` event log.

    Detects fiber degradation scheduling changes between consecutive snapshots
    and emits one event per change.  This conversion is **lossy**: exact event
    timestamps are replaced by the snapshot time at which the change was first
    observed.

    :param input_data: Data collection containing ``f_deg_time``, ``tsave``,
        and ``params``.
    :type input_data: DataCollectionType
    :return: List of structured arrays (one per simulation), each with fields
        ``("Simulation Time Elapsed", "Grid Location Index",
        "Fiber New Degrade Time")`` matching the v1.95.0 dtype.
    :rtype: list[np.ndarray]
    """
    initial_degrade_time = 9.9e100
    output = []
    for f_deg_time, tsave in zip(input_data["f_deg_time"], input_data["tsave"]):
        events = []
        n_snapshots, n_edges = f_deg_time.shape
        # Initial snapshot: edges already scheduled (not at sentinel)
        for idx in np.where(f_deg_time[0] < initial_degrade_time)[0]:
            events.append((tsave[0], idx + 1, f_deg_time[0, idx]))  # 1-based index
        # Subsequent snapshots: detect changes from previous snapshot
        for t_idx in range(1, n_snapshots):
            for idx in np.where(f_deg_time[t_idx] != f_deg_time[t_idx - 1])[0]:
                events.append((tsave[t_idx], idx + 1, f_deg_time[t_idx, idx]))
        f_deg_list_arr = np.array(
            events,
            dtype=np.dtype(
                [
                    ("Simulation Time Elapsed", np.float64),
                    ("Grid Location Index", np.int32),
                    ("Fiber New Degrade Time", np.float64),
                ]
            ),
        )
        output.append(f_deg_list_arr)
    return output


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


def _resolve_fortran_names(params):
    """Convert Fortran parameter names to Python names and apply transforms.

    Uses :meth:`~.parameters.Parameters.inverse_fortran_map` for the mapping.
    Keys that are already Python parameter names pass through unchanged
    (idempotent).  Unknown keys (not in inverse_map and not a Python name)
    also pass through.

    :param params: Parameters dict with ``{section_key: {name: value, ...}}``
        structure.  Non-dict sections pass through unchanged.
    :type params: dict
    :return: Parameters dict with Fortran names resolved to Python names.
    :rtype: dict
    """
    micro_inverse = MicroParameters.inverse_fortran_map()
    macro_inverse = MacroParameters.inverse_fortran_map(extra_cls=MicroParameters)
    micro_fields = {f.name for f in dataclasses.fields(MicroParameters)}
    macro_fields = {f.name for f in dataclasses.fields(MacroParameters)}
    all_fields = micro_fields | macro_fields

    out = {}
    for section_key, section in params.items():
        if not isinstance(section, dict):
            out[section_key] = section
            continue
        inverse = micro_inverse if "micro" in section_key else macro_inverse
        out[section_key] = {}
        for key, value in section.items():
            if key in all_fields:
                # Already a Python name — pass through unchanged
                out[section_key][key] = value
            elif key in inverse:
                py_name, transform, _ = inverse[key]
                out[section_key][py_name] = Parameters.apply_fortran_transform(
                    value, transform
                )
            else:
                # Unknown key — pass through for version-specific handling
                out[section_key][key] = value
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Legacy v1.90.0 parameter renames
# ─────────────────────────────────────────────────────────────────────────────

#: v1.90.0 micro log keys that differ from the v1.95.0 Python names.
_V190_MICRO_RENAMES = {"runs": "micro_simulations", "seed": "micro_seed"}
#: v1.90.0 micro log keys to silently discard.
_V190_MICRO_REMOVE = {"stats"}
#: v1.90.0 macro log keys that differ from the v1.95.0 Python names.
_V190_MACRO_RENAMES = {
    "total_trials": "macro_simulations",
    "log_lvl": "macro_log_lvl",
    "seed": "macro_seed",
}


def _convert_params_v190_to_v195(params):
    """Resolve Fortran names, rename v1.90.0 legacy keys, remove noise.

    :param params: Parameters dict with raw Fortran names.
    :type params: dict
    :return: Parameters dict with Python names (v1.95.0 convention).
    :rtype: dict
    """
    params = _resolve_fortran_names(params)
    out = {}
    for section_key, section in params.items():
        if not isinstance(section, dict):
            out[section_key] = section
            continue
        renames = _V190_MICRO_RENAMES if "micro" in section_key else _V190_MACRO_RENAMES
        removes = _V190_MICRO_REMOVE if "micro" in section_key else set()
        out[section_key] = {}
        for key, value in section.items():
            if key in removes:
                continue
            out[section_key][renames.get(key, key)] = value
    return out


def _convert_params_v195_to_v190(params):
    """Reverse rename Python names back to v1.90.0 legacy keys.

    :param params: Parameters dict with Python names (v1.95.0 convention).
    :type params: dict
    :return: Parameters dict with v1.90.0 keys where applicable.
    :rtype: dict
    """
    # Build reverse renames
    micro_reverse = {v: k for k, v in _V190_MICRO_RENAMES.items()}
    macro_reverse = {v: k for k, v in _V190_MACRO_RENAMES.items()}
    out = {}
    for section_key, section in params.items():
        if not isinstance(section, dict):
            out[section_key] = section
            continue
        reverse = micro_reverse if "micro" in section_key else macro_reverse
        out[section_key] = {}
        for key, value in section.items():
            out[section_key][reverse.get(key, key)] = value
    return out


def _convert_params_add_units(params: dict) -> dict:
    """Convert parameters from v1.95.0 (bare magnitudes) to v1.99.0 (with units).

    Resolves Fortran names first (for data read directly as v1.95.0), then
    wraps numeric values that have known units as Pint Quantity strings
    (e.g., ``0.000534`` → ``"0.000534 centimeter"``). Non-numeric values and
    parameters without known units pass through unchanged.

    :param params: Parameters dict with bare magnitude values
    :type params: dict
    :return: Parameters dict with dimensioned values as unit strings
    :rtype: dict
    """
    params = _resolve_fortran_names(params)
    units = Parameters.units()
    out = {}
    for section_key, section in params.items():
        if not isinstance(section, dict):
            out[section_key] = section
            continue
        out[section_key] = {}
        for key, value in section.items():
            if key in units and isinstance(value, (int, float)):
                out[section_key][key] = str(Q_(value, units[key]))
            else:
                out[section_key][key] = value
    return out


def _convert_params_strip_units(params: dict) -> dict:
    """Convert parameters from v1.99.0 (with units) to v1.95.0 (bare magnitudes).

    For each parameter section that is a dict, strips units from string values
    that have known units (from :meth:`Parameters.units`), converting them to
    bare magnitudes in the canonical unit (e.g., ``"0.000534 centimeter"`` →
    ``0.000534``). Non-string values and parameters without known units pass
    through unchanged.

    :param params: Parameters dict with dimensioned string values
    :type params: dict
    :return: Parameters dict with bare magnitude values
    :rtype: dict
    """
    units = Parameters.units()
    out = {}
    for section_key, section in params.items():
        if not isinstance(section, dict):
            out[section_key] = section
            continue
        out[section_key] = {}
        for key, value in section.items():
            if key in units and isinstance(value, str):
                out[section_key][key] = Q_(value).to(units[key]).magnitude
            else:
                out[section_key][key] = value
    return out


# Registry mapping (input_spec, output_spec) pairs to parameter converter functions.
# Used by convert_data() to transform parameter dictionaries between spec versions.
params_converters: dict[tuple[str, str], Callable] = {
    ("v1.95.0", "v1.99.0"): _convert_params_add_units,
    ("v1.99.0", "v1.95.0"): _convert_params_strip_units,
    ("v1.90.0", "v1.95.0"): _convert_params_v190_to_v195,
    ("v1.95.0", "v1.90.0"): _convert_params_v195_to_v190,
}


# Dictionary mapping (input_spec, output_spec) pairs to dataset converter functions
# Structure: {(input_version, output_version): {dataset_name: converter_function}}
#
# Converter functions take a DataCollectionType and return a DataSetType (or list thereof)
# They can be:
#   - Lambda functions for direct field mapping: lambda data: data["field_name"]
#   - Partial functions for grid-index conversion:
#       functools.partial(convert_structured_grid_fields, ...)
#       functools.partial(convert_location_snapshot, ...)
#   - Not-implemented stubs: _not_implemented("name", "in_spec", "out_spec")
#
# Routing: all conversions go through the v1.99.0 <-> v2.0.0 pair.
# To add an older spec (e.g., v1.0.0), add ("v1.0.0", "v1.99.0") and
# ("v1.99.0", "v1.0.0") entries. To add a newer spec (e.g., v3.0.0),
# add ("v2.0.0", "v3.0.0") and ("v3.0.0", "v2.0.0") entries.
data_converters: dict[
    tuple[str, str], dict[str, Callable[[DataCollectionType], DataSetType]]
] = {
    # ============================================================================
    # Convert from v1.90.0 to v1.95.0 (reconstruct event log from snapshots — lossy)
    # ============================================================================
    ("v1.90.0", "v1.95.0"): {
        name: (lambda data, n=name: data[n])
        for collection in dataspec["v1.90.0"].values()
        for name in collection.data
        if name != "f_deg_time"
    }
    | {
        "f_deg_list": convert_f_deg_time_to_f_deg_list,
    },
    # ============================================================================
    # Convert from v1.95.0 to v1.90.0 (replay event log into snapshots)
    # ============================================================================
    ("v1.95.0", "v1.90.0"): {
        name: (lambda data, n=name: data[n])
        for collection in dataspec["v1.95.0"].values()
        for name in collection.data
        if name != "f_deg_list"
    }
    | {
        "f_deg_time": convert_f_deg_list_to_f_deg_time,
    },
    # ============================================================================
    # Convert from v1.95.0 to v1.99.0 (identity — data datasets are identical)
    # ============================================================================
    ("v1.95.0", "v1.99.0"): {
        name: (lambda data, n=name: data[n])
        for collection in dataspec["v1.95.0"].values()
        for name in collection.data
    },
    # ============================================================================
    # Convert from v1.99.0 to v1.95.0 (identity — data datasets are identical)
    # ============================================================================
    ("v1.99.0", "v1.95.0"): {
        name: (lambda data, n=name: data[n])
        for collection in dataspec["v1.99.0"].values()
        for name in collection.data
    },
    # ============================================================================
    # Convert from v2.0.0 (HDF5 unified format) to v1.99.0 (Fortran file-based format)
    # ============================================================================
    ("v2.0.0", "v1.99.0"): {
        # ------------------------------------------------------------------------
        # Microscale output datasets (fully implemented - direct field mapping)
        # ------------------------------------------------------------------------
        "micro_log": lambda data: data["micro_log"],  # Direct mapping
        "firstPLi": lambda data: data["pli_first_time"],  # Direct mapping
        "lasttPA": lambda data: data["tpa_final_num"],  # Direct mapping
        "lyscomplete": lambda data: data["fiber_degraded"],  # Direct mapping
        "lysis": lambda data: data["sim_final_time"],  # Direct mapping
        "PLi": lambda data: data["pli_generated_num"],  # Direct mapping
        "tPA_time": lambda data: data["tpa_leaving_time"],  # Direct mapping
        "tPAPLiunbd": lambda data: data["tpa_unbound_by_pli"],  # Direct mapping
        "tPAunbind": lambda data: data["tpa_unbound_kinetic"],  # Direct mapping
        # ------------------------------------------------------------------------
        # Macroscale input datasets (fully implemented)
        # See dataspec.py and docs/usage/data_specification.rst for field descriptions
        # ------------------------------------------------------------------------
        "tPAleave": lambda data: data[
            "bin_edge_proportions"
        ],  # Proportion of tPA leaving each edge per bin
        "tsectPA": lambda data: data[
            "bin_edge_tpa_leaving_time"
        ],  # tPA leaving times per bin
        "lysismat": lambda data: np.where(
            data["binned_fiber_degrade_time"]
            == float("inf"),  # v2.0.0 uses infinity sentinel
            6_000,  # v1.99.0/Fortran uses 6,000 sentinel for incomplete degradation
            data["binned_fiber_degrade_time"],
        ),
        "lenlysisvect": lambda data: data["binned_fiber_degraded"]
        + 1,  # Count of degraded fibers per bin (Fortran uses 1-based counting)
        "neighbors": lambda data: np.reshape(
            data["edge_grid_neighbors"] + 1,  # Grid neighbor indices (Fortran: 1-based)
            (-1, 1),  # Flatten to column vector for Fortran I/O
        ),
        # ------------------------------------------------------------------------
        # Macroscale output datasets
        # Uses generic converters: convert_structured_grid_fields() for
        # structured arrays and convert_location_snapshot() for location data
        # ------------------------------------------------------------------------
        "macro_log": lambda data: data["macro_log"],  # Direct mapping
        "Nsave": lambda data: [
            np.int32(len(x) - 1) for x in data["snapshot_time"]
        ],  # Fortran Nsave counts save intervals, excludes initial time
        "tsave": lambda data: data["snapshot_time"],  # Direct mapping
        "f_deg_list": functools.partial(
            convert_structured_grid_fields,
            input_dataset="fiber_degrade_time",
            output_dataset="f_deg_list",
            input_spec="v2.0.0",
            output_spec="v1.99.0",
        ),
        "m_bind_t": functools.partial(
            convert_structured_grid_fields,
            input_dataset="tpa_bind_events",
            output_dataset="m_bind_t",
            input_spec="v2.0.0",
            output_spec="v1.99.0",
            field_offsets={"tPA Molecule Index": 1},  # 0-based → 1-based Fortran
        ),
        "m_loc": functools.partial(
            convert_location_snapshot,
            input_dataset="tpa_location_snapshot",
            input_spec="v2.0.0",
            output_spec="v1.99.0",
        ),
        # m_bound: reconstructed from tpa_bind_events + snapshot_time
        # by replaying the event log (see convert_bind_events_to_bound)
        "m_bound": functools.partial(
            convert_bind_events_to_bound,
            input_dataset="tpa_bind_events",
            snapshot_dataset="snapshot_time",
        ),
        "mfpt": lambda data: data["tpa_transit_time"],  # Direct mapping
    },
    # ============================================================================
    # Convert from v1.99.0 (Fortran file-based format) to v2.0.0 (HDF5 unified format)
    # ============================================================================
    ("v1.99.0", "v2.0.0"): {
        # ------------------------------------------------------------------------
        # Microscale output datasets (fully implemented - direct field mapping)
        # ------------------------------------------------------------------------
        "micro_log": lambda data: data["micro_log"],  # Direct mapping
        "pli_first_time": lambda data: data["firstPLi"],  # Direct mapping
        "tpa_final_num": lambda data: data["lasttPA"],  # Direct mapping
        "fiber_degraded": lambda data: data["lyscomplete"],  # Direct mapping
        "sim_final_time": lambda data: data["lysis"],  # Direct mapping
        "pli_generated_num": lambda data: data["PLi"],  # Direct mapping
        "tpa_leaving_time": lambda data: data["tPA_time"],  # Direct mapping
        "tpa_unbound_by_pli": lambda data: data["tPAPLiunbd"],  # Direct mapping
        "tpa_unbound_kinetic": lambda data: data["tPAunbind"],  # Direct mapping
        # ------------------------------------------------------------------------
        # Macroscale input datasets (For testing/round-trip conversion only!)
        # Generated by generate_macroscale_in() - don't store redundant copies in HDF5
        # See dataspec.py and docs/usage/data_specification.rst for field descriptions
        # ------------------------------------------------------------------------
        "bin_edge_proportions": lambda data: data[
            "tPAleave"
        ],  # Proportion of tPA leaving each edge per bin
        "bin_edge_tpa_leaving_time": lambda data: data[
            "tsectPA"
        ],  # tPA leaving times per bin
        "binned_fiber_degrade_time": lambda data: np.where(
            data["lysismat"] == 6_000,  # v1.99.0/Fortran sentinel
            float("inf"),  # v2.0.0 uses infinity for incomplete degradation
            data["lysismat"],
        ),
        "binned_fiber_degraded": lambda data: data["lenlysisvect"]
        - 1,  # Count of degraded fibers per bin (convert Fortran 1-based to 0-based)
        "edge_grid_neighbors": lambda data: np.reshape(
            data["neighbors"]
            - 1,  # Grid neighbor indices (convert Fortran 1-based to 0-based)
            (
                -1,
                8,
            ),  # Unflatten from column vector: one row per edge, 8 neighbors per row
        ),
        # ------------------------------------------------------------------------
        # Macroscale output datasets
        # Uses generic converters: convert_structured_grid_fields() for
        # structured arrays and convert_location_snapshot() for location data
        # ------------------------------------------------------------------------
        "macro_log": lambda data: data["macro_log"],  # Direct mapping
        "snapshot_time": lambda data: data["tsave"],  # Direct mapping
        "fiber_degrade_time": functools.partial(
            convert_structured_grid_fields,
            input_dataset="f_deg_list",
            output_dataset="fiber_degrade_time",
            input_spec="v1.99.0",
            output_spec="v2.0.0",
        ),
        "tpa_bind_events": functools.partial(
            convert_structured_grid_fields,
            input_dataset="m_bind_t",
            output_dataset="tpa_bind_events",
            input_spec="v1.99.0",
            output_spec="v2.0.0",
            field_offsets={"tPA Molecule Index": -1},  # 1-based Fortran → 0-based
        ),
        "tpa_location_snapshot": functools.partial(
            convert_location_snapshot,
            input_dataset="m_loc",
            input_spec="v1.99.0",
            output_spec="v2.0.0",
        ),
        "tpa_transit_time": lambda data: data["mfpt"],  # Direct mapping
    },
}


def _build_conversion_graph() -> dict[str, set[str]]:
    """Build an undirected adjacency dict from ``data_converters`` keys.

    An undirected edge ``(a, b)`` exists only when **both** ``(a, b)`` and
    ``(b, a)`` are present in ``data_converters``.  If only one direction
    exists, a :func:`warnings.warn` is emitted.

    :return: Adjacency dict mapping each version to its set of neighbours.
    :rtype: dict[str, set[str]]
    """
    directed_edges = set(data_converters.keys())
    graph: dict[str, set[str]] = {}

    # Collect all version strings that appear
    for a, b in directed_edges:
        graph.setdefault(a, set())
        graph.setdefault(b, set())

    for a, b in directed_edges:
        if (b, a) in directed_edges:
            graph[a].add(b)
            graph[b].add(a)
        else:
            warnings.warn(
                f"data_converters has ({a!r}, {b!r}) but not the reverse "
                f"({b!r}, {a!r}). This edge will be excluded from the "
                f"conversion graph.",
                stacklevel=2,
            )

    return graph


def _build_spanning_tree(
    graph: dict[str, set[str]], root: str
) -> dict[str, str | None]:
    """BFS spanning tree from *root*.

    :param graph: Undirected adjacency dict.
    :type graph: dict[str, set[str]]
    :param root: Starting vertex (must be a key in *graph*).
    :type root: str
    :return: Parent-pointer dict ``{node: parent}``.  Root's parent is
        ``None``.
    :rtype: dict[str, str | None]
    """
    parent: dict[str, str | None] = {root: None}
    queue: deque[str] = deque([root])

    while queue:
        node = queue.popleft()
        for neighbour in sorted(graph[node]):
            if neighbour not in parent:
                parent[neighbour] = node
                queue.append(neighbour)

    return parent


def _extract_tree_path(
    parent: dict[str, str | None], src: str, dst: str
) -> list[str]:
    """Find the unique path between *src* and *dst* in a spanning tree.

    Uses the lowest-common-ancestor (LCA) approach: walk both nodes up
    to the root, then splice the two half-paths at the LCA.

    :param parent: Parent-pointer dict produced by :func:`_build_spanning_tree`.
    :type parent: dict[str, str | None]
    :param src: Start version.
    :type src: str
    :param dst: End version.
    :type dst: str
    :return: Ordered list of versions from *src* to *dst* (inclusive).
    :rtype: list[str]
    """
    # Ancestors of src (in root-first order)
    ancestors_src: list[str] = []
    node = src
    while node is not None:
        ancestors_src.append(node)
        node = parent[node]
    ancestors_src.reverse()  # root → ... → src

    ancestors_dst: list[str] = []
    node = dst
    while node is not None:
        ancestors_dst.append(node)
        node = parent[node]
    ancestors_dst.reverse()  # root → ... → dst

    # Find LCA (last common prefix element)
    lca_idx = 0
    for i in range(min(len(ancestors_src), len(ancestors_dst))):
        if ancestors_src[i] == ancestors_dst[i]:
            lca_idx = i
        else:
            break

    # Path = src → ... → LCA → ... → dst
    path_src_to_lca = list(reversed(ancestors_src[lca_idx:]))  # src → LCA
    path_lca_to_dst = ancestors_dst[lca_idx + 1 :]  # LCA children → dst
    return path_src_to_lca + path_lca_to_dst


def _build_conversion_paths() -> dict[tuple[str, str], list[str]]:
    """Build all pairwise conversion paths from the ``data_converters`` graph.

    Picks ``min(graph.keys())`` as the BFS root for determinism, builds
    a spanning tree, and precomputes all reachable ``(src, dst)`` paths.

    Warns about versions that are unreachable from the root (disconnected
    graph components).

    :return: Dict mapping ``(src_version, dst_version)`` to a list of
        intermediate version strings (inclusive of both endpoints).
    :rtype: dict[tuple[str, str], list[str]]
    """
    graph = _build_conversion_graph()

    if not graph:
        return {}

    root = min(graph.keys())
    parent = _build_spanning_tree(graph, root)

    # Warn about unreachable versions
    unreachable = set(graph.keys()) - set(parent.keys())
    if unreachable:
        warnings.warn(
            f"The following spec versions are not reachable from the "
            f"spanning-tree root {root!r} and cannot participate in "
            f"multi-step conversion: {sorted(unreachable)}",
            stacklevel=2,
        )

    reachable = sorted(parent.keys())
    paths: dict[tuple[str, str], list[str]] = {}
    for src in reachable:
        for dst in reachable:
            if src != dst:
                paths[src, dst] = _extract_tree_path(parent, src, dst)

    return paths


conversion_paths: dict[tuple[str, str], list[str]] = _build_conversion_paths()
"""Pre-computed conversion paths between all reachable spec version pairs.

Each key is a ``(src_version, dst_version)`` tuple and each value is the
ordered list of versions to traverse (inclusive), e.g.
``("v1.95.0", "v2.0.0")`` → ``["v1.95.0", "v1.99.0", "v2.0.0"]``.
"""


def _convert_single_step(
    input_data: DataCollectionType,
    input_set_spec: str,
    output_set_spec: str,
) -> DataCollectionType:
    """Perform a single-step conversion between adjacent specification versions.

    This handles one hop in the conversion graph: it applies the converter
    functions registered in ``data_converters[(input_set_spec, output_set_spec)]``
    and converts parameters via ``params_converters`` if an entry exists.

    The caller (``convert_data``) is responsible for tag resolution and for
    chaining multiple single-step conversions together.

    :param input_data: Complete data collection to convert, including all
                       datasets and parameters
    :type input_data: DataCollectionType
    :param input_set_spec: Input specification version (e.g., "v1.99.0").
        Must already be resolved (no tags).
    :type input_set_spec: str
    :param output_set_spec: Output specification version (e.g., "v2.0.0").
        Must already be resolved (no tags).
    :type output_set_spec: str
    :return: Converted data collection meeting the output specification
    :rtype: DataCollectionType
    :raises ValueError: If input data has no parameters
    :raises NotImplementedError: If a required converter is not implemented
    :raises TypeError: If type conversion fails
    :raises OverflowError: If numerical values don't fit in target dtype
    """
    # Validate that parameters are present — data must always include parameters
    if "params" not in input_data or not input_data["params"]:
        raise ValueError(
            "Input data must include parameters (data['params']). "
            "Data without parameters is not valid."
        )

    # Initialize output dictionary and convert/copy parameters
    out_data = {}
    if (input_set_spec, output_set_spec) in params_converters:
        out_data["params"] = params_converters[input_set_spec, output_set_spec](
            input_data["params"]
        )
    else:
        out_data["params"] = input_data["params"]

    # Iterate through all data collections in the output specification
    # (e.g., microscale_out, macroscale_in, macroscale_out)
    for name, collection in dataspec[output_set_spec].items():
        # Check if this collection exists in the input data
        # A collection is considered to exist if at least one of its datasets is present
        # This allows conversion to work with partial data (e.g., only microscale_out
        # without macroscale_in/out), preventing NotImplementedError for missing collections
        collection_exists = True
        for dataset in dataspec[input_set_spec][name].data.keys():
            if not dataset in input_data:
                collection_exists = False
                break  # No need to check further if any dataset is missing

        # Skip this collection if it doesn't exist in the input data
        # This is expected when converting data that only contains some collections
        # (e.g., microscale output only, without macroscale data)
        if not collection_exists:
            continue

        # Process each dataset within the collection
        for dataset_needed in collection.data.keys():
            # Apply the appropriate converter function for this dataset
            out = data_converters[input_set_spec, output_set_spec][dataset_needed](
                input_data
            )

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
                        case (
                            "U" | "S" | "O"
                        ):  # String types (Unicode, byte string, object)
                            output_table = safe_np_string_conversion(
                                data_table, dtype=dt
                            )
                        case "f":  # Float - safe conversion not yet implemented
                            # TODO: Create a function that does the same thing as the safe_np_int_conversion, but for floats.
                            raise NotImplementedError("Not implemented yet")
                        case _:  # Other types - re-raise original error
                            raise e
                out_data[dataset_needed].append(output_table)

            # Unwrap list if simulations_combined=True (single array, not list of arrays)
            if collection.simulations_combined:
                out_data[dataset_needed] = out_data[dataset_needed][0]

    # Inject flag when converting from v1.90.0 (f_deg_list timestamps are approximated)
    if input_set_spec == "v1.90.0" and output_set_spec == "v1.95.0":
        out_data["params"][CONST.APPROX_F_DEG_LIST_ATTR] = True
    # Propagate flag through subsequent conversion steps
    # (params_converters may rebuild the params dict and drop unknown keys)
    elif input_data.get("params", {}).get(CONST.APPROX_F_DEG_LIST_ATTR):
        out_data["params"][CONST.APPROX_F_DEG_LIST_ATTR] = True

    return out_data


def convert_data(
    input_data: DataCollectionType,
    input_set_spec: str,
    output_set_spec: str,
) -> DataCollectionType:
    """Convert a complete data collection between specification versions.

    This is the main entry point for data conversion. It converts an entire
    data collection from one specification format to another by applying
    appropriate converter functions to each dataset. Multi-step conversions
    are handled automatically by chaining through intermediate versions
    using the precomputed :data:`conversion_paths`. The function:

    1. Resolves tag aliases (like "current", "fortran") to version numbers
    2. Short-circuits if input and output specs are identical
    3. Looks up the conversion path in ``conversion_paths``
    4. Chains through consecutive single-step conversions along the path
    5. Each step: applies converters, performs safe type conversions, handles
       both combined and per-simulation data structures

    The function gracefully handles partial data by checking if each collection
    exists before attempting conversion. This allows converting only microscale
    data without requiring macroscale data to be present.

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
    :raises ValueError: If no conversion path exists between the two specs
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

    Multi-step conversion (automatically routed)::

        >>> convert_data(data, "v1.95.0", "v2.0.0")
        # Automatically chains: v1.95.0 → v1.99.0 → v2.0.0

    Convert microscale results for macroscale input::

        >>> micro_results = load_microscale_output()
        >>> macro_input = generate_macroscale_in(micro_results)
        >>> fortran_input = convert_data(macro_input, "v2.0.0", "v1.99.0")

    Partial data conversion (only microscale, no macroscale)::

        >>> # Load only microscale output (no macroscale data)
        >>> micro_only = read_data_collection(path, [microscale_out_spec], [""])
        >>> # Convert successfully - macroscale collections are skipped automatically
        >>> hdf5_data = convert_data(micro_only, "v1.99.0", "v2.0.0")

    Notes
    -----
    - Tag aliases are automatically resolved before conversion
    - If input and output specs resolve to the same version, input data is
      returned as-is
    - Parameters are converted at each step if a ``params_converters`` entry exists
    - **Collections that don't exist in input data are skipped** (no error raised)
    - Type conversions use safe methods that check bounds to prevent overflow
    - Supports partial data conversion (e.g., microscale only, without macroscale)
    - TODO: Validate input_data against input_set_spec before conversion
    - TODO: Validate that input_set_spec and output_set_spec exist

    See Also
    --------
    :func:`generate_macroscale_in` : Generate macroscale input from microscale output
    :func:`safe_np_int_conversion` : Safe integer type conversion
    :func:`safe_np_bool_conversion` : Safe boolean type conversion
    :data:`data_converters` : Registry of converter functions
    :data:`conversion_paths` : Precomputed multi-step conversion routes
    """
    # Resolve any tag aliases to actual version numbers
    # E.g., "fortran" → "v1.99.0", "current" → "v2.0.0"
    while input_set_spec in tags:
        input_set_spec = tags[input_set_spec]
    while output_set_spec in tags:
        output_set_spec = tags[output_set_spec]

    # TODO: Check that `input_set_spec` and `output_set_spec` exist
    # TODO: Add code to check that `input_data` meets the specifications of `input_set_spec`.

    # Short-circuit: if input and output specs are the same, return data as-is
    if input_set_spec == output_set_spec:
        return input_data

    # Look up the precomputed conversion path
    path = conversion_paths.get((input_set_spec, output_set_spec))
    if path is None:
        raise ValueError(
            f"No conversion path from {input_set_spec!r} to "
            f"{output_set_spec!r}. Available paths: "
            f"{sorted(conversion_paths.keys())}"
        )

    # Chain through consecutive pairs in the path
    data = input_data
    for step_in, step_out in zip(path[:-1], path[1:]):
        data = _convert_single_step(data, step_in, step_out)

    return data
