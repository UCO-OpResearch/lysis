"""File I/O operations for lysis simulation data.

This module provides a unified interface for reading and writing simulation data
across multiple storage formats. It abstracts away the complexity of different
file types, allowing the same high-level functions to work with Fortran text files,
binary files, JSON parameter files, and HDF5 datasets.

Storage Types
-------------

The module supports five storage types, each with specialized reader/writer functions:

1. **FILE_TEXT**: Plain text files with delimited data (e.g., CSV, space-separated)
   - Uses numpy.loadtxt/savetxt
   - Supports custom delimiters

2. **FILE_BINARY**: Binary files with raw numerical data
   - Uses numpy.fromfile
   - Requires shape information from parameters for correct reshaping

3. **FILE_JSON**: JSON files for parameter storage
   - Uses Python json library
   - Stores nested parameter dictionaries

4. **HDF5_DATASET**: HDF5 datasets for large numerical arrays
   - Uses h5py library
   - Supports compression, chunking, and dynamic shapes

5. **HDF5_ATTR**: HDF5 attributes for parameter metadata
   - Stores parameters as HDF5 group attributes
   - More efficient than separate parameter files

Architecture
------------

The module uses a registry-based design pattern:

- ``data_readers``: Dictionary mapping storage types to reader functions
- ``data_writers``: Dictionary mapping storage types to writer functions

Each reader/writer function follows a common signature, allowing the high-level
``read_dataset()`` and ``write_dataset()`` functions to dispatch to the
appropriate handler based on the storage type specified in the ``DataSetSpec``.

Key Functions
-------------

**High-level API** (recommended for most use cases):
  - :func:`read_data_collection`: Read complete data collections with all datasets
  - :func:`write_data_collection`: Write complete data collections with all datasets

**Mid-level API** (for individual datasets):
  - :func:`read_dataset`: Read a single dataset
  - :func:`write_dataset`: Write a single dataset

**Low-level API** (internal use):
  - ``_read_file_text``, ``_read_file_binary``, ``_read_file_json``
  - ``_read_hdf5_attr``, ``_read_hdf5_dataset``
  - ``_write_file_text``, ``_write_hdf5_attr``, ``_write_hdf5_dataset``

Example Usage
-------------

Reading a complete data collection::

    from lysis.dataio.fileops import read_data_collection
    from lysis.dataio.dataspec import dataspec

    # Read v1.99.0 Fortran microscale output
    collections = list(dataspec["v1.99.0"].values())
    data = read_data_collection(
        path="/path/to/fortran/output",
        collections=[collections[0]],  # microscale_out only
        file_codes=[""]
    )

    # Access datasets
    lysis_times = data["lysis"]
    parameters = data["params"]

Writing a complete data collection::

    from lysis.dataio.fileops import write_data_collection
    from lysis.dataio.dataspec import dataspec

    # Write v2.0.0 HDF5 format
    collections = list(dataspec["v2.0.0"].values())
    write_data_collection(
        data=simulation_results,
        path="/path/to/output.h5",
        collections=collections,
        file_codes=[""]
    )

Reading/writing individual datasets::

    from lysis.dataio.fileops import read_dataset, write_dataset

    # Read a single dataset
    spec = dataspec["v1.99.0"]["microscale_out"].data["lysis"]
    lysis_data = read_dataset(
        path="/path/to/data",
        spec=spec,
        params=parameters,
        sim=0  # First simulation
    )

    # Write a single dataset
    write_dataset(
        data=lysis_data,
        path="/path/to/output",
        spec=spec,
        params=parameters,
        sim=0
    )

File Codes and Simulation Indexing
-----------------------------------

**file_code**: Optional string inserted into filenames for organizing outputs
  - Example: ``file_code="run1"`` → ``microscale_out_run1.h5``
  - Useful for batch processing or parameter sweeps

**sim**: Simulation index for per-simulation storage
  - When ``simulations_combined=False``, each simulation is stored separately
  - sim=0, 1, 2, ... indexes individual simulation files
  - When ``simulations_combined=True``, all simulations are in one file (sim is ignored)

Notes
-----

- All write operations validate data against specifications using ``check_dataset_spec()``
- HDF5 writes use gzip compression by default
- Binary file reading requires parameter-based shape resolution
- Missing files during multi-simulation reads are handled gracefully (stop iteration)
- File paths are constructed from ``spec.data_location`` with format string substitution

See Also
--------

:mod:`lysis.dataio.dataspec` : Data specification definitions
:mod:`lysis.dataio.dataconvert` : Data format conversion utilities
:mod:`lysis.config.constants` : Storage type constants
"""

import json
import os
import re
import warnings
from pathlib import Path

from typing import AnyStr, Callable

import numpy as np
import h5py

from ..config.constants import CONST, Q_
from pint import Quantity

from .dataspec import (
    DataCollectionSpec,
    DataSetSpec,
    DataCollectionType,
    BaseParamsType,
    CONVERSION_WARNINGS,
    parse_shape,
    check_dataset_spec,
)

__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


# ─────────────────────────────────────────────────────────────────────────────
# Fortran log parsing
# ─────────────────────────────────────────────────────────────────────────────


def _parse_fortran_kv(lines):
    """Extract ``key = value`` pairs from Fortran log lines.

    Matches two formats:

    * ``key = value`` -- standard Fortran output.
    * ``Setting key = value`` -- command-line argument echoes written by the
      Fortran code during startup.

    Lines that match neither pattern are silently skipped.  When the same key
    appears in both a ``Setting`` line and a later ``key = value`` line, the
    later value overwrites the earlier one.

    :param lines: Iterable of text lines from a Fortran log file.
    :return: ``{fortran_name_lower: raw_value_str}``.
    :rtype: dict[str, str]
    """
    kv_pattern = re.compile(r"^\s*(\w+)\s*=\s*(.+?)\s*$")
    setting_pattern = re.compile(r"^\s*Setting\s+(\w+)\s*=\s*(.+?)\s*$")
    result = {}
    for line in lines:
        m = setting_pattern.match(line) or kv_pattern.match(line)
        if m:
            result[m.group(1).lower()] = m.group(2).strip()
    return result


def _parse_log(path, stop_pattern):
    """Core Fortran log parser: extract ``{fortran_name_lower: float}`` pairs.

    Reads *path* line by line, stopping before the first line matching
    *stop_pattern*.  Extracts ``key = value`` pairs and keeps only those
    whose values are parseable as ``float``.  Non-numeric values (e.g.
    ``filetype=binary``) are silently skipped.

    No name resolution, transforms, units, or error-checking is performed.

    :param path: Path to the Fortran log file.
    :type path: str | Path
    :param stop_pattern: Compiled :mod:`re` pattern; parsing stops at the
        first line that matches.
    :return: ``{fortran_name_lower: float}`` for all numeric KV pairs.
    :rtype: dict[str, float]
    """
    path = Path(path)
    with path.open("r", errors="replace") as fh:
        lines = []
        for line in fh:
            if stop_pattern.search(line):
                break
            lines.append(line)

    raw = _parse_fortran_kv(lines)

    result = {}
    for fortran_lower, raw_str in raw.items():
        try:
            result[fortran_lower] = float(raw_str)
        except (ValueError, TypeError):
            pass  # Non-numeric; skip silently
    return result


def parse_micro_log(path, *, stop_pattern=None):
    """Parse a Fortran microscale log file into raw parameter values.

    Reads from the start of *path* until the first line matching
    *stop_pattern* (exclusive).  Returns ``{fortran_name_lower: float}``
    with no name resolution or transforms applied.

    :param path: Path to the micro log file.
    :type path: str | Path
    :param stop_pattern: Compiled regex or string pattern; parsing stops at
        the first matching line.  Defaults to ``r"^\\s*stats\\s*="``.
    :type stop_pattern: re.Pattern | str | None
    :return: ``{fortran_name_lower: float}`` for all numeric KV pairs.
    :rtype: dict[str, float]
    """
    if stop_pattern is None:
        stop_pattern = re.compile(r"^\s*stats\s*=")
    elif isinstance(stop_pattern, str):
        stop_pattern = re.compile(stop_pattern)
    return _parse_log(path, stop_pattern)


def parse_macro_log(path, *, stop_pattern=None):
    """Parse a Fortran macroscale log file into raw parameter values.

    Reads from the start of *path* until the first line beginning with
    ``After `` (exclusive).  Returns ``{fortran_name_lower: float}`` with no
    name resolution or transforms applied.

    :param path: Path to the macro log file.
    :type path: str | Path
    :param stop_pattern: Compiled regex or string pattern; parsing stops at
        the first matching line.  Defaults to ``r"^After\\s"``.
    :type stop_pattern: re.Pattern | str | None
    :return: ``{fortran_name_lower: float}`` for all numeric KV pairs.
    :rtype: dict[str, float]
    """
    if stop_pattern is None:
        stop_pattern = re.compile(r"^After\s")
    elif isinstance(stop_pattern, str):
        stop_pattern = re.compile(stop_pattern)
    return _parse_log(path, stop_pattern)


def parse_micro_file_code(file_code):
    """Extract microscale parameters encoded in a file code string.

    Parses a file code like ``_PLG2_tPA01_Q4`` and returns parameter values
    for any recognized fiber type code segment.  The fiber type codes and
    their associated parameter values are defined in
    :data:`~lysis.config.constants.FIBER_TYPES`.

    :param file_code: File code string (e.g. ``"_PLG2_tPA01_Q4"``).
    :type file_code: str
    :return: ``{python_name: value}`` for recognized fiber type parameters
        (``fiber_radius`` and ``nodes_in_micro_row``).
    :rtype: dict
    """
    from ..config.constants import FIBER_TYPES

    result = {}
    parts = [p for p in file_code.split("_") if p]

    for part in parts:
        if part in FIBER_TYPES:
            result.update(FIBER_TYPES[part])

    return result


def _not_implemented(*args, **kwargs):
    """Placeholder for storage types not yet implemented.

    :raises NotImplementedError: Always raised when called
    """
    raise NotImplementedError("This function is not yet implemented.")


def _read_file_text(
    path: AnyStr,
    spec: DataSetSpec,
    params: BaseParamsType = None,
    sim: int = None,
    file_code: str = "",
) -> np.ndarray:
    """Read a delimited text file (CSV, space-separated, etc.) into a NumPy array.

    Uses numpy.loadtxt to read the file with the delimiter specified in the spec.
    Common for Fortran output files in v1.99.0 format.

    :param path: Directory containing the file (without filename)
    :type path: AnyStr
    :param spec: Dataset specification containing filename pattern, dtype, and delimiter
    :type spec: DataSetSpec
    :param params: Simulation parameters (not used for text files, but kept for interface consistency)
    :type params: BaseParamsType, optional
    :param sim: Simulation index for per-simulation files (e.g., 0, 1, 2...)
    :type sim: int, optional
    :param file_code: Additional code to insert into filename (e.g., "_PLG2_tPA01_TB-xiii.dat")
    :type file_code: str, optional
    :return: Array containing the file data
    :rtype: np.ndarray
    """
    return np.loadtxt(
        os.path.join(path, spec.data_location.format(sim=sim, file_code=file_code)),
        dtype=spec.dtype,
        delimiter=spec.delimiter,
    )


def _read_file_binary(
    path: AnyStr,
    spec: DataSetSpec,
    params: BaseParamsType = None,
    sim: int = None,
    file_code: str = "",
) -> np.ndarray:
    """Read a binary file into a NumPy array with proper reshaping.

    Uses numpy.fromfile to read raw binary data, then reshapes it according to
    the shape specification. The shape may reference parameters (e.g., "micro_simulations")
    which are resolved using parse_shape().

    Common for Fortran binary output files in v1.99.0 format.

    :param path: Directory containing the file (without filename)
    :type path: AnyStr
    :param spec: Dataset specification containing filename pattern, dtype, and shape
    :type spec: DataSetSpec
    :param params: Simulation parameters used to resolve dynamic shapes
    :type params: BaseParamsType, optional
    :param sim: Simulation index for per-simulation files
    :type sim: int, optional
    :param file_code: Additional code to insert into filename
    :type file_code: str, optional
    :return: Reshaped array containing the file data
    :rtype: np.ndarray
    """
    dataset = np.fromfile(
        os.path.join(path, spec.data_location.format(sim=sim, file_code=file_code)),
        dtype=spec.dtype,
    )
    return dataset.reshape(parse_shape(spec.shape, params=params))


def _read_file_json(
    path: AnyStr,
    spec: DataSetSpec,
    params: BaseParamsType = None,
    sim: int = None,
    file_code: str = "",
) -> BaseParamsType:
    """Read a JSON file containing parameter dictionaries.

    Uses Python's json library to parse the file. Typically used for reading
    parameter files (micro_params, macro_params) in both v1.99.0 and v2.0.0 formats.

    :param path: Directory containing the file (without filename)
    :type path: AnyStr
    :param spec: Dataset specification containing filename pattern
    :type spec: DataSetSpec
    :param params: Not used for JSON files, kept for interface consistency
    :type params: BaseParamsType, optional
    :param sim: Simulation index for per-simulation files
    :type sim: int, optional
    :param file_code: Additional code to insert into filename
    :type file_code: str, optional
    :return: Dictionary containing the parsed JSON data
    :rtype: BaseParamsType
    """
    with open(
        os.path.join(path, spec.data_location.format(sim=sim, file_code=file_code)), "r"
    ) as file:
        # Use the JSON library to read in the parameters as a dictionary
        data = json.load(file)
    return data


_COLLECTION_LOG_PARSER = {
    ("v1.95.0", "microscale_out"): {
        "parser": parse_micro_log,
        "params_key": "micro_params",
    },
    ("v1.90.0", "microscale_out"): {
        "parser": parse_micro_log,
        "params_key": "micro_params",
        "stop_pattern": r"^\s*init_entry\s*=",
    },
}
"""Map ``(spec_version, collection_name)`` to parser config for FILE_PARSED storage."""


def _read_file_parsed(
    path: AnyStr,
    spec: DataSetSpec,
    params: BaseParamsType = None,
    sim: int = None,
    file_code: str = "",
) -> BaseParamsType:
    """Read parameters by parsing a Fortran log file.

    Delegates to the appropriate log parser based on the spec version and
    collection name.  The parser returns ``{fortran_name_lower: float}``.
    File-code parameters are merged in for micro logs.

    :param path: Directory containing the log file.
    :type path: AnyStr
    :param spec: Dataset specification; ``spec.version`` and ``spec.collection``
        select the parser and ``spec.data_location`` provides the filename
        template.
    :type spec: DataSetSpec
    :param params: Not used; kept for interface consistency.
    :type params: BaseParamsType, optional
    :param sim: Not used; kept for interface consistency.
    :type sim: int, optional
    :param file_code: Code inserted into the filename template.
    :type file_code: str, optional
    :return: ``{params_key: {fortran_name_lower: float, ...}}`` ready for
        merging into the data-collection params dict.
    :rtype: BaseParamsType
    :raises NotImplementedError: If ``spec.collection`` is not supported.
    """
    if (spec.version, spec.collection) not in _COLLECTION_LOG_PARSER:
        raise NotImplementedError(
            f"FILE_PARSED not supported for collection '{spec.collection}'"
        )

    config = _COLLECTION_LOG_PARSER[spec.version, spec.collection]
    parser = config["parser"]
    params_key = config["params_key"]
    filepath = os.path.join(
        path, spec.data_location.format(sim=sim, file_code=file_code)
    )

    kwargs = {}
    if "stop_pattern" in config:
        kwargs["stop_pattern"] = config["stop_pattern"]
    parsed = parser(filepath, **kwargs)

    # Fill in missing micro parameters from the file code (e.g. fiber type
    # codes like Q4 → fiber_radius and nodes_in_micro_row).
    if file_code and params_key == "micro_params":
        file_code_params = parse_micro_file_code(file_code)
        for k, v in file_code_params.items():
            if k not in parsed:
                # Convert Quantities to strings for consistency
                if isinstance(v, Quantity):
                    parsed[k] = str(v)
                else:
                    parsed[k] = v

    return {params_key: parsed}


def _validate_hdf5_version(file: h5py.File, path: AnyStr, version: str):
    """Check that an open HDF5 file has the expected dataspec_version attribute.

    Skips validation if ``version`` is empty (e.g., for standalone specs not
    registered in a :class:`~lysis.dataio.dataspec.DataSpec`).

    Also emits a :class:`UserWarning` when the file carries the
    :attr:`~lysis.config.constants.Const.CONVERTED_FROM_ATTR` attribute and the
    source version has an entry in
    :data:`~lysis.dataio.dataspec.CONVERSION_WARNINGS`.

    :param file: An already-open HDF5 file handle.
    :type file: h5py.File
    :param path: Path to the HDF5 file (used in error messages).
    :type path: AnyStr
    :param version: Expected version string (e.g., ``"v2.0.0"``). Empty string
        skips validation.
    :type version: str
    :raises ValueError: If the file has no dataspec_version attribute, or if the
        attribute does not match the expected version.
    """
    if not version:
        return
    found = file.attrs.get(CONST.DATASPEC_VERSION_ATTR)
    if found is None:
        raise ValueError(
            f"HDF5 file has no '{CONST.DATASPEC_VERSION_ATTR}' attribute: {path}"
        )
    if found != version:
        raise ValueError(
            f"Dataspec version mismatch in '{path}': "
            f"file has '{found}', expected '{version}'"
        )
    converted_from = file.attrs.get(CONST.CONVERTED_FROM_ATTR)
    if converted_from and converted_from in CONVERSION_WARNINGS:
        warnings.warn(CONVERSION_WARNINGS[converted_from], UserWarning, stacklevel=2)


def ensure_hdf5_version(path: AnyStr, version: str, data: dict = None):
    """Ensure an HDF5 file exists with the correct dataspec_version attribute.

    If the file does not exist, it is created with the version attribute written
    to the root group. If the file already exists, the version attribute is
    validated against the expected version.

    Skips all checks if ``version`` is empty.

    :param path: Path to the HDF5 file.
    :type path: AnyStr
    :param version: Expected version string (e.g., ``"v2.0.0"``). Empty string
        skips all checks.
    :type version: str
    :param data: Optional params dict (as passed to writer functions). When
        provided and the dict contains the
        :attr:`~lysis.config.constants.Const.CONVERTED_FROM_ATTR` key,
        that value is also written to the new file's root attributes.
    :type data: dict, optional
    :raises ValueError: If the existing file has no dataspec_version attribute,
        or if it does not match the expected version.
    """
    if not version:
        return
    if os.path.exists(path):
        with h5py.File(path, "r") as file:
            _validate_hdf5_version(file, path, version)
    else:
        with h5py.File(path, "w") as file:
            file.attrs[CONST.DATASPEC_VERSION_ATTR] = version
            converted_from = data.get(CONST.CONVERTED_FROM_ATTR) if data else None
            if converted_from:
                file.attrs[CONST.CONVERTED_FROM_ATTR] = converted_from


def _read_hdf5_attr(
    path: AnyStr,
    spec: DataSetSpec,
    params: BaseParamsType = None,
    sim: int = None,
    file_code: str = "",
) -> BaseParamsType:
    """Read HDF5 group attributes as parameter dictionaries.

    Reads all attributes from an HDF5 group and returns them as a nested dictionary.
    Used for reading parameters stored as HDF5 attributes in v2.0.0 format.

    The output structure is: {group_path: {attr_name: attr_value, ...}}

    :param path: Path to the HDF5 file (full file path, not just directory)
    :type path: AnyStr
    :param spec: Dataset specification containing the HDF5 group path
    :type spec: DataSetSpec
    :param params: Not used for HDF5 attributes, kept for interface consistency
    :type params: BaseParamsType, optional
    :param sim: Simulation index for per-simulation files
    :type sim: int, optional
    :param file_code: Additional code to insert into group path
    :type file_code: str, optional
    :return: Nested dictionary containing {group_path: {attribute_name: value}}
    :rtype: BaseParamsType
    """
    out = {}
    out[spec.data_location.format(sim=sim, file_code=file_code)] = {}
    with h5py.File(path, "r") as file:
        _validate_hdf5_version(file, path, spec.version)
        items = file[
            spec.data_location.format(sim=sim, file_code=file_code)
        ].attrs.items()
        for k, v in items:
            out[spec.data_location.format(sim=sim, file_code=file_code)][k] = v
    return out


def _read_hdf5_dataset(
    path: AnyStr,
    spec: DataSetSpec,
    params: BaseParamsType = None,
    sim: int = None,
    file_code: str = "",
) -> np.ndarray:
    """Read an HDF5 dataset into a NumPy array.

    Reads a complete HDF5 dataset using the [:] slice notation to load all data
    into memory. Used for reading numerical datasets in v2.0.0 HDF5 format.

    :param path: Path to the HDF5 file (full file path, not just directory)
    :type path: AnyStr
    :param spec: Dataset specification containing the HDF5 dataset path
    :type spec: DataSetSpec
    :param params: Not used for HDF5 datasets, kept for interface consistency
    :type params: BaseParamsType, optional
    :param sim: Simulation index for per-simulation files
    :type sim: int, optional
    :param file_code: Additional code to insert into dataset path
    :type file_code: str, optional
    :return: Array containing the HDF5 dataset data
    :rtype: np.ndarray
    """
    with h5py.File(path, "r") as file:
        _validate_hdf5_version(file, path, spec.version)
        table = file[spec.data_location.format(sim=sim, file_code=file_code)][:]
    return table


# Registry mapping storage types to their corresponding reader functions
# This enables the dispatcher pattern in read_dataset() - based on the
# dataset_storage_type field in a DataSetSpec, the appropriate reader
# function is automatically selected and called.
#
# Each reader function must accept: (path, spec, params, sim, file_code)
# and return either a NumPy array (for data) or dict (for parameters)
data_readers: dict[
    DataSetSpec,
    Callable[
        [AnyStr, DataSetSpec, BaseParamsType, int, str], np.ndarray | BaseParamsType
    ],
] = {
    CONST.DATASET_STORAGE_TYPE.FILE_TEXT: _read_file_text,  # Delimited text files (CSV, etc.)
    CONST.DATASET_STORAGE_TYPE.FILE_PARSED: _read_file_parsed,  # Parameters parsed from log files
    CONST.DATASET_STORAGE_TYPE.FILE_BINARY: _read_file_binary,  # Raw binary files (Fortran)
    CONST.DATASET_STORAGE_TYPE.FILE_JSON: _read_file_json,  # JSON parameter files
    CONST.DATASET_STORAGE_TYPE.HDF5_ATTR: _read_hdf5_attr,  # HDF5 group attributes (params)
    CONST.DATASET_STORAGE_TYPE.HDF5_DATASET: _read_hdf5_dataset,  # HDF5 datasets (numerical data)
}


def read_dataset(
    path: AnyStr,
    spec: DataSetSpec,
    params: BaseParamsType = None,
    sim: int = None,
    file_code: str = "",
) -> np.ndarray | BaseParamsType:
    """Read a single dataset using the appropriate reader for its storage type.

    This is the mid-level API for reading individual datasets. It automatically
    dispatches to the correct reader function (_read_file_text, _read_file_binary,
    _read_file_json, _read_hdf5_attr, or _read_hdf5_dataset) based on the
    dataset_storage_type specified in the DataSetSpec.

    For most use cases, prefer read_data_collection() which handles complete
    collections automatically.

    :param path: Directory containing files (for file-based storage) or path to
                 HDF5 file (for HDF5-based storage)
    :type path: AnyStr
    :param spec: Dataset specification defining storage type, location, dtype, and shape
    :type spec: DataSetSpec
    :param params: Simulation parameters used for resolving dynamic shapes in binary files
    :type params: BaseParamsType, optional
    :param sim: Simulation index for per-simulation storage (e.g., 0, 1, 2...).
                None for combined storage.
    :type sim: int, optional
    :param file_code: Optional string to insert into filename/path patterns.
                      Example: "_PLG2_tPA01_TB-xiii.dat" for Fortran files.
    :type file_code: str, optional
    :return: For numerical datasets: NumPy array. For parameter datasets: dictionary.
    :rtype: np.ndarray | BaseParamsType

    Examples
    --------
    Reading a single Fortran binary file::

        >>> from lysis.dataio.dataspec import dataspec
        >>> spec = dataspec["v1.99.0"]["microscale_out"].data["lysis"]
        >>> data = read_dataset(
        ...     path="/path/to/data",
        ...     spec=spec,
        ...     params=params,
        ...     file_code="_PLG2_tPA01_TB-xiii.dat"
        ... )

    Reading an HDF5 dataset::

        >>> spec = dataspec["v2.0.0"]["microscale_out"].data["sim_final_time"]
        >>> data = read_dataset(
        ...     path="/path/to/file.h5",
        ...     spec=spec
        ... )

    See Also
    --------
    :func:`read_data_collection` : Read complete data collections (recommended)
    """
    return data_readers[spec.dataset_storage_type](
        path, spec, params=params, sim=sim, file_code=file_code
    )


def read_data_collection(
    path: AnyStr,
    collections: list[DataCollectionSpec],
    file_codes: list[str],
) -> DataCollectionType:
    """Read complete data collections from disk.

    This is the high-level API for reading simulation data. It handles:
    - Reading parameters first (required for resolving dynamic shapes)
    - Determining whether simulations are stored separately or combined
    - Iterating through all datasets in each collection
    - Loading per-simulation files until no more are found
    - Merging parameters from multiple collections

    For Fortran data, parameters are returned as raw
    ``{fortran_name_lower: float}`` dicts.  Name resolution, transforms, and
    unit addition are handled downstream by
    :func:`~lysis.dataio.dataconvert.convert_data`.

    :param path: Directory containing files (for Fortran format) or path to
                 HDF5 file (for v2.0.0 format)
    :type path: AnyStr
    :param collections: List of collection specifications to read.
    :type collections: list[DataCollectionSpec]
    :param file_codes: List of file code strings, one per collection.
    :type file_codes: list[str]
    :return: Dictionary containing all loaded datasets plus a ``"params"`` key
        with merged parameters.
    :rtype: DataCollectionType
    :raises FileNotFoundError: If sim=0 file is not found (no data exists)
    :raises KeyError: If required HDF5 dataset/group is not found
    """
    # Initialize output dictionary with empty params
    data = {}
    data["params"] = {}

    # Process each collection in order
    for idx, collection in enumerate(collections):
        # Get the file code for this collection (use empty string if not provided)
        if len(file_codes) > idx:
            file_code = file_codes[idx]
        else:
            file_code = ""

        # Read parameters first - needed for resolving dynamic shapes in later datasets
        if collection.params is not None:
            if (
                collection.params.dataset_storage_type
                == CONST.DATASET_STORAGE_TYPE.FILE_PARSED
            ):
                params = _read_file_parsed(
                    path,
                    collection.params,
                    file_code=file_code,
                )
            else:
                params = read_dataset(
                    path, collection.params, params=None, file_code=file_code
                )
            # Merge parameters from this collection with previously loaded params.
            # Skip None values so old-spec macroscale collections don't erase
            # microscale parameters that were already loaded.
            data["params"] = data["params"] | {
                k: v
                for k, v in params.items()
                if v is not None or k not in data["params"]
            }

        # Read each dataset in the collection
        for name, spec in collection.data.items():
            if collection.simulations_combined is True:
                # All simulations stored together in a single file/dataset
                try:
                    data[name] = read_dataset(
                        path, spec, params=data["params"], file_code=file_code
                    )
                except (FileNotFoundError, KeyError):
                    if not spec.optional:
                        raise
            else:
                # Each simulation stored in a separate file/dataset
                data[name] = []
                sim = 0
                next_sim = True
                # Keep reading simulation files until we run out
                while next_sim:
                    try:
                        dataset = read_dataset(
                            path,
                            spec,
                            params=data["params"],
                            sim=sim,
                            file_code=file_code,
                        )
                    except (FileNotFoundError, KeyError) as e:
                        # If the first simulation file doesn't exist, that's an error
                        if sim == 0:
                            raise e
                        # Otherwise, we've simply read all available simulations
                        else:
                            next_sim = False
                    else:
                        # Successfully read this simulation, add to list
                        data[name].append(dataset)
                        sim += 1

    return data


def _write_file_text(
    data: np.ndarray,
    path: AnyStr,
    spec: DataSetSpec,
    params: BaseParamsType = None,
    sim: int = None,
    file_code: str = "",
    overwrite: bool = False,
):
    """
    Writes an array to disk as a text file

    Uses CONST.get_savetxt_format() to determine the appropriate format string
    for the data's dtype, ensuring proper precision and representation.

    :param data: The array to write.
    :type data: np.ndarray
    :param path: The folder in which to store the file.
        Note: The simulation folder and filename should not be included here as they will be included in the data specification.
    :type path: AnyStr
    :param spec: The specification for the data.
    :type spec: DataSetSpec
    :param params: The parameters matching the data, defaults to None
    :type params: dict[str, Any], optional
    :param sim: The index of the simulation this data is from, if stored individually, defaults to None.
        Format: {"micro_params": dict[str, float | int | Quantity | str], "macro_params": dict[str, float | int | Quantity | str]}
    :type sim: int, optional
    :param file_code: Any code that needs to be attached to the filename, defaults to ""
    :type file_code: str, optional
    :param overwrite: Accepted for interface uniformity; file-based writers always
        overwrite existing files, defaults to False
    :type overwrite: bool, optional
    :raises TypeError: Raised if the data does not match the specification.
    :raises ValueError: Raised if the data's dtype is not supported for text output.
    """
    if not check_dataset_spec(data, spec, params=params):
        raise TypeError(
            f"Data sent for writing does not meet the specification {spec}."
        )

    # Get the appropriate format string for this dtype
    fmt = CONST.get_savetxt_format(data.dtype)

    np.savetxt(
        os.path.join(path, spec.data_location.format(sim=sim, file_code=file_code)),
        data,
        fmt=fmt,
        delimiter=spec.delimiter if not spec.delimiter is None else " ",
    )


def _write_file_binary(
    data: np.ndarray,
    path: AnyStr,
    spec: DataSetSpec,
    params: BaseParamsType = None,
    sim: int = None,
    file_code: str = "",
    overwrite: bool = False,
):
    """Write an array to disk as a raw binary file.

    Inverse of :func:`_read_file_binary`. Uses :meth:`numpy.ndarray.tofile`
    to write the array in its native dtype without any header or delimiter.

    :param data: The array to write.
    :type data: np.ndarray
    :param path: The folder in which to store the file.
    :type path: AnyStr
    :param spec: The specification for the data.
    :type spec: DataSetSpec
    :param params: The parameters matching the data, defaults to None
    :type params: BaseParamsType, optional
    :param sim: The index of the simulation, defaults to None
    :type sim: int, optional
    :param file_code: Code to insert into the filename, defaults to ""
    :type file_code: str, optional
    :param overwrite: Accepted for interface uniformity; file-based writers always
        overwrite existing files, defaults to False
    :type overwrite: bool, optional
    :raises TypeError: Raised if the data does not meet the specification.
    """
    if not check_dataset_spec(data, spec, params=params):
        raise TypeError(
            f"Data sent for writing does not meet the specification {spec}."
        )
    data.tofile(
        os.path.join(path, spec.data_location.format(sim=sim, file_code=file_code))
    )


def _params_json_default(obj):
    """JSON encoder for parameter dicts: numpy scalars → native Python, others → str.

    The inverse of
    :meth:`~lysis.config.parameters.Parameters.to_quantity_or_number`: numpy
    wrapper types are unwrapped to plain Python ``int`` / ``float`` / ``bool``
    so ``json.dump`` writes them as JSON numbers rather than quoted strings.
    Pint :class:`~pint.Quantity` objects (and any other non-serialisable type)
    fall back to their ``str()`` representation, which preserves the unit.

    :param obj: Object to encode.
    :return: JSON-serialisable equivalent.
    """
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return str(obj)


def _write_file_json(
    data: BaseParamsType,
    path: AnyStr,
    spec: DataSetSpec,
    params: BaseParamsType = None,
    sim: int = None,
    file_code: str = "",
    overwrite: bool = False,
):
    """Write a parameter dictionary to a JSON file.

    Inverse of :func:`_read_file_json`. Writes the data dictionary to a
    human-readable JSON file with 4-space indentation.

    NumPy scalar types (``np.integer``, ``np.floating``, ``np.bool_``) are
    written as plain JSON numbers so they round-trip correctly when the file
    is read back.  Pint :class:`~pint.Quantity` values (and any other
    non-serialisable type) fall back to their ``str()`` representation,
    preserving the unit string.

    :param data: Parameter dictionary to write
    :type data: BaseParamsType
    :param path: Directory for the output file
    :type path: AnyStr
    :param spec: Dataset specification containing filename pattern
    :type spec: DataSetSpec
    :param params: Not used for JSON files, kept for interface consistency
    :type params: BaseParamsType, optional
    :param sim: Simulation index for per-simulation files
    :type sim: int, optional
    :param file_code: Additional code to insert into filename
    :type file_code: str, optional
    :param overwrite: Accepted for interface uniformity; file-based writers always
        overwrite existing files, defaults to False
    :type overwrite: bool, optional
    """
    with open(
        os.path.join(path, spec.data_location.format(sim=sim, file_code=file_code)), "w"
    ) as file:
        json.dump(data, file, indent=4, default=_params_json_default)


def _write_hdf5_dataset(
    data: np.ndarray,
    path: AnyStr,
    spec: DataSetSpec,
    params: BaseParamsType = None,
    sim: int = None,
    file_code: str = "",
    overwrite: bool = False,
):
    """
    Writes an array to an HDF5 file as a DataSet

    :param data: The array to write.
    :type data: np.ndarray
    :param path: The folder in which to store the file.
        Note: The simulation folder and filename should not be included here as they will be included in the data specification.
    :type path: AnyStr
    :param spec: The specification for the data.
    :type spec: DataSetSpec
    :param params: The parameters matching the data, defaults to None
    :type params: dict[str, Any], optional
    :param sim: The index of the simulation this data is from, if stored individually, defaults to None.
        Format: {"micro_params": dict[str, float | int | Quantity | str], "macro_params": dict[str, float | int | Quantity | str]}
    :type sim: int, optional
    :param file_code: Any code that needs to be attached to the filename, defaults to ""
    :type file_code: str, optional
    :param overwrite: When ``True`` and the dataset already exists (e.g. a
        zero-length placeholder created by
        :meth:`~lysis.dataio.datastore.DataStore.create`), delete it before
        writing, defaults to False
    :type overwrite: bool, optional
    :raises TypeError: Raised if the data does not match the specification.
    """
    ensure_hdf5_version(path, spec.version)
    if not check_dataset_spec(data, spec, params=params):
        raise TypeError(
            f"Data sent for writing does not meet the specification {spec, data.dtype, data.shape}."
        )
    maxshape = []
    for i in spec.shape:
        if i < 0:
            maxshape.append(None)
        else:
            maxshape.append(i)
    maxshape = tuple(maxshape)
    loc = spec.data_location.format(sim=sim, file_code=file_code)
    with h5py.File(path, "a") as file:
        if overwrite and loc in file:
            del file[loc]
        file.create_dataset(
            loc,
            maxshape=maxshape,
            compression="gzip",
            dtype=spec.dtype,
            data=data.astype(spec.dtype),
            # TODO: Add chunk calculation. Maybe not. Check if it calculates automatically from the size of the data
        )


def _write_hdf5_attr(
    data: BaseParamsType,
    path: AnyStr,
    spec: DataSetSpec,
    params: BaseParamsType = None,
    sim: int = None,
    file_code: str = "",
    overwrite: bool = False,
):
    """Write parameters as HDF5 group attributes.

    Creates or opens an HDF5 group and writes all parameters as attributes to that
    group. Used for storing micro_params and macro_params in v2.0.0 HDF5 format.

    The function converts the data location path to a parameter group name by
    replacing "_data" with "_params" (e.g., "micro_data" → "micro_params").

    :param data: Dictionary containing parameters to write. Should have structure
                 {param_group_name: {key: value, ...}}
    :type data: BaseParamsType
    :param path: Path to the HDF5 file (full file path)
    :type path: AnyStr
    :param spec: Dataset specification containing the HDF5 group path
    :type spec: DataSetSpec
    :param params: Not used for HDF5 attributes, kept for interface consistency
    :type params: BaseParamsType, optional
    :param sim: Simulation index for per-simulation files
    :type sim: int, optional
    :param file_code: Additional code to insert into group path
    :type file_code: str, optional
    :param overwrite: Accepted for interface uniformity; HDF5 attr writer uses
        ``require_group`` which is idempotent, defaults to False
    :type overwrite: bool, optional
    """
    ensure_hdf5_version(path, spec.version, data=data)
    with h5py.File(path, "a") as file:
        # Create or get the HDF5 group for parameters
        group = file.require_group(
            spec.data_location.format(sim=sim, file_code=file_code)
        )
        # Convert data location (e.g., "micro_data/...") to params group (e.g., "micro_params")
        param_group_name = spec.data_location.format(
            sim=sim, file_code=file_code
        ).replace("_data", "_params")
        # Write each parameter as a group attribute
        for k, v in data[param_group_name].items():
            group.attrs[k] = v


# Registry mapping storage types to their corresponding writer functions
# This enables the dispatcher pattern in write_dataset() - based on the
# dataset_storage_type field in a DataSetSpec, the appropriate writer
# function is automatically selected and called.
#
# Each writer function must accept: (data, path, spec, params, sim, file_code)
# and return None (writes occur as side effects to disk)
#
# All five storage types are now implemented.
data_writers: dict[
    DataSetSpec,
    Callable[
        [np.ndarray | BaseParamsType, AnyStr, DataSetSpec, BaseParamsType, int, str],
        None,
    ],
] = {
    CONST.DATASET_STORAGE_TYPE.FILE_TEXT: _write_file_text,  # Delimited text files (CSV, etc.)
    CONST.DATASET_STORAGE_TYPE.FILE_PARSED: _not_implemented,  # Parameters printed in log files
    CONST.DATASET_STORAGE_TYPE.FILE_BINARY: _write_file_binary,  # Raw binary files (Fortran output)
    CONST.DATASET_STORAGE_TYPE.FILE_JSON: _write_file_json,  # JSON parameter files
    CONST.DATASET_STORAGE_TYPE.HDF5_ATTR: _write_hdf5_attr,  # HDF5 group attributes (params)
    CONST.DATASET_STORAGE_TYPE.HDF5_DATASET: _write_hdf5_dataset,  # HDF5 datasets (numerical data)
    None: lambda *args, **kwargs: None,
}


def write_dataset(
    data: DataCollectionType,
    path: AnyStr,
    spec: DataSetSpec,
    params: BaseParamsType = None,
    sim: int = None,
    file_code: str = "",
    overwrite: bool = False,
) -> None:
    """Write a single dataset using the appropriate writer for its storage type.

    This is the mid-level API for writing individual datasets. It automatically
    dispatches to the correct writer function based on the dataset_storage_type
    specified in the DataSetSpec.

    For most use cases, prefer write_data_collection() which handles complete
    collections automatically.

    :param data: Data to write (array for numerical datasets, dict for parameters)
    :type data: DataCollectionType
    :param path: Directory for files (file-based storage) or HDF5 file path
    :type path: AnyStr
    :param spec: Dataset specification defining storage type, location, and dtype
    :type spec: DataSetSpec
    :param params: Simulation parameters used for validation
    :type params: BaseParamsType, optional
    :param sim: Simulation index for per-simulation storage
    :type sim: int, optional
    :param file_code: Optional string to insert into filename/path patterns
    :type file_code: str, optional
    :param overwrite: When ``True``, overwrite an existing HDF5 dataset (ignored
        for file-based writers which always overwrite), defaults to False
    :type overwrite: bool, optional

    See Also
    --------
    :func:`write_data_collection` : Write complete data collections (recommended)
    :func:`read_dataset` : Read individual datasets
    """
    data_writers[spec.dataset_storage_type](
        data, path, spec, params=params, sim=sim, file_code=file_code, overwrite=overwrite
    )


def write_data_collection(
    data: DataCollectionType,
    path: AnyStr,
    collections: list[DataCollectionSpec],
    file_codes: list[str],
):
    """Write complete data collections to disk.

    This is the high-level API for writing simulation data. It handles:
    - Writing parameters first
    - Determining whether simulations should be stored separately or combined
    - Iterating through all datasets in each collection
    - Writing per-simulation files for each simulation in the data

    **Recommended** for most use cases as it handles all the complexity of
    writing multi-collection, multi-simulation data.

    :param data: Dictionary containing all datasets plus a "params" key with
                 parameters. For per-simulation storage, datasets should be
                 lists of arrays. For combined storage, datasets should be
                 single arrays.
    :type data: DataCollectionType
    :param path: Directory for files (v1.99.0 Fortran format) or path to HDF5
                 file (v2.0.0 format)
    :type path: AnyStr
    :param collections: List of collection specifications to write. Order should
                        match how data was organized. Typically from
                        dataspec[version].values()
    :type collections: list[DataCollectionSpec]
    :param file_codes: List of file code strings, one per collection. Use [""]
                       for no file codes.
    :type file_codes: list[str]

    Examples
    --------
    Writing HDF5 v2.0.0 format::

        >>> from lysis.dataio.dataspec import dataspec
        >>> collections = list(dataspec["v2.0.0"].values())
        >>> write_data_collection(
        ...     data=simulation_results,
        ...     path="/path/to/output.h5",
        ...     collections=collections,
        ...     file_codes=[""]
        ... )

    Writing Fortran v1.99.0 format::

        >>> collections = [dataspec["v1.99.0"]["microscale_out"]]
        >>> write_data_collection(
        ...     data=fortran_data,
        ...     path="/path/to/output/directory",
        ...     collections=collections,
        ...     file_codes=["_PLG2_tPA01_TB-xiii.dat"]
        ... )

    See Also
    --------
    :func:`write_dataset` : Write individual datasets
    :func:`read_data_collection` : Read complete data collections
    """
    # Validate that parameters are present when required
    if "params" not in data or not data["params"]:
        for collection in collections:
            if collection.params is not None:
                raise ValueError(
                    "Data collection requires parameters but data['params'] "
                    "is missing or empty. Parameters must always accompany data."
                )

    # Process each collection in order
    for idx, collection in enumerate(collections):
        # Get the file code for this collection (use empty string if not provided)
        if len(file_codes) > idx:
            file_code = file_codes[idx]
        else:
            file_code = ""

        # Write parameters first
        if collection.params is not None:
            write_dataset(
                data["params"],
                path,
                collection.params,
                params=None,
                file_code=file_code,
            )

        # Write each dataset in the collection
        for name, spec in collection.data.items():
            if collection.simulations_combined is True:
                # All simulations in a single file/dataset
                write_dataset(
                    data[name], path, spec, params=data["params"], file_code=file_code
                )
            else:
                # Each simulation in a separate file/dataset
                for sim, table in enumerate(data[name]):
                    write_dataset(
                        table,
                        path,
                        spec,
                        params=data["params"],
                        sim=sim,
                        file_code=file_code,
                    )
