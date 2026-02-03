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

    from lysis.util.fileops import read_data_collection
    from lysis.util.dataspec import dataspec

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

    from lysis.util.fileops import write_data_collection
    from lysis.util.dataspec import dataspec

    # Write v2.0.0 HDF5 format
    collections = list(dataspec["v2.0.0"].values())
    write_data_collection(
        data=simulation_results,
        path="/path/to/output.h5",
        collections=collections,
        file_codes=[""]
    )

Reading/writing individual datasets::

    from lysis.util.fileops import read_dataset, write_dataset

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

:mod:`lysis.util.dataspec` : Data specification definitions
:mod:`lysis.util.dataconvert` : Data format conversion utilities
:mod:`lysis.util.constants` : Storage type constants
"""

import json
import os

from enum import Flag, auto, unique
from typing import Any, AnyStr, Mapping, Union, Callable

import numpy as np
import h5py

from .constants import CONST
from .dataspec import (
    DataCollectionSpec,
    DataSetSpec,
    DataCollectionType,
    UnitParamsType,
    BaseParamsType,
    dataspec,
    parse_shape,
    check_dataset_spec,
)

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2025, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


def _not_implemented(*args, **kwargs):
    raise NotImplementedError("This function is not yet implemented.")


def _read_file_text(
    path: AnyStr,
    spec: DataSetSpec,
    params: BaseParamsType = None,
    sim: int = None,
    file_code: str = "",
) -> np.ndarray:
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
    with open(
        os.path.join(path, spec.data_location.format(sim=sim, file_code=file_code)), "r"
    ) as file:
        # Use the JSON library to read in the parameters as a dictionary
        data = json.load(file)
    return data


def _read_hdf5_attr(
    path: AnyStr,
    spec: DataSetSpec,
    params: BaseParamsType = None,
    sim: int = None,
    file_code: str = "",
) -> BaseParamsType:
    out = {}
    out[spec.data_location.format(sim=sim, file_code=file_code)] = {}
    with h5py.File(path, "r") as file:
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
    with h5py.File(path, "r") as file:
        table = file[spec.data_location.format(sim=sim, file_code=file_code)][:]
    return table


data_readers: dict[
    DataSetSpec,
    Callable[
        [AnyStr, DataSetSpec, BaseParamsType, int, str], np.ndarray | BaseParamsType
    ],
] = {
    CONST.DATASET_STORAGE_TYPE.FILE_TEXT: _read_file_text,
    CONST.DATASET_STORAGE_TYPE.FILE_BINARY: _read_file_binary,
    CONST.DATASET_STORAGE_TYPE.FILE_JSON: _read_file_json,
    CONST.DATASET_STORAGE_TYPE.HDF5_ATTR: _read_hdf5_attr,
    CONST.DATASET_STORAGE_TYPE.HDF5_DATASET: _read_hdf5_dataset,
}


def read_dataset(
    path: AnyStr,
    spec: DataSetSpec,
    params: BaseParamsType = None,
    sim: int = None,
    file_code: str = "",
) -> np.ndarray | BaseParamsType:
    """
    Reads a single table

    This function calls the reader function matching the 
    type of this dataset from the `data_readers` dictionary

    :param path: _description_
    :type path: AnyStr
    :param spec: _description_
    :type spec: DataSetSpec
    :param params: _description_, defaults to None
    :type params: BaseParamsType, optional
    :param sim: _description_, defaults to None
    :type sim: int, optional
    :param file_code: _description_, defaults to ""
    :type file_code: str, optional
    :return: _description_
    :rtype: np.ndarray | BaseParamsType
    """
    return data_readers[spec.dataset_storage_type](
        path, spec, params=params, sim=sim, file_code=file_code
    )


def read_data_collection(
    path: AnyStr,
    collections: list[DataCollectionSpec],
    file_codes: list[str],
) -> DataCollectionType:
    """
    Iterates over the items in a data collection, calling the `read_dataset` function
    for each item.
    This function also determines whether simulations are stored together or separately
    and calls the appropriate functions.
    This function also handles parameter loading

    :param path: _description_
    :type path: AnyStr
    :param collections: _description_
    :type collections: list[DataCollectionSpec]
    :param file_codes: _description_
    :type file_codes: list[str]
    :raises e: _description_
    :return: _description_
    :rtype: DataCollectionType
    """
    data = {}
    data["params"] = {}
    for idx, collection in enumerate(collections):
        if len(file_codes) > idx:
            file_code = file_codes[idx]
        else:
            file_code = ""
        params = read_dataset(path, collection.params, params=None, file_code=file_code)
        data["params"] = data["params"] | params
        for name, spec in collection.data.items():
            if collection.simulations_combined is True:
                data[name] = read_dataset(
                    path, spec, params=params, file_code=file_code
                )
            else:
                data[name] = []
                sim = 0
                next_sim = True
                while next_sim:
                    try:
                        dataset = read_dataset(
                            path, spec, params=params, sim=sim, file_code=file_code
                        )
                    except (FileNotFoundError, KeyError) as e:
                        if sim == 0:
                            raise e
                        else:
                            next_sim = False
                    else:
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
):
    """
    Writes an array to disk as a text file

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
    :raises TypeError: Raised if the data does not match the specification.
    """
    if not check_dataset_spec(data, spec, params=params):
        raise TypeError(
            f"Data sent for writing does not meet the specification {spec}."
        )
    np.savetxt(
        os.path.join(path, spec.data_location.format(sim=sim, file_code=file_code)),
        data,
        delimiter=spec.delimiter,
    )


def _write_hdf5_dataset(
    data: np.ndarray,
    path: AnyStr,
    spec: DataSetSpec,
    params: BaseParamsType = None,
    sim: int = None,
    file_code: str = "",
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
    :raises TypeError: Raised if the data does not match the specification.
    """
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
    with h5py.File(path, "a") as file:
        file.create_dataset(
            spec.data_location.format(sim=sim, file_code=file_code),
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
):
    """
    _summary_

    :param data: _description_
    :type data: dict[str, Any]
    :param path: _description_
    :type path: AnyStr
    :param spec: _description_
    :type spec: DataSetSpec
    :param params: _description_, defaults to None
    :type params: BaseParamsType, optional
    :param sim: _description_, defaults to None
    :type sim: int, optional
    :param file_code: _description_, defaults to ""
    :type file_code: str, optional
    """
    with h5py.File(path, "w") as file:
        group = file.require_group(
            spec.data_location.format(sim=sim, file_code=file_code)
        )
        # Convert data location (e.g., "micro_data/...") to params group (e.g., "micro_params")
        param_group_name = spec.data_location.format(sim=sim, file_code=file_code).replace("_data", "_params")
        for k, v in data[param_group_name].items():
            group.attrs[k] = v


data_writers: dict[
    DataSetSpec,
    Callable[
        [np.ndarray | BaseParamsType, AnyStr, DataSetSpec, BaseParamsType, int, str],
        None,
    ],
] = {
    CONST.DATASET_STORAGE_TYPE.FILE_TEXT: _write_file_text,
    CONST.DATASET_STORAGE_TYPE.FILE_JSON: _not_implemented,
    CONST.DATASET_STORAGE_TYPE.HDF5_ATTR: _write_hdf5_attr,
    CONST.DATASET_STORAGE_TYPE.HDF5_DATASET: _write_hdf5_dataset,
}



def write_dataset(
    data: DataCollectionType,
    path: AnyStr,
    spec: DataSetSpec,
    params: BaseParamsType = None,
    sim: int = None,
    file_code: str = "",
) -> None:
    data_writers[spec.dataset_storage_type](
        data, path, spec, params=params, sim=sim, file_code=file_code
    )


def write_data_collection(
    data: DataCollectionType,
    path: AnyStr,
    collections: list[DataCollectionSpec],
    file_codes: list[str],
):
    for idx, collection in enumerate(collections):
        if len(file_codes) > idx:
            file_code = file_codes[idx]
        else:
            file_code = ""
        write_dataset(
            data["params"], path, collection.params, params=None, file_code=file_code
        )
        for name, spec in collection.data.items():
            if collection.simulations_combined is True:
                write_dataset(
                    data[name], path, spec, params=data["params"], file_code=file_code
                )
            else:
                for sim, table in enumerate(data[name]):
                    write_dataset(
                        table,
                        path,
                        spec,
                        params=data["params"],
                        sim=sim,
                        file_code=file_code,
                    )
