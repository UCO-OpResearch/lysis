import os

from enum import Flag, auto, unique
from typing import Any, AnyStr, List, Mapping, Union

import numpy as np
import h5py

from .constants import CONST
from .dataspec import DataCollectionSpec, DataSetSpec, dataspec
from .parameters import read_param_file

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2025, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"




def read_file_text(
    path: AnyStr,
    spec: DataSetSpec,
    params: dict[str, Any] = None,
    sim: int = None,
    file_code: str = "",
) -> np.ndarray:
    return np.loadtxt(
        os.path.join(path, spec.data_location.format(sim=sim, file_code=file_code)),
        dtype=spec.dtype,
        delimiter=spec.delimiter,
    )


def read_file_binary(
    path: AnyStr,
    spec: DataSetSpec,
    params: dict[str, Any] = None,
    sim: int = None,
    file_code: str = "",
) -> np.ndarray:
    dataset = np.fromfile(
        os.path.join(path, spec.data_location.format(sim=sim, file_code=file_code)),
        dtype=spec.dtype,
    )
    shape = []
    for i in spec.shape:
        if isinstance(i, int):
            shape.append(i)
        elif isinstance(i, str):
            shape.append(params[i])
        else:
            raise RuntimeError("Incorrect shape format {i}.")
    return dataset.reshape(tuple(shape))


def read_file_json(
    path: AnyStr,
    spec: DataSetSpec,
    params: dict[str, Any] = None,
    sim: int = None,
    file_code: str = "",
) -> dict[str, Any]:
    micro, macro = read_param_file(
        os.path.join(path, spec.data_location.format(sim=sim, file_code=file_code))
    )
    return micro | macro


def read_file_hdf5(
    path: AnyStr,
    spec: DataSetSpec,
    params: dict[str, Any] = None,
    sim: int = None,
    file_code: str = "",
):
    raise NotImplementedError


def read_hdf5_attr(
    path: AnyStr,
    spec: DataSetSpec,
    params: dict[str, Any] = None,
    sim: int = None,
    file_code: str = "",
):
    raise NotImplementedError


def read_hdf5_group(
    path: AnyStr,
    spec: DataSetSpec,
    params: dict[str, Any] = None,
    sim: int = None,
    file_code: str = "",
):
    raise NotImplementedError


def read_hdf5_dataset(
    path: AnyStr,
    spec: DataSetSpec,
    params: dict[str, Any] = None,
    sim: int = None,
    file_code: str = "",
):
    raise NotImplementedError


data_readers = {
    CONST.DATASET_TYPE.FILE_TEXT: read_file_text,
    CONST.DATASET_TYPE.FILE_BINARY: read_file_binary,
    CONST.DATASET_TYPE.FILE_JSON: read_file_json,
    CONST.DATASET_TYPE.FILE_HDF5: read_file_hdf5,
    CONST.DATASET_TYPE.HDF5_ATTR: read_hdf5_attr,
    CONST.DATASET_TYPE.HDF5_GROUP: read_hdf5_group,
    CONST.DATASET_TYPE.HDF5_DATASET: read_hdf5_dataset,
}



def read_data_set(
    path: AnyStr,
    spec: DataSetSpec,
    params: dict[str, Any] = None,
    sim: int = None,
    file_code: str = "",
) -> np.ndarray | dict[str, Any]:
    return data_readers[spec.dataset_type](
        path, spec, params=params, sim=sim, file_code=file_code
    )


def read_data_collection(
    path: AnyStr,
    collections: list[DataCollectionSpec],
    file_codes: list[str],
) -> dict[str, np.ndarray] | dict[str, list[np.ndarray]]:
    data = {}
    data["params"] = {}
    for idx, collection in enumerate(collections):
        if len(file_codes) > idx:
            file_code = file_codes[idx]
        else:
            file_code = ""
        params = read_data_set(
            path, collection.params, params=None, file_code=file_code
        )
        data["params"] = data["params"] | params
        for name, spec in collection.data.items():
            if collection.simulations_combined is True:
                data[name] = read_data_set(
                    path, spec, params=params, file_code=file_code
                )
            else:
                data[name] = []
                sim = 0
                next_sim = True
                while next_sim:
                    try:
                        dataset = read_data_set(
                            path, spec, params=params, sim=sim, file_code=file_code
                        )
                    except FileNotFoundError as e:
                        if sim == 0:
                            raise e
                        else:
                            next_sim = False
                    else:
                        data[name].append(dataset)
                        sim += 1
    return data
