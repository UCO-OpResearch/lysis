import os

from collections.abc import Callable
from dataclasses import asdict, dataclass, field
from typing import Any, AnyStr, NewType, Union

import h5py
import numpy as np

from pint import Quantity

from .constants import CONST, DataSetStorageType

__author__ = "Bradley Paynter"
__copyright__ = "Copyright 2025, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


BaseParamsType = NewType("BaseParamsType", dict[str, dict[str, int | float | str]])
UnitParamsType = NewType(
    "UnitParamsType", dict[str, dict[str, int | float | str | Quantity]]
)
DataSetType = Union[np.ndarray | list[np.ndarray] | BaseParamsType]
DataCollectionType = NewType("DataCollectionType", dict[str, DataSetType])


@dataclass(frozen=True)
class DataSetSpec:
    dataset_storage_type: DataSetStorageType
    dtype: np.dtype
    data_location: str | None = None
    shape: tuple[int, ...] = (-1,)
    delimiter: str | None = None


@dataclass(frozen=True)
class DataCollectionSpec:
    simulations_combined: bool
    params: DataSetSpec
    data: dict[str, DataSetSpec]


dataspec: dict[str, dict[str, DataCollectionSpec]] = {
    "v1.99.0": {
        "microscale_out": DataCollectionSpec(
            simulations_combined=True,
            params=DataSetSpec(
                data_location="params.json",
                dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_JSON,
                dtype=Quantity,
            ),
            data={
                "micro_log": DataSetSpec(
                    data_location="micro{file_code}.txt",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
                    dtype=str,
                    delimiter="\u0000",
                ),
                "micro_code": DataSetSpec(
                    data_location="micro_rates.f90",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
                    dtype=str,
                    delimiter="\u0000",
                ),
                "firstPLi": DataSetSpec(
                    data_location="firstPLi{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.float64,
                ),
                "lasttPA": DataSetSpec(
                    data_location="lasttPA{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.int32,
                ),
                "lyscomplete": DataSetSpec(
                    data_location="lyscomplete{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.int32,
                ),
                "lysis": DataSetSpec(
                    data_location="lysis{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.float64,
                ),
                "PLi": DataSetSpec(
                    data_location="PLi{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.int32,
                ),
                "tPA_time": DataSetSpec(
                    data_location="tPA_time{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.float64,
                ),
                "tPAPLiunbd": DataSetSpec(
                    data_location="tPAPLiunbd{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.int32,
                ),
                "tPAunbind": DataSetSpec(
                    data_location="tPAunbind{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.int32,
                ),
            },
        ),
        "macroscale_in": DataCollectionSpec(
            simulations_combined=True,
            params=None,
            data={
                "tPAleave": DataSetSpec(
                    data_location="tPAleave{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
                    dtype=np.float64,
                    shape=(101,),
                ),
                "tsectPA": DataSetSpec(
                    data_location="tsectPA{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
                    dtype=np.float64,
                    shape=(101,),
                ),
                "lysismat": DataSetSpec(
                    data_location="lysismat{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
                    dtype=np.float64,
                    shape=(-1, 100),
                ),
                "lenlysisvect": DataSetSpec(
                    data_location="lenlysisvect{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
                    dtype=np.int32,
                    shape=(100,),
                ),
                "neighbors": DataSetSpec(
                    data_location="neighbors{file_code}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
                    dtype=np.int32,
                    shape=(-1,),
                ),
            },
        ),
        "macroscale_out": DataCollectionSpec(
            simulations_combined=False,
            params=DataSetSpec(
                data_location="params.json",
                dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_JSON,
                dtype=Quantity,
            ),
            data={
                "macro_log": DataSetSpec(
                    data_location="{sim:02}/macro{file_code}_{sim:02}.txt",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
                    dtype=str,
                    delimiter="\u0000",
                ),
                "Nsave": DataSetSpec(
                    data_location="{sim:02}/Nsave{file_code}_{sim:02}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.int32,
                    shape=(),
                ),
                "tsave": DataSetSpec(
                    data_location="{sim:02}/tsave{file_code}_{sim:02}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.float64,
                ),
                "f_deg_list": DataSetSpec(
                    data_location="{sim:02}/f_deg_list{file_code}_{sim:02}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_TEXT,
                    dtype=np.dtype(
                        [
                            ("Simulation Time Elapsed", np.float64),
                            ("Grid Location Index", np.int32),
                            ("Fiber New Degrade Time", np.float64),
                        ]
                    ),
                    delimiter=",",
                ),
                "m_bind_t": DataSetSpec(
                    data_location="{sim:02}/m_bind_t{file_code}_{sim:02}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.dtype(
                        [
                            ("Simulation Time Elapsed", np.float64),
                            ("tPA Molecule Index", np.int32),
                            ("Molecule New Status", np.int32),
                            ("Grid Location Index", np.int32),
                        ]
                    ),
                    delimiter=",",
                ),
                "m_loc": DataSetSpec(
                    data_location="{sim:02}/m_loc{file_code}_{sim:02}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.int32,
                    shape=(-1, "macro_params.total_molecules"),
                ),
                "m_bound": DataSetSpec(
                    data_location="{sim:02}/m_bound{file_code}_{sim:02}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.int32,
                    shape=(-1, "macro_params.total_molecules"),
                ),
                "mfpt": DataSetSpec(
                    data_location="{sim:02}/mfpt{file_code}_{sim:02}.dat",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.FILE_BINARY,
                    dtype=np.float64,
                ),
            },
        ),
    },
    "v2.0.0": {
        "microscale_out": DataCollectionSpec(
            simulations_combined=True,
            params=DataSetSpec(
                data_location="micro_data",
                dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_ATTR,
                dtype=Quantity,
            ),
            data={
                "micro_log": DataSetSpec(
                    data_location="log_files/micro_log",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=h5py.string_dtype(),
                ),
                "pli_first_time": DataSetSpec(
                    data_location="micro_data/pli_first_time",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.float64,
                ),
                "tpa_final_num": DataSetSpec(
                    data_location="micro_data/tpa_final_num",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.uint8,
                ),
                "fiber_degraded": DataSetSpec(
                    data_location="micro_data/fiber_degraded",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.bool,
                ),
                "sim_final_time": DataSetSpec(
                    data_location="micro_data/sim_final_time",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.float64,
                ),
                "pli_generated_num": DataSetSpec(
                    data_location="micro_data/pli_generated_num",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.int16,
                ),
                "tpa_leaving_time": DataSetSpec(
                    data_location="micro_data/tpa_leaving_time",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.float64,
                ),
                "tpa_unbound_by_pli": DataSetSpec(
                    data_location="micro_data/tpa_unbound_by_pli",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.bool,
                ),
                "tpa_unbound_kinetic": DataSetSpec(
                    data_location="micro_data/tpa_unbound_kinetic",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.bool,
                ),
            },
        ),
        "macroscale_in": DataCollectionSpec(
            simulations_combined=False,
            params=None,
            data={
                "bin_edge_proportions": DataSetSpec(
                    data_location=None,
                    dataset_storage_type=None,
                    dtype=np.float64,
                    shape=(101,),
                ),
                "bin_edge_tpa_leaving_time": DataSetSpec(
                    data_location=None,
                    dataset_storage_type=None,
                    dtype=np.float64,
                    shape=(101,),
                ),
                "binned_fiber_degrade_time": DataSetSpec(
                    data_location=None,
                    dataset_storage_type=None,
                    dtype=np.float64,
                    shape=(-1, 100),
                ),
                "binned_fiber_degraded": DataSetSpec(
                    data_location=None,
                    dataset_storage_type=None,
                    dtype=np.uint16,
                    shape=(100,),
                ),
                "edge_grid_neighbors": DataSetSpec(
                    data_location=None,
                    dataset_storage_type=None,
                    dtype=np.uint32,
                    shape=(-1, 8),
                ),
            },
        ),
        "macroscale_out": DataCollectionSpec(
            simulations_combined=False,
            params=DataSetSpec(
                data_location="macro_data",
                dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_ATTR,
                dtype=Quantity,
            ),
            data={
                "macro_log": DataSetSpec(
                    data_location="log_files/macro_log__sim_{sim:02}",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=h5py.string_dtype(),
                ),
                "snapshot_time": DataSetSpec(
                    data_location="macro_data/sim_{sim:02}/snapshot_time",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.float64,
                ),
                "fiber_degrade_time": DataSetSpec(  # Done
                    data_location="macro_data/sim_{sim:02}/fiber_degrade_time",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.dtype(
                        [
                            ("Simulation Time Elapsed", np.float64),
                            ("Grid Location Row", np.uint32),
                            ("Grid Location Rank", np.uint32),
                            ("Fiber New Degrade Time", np.float64),
                        ]
                    ),
                ),
                "tpa_bind_events": DataSetSpec(  # Still need to do
                    data_location="macro_data/sim_{sim:02}/tpa_bind_events",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.dtype(
                        [
                            ("Simulation Time Elapsed", np.float64),
                            ("tPA Molecule Index", np.int64),
                            (
                                "Molecule New Status",
                                h5py.enum_dtype(
                                    {i.name: i.value for i in CONST.MOL_STATUS},
                                    basetype="u1",
                                ),
                            ),
                            ("Grid Location Row", np.uint32),
                            ("Grid Location Rank", np.uint32),
                        ]
                    ),
                ),
                "tpa_location_snapshot": DataSetSpec(  # Still need to do
                    data_location="macro_data/sim_{sim:02}/tpa_location_snapshot",
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.int32,
                    shape=(-1, 2, -1),
                ),
                "tpa_transit_time": DataSetSpec(
                    data_location="macro_data/sim_{sim:02}/tpa_transit_time",  # Still need to do
                    dataset_storage_type=CONST.DATASET_STORAGE_TYPE.HDF5_DATASET,
                    dtype=np.float64,
                ),
            },
        ),
    },
}

# Define tags
tags = {
    "fortran": "v1.99.0",
    "hdf5": "v2.0.0",
    "current": "hdf5",
}

# Add tags to data spec
for k, v in tags.items():
    dataspec[k] = dataspec[v]


def parse_shape(
    shape: tuple[int | str, ...], params: BaseParamsType = None
) -> tuple[int, ...]:
    """
    Fills in any unknown values in the specification shape from the params dictionary.

    :param shape: The shape from the specification.
        May contain integers, which are passed through unchanged,
        or strings, which are used as keys in the params dictionary to look up the correct value.
    :type shape: tuple[int  |  str, ...]
    :param params: The dictionary containing parameters necessary for calculating shape, defaults to None.
    :type params: BaseParamsType, optional
    :raises RuntimeError: Raised if the shape cannot be parsed correctly.
    :return: Returns a shape appropriate for use in numpy methods.
    :rtype: tuple[int, ...]
    """
    parsed_shape = []
    for i in shape:
        if isinstance(i, int):
            parsed_shape.append(i)
        elif isinstance(i, str):
            parts = i.split(".")
            parsed_shape.append(params[parts[0]][parts[1]])
        else:
            raise RuntimeError("Incorrect shape format {i}.")
    return tuple(parsed_shape)


def check_dataset_spec(
    data: np.ndarray, spec: DataSetSpec, params: BaseParamsType = None
) -> bool:
    """
    Checks whether or not an array of data meets the given specification.

    :param data: The array of data to be checked.
    :type data: np.ndarray
    :param spec: The specification to check the data against.
    :type spec: DataSetSpec
    :param params: A dictionary of parameters matching the data and specifications, defaults to None.
    :type params: BaseParamsType, optional
    :return: True if the data matches the specification, False else.
    :rtype: bool
    """
    if not np.can_cast(data.dtype, spec.dtype):
        return False
    for idx, i in enumerate(parse_shape(spec.shape, params=params)):
        if i < 0:
            continue
        elif i != data.shape[idx]:
            return False
    return True
