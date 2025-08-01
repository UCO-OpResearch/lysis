import os

from collections.abc import Callable
from dataclasses import asdict, dataclass, field
from typing import Any, AnyStr

import h5py
import numpy as np

from pint import Quantity

from .constants import CONST, DataSetType
from .parameters import read_param_file

__author__ = "Bradley Paynter"
__copyright__ = "Copyright 2025, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


@dataclass(frozen=True)
class DataSetSpec:
    dataset_type: DataSetType
    dtype: np.dtype
    data_location: str = None
    shape: tuple[int] = (-1,)
    delimiter: str = None


@dataclass(frozen=True)
class DataCollectionSpec:
    simulations_combined: bool
    params: DataSetSpec
    data: dict[str, DataSetSpec]


@dataclass(frozen=True)
class DataSpec:
    microscale_out: DataCollectionSpec
    macroscale_in: DataCollectionSpec
    macroscale_out: DataCollectionSpec


dataspec: dict[str, DataSpec] = {
    "v1.99.0": DataSpec(
        microscale_out=DataCollectionSpec(
            simulations_combined=True,
            params=DataSetSpec(
                data_location="params.json",
                dataset_type=CONST.DATASET_TYPE.FILE_JSON,
                dtype=Quantity,
            ),
            data={
                "micro": DataSetSpec(
                    data_location="micro{file_code}.txt",
                    dataset_type=CONST.DATASET_TYPE.FILE_TEXT,
                    dtype=str,
                    delimiter="\u0000",
                ),
                "micro_rates": DataSetSpec(
                    data_location="micro_rates.f90",
                    dataset_type=CONST.DATASET_TYPE.FILE_TEXT,
                    dtype=str,
                    delimiter="\u0000",
                ),
                "firstPLi": DataSetSpec(
                    data_location="firstPLi{file_code}.dat",
                    dataset_type=CONST.DATASET_TYPE.FILE_BINARY,
                    dtype=np.float64,
                ),
                "lasttPA": DataSetSpec(
                    data_location="lasttPA{file_code}.dat",
                    dataset_type=CONST.DATASET_TYPE.FILE_BINARY,
                    dtype=np.int32,
                ),
                "lyscomplete": DataSetSpec(
                    data_location="lyscomplete{file_code}.dat",
                    dataset_type=CONST.DATASET_TYPE.FILE_BINARY,
                    dtype=np.int32,
                ),
                "lysis": DataSetSpec(
                    data_location="lysis{file_code}.dat",
                    dataset_type=CONST.DATASET_TYPE.FILE_BINARY,
                    dtype=np.float64,
                ),
                "PLi": DataSetSpec(
                    data_location="PLi{file_code}.dat",
                    dataset_type=CONST.DATASET_TYPE.FILE_BINARY,
                    dtype=np.int32,
                ),
                "tPA_time": DataSetSpec(
                    data_location="tPA_time{file_code}.dat",
                    dataset_type=CONST.DATASET_TYPE.FILE_BINARY,
                    dtype=np.float64,
                ),
                "tPAPLiunbd": DataSetSpec(
                    data_location="tPAPLiunbd{file_code}.dat",
                    dataset_type=CONST.DATASET_TYPE.FILE_BINARY,
                    dtype=np.int32,
                ),
                "tPAunbind": DataSetSpec(
                    data_location="tPAunbind{file_code}.dat",
                    dataset_type=CONST.DATASET_TYPE.FILE_BINARY,
                    dtype=np.int32,
                ),
            },
        ),
        macroscale_in=DataCollectionSpec(
            simulations_combined=True,
            params=None,
            data={
                "tPAleave": DataSetSpec(
                    data_location="tPAleave{file_code}.dat",
                    dataset_type=CONST.DATASET_TYPE.FILE_TEXT,
                    dtype=np.float64,
                    shape=(101,),
                ),
                "tsectPA": DataSetSpec(
                    data_location="tsectPA{file_code}.dat",
                    dataset_type=CONST.DATASET_TYPE.FILE_TEXT,
                    dtype=np.float64,
                    shape=(101,),
                ),
                "lysismat": DataSetSpec(
                    data_location="lysismat{file_code}.dat",
                    dataset_type=CONST.DATASET_TYPE.FILE_TEXT,
                    dtype=np.float64,
                    shape=(-1, 100),
                ),
                "lenlysisvect": DataSetSpec(
                    data_location="lenlysisvect{file_code}.dat",
                    dataset_type=CONST.DATASET_TYPE.FILE_TEXT,
                    dtype=np.int32,
                    shape=(100,),
                ),
                "neighbors": DataSetSpec(
                    data_location="neighbors{file_code}.dat",
                    dataset_type=CONST.DATASET_TYPE.FILE_TEXT,
                    dtype=np.int32,
                    shape=(8, -1),
                ),
            },
        ),
        macroscale_out=DataCollectionSpec(
            simulations_combined=False,
            params=DataSetSpec(
                data_location="params.json",
                dataset_type=CONST.DATASET_TYPE.FILE_JSON,
                dtype=Quantity,
            ),
            data={
                "macro": DataSetSpec(
                    data_location="{sim:02}/macro{file_code}_{sim:02}.txt",
                    dataset_type=CONST.DATASET_TYPE.FILE_TEXT,
                    dtype=str,
                    delimiter="\u0000",
                ),
                "Nsave": DataSetSpec(
                    data_location="{sim:02}/Nsave{file_code}_{sim:02}.dat",
                    dataset_type=CONST.DATASET_TYPE.FILE_BINARY,
                    dtype=np.int32,
                    shape=(),
                ),
                "tsave": DataSetSpec(
                    data_location="{sim:02}/tsave{file_code}_{sim:02}.dat",
                    dataset_type=CONST.DATASET_TYPE.FILE_BINARY,
                    dtype=np.float64,
                ),
                "f_deg_list": DataSetSpec(
                    data_location="{sim:02}/f_deg_list{file_code}_{sim:02}.dat",
                    dataset_type=CONST.DATASET_TYPE.FILE_TEXT,
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
                    dataset_type=CONST.DATASET_TYPE.FILE_BINARY,
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
                    dataset_type=CONST.DATASET_TYPE.FILE_BINARY,
                    dtype=np.int32,
                    shape=("total_molecules", -1),
                ),
                "m_bound": DataSetSpec(
                    data_location="{sim:02}/m_bound{file_code}_{sim:02}.dat",
                    dataset_type=CONST.DATASET_TYPE.FILE_BINARY,
                    dtype=np.int32,
                    shape=("total_molecules", -1),
                ),
                "mfpt": DataSetSpec(
                    data_location="{sim:02}/mfpt{file_code}_{sim:02}.dat",
                    dataset_type=CONST.DATASET_TYPE.FILE_BINARY,
                    dtype=np.float64,
                ),
            },
        ),
    ),
    "v2.0.0": DataSpec(
        microscale_out=DataCollectionSpec(
            simulations_combined=True,
            params=DataSetSpec(
                data_location="micro_data",
                dataset_type=CONST.DATASET_TYPE.HDF5_ATTR,
                dtype=Quantity,
            ),
            data={
                "micro_log": DataSetSpec(
                    data_location="log_files/micro_log",
                    dataset_type=CONST.DATASET_TYPE.HDF5_DATASET,
                    dtype=str,
                ),
                "pli_first_time": DataSetSpec(
                    data_location="micro_data/pli_first_time",
                    dataset_type=CONST.DATASET_TYPE.HDF5_DATASET,
                    dtype=np.float64,
                ),
                "tpa_final_num": DataSetSpec(
                    data_location="micro_data/tpa_final_num",
                    dataset_type=CONST.DATASET_TYPE.HDF5_DATASET,
                    dtype=np.uint8,
                ),
                "fiber_degraded": DataSetSpec(
                    data_location="micro_data/fiber_degraded",
                    dataset_type=CONST.DATASET_TYPE.HDF5_DATASET,
                    dtype=np.bool,
                ),
                "sim_final_time": DataSetSpec(
                    data_location="micro_data/sim_final_time",
                    dataset_type=CONST.DATASET_TYPE.HDF5_DATASET,
                    dtype=np.float64,
                ),
                "pli_generated_num": DataSetSpec(
                    data_location="micro_data/pli_generated_num",
                    dataset_type=CONST.DATASET_TYPE.HDF5_DATASET,
                    dtype=np.int16,
                ),
                "tpa_leaving_time": DataSetSpec(
                    data_location="micro_data/tpa_leaving_time",
                    dataset_type=CONST.DATASET_TYPE.HDF5_DATASET,
                    dtype=np.float64,
                ),
                "tpa_unbound_by_pli": DataSetSpec(
                    data_location="micro_data/tpa_unbound_by_pli",
                    dataset_type=CONST.DATASET_TYPE.HDF5_DATASET,
                    dtype=np.bool,
                ),
                "tpa_unbound_kinetic": DataSetSpec(
                    data_location="micro_data/tpa_unbound_kinetic",
                    dataset_type=CONST.DATASET_TYPE.HDF5_DATASET,
                    dtype=np.bool,
                ),
            },
        ),
        macroscale_in=None,
        macroscale_out=DataCollectionSpec(
            simulations_combined=False,
            params=DataSetSpec(
                data_location="macro_data",
                dataset_type=CONST.DATASET_TYPE.HDF5_ATTR,
                dtype=Quantity,
            ),
            data={
                "macro_log": DataSetSpec(
                    data_location="log_files/macro_log__sim_{sim:02}",
                    dataset_type=CONST.DATASET_TYPE.HDF5_DATASET,
                    dtype=str,
                ),
                "snapshot_time": DataSetSpec(
                    data_location="macro_data/sim_{sim:02}/snapshot_time",
                    dataset_type=CONST.DATASET_TYPE.HDF5_DATASET,
                    dtype=np.float64,
                ),
                "fiber_degrade_time": DataSetSpec(
                    data_location="macro_data/sim_{sim:02}/fiber_degrade_time",
                    dataset_type=CONST.DATASET_TYPE.HDF5_DATASET,
                    dtype=np.dtype(
                        [
                            ("Simulation Time Elapsed", np.float64),
                            ("Grid Location Row", np.uint32),
                            ("Grid Location Rank", np.uint32),
                            ("Fiber New Degrade Time", np.float64),
                        ]
                    ),
                ),
                "tpa_bind_events": DataSetSpec(
                    data_location="macro_data/sim_{sim:02}/tpa_bind_events",
                    dataset_type=CONST.DATASET_TYPE.HDF5_DATASET,
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
                "tpa_location_snapshot": DataSetSpec(
                    data_location="macro_data/sim_{sim:02}/tpa_location_snapshot",
                    dataset_type=CONST.DATASET_TYPE.HDF5_DATASET,
                    dtype=np.int32,
                    shape=(-1, 2, -1),
                ),
                "tpa_transit_time": DataSetSpec(
                    data_location="macro_data/sim_{sim:02}/tpa_transit_time",
                    dataset_type=CONST.DATASET_TYPE.HDF5_DATASET,
                    dtype=np.float64,
                ),
            },
        ),
    ),
}

dataspec["fortran"] = dataspec["v1.99.0"]
dataspec["hdf5"] = dataspec["v2.0.0"]
dataspec["current"] = dataspec["hdf5"]

# Format: data_converters["From DataSpec", "To DataSpec"]["Name in new DataSpec"] = function(read_data_collection("From DataSpec")) |-> "New array"
data_converters: dict[
    tuple[str, str], dict[str, Callable[[dict[str, np.ndarray]], np.ndarray]]
] = {
    ("v1.99.0", "v2.0.0"): {"log_files/micro_log": lambda x: x[""]},
}

