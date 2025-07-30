from dataclasses import asdict, dataclass, field

import h5py
import numpy as np

from .constants import CONST

__author__ = "Bradley Paynter"
__copyright__ = "Copyright 2025, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


@dataclass
class DataCollection:
    simulations_combined: bool
    text: list[str]
    tables: list[str]
    dtype: dict[str, np.dtype]
    shape: dict[str, tuple[int]]
    delimiter: dict[str, str]


@dataclass
class DataSpec:
    microscale_out: DataCollection
    macroscale_in: DataCollection
    macroscale_out: DataCollection


dataspec: dict[str, DataSpec] = {
    "v1.99.0": DataSpec(
        microscale_out=DataCollection(
            simulations_combined=True,
            text=["micro.txt", "micro_rates.f90"],
            tables=[
                "firstPLi",
                "lasttPA",
                "lyscomplete",
                "lysis",
                "PLi",
                "tPA_time",
                "tPAPLiunbd",
                "tPAunbind",
            ],
            dtype={
                "firstPLi": np.float64,
                "lasttPA": np.int32,
                "lyscomplete": np.int32,
                "lysis": np.float64,
                "PLi": np.int32,
                "tPA_time": np.float64,
                "tPAPLiunbd": np.int32,
                "tPAunbind": np.int32,
            },
            shape={
                "firstPLi": (-1,),
                "lasttPA": (-1,),
                "lyscomplete": (-1,),
                "lysis": (-1,),
                "PLi": (-1,),
                "tPA_time": (-1,),
                "tPAPLiunbd": (-1,),
                "tPAunbind": (-1,),
            },
            delimiter={
                "firstPLi": None,
                "lasttPA": None,
                "lyscomplete": None,
                "lysis": None,
                "PLi": None,
                "tPA_time": None,
                "tPAPLiunbd": None,
                "tPAunbind": None,
            },
        ),
        macroscale_in=DataCollection(
            simulations_combined=True,
            text=[],
            tables=["tPAleave", "tsectPA", "lysismat", "lenlysisvect", "neighbors"],
            dtype={
                "tPAleave": np.float64,
                "tsectPA": np.float64,
                "lysismat": np.float64,
                "lenlysisvect": np.int32,
                "neighbors": np.int32,
            },
            shape={
                "tPAleave": (101,),
                "tsectPA": (101,),
                "lysismat": (-1, 100),
                "lenlysisvect": (100,),
                "neighbors": (8, -1),
            },
            delimiter={
                "tPAleave": " ",
                "tsectPA": " ",
                "lysismat": " ",
                "lenlysisvect": " ",
                "neighbors": " ",
            },
        ),
        macroscale_out=DataCollection(
            simulations_combined=False,
            text=["macro.txt"],
            tables=[
                "Nsave",
                "tsave",
                "f_deg_list",
                "m_bind_t",
                "m_loc",
                "m_bound",
                "mfpt",
            ],
            dtype={
                "f_deg_list": np.dtype(
                    [
                        ("Simulation Time Elapsed", np.float64),
                        ("Grid Location Index", np.int32),
                        ("Fiber New Degrade Time", np.float64),
                    ]
                ),
                "m_bind_t": np.dtype(
                    [
                        ("Simulation Time Elapsed", np.float64),
                        ("tPA Molecule Index", np.int32),
                        ("Molecule New Status", np.int32),
                        ("Grid Location Index", np.int32),
                    ]
                ),
                "m_loc": np.int32,
                "m_bound": np.int32,
                "mfpt": np.float64,
                "tsave": np.float64,
                "Nsave": np.int32,
            },
            shape={
                "f_deg_list": (-1,),
                "m_bind_t": (-1,),
                "m_loc": (-1, "tsave"),
                "m_bound": (-1, "tsave"),
                "mfpt": (-1,),
                "tsave": (-1,),
                "Nsave": None,
            },
            delimiter={
                "f_deg_list": ",",
                "m_bind_t": ",",
                "m_loc": None,
                "m_bound": None,
                "mfpt": None,
                "tsave": None,
                "Nsave": None,
            },
        ),
    ),
    "v2.0.0": DataSpec(
        microscale_out=DataCollection(
            simulations_combined=True,
            text=["log_files/micro_log"],
            tables=[
                "pli_first_time",
                "tpa_final_num",
                "fiber_degraded",
                "sim_final_time",
                "pli_generated_num",
                "tpa_leaving_time",
                "tpa_unbound_by_pli",
                "tpa_unbound_kinetic",
            ],
            dtype={
                "pli_first_time": np.float64,
                "tpa_final_num": np.uint8,
                "fiber_degraded": np.bool,
                "sim_final_time": np.float64,
                "pli_generated_num": np.int16,
                "tpa_leaving_time": np.float64,
                "tpa_unbound_by_pli": np.bool,
                "tpa_unbound_kinetic": np.bool,
            },
            shape={
                "pli_first_time": (-1,),
                "tpa_final_num": (-1,),
                "fiber_degraded": (-1,),
                "sim_final_time": (-1,),
                "pli_generated_num": (-1,),
                "tpa_leaving_time": (-1,),
                "tpa_unbound_by_pli": (-1,),
                "tpa_unbound_kinetic": (-1,),
            },
            delimiter=None,
        ),
        macroscale_in=None,
        macroscale_out=DataCollection(
            simulations_combined=False,
            text=["log_files/macro_log__sim_{sim:02}"],
            tables=[
                "f_degfiber_degrade_time_list",
                "tpa_bind_events",
                "tpa_location_snapshot",
                "m_bound",
                "mfpt",
                "tsave",
            ],
            dtype={
                "fiber_degrade_time": np.dtype(
                    [
                        ("Simulation Time Elapsed", np.float64),
                        ("Grid Location Row", np.uint32),
                        ("Grid Location Rank", np.uint32),
                        ("Fiber New Degrade Time", np.float64),
                    ]
                ),
                "tpa_bind_events": np.dtype(
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
                "tpa_location_snapshot": np.int32,
                "m_bound": np.int32,
                "mfpt": np.float64,
                "tsave": np.float64,
            },
            shape={
                "fiber_degrade_time": (-1,),
                "tpa_bind_events": (-1,),
                "tpa_location_snapshot": (-1, 2, "tsave"),
                "m_bound": (-1, "tsave"),
                "mfpt": (-1,),
                "tsave": (-1,),
            },
            delimiter=None,
        ),
    ),
}

dataspec["fortran"] = dataspec["v1.99.0"]
dataspec["hdf5"] = dataspec["v2.0.0"]
dataspec["current"] = dataspec["hdf5"]
