# *_* coding: utf-8 *_*
#
# Clot Lysis Simulation
# Copyright (C) 2023  Bradley Paynter & Brittany Bannish
#
# constants.py
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

"""Constants used throughout the package

"""

from enum import Enum, IntEnum, unique
from typing import AnyStr

import numpy as np
import pandas as pd

from .data_file import DataType

_cupy_loaded = True
try:
    import cupy as cp
except ImportError:
    _cupy_loaded = False
    cp = None


default_filenames = {
    "unbinding_time": "tsectPA.dat",  # Fortran: tsec1
    # 'leaving_time': "tPAleave.dat",  # Fortran: CDFtPA
    "lysis_time": "lysismat.dat",  # Fortran: lysismat
    "total_lyses": "lenlysisvect.dat",  # Fortran: lenlysismat
    "degradation_state": "deg.p.npy",  # Fortran: degnext
    "molecule_location": "m_loc.p.npy",
    "molecule_state": "m_bound.p.npy",
    "save_time": "tsave.p.npy",  # Fortran: tsave
}


def _save_pandas_json(filepath: AnyStr, data: pd.DataFrame, *args, **kwargs):
    return data.to_json(filepath, *args, **kwargs)


if _cupy_loaded:

    def _save_cupy(filepath: AnyStr, data: cp.ndarray, *args, **kwargs):
        return np.save(filepath, cp.asnumpy(data), *args, **kwargs)

else:
    _save_cupy = None
    

class DataTypes:
    numpy = DataType(name="NumPy", extension=".npy", load=np.load, save=np.save)
    fortran = DataType(name="NumPy", extension=".dat", load=np.fromfile, save=None)
    matlab = DataType(name="Matlab", extension=".dat", load=np.loadtxt, save=None)
    pandas_json = DataType(
        name="Pandas JSON", extension=".json", load=pd.read_json, save=_save_pandas_json
    )
    if _cupy_loaded:
        cupy = DataType(
            name="CuPy",
            extension=".npy",
            load=lambda *a, **kw: cp.array(np.load(*a, **kw)),
            save=_save_cupy,
        )


class Const:
    def __init__(self):
        self.BOUND = BoundaryDirection
        self.BOUND_COND = BoundaryCondition
        self.DIR = FiberDirection
        self.NEIGHBORHOOD = Neighbors()


class Neighbors:
    def __init__(self):
        self.X = ((-1, -1, 0, 0, 0, 0, 0, 0), (-2, 1, -2, 1, -1, -1, 2, 2))
        self.Y = ((0, 0, 1, 1, 0, 1, 0, 1), (1, 1, 1, 1, -1, -1, 2, 2))
        self.Z = ((-1, -1, 0, 0, 0, 0, 0, 0), (-1, -1, -1, -1, -2, -2, 1, 1))
        self.TOP_REFL = (0, 0, -1, -1, 0, 0, 0, 0)
        self.BOTTOM_REFL = (1, 1, 0, 0, 0, 0, 0, 0)
        self.LEFT_REFL = (0, 0, 0, 0, 3, 3, 0, 0)
        self.RIGHT_REFL = (0, 0, 0, 0, 0, 0, -3, -3)


@unique
class RandomDraw(IntEnum):
    BINDING_TIME_WHEN_UNBINDING = 0
    BINDING_TIME_WHEN_MOVING = 1
    MICRO_UNBIND = 2
    MOVE = 3
    UNBINDING_TIME = 4
    LYSIS_TIME = 5
    CONFLICT_RESOLUTION = 6
    RESTRICTED_MOVE = 7


@unique
class BoundaryDirection(IntEnum):
    TOP = 0
    BOTTOM = 1
    LEFT = 2
    RIGHT = 3
    FRONT = 4
    BACK = 5


@unique
class BoundaryCondition(Enum):
    REFLECTING = 0
    PERIODIC = 1
    CONTINUING = 2


@unique
class FiberDirection(Enum):
    UP = 1
    DOWN = -1
    LEFT = 2
    RIGHT = -2
    OUT = 3
    IN = -3
