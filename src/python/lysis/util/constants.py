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

import secrets
import time

from enum import Enum, IntEnum, unique
from math import log, floor

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

tokens = "23456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz"
_last_v8_timestamp = None


def uuid8code():
    num = uuid8()
    length = floor(128 * log(2) / log(len(tokens))) + 1
    result = ""
    while num > 0:
        num, remainder = divmod(num, len(tokens))
        result += tokens[remainder]
    for i in range(length - len(result)):
        result += "0"
    result = result[::-1]
    result = "-".join([result[:4], result[4:8], result[8:12], result[12:16], result[16:]])
    return result


def uuid8():
    r"""UUID version 8 features a time-ordered value field derived from the
    widely implemented and well known Unix Epoch timestamp source, the
    number of nanoseconds seconds since midnight 1 Jan 1970 UTC, leap
    seconds excluded.

    Copied from https://github.com/oittaa/uuid6-python under MIT License
    """

    global _last_v8_timestamp

    nanoseconds = time.time_ns()
    if _last_v8_timestamp is not None and nanoseconds <= _last_v8_timestamp:
        nanoseconds = _last_v8_timestamp + 1
    _last_v8_timestamp = nanoseconds
    timestamp_ms, timestamp_ns = divmod(nanoseconds, 10**6)
    subsec = timestamp_ns * 2**20 // 10**6
    subsec_a = subsec >> 8
    subsec_b = subsec & 0xFF
    uuid_int = (timestamp_ms & 0xFFFFFFFFFFFF) << 80
    uuid_int |= subsec_a << 64
    uuid_int |= subsec_b << 54
    uuid_int |= secrets.randbits(54)
    return uuid_int

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
