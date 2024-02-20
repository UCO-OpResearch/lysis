# *_* coding: utf-8 *_*
#
# Clot Lysis Simulation
# Copyright (C) 2024  Bradley Paynter & Brittany Bannish
#
# DegStruct.py
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

"""Code for wrapping an IntervalTree so that it appears as an ndarray."""

from typing import AnyStr

import numpy as np
import pandas as pd

from ..util import Experiment


class DegStruct(object):
    """A Degradation Structure that takes a deg_list, stores it in an IntervalTree,
    and serves it as if it were a deg_time array


    Args:
        deg_list_file:
        save_t:
    """

    def __init__(self, e: Experiment, deg_list_file: AnyStr, save_t: np.ndarray):
        self.e = e
        self.deg_list = pd.read_csv(
            deg_list_file,
            names=["Simulation Time", "Location Index", "Degrade Time"],
        )
        self.save_t = save_t

    def __getitem__(self, item):
        if isinstance(item, tuple):
            passthrough = item[1:]
            item = item[0]
        else:
            passthrough = None
        if isinstance(item, slice):
            raise NotImplementedError("Slicing not supported.")
        valid_times = self.deg_list[
            self.deg_list["Simulation Time"] <= self.save_t[item]
        ]
        valid_times = valid_times.drop_duplicates(
            subset=["Location Index"], keep="last"
        )
        valid_times = valid_times.sort_values("Location Index")
        valid_times = valid_times.set_index("Location Index")
        valid_times = valid_times.reindex(
            pd.Index(range(self.e.macro_params.total_fibers), name="Location Index"),
            axis="index",
        )
        valid_times = valid_times.fillna(9.9e100)
        return_times = valid_times["Degrade Time"].to_numpy()
        if passthrough is None:
            return return_times
        else:
            return return_times[passthrough]

    def __eq__(self, other):
        print("I've been called!")
        for i in range(len(self.save_t)):
            if np.any(self[i] != other[i]):
                return False
        return True

    #
    # def __ne__(self, other: np.ndarray) -> bool:
    #     for i in range(len(self.save_t)):
    #         if np.all(self[i] == other[i]):
    #             return False
    #     return True
