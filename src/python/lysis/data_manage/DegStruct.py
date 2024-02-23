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

"""Code for wrapping a deg_list so that it appears as an ndarray."""

from typing import AnyStr

import numpy as np
import pandas as pd

from ..util import Experiment


class DegStruct(object):
    """A Degradation Structure that takes a deg_list and serves it as if it were a deg_time array


    Args:
        deg_list_file:
        save_t:
    """

    def __init__(self, e: Experiment, deg_list_file: AnyStr, save_t: np.ndarray):
        # Store the Experiment Parameters
        self.e = e
        # Read in the deg_list to a DataFrame
        self.deg_list = pd.read_csv(
            deg_list_file,
            names=["Simulation Time", "Location Index", "Degrade Time"],
        )
        # Convert Fortran indices to Python
        self.deg_list["Location Index"] -= 1
        # Store the list of save times
        self.save_t = save_t
        self.shape = (self.save_t.shape[0], self.e.macro_params.total_edges)

    def __getitem__(self, item: int) -> np.ndarray:
        """

        Args:
            item: The index in the time sequence for which you want to get an array of deg_times

        Returns:
            np.ndarray:

        """
        # Check if a specific piece is being requested
        # TODO(bpaynter): Use this more intelligently to only process the locations being requested
        if isinstance(item, tuple):
            passthrough = item[1:]
            item = item[0]
        else:
            passthrough = None
        # Ensure that we are getting an appropriate index
        # TODO(bpaynter): Add support for slices and iterables
        if isinstance(item, slice):
            raise NotImplementedError("Slicing not supported.")
        # Filter only deg_times that would have been set by this point
        valid_times = self.deg_list[
            self.deg_list["Simulation Time"]
            <= self.save_t[item] + self.e.macro_params.time_step / 2
        ]
        # Filter out fibers with multiple binds (keep the one with the earliest deg time
        # NOTE: This currently assumes that we only record updated deg_times in the list
        valid_times = valid_times.drop_duplicates(
            subset=["Location Index"], keep="last"
        )
        # Sort by location and fill in locations with no deg_time
        valid_times = valid_times.sort_values("Location Index")
        valid_times = valid_times.set_index("Location Index")
        valid_times = valid_times.reindex(
            pd.Index(range(self.e.macro_params.total_edges), name="Location Index"),
            axis="index",
        )
        valid_times = valid_times.fillna(9.9e100)
        # TODO(bpaynter): Add entries for empty rows
        # Pull out the column of deg_times as a numpy array
        return_times = valid_times["Degrade Time"].to_numpy()
        # Return the requested times
        if passthrough is None:
            return return_times
        else:
            return return_times[passthrough]

    def __eq__(self, other: np.ndarray) -> bool:
        """Compare the DegStruct contents to an ndarray.
        Ususally for the purpose of comparing to an older format f_deg_time array

        Args:
            other: The f_deg_time array we wish to compare against

        Returns: True if the f_deg times are all within half a timestep of each other, False otherwise.

        """
        for i in range(len(self.save_t)):
            # Due to issues of float precision between data formats,
            # we just ensure that times are less than half a timestep apart
            if np.any(self[i] - other[i] > self.e.macro_params.time_step / 2):
                return False
        return True

    #
    # def __ne__(self, other: np.ndarray) -> bool:
    #     for i in range(len(self.save_t)):
    #         if np.all(self[i] == other[i]):
    #             return False
    #     return True
