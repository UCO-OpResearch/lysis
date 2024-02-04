# *_* coding: utf-8 *_*
#
# Clot Lysis Simulation
# Copyright (C) 2024  Bradley Paynter & Brittany Bannish
#
# datapackage_fortran.py
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

"""Code for reading data files produced by the 
Fortran microscale and macroscale models"""

import os
import sys

from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd

from ..util import Experiment


@dataclass(eq=False)
class FortranData:
    """Class containing the data output by
    the Fortran macroscale and microscale models."""

    n_saves: list[int]
    f_deg_t: list[np.ndarray]
    save_t: list[np.ndarray]
    m_firstpass_t: list[np.ndarray]
    m_location: list[np.ndarray]
    m_bound: list[np.ndarray]

    def __eq__(self, other):
        a = asdict(self)
        b = asdict(other)
        if len(a) != len(b):
            return False
        for key, data in a.items():
            if b[key] is None:
                return False
            elif len(b[key]) != len(data):
                return False
            else:
                for i in range(len(data)):
                    if b[key][i].dtype == np.float64:
                        if np.any(abs(b[key][i] - data[i]) > 1e-8):
                            return False
                    elif np.any(b[key][i] != data[i]):
                        return False
        return True


# def widen_t_degrade()


# TODO(bpaynter): Reorganize using `Datastore` once the new `Experiment` framework is implemented.
def read_data(e: Experiment, file_code: str) -> FortranData:
    n_saves = [
        np.fromfile(
            os.path.join(
                e.os_path,
                f"{sim:02}",
                f"Nsave{file_code[:-4]}_{sim:02}{file_code[-4:]}",
            ),
            dtype=np.int32,
        )[0]
        for sim in range(e.macro_params.total_trials)
    ]
    n_saves = [n + 1 for n in n_saves]

    tsave = [
        np.fromfile(
            os.path.join(
                e.os_path,
                f"{sim:02}",
                f"tsave{file_code[:-4]}_{sim:02}{file_code[-4:]}",
            )
        )
        for sim in range(e.macro_params.total_trials)
    ]

    mfpt = [
        np.fromfile(
            os.path.join(
                e.os_path,
                f"{sim:02}",
                f"mfpt{file_code[:-4]}_{sim:02}{file_code[-4:]}",
            )
        )
        for sim in range(e.macro_params.total_trials)
    ]

    mol_location = []
    mol_status = []
    for sim in range(e.macro_params.total_trials):
        raw_mol_location = np.fromfile(
            os.path.join(
                e.os_path,
                f"{sim:02}",
                f"m_loc{file_code[:-4]}_{sim:02}{file_code[-4:]}",
            ),
            dtype=np.int32,
        )
        mol_location.append(
            raw_mol_location.reshape(n_saves[sim], e.macro_params.total_molecules) - 1
        )

        raw_mol_status = np.fromfile(
            os.path.join(
                e.os_path,
                f"{sim:02}",
                f"m_bound{file_code[:-4]}_{sim:02}{file_code[-4:]}",
            ),
            dtype=np.int32,
        )
        raw_mol_status = raw_mol_status.astype(np.bool_)
        mol_status.append(
            raw_mol_status.reshape(n_saves[sim], e.macro_params.total_molecules)
        )

    mapped_deg = []
    # # Code for storing the processed "wide" deg array if time is more important than storage
    # filename = os.path.join(e.os_path, f"f_deg_time{file_code[:-4]}.npz")
    # if os.path.isfile(filename):
    #     deg_arrays = np.load(filename)
    #     mapped_deg = [deg for deg in deg_arrays.values()]
    # else:
    for sim in range(e.macro_params.total_trials):
        # If an old format (wide) deg array exists
        if os.path.isfile(
            os.path.join(
                e.os_path,
                f"{sim:02}",
                f"f_deg_time{file_code[:-4]}_{sim:02}{file_code[-4:]}",
            )
        ):
            # Read it in
            raw_deg = np.fromfile(
                os.path.join(
                    e.os_path,
                    f"{sim:02}",
                    f"f_deg_time{file_code[:-4]}_{sim:02}{file_code[-4:]}",
                )
            )
            # Reshape it appropriately and append it to the list
            mapped_deg.append(raw_deg.reshape(n_saves[sim], e.macro_params.total_edges))
        else:
            # Assume a new format (long) deg array exists
            # Read in the file to a Pandas DataFrame
            raw_fiber_events = pd.read_csv(
                os.path.join(
                    e.os_path,
                    f"{sim:02}",
                    f"f_deg_list{file_code[:-4]}_{sim:02}{file_code[-4:]}",
                ),
                names=["Simulation Time", "Fiber Index", "Degrade Time"],
            )

            # Bin the events by save interval
            raw_fiber_events["Save Interval"] = pd.cut(
                raw_fiber_events["Simulation Time"],
                tsave[sim] + (e.macro_params.time_step / 2),
                labels=range(1, len(tsave[sim])),
            )
            # Drop any duplicates (more than one binding in a single save interval)
            raw_fiber_events = raw_fiber_events.drop_duplicates(
                subset=["Fiber Index", "Save Interval"], keep="last"
            )
            # Pivot the "long" table into a "wide" table
            deg_df = raw_fiber_events.pivot(
                index="Fiber Index", columns="Save Interval", values="Degrade Time"
            )
            # Add an initial value for each fiber (essentially infinite)
            deg_df[0] = 9.9e100
            # Add any save intervals that had no changes to degradation
            for t in range(len(tsave[sim])):
                if t not in deg_df.columns:
                    deg_df[t] = np.nan
            # Sort by save interval
            deg_df = deg_df.sort_index(axis=1)
            # Copy the degrade time through any save intervals where the degrade time didn't change
            deg_df = deg_df.ffill(axis=1)
            # Convert the DataFrame to a NumPy array and transpose
            deg = deg_df.to_numpy().T
            # Append to the list
            mapped_deg.append(deg)
        # # More code for storing the processed "wide" deg array
        # deg_arrays = {f"f_deg_time_{sim:02}": deg for (sim, deg) in enumerate(mapped_deg)}
        # np.savez_compressed(os.path.join(e.os_path, f"f_deg_time{file_code[:-4]}"), **deg_arrays)

    return FortranData(
        n_saves=n_saves,
        f_deg_t=mapped_deg,
        save_t=tsave,
        m_firstpass_t=mfpt,
        m_location=mol_location,
        m_bound=mol_status,
    )


def compare_data(
    e_new: Experiment, e_old: Experiment, file_code: str, new_data: FortranData = None
):
    if new_data is None:
        new_data = read_data(e_new, file_code)
    old_data = read_data(e_old, file_code)
    return new_data == old_data


def get_parameters(scenarios=dict[str, Experiment]) -> pd.DataFrame:
    # Get default parameters
    e = Experiment(
        os.path.join("..", "..", "..", "data"), experiment_code="0000-00-00-0000"
    )
    e.initialize_macro_param()
    df = pd.DataFrame([e.macro_params])
    index = pd.Index(scenarios.keys(), name="Scenario")
    parameters_df = pd.DataFrame(index=index, columns=df.columns)
    for scen, exp in scenarios.items():
        # Read parameters into a dataframe
        parameters_df.loc[scen] = pd.DataFrame([exp.macro_params]).loc[0]
    return parameters_df
