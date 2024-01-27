# *_* coding: utf-8 *_*
#
# Clot Lysis Simulation
# Copyright (C) 2024  Bradley Paynter & Brittany Bannish
#
# datapackage-fortran.py
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

import numpy as np

from ..util import Experiment


# TODO(bpaynter): Reorganize using `Datastore` once the new `Experiment` framework is implemented.
def read_data(e: Experiment, file_code: str) -> dict[str, list[np.ndarray]]:
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

    mfpt = np.asarray(
        [
            np.fromfile(
                os.path.join(
                    e.os_path,
                    f"{sim:02}",
                    f"mfpt{file_code[:-4]}_{sim:02}{file_code[-4:]}",
                )
            )
            for sim in range(e.macro_params.total_trials)
        ]
    )

    mol_location = []
    mol_status = []
    mapped_deg = []
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

        raw_deg = np.fromfile(
            os.path.join(
                e.os_path,
                f"{sim:02}",
                f"f_deg_time{file_code[:-4]}_{sim:02}{file_code[-4:]}",
            )
        )
        mapped_deg.append(raw_deg.reshape(n_saves[sim], e.macro_params.total_edges))

    return {
        "n_saves": n_saves,
        "f_deg_t": mapped_deg,
        "save_t": tsave,
        "m_firstpass_t": mfpt,
        "m_location": mol_location,
        "m_status": mol_status,
    }


def compare_data(e_new: Experiment, e_old: Experiment, file_code: str):
    new_data = read_data(e_new, file_code)
    old_data = read_data(e_old, file_code)
    if np.count_nonzero(n_save != n_save_old) > 0:
        return False
    for i in range(e.macro_params.total_trials):
        if np.count_nonzero(deg[i] != deg_old[i]) > 0:
            return False
        if np.count_nonzero(tsave[i] != tsave_old[i]) > 0:
            return False
        if np.count_nonzero(mfpt[i] != mfpt_old[i]) > 0:
            return False
        if np.count_nonzero(mol_location[i] != mol_location_old[i]) > 0:
            return False
        if np.count_nonzero(mol_status[i] != mol_status_old[i]) > 0:
            return False
    return True
