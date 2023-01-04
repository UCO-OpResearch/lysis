import math

import numpy as np

from tqdm import tqdm

from .util import Experiment
from .edge_grid import EdgeGrid
from .molecule import Molecule

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2022, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


class MacroscaleRun:
    def __init__(self, exp: Experiment):
        self.exp = exp
        assert self.exp.macro_params is not None


def run():
    rng = np.random.default_rng(seed=24375893)
    exp = Experiment(r'..\..\data', experiment_code='2022-12-27-1100')
    exp.initialize_macro_param()
    edge_grid = EdgeGrid(exp)
    molecule_start = (rng.integers(exp.macro_params.empty_rows,
                                   size=exp.macro_params.total_molecules),
                      rng.integers(exp.macro_params.full_row,
                                   size=exp.macro_params.total_molecules))
    molecules = [Molecule(location_i=molecule_start[0][i],
                          location_j=molecule_start[1][i]) for i in range(exp.macro_params.total_molecules)]

    for idx, m in enumerate(molecules):
        edge_grid.molecules[m.location_i, m.location_j] = idx

    for ts in tqdm(range(exp.macro_params.total_time_steps)):
        current_time = ts * exp.macro_params.time_step

        for idx, m in enumerate(molecules):
            if m.bound:
                if m.location_i >= exp.macro_params.empty_rows:
                    if edge_grid.fiber_status[m.location_i, m.location_j] < current_time:
                        m.bound = False
                        m.unbound_by_degradation = True
                        m.waiting_time = exp.macro_params.average_bind_time - exp.macro_params.time_step/2
                        m.leaving_time = float('int')
                    elif m.leaving_time < current_time:
                        m.bound = False
                        binding_time_interval = -math.log(rng.random())/(exp.macro_params.binding_rate
                                                                         * exp.macro_params.binding_sites)
                        m.binding_time = current_time + binding_time_interval - exp.macro_params.time_step/2

