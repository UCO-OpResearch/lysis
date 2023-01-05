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

        self.rng = np.random.default_rng(seed=abs(exp.macro_params.seed))

        self.edge_grid = EdgeGrid(exp)

        molecule_start = (self.rng.integers(exp.macro_params.empty_rows,
                                            size=exp.macro_params.total_molecules),
                          self.rng.integers(exp.macro_params.full_row,
                                            size=exp.macro_params.total_molecules))

        self.molecules = [Molecule(index=i,
                                   location_i=molecule_start[0][i],
                                   location_j=molecule_start[1][i]) for i in range(exp.macro_params.total_molecules)]
        # for m in self.molecules:
        #     self.edge_grid.molecules[m.location_i, m.location_j].append(m.index)

        self.total_macro_unbinds = 0
        self.total_micro_unbinds = 0
        self.total_binds = 0
        self.independent_binds = 0

    def unbind_by_degradation(self, m, current_time):
        m.bound = False
        m.unbound_by_degradation = current_time
        self.total_macro_unbinds += 1
        m.waiting_time = (current_time
                          + self.exp.macro_params.average_bind_time
                          - self.exp.macro_params.time_step / 2)
        m.binding_time = float('int')
        m.leaving_time = 0

    def unbind_by_time(self, m: Molecule, current_time: float):
        m.bound = False
        m.binding_time = self.find_binding_time(current_time)
        if self.rng.random() <= self.exp.macro_params.forced_unbind:
            m.waiting_time = (current_time
                              + self.exp.macro_params.average_bind_time
                              - self.exp.macro_params.time_step / 2)
            m.binding_time = float('inf')
            self.total_micro_unbinds += 1

    def bind(self, m: Molecule, current_time: float):
        m.bound = True
        m.waiting_time = 0
        self.total_binds += 1
        r = self.rng.random()
        m.binding_time = self.find_unbinding_time(r, current_time)
        lysis_time = self.find_lysis_time(r, current_time)
        if lysis_time is not None:
            current_lysis_time = self.edge_grid.fiber_status[m.location_i, m.location_j]
            if current_lysis_time < float('inf'):
                self.edge_grid.fiber_status[m.location_i, m.location_j] = lysis_time
                self.independent_binds += 1
            else:
                self.edge_grid.fiber_status[m.location_i, m.location_j] = min(current_lysis_time, lysis_time)

    def move_to_empty_edge(self, m: Molecule, current_time: float):
        neighborhood = [self.edge_grid.neighbor(m.location_i, m.location_j, k) for k in range(8)]
        neighborhood.append((m.location_i, m.location_j))
        while len(neighborhood) > 0:
            neighborhood_index = self.rng.integers(len(neighborhood))
            if self.edge_grid.fiber_status[neighborhood[neighborhood_index]] < current_time:
                neighborhood.pop(neighborhood_index)
            else:
                # self.edge_grid.move_molecule(m.index, (m.location_i, m.location_j), neighborhood[neighborhood_index])
                m.location_i, m.location_j = neighborhood[neighborhood_index]
                neighborhood = []

    def move(self, m: Molecule, current_time: float):
        if m.waiting_time > current_time and m.unbound_by_degradation == current_time:
            self.move_to_empty_edge(m, current_time)
        else:
            destination = self.edge_grid.neighbor(m.location_i, m.location_j, self.rng.integers(8))
            # self.edge_grid.move_molecule(m.index, (m.location_i, m.location_j), destination)
            m.location_i, m.location_j = destination
            m.binding_time = self.find_binding_time(current_time)
        if m.location_i == self.exp.macro_params.rows-1 and current_time < m.time_to_reach_back_row:
            m.time_to_reach_back_row = current_time

    def find_binding_time(self, current_time: float) -> float:
        binding_time_interval = -math.log(self.rng.random()) / (self.exp.macro_params.binding_rate
                                                                * self.exp.macro_params.binding_sites)
        return current_time + binding_time_interval - self.exp.macro_params.time_step / 2

    def find_unbinding_time(self, r: float, current_time: float) -> float:
        """Find the experiment time at which the tPA molecule will unbind

        If we think of the binding_time matrix as a function f(x) where binding_time[100*x] = f(x)
        We want to draw uniformly from the range of f, that is, we want the unbinding timestep to be
        f(r) where r~U(0,1).

        If -- magically -- 100 * r is an integer i, this is easy because we just look up binding_time[i],
        but what if it isn't a perfect integer?

        Then f(r) lies in the interval ( f(floor(100*r)), f(ceil(100*r)) ) so we interpolate it linearly
        That is, define a linear function g(x) such that
        g(floor(100*r)) = f(floor(100*r)) and g(ceil(100*r)) = f(ceil(100*r))
        and then use g(r) to approximate f(r).

        Another way to see this is,
        if i is an integer such that 100i < 100r < 100(i+1) and
        lambda is such that 100i + lambda = 100r,
        then f(r) ~ (1-lambda)*f(i/100) + lambda*f( (i+1)/100 )"""

        lam = 100 * r - math.floor(100 * r)
        y1 = self.exp.data.unbinding_time[math.floor(100 * r)]
        y2 = self.exp.data.unbinding_time[math.ceil(100 * r)]
        step = (1 - lam) * y1 + lam * y2
        return current_time + step - self.exp.macro_params.time_step / 2

    def find_lysis_time(self, r: float, current_time: float) -> float | None:
        lysis_r = self.rng.random()
        lysis_bin = math.floor(self.exp.macro_params.microscale_runs / 100 * lysis_r)
        if lysis_bin < self.exp.data.total_lyses[math.floor(100 * r)] - 2:
            if lysis_bin == self.exp.data.total_lyses[math.floor(100 * r)] - 3:
                lysis_step = self.exp.data.lysis_time[math.floor(100 * r), lysis_bin]
            else:
                lam = (self.exp.macro_params.microscale_runs / 100 * lysis_r
                       - math.floor(self.exp.macro_params.microscale_runs / 100 * lysis_r))
                y1 = self.exp.data.lysis_time[math.floor(100 * r), lysis_bin]
                y2 = self.exp.data.lysis_time[math.floor(100 * r), lysis_bin + 1]
                lysis_step = (1 - lam) * y1 + lam * y2
            return current_time + lysis_step - self.exp.macro_params.time_step / 2
        else:
            return None

    def run(self):
        for ts in tqdm(range(self.exp.macro_params.total_time_steps)):
            current_time = ts * self.exp.macro_params.time_step

            for m in self.molecules:
                if m.bound:
                    if self.edge_grid.fiber_status[m.location_i, m.location_j] < current_time:
                        self.unbind_by_degradation(m, current_time)
                    elif m.leaving_time < current_time:
                        self.unbind_by_time(m, current_time)

                if not m.bound:
                    should_bind = (m.binding_time < current_time < self.edge_grid.fiber_status[m.location_i,
                                                                                               m.location_j]
                                   and m.waiting_time < current_time)
                    should_move = self.rng.random() >= self.exp.macro_params.moving_probability
                    if should_move and should_bind:
                        should_bind = (self.rng.random() <= ((current_time - m.binding_time)
                                                             / self.exp.macro_params.time_step))
                        should_move = ~should_bind
                    if should_bind:
                        self.bind(m, current_time)
                    elif should_move:
                        self.move(m, current_time)
