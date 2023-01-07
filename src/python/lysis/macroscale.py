import math

import numpy as np

from tqdm import tqdm

from .util import Experiment
from .edge_grid import EdgeGrid

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

        # try:
        #     self.rng = KissRandomGenerator(seed=abs(exp.macro_params.seed))
        # except OSError:
        #     print("KISS Random Generator Error. Using NumPy instead.")
        self.rng = np.random.default_rng(seed=abs(exp.macro_params.seed))

        self.binding_time_factory = self._BindingTimeFactory(self.exp, self.rng)

        # self.edge_grid = EdgeGrid(exp)
        self.fiber_status = np.full((self.exp.macro_params.rows, self.exp.macro_params.full_row),
                                    float('inf'), dtype=np.double)
        self.fiber_status[:self.exp.macro_params.empty_rows] = 0

        self.neighbors_i, self.neighbors_j = EdgeGrid.generate_neighborhood_structure(exp)

        self.location_i = self.rng.integers(exp.macro_params.empty_rows,
                                            size=exp.macro_params.total_molecules)
        self.location_j = self.rng.integers(exp.macro_params.full_row,
                                            size=exp.macro_params.total_molecules)
        self.m_fiber_status = None

        self.bound = np.full(exp.macro_params.total_molecules, False, dtype=bool)
        # self.leaving_time = np.full(exp.macro_params.total_molecules, float('inf'), dtype=float)
        self.waiting_time = np.full(exp.macro_params.total_molecules, 0, dtype=float)
        self.binding_time = np.full(exp.macro_params.total_molecules, float('inf'), dtype=float)
        self.unbound_by_degradation = np.full(exp.macro_params.total_molecules, 0, dtype=float)
        self.time_to_reach_back_row = np.full(exp.macro_params.total_molecules, float('inf'), dtype=float)

        self.total_macro_unbinds = 0
        self.total_micro_unbinds = 0
        self.total_binds = 0
        self.independent_binds = 0
        self.total_moves = 0

    def unbind_by_degradation(self, m: np.ndarray, current_time: float):
        count = np.count_nonzero(m)
        if count == 0:
            return None
        self.bound[m] = False
        self.unbound_by_degradation[m] = current_time
        self.total_macro_unbinds += 1
        self.waiting_time[m] = (current_time
                                + self.exp.macro_params.average_bind_time
                                - self.exp.macro_params.time_step / 2)
        self.binding_time[m] = float('inf')

    def unbind_by_time(self, m: np.ndarray, current_time: float):
        count = np.count_nonzero(m)
        if count == 0:
            return None
        self.bound[m] = False
        forced = np.full(self.exp.macro_params.total_molecules, False)
        forced[m] = self.rng.random(count) <= self.exp.macro_params.forced_unbind
        self.waiting_time[m & forced] = (current_time
                                         + self.exp.macro_params.average_bind_time
                                         - self.exp.macro_params.time_step / 2)
        self.binding_time[m & forced] = float('inf')
        self.total_micro_unbinds += np.count_nonzero(m & (forced <= self.exp.macro_params.forced_unbind))
        self.binding_time[m & ~forced] = self.find_binding_time(current_time, np.count_nonzero(m & ~forced))

    def bind(self, m: np.ndarray, current_time: float):
        count = np.count_nonzero(m)
        self.bound[m] = True
        self.waiting_time[m] = 0
        self.total_binds += count

        r = self.rng.random(count)
        self.binding_time[m] = self.find_unbinding_time(r, current_time, count)

        lysis_time = self.find_lysis_time(r, current_time, count)
        lysis = np.full(self.exp.macro_params.total_molecules, False)
        lysis[m] = lysis_time < float('inf')

        new_lysis = np.full(self.exp.macro_params.total_molecules, False)
        new_lysis[m] = lysis[m] & (self.m_fiber_status[m] == float('inf'))
        new_lysis_fibers = np.concatenate((self.location_i[new_lysis],
                                           self.location_j[new_lysis])).reshape((2, -1))
        self.independent_binds += np.unique(new_lysis_fibers, axis=1).shape[1]

        self.fiber_status[self.location_i[lysis],
                          self.location_j[lysis]] = np.fmin(self.fiber_status[self.location_i[lysis],
                                                                              self.location_j[lysis]],
                                                            lysis_time[lysis[m]])

    def move_to_empty_edge(self, m: np.ndarray, move_chance: np.ndarray, current_time: float):
        for n in np.arange(self.exp.macro_params.total_molecules)[m]:
            neighborhood = [(self.neighbors_i[self.location_i[n], self.location_j[n], k],
                             self.neighbors_j[self.location_i[n], self.location_j[n], k]) for k in range(8)]
            neighborhood.append((self.location_i[n], self.location_j[n]))
            while len(neighborhood) > 0:
                neighborhood_index = int(len(neighborhood) * move_chance[n] / self.exp.macro_params.moving_probability)
                if self.fiber_status[neighborhood[neighborhood_index]] < current_time:
                    neighborhood.pop(neighborhood_index)
                else:
                    self.location_i[n], self.location_j[n] = neighborhood[neighborhood_index]
                    neighborhood = []
                    self.total_moves += 1

    def move(self, m: np.ndarray, move_chance: np.ndarray, current_time: float):
        still_stuck_to_fiber = np.full(self.exp.macro_params.total_molecules, False)
        still_stuck_to_fiber[m] = ((self.waiting_time[m] > current_time)
                                   & (self.unbound_by_degradation[m] == current_time))
        self.move_to_empty_edge(still_stuck_to_fiber, move_chance, current_time)

        np.logical_and(m, ~still_stuck_to_fiber, out=still_stuck_to_fiber)
        free_to_move = still_stuck_to_fiber
        # neighbor = np.empty(np.count_nonzero(free_to_move), dtype=int)
        neighbor = move_chance[free_to_move] * (8 / self.exp.macro_params.moving_probability)
        neighbor = neighbor.astype(int)
        destination_i = self.neighbors_i[self.location_i[free_to_move],
                                         self.location_j[free_to_move],
                                         neighbor]
        self.location_j[free_to_move] = self.neighbors_j[self.location_i[free_to_move],
                                                         self.location_j[free_to_move],
                                                         neighbor]
        self.location_i[free_to_move] = destination_i

        self.total_moves = np.count_nonzero(free_to_move)

        self.binding_time[free_to_move] = self.find_binding_time(current_time, np.count_nonzero(free_to_move))

        # reached_rear = np.full(self.exp.macro_params.total_molecules, False)
        free_to_move[m] = ((self.location_i[m] == self.exp.macro_params.rows-1)
                           & (current_time < self.time_to_reach_back_row[m]))
        self.time_to_reach_back_row[free_to_move] = current_time

    class _BindingTimeFactory:
        def __init__(self, exp: Experiment, rng: np.random.Generator):
            self.period = max(1000000, 10 * exp.macro_params.total_molecules)
            self.exp = exp
            self.rng = rng
            self.pointer = self.period
            self.list = np.empty((self.period,), dtype=float)

        def fill_list(self):
            self.rng.random(out=self.list)
            np.log(self.list, out=self.list)
            denominator = self.exp.macro_params.binding_rate * self.exp.macro_params.binding_sites
            np.divide(self.list, denominator, out=self.list)
            np.subtract(-self.exp.macro_params.time_step / 2, self.list, out=self.list)
            # print(f"Refilled BindingTimeFactory list.")

        def next(self, count: int = 1):
            self.pointer += 1
            if self.pointer >= self.period:
                self.fill_list()
                self.pointer = 0
            if count == 1:
                return self.list[self.pointer]
            else:
                out = np.empty(count, dtype=float)
                if count + self.pointer <= self.period:
                    out = self.list[self.pointer:self.pointer+count]
                    self.pointer += count - 1
                else:
                    out[:self.period - self.pointer] = self.list[self.pointer:]
                    self.fill_list()
                    out[self.period - self.pointer:] = self.list[:self.pointer + count - self.period]
                    self.pointer += count - self.period - 1
                return out


    def find_binding_time(self, current_time: float, count: int = 1) -> float:
        # binding_time_interval = -math.log(self.rng.random()) / (self.exp.macro_params.binding_rate
        #                                                         * self.exp.macro_params.binding_sites)
        # return current_time + binding_time_interval - self.exp.macro_params.time_step / 2
        return current_time + self.binding_time_factory.next(count)

    def find_unbinding_time(self, r: np.ndarray, current_time: float, count: int) -> np.ndarray:
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
        # TODO(bpaynter): Use np.interp
        out = np.empty(count, dtype=float)
        for i in np.arange(count):
            lam = 100 * r[i] - math.floor(100 * r[i])
            y1 = self.exp.data.unbinding_time[math.floor(100 * r[i])]
            y2 = self.exp.data.unbinding_time[math.ceil(100 * r[i])]
            step = (1 - lam) * y1 + lam * y2
            out[i] = current_time + step - self.exp.macro_params.time_step / 2
        return out

    def find_lysis_time(self, r: np.ndarray, current_time: float, count: int) -> np.ndarray:
        lysis_r = self.rng.random(count)
        out = np.full(count, float('inf'))
        for i in np.arange(count):
            lysis_bin = math.floor(self.exp.macro_params.microscale_runs / 100 * lysis_r[i])
            if lysis_bin < self.exp.data.total_lyses[math.floor(100 * r[i])] - 2:
                if lysis_bin == self.exp.data.total_lyses[math.floor(100 * r[i])] - 3:
                    lysis_step = self.exp.data.lysis_time[math.floor(100 * r[i]), lysis_bin]
                else:
                    lam = (self.exp.macro_params.microscale_runs / 100 * lysis_r[i]
                           - math.floor(self.exp.macro_params.microscale_runs / 100 * lysis_r[i]))
                    y1 = self.exp.data.lysis_time[math.floor(100 * r[i]), lysis_bin]
                    y2 = self.exp.data.lysis_time[math.floor(100 * r[i]), lysis_bin + 1]
                    lysis_step = (1 - lam) * y1 + lam * y2
                out[i] = current_time + lysis_step - self.exp.macro_params.time_step / 2
        return out

    def run(self):
        for ts in tqdm(np.arange(self.exp.macro_params.total_time_steps)):
            current_time = ts * self.exp.macro_params.time_step

            self.m_fiber_status = self.fiber_status[self.location_i, self.location_j]

            self.unbind_by_degradation(self.bound & (self.m_fiber_status < current_time), current_time)
            self.unbind_by_time(self.bound & (self.binding_time < current_time), current_time)

            should_bind = (~self.bound
                           & (self.binding_time < current_time)
                           & (self.m_fiber_status > current_time)
                           & (self.waiting_time < current_time))
            move_chance = np.ones(self.exp.macro_params.total_molecules)
            move_chance[~self.bound] = self.rng.random(np.count_nonzero(~self.bound))
            should_move = move_chance < self.exp.macro_params.moving_probability
            conflict = should_bind & should_move
            threshold = np.full(self.exp.macro_params.total_molecules, -1)
            threshold[conflict] = ((current_time - self.binding_time[conflict])
                                                             / self.exp.macro_params.time_step)
            should_bind[conflict] = self.rng.random(np.count_nonzero(conflict)) <= threshold[conflict]
            should_move[conflict] = ~should_bind[conflict]

            self.bind(should_bind, current_time)
            self.move(should_move, move_chance, current_time)

        print(f"Total binds: {self.total_binds}")
        print(f"Total moves: {self.total_moves}")


