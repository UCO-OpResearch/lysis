from functools import partial

import numpy as cp
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
        self.rng = cp.random.default_rng(seed=abs(exp.macro_params.seed))

        self.binding_time_factory = self._BindingTimeFactory(self.exp, self.rng)

        edge_lookup = partial(np.ravel_multi_index, dims=(exp.macro_params.rows, exp.macro_params.full_row))

        # self.edge_grid = EdgeGrid(exp)
        self.fiber_status = cp.full(self.exp.macro_params.rows * self.exp.macro_params.full_row,
                                    float('inf'), dtype=cp.float_)
        for i, j in np.ndindex(self.exp.macro_params.empty_rows, self.exp.macro_params.full_row):
            self.fiber_status[edge_lookup((i, j))] = 0
        for j in range(exp.macro_params.cols):
            self.fiber_status[edge_lookup((self.exp.macro_params.rows-1, 3*j,))] = 0

        self.neighbors = cp.array(EdgeGrid.generate_neighborhood_structure(exp))

        location_i = cp.asnumpy(self.rng.integers(exp.macro_params.empty_rows,
                                                  size=exp.macro_params.total_molecules,
                                                  dtype=cp.short))
        location_j = cp.asnumpy(self.rng.integers(exp.macro_params.full_row,
                                                  size=exp.macro_params.total_molecules,
                                                  dtype=cp.short))

        self.location = cp.array(edge_lookup((location_i, location_j)))

        self.m_fiber_status = None
        self.had_macro_degrade = False

        self.bound = cp.full(exp.macro_params.total_molecules, False, dtype=cp.bool_)
        # self.leaving_time = cp.full(exp.macro_params.total_molecules, float('inf'), dtype=cp.float_)
        self.waiting_time = cp.full(exp.macro_params.total_molecules, 0, dtype=cp.float_)
        self.binding_time = cp.full(exp.macro_params.total_molecules, float('inf'), dtype=cp.float_)
        self.unbound_by_degradation = cp.full(exp.macro_params.total_molecules, 0, dtype=cp.float_)
        self.time_to_reach_back_row = cp.full(exp.macro_params.total_molecules, float('inf'), dtype=cp.float_)
        self.xp = cp.arange(101)

        self.total_macro_unbinds = 0
        self.total_micro_unbinds = 0
        self.total_binds = 0
        self.independent_binds = 0
        self.total_moves = 0
        self.timesteps_with_fiber_changes = 0

    def unbind_by_degradation(self, m: cp.ndarray, current_time: float):
        count = int(cp.count_nonzero(m))
        if count == 0:
            return None
        self.bound[m] = False
        self.unbound_by_degradation[m] = current_time
        self.total_macro_unbinds += 1
        self.waiting_time[m] = (current_time
                                + self.exp.macro_params.average_bind_time
                                - self.exp.macro_params.time_step / 2)
        self.binding_time[m] = float('inf')
        self.had_macro_degrade = True

    def unbind_by_time(self, m: cp.ndarray, current_time: float):
        count = int(cp.count_nonzero(m))
        if count == 0:
            return None
        self.bound[m] = False
        forced = cp.full(self.exp.macro_params.total_molecules, False, dtype=cp.bool_)
        forced[m] = self.rng.random(count) <= self.exp.macro_params.forced_unbind
        self.waiting_time[m & forced] = (current_time
                                         + self.exp.macro_params.average_bind_time
                                         - self.exp.macro_params.time_step / 2)
        self.binding_time[m & forced] = float('inf')
        self.total_micro_unbinds += int(cp.count_nonzero(m & (forced <= self.exp.macro_params.forced_unbind)))
        self.binding_time[m & ~forced] = self.find_binding_time(current_time, int(cp.count_nonzero(m & ~forced)))

    def bind(self, m: cp.ndarray, current_time: float):
        count = int(cp.count_nonzero(m))
        if count == 0:
            return
        self.bound[m] = True
        self.waiting_time[m] = 0
        self.total_binds += count

        r: cp.ndarray = self.rng.random(count) * 100
        self.binding_time[m] = self.find_unbinding_time(r, current_time)

        lysis_time = self.find_lysis_time(r, current_time, count)
        lysis = cp.full(self.exp.macro_params.total_molecules, False, dtype=cp.bool_)
        lysis[m] = lysis_time < float('inf')

        new_lysis = cp.full(self.exp.macro_params.total_molecules, False, dtype=cp.bool_)
        new_lysis[m] = lysis[m] & (self.m_fiber_status[m] == float('inf'))
        new_lysis_fibers = self.location[new_lysis]
        self.independent_binds += len(cp.unique(new_lysis_fibers))

        if int(cp.count_nonzero(self.fiber_status[self.location[m]] > lysis_time)) > 0:
            self.timesteps_with_fiber_changes += 1

        self.fiber_status[self.location[m]] = cp.fmin(self.fiber_status[self.location[m]],
                                                      lysis_time)

    def move_to_empty_edge(self, m: cp.ndarray, move_chance: cp.ndarray, current_time: float):
        for n in cp.arange(self.exp.macro_params.total_molecules)[m]:
            neighborhood = [(self.neighbors[self.location[n], k]) for k in range(8)]
            # neighborhood.append(self.location[n])
            while len(neighborhood) > 0:
                neighborhood_index = int((len(neighborhood)+1)
                                         * move_chance[n]
                                         / self.exp.macro_params.moving_probability)
                if neighborhood_index == len(neighborhood):
                    break
                elif self.fiber_status[neighborhood[neighborhood_index]] < current_time:
                    neighborhood.pop(neighborhood_index)
                else:
                    self.location[n] = neighborhood[neighborhood_index]
                    self.total_moves += 1
                    break

    def find_still_stuck(self, m: cp.ndarray, move_chance: cp.ndarray, current_time: float) -> cp.ndarray:
        still_stuck_to_fiber = cp.full(self.exp.macro_params.total_molecules, False, dtype=cp.bool_)
        if self.had_macro_degrade:
            still_stuck_to_fiber[m] = self.waiting_time[m] > current_time
            just_degraded = cp.full(self.exp.macro_params.total_molecules, False, dtype=cp.bool_)
            just_degraded[still_stuck_to_fiber] = self.unbound_by_degradation[still_stuck_to_fiber] == current_time
            cp.logical_and(still_stuck_to_fiber, just_degraded, out=still_stuck_to_fiber)
            self.move_to_empty_edge(still_stuck_to_fiber, move_chance, current_time)
        return still_stuck_to_fiber

    @staticmethod
    def find_free_to_move(m: cp.ndarray, still_stuck_to_fiber: cp.ndarray) -> cp.ndarray:
        cp.logical_and(m, ~still_stuck_to_fiber, out=still_stuck_to_fiber)
        free_to_move = still_stuck_to_fiber
        return free_to_move

    def find_which_neighbor(self, free_to_move: cp.ndarray, move_chance: cp.ndarray) -> cp.ndarray:
        # neighbor = cp.empty(int(cp.count_nonzero(free_to_move)), dtype=cp.byte)
        neighbor = move_chance[free_to_move] * (8 / self.exp.macro_params.moving_probability)
        neighbor = neighbor.astype(int)
        return neighbor

    def actual_move(self, free_to_move: cp.ndarray, neighbor: cp.ndarray):
        self.location[free_to_move] = self.neighbors[self.location[free_to_move], neighbor]

        self.total_moves += int(cp.count_nonzero(free_to_move))

    def check_back_row(self, m: cp.ndarray, current_time: float):
        reached_rear = cp.full(self.exp.macro_params.total_molecules, False)
        reached_rear[m] = (self.location[m] > (self.exp.macro_params.rows - 2)*self.exp.macro_params.full_row)
        first_time = cp.full(self.exp.macro_params.total_molecules, False)
        first_time[reached_rear] = current_time < self.time_to_reach_back_row[reached_rear]
        cp.logical_and(reached_rear, first_time, out=reached_rear)
        self.time_to_reach_back_row[reached_rear] = current_time

    def move(self, m: cp.ndarray, move_chance: cp.ndarray, current_time: float):

        still_stuck_to_fiber = self.find_still_stuck(m, move_chance, current_time)

        free_to_move = MacroscaleRun.find_free_to_move(m, still_stuck_to_fiber)

        neighbor = self.find_which_neighbor(free_to_move, move_chance)

        self.actual_move(free_to_move, neighbor)

        cp.logical_and(free_to_move, self.m_fiber_status >= current_time, out=free_to_move)

        self.binding_time[free_to_move] = self.find_binding_time(current_time, int(cp.count_nonzero(free_to_move)))

        self.check_back_row(m, current_time)

    class _BindingTimeFactory:
        def __init__(self, exp: Experiment, rng: cp.random.Generator):
            self.period = max(1000000, 10 * exp.macro_params.total_molecules)
            self.exp = exp
            self.rng = rng
            self.pointer = self.period
            self.list = cp.empty((self.period,), dtype=cp.double)

        def fill_list(self):
            self.rng.random(dtype=cp.double, out=self.list)
            cp.log(self.list, out=self.list)
            denominator = self.exp.macro_params.binding_rate * self.exp.macro_params.binding_sites
            cp.divide(self.list, denominator, out=self.list)
            cp.subtract(-self.exp.macro_params.time_step / 2, self.list, out=self.list)
            # print(f"Refilled BindingTimeFactory list.")

        def next(self, count: int = 1):
            self.pointer += 1
            if self.pointer >= self.period:
                self.fill_list()
                self.pointer = 0
            if count == 1:
                return self.list[self.pointer]
            else:
                out = cp.empty(count, dtype=cp.double)
                if count + self.pointer <= self.period:
                    out = self.list[self.pointer:self.pointer+count]
                    self.pointer += count - 1
                else:
                    out[:self.period - self.pointer] = self.list[self.pointer:]
                    self.fill_list()
                    out[self.period - self.pointer:] = self.list[:self.pointer + count - self.period]
                    self.pointer += count - self.period - 1
                return out


    def find_binding_time(self, current_time: float, count: int = 1) -> cp.ndarray | None:
        if count > 0:
            return current_time + self.binding_time_factory.next(count)
        else:
            return None

    def find_unbinding_time(self, r: cp.ndarray, current_time: float) -> cp.ndarray:
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
        interp = cp.interp(r, self.xp, self.exp.data.unbinding_time)
        return interp + (current_time - self.exp.macro_params.time_step / 2)

    def find_lysis_time(self, r: cp.ndarray, current_time: float, count: int) -> cp.ndarray:
        lysis_r = self.rng.random(count, dtype=cp.double) * (self.exp.macro_params.microscale_runs / 100)
        interp = cp.full(count, float('inf'), dtype=cp.double)
        r_idx = r.astype(int)
        total_lyses = self.exp.data.total_lyses[r_idx] - 1
        lysis_happens = lysis_r < total_lyses
        interp[~lysis_happens] = float('inf')
        for i in cp.arange(count)[lysis_happens]:
            interp[i] = cp.interp(lysis_r[i],
                                  cp.arange(total_lyses[i]),
                                  self.exp.data.lysis_time[r_idx[i], :total_lyses[i]])
        interp = interp + (current_time - self.exp.macro_params.time_step / 2)

        return interp

    def find_molecule_fiber_status(self) -> cp.ndarray:
        return self.fiber_status[self.location]

    def find_should_bind(self, current_time: float):
        should_bind = (~self.bound
                       & (self.binding_time < current_time)
                       & (self.m_fiber_status > current_time)
                       & (self.waiting_time < current_time))
        return should_bind

    def run(self):
        for ts in tqdm(cp.arange(self.exp.macro_params.total_time_steps)):
            current_time = ts * self.exp.macro_params.time_step

            self.m_fiber_status = self.find_molecule_fiber_status()
            self.had_macro_degrade = False

            self.unbind_by_degradation(self.bound & (self.m_fiber_status < current_time), current_time)
            self.unbind_by_time(self.bound & (self.binding_time < current_time), current_time)

            should_bind = self.find_should_bind(current_time)

            move_chance = cp.ones(self.exp.macro_params.total_molecules, dtype=cp.double)
            move_chance[~self.bound] = self.rng.random(int(cp.count_nonzero(~self.bound)), dtype=cp.double)
            should_move = move_chance < self.exp.macro_params.moving_probability
            conflict = should_bind & should_move
            threshold = cp.full(self.exp.macro_params.total_molecules, -1, dtype=cp.double)
            threshold[conflict] = ((current_time - self.binding_time[conflict])
                                   / self.exp.macro_params.time_step)
            should_bind[conflict] = self.rng.random(int(cp.count_nonzero(conflict))) <= threshold[conflict]
            should_move[conflict] = ~should_bind[conflict]

            self.bind(should_bind, current_time)
            self.move(should_move, move_chance, current_time)

            if ts % 10000 == 10000-1:
                unlysed_fibers = int(cp.count_nonzero(self.fiber_status > current_time))

                if unlysed_fibers == 0:
                    print()
                    print(f"All fibers degraded after {current_time:.2f} sec. Terminating")
                    # break
                else:
                    degraded_fiber_percent = 100 - unlysed_fibers / self.exp.macro_params.total_fibers * 100
                    print()
                    print(f"After {current_time:.2f} sec, {self.exp.macro_params.total_fibers - unlysed_fibers:,} "
                          f"fibers are degraded.")
                    print(f"After {current_time:.2f} sec, {degraded_fiber_percent:.1f}% of fibers are degraded.")

        print(f"Total binds: {self.total_binds:,}")
        print(f"Timesteps with changes to degrade time: {self.timesteps_with_fiber_changes:,}")
        print(f"Total moves: {self.total_moves:,}")
        print(f"Total macro unbinds: {self.total_macro_unbinds:,}")
        print(f"Total micro unbinds: {self.total_micro_unbinds:,}")
