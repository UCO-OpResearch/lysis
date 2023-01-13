from functools import partial

import numpy as np
from tqdm.auto import tqdm

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

        self.rng = np.random.default_rng(seed=abs(exp.macro_params.seed))

        self.binding_time_factory = self._BindingTimeFactory(self.exp, self.rng)

        edge_lookup = partial(np.ravel_multi_index, dims=(exp.macro_params.rows, exp.macro_params.full_row))

        self.fiber_status = np.full(self.exp.macro_params.rows * self.exp.macro_params.full_row,
                                    float('inf'), dtype=np.float_)
        for i, j in np.ndindex(self.exp.macro_params.empty_rows, self.exp.macro_params.full_row):
            self.fiber_status[edge_lookup((i, j))] = 0
        for j in range(exp.macro_params.cols):
            self.fiber_status[edge_lookup((self.exp.macro_params.rows-1, 3*j,))] = 0

        self.neighbors = EdgeGrid.generate_neighborhood_structure(exp)

        location_i = self.rng.integers(exp.macro_params.empty_rows,
                                       size=exp.macro_params.total_molecules,
                                       dtype=np.short)
        location_j = self.rng.integers(exp.macro_params.full_row,
                                       size=exp.macro_params.total_molecules,
                                       dtype=np.short)

        self.location = edge_lookup((location_i, location_j))

        self.m_fiber_status = None

        self.bound = np.full(exp.macro_params.total_molecules, False, dtype=np.bool_)
        # self.leaving_time = np.full(exp.macro_params.total_molecules, float('inf'), dtype=np.float_)
        self.waiting_time = np.full(exp.macro_params.total_molecules, 0, dtype=np.float_)
        self.binding_time = np.full(exp.macro_params.total_molecules, float('inf'), dtype=np.float_)
        self.unbound_by_degradation = np.full(exp.macro_params.total_molecules, 0, dtype=np.bool_)
        self.time_to_reach_back_row = np.full(exp.macro_params.total_molecules, float('inf'), dtype=np.float_)
        self.reached_back_row = np.full(exp.macro_params.total_molecules, False, dtype=np.bool_)
        self.xp = np.arange(101)

        self.total_macro_unbinds = 0
        self.total_micro_unbinds = 0
        self.total_binds = 0
        self.independent_binds = 0
        self.total_regular_moves = 0
        self.total_restricted_moves = 0
        self.timesteps_with_fiber_changes = 0
        self.number_reached_back_row = 0

    def unbind_by_degradation(self, m: np.ndarray, current_time: float):
        count = np.count_nonzero(m)
        if count == 0:
            return None
        self.bound[m] = False
        self.unbound_by_degradation[m] = True
        self.total_macro_unbinds += count
        self.waiting_time[m] = (current_time
                                + self.exp.macro_params.average_bind_time
                                - self.exp.macro_params.time_step / 2)
        self.binding_time[m] = float('inf')

    def unbind_by_time(self, m: np.ndarray, current_time: float):
        count = np.count_nonzero(m)
        if count == 0:
            return None
        self.bound[m] = False
        self.unbound_by_degradation[m] = False
        forced = np.full(self.exp.macro_params.total_molecules, False, dtype=np.bool_)
        forced[m] = self.rng.random(count) <= self.exp.macro_params.forced_unbind
        self.waiting_time[m & forced] = (current_time
                                         + self.exp.macro_params.average_bind_time
                                         - self.exp.macro_params.time_step / 2)
        self.binding_time[m & forced] = float('inf')
        self.total_micro_unbinds += np.count_nonzero(forced)
        unbinding_times_needed = np.count_nonzero(m & ~forced)
        if unbinding_times_needed > 0:
            self.binding_time[m & ~forced] = current_time + self.binding_time_factory.next(unbinding_times_needed)

    def bind(self, m: np.ndarray, current_time: float):
        count = np.count_nonzero(m)
        if count == 0:
            return
        self.bound[m] = True
        self.waiting_time[m] = 0
        self.total_binds += count

        r: np.ndarray = self.rng.random(count) * 100
        self.binding_time[m] = self.find_unbinding_time(r, current_time)

        lysis_time = self.find_lysis_time(r, current_time, count)
        lysis = np.full(self.exp.macro_params.total_molecules, False, dtype=np.bool_)
        lysis[m] = lysis_time < float('inf')

        new_lysis = np.full(self.exp.macro_params.total_molecules, False, dtype=np.bool_)
        new_lysis[m] = lysis[m] & (self.m_fiber_status[m] == float('inf'))
        new_lysis_fibers = self.location[new_lysis]
        self.independent_binds += len(np.unique(new_lysis_fibers))

        if int(np.count_nonzero(self.fiber_status[self.location[m]] > lysis_time)) > 0:
            self.timesteps_with_fiber_changes += 1

        self.fiber_status[self.location[m]] = np.fmin(self.fiber_status[self.location[m]],
                                                      lysis_time)

    def get_neighbors(self, m: np.ndarray, count: int):
        return np.append(self.neighbors[self.location[m]],
                         self.location[m].reshape(count, 1),
                         axis=1)

    def select_neighbor(self, count: int, potential_neighbors: int):
        # neighborhood_index = np.full(self.exp.macro_params.total_molecules, 0, dtype=int)

        return self.rng.integers(potential_neighbors, size=count)

        # return neighborhood_index.astype(int, copy=True)

    def check_neighbor_empty(self,
                             count: int,
                             need_neighbors: np.ndarray,
                             neighborhoods: np.ndarray,
                             neighborhood_index: np.ndarray,
                             current_time: float):
        valid_neighbors_short = np.full(count, False, dtype=np.bool_)

        valid_neighbors_short[need_neighbors] = (
                self.fiber_status[neighborhoods[need_neighbors,
                                                neighborhood_index[need_neighbors]]] < current_time
        )

        return valid_neighbors_short

    def reset_neighbors(self,
                        neighborhoods: np.ndarray,
                        need_neighbors: np.ndarray,
                        neighborhood_index: np.ndarray,
                        potential_neighbors: int):
        neighborhoods[need_neighbors,
                      neighborhood_index[need_neighbors]] = neighborhoods[need_neighbors, potential_neighbors]

        return neighborhoods

    def move_to_empty_edge(self, m: np.ndarray, current_time: float):
        count = np.count_nonzero(m)
        if count == 0:
            return

        self.total_restricted_moves += count

        neighborhoods = self.get_neighbors(m, count)

        need_neighbors = np.full(count, True, dtype=np.bool_)
        valid_neighbors = np.full(self.exp.macro_params.total_molecules, False, dtype=np.bool_)
        potential_neighbors = 9

        while np.count_nonzero(need_neighbors) > 0:
            if potential_neighbors == 0:
                raise RuntimeError(f"No valid neighbors!")

            neighborhood_index = self.select_neighbor(count, potential_neighbors)

            valid_neighbors_short = self.check_neighbor_empty(count,
                                                              need_neighbors,
                                                              neighborhoods,
                                                              neighborhood_index,
                                                              current_time)

            valid_neighbors[m] = valid_neighbors_short

            self.location[valid_neighbors] = neighborhoods[valid_neighbors_short,
                                                           neighborhood_index[valid_neighbors_short]]

            potential_neighbors -= 1

            neighborhoods = self.reset_neighbors(neighborhoods,
                                                 need_neighbors,
                                                 neighborhood_index,
                                                 potential_neighbors)

            need_neighbors = need_neighbors & ~valid_neighbors_short

        if np.count_nonzero(self.fiber_status[self.location[m]] >= current_time):
            raise RuntimeError(f"We moved to a non-empty edge!")

        # for n in np.arange(self.exp.macro_params.total_molecules)[m]:
        #     # print(f"Time: {current_time:.2f}; Molecule: {n:,}")
        #     neighborhood = list(range(8))
        #     while len(neighborhood) > 0:
        #         neighborhood_index = int((len(neighborhood)+1)
        #                                  * move_chance[n]
        #                                  / self.exp.macro_params.moving_probability)
        #         if neighborhood_index == len(neighborhood):
        #             break
        #         neighbor = self.neighbors[self.location[n], neighborhood[neighborhood_index]]
        #         if self.fiber_status[neighbor] < current_time:
        #             neighborhood.pop(neighborhood_index)
        #         else:
        #             self.location[n] = neighborhood[neighborhood_index]
        #             self.total_moves += 1
        #             break

    def find_still_stuck(self, m: np.ndarray, current_time: float):
        # still_stuck_to_fiber = np.full(self.exp.macro_params.total_molecules, False, dtype=np.bool_)
        still_stuck_to_fiber = self.waiting_time > current_time
        still_stuck_to_fiber = np.logical_and(still_stuck_to_fiber,
                                              self.unbound_by_degradation,
                                              out=still_stuck_to_fiber)
        still_stuck_to_fiber = np.logical_and(still_stuck_to_fiber,
                                              m,
                                              out=still_stuck_to_fiber)
        # count = np.count_nonzero(still_stuck_to_fiber)
        # if count > 0:
        #     print(f"Time: {current_time:.2f}; Molecules stuck to macro fibers: {count:,}")
        return still_stuck_to_fiber

    def actual_move(self, free_to_move: np.ndarray, current_time: float):
        # neighbor = move_chance[free_to_move] * (8 / self.exp.macro_params.moving_probability)
        # neighbor = neighbor.astype(int)
        neighbor = self.rng.integers(8, size=np.count_nonzero(free_to_move))

        self.location[free_to_move] = self.neighbors[self.location[free_to_move], neighbor]
        self.total_regular_moves += np.count_nonzero(free_to_move)

        free_to_move = np.logical_and(free_to_move, self.m_fiber_status >= current_time, out=free_to_move)

        move_to_fiber = np.count_nonzero(free_to_move)
        if move_to_fiber > 0:
            self.binding_time[free_to_move] = current_time + self.binding_time_factory.next(move_to_fiber)

    def find_reached_rear(self, m: np.ndarray, current_time: float):
        reached_rear = np.full(self.exp.macro_params.total_molecules, False)
        reached_rear[m] = (self.location[m] > (self.exp.macro_params.rows - 2) * self.exp.macro_params.full_row)
        # first_time = np.full(self.exp.macro_params.total_molecules, False)
        first_time = current_time < self.time_to_reach_back_row
        reached_rear = np.logical_and(reached_rear, first_time, out=reached_rear)
        self.time_to_reach_back_row[reached_rear] = current_time
        self.number_reached_back_row += np.count_nonzero(reached_rear)

        # if self.reached_back_row == self.exp.macro_params.total_molecules:
        #     return
        # first_time = (~self.reached_back_row
        #               & self.location > (self.exp.macro_params.rows - 2) * self.exp.macro_params.full_row)
        # self.time_to_reach_back_row[first_time] = current_time
        # self.reached_back_row[first_time] = True
        # self.number_reached_back_row += np.count_nonzero(first_time)

    def move(self, m: np.ndarray, current_time: float):
        still_stuck_to_fiber = self.find_still_stuck(m, current_time)
        self.move_to_empty_edge(still_stuck_to_fiber, current_time)

        still_stuck_to_fiber = np.logical_and(m, ~still_stuck_to_fiber, out=still_stuck_to_fiber)
        free_to_move = still_stuck_to_fiber

        self.actual_move(free_to_move, current_time)

        self.find_reached_rear(m, current_time)


    class _BindingTimeFactory:
        def __init__(self, exp: Experiment, rng: np.random.Generator):
            self.period = max(1000000, 10 * exp.macro_params.total_molecules)
            self.exp = exp
            self.rng = rng
            self.pointer = self.period
            self.list = np.empty((self.period,), dtype=np.double)

        def fill_list(self):
            self.list = self.rng.random(dtype=np.double, out=self.list)
            self.list = np.log(self.list, out=self.list)
            denominator = self.exp.macro_params.binding_rate * self.exp.macro_params.binding_sites
            self.list = np.divide(self.list, denominator, out=self.list)
            self.list = np.subtract(-self.exp.macro_params.time_step / 2, self.list, out=self.list)
            # print(f"Refilled BindingTimeFactory list.")

        def next(self, count: int = 1):
            self.pointer += 1
            if self.pointer >= self.period:
                self.fill_list()
                self.pointer = 0
            if count == 1:
                return self.list[self.pointer]
            else:
                out = np.empty(count, dtype=np.double)
                if count + self.pointer <= self.period:
                    out = self.list[self.pointer:self.pointer+count]
                    self.pointer += count - 1
                else:
                    out[:self.period - self.pointer] = self.list[self.pointer:]
                    self.fill_list()
                    out[self.period - self.pointer:] = self.list[:self.pointer + count - self.period]
                    self.pointer += count - self.period - 1
                return out


    def find_unbinding_time(self, r: np.ndarray, current_time: float) -> np.ndarray:
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
        interp = np.interp(r, self.xp, self.exp.data.unbinding_time)
        return interp + (current_time - self.exp.macro_params.time_step / 2)

    def find_lysis_time(self, r: np.ndarray, current_time: float, count: int) -> np.ndarray:
        lysis_r = self.rng.random(count, dtype=np.double) * (self.exp.macro_params.microscale_runs / 100)
        interp = np.full(count, float('inf'), dtype=np.double)
        r_idx = r.astype(int)
        total_lyses = self.exp.data.total_lyses[r_idx] - 1
        lysis_happens = lysis_r < total_lyses
        interp[~lysis_happens] = float('inf')
        for i in np.arange(count)[lysis_happens]:
            interp[i] = np.interp(lysis_r[i],
                                  np.arange(total_lyses[i]),
                                  self.exp.data.lysis_time[r_idx[i], :total_lyses[i]])
        interp = interp + (current_time - self.exp.macro_params.time_step / 2)

        return interp

    def run(self):
        for ts in tqdm(np.arange(self.exp.macro_params.total_time_steps)):
            current_time = ts * self.exp.macro_params.time_step

            self.m_fiber_status = self.fiber_status[self.location]

            self.unbind_by_degradation(self.bound & (self.m_fiber_status < current_time), current_time)
            self.unbind_by_time(self.bound & (self.binding_time < current_time), current_time)

            should_bind = (~self.bound
                           & (self.binding_time < current_time)
                           & (self.m_fiber_status > current_time)
                           & (self.waiting_time < current_time))

            move_chance = np.ones(self.exp.macro_params.total_molecules, dtype=np.double)
            move_chance[~self.bound] = self.rng.random(np.count_nonzero(~self.bound), dtype=np.double)
            should_move = move_chance < self.exp.macro_params.moving_probability
            conflict = should_bind & should_move
            # threshold = np.full(self.exp.macro_params.total_molecules, -1, dtype=np.double)
            threshold = ((current_time - self.binding_time[conflict])
                         / self.exp.macro_params.time_step)
            should_bind[conflict] = self.rng.random(np.count_nonzero(conflict)) <= threshold
            should_move[conflict] = ~should_bind[conflict]

            self.bind(should_bind, current_time)
            self.move(should_move, current_time)

            if ts % 100000 == 100000-1:
                unlysed_fibers = np.count_nonzero(self.fiber_status > current_time)

                if unlysed_fibers == 0:
                    print()
                    print(f"All fibers degraded after {current_time:.2f} sec. Terminating")
                    # break
                else:
                    unlysed_fiber_percent = 100 - unlysed_fibers / self.exp.macro_params.total_fibers * 100
                    reached_back_row_percent = self.number_reached_back_row / self.exp.macro_params.total_molecules * 100
                    print()
                    print(f"After {current_time:.2f} sec, {self.exp.macro_params.total_fibers - unlysed_fibers:,} "
                          f"fibers are degraded ({unlysed_fiber_percent:.1f}% of total) and "
                          f"{self.number_reached_back_row:,} molecules have reached the back row "
                          f"({reached_back_row_percent:.2f}% of total).")

        print(f"Total binds: {self.total_binds:,}")
        print(f"Timesteps with changes to degrade time: {self.timesteps_with_fiber_changes:,}")
        print(f"Total regular moves: {self.total_regular_moves:,}")
        print(f"Total restricted moves: {self.total_restricted_moves:,}")
        print(f"Total macro unbinds: {self.total_macro_unbinds:,}")
        print(f"Total micro unbinds: {self.total_micro_unbinds:,}")
        print(f"Molecules which reached the back row: {self.number_reached_back_row:,}")
