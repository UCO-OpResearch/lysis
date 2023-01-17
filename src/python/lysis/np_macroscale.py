import logging
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

        self.logger = logging.getLogger(__name__)
        self.logger.debug(f"Initializing MacroscaleRun")
        self.rng = np.random.default_rng(seed=abs(exp.macro_params.seed))

        self.binding_time_factory = self._BindingTimeFactory(self.exp, self.rng)

        edge_lookup = partial(np.ravel_multi_index, dims=(exp.macro_params.rows, exp.macro_params.full_row))

        self.fiber_status = np.full(self.exp.macro_params.rows * self.exp.macro_params.full_row,
                                    float('inf'), dtype=np.float_)
        for i, j in np.ndindex(self.exp.macro_params.empty_rows, self.exp.macro_params.full_row):
            self.fiber_status[edge_lookup((i, j))] = 0
        for j in range(exp.macro_params.cols):
            self.fiber_status[edge_lookup((self.exp.macro_params.rows-1, 3*j,))] = 0

        self.logger.debug(f"Precalculating neighbors.")
        self.neighbors = EdgeGrid.generate_neighborhood_structure(exp)

        self.logger.info(f"Placing molecules on empty edges.")
        location_i = self.rng.integers(exp.macro_params.empty_rows,
                                       size=exp.macro_params.total_molecules,
                                       dtype=np.short)
        location_j = self.rng.integers(exp.macro_params.full_row,
                                       size=exp.macro_params.total_molecules,
                                       dtype=np.short)

        self.location = edge_lookup((location_i, location_j))

        self.m_fiber_status = None

        self.bound = np.full(exp.macro_params.total_molecules, False, dtype=np.bool_)
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

        self.current_save_interval = 0

        self.data.saved_degradation_state = np.empty(())

        self.logger.debug(f"Initialization complete.")


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


    def unbind_by_degradation(self, m: np.ndarray, current_time: float):
        count = np.count_nonzero(m)
        if count == 0:
            return None
        self.bound = self.bound & ~m
        self.unbound_by_degradation = self.unbound_by_degradation | m
        self.total_macro_unbinds += count
        self.waiting_time[m] = (current_time
                                + self.exp.macro_params.average_bind_time
                                - self.exp.macro_params.time_step / 2)
        self.binding_time[m] = float('inf')

    def unbind_by_time(self, m: np.ndarray, current_time: float):
        count = np.count_nonzero(m)
        if count == 0:
            return None
        self.bound = self.bound & ~m
        self.unbound_by_degradation = self.unbound_by_degradation & ~m
        forced = np.full(self.exp.macro_params.total_molecules, False, dtype=np.bool_)
        forced[m] = (self.rng.random(count) <= self.exp.macro_params.forced_unbind)
        self.waiting_time[forced] = (current_time
                                     + self.exp.macro_params.average_bind_time
                                     - self.exp.macro_params.time_step / 2)
        self.binding_time[forced] = float('inf')
        num_forced = np.count_nonzero(forced)
        self.total_micro_unbinds += num_forced
        #
        if num_forced < count:
            self.binding_time[m & ~forced] = current_time + self.binding_time_factory.next(count - num_forced)

    def find_unbinding_time(self, unbinding_time_bin: np.ndarray, current_time: float) -> np.ndarray:
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
        interp = np.interp(unbinding_time_bin, self.xp, self.exp.data.unbinding_time)
        return interp + (current_time - self.exp.macro_params.time_step / 2)

    def find_lysis_time(self, unbinding_time_bin: np.ndarray, current_time: float, count: int) -> np.ndarray:
        lysis_time_bin = self.rng.random(count, dtype=np.double) * (self.exp.macro_params.microscale_runs / 100)
        interp = np.full(count, float('inf'), dtype=np.double)
        unbinding_time_bin = unbinding_time_bin.astype(int)
        total_lyses = self.exp.data.total_lyses[unbinding_time_bin] - 1
        lysis_happens = lysis_time_bin < total_lyses
        interp[~lysis_happens] = float('inf')

        # TODO(bpaynter): Should be able to improve this with a 2D interpolation or a custom numpy kernel
        for i in np.arange(count)[lysis_happens]:
            interp[i] = np.interp(lysis_time_bin[i],
                                  np.arange(total_lyses[i]),
                                  self.exp.data.lysis_time[unbinding_time_bin[i], :total_lyses[i]])
        return interp + (current_time - self.exp.macro_params.time_step / 2)

    def bind(self, m: np.ndarray, current_time: float):
        count = np.count_nonzero(m)
        if count == 0:
            return

        self.bound = self.bound | m
        self.waiting_time[m] = 0
        self.total_binds += count

        unbinding_time_bin: np.ndarray = self.rng.random(count) * 100
        self.binding_time[m] = self.find_unbinding_time(unbinding_time_bin, current_time)

        lysis_time = self.find_lysis_time(unbinding_time_bin, current_time, count)
        locations = self.location[m]
        if np.count_nonzero(self.fiber_status[locations] > lysis_time) > 0:
            self.timesteps_with_fiber_changes += 1

        self.fiber_status[locations] = np.fmin(self.fiber_status[locations],
                                               lysis_time)

    def move_to_empty_edge(self, m: np.ndarray, current_time: float):
        count = np.count_nonzero(m)
        if count == 0:
            return

        self.total_restricted_moves += count

        current_locations = self.location[m]

        neighborhoods = self.neighbors[current_locations]
        valid_neighbors = self.fiber_status[neighborhoods] < current_time
        valid_neighborhood_index = np.argsort(~valid_neighbors, axis=1)
        valid_neighborhoods = np.take_along_axis(neighborhoods, valid_neighborhood_index, axis=1)
        valid_neighborhoods = np.append(current_locations.reshape(count, 1),
                                        valid_neighborhoods,
                                        axis=1)
        num_valid_neighbors = np.count_nonzero(valid_neighbors, axis=1)
        neighbor = self.rng.random(count) * (num_valid_neighbors + 1)
        neighbor = neighbor.astype(int, copy=False)

        self.location[m] = valid_neighborhoods[np.full(count, True), neighbor]

    def find_still_stuck(self, m: np.ndarray, current_time: float):
        return ((self.waiting_time > current_time)
                & self.unbound_by_degradation
                & m)

    def unrestricted_move(self, free_to_move: np.ndarray, current_time: float):
        neighbor = self.rng.integers(8, size=np.count_nonzero(free_to_move))

        self.location[free_to_move] = self.neighbors[self.location[free_to_move], neighbor]
        self.total_regular_moves += np.count_nonzero(free_to_move)

        move_to_fiber = free_to_move & (self.m_fiber_status >= current_time)

        num_move_to_fiber = np.count_nonzero(move_to_fiber)
        if num_move_to_fiber > 0:
            self.binding_time[move_to_fiber] = current_time + self.binding_time_factory.next(num_move_to_fiber)

    def move(self, m: np.ndarray, current_time: float):
        still_stuck_to_fiber = self.find_still_stuck(m, current_time)
        self.move_to_empty_edge(still_stuck_to_fiber, current_time)

        free_to_move = m & ~still_stuck_to_fiber

        self.unrestricted_move(free_to_move, current_time)

        if self.number_reached_back_row < self.exp.macro_params.total_molecules:
            first_time = (~self.reached_back_row
                          & (self.location > (self.exp.macro_params.rows - 1) * self.exp.macro_params.full_row - 1))
            self.time_to_reach_back_row[first_time] = current_time
            self.reached_back_row = self.reached_back_row | first_time
            self.number_reached_back_row += np.count_nonzero(first_time)

    def run(self):
        for ts in tqdm(np.arange(self.exp.macro_params.total_time_steps), mininterval=2):
            current_time = ts * self.exp.macro_params.time_step

            self.m_fiber_status = self.fiber_status[self.location]

            self.unbind_by_degradation(self.bound & (self.m_fiber_status < current_time), current_time)
            self.unbind_by_time(self.bound & (self.binding_time < current_time), current_time)

            should_bind = (~self.bound
                           & (self.binding_time < current_time)
                           & (self.m_fiber_status > current_time)
                           & (self.waiting_time < current_time))

            move_chance = self.rng.random(self.exp.macro_params.total_molecules, dtype=np.double)
            should_move = (move_chance < self.exp.macro_params.moving_probability) & ~self.bound

            # TODO(bpaynter): Since these are all boolean vectors,
            #                 I should be able to use & and | instead of the lookup masks
            conflict = should_bind & should_move
            threshold = ((current_time - self.binding_time[conflict])
                         / self.exp.macro_params.time_step)
            should_bind[conflict] = self.rng.random(np.count_nonzero(conflict)) <= threshold
            should_move[conflict] = ~should_bind[conflict]

            self.bind(should_bind, current_time)
            self.move(should_move, current_time)

            if ts % 100000 == 100000-1:
                unlysed_fibers = np.count_nonzero(self.fiber_status > current_time)

                if unlysed_fibers == 0:
                    self.logger.info(f"All fibers degraded after {current_time:.2f} sec. Terminating")
                    # break
                else:
                    unlysed_fiber_percent = 100 - unlysed_fibers / self.exp.macro_params.total_fibers * 100
                    reached_back_row_percent = (self.number_reached_back_row
                                                / self.exp.macro_params.total_molecules
                                                * 100)
                    self.logger.info(f"After {current_time:.2f} sec, "
                                     f"{self.exp.macro_params.total_fibers - unlysed_fibers:,} "
                                     f"fibers are degraded ({unlysed_fiber_percent:.1f}% of total) and "
                                     f"{self.number_reached_back_row:,} molecules have reached the back row "
                                     f"({reached_back_row_percent:.1f}% of total).")

        self.logger.info(f"Total binds: {self.total_binds:,}")
        self.logger.info(f"Timesteps with changes to degrade time: {self.timesteps_with_fiber_changes:,}")
        self.logger.info(f"Total regular moves: {self.total_regular_moves:,}")
        self.logger.info(f"Total restricted moves: {self.total_restricted_moves:,}")
        self.logger.info(f"Total macro unbinds: {self.total_macro_unbinds:,}")
        self.logger.info(f"Total micro unbinds: {self.total_micro_unbinds:,}")
        self.logger.info(f"Molecules which reached the back row: {self.number_reached_back_row:,}")
