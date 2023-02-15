import logging
import os
from functools import partial

import numpy as np
from tqdm.auto import tqdm

from .util import Experiment, KissRandomGenerator, RandomDraw
from .edge_grid import EdgeGrid, from_fortran_edge_index, to_fortran_edge_index

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2022, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


class MacroscaleRun:
    def __init__(self, exp: Experiment, instance: int = None, seed: int = None):
        self.exp = exp
        assert self.exp.macro_params is not None

        self.logger = logging.getLogger(__name__)
        self.logger.debug(f"Initializing MacroscaleRun")
        if seed is None:
            seed = exp.macro_params.seed
        if exp.macro_params.duplicate_fortran:
            self.rng = KissRandomGenerator(seed)
        else:
            self.rng = np.random.default_rng(seed=abs(seed))

        self.binding_time_factory = self._BindingTimeFactory(self.exp, self.rng)

        self.edge_lookup = partial(
            np.ravel_multi_index,
            dims=(exp.macro_params.rows, exp.macro_params.full_row),
        )

        self.fiber_status = np.full(
            self.exp.macro_params.rows * self.exp.macro_params.full_row,
            float("inf"),
            dtype=np.float_,
        )
        self.real_fiber = np.full(
            self.exp.macro_params.rows * self.exp.macro_params.full_row,
            True,
            dtype=np.bool_,
        )
        for i, j in np.ndindex(
            self.exp.macro_params.empty_rows, self.exp.macro_params.full_row
        ):
            self.real_fiber[self.edge_lookup((i, j))] = False
        for j in range(exp.macro_params.cols):
            self.real_fiber[
                self.edge_lookup(
                    (
                        self.exp.macro_params.rows - 1,
                        3 * j,
                    )
                )
            ] = False
        self.fiber_status[~self.real_fiber] = 0

        self.logger.debug(f"Precalculating neighbors.")
        self.neighbors = EdgeGrid.generate_neighborhood_structure(exp)

        if exp.macro_params.duplicate_fortran:
            perm = np.array(
                [
                    from_fortran_edge_index(
                        k, exp.macro_params.rows, exp.macro_params.cols
                    )
                    for k in range(exp.macro_params.total_edges)
                ]
            )
            perm = np.ravel_multi_index(
                (perm[:, 0], perm[:, 1]),
                dims=(exp.macro_params.rows, exp.macro_params.full_row),
            )
            inv_perm = []
            for k in range(exp.macro_params.rows * exp.macro_params.full_row):
                p = np.argwhere(perm == k)
                if p.size > 0:
                    inv_perm.append(p[0, 0])
                else:
                    inv_perm.append(-1)
            inv_perm = np.array(inv_perm)
            sort = np.argsort(inv_perm[self.neighbors[perm]], axis=1)
            self.neighbors[perm] = np.take_along_axis(
                self.neighbors[perm], sort, axis=1
            )

        self.logger.info(f"Placing molecules on empty edges.")
        if exp.macro_params.duplicate_fortran:
            location = self.rng.random(exp.macro_params.total_molecules)
            location = (
                exp.macro_params.empty_rows * exp.macro_params.full_row * location
            )
            location = location.astype(int, copy=False)
            location_i = np.empty(exp.macro_params.total_molecules, dtype=np.int_)
            location_j = np.empty(exp.macro_params.total_molecules, dtype=np.int_)
            for m in range(len(location)):
                location_i[m], location_j[m] = from_fortran_edge_index(
                    location[m], exp.macro_params.rows, exp.macro_params.cols
                )
        else:
            location_i = self.rng.integers(
                exp.macro_params.empty_rows,
                size=exp.macro_params.total_molecules,
                dtype=np.short,
            )
            location_j = self.rng.integers(
                exp.macro_params.full_row,
                size=exp.macro_params.total_molecules,
                dtype=np.short,
            )

        self.location = self.edge_lookup((location_i, location_j))

        m_location_ij = np.unravel_index(
            self.location, (self.exp.macro_params.rows, self.exp.macro_params.full_row)
        )

        f_location = np.array(
            [
                to_fortran_edge_index(
                    m_location_ij[0][k],
                    m_location_ij[1][k],
                    self.exp.macro_params.rows,
                    self.exp.macro_params.cols,
                )
                for k in range(self.exp.macro_params.total_molecules)
            ]
        )

        self.m_fiber_status = None

        self.bound = np.full(exp.macro_params.total_molecules, False, dtype=np.bool_)
        self.waiting_time = np.full(
            exp.macro_params.total_molecules, 0, dtype=np.float_
        )
        self.binding_time = np.full(
            exp.macro_params.total_molecules, float("inf"), dtype=np.float_
        )
        self.unbound_by_degradation = np.full(
            exp.macro_params.total_molecules, 0, dtype=np.bool_
        )
        self.time_to_reach_back_row = np.full(
            exp.macro_params.total_molecules, float("inf"), dtype=np.float_
        )
        self.reached_back_row = np.full(
            exp.macro_params.total_molecules, False, dtype=np.bool_
        )
        self.xp = np.arange(101)
        self.random_numbers = None

        self.total_macro_unbinds = 0
        self.total_micro_unbinds = 0
        self.total_binds = 0
        self.independent_binds = 0
        self.total_regular_moves = 0
        self.total_restricted_moves = 0
        self.timesteps_with_fiber_changes = 0
        self.number_reached_back_row = 0
        self.last_degrade_time = float("inf")

        self.current_save_interval = 0

        self.exp.data.degradation_state = np.empty(
            (
                self.exp.macro_params.number_of_saves,
                self.exp.macro_params.rows * self.exp.macro_params.full_row,
            ),
            dtype=np.float_,
        )
        self.exp.data.molecule_location = np.empty(
            (
                self.exp.macro_params.number_of_saves,
                self.exp.macro_params.total_molecules,
            ),
            dtype=np.int_,
        )
        self.exp.data.molecule_state = np.empty(
            (
                self.exp.macro_params.number_of_saves,
                self.exp.macro_params.total_molecules,
            ),
            dtype=np.bool_,
        )
        self.exp.data.save_time = np.empty(
            (self.exp.macro_params.number_of_saves,),
            dtype=np.float_,
        )

        self.degrade_time_tracker = []

        # self.m_tracker_index = 21
        # self.m_tracker_last = self.location[self.m_tracker_index]
        # self.m_tracker = [(0, self.location[self.m_tracker_index])]

        self.logger.debug(f"Initialization complete.")

    class _BindingTimeFactory:
        def __init__(self, exp: Experiment, rng: np.random.Generator):
            self.period = max(1000000, 10 * exp.macro_params.total_molecules)
            self.exp = exp
            self.rng = rng
            self.pointer = self.period
            self.list = np.empty((self.period,), dtype=np.double)

        def fill_list(
            self,
            random_numbers: np.ndarray = None,
        ) -> np.ndarray | None:
            if random_numbers is None:
                binding_time_list = self.rng.random(
                    len(self.list)
                )  # dtype=np.double, out=self.list)
                return_list = False
            else:
                binding_time_list = random_numbers
                return_list = True
            binding_time_list = np.log(binding_time_list, out=binding_time_list)
            denominator = (
                self.exp.macro_params.binding_rate * self.exp.macro_params.binding_sites
            )
            binding_time_list = np.divide(
                binding_time_list, denominator, out=binding_time_list
            )
            binding_time_list = np.subtract(
                -self.exp.macro_params.time_step / 2,
                binding_time_list,
                out=binding_time_list,
            )
            if return_list:
                return binding_time_list
            else:
                self.list = binding_time_list

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
                    out = self.list[self.pointer : self.pointer + count]
                    self.pointer += count - 1
                else:
                    out[: self.period - self.pointer] = self.list[self.pointer :]
                    self.fill_list()
                    out[self.period - self.pointer :] = self.list[
                        : self.pointer + count - self.period
                    ]
                    self.pointer += count - self.period - 1
                return out

    def unbind_by_degradation(self, m: np.ndarray, current_time: float):
        count = np.count_nonzero(m)
        if count == 0:
            return None
        self.bound = self.bound & ~m
        self.unbound_by_degradation = self.unbound_by_degradation | m
        self.total_macro_unbinds += count
        self.waiting_time[m] = (
            current_time
            + self.exp.macro_params.average_bind_time
            - self.exp.macro_params.time_step / 2
        )
        self.binding_time[m] = float("inf")

    def unbind_by_time(self, m: np.ndarray, current_time: float):
        count = np.count_nonzero(m)
        if count == 0:
            return None
        self.bound = self.bound & ~m
        self.unbound_by_degradation = self.unbound_by_degradation & ~m
        forced = np.full(self.exp.macro_params.total_molecules, False, dtype=np.bool_)
        if self.exp.macro_params.duplicate_fortran:
            forced[m] = (
                self.random_numbers[RandomDraw.MICRO_UNBIND][m]
                <= self.exp.macro_params.forced_unbind
            )
        else:
            forced[m] = self.rng.random(count) <= self.exp.macro_params.forced_unbind
        self.waiting_time[forced] = (
            current_time
            + self.exp.macro_params.average_bind_time
            - self.exp.macro_params.time_step / 2
        )
        self.binding_time[forced] = float("inf")
        num_forced = np.count_nonzero(forced)
        self.total_micro_unbinds += num_forced
        if num_forced < count:
            if self.exp.macro_params.duplicate_fortran:
                self.binding_time[
                    m & ~forced
                ] = current_time + self.binding_time_factory.fill_list(
                    self.random_numbers[RandomDraw.BINDING_TIME_WHEN_UNBINDING][
                        m & ~forced
                    ]
                )
                # c = np.empty((3 * (count - num_forced),), dtype=np.float_)
                # c[0::3] = self.random_numbers[RandomDraw.BINDING_TIME_WHEN_UNBINDING][
                #     m & ~forced
                # ]
                # c[1::3] = current_time
                # c[2::3] = self.binding_time[m & ~forced]
                # self.bind_info_file.write(np.array2string(c)[1:-1] + os.linesep)
            else:
                self.binding_time[
                    m & ~forced
                ] = current_time + self.binding_time_factory.next(count - num_forced)

    def find_unbinding_time(
        self, unbinding_time_bin: np.ndarray, current_time: float
    ) -> np.ndarray:
        """Find the experiment time at which the tPA molecule will unbind

        If we think of the binding_time matrix as a function f(x) where
        binding_time[100*x] = f(x), we want to draw uniformly from the range of
        f, that is, we want the unbinding timestep to be f(r) where r~U(0,1).

        If -- magically -- 100 * r is an integer i, this is easy because we
        just look up binding_time[i], but what if it isn't a perfect integer?

        Then f(r) lies in the interval ( f(floor(100*r)), f(ceil(100*r)) ) so
        we interpolate it linearly. That is, define a linear function g(x) such
        that
        g(floor(100*r)) = f(floor(100*r)) and g(ceil(100*r)) = f(ceil(100*r))
        and then use g(r) to approximate f(r).

        Another way to see this is,
        if i is an integer such that 100i < 100r < 100(i+1) and
        lambda is such that 100i + lambda = 100r,
        then f(r) ~ (1-lambda)*f(i/100) + lambda*f( (i+1)/100 )"""
        interp = np.interp(unbinding_time_bin, self.xp, self.exp.data.unbinding_time)
        return interp + (current_time - self.exp.macro_params.time_step / 2)

    def find_lysis_time(
        self,
        m: np.ndarray,
        unbinding_time_bin: np.ndarray,
        current_time: float,
        count: int,
    ) -> np.ndarray:
        if self.exp.macro_params.duplicate_fortran:
            lysis_time_bin = self.random_numbers[RandomDraw.LYSIS_TIME][m]
        else:
            lysis_time_bin = self.rng.random(count)
        lysis_time_bin = lysis_time_bin * (self.exp.macro_params.microscale_runs / 100)
        interp = np.full(count, float("inf"), dtype=np.double)
        unbinding_time_bin = unbinding_time_bin.astype(int)
        total_lyses = self.exp.data.total_lyses[unbinding_time_bin] - 1
        lysis_happens = lysis_time_bin < total_lyses
        interp[~lysis_happens] = float("inf")

        mol_index = np.arange(self.exp.macro_params.total_molecules)
        mol_index = mol_index[m]

        # TODO(bpaynter): Should be able to improve this with a 2D
        #                 interpolation or a custom numpy kernel
        for i in np.arange(count)[lysis_happens]:
            interp[i] = np.interp(
                lysis_time_bin[i],
                np.arange(total_lyses[i]),
                self.exp.data.lysis_time[unbinding_time_bin[i], : total_lyses[i]],
            )
            if (
                interp[i] + (current_time - self.exp.macro_params.time_step / 2)
                < self.fiber_status[self.location[mol_index[i]]]
            ):
                self.degrade_time_tracker.append(
                    (
                        current_time,
                        self.location[mol_index[i]],
                        self.random_numbers[RandomDraw.LYSIS_TIME][mol_index[i]],
                        interp[i],
                        interp[i]
                        + (current_time - self.exp.macro_params.time_step / 2),
                    )
                )
        return interp + (current_time - self.exp.macro_params.time_step / 2)

    def bind(self, m: np.ndarray, current_time: float):
        count = np.count_nonzero(m)
        if count == 0:
            return

        self.bound = self.bound | m
        self.waiting_time[m] = 0
        self.total_binds += count

        if self.exp.macro_params.duplicate_fortran:
            unbinding_time_bin = self.random_numbers[RandomDraw.UNBINDING_TIME][m] * 100
        else:
            unbinding_time_bin = self.rng.random(count) * 100
        self.binding_time[m] = self.find_unbinding_time(
            unbinding_time_bin, current_time
        )

        lysis_time = self.find_lysis_time(m, unbinding_time_bin, current_time, count)
        locations = self.location[m]
        for i in range(count):
            if locations[i] == 7794:
                if self.fiber_status[locations[i]] > lysis_time[i]:
                    pass
            if lysis_time[i] < float("inf"):
                self.fiber_status[locations[i]] = min(
                    self.fiber_status[locations[i]], lysis_time[i]
                )
        # NOTE:
        # This line should NOT be vectorized as below. If two (or more) molecules bind
        # to the same fiber on the same timestep, the second molecule's lysis time will
        # override the first one's, even if the first one's lysis time is lower!
        #
        # self.fiber_status[locations] = np.fmin(self.fiber_status[locations], lysis_time)

    def move_to_empty_edge(self, m: np.ndarray, current_time: float):
        count = np.count_nonzero(m)
        if count == 0:
            return

        self.total_restricted_moves += count

        current_locations = self.location[m]

        neighborhoods = self.neighbors[current_locations]
        valid_neighbors = self.fiber_status[neighborhoods] < current_time
        valid_neighborhood_index = np.argsort(~valid_neighbors, axis=1)
        valid_neighborhoods = np.take_along_axis(
            neighborhoods, valid_neighborhood_index, axis=1
        )
        valid_neighborhoods = np.append(
            current_locations.reshape(count, 1), valid_neighborhoods, axis=1
        )
        num_valid_neighbors = np.count_nonzero(valid_neighbors, axis=1)
        if self.exp.macro_params.duplicate_fortran:
            neighbor = self.random_numbers[RandomDraw.RESTRICTED_MOVE][m] * (
                num_valid_neighbors + 1
            )
        else:
            neighbor = self.rng.random(count) * (num_valid_neighbors + 1)
        neighbor = neighbor.astype(int, copy=False)

        self.location[m] = valid_neighborhoods[np.full(count, True), neighbor]

    def find_still_stuck(self, m: np.ndarray, current_time: float):
        return (self.waiting_time > current_time) & self.unbound_by_degradation & m

    def unrestricted_move(self, free_to_move: np.ndarray, current_time: float):
        if self.exp.macro_params.duplicate_fortran:
            neighbor = self.random_numbers[RandomDraw.MOVE][free_to_move]
            neighbor = neighbor - (1 - self.exp.macro_params.moving_probability)
            neighbor = neighbor / self.exp.macro_params.moving_probability
            neighbor = neighbor * 8
            neighbor = neighbor.astype(int, copy=False)
        else:
            neighbor = self.rng.integers(8, size=np.count_nonzero(free_to_move))

        self.location[free_to_move] = self.neighbors[
            self.location[free_to_move], neighbor
        ]
        self.total_regular_moves += np.count_nonzero(free_to_move)

        # move_to_fiber = free_to_move & (self.m_fiber_status >= current_time)

        num_move_to_fiber = np.count_nonzero(free_to_move)
        if num_move_to_fiber > 0:
            if self.exp.macro_params.duplicate_fortran:
                self.binding_time[
                    free_to_move
                ] = current_time + self.binding_time_factory.fill_list(
                    self.random_numbers[RandomDraw.BINDING_TIME_WHEN_MOVING][
                        free_to_move
                    ]
                )
                # c = np.empty((3 * num_move_to_fiber,), dtype=np.float_)
                # c[0::3] = self.random_numbers[RandomDraw.BINDING_TIME_WHEN_MOVING][
                #     move_to_fiber
                # ]
                # c[1::3] = current_time
                # c[2::3] = self.binding_time[move_to_fiber]
                # self.bind_info_file.write(np.array2string(c)[1:-1] + os.linesep)
            else:
                self.binding_time[
                    free_to_move
                ] = current_time + self.binding_time_factory.next(num_move_to_fiber)

    def move(self, m: np.ndarray, current_time: float):
        still_stuck_to_fiber = self.find_still_stuck(m, current_time)
        self.move_to_empty_edge(still_stuck_to_fiber, current_time)

        free_to_move = m & ~still_stuck_to_fiber

        self.unrestricted_move(free_to_move, current_time)

        if self.number_reached_back_row < self.exp.macro_params.total_molecules:
            first_time = ~self.reached_back_row & (
                self.location
                > (self.exp.macro_params.rows - 1) * self.exp.macro_params.full_row - 1
            )
            self.time_to_reach_back_row[first_time] = current_time
            self.reached_back_row = self.reached_back_row | first_time
            self.number_reached_back_row += np.count_nonzero(first_time)

    def save_data(self, current_time):
        # self.logger.info(f"Saving data at time {current_time:.2f} sec.")
        self.exp.data.degradation_state[self.current_save_interval] = self.fiber_status
        self.exp.data.molecule_location[self.current_save_interval] = self.location
        self.exp.data.molecule_state[self.current_save_interval] = self.bound
        self.exp.data.save_time[self.current_save_interval] = current_time
        self.current_save_interval += 1

    def record_data_to_disk(self):
        self.logger.info(f"Saving data to disk.")
        self.exp.data.save_to_disk("degradation_state")
        self.exp.data.save_to_disk("molecule_location")
        self.exp.data.save_to_disk("molecule_state")
        self.exp.data.save_to_disk("save_time")

    def run(self):
        self.save_data(0)
        for ts in tqdm(
            np.arange(self.exp.macro_params.total_time_steps), mininterval=2
        ):
            current_time = ts * self.exp.macro_params.time_step
            if self.exp.macro_params.duplicate_fortran:
                current_time += self.exp.macro_params.time_step

            if ts == 23:
                pass

            if self.exp.macro_params.duplicate_fortran:
                self.random_numbers = np.empty(
                    (8, self.exp.macro_params.total_molecules), np.float_
                )
                for i in range(8):
                    self.random_numbers[i] = self.rng.random(
                        self.exp.macro_params.total_molecules
                    )

            self.m_fiber_status = self.fiber_status[self.location]

            self.unbind_by_degradation(
                self.bound & (self.m_fiber_status < current_time), current_time
            )
            self.unbind_by_time(
                self.bound & (self.binding_time < current_time), current_time
            )

            should_bind = (
                ~self.bound
                & (self.binding_time < current_time)
                & (self.m_fiber_status > current_time)
                & (self.waiting_time < current_time)
            )
            if self.exp.macro_params.duplicate_fortran:
                move_chance = self.random_numbers[RandomDraw.MOVE]
                should_move = (
                    move_chance > 1 - self.exp.macro_params.moving_probability
                ) & ~self.bound
            else:
                move_chance = self.rng.random(self.exp.macro_params.total_molecules)
                should_move = (
                    move_chance < self.exp.macro_params.moving_probability
                ) & ~self.bound

            # TODO(bpaynter): Since these are all boolean vectors, I should be
            #                 able to use & and | instead of the lookup masks
            conflict = should_bind & should_move
            threshold = (
                current_time - self.binding_time[conflict]
            ) / self.exp.macro_params.time_step
            if self.exp.macro_params.duplicate_fortran:
                should_bind[conflict] = (
                    self.random_numbers[RandomDraw.CONFLICT_RESOLUTION][conflict]
                    <= threshold
                )
            else:
                should_bind[conflict] = (
                    self.rng.random(np.count_nonzero(conflict)) <= threshold
                )
            should_move[conflict] = ~should_bind[conflict]

            self.bind(should_bind, current_time)
            self.move(should_move, current_time)

            if (
                current_time
                >= self.exp.macro_params.save_interval * self.current_save_interval
            ):
                self.save_data(current_time)

            if ts % 100000 == 100000 - 1:
                unlysed_fibers = np.count_nonzero(self.fiber_status > current_time)

                if unlysed_fibers == 0:
                    self.logger.info(
                        f"All fibers degraded after {current_time:.2f} sec. Terminating"
                    )
                    self.last_degrade_time = np.max(self.fiber_status)
                    # break
                else:
                    unlysed_fiber_percent = (
                        100 - unlysed_fibers / self.exp.macro_params.total_fibers * 100
                    )
                    reached_back_row_percent = (
                        self.number_reached_back_row
                        / self.exp.macro_params.total_molecules
                        * 100
                    )
                    self.logger.info(
                        f"After {current_time:.2f} sec, "
                        f"{self.exp.macro_params.total_fibers - unlysed_fibers:,} "
                        f"fibers are degraded ({unlysed_fiber_percent:.1f}% of total) "
                        f"and {self.number_reached_back_row:,} molecules have reached "
                        f"the back row {reached_back_row_percent:.1f}% of total)."
                    )

            # m_location_ij = np.unravel_index(
            #     self.location,
            #     (self.exp.macro_params.rows, self.exp.macro_params.full_row),
            # )
            # f_location = np.array(
            #     [
            #         to_fortran_edge_index(
            #             m_location_ij[0][k],
            #             m_location_ij[1][k],
            #             self.exp.macro_params.rows,
            #             self.exp.macro_params.cols,
            #         )
            #         for k in range(self.exp.macro_params.total_molecules)
            #     ]
            # )
            #
            # if self.location[self.m_tracker_index] != self.m_tracker_last:
            #     self.m_tracker.append((ts, self.location[self.m_tracker_index]))
            #     self.m_tracker_last = self.location[self.m_tracker_index]
            # pass

        degrade_time_tracker = np.array(self.degrade_time_tracker)
        np.save(
            os.path.join(self.exp.os_path, "deg_tracker.p.npy"), degrade_time_tracker
        )
        self.save_data(self.exp.macro_params.total_time)
        self.record_data_to_disk()

        self.logger.info(f"Total binds: {self.total_binds:,}")
        self.logger.info(
            f"Timesteps with changes to degrade time: "
            f"{self.timesteps_with_fiber_changes:,}"
        )
        self.logger.info(f"Total regular moves: {self.total_regular_moves:,}")
        self.logger.info(f"Total restricted moves: {self.total_restricted_moves:,}")
        self.logger.info(f"Total macro unbinds: {self.total_macro_unbinds:,}")
        self.logger.info(f"Total micro unbinds: {self.total_micro_unbinds:,}")
        self.logger.info(
            f"Molecules which reached the back row: {self.number_reached_back_row:,}"
        )
        self.logger.info(f"Last fiber degraded at: {self.last_degrade_time:2f} sec")
