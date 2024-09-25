import logging
import os
from functools import partial

import numpy as np
from tqdm.auto import tqdm

from .util import (
    Run,
    KissRandomGenerator,
    RandomDraw,
    EdgeGrid,
    from_fortran_edge_index,
    to_fortran_edge_index,
)


class MacroscaleSim:
    def __init__(self, run: Run, instance: int = None, seed: int = None):
        self.run = run
        assert self.run.macro_params is not None

        self.logger = logging.getLogger(__name__)
        self.logger.debug(f"Initializing MacroscaleSim")
        if seed is None:
            seed = run.macro_params.seed
        if run.macro_params.duplicate_fortran:
            self.rng = KissRandomGenerator(seed)
        else:
            self.rng = np.random.default_rng(seed=abs(seed))

        self.binding_time_factory = self._BindingTimeFactory(self.run, self.rng)

        self.edge_lookup = partial(
            np.ravel_multi_index,
            dims=(run.macro_params.rows, run.macro_params.full_row),
        )

        self.fiber_status = np.full(
            self.run.macro_params.rows * self.run.macro_params.full_row,
            float("inf"),
            dtype=np.float_,
        )
        self.real_fiber = np.full(
            self.run.macro_params.rows * self.run.macro_params.full_row,
            True,
            dtype=np.bool_,
        )
        for i, j in np.ndindex(
            self.run.macro_params.empty_rows, self.run.macro_params.full_row
        ):
            self.real_fiber[self.edge_lookup((i, j))] = False
        for j in range(run.macro_params.cols):
            self.real_fiber[
                self.edge_lookup(
                    (
                        self.run.macro_params.rows - 1,
                        3 * j,
                    )
                )
            ] = False
        self.fiber_status[~self.real_fiber] = 0

        self.logger.debug(f"Precalculating neighbors.")
        self.neighbors = EdgeGrid.generate_neighborhood_structure(run)

        if run.macro_params.duplicate_fortran:
            perm = np.array(
                [
                    from_fortran_edge_index(
                        k, run.macro_params.rows, run.macro_params.cols
                    )
                    for k in range(run.macro_params.total_edges)
                ]
            )
            perm = np.ravel_multi_index(
                (perm[:, 0], perm[:, 1]),
                dims=(run.macro_params.rows, run.macro_params.full_row),
            )
            inv_perm = []
            for k in range(run.macro_params.rows * run.macro_params.full_row):
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
        if run.macro_params.duplicate_fortran:
            location = self.rng.random(run.macro_params.total_molecules)
            location = (
                run.macro_params.empty_rows * run.macro_params.full_row * location
            )
            location = location.astype(int, copy=False)
            location_i = np.empty(run.macro_params.total_molecules, dtype=np.int_)
            location_j = np.empty(run.macro_params.total_molecules, dtype=np.int_)
            for m in range(len(location)):
                location_i[m], location_j[m] = from_fortran_edge_index(
                    location[m], run.macro_params.rows, run.macro_params.cols
                )
        else:
            location_i = self.rng.integers(
                run.macro_params.empty_rows,
                size=run.macro_params.total_molecules,
                dtype=np.short,
            )
            location_j = self.rng.integers(
                run.macro_params.full_row,
                size=run.macro_params.total_molecules,
                dtype=np.short,
            )

        self.location = self.edge_lookup((location_i, location_j))

        m_location_ij = np.unravel_index(
            self.location, (self.run.macro_params.rows, self.run.macro_params.full_row)
        )

        f_location = np.array(
            [
                to_fortran_edge_index(
                    m_location_ij[0][k],
                    m_location_ij[1][k],
                    self.run.macro_params.rows,
                    self.run.macro_params.cols,
                )
                for k in range(self.run.macro_params.total_molecules)
            ]
        )

        self.m_fiber_status = None

        self.bound = np.full(run.macro_params.total_molecules, False, dtype=np.bool_)
        self.waiting_time = np.full(
            run.macro_params.total_molecules, 0, dtype=np.float_
        )
        self.binding_time = np.full(
            run.macro_params.total_molecules, float("inf"), dtype=np.float_
        )
        self.unbound_by_degradation = np.full(
            run.macro_params.total_molecules, 0, dtype=np.bool_
        )
        self.time_to_reach_back_row = np.full(
            run.macro_params.total_molecules, float("inf"), dtype=np.float_
        )
        self.reached_back_row = np.full(
            run.macro_params.total_molecules, False, dtype=np.bool_
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

        self.run.data.degradation_state = np.empty(
            (
                self.run.macro_params.number_of_saves,
                self.run.macro_params.rows * self.run.macro_params.full_row,
            ),
            dtype=np.float_,
        )
        self.run.data.molecule_location = np.empty(
            (
                self.run.macro_params.number_of_saves,
                self.run.macro_params.total_molecules,
            ),
            dtype=np.int_,
        )
        self.run.data.molecule_state = np.empty(
            (
                self.run.macro_params.number_of_saves,
                self.run.macro_params.total_molecules,
            ),
            dtype=np.bool_,
        )
        self.run.data.save_time = np.empty(
            (self.run.macro_params.number_of_saves,),
            dtype=np.float_,
        )

        self.logger.debug(f"Initialization complete.")

    class _BindingTimeFactory:
        """Provides binding times from an internal array which is automatically
        refreshed as necessary

        This class is used internally to provide binding times to moving
        molecules. It pre-generates an array of fixed size and then serves from
        that array as requested, refilling the array when empty.
        Since generating one large array is far more efficient than generating
        many smaller ones, we are able to achieve better amortized performance.

        Args:
            run: The Run object that this model is a part of
            rng: A generator of random numbers. This should be the same
                    generator object that the rest of the model is using to
                    prevent overlap issues.
        """

        def __init__(self, run: Run, rng: np.random.Generator):
            # The size of the internal list
            self.period = max(1000000, 10 * run.macro_params.total_molecules)
            # The run this object is a part of
            self.run = run
            # The random number generator to be used
            self.rng = rng
            # Initialize the list. It starts empty and will be filled the first
            # time `next()` is called.
            self.list = np.empty((self.period,), dtype=np.double)
            # This points to the next binding time in the list to be served
            # We initialize it here at the end of the list as the list is
            # currently empty
            self.pointer = self.period

        def fill_list(
            self,
            random_numbers: np.ndarray = None,
        ) -> np.ndarray | None:
            """Fills an array with binding times.

            This method, when called with no arguments, will refill the array
            of binding times stored internally by this object. In this
            situation the method will return nothing.
            If an array of random numbers is passed, the method will return an
            array of binding times calculated using those random numbers. The
            returned array will have the same number of entries as the array
            passed. This usage will NOT benefit from the amortized performance
            acceleration.
            NOTE: Binding times generated by this object do NOT include the
            current time.

            Args:
                random_numbers: The array of random numbers that should be used
                    to generate new binding times. This argument should only be
                    used if NOT using the internal array

            Returns: An array of binding times, if an array was passed.
                Else, nothing.
            """
            # Check to see if we were passed any random numbers
            if random_numbers is None:
                # If not, we need to generate random numbers
                binding_time_list = self.rng.random(len(self.list))
                # Set a flag so that we know to use the internal list
                return_list = False
            else:
                # If we were passed random numbers, we use those
                binding_time_list = random_numbers
                # And set a flag so that we know to return the times we generate
                return_list = True

            # Binding time is calculated as:
            # (current_time) - ln(X)/(binding_rate * binding_sites) - time_step/2
            # where X ~ U(0,1)
            # In general, we don't know the current time when the binding time
            # is generated, so we depend on the calling code to add that.

            # To save processing time, we do as many operations 'in-place' as
            # possible.
            binding_time_list = np.log(binding_time_list, out=binding_time_list)
            # Since the denominator is the same for all terms, it is better to
            # calculate this once and re-use it.
            denominator = (
                self.run.macro_params.binding_rate * self.run.macro_params.binding_sites
            )
            # Use the vectorized, in-place divide function for the middle term
            binding_time_list = np.divide(
                binding_time_list, denominator, out=binding_time_list
            )
            # Using vectorized, in-place subtract to combine the last two terms
            # Using subtract so that we don't have to multiply
            # each 'ln(X)' term by -1. Don't know if this makes a performance
            # difference or not.
            binding_time_list = np.subtract(
                -self.run.macro_params.time_step / 2,
                binding_time_list,
                out=binding_time_list,
            )
            # Either return the binding times we calculated or store them in
            # the internal list.
            if return_list:
                return binding_time_list
            else:
                self.list = binding_time_list

        def next(self, count: int = 1):
            """Provide an array of binding times.

            If the internal array does not contain enough remaining binding
            times, it will automatically refill the array and provide from
            there. This means that the running time of this method will be
            inconsistent (fast most of the time, really slow whenever the list
            empty), but the overall (amortized) running time is fast.

            Args:
                count: the length of the array being requested

            Returns: An array of binding times of length 'count'
            """
            # Allocate the array for the output
            out = np.empty(count, dtype=np.double)
            # If there are enough binding times left in the internal array to
            # provide what is being requested
            if count + self.pointer <= self.period:
                # Fill the output array from the internal array
                out = self.list[self.pointer : self.pointer + count]
                # Update the pointer
                self.pointer += count
            # Otherwise, if there are not enough stored binding times
            else:
                # Put the remaining binding times at the start of the output
                # array
                out[: self.period - self.pointer] = self.list[self.pointer :]
                # Refill the internal array
                self.fill_list()
                # Fill the remainder of the output array from the new list
                out[self.period - self.pointer :] = self.list[
                    : self.pointer + count - self.period
                ]
                # Update the pointer
                self.pointer += count - self.period
            return out

    def unbind_by_degradation(self, m: np.ndarray, current_time: float):
        """Deals with molecules bound to a fiber that just degraded

        Args:
            m: A boolean array of molecules. Used as a mask to apply this
                function to.
            current_time: The current model time
        """
        # Count how many molecules need to be unbound
        # that is, how many elements of `m` are `True`
        count = np.count_nonzero(m)
        # If none, then return without any more work
        if count == 0:
            return
        # Set the bound status
        # True if we are currently bound and not selected for unbinding
        self.bound = self.bound & ~m
        # Set the 'macro unbind' flag
        # True if we are currently flagged, or were selected for unbinding
        self.unbound_by_degradation = self.unbound_by_degradation | m
        # Update the 'macro unbind' counter
        self.total_macro_unbinds += count
        # Set the waiting time of selected molecules to the current time,
        # plus the average bind time, minus half a timestep
        self.waiting_time[m] = (
            current_time
            + self.run.macro_params.average_bound_time
            - self.run.macro_params.time_step / 2
        )
        # Set the binding/unbinding time of the selected molecules to infinite
        # The fiber in their current location just degraded, so they can never
        # bind to it.
        self.binding_time[m] = float("inf")

    def unbind_by_time(self, m: np.ndarray, current_time: float):
        """Deals with bound molecules whose 'unbind' timer just expired.

        Args:
            m: A boolean array of molecules. Used as a mask to apply this
                function to.
            current_time: The current model time
        """
        # Count how many molecules need to be unbound
        # that is, how many elements of `m` are `True`
        count = np.count_nonzero(m)
        # If none, then return without any more work
        if count == 0:
            return
        self.bound = self.bound & ~m
        self.unbound_by_degradation = self.unbound_by_degradation & ~m
        forced = np.full(self.run.macro_params.total_molecules, False, dtype=np.bool_)
        if self.run.macro_params.duplicate_fortran:
            forced[m] = (
                self.random_numbers[RandomDraw.MICRO_UNBIND][m]
                <= self.run.macro_params.forced_unbind
            )
        else:
            forced[m] = self.rng.random(count) <= self.run.macro_params.forced_unbind
        self.waiting_time[forced] = (
            current_time
            + self.run.macro_params.average_bound_time
            - self.run.macro_params.time_step / 2
        )
        self.binding_time[forced] = float("inf")
        num_forced = np.count_nonzero(forced)
        self.total_micro_unbinds += num_forced
        if num_forced < count:
            if self.run.macro_params.duplicate_fortran:
                self.binding_time[m & ~forced] = (
                    current_time
                    + self.binding_time_factory.fill_list(
                        self.random_numbers[RandomDraw.BINDING_TIME_WHEN_UNBINDING][
                            m & ~forced
                        ]
                    )
                )
            else:
                self.binding_time[m & ~forced] = (
                    current_time + self.binding_time_factory.next(count - num_forced)
                )

    def find_unbinding_time(
        self, unbinding_time_bin: np.ndarray, current_time: float
    ) -> np.ndarray:
        """Find the simulation time at which the tPA molecule will unbind

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
        interp = np.interp(unbinding_time_bin, self.xp, self.run.data.unbinding_time)
        return interp + (current_time - self.run.macro_params.time_step / 2)

    def find_lysis_time(
        self,
        m: np.ndarray,
        unbinding_time_bin: np.ndarray,
        current_time: float,
        count: int,
    ) -> np.ndarray:
        if self.run.macro_params.duplicate_fortran:
            lysis_time_bin = self.random_numbers[RandomDraw.LYSIS_TIME][m]
        else:
            lysis_time_bin = self.rng.random(count)
        lysis_time_bin = lysis_time_bin * (self.run.macro_params.microscale_runs / 100)
        interp = np.full(count, float("inf"), dtype=np.double)
        unbinding_time_bin = unbinding_time_bin.astype(int)
        total_lyses = self.run.data.total_lyses[unbinding_time_bin] - 1
        lysis_happens = lysis_time_bin < total_lyses
        interp[~lysis_happens] = float("inf")

        # TODO(bpaynter): Should be able to improve this with a 2D
        #                 interpolation or a custom numpy kernel
        for i in np.arange(count)[lysis_happens]:
            interp[i] = np.interp(
                lysis_time_bin[i],
                np.arange(total_lyses[i]),
                self.run.data.lysis_time[unbinding_time_bin[i], : total_lyses[i]],
            )
        return interp + (current_time - self.run.macro_params.time_step / 2)

    def bind(self, m: np.ndarray, current_time: float):
        count = np.count_nonzero(m)
        if count == 0:
            return

        self.bound = self.bound | m
        self.waiting_time[m] = 0
        self.total_binds += count

        if self.run.macro_params.duplicate_fortran:
            unbinding_time_bin = self.random_numbers[RandomDraw.UNBINDING_TIME][m] * 100
        else:
            unbinding_time_bin = self.rng.random(count) * 100
        self.binding_time[m] = self.find_unbinding_time(
            unbinding_time_bin, current_time
        )

        lysis_time = self.find_lysis_time(m, unbinding_time_bin, current_time, count)
        locations = self.location[m]
        for i in range(count):
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
        if self.run.macro_params.duplicate_fortran:
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
        if self.run.macro_params.duplicate_fortran:
            neighbor = self.random_numbers[RandomDraw.MOVE][free_to_move]
            neighbor = neighbor - (1 - self.run.macro_params.moving_probability)
            neighbor = neighbor / self.run.macro_params.moving_probability
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
            if self.run.macro_params.duplicate_fortran:
                self.binding_time[free_to_move] = (
                    current_time
                    + self.binding_time_factory.fill_list(
                        self.random_numbers[RandomDraw.BINDING_TIME_WHEN_MOVING][
                            free_to_move
                        ]
                    )
                )
            else:
                self.binding_time[free_to_move] = (
                    current_time + self.binding_time_factory.next(num_move_to_fiber)
                )

    def move(self, m: np.ndarray, current_time: float):
        still_stuck_to_fiber = self.find_still_stuck(m, current_time)
        self.move_to_empty_edge(still_stuck_to_fiber, current_time)

        free_to_move = m & ~still_stuck_to_fiber

        self.unrestricted_move(free_to_move, current_time)

        if self.number_reached_back_row < self.run.macro_params.total_molecules:
            first_time = ~self.reached_back_row & (
                self.location
                > (self.run.macro_params.rows - 1) * self.run.macro_params.full_row - 1
            )
            self.time_to_reach_back_row[first_time] = current_time
            self.reached_back_row = self.reached_back_row | first_time
            self.number_reached_back_row += np.count_nonzero(first_time)

    def save_data(self, current_time):
        self.run.data.degradation_state[self.current_save_interval] = self.fiber_status
        self.run.data.molecule_location[self.current_save_interval] = self.location
        self.run.data.molecule_state[self.current_save_interval] = self.bound
        self.run.data.save_time[self.current_save_interval] = current_time
        self.current_save_interval += 1

    def record_data_to_disk(self):
        self.logger.info(f"Saving data to disk.")
        self.run.data.save_to_disk("degradation_state")
        self.run.data.save_to_disk("molecule_location")
        self.run.data.save_to_disk("molecule_state")
        self.run.data.save_to_disk("save_time")

    def go(self):
        self.save_data(0)
        for ts in tqdm(
            np.arange(self.run.macro_params.total_time_steps), mininterval=2
        ):
            current_time = ts * self.run.macro_params.time_step
            if self.run.macro_params.duplicate_fortran:
                current_time += self.run.macro_params.time_step

            if ts == 23:
                pass

            if self.run.macro_params.duplicate_fortran:
                self.random_numbers = np.empty(
                    (8, self.run.macro_params.total_molecules), np.float_
                )
                for i in range(8):
                    self.random_numbers[i] = self.rng.random(
                        self.run.macro_params.total_molecules
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
            if self.run.macro_params.duplicate_fortran:
                move_chance = self.random_numbers[RandomDraw.MOVE]
                should_move = (
                    move_chance > 1 - self.run.macro_params.moving_probability
                ) & ~self.bound
            else:
                move_chance = self.rng.random(self.run.macro_params.total_molecules)
                should_move = (
                    move_chance < self.run.macro_params.moving_probability
                ) & ~self.bound

            # TODO(bpaynter): Since these are all boolean vectors, I should be
            #                 able to use & and | instead of the lookup masks
            conflict = should_bind & should_move
            threshold = (
                current_time - self.binding_time[conflict]
            ) / self.run.macro_params.time_step
            if self.run.macro_params.duplicate_fortran:
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
                >= self.run.macro_params.save_interval * self.current_save_interval
            ):
                self.save_data(current_time)

            if ts % 100000 == 100000 - 1:
                unlysed_fibers = np.count_nonzero(self.fiber_status > current_time)

                if unlysed_fibers == 0:
                    self.logger.info(
                        f"All fibers degraded after {current_time:.2f} sec. Terminating"
                    )
                    self.last_degrade_time = np.max(self.fiber_status)
                    break
                else:
                    unlysed_fiber_percent = (
                        100 - unlysed_fibers / self.run.macro_params.total_fibers * 100
                    )
                    reached_back_row_percent = (
                        self.number_reached_back_row
                        / self.run.macro_params.total_molecules
                        * 100
                    )
                    self.logger.info(
                        f"After {current_time:.2f} sec, "
                        f"{self.run.macro_params.total_fibers - unlysed_fibers:,} "
                        f"fibers are degraded ({unlysed_fiber_percent:.1f}% of total) "
                        f"and {self.number_reached_back_row:,} molecules have reached "
                        f"the back row {reached_back_row_percent:.1f}% of total)."
                    )

        self.save_data(self.run.macro_params.total_time)
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
