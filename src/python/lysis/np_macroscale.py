# *_* coding: utf-8 *_*
#
# Clot Lysis Simulation
# Copyright (C) 2023  Bradley Paynter & Brittany Bannish
#
# np_macroscale.py
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
"""Macroscale clot lysis simulation"""

import logging
from functools import partial

import numpy as np
from tqdm.auto import tqdm

from .util import Experiment, KissRandomGenerator, RandomDraw
from .edge_grid import EdgeGrid, from_fortran_edge_index, to_fortran_edge_index


class MacroscaleSimulation:
    def __init__(self, exp: Experiment, instance: int = None, seed: int = None):
        self.exp = exp
        assert self.exp.params is not None

        self.logger = logging.getLogger(__name__)
        self.logger.debug(f"Initializing MacroscaleRun")
        if seed is None:
            seed = exp.params.seed
        if exp.params.duplicate_fortran:
            self.rng = KissRandomGenerator(seed)
        else:
            self.rng = np.random.default_rng(seed=abs(seed))

        self.binding_time_factory = self._BindingTimeFactory(self.exp, self.rng)

        self.edge_lookup = partial(
            np.ravel_multi_index,
            dims=(exp.params.rows, exp.params.full_row),
        )

        self.fiber_status = np.full(
            self.exp.params.rows * self.exp.params.full_row,
            float("inf"),
            dtype=np.float_,
        )
        self.real_fiber = np.full(
            self.exp.params.rows * self.exp.params.full_row,
            True,
            dtype=np.bool_,
        )
        for i, j in np.ndindex(self.exp.params.empty_rows, self.exp.params.full_row):
            self.real_fiber[self.edge_lookup((i, j))] = False
        for j in range(exp.params.cols):
            self.real_fiber[
                self.edge_lookup(
                    (
                        self.exp.params.rows - 1,
                        3 * j,
                    )
                )
            ] = False
        self.fiber_status[~self.real_fiber] = 0

        self.logger.debug(f"Precalculating neighbors.")
        self.neighbors = EdgeGrid.generate_neighborhood_structure(exp)

        if exp.params.duplicate_fortran:
            perm = np.array(
                [
                    from_fortran_edge_index(k, exp.params.rows, exp.params.cols)
                    for k in range(exp.params.total_edges)
                ]
            )
            perm = np.ravel_multi_index(
                (perm[:, 0], perm[:, 1]),
                dims=(exp.params.rows, exp.params.full_row),
            )
            inv_perm = []
            for k in range(exp.params.rows * exp.params.full_row):
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
        if exp.params.duplicate_fortran:
            location = self.rng.random(exp.params.total_molecules)
            location = exp.params.empty_rows * exp.params.full_row * location
            location = location.astype(int, copy=False)
            location_i = np.empty(exp.params.total_molecules, dtype=np.int_)
            location_j = np.empty(exp.params.total_molecules, dtype=np.int_)
            for m in range(len(location)):
                location_i[m], location_j[m] = from_fortran_edge_index(
                    location[m], exp.params.rows, exp.params.cols
                )
        else:
            location_i = self.rng.integers(
                exp.params.empty_rows,
                size=exp.params.total_molecules,
                dtype=np.short,
            )
            location_j = self.rng.integers(
                exp.params.full_row,
                size=exp.params.total_molecules,
                dtype=np.short,
            )

        self.location = self.edge_lookup((location_i, location_j))

        m_location_ij = np.unravel_index(
            self.location, (self.exp.params.rows, self.exp.params.full_row)
        )

        f_location = np.array(
            [
                to_fortran_edge_index(
                    m_location_ij[0][k],
                    m_location_ij[1][k],
                    self.exp.params.rows,
                    self.exp.params.cols,
                )
                for k in range(self.exp.params.total_molecules)
            ]
        )

        self.m_fiber_status = None

        self.bound = np.full(exp.params.total_molecules, False, dtype=np.bool_)
        self.waiting_time = np.full(exp.params.total_molecules, 0, dtype=np.float_)
        self.binding_time = np.full(
            exp.params.total_molecules, float("inf"), dtype=np.float_
        )
        self.unbound_by_degradation = np.full(
            exp.params.total_molecules, 0, dtype=np.bool_
        )
        self.time_to_reach_back_row = np.full(
            exp.params.total_molecules, float("inf"), dtype=np.float_
        )
        self.reached_back_row = np.full(
            exp.params.total_molecules, False, dtype=np.bool_
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
                self.exp.params.number_of_saves,
                self.exp.params.rows * self.exp.params.full_row,
            ),
            dtype=np.float_,
        )
        self.exp.data.molecule_location = np.empty(
            (
                self.exp.params.number_of_saves,
                self.exp.params.total_molecules,
            ),
            dtype=np.int_,
        )
        self.exp.data.molecule_state = np.empty(
            (
                self.exp.params.number_of_saves,
                self.exp.params.total_molecules,
            ),
            dtype=np.bool_,
        )
        self.exp.data.save_time = np.empty(
            (self.exp.params.number_of_saves,),
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
            exp: The Experiment object that this model is a part of
            rng: A generator of random numbers. This should be the same
                    generator object that the rest of the model is using to
                    prevent overlap issues.
        """

        def __init__(self, exp: Experiment, rng: np.random.Generator):
            # The size of the internal list
            self.period = max(1000000, 10 * exp.params.total_molecules)
            # The experiment this object is a part of
            self.exp = exp
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
            denominator = self.exp.params.binding_rate * self.exp.params.binding_sites
            # Use the vectorized, in-place divide function for the middle term
            binding_time_list = np.divide(
                binding_time_list, denominator, out=binding_time_list
            )
            # Using vectorized, in-place subtract to combine the last two terms
            # Using subtract so that we don't have to multiply
            # each 'ln(X)' term by -1. Don't know if this makes a performance
            # difference or not.
            binding_time_list = np.subtract(
                -self.exp.params.time_step / 2,
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
            m: A boolean array of molecules. Used as a mask when applying this
                function.
            current_time: The current simulation time.
        """
        # Count how many molecules need to be unbound.
        # That is, how many elements of `m` are `True`
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
            + self.exp.params.average_bind_time
            - self.exp.params.time_step / 2
        )
        # Set the binding/unbinding time of the selected molecules to infinite
        # The fiber in their current location just degraded, so they can never
        # bind to it.
        self.binding_time[m] = float("inf")

    def unbind_by_time(self, m: np.ndarray, current_time: float):
        """Deals with bound molecules whose 'unbind' timer just expired.

        Args:
            m: A boolean array of molecules. Used as a mask when applying this
                function.
            current_time: The current simulation time.
        """
        # Count how many molecules need to be unbound.
        # That is, how many elements of `m` are `True`
        count = np.count_nonzero(m)
        # If none, then return without any more work
        if count == 0:
            return
        # Set the bound status
        # True if we are currently bound and not selected for unbinding
        self.bound = self.bound & ~m
        # Set the 'macro unbind' flag to False
        self.unbound_by_degradation = self.unbound_by_degradation & ~m

        # Deal with Micro unbinding
        # Allocate space
        forced = np.full(self.exp.params.total_molecules, False, dtype=np.bool_)
        # Draw random numbers
        if self.exp.params.duplicate_fortran:
            random = self.random_numbers[RandomDraw.MICRO_UNBIND][m]
        else:
            random = self.rng.random(count)
        # Determine which molecules were forced to unbind
        forced[m] = random <= self.exp.params.forced_unbind
        # Calculate waiting time
        self.waiting_time[forced] = (
            current_time
            + self.exp.params.average_bind_time
            - self.exp.params.time_step / 2
        )
        # Micro unbound molecules cannot bind
        self.binding_time[forced] = float("inf")
        num_forced = np.count_nonzero(forced)
        self.total_micro_unbinds += num_forced

        # If a molecule unbound from an intact fiber and is not micro-unbound
        # then it could rebind at some point, so we should calculate a binding
        # time.
        if num_forced < count:
            # Determine the amount of simulation time from now that the
            # molecule will bind
            if self.exp.params.duplicate_fortran:
                binding_times = self.binding_time_factory.fill_list(
                    self.random_numbers[RandomDraw.BINDING_TIME_WHEN_UNBINDING][
                        m & ~forced
                    ]
                )
            else:
                binding_times = self.binding_time_factory.next(count - num_forced)
            # Add the binding time to the current time and store.
            self.binding_time[m & ~forced] = current_time + binding_times

    def find_unbinding_time(
        self, unbinding_time_bin: np.ndarray, current_time: float
    ) -> np.ndarray:
        """Find the simulation time at which the tPA molecule will unbind

        If we think of the `binding_time` array as a function :math:`f(x)` where
        `binding_time[100*x] =` :math:`f(x)`, we want to draw uniformly from
        the range of :math:`f`. That is, we want the unbinding timestep to be
        :math:`f(r)` where :math:`r\\sim U(0,1)`.

        If :math:`100r` is an integer :math:`i`, this is easy because we
        just look up `binding_time[i]`, but what if it isn't a perfect integer?

        Then :math:`f(r)` lies in the interval
        :math:`\\Bigl(f\\bigl(\\lfloor 100r \\rfloor\\bigr),
        f\\bigl(\\lceil 100r \\rceil\\bigr)\\Bigr)` so we
        interpolate it linearly. That is, define a linear function :math:`g(x)`
        such that
        :math:`g\\bigl(\\lfloor 100r \\rfloor\\bigr) =
        f\\bigl(\\lfloor 100r \\rfloor\\bigr)` and
        :math:`g\\bigl(\\lceil 100r \\rceil\\bigr) =
        f\\bigl(\\lceil 100r \\rceil\\bigr)`.
        Then use :math:`g(r)` to approximate :math:`f(r)`.

        Another way to see this is, if :math:`i` is an integer such that
        :math:`100i < 100r < 100(i+1)` and
        :math:`\\lambda` is such that :math:`100i + \\lambda = 100r`,
        then
        :math:`f(r) \\approx (1-\\lambda)f(\\frac{i}{100}) +
        \\lambda f(\\frac{i+1}{100})`.

        Args:
            unbinding_time_bin: Which bin of unbinding times are we drawing from?
                Entries should be floats, drawn uniformly from (0, 100).
            current_time: The current simulation time.

        Returns: An array of unbinding times, randomly drawn from the provided
            microscale data.
        """
        # Use NumPy's interp function to do the interpolation of the microscale
        # unbinding data
        interp = np.interp(unbinding_time_bin, self.xp, self.exp.data.unbinding_time)
        # Add the current simulation time and a half timestep. Then return
        return interp + (current_time - self.exp.params.time_step / 2)

    def find_lysis_time(
        self,
        m: np.ndarray,
        unbinding_time_bin: np.ndarray,
        current_time: float,
        count: int,
    ) -> np.ndarray:
        """Find the (possibly infinite) simulation time at which the molecule
        will degrade the fiber to which it just bound.

        Args:
            m: A boolean array of molecules. Used as a mask when applying this
                function.
            unbinding_time_bin: Which bin of unbinding times are we drawing from?
                This should be the same array passed to the
                `find_unbinding_time` method.
                Entries should be floats, drawn uniformly from (0, 100).
            current_time: The current simulation time.
            count: The number of molecules for which we are finding lysis times.

        Returns: An array of lysis times, randomly drawn from the provided
            microscale data.

        """
        # Draw random numbers to determine where we will interpolate the lysis
        # time from.
        # `unbinding_time_bin` determines the row of the `lysis_time` matrix,
        # while `lysis_time_bin` determines the column.
        if self.exp.params.duplicate_fortran:
            lysis_time_bin = self.random_numbers[RandomDraw.LYSIS_TIME][m]
        else:
            lysis_time_bin = self.rng.random(count)
        lysis_time_bin = lysis_time_bin * (self.exp.params.microscale_runs / 100)
        # Allocate memory for the lysis times
        interp = np.full(count, float("inf"), dtype=np.double)
        # Floor the `unbinding_time_bin` so that we can use it indexes the
        # appropriate row of the lysis matrix
        unbinding_time_bin = unbinding_time_bin.astype(int)
        # Determine how many entries in the current row of `lysis_time` have
        # finite entries (infinite entries are stored as `6,000` for legacy
        # reasons).
        # The `total_lyses` array stores the address of the first infinite entry
        # so we need to subtract one to find the number of finite entries.
        total_lyses = self.exp.data.total_lyses[unbinding_time_bin] - 1
        # Determine which bindings actually result in lysis (non-infinite)
        lysis_happens = lysis_time_bin < total_lyses
        # Set the other entries to infinite lysis times.
        interp[~lysis_happens] = float("inf")

        # TODO(bpaynter): Should be able to improve this with a 2D
        #                 interpolation or a custom numpy kernel
        # For each molecule where lysis takes place
        for i in np.arange(count)[lysis_happens]:
            # Interpolate the `lysis_time_bin` into the appropriate row of the
            # `lysis_time` matrix using the NumPy `interp` function.
            interp[i] = np.interp(
                lysis_time_bin[i],
                np.arange(total_lyses[i]),
                self.exp.data.lysis_time[unbinding_time_bin[i], : total_lyses[i]],
            )
        # Add the current time and a half timestep. Then return.
        return interp + (current_time - self.exp.params.time_step / 2)

    def bind(self, m: np.ndarray, current_time: float):
        """Deals with molecules that just bound to a fiber.

        Args:
            m: A boolean array of molecules. Used as a mask when applying this
                function.
            current_time: The current simulation time.
        """
        # Count how many molecules need to be bound.
        # That is, how many elements of `m` are `True`
        count = np.count_nonzero(m)
        # If none, then return without any more work
        if count == 0:
            return

        # Set the molecule status to 'bound'
        # This is true if the molecule was already bound, or if it is binding
        # on this timestep
        self.bound = self.bound | m
        # Reset the molecule's waiting time
        self.waiting_time[m] = 0
        # Add to the total number of binds for the simulation
        self.total_binds += count

        # Draw random numbers for the unbinding time
        if self.exp.params.duplicate_fortran:
            unbinding_time_bin = self.random_numbers[RandomDraw.UNBINDING_TIME][m] * 100
        else:
            unbinding_time_bin = self.rng.random(count) * 100
        # Determine the time at which the molecule would unbind itself
        self.binding_time[m] = self.find_unbinding_time(
            unbinding_time_bin, current_time
        )
        # Determine the (possibly infinite) simulation time at which the
        # molecule would degrade this fiber
        lysis_time = self.find_lysis_time(m, unbinding_time_bin, current_time, count)
        # Find the fibers to which the molecules have just been bound
        locations = self.location[m]
        # For each molecule
        for i in range(count):
            # Check if the molecule would degrade the fiber
            if lysis_time[i] < float("inf"):
                # Set the fiber's status (degrade time) if the new time is
                # smaller than its previously stored time.
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
        """Deal with molecules that need to move this timestep, but are still
        bound to a fiber particle due to macro- or micro-unbinding.

        Also known as 'restricted move'.

        Here, a molecule cannot move to an edge containing an intact fiber.
        The new location for the molecule is drawn uniformly from a set
        containing the molecule's current location and any neighboring edge
        without a fiber (whether always empty or empty due to degradation).

        Args:
            m: A boolean array of molecules. Used as a mask when applying this
                function.
            current_time: The current simulation time.
        """
        # Count how many molecules need to be moved.
        # That is, how many elements of `m` are `True`
        count = np.count_nonzero(m)
        # If none, then return without any more work
        if count == 0:
            return

        # Add to the total restricted moves for the simulation
        self.total_restricted_moves += count
        # Get the current location of the molecules
        current_locations = self.location[m]

        # Get the full neighborhoods of the molecule's current location
        neighborhoods = self.neighbors[current_locations]
        # Neighbors that are valid for movement were either empty from the
        # start of the simulation (`fiber_status` = 0) or have degraded by this
        # point of the simulation (`fiber_status` < current_time)
        valid_neighbors = self.fiber_status[neighborhoods] < current_time
        # We then determine the indices of the valid neighbors
        # (True sorts before False)
        valid_neighborhood_index = np.argsort(~valid_neighbors, axis=1)
        # And sort so that the valid neighbors are at the front of each
        # molecule's neighborhood
        valid_neighborhoods = np.take_along_axis(
            neighborhoods, valid_neighborhood_index, axis=1
        )
        # Add the molecules' current locations to the front of their array of
        # valid neighbors
        valid_neighborhoods = np.append(
            current_locations.reshape(count, 1), valid_neighborhoods, axis=1
        )
        # Determine the number of valid neighbors for each molecule
        # (not counting its current location)
        num_valid_neighbors = np.count_nonzero(valid_neighbors, axis=1)
        # Draw random numbers to determine which neighbor
        if self.exp.params.duplicate_fortran:
            neighbor = self.random_numbers[RandomDraw.RESTRICTED_MOVE][m]
        else:
            neighbor = self.rng.random(count)
        # Convert the random float into an integer
        neighbor = (neighbor * (num_valid_neighbors + 1)).astype(int, copy=False)
        # Set the molecules' locations to the randomly determined neighbor.
        self.location[m] = valid_neighborhoods[np.full(count, True), neighbor]

    def unrestricted_move(self, m: np.ndarray, current_time: float):
        """Deals with molecules that need to move this timestep and are free to
        move to any neighboring edge.

        Args:
            m: A boolean array of molecules. Used as a mask when applying this
                function.
            current_time: The current simulation time.
        """
        # Count how many molecules need to be moved.
        # That is, how many elements of `free_to_move` are `True`
        count = np.count_nonzero(m)
        # If none, then return without any more work
        if count == 0:
            return

        # Draw random numbers to determine the neighbor we should move to
        if self.exp.params.duplicate_fortran:
            # The Fortran code uses the same random number that determines IF a
            # molecule should move to determine WHERE a molecule should move
            neighbor = self.random_numbers[RandomDraw.MOVE][m]
            # So all of these numbers will be in the interval (.8, 1), so they
            # first need to be scaled back to (0, 1)
            neighbor = neighbor - (1 - self.exp.params.moving_probability)
            neighbor = neighbor / self.exp.params.moving_probability
            # Then to (0, 8)
            neighbor = neighbor * 8
            # Then converted to an integer in {1..8}
            neighbor = neighbor.astype(int, copy=False)
        else:
            # All neighborhoods have eight entries, so choose one uniformly
            neighbor = self.rng.integers(8, size=count)

        # Set the molecules' new location to that of the selected neighbor.
        self.location[m] = self.neighbors[self.location[m], neighbor]
        # Add to the total moves for this simulation
        self.total_regular_moves += count

        # TODO(bpaynter): It might be quicker to only assign binding times
        #   to those molecules whose new location contains a non-degraded fiber

        # Generate binding times for each molecule
        if self.exp.params.duplicate_fortran:
            binding_time = self.binding_time_factory.fill_list(
                self.random_numbers[RandomDraw.BINDING_TIME_WHEN_MOVING][m]
            )
        else:
            binding_time = self.binding_time_factory.next(count)
        # Add the binding times to the current time and store
        self.binding_time[m] = current_time + binding_time

    def move(self, m: np.ndarray, current_time: float):
        """Deals with all molecules that need to move on this timestep

        Args:
            m: A boolean array of molecules. Used as a mask when applying this
                function.
            current_time: The current simulation time.
        """
        # Determine which molecules are restricted in their movement. These are
        # molecules that have been macro- or micro-unbound.
        still_stuck_to_fiber = (
            (self.waiting_time > current_time) & self.unbound_by_degradation & m
        )
        # Process those molecules that have restricted movement
        self.move_to_empty_edge(still_stuck_to_fiber, current_time)

        # Determine the molecules which are free to move to any neighbor
        free_to_move = m & ~still_stuck_to_fiber
        # Process those molecules that are unrestricted in their movement.
        self.unrestricted_move(free_to_move, current_time)

        # If there are still molecules which have not yet passed through the clot
        if self.number_reached_back_row < self.exp.params.total_molecules:
            back_row = (self.exp.params.rows - 1) * self.exp.params.full_row - 1
            # Determine those molecules which are currently on the back row of
            # the grid and have not been counted before
            first_time = (self.location > back_row) & ~self.reached_back_row
            # Count how many molecules have passed through for the first time
            # on this timestep
            count = np.count_nonzero(first_time)
            if count > 0:
                # Store the 'first passage time' for these molecules
                self.time_to_reach_back_row[first_time] = current_time
                # Flag these molecules as having passed through the clot
                self.reached_back_row = self.reached_back_row | first_time
                # Add the number to the total for the simulation
                self.number_reached_back_row += count

    def save_data(self, current_time):
        """Saves the current simulation state.

        Args:
            current_time: The current simulation time.
        """
        self.exp.data.degradation_state[self.current_save_interval] = self.fiber_status
        self.exp.data.molecule_location[self.current_save_interval] = self.location
        self.exp.data.molecule_state[self.current_save_interval] = self.bound
        self.exp.data.save_time[self.current_save_interval] = current_time
        self.current_save_interval += 1

    def record_data_to_disk(self):
        """Records the saved data to disk."""
        self.logger.info(f"Saving data to disk.")
        self.exp.data.save_to_disk("degradation_state")
        self.exp.data.save_to_disk("molecule_location")
        self.exp.data.save_to_disk("molecule_state")
        self.exp.data.save_to_disk("save_time")

    def run(self):
        """Runs the simulation."""
        # Save the initial state of the simulation.
        self.save_data(0)
        # Start the main loop of the simulation, looping over all timesteps.
        for ts in tqdm(np.arange(self.exp.params.total_time_steps), mininterval=2):
            # We process each timestep at the end of that time interval.
            # Therefor, the current time is the number of timesteps,
            # times the length of one timestep,
            # plus the length of the current timestep.
            current_time = (ts + 1) * self.exp.params.time_step

            # If we are matching the 'rng-array' Fortran code step-by-step,
            # then we need to roll all random numbers that might be used for
            # this timestep at the beginning.
            if self.exp.params.duplicate_fortran:
                # Allocate space for eight random numbers per molecule.
                # See lysis.util.constants.RandomDraw for their use-codes
                self.random_numbers = np.empty(
                    (8, self.exp.params.total_molecules), np.float_
                )
                # Draw random numbers
                for i in range(8):
                    self.random_numbers[i] = self.rng.random(
                        self.exp.params.total_molecules
                    )

            # Store the fiber status of each molecule's current location.
            # This sort of complex indexing function is resource intensive and
            # reducing the number of calls of this type will improve running
            # time.
            self.m_fiber_status = self.fiber_status[self.location]

            # Find any molecules bound to a degraded fiber and deal with them.
            self.unbind_by_degradation(
                self.bound & (self.m_fiber_status < current_time), current_time
            )
            # Find any bound molecules whose binding/unbinding timer has
            # expired and deal with them
            self.unbind_by_time(
                self.bound & (self.binding_time < current_time), current_time
            )

            # Determine which molecules are ready to bind on this timestep.
            should_bind = (
                ~self.bound  # Not currently bound
                & (self.binding_time < current_time)  # Binding time passed
                & (self.m_fiber_status > current_time)  # Fiber not degraded
                & (self.waiting_time < current_time)  # Not waiting
            )

            # Draw random numbers to determine if molecules should move on this
            # timestep
            if self.exp.params.duplicate_fortran:
                move_chance = self.random_numbers[RandomDraw.MOVE]
            else:
                move_chance = self.rng.random(self.exp.params.total_molecules)
            # Molecules should move if their moving chance is rolled and they
            # are not currently bound
            # NOTE: It is faster to roll for all molecules and mask out those
            # that are bound, than it is to look up those that are unbound and
            # just roll for them
            should_move = (
                move_chance > 1 - self.exp.params.moving_probability
            ) & ~self.bound

            # TODO(bpaynter): Since these are all boolean vectors, I should be
            #                 able to use & and | instead of the lookup masks
            # Find those molecules that are scheduled to both bind and move on
            # this timestep
            conflict = should_bind & should_move
            # Determine how far into this timestep the molecule will bind
            threshold = (
                current_time - self.binding_time[conflict]
            ) / self.exp.params.time_step

            # Draw random numbers to resolve move & bind conflicts
            if self.exp.params.duplicate_fortran:
                random = self.random_numbers[RandomDraw.CONFLICT_RESOLUTION][conflict]
            else:
                random = self.rng.random(np.count_nonzero(conflict))
            # The molecule binds if its move time (uniformly drawn) is earlier(??)
            # than its binding time.
            should_bind[conflict] = random <= threshold
            # If it has been chosen to bind, it should not move
            should_move[conflict] = ~should_bind[conflict]

            # Deal with all molecules that will bind on this timestep
            self.bind(should_bind, current_time)
            # Deal with all molecules that should move on this timestep
            self.move(should_move, current_time)

            # If a full save interval has passed, then save the current
            # simulation state.
            if (
                current_time
                >= self.exp.params.save_interval * self.current_save_interval
            ):
                self.save_data(current_time)

            # If 100,000 timesteps have passed,
            # output status and check termination criteria
            if ts % 100_000 == 100_000 - 1:
                # Count the number of non-degraded fibers
                unlysed_fibers = np.count_nonzero(self.fiber_status > current_time)
                # Determine the percentage of fibers that have degraded
                lysed_fiber_percent = (
                    100 - unlysed_fibers / self.exp.params.total_fibers * 100
                )
                # Determine the percentage of molecules that have passed
                # through the clot
                reached_back_row_percent = (
                    self.number_reached_back_row / self.exp.params.total_molecules * 100
                )
                # Output status to log
                self.logger.info(
                    f"After {current_time:.2f} sec, "
                    f"{self.exp.params.total_fibers - unlysed_fibers:,} "
                    f"fibers are degraded ({lysed_fiber_percent:.1f}% of total) "
                    f"and {self.number_reached_back_row:,} molecules have reached "
                    f"the back row {reached_back_row_percent:.1f}% of total)."
                )
                # Determine if the simulation should terminate
                if unlysed_fibers == 0 and reached_back_row_percent > 95:
                    self.last_degrade_time = np.max(self.fiber_status)
                    break

        # Save the final status of the simulation
        self.save_data(self.exp.params.total_time)
        # Store simulation data to disk
        self.record_data_to_disk()

        # Output final simulation metrics to log
        self.logger.info(f"Total binds: {self.total_binds:,}")
        self.logger.info(f"Total regular moves: {self.total_regular_moves:,}")
        self.logger.info(f"Total restricted moves: {self.total_restricted_moves:,}")
        self.logger.info(f"Total macro unbinds: {self.total_macro_unbinds:,}")
        self.logger.info(f"Total micro unbinds: {self.total_micro_unbinds:,}")
        self.logger.info(
            f"Molecules which reached the back row: {self.number_reached_back_row:,}"
        )
        self.logger.info(f"Last fiber degraded at: {self.last_degrade_time:2f} sec")
