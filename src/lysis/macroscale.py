"""NumPy-based macroscale simulation of tPA-mediated fibrin clot lysis.

This module implements the macroscale component of the multi-scale clot lysis
simulation. It models the diffusion of tissue plasminogen activator (tPA) molecules
through a fibrin network and their binding to and degradation of individual fibers.

The macroscale simulation operates on a rectilinear edge grid where each edge
represents a fibrin fiber. tPA molecules perform a random walk through the grid,
binding to fibers, causing degradation, and unbinding according to statistics
derived from microscale simulations.

Simulation Algorithm
--------------------

Each timestep involves unbinding, movement, binding, conflict resolution, fiber
degradation, and data collection. See :meth:`MacroscaleSim.go` for the detailed
step-by-step algorithm.

Operational Modes
-----------------

The simulation supports two modes:

**Native NumPy mode** (default):
    - Optimized for performance using NumPy's Mersenne Twister RNG
    - Generates random numbers on-demand as needed
    - Vectorized operations for maximum speed
    - Recommended for production simulations

**Fortran duplication mode** (``duplicate_fortran=True``):
    - Replicates the original Fortran implementation exactly
    - Uses KISS RNG for identical random number sequences
    - Pre-generates all random numbers per timestep in Fortran order
    - Used for validation and comparison with legacy code

Key Features
------------

- **Vectorization**: Extensive use of NumPy array operations for performance
- **Pre-computed neighbors**: EdgeGrid neighbor structure calculated once at startup
- **Batch binding time generation**: _BindingTimeFactory amortizes RNG costs
- **Microscale integration**: Uses binding/unbinding/lysis statistics from microscale
- **Restricted movement**: Molecules unbound by degradation temporarily restricted
- **Progress tracking**: Monitors degradation progress and molecule migration

Main Classes
------------

**MacroscaleSim**:
    The primary simulation class. Initializes the grid, places molecules, and
    executes the main simulation loop. Contains all state arrays (fiber_status,
    molecule locations, binding states) and simulation logic.

**_BindingTimeFactory** (nested):
    Factory class for efficient generation of binding times. Pre-generates large
    batches of binding times and serves them on demand, dramatically improving
    performance over generating times individually.

Molecule States
---------------

Molecules can be in one of several states:

- **Unbound**: Free to move via random walk, may bind to intact fibers
- **Bound**: Attached to a fiber, cannot move, may cause lysis
- **Macro-unbound**: Forcibly unbound by fiber degradation, restricted movement
- **Micro-unbound**: Forcibly unbound during natural unbinding, restricted movement

Restricted movement means the molecule can only move to neighboring edges that
have already degraded, or remain in place.

Data Storage
------------

Periodic snapshots of fiber and molecule state are saved at regular intervals.
See :meth:`MacroscaleSim.save_data` and :meth:`MacroscaleSim.record_data_to_disk`
for details on which arrays are stored.

Performance Considerations
--------------------------

The simulation is optimized for large-scale problems:

- Neighbor structure pre-computation: O(1) lookups vs O(1) calculations per timestep
- Vectorized operations: Process all molecules simultaneously when possible
- Batch RNG: Generate binding times in large batches (1M+ at a time)
- In-place operations: Minimize memory allocations during time-stepping
- Linear scaling with molecule count and timestep count

Example Usage
-------------

Basic simulation setup::

    >>> from lysis.config.run import Run
    >>> from lysis.macroscale import MacroscaleSim
    >>>
    >>> # Create run with parameters
    >>> run = Run(...)
    >>> run.macro_params.rows = 100
    >>> run.macro_params.cols = 100
    >>> run.macro_params.total_molecules = 10000
    >>> run.macro_params.total_time_steps = 1000000
    >>>
    >>> # Initialize and run simulation
    >>> sim = MacroscaleSim(run)
    >>> sim.go()
    >>>
    >>> # Results are stored in run.data
    >>> print(f"Final degradation: {sim.fiber_status}")

Fortran compatibility mode for validation::

    >>> run.macro_params.duplicate_fortran = True
    >>> run.macro_params.macro_seed = 12345
    >>> sim = MacroscaleSim(run)
    >>> sim.go()  # Will exactly match Fortran output

Access saved data::

    >>> # Data is saved to disk after simulation
    >>> sim_view = run.data.macroscale_out[0]
    >>> times = sim_view.snapshot_time[:]
    >>> locations = sim_view.tpa_location_snapshot[:]
    >>> events = sim_view.tpa_bind_events[:]

Algorithm Details
-----------------

**Binding time calculation**:
    Binding times are drawn from an exponential distribution:
    t_bind = -ln(X) / (rate * sites) - dt/2
    where X ~ Uniform(0,1), rate is the binding rate, and sites is the number
    of binding sites per fiber.

**Unbinding time calculation**:
    Unbinding times are drawn from microscale simulation data. A random uniform
    value selects a bin [0, 100], and linear interpolation provides the
    unbinding time within that bin.

**Lysis time calculation**:
    Lysis times (when bound molecules cause fiber degradation) are also drawn
    from microscale data, indexed by the unbinding time bin. Not all bindings
    result in lysis.

**Conflict resolution**:
    When both binding and movement are scheduled in the same timestep, the
    probability of binding increases linearly with how overdue the binding is:
    P(bind) = (t_current - t_bind) / dt

Limitations and Assumptions
----------------------------

- Assumes well-mixed system at the microscale (binding/unbinding rates uniform)
- Molecules do not interact with each other (independent random walks)
- Fibers degrade instantly and completely (no partial degradation)
- Boundary conditions are reflecting or periodic (no molecule creation/destruction)
- Grid is rectilinear (not truly hexagonal, though edge structure mimics it)

See Also
--------

- ``geometry.edge_grid``: Rectilinear grid structure and neighbor calculations
- ``config.parameters``: MacroParameters configuration class
- ``config.run``: Run container for parameters and data

References
----------

.. TODO:: Add published references for the multi-scale fibrinolysis model.

Notes
------

The simulation can terminate early if all fibers degrade before the configured
end time. Check ``sim.last_degrade_time`` for the actual termination time.

Progress is logged every 100K timesteps, showing degradation percentage and
molecules reaching the back row.
"""

import itertools
import logging
import os
from functools import partial

import numpy as np
from tqdm.auto import tqdm

from .config.constants import MolStatus, RandomDraw
from .config.run import Run
from .dataio.dataspec import dataspec
from .geometry.edge_grid import EdgeGrid, from_fortran_edge_index, to_fortran_edge_index
from .tools.kiss import KissRandomGenerator

_macro_out_spec = dataspec["v2.0.0"]["macroscale_out"]

#: NumPy structured dtype for tPA bind/unbind event logging (from dataspec).
_BIND_EVENT_DTYPE = _macro_out_spec.data["tpa_bind_events"].dtype

#: NumPy structured dtype for fiber degrade-time event logging (from dataspec).
_FIBER_DEGRADE_DTYPE = _macro_out_spec.data["fiber_degrade_time"].dtype


class MacroscaleSim:
    """Macroscale simulation of tPA-mediated fibrin clot lysis.

    This class implements a stochastic simulation of tissue plasminogen activator
    (tPA) molecules moving through and degrading a fibrin network represented as
    a hexagonal edge grid. The simulation tracks molecule binding/unbinding,
    movement via random walk, and fiber degradation over time.

    The simulation supports two modes:
    - Native NumPy mode: Optimized for performance using NumPy's RNG
    - Fortran duplication mode: Replicates the original Fortran implementation
    step-by-step for validation purposes

    Key features:
    - Pre-computed neighborhood structure for efficient lookups
    - Batch generation of binding times for performance
    - Integration with microscale simulation data for binding/lysis statistics
    - Periodic data snapshots for analysis
    """

    def __init__(self, run: Run, sim_number: int = 0):
        """Initialize a macroscale simulation.

        Sets up the simulation grid, places molecules, initializes fiber states,
        pre-calculates neighbor relationships, and allocates data storage arrays.

        :param run: The Run object containing simulation parameters and data storage.
            Must have macro_params defined.
        :type run: Run
        :param sim_number: The simulation index within the Run (0-based). Used to
            select the correct per-simulation datasets in macroscale_out.
        :type sim_number: int
        :raises AssertionError: If run.macro_params is None
        """
        self.run = run
        self.sim_number = sim_number
        # Check that the run has macro parameters
        assert self.run.macro_params is not None

        # Initialize logger
        self.logger = logging.getLogger(__name__)
        self.logger.debug(f"Initializing MacroscaleSim")

        # Initialize pseudorandom number generator
        if run.macro_params.duplicate_fortran:
            # If we are copying the Fortran code step-by-step,
            # we want to use the same RNG which is the KISS C code
            self.rng = KissRandomGenerator(run.macro_params.macro_seed)
        else:
            # Use NumPy's default RNG (Mersenne Twister) for native mode
            self.rng = np.random.default_rng(seed=abs(run.macro_params.macro_seed))

        # Initialize the Binding Time Factory which will generate random binding times for the simulation
        # We use the factory to generate these times in batches, giving us faster amortized run times
        self.binding_time_factory = self._BindingTimeFactory(self.run, self.rng)

        # We humans would rather interact with the edge grid in 2-D co-ordinates,
        # but it is faster to work with a 1-D array,
        # so we use the NumPy function ravel_multi_index to do that for us.
        # We use functools.partial to pre-fill the grid dimensions
        self.edge_lookup = partial(
            np.ravel_multi_index,
            dims=(run.macro_params.rows, run.macro_params.full_row),
        )

        # Initialize the status of all fibers to intact and never degrading.
        # This is achieved by setting all degrade times to infinity
        self.fiber_status = np.full(
            self.run.macro_params.rows * self.run.macro_params.full_row,
            float("inf"),
            dtype=np.float64,
        )
        # Determine which edges represent real fibers vs. non-existent edges.
        # Three types of non-existent edges:
        # 1. Empty rows at the front (where tPA starts external to clot)
        # 2. Missing vertical edges on the top row (hexagonal grid boundary)
        #
        # Start by initializing all edges as real fibers
        self.real_fiber = np.full(
            self.run.macro_params.rows * self.run.macro_params.full_row,
            True,
            dtype=np.bool_,
        )
        # Mark empty rows at the front as non-existent
        for i, j in np.ndindex(
            self.run.macro_params.empty_rows, self.run.macro_params.full_row
        ):
            self.real_fiber[self.edge_lookup((i, j))] = False
        # Find the non-existent vertical edges on the top row
        for j in range(run.macro_params.cols):
            self.real_fiber[
                self.edge_lookup(
                    (
                        self.run.macro_params.rows - 1,
                        3 * j,
                    )
                )
            ] = False
        # Set non-existent edges to already degraded (degrade time = 0) so molecules
        # can move through them freely
        self.fiber_status[~self.real_fiber] = 0

        # Precalculate the neighbors of all edges.
        # This is done once at the start of the simulation and generates significant
        # run-time savings.
        self.logger.debug(f"Precalculating neighbors.")
        if run.macro_params.duplicate_fortran:
            self.neighbors = EdgeGrid.generate_fortran_neighborhood_structure(run)
        else:
            self.neighbors = EdgeGrid.generate_neighborhood_structure(run)

        self.logger.info(f"Placing molecules on empty edges.")
        if run.macro_params.duplicate_fortran:
            # Fortran mode: Replicate the original Fortran placement algorithm
            # Generate random positions as flat indices in [0, empty_rows * full_row)
            location = self.rng.random(run.macro_params.total_molecules)
            location = (
                run.macro_params.empty_rows * run.macro_params.full_row * location
            )
            location = location.astype(int, copy=False)
            # Convert flat indices to (i,j) coordinates using Fortran indexing
            location_i = np.empty(run.macro_params.total_molecules, dtype=np.int_)
            location_j = np.empty(run.macro_params.total_molecules, dtype=np.int_)
            for m in range(len(location)):
                location_i[m], location_j[m] = from_fortran_edge_index(
                    location[m], run.macro_params.rows, run.macro_params.cols
                )
        else:
            # Native mode: Generate random (i,j) coordinates directly
            # More efficient and readable than the Fortran approach
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

        # Convert (i,j) coordinates to flat indices for efficient array operations
        self.location = self.edge_lookup((location_i, location_j))

        # Cache for fiber status at molecule locations (populated each timestep)
        self.m_fiber_status = None

        # Molecule binding state: True if currently bound to a fiber
        self.bound = np.full(run.macro_params.total_molecules, False, dtype=np.bool_)
        # Waiting time: molecules must wait until this time before moving freely after macro unbind
        self.waiting_time = np.full(
            run.macro_params.total_molecules, 0, dtype=np.float64
        )
        # Binding time: when the molecule will bind (if unbound) or unbind (if bound)
        self.binding_time = np.full(
            run.macro_params.total_molecules, float("inf"), dtype=np.float64
        )
        # Flag: True if molecule was unbound due to fiber degradation (macro unbind)
        self.unbound_by_degradation = np.full(
            run.macro_params.total_molecules, 0, dtype=np.bool_
        )
        # Track when each molecule first reaches the back row of the clot
        self.time_to_reach_back_row = np.full(
            run.macro_params.total_molecules, float("inf"), dtype=np.float64
        )
        # Flag: True if molecule has reached the back row at least once
        self.reached_back_row = np.full(
            run.macro_params.total_molecules, False, dtype=np.bool_
        )
        # X-axis values [0, 1, 2, ..., 100] for interpolating microscale data
        self.xp = np.arange(101)
        # Storage for pre-generated random numbers (Fortran mode only)
        self.random_numbers = None

        # Statistics: count of forced unbinds due to fiber degradation
        self.total_macro_unbinds = 0
        # Statistics: count of forced unbinds during natural unbinding
        self.total_micro_unbinds = 0
        # Statistics: count of binding events
        self.total_binds = 0
        # Statistics: (currently unused)
        self.independent_binds = 0
        # Statistics: count of unrestricted moves (normal random walk)
        self.total_regular_moves = 0
        # Statistics: count of restricted moves (to degraded edges only)
        self.total_restricted_moves = 0
        # Statistics: count of timesteps where fiber degradation times changed
        self.timesteps_with_fiber_changes = 0
        # Track how many molecules have reached the back row
        self.number_reached_back_row = 0
        # Track when the last fiber degrades (for early termination)
        self.last_degrade_time = float("inf")

        # Index for the next save operation
        self.current_save_interval = 0

        # Pre-allocate arrays to store simulation snapshots at save intervals.
        # The capacity is only an initial estimate: save_data() doubles these
        # buffers on demand, so run-to-completion simulations (total_time == 0,
        # number_of_saves == 1) whose length is unknown still work.  Start at
        # a minimum of 1 so the doubling in _grow_snapshot_buffers() can grow.
        initial_saves = max(1, self.run.macro_params.number_of_saves)
        # Stores molecule locations as (row, rank) at each save point
        self.tpa_location_snapshot = np.empty(
            (
                initial_saves,
                2,
                self.run.macro_params.total_molecules,
            ),
            dtype=np.uint32,
        )
        # Stores the simulation time corresponding to each save point
        self.snapshot_time = np.empty(
            (initial_saves,),
            dtype=np.float64,
        )

        # Event accumulation lists — each element is a structured numpy array
        self._bind_events = []
        self._fiber_degrade_events = []

        self.logger.debug(f"Initialization complete.")

    class _BindingTimeFactory:
        """Factory for efficient batch generation of molecule binding times.

        This class pre-generates binding times in large batches to amortize the
        cost of random number generation. It maintains an internal buffer that is
        automatically refilled when exhausted, providing better performance than
        generating binding times one at a time or in small batches.

        Binding times are calculated as:
            binding_time = -ln(X) / (binding_rate * binding_sites) - time_step/2
        where X ~ Uniform(0,1).

        Note: The current time is NOT included in generated binding times and must
        be added by the caller.

        :param run: The Run object containing simulation parameters
        :type run: Run
        :param rng: Random number generator shared with the main simulation to
            prevent statistical overlap
        :type rng: np.random.Generator or KissRandomGenerator
        """

        def __init__(self, run: Run, rng: np.random.Generator):
            # Buffer size: 1M minimum, or 10x molecule count (to amortize refill costs)
            self.period = max(1_000_000, 10 * run.macro_params.total_molecules)
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
            """Generate binding times from random numbers.

            This method operates in two modes:

            1. Internal mode (random_numbers=None): Refills the internal buffer with
                new binding times. This is the standard mode that provides amortized
                performance benefits. Returns None.

            2. External mode (random_numbers provided): Calculates binding times for
                the provided random numbers without affecting the internal buffer.
                Returns the calculated binding times. Does NOT benefit from amortization.

            Binding times are calculated as:
                binding_time = -ln(X) / (binding_rate * binding_sites) - time_step/2
            where X values come from either the RNG or the provided array.

            Note: Generated binding times do NOT include the current time. The caller
            must add current_time to get absolute binding times.

            :param random_numbers: Optional array of random numbers in [0,1) to convert
                to binding times. If None, generates new random numbers and refills
                the internal buffer.
            :type random_numbers: np.ndarray, optional
            :return: Array of binding times if random_numbers was provided, otherwise None
            :rtype: np.ndarray | None
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

            # Calculate binding time = -ln(X) / (rate * sites) - timestep/2
            # All operations done in-place for memory efficiency
            binding_time_list = np.log(binding_time_list, out=binding_time_list)
            denominator = (
                self.run.macro_params.micro_params.bind_rate_tPA
                * self.run.macro_params.micro_params.binding_sites
            ).magnitude
            binding_time_list = np.divide(
                binding_time_list, denominator, out=binding_time_list
            )
            # Use subtract instead of negating then adding (equivalent but cleaner)
            binding_time_list = np.subtract(
                -self.run.macro_params.time_step.magnitude / 2,
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
            """Provide the next batch of binding times from the internal buffer.

            Serves binding times from the internal buffer, automatically refilling it
            when necessary. The running time is inconsistent (fast most of the time,
            slow when refilling), but amortized performance is excellent.

            This is the preferred method for obtaining binding times during simulation
            as it provides optimal performance through batching.

            :param count: Number of binding times requested
            :type count: int
            :return: Array of binding times (excluding current time, which must be added by caller)
            :rtype: np.ndarray
            :raises: May raise memory errors if count exceeds available memory
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
                # Not enough binding times left in buffer - need to span across refill
                remaining = self.period - self.pointer
                # Copy remaining times from old buffer
                out[:remaining] = self.list[self.pointer :]
                # Refill the internal buffer
                self.fill_list()
                # Copy additional needed times from newly filled buffer
                needed_from_new = count - remaining
                out[remaining:] = self.list[:needed_from_new]
                # Update pointer to position after times we just served
                self.pointer = needed_from_new
            return out

    def unbind_by_degradation(self, m: np.ndarray, current_time: float):
        """Unbind molecules whose fibers have just degraded (macro unbind).

        When a fiber degrades, all bound molecules are forcibly unbound and enter
        a waiting period before they can move freely. During this waiting period,
        they can only move to already-degraded neighboring edges (restricted movement).

        This method:
        - Updates bound status to False for affected molecules
        - Sets the unbound_by_degradation flag to True
        - Calculates waiting time = current_time + average_bound_time - time_step/2
        - Sets binding_time to infinity (cannot rebind to degraded fiber)
        - Updates the total_macro_unbinds counter

        :param m: Boolean array mask of molecules to unbind. Typically molecules
            where bound=True and fiber_status[location] < current_time
        :type m: np.ndarray
        :param current_time: The current simulation time
        :type current_time: float
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
            + self.run.macro_params.average_bound_time.magnitude
            - self.run.macro_params.time_step.magnitude / 2
        )
        # Set the binding/unbinding time of the selected molecules to infinite
        # The fiber in their current location just degraded, so they can never
        # bind to it.
        self.binding_time[m] = float("inf")

        # Record MACRO_UNBOUND (unbinding by degradation) events
        mol_ids = np.where(m)[0]
        locations = self.location[m]
        rows, ranks = np.unravel_index(
            locations,
            (self.run.macro_params.rows, self.run.macro_params.full_row),
        )
        events = np.empty(count, dtype=_BIND_EVENT_DTYPE)
        events["Simulation Time Elapsed"] = current_time
        events["tPA Molecule Index"] = mol_ids
        events["Molecule New Status"] = MolStatus.MACRO_UNBOUND.value
        events["Grid Location Row"] = rows
        events["Grid Location Rank"] = ranks
        self._bind_events.append(events)

    def unbind_by_time(self, m: np.ndarray, current_time: float):
        """Unbind molecules whose binding duration has expired.

        When a molecule's binding_time is reached, it unbinds stochastically.
        With probability ``forced_unbind``, it becomes MICRO_UNBOUND (forced
        unbinding): the molecule has unrestricted movement but cannot rebind
        until its waiting period expires. Otherwise, it becomes UNBOUND
        (kinetic unbinding): the molecule can immediately rebind with a newly
        generated binding time.

        This method:

        - Updates bound status to False for affected molecules
        - Clears the unbound_by_degradation flag
        - For forced unbinds (MICRO_UNBOUND): sets waiting time and
          binding_time to infinity
        - For kinetic unbinds (UNBOUND): generates new binding times for
          immediate rebinding
        - Updates the total_micro_unbinds counter

        :param m: Boolean array mask of molecules to unbind. Typically molecules
            where bound=True and binding_time < current_time
        :type m: np.ndarray
        :param current_time: The current simulation time
        :type current_time: float
        """
        # Count how many molecules need to be unbound
        # that is, how many elements of `m` are `True`
        count = np.count_nonzero(m)
        # If none, then return without any more work
        if count == 0:
            return
        # Update binding state: no longer bound
        self.bound = self.bound & ~m
        # Clear the unbound_by_degradation flag (micro-unbound molecules have
        # unrestricted movement, unlike macro-unbound)
        self.unbound_by_degradation = self.unbound_by_degradation & ~m
        # Determine which molecules experience forced unbinding (MICRO_UNBOUND)
        # vs kinetic unbinding (UNBOUND)
        forced = np.full(self.run.macro_params.total_molecules, False, dtype=np.bool_)
        if self.run.macro_params.duplicate_fortran:
            # Use pre-generated random numbers in Fortran mode
            forced[m] = (
                self.random_numbers[RandomDraw.MICRO_UNBIND][m]
                <= self.run.macro_params.forced_unbind
            )
        else:
            # Generate fresh random numbers in native mode
            forced[m] = self.rng.random(count) <= self.run.macro_params.forced_unbind
        # Forced unbinds (MICRO_UNBOUND): cannot rebind until waiting period expires
        self.waiting_time[forced] = (
            current_time
            + self.run.macro_params.average_bound_time.magnitude
            - self.run.macro_params.time_step.magnitude / 2
        )
        self.binding_time[forced] = float("inf")
        num_forced = np.count_nonzero(forced)
        self.total_micro_unbinds += num_forced
        # Kinetic unbinds (UNBOUND): molecules can immediately attempt to rebind
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

        # Record MICRO_UNBOUND (forced unbinding) events
        if num_forced > 0:
            mol_ids = np.where(forced)[0]
            locations = self.location[forced]
            rows, ranks = np.unravel_index(
                locations,
                (self.run.macro_params.rows, self.run.macro_params.full_row),
            )
            events = np.empty(num_forced, dtype=_BIND_EVENT_DTYPE)
            events["Simulation Time Elapsed"] = current_time
            events["tPA Molecule Index"] = mol_ids
            events["Molecule New Status"] = MolStatus.MICRO_UNBOUND.value
            events["Grid Location Row"] = rows
            events["Grid Location Rank"] = ranks
            self._bind_events.append(events)

        # Record UNBOUND (kinetic unbinding) events
        non_forced = m & ~forced
        num_non_forced = count - num_forced
        if num_non_forced > 0:
            mol_ids = np.where(non_forced)[0]
            locations = self.location[non_forced]
            rows, ranks = np.unravel_index(
                locations,
                (self.run.macro_params.rows, self.run.macro_params.full_row),
            )
            events = np.empty(num_non_forced, dtype=_BIND_EVENT_DTYPE)
            events["Simulation Time Elapsed"] = current_time
            events["tPA Molecule Index"] = mol_ids
            events["Molecule New Status"] = MolStatus.UNBOUND.value
            events["Grid Location Row"] = rows
            events["Grid Location Rank"] = ranks
            self._bind_events.append(events)

    def expire_waiting_period(self, current_time: float):
        """Record UNBOUND events for molecules whose waiting period has expired.

        After a molecule is macro-unbound (fiber degraded) or micro-unbound
        (forced unbinding), it enters a waiting period during which it cannot
        rebind. When the waiting period expires, the Fortran code records an
        UNBOUND event (macro_diffuse_into_and_along__external.f90, lines
        1021-1039). This method replicates that behavior.

        For each expired molecule, this method:
        - Clears the ``unbound_by_degradation`` flag (no longer restricted)
        - Resets ``waiting_time`` to 0 (so the event fires exactly once)
        - Records an UNBOUND event in ``_bind_events``

        :param current_time: The current simulation time
        :type current_time: float
        """
        expired = (
            ~self.bound
            & (self.waiting_time > 0)
            & (self.waiting_time <= current_time)
        )
        count = np.count_nonzero(expired)
        if count == 0:
            return

        # Clear stale state
        self.unbound_by_degradation = self.unbound_by_degradation & ~expired
        self.waiting_time[expired] = 0

        # Record UNBOUND events
        mol_ids = np.where(expired)[0]
        locations = self.location[expired]
        rows, ranks = np.unravel_index(
            locations,
            (self.run.macro_params.rows, self.run.macro_params.full_row),
        )
        events = np.empty(count, dtype=_BIND_EVENT_DTYPE)
        events["Simulation Time Elapsed"] = current_time
        events["tPA Molecule Index"] = mol_ids
        events["Molecule New Status"] = MolStatus.UNBOUND.value
        events["Grid Location Row"] = rows
        events["Grid Location Rank"] = ranks
        self._bind_events.append(events)

    def find_unbinding_time(
        self, unbinding_time_bin: np.ndarray, current_time: float
    ) -> np.ndarray:
        """Calculate simulation time at which bound molecules will unbind.

        Uses microscale simulation data to determine unbinding times. The microscale
        provides a lookup table (unbinding_time) indexed by bins 0-100. This method
        samples uniformly from that distribution using linear interpolation.

        Mathematical approach:
        Think of unbinding_time as a function f(x) where unbinding_time[100*x] = f(x).
        We want to draw uniformly from the range of f, i.e., find f(r) where r~U(0,1).

        If 100*r is an integer i, we simply return unbinding_time[i].

        If 100*r is not an integer, then f(r) lies in the interval::

            (f(floor(100*r)), f(ceil(100*r)))

        We define a linear function g(x) such that::

            g(floor(100*r)) = f(floor(100*r))
            g(ceil(100*r)) = f(ceil(100*r))

        and use g(r) to approximate f(r).

        Equivalently, if i is an integer such that 100i < 100r < 100(i+1) and
        λ is defined such that 100i + λ = 100r, then::

            f(r) ≈ (1-λ)*f(i/100) + λ*f((i+1)/100)

        :param unbinding_time_bin: Array of values in [0, 100] representing which bin
            of the microscale unbinding time distribution to sample from
        :type unbinding_time_bin: np.ndarray
        :param current_time: The current simulation time to add to relative unbinding times
        :type current_time: float
        :return: Array of absolute unbinding times (current_time + relative_time - time_step/2)
        :rtype: np.ndarray
        """
        interp = np.interp(
            unbinding_time_bin,
            self.xp,
            self.run.data.macroscale_in.bin_edge_tpa_leaving_time,
        )
        return interp + (current_time - self.run.macro_params.time_step.magnitude / 2)

    def find_lysis_time(
        self,
        m: np.ndarray,
        unbinding_time_bin: np.ndarray,
        current_time: float,
        count: int,
    ) -> np.ndarray:
        """Calculate the time at which a bound molecule causes fiber lysis.

        Uses microscale simulation data to determine if and when a bound molecule
        will cause its fiber to degrade. Interpolates from the lysis_time data
        stored in the Run object.

        :param m: Boolean array mask indicating which molecules to calculate lysis times for
        :type m: np.ndarray
        :param unbinding_time_bin: The unbinding time bin index for each molecule (0-100)
        :type unbinding_time_bin: np.ndarray
        :param current_time: The current simulation time
        :type current_time: float
        :param count: The number of molecules being processed (should equal np.count_nonzero(m))
        :type count: int
        :return: Array of lysis times for each molecule. Returns infinity if lysis does not occur.
        :rtype: np.ndarray
        """
        # Draw random lysis event from uniform distribution
        if self.run.macro_params.duplicate_fortran:
            lysis_time_bin = self.random_numbers[RandomDraw.LYSIS_TIME][m]
        else:
            lysis_time_bin = self.rng.random(count)
        # Scale to match the number of microscale simulation runs
        lysis_time_bin = lysis_time_bin * (
            self.run.macro_params.micro_params.micro_simulations / 100
        )
        # Initialize all lysis times to infinity (no lysis by default)
        interp = np.full(count, float("inf"), dtype=np.double)
        # Convert unbinding bins to integer indices
        unbinding_time_bin = unbinding_time_bin.astype(int)
        # Get the number of lysis events recorded for each unbinding bin
        total_lyses = self.run.data.macroscale_in.binned_fiber_degraded[
            unbinding_time_bin
        ]
        # Check if lysis actually occurs (random draw < number of recorded lyses)
        lysis_happens = lysis_time_bin < total_lyses
        # Keep infinity for molecules that don't cause lysis
        # TODO: This line could probably be removed. Test!
        interp[~lysis_happens] = float("inf")

        # TODO(bpaynter): Performance bottleneck - this loop could potentially be
        #                 improved with 2D interpolation (scipy.interpolate.interp2d)
        #                 or a custom Numba/Cython kernel. Estimated 10-20% speedup
        #                 possible. Low priority as this runs infrequently.
        for i in np.arange(count)[lysis_happens]:
            interp[i] = np.interp(
                lysis_time_bin[i],
                np.arange(total_lyses[i]),
                self.run.data.macroscale_in.binned_fiber_degrade_time[
                    : total_lyses[i], unbinding_time_bin[i]
                ],
            )
        return interp + (current_time - self.run.macro_params.time_step.magnitude / 2)

    def bind(self, m: np.ndarray, current_time: float):
        """Bind molecules to fibers at their current locations.

        Updates molecule state to bound, calculates unbinding times from microscale
        data, and determines lysis times. If lysis occurs, updates the fiber_status
        to reflect when the fiber will degrade. Handles multiple molecules binding
        to the same fiber by taking the minimum lysis time.

        :param m: Boolean array mask indicating which molecules should bind
        :type m: np.ndarray
        :param current_time: The current simulation time
        :type current_time: float
        """
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

        # Record fiber degrade events where lysis_time actually changes fiber_status
        fiber_degrade_batch = []
        for i in range(count):
            if lysis_time[i] < float("inf"):
                new_val = min(self.fiber_status[locations[i]], lysis_time[i])
                if new_val < self.fiber_status[locations[i]]:
                    self.fiber_status[locations[i]] = new_val
                    fiber_degrade_batch.append(
                        (current_time, *np.unravel_index(
                            int(locations[i]),
                            (self.run.macro_params.rows,
                             self.run.macro_params.full_row),
                        ), new_val)
                    )
                # fiber_status unchanged if new_val == old_val (same min)
        # CRITICAL - DO NOT VECTORIZE:
        # We must use a loop here to handle the case where multiple molecules bind to
        # the same fiber in the same timestep. The vectorized approach below would allow
        # the last molecule's lysis time to overwrite earlier ones, potentially replacing
        # a lower (earlier) lysis time with a higher (later) one. Using min() in a loop
        # ensures we always keep the earliest lysis time for each fiber.
        #
        # INCORRECT (race condition):
        #   self.fiber_status[locations] = np.fmin(self.fiber_status[locations], lysis_time)
        # TODO: Investigate whether this could be prevented by first sorting by
        #       decreasing lysis time since last write wins.

        # Record tpa_bind_events for all molecules that bound
        mol_ids = np.where(m)[0]
        rows, ranks = np.unravel_index(
            locations,
            (self.run.macro_params.rows, self.run.macro_params.full_row),
        )
        events = np.empty(count, dtype=_BIND_EVENT_DTYPE)
        events["Simulation Time Elapsed"] = current_time
        events["tPA Molecule Index"] = mol_ids
        events["Molecule New Status"] = MolStatus.BOUND.value
        events["Grid Location Row"] = rows
        events["Grid Location Rank"] = ranks
        self._bind_events.append(events)

        # Record fiber degrade events
        if fiber_degrade_batch:
            fd_events = np.array(fiber_degrade_batch, dtype=_FIBER_DEGRADE_DTYPE)
            self._fiber_degrade_events.append(fd_events)

    def move_to_empty_edge(self, m: np.ndarray, current_time: float):
        """Move molecules to a random empty (degraded) neighboring edge.

        This restricted movement is used for molecules that are stuck after their
        fiber degraded. They can only move to neighboring edges that have already
        degraded (fiber_status < current_time), or stay in place if no such
        neighbors exist.

        :param m: Boolean array mask indicating which molecules should move
        :type m: np.ndarray
        :param current_time: The current simulation time
        :type current_time: float
        """
        count = np.count_nonzero(m)
        if count == 0:
            return

        self.total_restricted_moves += count

        # Get current locations of molecules that need to move
        current_locations = self.location[m]

        # Get all neighbors for each moving molecule
        neighborhoods = self.neighbors[current_locations]
        # Identify which neighbors are degraded (valid for restricted movement)
        valid_neighbors = self.fiber_status[neighborhoods] < current_time
        # Sort so valid neighbors come first (uses trick: ~bool sorts False before True)
        valid_neighborhood_index = np.argsort(~valid_neighbors, axis=1)
        valid_neighborhoods = np.take_along_axis(
            neighborhoods, valid_neighborhood_index, axis=1
        )
        # Add current location as option 0 (molecule can stay in place)
        valid_neighborhoods = np.append(
            current_locations.reshape(count, 1), valid_neighborhoods, axis=1
        )
        # Count how many valid neighbors each molecule has
        num_valid_neighbors = np.count_nonzero(valid_neighbors, axis=1)
        # Randomly select one of the valid options (including staying put)
        if self.run.macro_params.duplicate_fortran:
            neighbor = self.random_numbers[RandomDraw.RESTRICTED_MOVE][m] * (
                num_valid_neighbors + 1
            )
        else:
            neighbor = self.rng.random(count) * (num_valid_neighbors + 1)
        neighbor = neighbor.astype(int, copy=False)

        # Update molecule locations using the selected neighbor indices
        self.location[m] = valid_neighborhoods[np.full(count, True), neighbor]

    def find_still_stuck(self, m: np.ndarray, current_time: float) -> np.ndarray:
        """Identify molecules that are still waiting after being unbound by degradation.

        When a molecule is unbound by fiber degradation (macro unbind), it must
        wait for a period before it can move freely. This method identifies
        molecules that are still in this waiting period.

        :param m: Boolean array mask of molecules to check
        :type m: np.ndarray
        :param current_time: The current simulation time
        :type current_time: float
        :return: Boolean array indicating which molecules are still stuck waiting
        :rtype: np.ndarray
        """
        return (self.waiting_time > current_time) & self.unbound_by_degradation & m

    def unrestricted_move(self, free_to_move: np.ndarray, current_time: float):
        """Move molecules randomly to one of their eight neighboring edges.

        Performs unrestricted random walk movement for unbound molecules that are
        free to move. Each molecule moves to one of its 8 neighbors with equal
        probability. After moving, calculates new binding times for molecules that
        moved to intact fibers.

        :param free_to_move: Boolean array mask indicating which molecules can move freely
        :type free_to_move: np.ndarray
        :param current_time: The current simulation time
        :type current_time: float
        """
        if self.run.macro_params.duplicate_fortran:
            # Fortran mode: Transform random [0,1) to neighbor index [0,7]
            # First subtract (1 - moving_probability) to account for movement check
            neighbor = self.random_numbers[RandomDraw.MOVE][free_to_move]
            neighbor = neighbor - (1 - self.run.macro_params.moving_probability)
            # Normalize to [0,1) range
            neighbor = neighbor / self.run.macro_params.moving_probability
            # Scale to [0,8) range
            neighbor = neighbor * 8
            neighbor = neighbor.astype(int, copy=False)
        else:
            # Native mode: Directly generate random neighbor index [0,7]
            neighbor = self.rng.integers(8, size=np.count_nonzero(free_to_move))

        # Update molecule locations using pre-computed neighbor arrays
        self.location[free_to_move] = self.neighbors[
            self.location[free_to_move], neighbor
        ]
        self.total_regular_moves += np.count_nonzero(free_to_move)

        # TODO: Determine whether this code can be used if the 'duplicate_fortran'
        # options are removed.
        #
        # move_to_fiber = free_to_move & (self.m_fiber_status >= current_time)

        # Generate new binding times for all molecules that moved
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
        """Move molecules according to their state and constraints.

        Orchestrates molecule movement by first moving molecules that are still
        stuck (restricted movement to degraded edges only), then moving molecules
        that are free to move (unrestricted random walk). Also tracks when
        molecules reach the back row of the clot for the first time.

        :param m: Boolean array mask indicating which molecules should attempt to move
        :type m: np.ndarray
        :param current_time: The current simulation time
        :type current_time: float
        """
        # Find macro-unbound molecules still in their waiting period
        # (restricted movement — can only move to degraded edges)
        still_stuck_to_fiber = self.find_still_stuck(m, current_time)
        self.move_to_empty_edge(still_stuck_to_fiber, current_time)

        # Find those molecules that are unrestricted
        # and have them move to a random neighbor
        free_to_move = m & ~still_stuck_to_fiber
        self.unrestricted_move(free_to_move, current_time)

        # Only check if some molecules haven't reached the back row yet
        if self.number_reached_back_row < self.run.macro_params.total_molecules:
            # Identify molecules reaching the back row for the first time
            # Back row starts at index: (rows-1) * full_row
            first_time = ~self.reached_back_row & (
                self.location
                > (self.run.macro_params.rows - 1) * self.run.macro_params.full_row - 1
            )
            # Record the time they reached the back row
            self.time_to_reach_back_row[first_time] = current_time
            # Mark them as having reached the back row
            self.reached_back_row = self.reached_back_row | first_time
            # Update the counter
            self.number_reached_back_row += np.count_nonzero(first_time)

    def _grow_snapshot_buffers(self):
        """Double the capacity of the in-memory snapshot buffers.

        The number of saves is not known in advance for run-to-completion
        simulations (``total_time == 0``), so the snapshot buffers start small
        and are doubled on demand (amortized O(1) appends).
        :meth:`record_data_to_disk` later trims to the actual count on write.
        """
        old_capacity = self.snapshot_time.shape[0]
        new_capacity = max(1, old_capacity * 2)

        new_times = np.empty((new_capacity,), dtype=self.snapshot_time.dtype)
        new_times[:old_capacity] = self.snapshot_time
        self.snapshot_time = new_times

        new_locations = np.empty(
            (new_capacity, 2, self.run.macro_params.total_molecules),
            dtype=self.tpa_location_snapshot.dtype,
        )
        new_locations[:old_capacity] = self.tpa_location_snapshot
        self.tpa_location_snapshot = new_locations

    def save_data(self, current_time):
        """Save the current simulation state as a v2.0.0 snapshot.

        Stores molecule locations as (row, rank) coordinates and the simulation
        time at the next save interval index. Increments the save interval counter.

        :param current_time: The current simulation time to record
        :type current_time: float
        """
        # Grow the buffers if they are full (run-to-completion runs do not know
        # their length up front).
        if self.current_save_interval >= self.snapshot_time.shape[0]:
            self._grow_snapshot_buffers()
        rows, ranks = np.unravel_index(
            self.location,
            (self.run.macro_params.rows, self.run.macro_params.full_row),
        )
        self.tpa_location_snapshot[self.current_save_interval, 0, :] = rows
        self.tpa_location_snapshot[self.current_save_interval, 1, :] = ranks
        self.snapshot_time[self.current_save_interval] = current_time
        self.current_save_interval += 1

    def record_data_to_disk(self):
        """Write all collected simulation data to HDF5 via the DataStore.

        Writes snapshot_time, tpa_location_snapshot, tpa_bind_events,
        fiber_degrade_time, and tpa_transit_time datasets to the per-simulation
        macroscale_out group. Called at the end of the simulation.
        """
        sim_view = self.run.data.macroscale_out[self.sim_number]
        n = self.current_save_interval

        # snapshot_time
        ds = sim_view.snapshot_time
        ds.resize((n,))
        ds[:] = self.snapshot_time[:n]

        # tpa_location_snapshot
        ds = sim_view.tpa_location_snapshot
        ds.resize((n, 2, self.run.macro_params.total_molecules))
        ds[:] = self.tpa_location_snapshot[:n]

        # tpa_bind_events
        if self._bind_events:
            all_events = np.concatenate(self._bind_events)
            ds = sim_view.tpa_bind_events
            ds.resize((len(all_events),))
            ds[:] = all_events

        # fiber_degrade_time
        if self._fiber_degrade_events:
            all_events = np.concatenate(self._fiber_degrade_events)
            ds = sim_view.fiber_degrade_time
            ds.resize((len(all_events),))
            ds[:] = all_events

        # tpa_transit_time
        reached = self.time_to_reach_back_row[self.reached_back_row]
        if len(reached) > 0:
            ds = sim_view.tpa_transit_time
            ds.resize((len(reached),))
            ds[:] = reached

        self.logger.info(
            f"Wrote macroscale_out for sim {self.sim_number}: "
            f"{n} snapshots, {sum(len(e) for e in self._bind_events)} bind events, "
            f"{sum(len(e) for e in self._fiber_degrade_events)} fiber degrade events, "
            f"{len(reached)} transit times."
        )

    def go(self):
        """Execute the main simulation loop.

        Runs the macroscale simulation for the configured number of timesteps.
        At each timestep:

        1. Unbind molecules whose fibers degraded or whose binding time expired
        2. Determine which molecules should bind or move
        3. Resolve conflicts when both binding and movement are possible
        4. Execute binding and movement operations
        5. Save data at configured intervals
        6. Check termination conditions (all fibers degraded)

        Terminates early if all fibers have degraded. Saves final state and
        statistics at completion.
        """
        # Save initial state at t=0
        self.save_data(0)

        # total_time == 0 is a sentinel meaning "run until all fibers degrade,
        # no matter how long that takes" rather than stopping after a fixed
        # number of timesteps.  In that mode the step count is unknown, so we
        # iterate without bound and rely on the all-degraded check (performed
        # at each save point below) to terminate.
        run_to_completion = self.run.macro_params.total_time_steps == 0
        if run_to_completion:
            step_iterator = itertools.count()
        else:
            step_iterator = np.arange(self.run.macro_params.total_time_steps)

        current_time = 0.0
        for ts in tqdm(step_iterator, mininterval=2):
            # Calculate current simulation time
            current_time = ts * self.run.macro_params.time_step.magnitude
            if self.run.macro_params.duplicate_fortran:
                # Fortran starts at timestep 1 instead of 0
                current_time += self.run.macro_params.time_step.magnitude

            # The Fortran code processes all events for one molecule at a time
            # This code processes each event for all molecules at the same time
            # This means that, in order to follow the Fortran code exactly,
            # we need to pre-generate all of the random numbers we need for the
            # timestep and use them in the order Fortran would
            # (by molecule or column, instead of by event or row).
            if self.run.macro_params.duplicate_fortran:
                self.random_numbers = np.empty(
                    (8, self.run.macro_params.total_molecules), np.float64
                )
                for i in range(8):
                    self.random_numbers[i] = self.rng.random(
                        self.run.macro_params.total_molecules
                    )

            # Cache fiber status at each molecule's location for efficiency
            self.m_fiber_status = self.fiber_status[self.location]

            # Unbind molecules whose fibers have degraded
            self.unbind_by_degradation(
                self.bound & (self.m_fiber_status < current_time), current_time
            )
            # Unbind molecules whose binding time has expired
            self.unbind_by_time(
                self.bound & (self.binding_time < current_time), current_time
            )

            # Record waiting period expirations
            # (MACRO_UNBOUND -> UNBOUND and MICRO_UNBOUND -> UNBOUND)
            self.expire_waiting_period(current_time)

            # Determine which molecules should bind:
            # - Not currently bound
            # - Binding time has arrived
            # - Fiber is still intact
            # - Waiting period has expired (applies to macro- and micro-unbound)
            should_bind = (
                ~self.bound
                & (self.binding_time < current_time)
                & (self.m_fiber_status > current_time)
                & (self.waiting_time < current_time)
            )
            # Determine which molecules should move
            if self.run.macro_params.duplicate_fortran:
                move_chance = self.random_numbers[RandomDraw.MOVE]
                # Fortran mode: same probability of moving,
                #       but which molecules are selected is different,
                #       e.g., move if x ∈ (0.8, 1.0) instead of move if x ∈ (0.0, 0.2)
                should_move = (
                    move_chance > 1 - self.run.macro_params.moving_probability
                ) & ~self.bound
            else:
                # Native mode: simpler, but equal probability comparison
                move_chance = self.rng.random(self.run.macro_params.total_molecules)
                should_move = (
                    move_chance < self.run.macro_params.moving_probability
                ) & ~self.bound

            # Resolve conflicts when both binding and moving are scheduled
            conflict = should_bind & should_move
            # The molecule has been selected to move in this timestep, but the molecule
            # is also due to bind sometime during this timestep since
            #        current_time - timestep <= binding_time < current_time.
            # The question is, which happens first, binding or moving.
            # Since the molecule is due to move, we choose a random moment in the
            # interval for that to happen.
            # If the move is scheduled for after the bind, the molecule binds and the
            # move is cancelled.
            # If the move is scheduled for before the bind, the molecule moves and a
            # new binding_time is chosen.
            # All of this is normalized to a proportion to simplify calculation.
            #
            # Example:
            #
            # timestep = 0.042 sec
            # current_time = 253.383 sec
            # current interval: [253.341, 253.383)
            # binding_time = 253.368 sec
            # threshold = (253.368 - 253.341) / 0.042 = 0.643
            #
            # This means that the binding will happen 64.3% of the way through the
            # current time interval. Since moving can happen at any time during the
            # interval, there is a 64.3% chance that moving will happen before binding.
            # So we draw a random number, to determine how far into the interval moving
            # happens. If that number is in the range [0, 0.643) then moving happens
            # (because it happens before binding could take place). If that number is in
            # the range [0.643, 1.0] then binding happens (because it happens before
            # moving could take place).
            #
            # Simulation Time:        (current_time - timestep)          (current_time)
            #     (seconds)     0            253.341                  |      253.383
            #                   |----- ... -----|-------------------------------|---...
            #                                   |------threshold------|
            #                                            0.643        |
            #                                                    binding_time
            #                                                      253.368
            #                                                         |
            #                                   |--------move---------|---bind--|
            # Random number:                    0                   0.643       1
            #
            # NOTE: The Fortran code does this same calculation, but reverses the
            #       interval.
            #       I.e., move if x ∈ (0.357, 1.0) instead of bind if x ∈ [0.643, 1.0).
            #
            if self.run.macro_params.duplicate_fortran:
                # Determine the threshold for binding (how far from the end of the
                # interval does binding happen?)
                threshold = (
                    current_time - self.binding_time[conflict]
                ) / self.run.macro_params.time_step.magnitude
                # Bind if moving happens in that ending part
                should_bind[conflict] = (
                    self.random_numbers[RandomDraw.CONFLICT_RESOLUTION][conflict]
                    <= threshold
                )
            else:
                # Determine the threshold for moving (how far from the beginning of the
                # interval does binding happen?)
                threshold = (
                    self.binding_time[conflict]
                    - (current_time - self.run.macro_params.time_step.magnitude)
                ) / self.run.macro_params.time_step.magnitude
                # Bind if movement would happen after binding.
                should_bind[conflict] = (
                    self.rng.random(np.count_nonzero(conflict)) >= threshold
                )
            # Whatever doesn't bind will move
            should_move[conflict] = ~should_bind[conflict]

            # Execute binding and movement operations
            self.bind(should_bind, current_time)
            self.move(should_move, current_time)

            # Save simulation state at regular intervals.  We pause here to
            # write data anyway, so this is also where we check progress and
            # the termination condition (all fibers degraded) — the only stop
            # condition for run-to-completion runs, and an early exit for
            # fixed-length runs.  Checking every timestep would be too
            # expensive; the save cadence is a natural, cheap place to do it.
            if (
                current_time
                >= self.run.macro_params.save_interval.magnitude * self.current_save_interval
            ):
                self.save_data(current_time)

                # Count how many fibers are still intact
                unlysed_fibers = np.count_nonzero(self.fiber_status > current_time)

                if unlysed_fibers == 0:
                    # All fibers degraded — terminate.
                    self.logger.info(
                        f"All fibers degraded after {current_time:.2f} sec. Terminating"
                    )
                    self.last_degrade_time = np.max(self.fiber_status)
                    break
                else:
                    # Log progress statistics
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
                        f"the back row ({reached_back_row_percent:.1f}% of total)."
                    )

        # Save final state at the actual last simulated time (which need not
        # land exactly on total_time) and write all data to disk.
        self.save_data(current_time)
        self.record_data_to_disk()

        # Log final statistics
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
