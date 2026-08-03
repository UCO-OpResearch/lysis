"""MacroscaleSim: the orchestrator.

Owns all simulation state and the per-timestep loop. Delegates every
mode-dependent decision to its injected strategies/services; the only
duplicate_fortran checks left here are for molecule placement, neighborhood
structure, and the current-time timestep offset -- structural differences
that don't belong to Move/Bind/Unbind/ConflictResolution (see the start of
this refactor for why).
"""

import itertools
import logging
import os
from functools import partial

import numpy as np
from tqdm.auto import tqdm

from ..config.constants import CONST
from ..config.run import Run
from ..geometry.edge_grid import EdgeGrid, from_fortran_edge_index


class MacroscaleSim:
    def __init__(
        self,
        run: Run,
        sim_number: int,
        draws,
        binding_time_factory,
        move_strategy,
        bind_strategy,
        unbind_strategy,
        conflict_strategy,
    ):
        self.run = run
        self.sim_number = sim_number
        assert self.run.macro_params is not None

        self.logger = logging.getLogger(__name__)
        self.logger.debug("Initializing MacroscaleSim")

        self.draws = draws
        self.binding_time_factory = binding_time_factory
        self.move_strategy = move_strategy
        self.bind_strategy = bind_strategy
        self.unbind_strategy = unbind_strategy
        self.conflict_strategy = conflict_strategy

        self.edge_lookup = partial(
            np.ravel_multi_index,
            dims=(run.macro_params.rows, run.macro_params.full_row),
        )

        self.fiber_status = np.full(
            self.run.macro_params.rows * self.run.macro_params.full_row,
            float("inf"),
            dtype=np.float64,
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
                self.edge_lookup((self.run.macro_params.rows - 1, 3 * j))
            ] = False
        self.fiber_status[~self.real_fiber] = 0

        self.logger.debug("Precalculating neighbors.")
        if run.macro_params.duplicate_fortran:
            self.neighbors = EdgeGrid.generate_fortran_neighborhood_structure(run)
        else:
            self.neighbors = EdgeGrid.generate_neighborhood_structure(run)

        self.logger.info("Placing molecules on empty edges.")
        if run.macro_params.duplicate_fortran:
            location = self.draws.rng.random(run.macro_params.total_molecules)
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
            location_i = self.draws.rng.integers(
                run.macro_params.empty_rows,
                size=run.macro_params.total_molecules,
                dtype=np.short,
            )
            location_j = self.draws.rng.integers(
                run.macro_params.full_row,
                size=run.macro_params.total_molecules,
                dtype=np.short,
            )
        self.location = self.edge_lookup((location_i, location_j))

        self.m_fiber_status = None

        self.bound = np.full(run.macro_params.total_molecules, False, dtype=np.bool_)
        self.waiting_time = np.full(
            run.macro_params.total_molecules, 0, dtype=np.float64
        )
        self.binding_time = np.full(
            run.macro_params.total_molecules, float("inf"), dtype=np.float64
        )
        self.unbound_by_degradation = np.full(
            run.macro_params.total_molecules, 0, dtype=np.bool_
        )
        self.time_to_reach_back_row = np.full(
            run.macro_params.total_molecules, float("inf"), dtype=np.float64
        )
        self.reached_back_row = np.full(
            run.macro_params.total_molecules, False, dtype=np.bool_
        )
        self.xp = np.arange(CONST.TPA_LEAVE_TIME_BINS + 1)

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

        initial_saves = max(1, self.run.macro_params.number_of_saves)
        self.tpa_location_snapshot = np.empty(
            (initial_saves, 2, self.run.macro_params.total_molecules),
            dtype=np.uint32,
        )
        self.snapshot_time = np.empty((initial_saves,), dtype=np.float64)

        self._bind_events = []
        self._fiber_degrade_events = []

        self.logger.debug("Initialization complete.")

    def _grow_snapshot_buffers(self):
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
        if self.current_save_interval >= self.snapshot_time.shape[0]:
            self._grow_snapshot_buffers()
        rows, ranks = np.unravel_index(
            self.location, (self.run.macro_params.rows, self.run.macro_params.full_row)
        )
        self.tpa_location_snapshot[self.current_save_interval, 0, :] = rows
        self.tpa_location_snapshot[self.current_save_interval, 1, :] = ranks
        self.snapshot_time[self.current_save_interval] = current_time
        self.current_save_interval += 1

    def record_data_to_disk(self):
        sim_view = self.run.data.macroscale_out[self.sim_number]
        n = self.current_save_interval

        ds = sim_view.snapshot_time
        ds.resize((n,))
        ds[:] = self.snapshot_time[:n]

        ds = sim_view.tpa_location_snapshot
        ds.resize((n, 2, self.run.macro_params.total_molecules))
        ds[:] = self.tpa_location_snapshot[:n]

        if self._bind_events:
            all_events = np.concatenate(self._bind_events)
            ds = sim_view.tpa_bind_events
            ds.resize((len(all_events),))
            ds[:] = all_events

        if self._fiber_degrade_events:
            all_events = np.concatenate(self._fiber_degrade_events)
            ds = sim_view.fiber_degrade_time
            ds.resize((len(all_events),))
            ds[:] = all_events

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
        self.save_data(0)

        run_to_completion = self.run.macro_params.total_time_steps == 0
        if run_to_completion:
            step_iterator = itertools.count()
        else:
            step_iterator = np.arange(self.run.macro_params.total_time_steps)

        current_time = 0.0
        for ts in tqdm(step_iterator, mininterval=2):
            current_time = ts * self.run.macro_params.time_step.magnitude
            if self.run.macro_params.duplicate_fortran:
                current_time += self.run.macro_params.time_step.magnitude

            self.draws.begin_timestep(self.run.macro_params.total_molecules)

            self.m_fiber_status = self.fiber_status[self.location]

            self.unbind_strategy.unbind_by_degradation(
                self, self.bound & (self.m_fiber_status < current_time), current_time
            )
            self.unbind_strategy.unbind_by_time(
                self, self.bound & (self.binding_time < current_time), current_time
            )
            self.unbind_strategy.expire_waiting_period(self, current_time)

            should_bind = (
                ~self.bound
                & (self.binding_time < current_time)
                & (self.m_fiber_status > current_time)
                & (self.waiting_time < current_time)
            )
            should_move = self.move_strategy.should_move(
                self.run.macro_params.total_molecules, self.bound
            )

            conflict = should_bind & should_move
            should_bind[conflict] = self.conflict_strategy.resolve_conflict(
                conflict, current_time, self.binding_time
            )
            should_move[conflict] = ~should_bind[conflict]

            self.bind_strategy.bind(self, should_bind, current_time)
            self.move_strategy.move(self, should_move, current_time)

            if (
                current_time
                >= self.run.macro_params.save_interval.magnitude
                * self.current_save_interval
            ):
                self.save_data(current_time)

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
                        f"the back row ({reached_back_row_percent:.1f}% of total)."
                    )

        self.save_data(current_time)
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