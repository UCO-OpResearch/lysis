"""NativeMoveStrategy / FortranMoveStrategy.

find_still_stuck/move_to_empty_edge/unrestricted_move/move are currently
identical in orchestration logic between the two classes -- only the neighbor-
selection formulas inlined inside move_to_empty_edge and unrestricted_move
differ. This duplication is intentional: MoveStrategy has no default
implementations, so future variants are forced to explicitly decide
whether to reuse this logic or write their own, rather than silently
inheriting it.

state is expected to be whatever owns simulation state (location, neighbors,
fiber_status, counters, binding_time_factory, etc.) -- in practice, the
orchestrator.
"""

import numpy as np

from ...config.constants import RandomDraw
from ..services.random_draw import RandomDrawSource
from .base import MoveStrategy


class NativeMoveStrategy(MoveStrategy):
    def __init__(self, draws: RandomDrawSource, moving_probability: float):
        self.draws = draws
        self.moving_probability = moving_probability

    def should_move(self, total_molecules: int, bound: np.ndarray) -> np.ndarray:
        move_chance = self.draws.draw_unmasked(RandomDraw.MOVE, total_molecules)
        return (move_chance < self.moving_probability) & ~bound

    def find_still_stuck(
        self, state, m: np.ndarray, current_time: float
    ) -> np.ndarray:
        return (state.waiting_time > current_time) & state.unbound_by_degradation & m

    def move_to_empty_edge(self, state, m: np.ndarray, current_time: float) -> None:
        count = np.count_nonzero(m)
        if count == 0:
            return

        state.total_restricted_moves += count
        current_locations = state.location[m]
        neighborhoods = state.neighbors[current_locations]
        valid_neighbors = state.fiber_status[neighborhoods] < current_time
        valid_neighborhood_index = np.argsort(~valid_neighbors, axis=1)
        valid_neighborhoods = np.take_along_axis(
            neighborhoods, valid_neighborhood_index, axis=1
        )
        valid_neighborhoods = np.append(
            current_locations.reshape(count, 1), valid_neighborhoods, axis=1
        )
        num_valid_neighbors = np.count_nonzero(valid_neighbors, axis=1)

        # pick_restricted_neighbor, inlined
        neighbor = self.draws.draw_masked(RandomDraw.RESTRICTED_MOVE, m) * (
            num_valid_neighbors + 1
        )
        neighbor = neighbor.astype(int, copy=False)

        state.location[m] = valid_neighborhoods[np.full(count, True), neighbor]

    def unrestricted_move(
        self, state, free_to_move: np.ndarray, current_time: float
    ) -> None:
        # pick_neighbor, inlined
        count = int(np.count_nonzero(free_to_move))
        neighbor = self.draws.integers(8, count)

        state.location[free_to_move] = state.neighbors[
            state.location[free_to_move], neighbor
        ]
        state.total_regular_moves += np.count_nonzero(free_to_move)

        num_move_to_fiber = np.count_nonzero(free_to_move)
        if num_move_to_fiber > 0:
            # binding_time_factory: still the pending _BindingTimeFactory
            # question from earlier -- ported verbatim, not re-architected.
            state.binding_time[free_to_move] = (
                current_time + state.binding_time_factory.next(num_move_to_fiber)
            )

    def move(self, state, m: np.ndarray, current_time: float) -> None:
        still_stuck_to_fiber = self.find_still_stuck(state, m, current_time)
        self.move_to_empty_edge(state, still_stuck_to_fiber, current_time)

        free_to_move = m & ~still_stuck_to_fiber
        self.unrestricted_move(state, free_to_move, current_time)

        if state.number_reached_back_row < state.run.macro_params.total_molecules:
            first_time = ~state.reached_back_row & (
                state.location
                > (state.run.macro_params.rows - 1) * state.run.macro_params.full_row
                - 1
            )
            state.time_to_reach_back_row[first_time] = current_time
            state.reached_back_row = state.reached_back_row | first_time
            state.number_reached_back_row += np.count_nonzero(first_time)


class FortranMoveStrategy(MoveStrategy):
    def __init__(self, draws: RandomDrawSource, moving_probability: float):
        self.draws = draws
        self.moving_probability = moving_probability

    def should_move(self, total_molecules: int, bound: np.ndarray) -> np.ndarray:
        move_chance = self.draws.draw_unmasked(RandomDraw.MOVE, total_molecules)
        return (move_chance > 1 - self.moving_probability) & ~bound

    def find_still_stuck(
        self, state, m: np.ndarray, current_time: float
    ) -> np.ndarray:
        return (state.waiting_time > current_time) & state.unbound_by_degradation & m

    def move_to_empty_edge(self, state, m: np.ndarray, current_time: float) -> None:
        count = np.count_nonzero(m)
        if count == 0:
            return

        state.total_restricted_moves += count
        current_locations = state.location[m]
        neighborhoods = state.neighbors[current_locations]
        valid_neighbors = state.fiber_status[neighborhoods] < current_time
        valid_neighborhood_index = np.argsort(~valid_neighbors, axis=1)
        valid_neighborhoods = np.take_along_axis(
            neighborhoods, valid_neighborhood_index, axis=1
        )
        valid_neighborhoods = np.append(
            current_locations.reshape(count, 1), valid_neighborhoods, axis=1
        )
        num_valid_neighbors = np.count_nonzero(valid_neighbors, axis=1)

        # pick_restricted_neighbor, inlined -- identical formula to Native's
        # version, only the draw source differs. Not a bug.
        neighbor = self.draws.draw_masked(RandomDraw.RESTRICTED_MOVE, m) * (
            num_valid_neighbors + 1
        )
        neighbor = neighbor.astype(int, copy=False)

        state.location[m] = valid_neighborhoods[np.full(count, True), neighbor]

    def unrestricted_move(
        self, state, free_to_move: np.ndarray, current_time: float
    ) -> None:
        # pick_neighbor, inlined -- reuses the same RandomDraw.MOVE draw
        # already consumed by should_move; Fortran folds the move-probability
        # check and the neighbor choice into one value per molecule.
        neighbor = self.draws.draw_masked(RandomDraw.MOVE, free_to_move)
        neighbor = neighbor - (1 - self.moving_probability)
        neighbor = neighbor / self.moving_probability
        neighbor = neighbor * 8
        neighbor = neighbor.astype(int, copy=False)

        state.location[free_to_move] = state.neighbors[
            state.location[free_to_move], neighbor
        ]
        state.total_regular_moves += np.count_nonzero(free_to_move)

        num_move_to_fiber = np.count_nonzero(free_to_move)
        if num_move_to_fiber > 0:
            draw = self.draws.draw_masked(
                RandomDraw.BINDING_TIME_WHEN_MOVING, free_to_move
            )
            state.binding_time[free_to_move] = (
                current_time + state.binding_time_factory.fill_list(draw)
            )

    def move(self, state, m: np.ndarray, current_time: float) -> None:
        still_stuck_to_fiber = self.find_still_stuck(state, m, current_time)
        self.move_to_empty_edge(state, still_stuck_to_fiber, current_time)

        free_to_move = m & ~still_stuck_to_fiber
        self.unrestricted_move(state, free_to_move, current_time)

        if state.number_reached_back_row < state.run.macro_params.total_molecules:
            first_time = ~state.reached_back_row & (
                state.location
                > (state.run.macro_params.rows - 1) * state.run.macro_params.full_row
                - 1
            )
            state.time_to_reach_back_row[first_time] = current_time
            state.reached_back_row = state.reached_back_row | first_time
            state.number_reached_back_row += np.count_nonzero(first_time)