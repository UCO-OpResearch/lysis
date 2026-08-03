"""NativeConflictResolutionStrategy / FortranConflictResolutionStrategy.

Decides, for molecules eligible for both bind and move in the same
timestep, which one wins -- based on where in the current timestep interval
binding_time falls versus a random draw. Native and Fortran define the
threshold over opposite ends of the interval, so this is a genuine
behavioral divergence, not just a draw-source difference.
"""

import numpy as np

from ...config.constants import RandomDraw
from ..services.random_draw import RandomDrawSource
from .base import ConflictResolutionStrategy


class NativeConflictResolutionStrategy(ConflictResolutionStrategy):
    def __init__(self, draws: RandomDrawSource, time_step: float):
        self.draws = draws
        self.time_step = time_step

    def resolve_conflict(
        self, conflict: np.ndarray, current_time: float, binding_time: np.ndarray
    ) -> np.ndarray:
        threshold = (
            binding_time[conflict] - (current_time - self.time_step)
        ) / self.time_step
        draw = self.draws.draw_masked(RandomDraw.CONFLICT_RESOLUTION, conflict)
        return draw >= threshold


class FortranConflictResolutionStrategy(ConflictResolutionStrategy):
    def __init__(self, draws: RandomDrawSource, time_step: float):
        self.draws = draws
        self.time_step = time_step

    def resolve_conflict(
        self, conflict: np.ndarray, current_time: float, binding_time: np.ndarray
    ) -> np.ndarray:
        threshold = (current_time - binding_time[conflict]) / self.time_step
        draw = self.draws.draw_masked(RandomDraw.CONFLICT_RESOLUTION, conflict)
        return draw <= threshold