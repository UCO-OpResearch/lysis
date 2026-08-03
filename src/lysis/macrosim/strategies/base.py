"""Abstract base classes for every strategy category."""

import abc

import numpy as np


class MoveStrategy(abc.ABC):
    """Decides whether molecules move, and which neighbor they move to."""

    @abc.abstractmethod
    def should_move(self, total_molecules: int, bound: np.ndarray) -> np.ndarray:
        """Boolean mask: which molecules attempt an unrestricted move this timestep."""

    @abc.abstractmethod
    def find_still_stuck(
        self, state, m: np.ndarray, current_time: float
    ) -> np.ndarray: ...

    @abc.abstractmethod
    def move_to_empty_edge(
        self, state, m: np.ndarray, current_time: float
    ) -> None: ...

    @abc.abstractmethod
    def unrestricted_move(
        self, state, free_to_move: np.ndarray, current_time: float
    ) -> None: ...

    @abc.abstractmethod
    def move(self, state, m: np.ndarray, current_time: float) -> None: ...


class BindStrategy(abc.ABC):
    """Binds molecules to fibers and computes unbinding/lysis times."""

    @abc.abstractmethod
    def bind(self, state, m: np.ndarray, current_time: float) -> None:
        """Bind molecules in mask m: update bound state, binding_time,
        fiber_status (lysis), and record bind / fiber-degrade events."""


class UnbindStrategy(abc.ABC):
    """Handles all three unbind pathways."""

    @abc.abstractmethod
    def unbind_by_degradation(self, m: np.ndarray, current_time: float) -> None: ...

    @abc.abstractmethod
    def unbind_by_time(self, m: np.ndarray, current_time: float) -> None: ...

    @abc.abstractmethod
    def expire_waiting_period(self, current_time: float) -> None: ...

class ConflictResolutionStrategy(abc.ABC):
    """Decides, for molecules eligible for both bind and move in the same
    timestep, which one wins."""

    @abc.abstractmethod
    def resolve_conflict(
        self, conflict: np.ndarray, current_time: float, binding_time: np.ndarray
    ) -> np.ndarray:
        """Return the overridden should-bind value for each True entry in
        conflict. Native/Fortran differ here."""