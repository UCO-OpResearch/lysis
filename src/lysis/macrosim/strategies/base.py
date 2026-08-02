"""Abstract base classes for every strategy category.

STATUS: placeholder.

Each ABC's abstract methods are the real "contract" that keeps this codebase
coupled-on-purpose: a new subclass that doesn't implement one raises at
instantiation time, not three timesteps into a run.

Target shape:

    class MoveStrategy(abc.ABC):
        @abc.abstractmethod
        def should_move(self, total_molecules, bound): ...
        @abc.abstractmethod
        def pick_neighbor(self, free_to_move): ...

    class BindStrategy(abc.ABC):
        ...

    class UnbindStrategy(abc.ABC):
        ...

"""
import abc

import numpy as np


class MoveStrategy(abc.ABC):
    """Decides whether molecules move, and which neighbor they move to."""

    @abc.abstractmethod
    def should_move(self, total_molecules: int, bound: np.ndarray) -> np.ndarray:
        """Boolean mask: which molecules attempt an unrestricted move this timestep."""

    @abc.abstractmethod
    def pick_neighbor(self, free_to_move: np.ndarray) -> np.ndarray:
        """Integer neighbor index in [0, 8) for each True entry in free_to_move."""

    @abc.abstractmethod
    def pick_restricted_neighbor(
        self, m: np.ndarray, num_valid_neighbors: np.ndarray
    ) -> np.ndarray:
        """Integer index selecting among valid degraded neighbors (plus
        'stay in place') for restricted movement, per molecule in m."""


class BindStrategy(abc.ABC):
    """Binds molecules to fibers and computes unbinding/lysis times."""

    @abc.abstractmethod
    def bind(self, m: np.ndarray, current_time: float) -> None:
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