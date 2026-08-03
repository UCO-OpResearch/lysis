"""RandomDrawSource: unifies where random draws come from in Native vs
Fortran mode. Every strategy calls .draw_masked()/.draw_unmasked()/
.integers() and never touches self.rng or self.random_numbers directly.
"""

import abc

import numpy as np

from ...config.constants import RandomDraw
from ...tools.kiss import KissRandomGenerator


class RandomDrawSource(abc.ABC):
    @abc.abstractmethod
    def begin_timestep(self, total_molecules: int) -> None:
        """Called once per timestep, before any draws are requested."""

    @abc.abstractmethod
    def draw_masked(self, kind: RandomDraw, mask: np.ndarray) -> np.ndarray:
        """One float64 draw in [0, 1) per True entry in mask."""

    @abc.abstractmethod
    def draw_unmasked(self, kind: RandomDraw, n: int) -> np.ndarray:
        """n float64 draws in [0, 1), for the whole population (e.g. should_move)."""

    @abc.abstractmethod
    def integers(self, high: int, n: int) -> np.ndarray:
        """n integers in [0, high). Native-only -- see FortranDrawSource."""


class NativeDrawSource(RandomDrawSource):
    def __init__(self, seed: int):
        self.rng = np.random.default_rng(seed=seed)

    def begin_timestep(self, total_molecules: int) -> None:
        pass  # native mode draws lazily, nothing to pre-generate

    def draw_masked(self, kind: RandomDraw, mask: np.ndarray) -> np.ndarray:
        return self.rng.random(int(np.count_nonzero(mask)))

    def draw_unmasked(self, kind: RandomDraw, n: int) -> np.ndarray:
        return self.rng.random(n)

    def integers(self, high: int, n: int) -> np.ndarray:
        return self.rng.integers(high, size=n)


class FortranDrawSource(RandomDrawSource):
    def __init__(self, seed: int):
        self.rng = KissRandomGenerator(seed & 0xFFFFFFFF)
        self.random_numbers: np.ndarray | None = None

    def begin_timestep(self, total_molecules: int) -> None:
        self.random_numbers = np.empty((8, total_molecules), dtype=np.float64)
        for i in range(8):
            self.random_numbers[i] = self.rng.random(total_molecules)

    def draw_masked(self, kind: RandomDraw, mask: np.ndarray) -> np.ndarray:
        return self.random_numbers[kind][mask]

    def draw_unmasked(self, kind: RandomDraw, n: int) -> np.ndarray:
        return self.random_numbers[kind]

    def integers(self, high: int, n: int) -> np.ndarray:
        raise NotImplementedError(
            "FortranDrawSource has no integer-draw path -- Fortran mode "
            "combines the move-probability check and neighbor pick into a "
            "single RandomDraw.MOVE draw instead. Reaching this means "
            "something is calling integers() on the wrong mode."
        )


def make_draw_source(duplicate_fortran: bool, seed: int) -> RandomDrawSource:
    return FortranDrawSource(seed) if duplicate_fortran else NativeDrawSource(seed)