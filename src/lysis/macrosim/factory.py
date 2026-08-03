"""Composition root -- the only module that reads duplicate_fortran."""

from .services.binding_time import BindingTimeFactory
from .services.random_draw import FortranDrawSource, NativeDrawSource, RandomDrawSource
from .strategies.bind import DefaultBindStrategy
from .strategies.conflict_resolution import (
    FortranConflictResolutionStrategy,
    NativeConflictResolutionStrategy,
)
from .strategies.move import FortranMoveStrategy, NativeMoveStrategy
from .strategies.unbind import FortranUnbindStrategy, NativeUnbindStrategy
from .macroscale_sim import MacroscaleSim

# Adding a new variant to any category is just one new entry here -- nothing
# else in this file needs to change.
_DRAW_SOURCES = {
    "native": NativeDrawSource,
    "fortran": FortranDrawSource,
}

_MOVE_STRATEGIES = {
    "native": NativeMoveStrategy,
    "fortran": FortranMoveStrategy,
}

_UNBIND_STRATEGIES = {
    "native": NativeUnbindStrategy,
    "fortran": FortranUnbindStrategy,
}

_CONFLICT_STRATEGIES = {
    "native": NativeConflictResolutionStrategy,
    "fortran": FortranConflictResolutionStrategy,
}

def make_draw_source(mode: str, seed: int) -> RandomDrawSource:
    return _DRAW_SOURCES[mode](seed)

def make_move_strategy(mode: str, draws, moving_probability: float):
    return _MOVE_STRATEGIES[mode](draws, moving_probability)


def make_bind_strategy(draws):
    # No registry -- there's currently only one variant, and nothing selects
    # between variants of Bind at all (see DefaultBindStrategy's docstring).
    return DefaultBindStrategy(draws)


def make_unbind_strategy(mode: str, draws):
    return _UNBIND_STRATEGIES[mode](draws)


def make_conflict_strategy(mode: str, draws, time_step: float):
    return _CONFLICT_STRATEGIES[mode](draws, time_step)


def build_simulation(run, sim_number=0) -> MacroscaleSim:
    # Placeholder until mode piping has a source -- for now, just read the run's duplicate_fortran flag.
    mode = "fortran" if run.macro_params.duplicate_fortran else "native"
    seed = run.macro_params.seed_as_int()

    draws = make_draw_source(mode, seed)
    binding_time_factory = BindingTimeFactory(run, draws.rng)

    return MacroscaleSim(
        run=run,
        sim_number=sim_number,
        draws=draws,
        binding_time_factory=binding_time_factory,
        move_strategy=make_move_strategy(
            mode, draws, run.macro_params.moving_probability
        ),
        bind_strategy=make_bind_strategy(draws),
        unbind_strategy=make_unbind_strategy(mode, draws),
        conflict_strategy=make_conflict_strategy(
            mode, draws, run.macro_params.time_step.magnitude
        ),
    )