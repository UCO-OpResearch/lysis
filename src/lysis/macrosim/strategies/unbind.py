"""NativeUnbindStrategy / FortranUnbindStrategy.

unbind_by_degradation and expire_waiting_period have no random draws at all
-- their bodies are identical between the two classes, so future variants 
are forced to decide for themselves rather than silently inheriting.

unbind_by_time is where the real divergence lives: the forced-unbind draw
(MICRO_UNBIND) is identical in both modes via draw_masked, but the
kinetic-rebind binding time goes through binding_time_factory.next() (Native)
vs .fill_list(draw) (Fortran).

state is expected to be whatever owns simulation state -- in practice, the
orchestrator.
"""

import numpy as np

from ...config.constants import CONST, MolStatus, RandomDraw
from ...dataio.dataspec import dataspec
from ..services.random_draw import RandomDrawSource
from .base import UnbindStrategy

_macro_out_spec = dataspec["v2.0.0"]["macroscale_out"]
_BIND_EVENT_DTYPE = _macro_out_spec.data["tpa_bind_events"].dtype


class NativeUnbindStrategy(UnbindStrategy):
    def __init__(self, draws: RandomDrawSource):
        self.draws = draws

    def unbind_by_degradation(
        self, state, m: np.ndarray, current_time: float
    ) -> None:
        count = np.count_nonzero(m)
        if count == 0:
            return

        state.bound = state.bound & ~m
        state.unbound_by_degradation = state.unbound_by_degradation | m
        state.total_macro_unbinds += count
        state.waiting_time[m] = (
            current_time
            + state.run.macro_params.average_bound_time.magnitude
            - state.run.macro_params.time_step.magnitude / 2
        )
        state.binding_time[m] = float("inf")

        mol_ids = np.where(m)[0]
        locations = state.location[m]
        rows, ranks = np.unravel_index(
            locations, (state.run.macro_params.rows, state.run.macro_params.full_row)
        )
        events = np.empty(count, dtype=_BIND_EVENT_DTYPE)
        events["Simulation Time Elapsed"] = current_time
        events["tPA Molecule Index"] = mol_ids
        events["Molecule New Status"] = MolStatus.MACRO_UNBOUND.value
        events["Grid Location Row"] = rows
        events["Grid Location Rank"] = ranks
        state._bind_events.append(events)

    def unbind_by_time(self, state, m: np.ndarray, current_time: float) -> None:
        count = np.count_nonzero(m)
        if count == 0:
            return

        state.bound = state.bound & ~m
        state.unbound_by_degradation = state.unbound_by_degradation & ~m

        forced = np.full(state.run.macro_params.total_molecules, False, dtype=np.bool_)
        forced[m] = (
            self.draws.draw_masked(RandomDraw.MICRO_UNBIND, m)
            <= state.run.macro_params.forced_unbind
        )

        state.waiting_time[forced] = (
            current_time
            + state.run.macro_params.average_bound_time.magnitude
            - state.run.macro_params.time_step.magnitude / 2
        )
        state.binding_time[forced] = float("inf")
        num_forced = np.count_nonzero(forced)
        state.total_micro_unbinds += num_forced

        if num_forced < count:
            state.binding_time[m & ~forced] = (
                current_time + state.binding_time_factory.next(count - num_forced)
            )

        if num_forced > 0:
            mol_ids = np.where(forced)[0]
            locations = state.location[forced]
            rows, ranks = np.unravel_index(
                locations,
                (state.run.macro_params.rows, state.run.macro_params.full_row),
            )
            events = np.empty(num_forced, dtype=_BIND_EVENT_DTYPE)
            events["Simulation Time Elapsed"] = current_time
            events["tPA Molecule Index"] = mol_ids
            events["Molecule New Status"] = MolStatus.MICRO_UNBOUND.value
            events["Grid Location Row"] = rows
            events["Grid Location Rank"] = ranks
            state._bind_events.append(events)

        non_forced = m & ~forced
        num_non_forced = count - num_forced
        if num_non_forced > 0:
            mol_ids = np.where(non_forced)[0]
            locations = state.location[non_forced]
            rows, ranks = np.unravel_index(
                locations,
                (state.run.macro_params.rows, state.run.macro_params.full_row),
            )
            events = np.empty(num_non_forced, dtype=_BIND_EVENT_DTYPE)
            events["Simulation Time Elapsed"] = current_time
            events["tPA Molecule Index"] = mol_ids
            events["Molecule New Status"] = MolStatus.UNBOUND.value
            events["Grid Location Row"] = rows
            events["Grid Location Rank"] = ranks
            state._bind_events.append(events)

    def expire_waiting_period(self, state, current_time: float) -> None:
        expired = (
            ~state.bound
            & (state.waiting_time > 0)
            & (state.waiting_time <= current_time)
        )
        count = np.count_nonzero(expired)
        if count == 0:
            return

        state.unbound_by_degradation = state.unbound_by_degradation & ~expired
        state.waiting_time[expired] = 0

        mol_ids = np.where(expired)[0]
        locations = state.location[expired]
        rows, ranks = np.unravel_index(
            locations, (state.run.macro_params.rows, state.run.macro_params.full_row)
        )
        events = np.empty(count, dtype=_BIND_EVENT_DTYPE)
        events["Simulation Time Elapsed"] = current_time
        events["tPA Molecule Index"] = mol_ids
        events["Molecule New Status"] = MolStatus.UNBOUND.value
        events["Grid Location Row"] = rows
        events["Grid Location Rank"] = ranks
        state._bind_events.append(events)


class FortranUnbindStrategy(UnbindStrategy):
    def __init__(self, draws: RandomDrawSource):
        self.draws = draws

    def unbind_by_degradation(
        self, state, m: np.ndarray, current_time: float
    ) -> None:
        # Identical to NativeUnbindStrategy's version -- no randomness
        # involved here at all. Duplicated deliberately, not a bug.
        count = np.count_nonzero(m)
        if count == 0:
            return

        state.bound = state.bound & ~m
        state.unbound_by_degradation = state.unbound_by_degradation | m
        state.total_macro_unbinds += count
        state.waiting_time[m] = (
            current_time
            + state.run.macro_params.average_bound_time.magnitude
            - state.run.macro_params.time_step.magnitude / 2
        )
        state.binding_time[m] = float("inf")

        mol_ids = np.where(m)[0]
        locations = state.location[m]
        rows, ranks = np.unravel_index(
            locations, (state.run.macro_params.rows, state.run.macro_params.full_row)
        )
        events = np.empty(count, dtype=_BIND_EVENT_DTYPE)
        events["Simulation Time Elapsed"] = current_time
        events["tPA Molecule Index"] = mol_ids
        events["Molecule New Status"] = MolStatus.MACRO_UNBOUND.value
        events["Grid Location Row"] = rows
        events["Grid Location Rank"] = ranks
        state._bind_events.append(events)

    def unbind_by_time(self, state, m: np.ndarray, current_time: float) -> None:
        count = np.count_nonzero(m)
        if count == 0:
            return

        state.bound = state.bound & ~m
        state.unbound_by_degradation = state.unbound_by_degradation & ~m

        forced = np.full(state.run.macro_params.total_molecules, False, dtype=np.bool_)
        forced[m] = (
            self.draws.draw_masked(RandomDraw.MICRO_UNBIND, m)
            <= state.run.macro_params.forced_unbind
        )

        state.waiting_time[forced] = (
            current_time
            + state.run.macro_params.average_bound_time.magnitude
            - state.run.macro_params.time_step.magnitude / 2
        )
        state.binding_time[forced] = float("inf")
        num_forced = np.count_nonzero(forced)
        state.total_micro_unbinds += num_forced

        if num_forced < count:
            draw = self.draws.draw_masked(
                RandomDraw.BINDING_TIME_WHEN_UNBINDING, m & ~forced
            )
            state.binding_time[m & ~forced] = (
                current_time + state.binding_time_factory.fill_list(draw)
            )

        if num_forced > 0:
            mol_ids = np.where(forced)[0]
            locations = state.location[forced]
            rows, ranks = np.unravel_index(
                locations,
                (state.run.macro_params.rows, state.run.macro_params.full_row),
            )
            events = np.empty(num_forced, dtype=_BIND_EVENT_DTYPE)
            events["Simulation Time Elapsed"] = current_time
            events["tPA Molecule Index"] = mol_ids
            events["Molecule New Status"] = MolStatus.MICRO_UNBOUND.value
            events["Grid Location Row"] = rows
            events["Grid Location Rank"] = ranks
            state._bind_events.append(events)

        non_forced = m & ~forced
        num_non_forced = count - num_forced
        if num_non_forced > 0:
            mol_ids = np.where(non_forced)[0]
            locations = state.location[non_forced]
            rows, ranks = np.unravel_index(
                locations,
                (state.run.macro_params.rows, state.run.macro_params.full_row),
            )
            events = np.empty(num_non_forced, dtype=_BIND_EVENT_DTYPE)
            events["Simulation Time Elapsed"] = current_time
            events["tPA Molecule Index"] = mol_ids
            events["Molecule New Status"] = MolStatus.UNBOUND.value
            events["Grid Location Row"] = rows
            events["Grid Location Rank"] = ranks
            state._bind_events.append(events)

    def expire_waiting_period(self, state, current_time: float) -> None:
        # Identical to NativeUnbindStrategy's version -- no randomness here.
        expired = (
            ~state.bound
            & (state.waiting_time > 0)
            & (state.waiting_time <= current_time)
        )
        count = np.count_nonzero(expired)
        if count == 0:
            return

        state.unbound_by_degradation = state.unbound_by_degradation & ~expired
        state.waiting_time[expired] = 0

        mol_ids = np.where(expired)[0]
        locations = state.location[expired]
        rows, ranks = np.unravel_index(
            locations, (state.run.macro_params.rows, state.run.macro_params.full_row)
        )
        events = np.empty(count, dtype=_BIND_EVENT_DTYPE)
        events["Simulation Time Elapsed"] = current_time
        events["tPA Molecule Index"] = mol_ids
        events["Molecule New Status"] = MolStatus.UNBOUND.value
        events["Grid Location Row"] = rows
        events["Grid Location Rank"] = ranks
        state._bind_events.append(events)