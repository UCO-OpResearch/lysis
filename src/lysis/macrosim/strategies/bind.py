"""DefaultBindStrategy.

find_unbinding_time and find_lysis_time from the original are inlined
directly into bind().

state is expected to be whatever owns simulation state -- in practice, the
orchestrator.
"""

import numpy as np

from ...config.constants import CONST, MolStatus, RandomDraw
from ...dataio.dataspec import dataspec
from ..services.random_draw import RandomDrawSource
from .base import BindStrategy

_macro_out_spec = dataspec["v2.0.0"]["macroscale_out"]
_BIND_EVENT_DTYPE = _macro_out_spec.data["tpa_bind_events"].dtype
_FIBER_DEGRADE_DTYPE = _macro_out_spec.data["fiber_degrade_time"].dtype


class DefaultBindStrategy(BindStrategy):
    def __init__(self, draws: RandomDrawSource):
        self.draws = draws

    def bind(self, state, m: np.ndarray, current_time: float) -> None:
        count = np.count_nonzero(m)
        if count == 0:
            return

        state.bound = state.bound | m
        state.waiting_time[m] = 0
        state.total_binds += count

        # find_unbinding_time, inlined
        unbinding_time_bin = (
            self.draws.draw_masked(RandomDraw.UNBINDING_TIME, m)
            * CONST.TPA_LEAVE_TIME_BINS
        )
        state.binding_time[m] = np.interp(
            unbinding_time_bin,
            state.xp,
            state.run.data.macroscale_in.bin_edge_tpa_leaving_time,
        ) + (current_time - state.run.macro_params.time_step.magnitude / 2)

        # find_lysis_time, inlined
        lysis_time_bin = self.draws.draw_masked(RandomDraw.LYSIS_TIME, m)
        lysis_time_bin = lysis_time_bin * (
            state.run.macro_params.micro_params.micro_simulations
            / CONST.TPA_LEAVE_TIME_BINS
        )
        lysis_interp = np.full(count, float("inf"), dtype=np.double)
        unbinding_time_bin_int = unbinding_time_bin.astype(int)
        total_lyses = state.run.data.macroscale_in.binned_fiber_degraded[
            unbinding_time_bin_int
        ]
        lysis_happens = lysis_time_bin < total_lyses
        # TODO: This line could probably be removed. Test!
        lysis_interp[~lysis_happens] = float("inf")
        # TODO(bpaynter): Performance bottleneck - this loop could potentially be
        #                 improved with 2D interpolation (scipy.interpolate.interp2d)
        #                 or a custom Numba/Cython kernel. Estimated 10-20% speedup
        #                 possible. Low priority as this runs infrequently.
        for i in np.arange(count)[lysis_happens]:
            lysis_interp[i] = np.interp(
                lysis_time_bin[i],
                np.arange(total_lyses[i]),
                state.run.data.macroscale_in.binned_fiber_degrade_time[
                    : total_lyses[i], unbinding_time_bin_int[i]
                ],
            )
        lysis_time = lysis_interp + (
            current_time - state.run.macro_params.time_step.magnitude / 2
        )

        locations = state.location[m]

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
        fiber_degrade_batch = []
        for i in range(count):
            if lysis_time[i] < float("inf"):
                new_val = min(state.fiber_status[locations[i]], lysis_time[i])
                if new_val < state.fiber_status[locations[i]]:
                    state.fiber_status[locations[i]] = new_val
                    fiber_degrade_batch.append(
                        (
                            current_time,
                            *np.unravel_index(
                                int(locations[i]),
                                (
                                    state.run.macro_params.rows,
                                    state.run.macro_params.full_row,
                                ),
                            ),
                            new_val,
                        )
                    )

        mol_ids = np.where(m)[0]
        rows, ranks = np.unravel_index(
            locations, (state.run.macro_params.rows, state.run.macro_params.full_row)
        )
        events = np.empty(count, dtype=_BIND_EVENT_DTYPE)
        events["Simulation Time Elapsed"] = current_time
        events["tPA Molecule Index"] = mol_ids
        events["Molecule New Status"] = MolStatus.BOUND.value
        events["Grid Location Row"] = rows
        events["Grid Location Rank"] = ranks
        state._bind_events.append(events)

        if fiber_degrade_batch:
            fd_events = np.array(fiber_degrade_batch, dtype=_FIBER_DEGRADE_DTYPE)
            state._fiber_degrade_events.append(fd_events)