"""Forward-replay cursor over a fiber degrade-time event log.

Provides :class:`FiberReplayCursor`, which wraps the structured-array
event log stored in ``run.data.macroscale_out[sim].fiber_degrade_time``
and exposes 2-D ``(row, rank)`` state that can be advanced forward through
simulation time without materialising the full snapshot array at every
save point.
"""

__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from lysis.config.run import Run

__all__ = ["FiberReplayCursor"]


class FiberReplayCursor:
    """Forward-replay cursor for fiber degradation events.

    Reads the ``fiber_degrade_time`` event log from
    ``run.data.macroscale_out[sim]`` and grid parameters from
    ``run.macro_params``.  The internal state array uses native
    ``(row, rank)`` indexing — no Fortran edge-index conversion is
    needed.

    :param run: Run object supplying macro grid parameters and data.
    :type run: Run
    :param sim: Zero-based simulation index.
    :type sim: int
    """

    def __init__(self, run: "Run", sim: int) -> None:
        mp = run.macro_params
        rows = mp.rows
        full_row = mp.full_row
        xz_row = mp.xz_row
        empty_rows = mp.empty_rows

        # Read event log
        event_log = run.data.macroscale_out[sim].fiber_degrade_time[:]

        if event_log.size == 0:
            self._event_times = np.empty(0, dtype=np.float64)
            self._event_rows = np.empty(0, dtype=np.uint32)
            self._event_ranks = np.empty(0, dtype=np.uint32)
            self._degrade_times = np.empty(0, dtype=np.float64)
        else:
            self._event_times = event_log["Simulation Time Elapsed"]
            self._event_rows = event_log["Grid Location Row"]
            self._event_ranks = event_log["Grid Location Rank"]
            self._degrade_times = event_log["Fiber New Degrade Time"]

        # Build initial state: shape (rows, full_row)
        state = np.full((rows, full_row), np.inf, dtype=np.float64)
        # Empty rows are always "degraded"
        state[:empty_rows, :] = 0.0
        # Last row: xz edges are inf (already set), y-edge padding is NaN
        state[rows - 1, xz_row:] = np.nan

        self._initial_state = state.copy()
        self._state = state
        self._cursor_pos = 0
        self._current_time = 0.0
        self._rows = rows
        self._full_row = full_row

    def advance_to(self, t: float) -> None:
        """Apply all events with time <= *t*.  Forward-only.

        :param t: Target time in seconds.
        :type t: float
        :raises ValueError: If *t* < :attr:`current_time`.
        """
        if t < self._current_time:
            raise ValueError(
                f"Cannot advance backward: current_time={self._current_time}, "
                f"requested t={t}"
            )
        new_cutoff = int(np.searchsorted(self._event_times, t, side="right"))
        if new_cutoff > self._cursor_pos:
            s = slice(self._cursor_pos, new_cutoff)
            self._state[self._event_rows[s], self._event_ranks[s]] = (
                self._degrade_times[s]
            )
            self._cursor_pos = new_cutoff
        self._current_time = t

    def reset(self) -> None:
        """Reset to initial state (time 0)."""
        self._state[:] = self._initial_state
        self._cursor_pos = 0
        self._current_time = 0.0

    @property
    def state(self) -> np.ndarray:
        """2-D degrade-time array, shape ``(rows, full_row)``.  Read-only view."""
        view = self._state.view()
        view.flags.writeable = False
        return view

    @property
    def current_time(self) -> float:
        """Time the cursor has been advanced to."""
        return self._current_time
