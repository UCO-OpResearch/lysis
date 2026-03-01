"""Unit tests for :mod:`lysis.analysis.fiber_replay`.

Tests cover:

* FiberReplayCursor.__init__  — initial state shape, empty/fiber/padding values
* advance_to                  — forward replay through events, state correctness
* advance backward            — raises ValueError
* reset                       — returns to initial state
* state read-only             — view is not writeable
* empty event log             — cursor works with zero events
"""

import numpy as np
import pytest

from lysis.analysis.fiber_replay import FiberReplayCursor

from .conftest import _StubData, _StubMacroOut, _StubRun, _StubSimView, _build_event_log


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_run_and_cursor(small_macro, sim=0):
    """Build a _StubRun with synthetic event log and return (run, cursor)."""
    rows = small_macro.rows
    cols = small_macro.cols
    n_sims = small_macro.macro_simulations
    snapshot_time = np.array([0.0, 75.0, 150.0, 225.0, 300.0])
    event_log = _build_event_log(rows, cols)
    views = [
        _StubSimView(snapshot_time.copy(), event_log.copy())
        for _ in range(n_sims)
    ]
    data = _StubData(_StubMacroOut(views))
    run = _StubRun(small_macro, data=data)
    cursor = FiberReplayCursor(run, sim)
    return run, cursor


# ===========================================================================
# TestInit
# ===========================================================================


class TestInit:
    """Tests for FiberReplayCursor.__init__."""

    def test_state_shape(self, small_macro):
        """State has shape (rows, full_row)."""
        _, cursor = _make_run_and_cursor(small_macro)
        assert cursor.state.shape == (small_macro.rows, small_macro.full_row)

    def test_empty_rows_are_zero(self, small_macro):
        """Empty rows (row 0) are all 0.0."""
        _, cursor = _make_run_and_cursor(small_macro)
        np.testing.assert_array_equal(
            cursor.state[: small_macro.empty_rows, :], 0.0
        )

    def test_fiber_rows_are_inf(self, small_macro):
        """Fiber rows (rows 1 to rows-2) are all np.inf."""
        _, cursor = _make_run_and_cursor(small_macro)
        fibrous = cursor.state[small_macro.empty_rows : small_macro.rows - 1, :]
        assert np.all(np.isinf(fibrous))
        assert np.all(fibrous > 0)

    def test_last_row_xz_region_is_inf(self, small_macro):
        """Last row, first xz_row positions, are np.inf."""
        _, cursor = _make_run_and_cursor(small_macro)
        last_row_xz = cursor.state[small_macro.rows - 1, : small_macro.xz_row]
        assert np.all(np.isinf(last_row_xz))

    def test_last_row_padding_is_nan(self, small_macro):
        """Last row, positions xz_row onwards, are NaN."""
        _, cursor = _make_run_and_cursor(small_macro)
        padding = cursor.state[small_macro.rows - 1, small_macro.xz_row :]
        assert np.all(np.isnan(padding))

    def test_current_time_is_zero(self, small_macro):
        """Initial current_time is 0."""
        _, cursor = _make_run_and_cursor(small_macro)
        assert cursor.current_time == 0.0


# ===========================================================================
# TestAdvanceTo
# ===========================================================================


class TestAdvanceTo:
    """Tests for FiberReplayCursor.advance_to."""

    def test_no_change_at_time_zero(self, small_macro):
        """Advancing to t=0 does not alter state (no events at t <= 0)."""
        _, cursor = _make_run_and_cursor(small_macro)
        initial = cursor.state.copy()
        cursor.advance_to(0.0)
        np.testing.assert_array_equal(cursor.state, initial)

    def test_first_batch_applied(self, small_macro):
        """After advance_to(50), 8 fibers in row 1 have degrade_time=50."""
        _, cursor = _make_run_and_cursor(small_macro)
        cursor.advance_to(50.0)
        row1 = cursor.state[1, :]
        assert np.count_nonzero(row1 == 50.0) == 8

    def test_state_after_75(self, small_macro):
        """After advance_to(75), only first batch events (t=50) applied."""
        _, cursor = _make_run_and_cursor(small_macro)
        cursor.advance_to(75.0)
        # Row 1: all 50
        assert np.all(cursor.state[1, :] == 50.0)
        # Row 2: still inf (events at t=100 not yet applied)
        assert np.all(np.isinf(cursor.state[2, :]))

    def test_state_after_150(self, small_macro):
        """After advance_to(150), first two batches applied."""
        _, cursor = _make_run_and_cursor(small_macro)
        cursor.advance_to(150.0)
        # Row 1: all 50
        assert np.all(cursor.state[1, :] == 50.0)
        # Row 2: 7 positions with degrade_time=100, 1 still inf (edge (2,6) at t=175)
        row2 = cursor.state[2, :]
        assert np.count_nonzero(row2 == 100.0) == 7
        assert np.count_nonzero(np.isinf(row2)) == 1

    def test_state_after_300(self, small_macro):
        """After advance_to(300), all events applied."""
        _, cursor = _make_run_and_cursor(small_macro)
        cursor.advance_to(300.0)
        # No inf remaining in rows 1-3
        for row in range(1, 4):
            assert not np.any(np.isinf(cursor.state[row, :]))

    def test_advance_updates_current_time(self, small_macro):
        """current_time tracks the latest advance_to target."""
        _, cursor = _make_run_and_cursor(small_macro)
        cursor.advance_to(75.0)
        assert cursor.current_time == 75.0
        cursor.advance_to(200.0)
        assert cursor.current_time == 200.0

    def test_advance_to_same_time_is_noop(self, small_macro):
        """Advancing to current_time again does not raise or change state."""
        _, cursor = _make_run_and_cursor(small_macro)
        cursor.advance_to(75.0)
        state_before = cursor.state.copy()
        cursor.advance_to(75.0)
        np.testing.assert_array_equal(cursor.state, state_before)


# ===========================================================================
# TestAdvanceBackwardRaises
# ===========================================================================


class TestAdvanceBackwardRaises:
    """Tests that advance_to raises ValueError for backward time."""

    def test_backward_raises(self, small_macro):
        """advance_to(t) with t < current_time raises ValueError."""
        _, cursor = _make_run_and_cursor(small_macro)
        cursor.advance_to(100.0)
        with pytest.raises(ValueError, match="Cannot advance backward"):
            cursor.advance_to(50.0)


# ===========================================================================
# TestReset
# ===========================================================================


class TestReset:
    """Tests for FiberReplayCursor.reset."""

    def test_reset_restores_initial_state(self, small_macro):
        """After advance + reset, state matches initial state."""
        _, cursor = _make_run_and_cursor(small_macro)
        initial = cursor.state.copy()
        cursor.advance_to(300.0)
        cursor.reset()
        np.testing.assert_array_equal(cursor.state, initial)

    def test_reset_restores_current_time(self, small_macro):
        """After reset, current_time is 0."""
        _, cursor = _make_run_and_cursor(small_macro)
        cursor.advance_to(150.0)
        cursor.reset()
        assert cursor.current_time == 0.0

    def test_can_advance_after_reset(self, small_macro):
        """After reset, can advance again and get the same results."""
        _, cursor = _make_run_and_cursor(small_macro)
        cursor.advance_to(75.0)
        state_at_75 = cursor.state.copy()
        cursor.reset()
        cursor.advance_to(75.0)
        np.testing.assert_array_equal(cursor.state, state_at_75)


# ===========================================================================
# TestStateReadOnly
# ===========================================================================


class TestStateReadOnly:
    """Tests that state property returns a read-only view."""

    def test_state_not_writeable(self, small_macro):
        """Writing to the state array raises ValueError."""
        _, cursor = _make_run_and_cursor(small_macro)
        with pytest.raises(ValueError):
            cursor.state[0, 0] = 999.0


# ===========================================================================
# TestEmptyEventLog
# ===========================================================================


class TestEmptyEventLog:
    """Tests that FiberReplayCursor works with an empty event log."""

    def test_empty_log_initial_state(self, small_macro):
        """With no events, state is all initial values."""
        n_sims = small_macro.macro_simulations
        snapshot_time = np.array([0.0, 75.0, 150.0, 225.0, 300.0])
        empty_log = np.array([], dtype=np.dtype(
            [
                ("Simulation Time Elapsed", np.float64),
                ("Grid Location Row", np.uint32),
                ("Grid Location Rank", np.uint32),
                ("Fiber New Degrade Time", np.float64),
            ]
        ))
        views = [
            _StubSimView(snapshot_time.copy(), empty_log.copy())
            for _ in range(n_sims)
        ]
        data = _StubData(_StubMacroOut(views))
        run = _StubRun(small_macro, data=data)
        cursor = FiberReplayCursor(run, 0)
        assert cursor.state.shape == (small_macro.rows, small_macro.full_row)
        assert cursor.current_time == 0.0

    def test_empty_log_advance_does_not_change_state(self, small_macro):
        """With no events, advancing does not change state."""
        n_sims = small_macro.macro_simulations
        snapshot_time = np.array([0.0, 75.0, 150.0, 225.0, 300.0])
        empty_log = np.array([], dtype=np.dtype(
            [
                ("Simulation Time Elapsed", np.float64),
                ("Grid Location Row", np.uint32),
                ("Grid Location Rank", np.uint32),
                ("Fiber New Degrade Time", np.float64),
            ]
        ))
        views = [
            _StubSimView(snapshot_time.copy(), empty_log.copy())
            for _ in range(n_sims)
        ]
        data = _StubData(_StubMacroOut(views))
        run = _StubRun(small_macro, data=data)
        cursor = FiberReplayCursor(run, 0)
        initial = cursor.state.copy()
        cursor.advance_to(300.0)
        np.testing.assert_array_equal(cursor.state, initial)
        assert cursor.current_time == 300.0
