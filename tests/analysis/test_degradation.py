"""Unit tests for :mod:`lysis.analysis.degradation`.

Tests cover:

* find_degraded_fraction          — fraction, shape, empty-edge exclusion
* find_degradation_marker_frames  — shape, validity, monotonicity, exact frames
* find_degradation_marker_times   — shape, units (minutes), exact times
* degradation_rates               — shape, sign, exact rate value
* mean_degradation_rate           — shape, sign, lag time
* calculate_time_row_exposed      — shape, row-0 zeroed, units, monotonicity
* find_degradation_fronts         — structure, y-distances, time ordering
* mean_front_velocity             — return type, positive velocity
* find_row_deg_fraction           — shape, value range, all-degraded case
* find_front                      — shape, range, all-degraded case
* fiber_degradation_linear_extrapolation — shape, range, initial/final values
* plot_degradation_percent        — returns Figure with axes
* plot_front_degradation          — returns Figure with axes
* get_unbind_amounts              — parses log files correctly
* get_processing_time             — parses log files correctly
* get_total_binds                 — parses log files correctly

All synthetic data is defined in conftest.py.  Key known values:

    degraded_fraction[sim] ≈ [0, 8/29, 15/29, 22/29, 1]
    marker_frames[sim, :]  =  [0, 1, 2, 3, 4]
    marker_times[sim, :]   =  [0, 1.25, 2.5, 3.75, 5.0] minutes
    degradation_rate (25%-75%) ≈ 14/29 / 150 * 60 ≈ 0.1931 fraction/min
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.figure import Figure

from lysis.analysis.degradation import (
    calculate_time_row_exposed,
    degradation_rates,
    fiber_degradation_linear_extrapolation,
    find_degradation_fronts,
    find_degradation_marker_frames,
    find_degradation_marker_times,
    find_degraded_fraction,
    find_front,
    find_row_deg_fraction,
    get_processing_time,
    get_total_binds,
    get_unbind_amounts,
    mean_degradation_rate,
    mean_front_velocity,
    plot_degradation_percent,
    plot_front_degradation,
)

# Re-import config constants used for assertions
from .conftest import FRONT_THRESHOLD, PERCENT_MARKERS, SLOPE_PAIRS


# ===========================================================================
# find_degraded_fraction
# ===========================================================================


class TestFindDegradedFraction:
    """Tests for :func:`find_degraded_fraction`."""

    def test_returns_list_of_length_n_sims(self, stub_run):
        """Return value is a list with one array per simulation."""
        result = find_degraded_fraction(stub_run)
        assert isinstance(result, list)
        assert len(result) == stub_run.macro_params.macro_simulations

    def test_each_element_has_shape_n_saves(self, stub_run):
        """Each element has shape (n_saves,)."""
        result = find_degraded_fraction(stub_run)
        n_saves = len(stub_run.data.macroscale_out[0].snapshot_time[:])
        for arr in result:
            assert arr.shape == (n_saves,)

    def test_initial_fraction_is_zero(self, degraded_fraction):
        """Degraded fraction is 0 at the first save point (no fiber has degraded yet)."""
        for arr in degraded_fraction:
            assert arr[0] == pytest.approx(0.0)

    def test_final_fraction_is_one(self, degraded_fraction):
        """Degraded fraction reaches 1 at the last save point (all fibers degraded)."""
        for arr in degraded_fraction:
            assert arr[-1] == pytest.approx(1.0)

    def test_values_non_decreasing(self, degraded_fraction):
        """Degraded fraction must not decrease over time."""
        for arr in degraded_fraction:
            assert np.all(np.diff(arr) >= -1e-12)

    def test_values_in_zero_to_one(self, degraded_fraction):
        """All degraded-fraction values lie in [0, 1]."""
        for arr in degraded_fraction:
            assert np.all(arr >= -1e-12)
            assert np.all(arr <= 1.0 + 1e-12)

    def test_known_fraction_at_frame_one(self, degraded_fraction):
        """At save-point 1 exactly 8 of 29 fibers are degraded."""
        expected = 8 / 29
        for arr in degraded_fraction:
            assert arr[1] == pytest.approx(expected)

    def test_known_fraction_at_frame_two(self, degraded_fraction):
        """At save-point 2 exactly 15 of 29 fibers are degraded."""
        expected = 15 / 29
        for arr in degraded_fraction:
            assert arr[2] == pytest.approx(expected)


# ===========================================================================
# find_degradation_marker_frames
# ===========================================================================


class TestFindDegradationMarkerFrames:
    """Tests for :func:`find_degradation_marker_frames`."""

    def test_returns_array_of_correct_shape(self, stub_run):
        """Returns an integer array of shape (n_sims, n_markers)."""
        result = find_degradation_marker_frames(stub_run, PERCENT_MARKERS)
        n_sims = stub_run.macro_params.macro_simulations
        assert result.shape == (n_sims, len(PERCENT_MARKERS))

    def test_frame_indices_are_valid(self, stub_run):
        """All frame indices are non-negative and within the save range."""
        n_saves = len(stub_run.data.macroscale_out[0].snapshot_time[:])
        result = find_degradation_marker_frames(stub_run, PERCENT_MARKERS)
        assert np.all(result >= 0)
        assert np.all(result < n_saves)

    def test_marker_frames_are_non_decreasing_per_sim(self, stub_run):
        """Frame indices must be non-decreasing for each simulation."""
        result = find_degradation_marker_frames(stub_run, PERCENT_MARKERS)
        for sim_frames in result:
            assert np.all(np.diff(sim_frames) >= 0)

    def test_known_marker_frames(self, stub_run):
        """With the synthetic data, every milestone is first hit at frames [0,1,2,3,4]."""
        result = find_degradation_marker_frames(stub_run, PERCENT_MARKERS)
        expected = np.array([0, 1, 2, 3, 4])
        for sim_frames in result:
            np.testing.assert_array_equal(sim_frames, expected)

    def test_zero_threshold_returns_frame_zero(self, stub_run):
        """A threshold of 0.0 is always satisfied at frame 0."""
        result = find_degradation_marker_frames(stub_run, [0.0])
        assert np.all(result[:, 0] == 0)


# ===========================================================================
# find_degradation_marker_times
# ===========================================================================


class TestFindDegradationMarkerTimes:
    """Tests for :func:`find_degradation_marker_times`."""

    def test_returns_array_of_correct_shape(self, stub_run):
        """Returns a float array of shape (n_sims, n_markers)."""
        result = find_degradation_marker_times(stub_run, PERCENT_MARKERS)
        n_sims = stub_run.macro_params.macro_simulations
        assert result.shape == (n_sims, len(PERCENT_MARKERS))

    def test_values_are_in_minutes(self, stub_run, marker_times):
        """All times must be strictly less than max tsave in seconds."""
        max_tsave_sec = max(
            stub_run.data.macroscale_out[sim].snapshot_time[:][-1]
            for sim in range(stub_run.macro_params.macro_simulations)
        )
        # If units were seconds, values would be near max_tsave_sec (e.g. 300),
        # but in minutes they must be much smaller.
        assert np.all(marker_times <= max_tsave_sec / 60 + 1e-9)

    def test_known_marker_times_in_minutes(self, marker_times):
        """With tsave=[0,75,150,225,300]s and frames=[0,1,2,3,4] the times
        must be [0, 1.25, 2.5, 3.75, 5.0] minutes."""
        expected = np.array([0.0, 1.25, 2.5, 3.75, 5.0])
        for sim_times in marker_times:
            np.testing.assert_allclose(sim_times, expected)

    def test_times_are_non_decreasing(self, marker_times):
        """Marker times must be non-decreasing for each simulation."""
        for sim_times in marker_times:
            assert np.all(np.diff(sim_times) >= -1e-12)


# ===========================================================================
# degradation_rates
# ===========================================================================


class TestDegradationRates:
    """Tests for :func:`degradation_rates`."""

    def test_returns_array_of_correct_shape(self, stub_run):
        """Returns a float array of shape (n_sims, n_slope_pairs)."""
        result = degradation_rates(stub_run, SLOPE_PAIRS, PERCENT_MARKERS)
        n_sims = stub_run.macro_params.macro_simulations
        assert result.shape == (n_sims, len(SLOPE_PAIRS))

    def test_rates_are_positive(self, rates):
        """Degradation rates are positive (fraction can only increase)."""
        assert np.all(rates > 0)

    def test_known_rate_value(self, rates):
        """25%-to-75% rate = (14/29) / 150s * 60 ≈ 0.1931 fraction/min."""
        expected = (14 / 29) / 150 * 60
        for sim_rates in rates:
            assert sim_rates[0] == pytest.approx(expected)

    def test_multiple_slope_pairs(self, stub_run):
        """Accepts multiple slope pairs and returns one column per pair."""
        pairs = [(0.0, 0.5), (0.25, 1.0)]
        result = degradation_rates(stub_run, pairs, PERCENT_MARKERS)
        assert result.shape == (stub_run.macro_params.macro_simulations, 2)


# ===========================================================================
# mean_degradation_rate
# ===========================================================================


class TestMeanDegradationRate:
    """Tests for :func:`mean_degradation_rate`."""

    def test_returns_three_arrays_of_length_n_sims(self, stub_run):
        """Returns a 3-tuple; each element is a 1-D array of length n_sims."""
        n_sims = stub_run.macro_params.macro_simulations
        dr, off, lag = mean_degradation_rate(stub_run)
        assert dr.shape == (n_sims,)
        assert off.shape == (n_sims,)
        assert lag.shape == (n_sims,)

    def test_degradation_rate_is_positive(self, stub_run):
        """Fitted slope must be positive (clot degrades over time)."""
        dr, _, _ = mean_degradation_rate(stub_run)
        assert np.all(dr > 0)

    def test_deg_start_time_is_non_negative(self, stub_run):
        """Lysis lag time must be non-negative."""
        _, _, lag = mean_degradation_rate(stub_run)
        assert np.all(lag >= 0)

    def test_deg_start_time_in_minutes(self, stub_run):
        """Lysis lag time must be less than the total simulation time in minutes."""
        total_min = stub_run.data.macroscale_out[0].snapshot_time[:][-1] / 60
        _, _, lag = mean_degradation_rate(stub_run)
        assert np.all(lag <= total_min)


# ===========================================================================
# calculate_time_row_exposed
# ===========================================================================


class TestCalculateTimeRowExposed:
    """Tests for :func:`calculate_time_row_exposed`."""

    def test_returns_array_of_correct_shape(self, stub_run):
        """Returns array of shape (n_sims, rows-1, cols)."""
        rows = stub_run.macro_params.rows
        cols = stub_run.macro_params.cols
        n_sims = stub_run.macro_params.macro_simulations
        result = calculate_time_row_exposed(stub_run)
        assert result.shape == (n_sims, rows - 1, cols)

    def test_row_zero_is_always_zero(self, exposed_time):
        """Row 0 is always exposed at t=0 (no fibrin above it)."""
        np.testing.assert_array_equal(exposed_time[:, 0, :], 0.0)

    def test_values_are_non_negative(self, exposed_time):
        """All exposure times must be >= 0."""
        assert np.all(exposed_time >= 0.0)

    def test_values_are_in_minutes(self, stub_run, exposed_time):
        """Exposure times must be in minutes (< total simulation time in seconds)."""
        max_tsave_sec = max(
            stub_run.data.macroscale_out[sim].snapshot_time[:][-1]
            for sim in range(stub_run.macro_params.macro_simulations)
        )
        # Values in minutes must be less than max tsave in seconds
        assert np.all(exposed_time <= max_tsave_sec)

    def test_exposure_non_decreasing_across_rows(self, exposed_time):
        """Exposure time of each row must be >= that of the preceding row."""
        # Compare consecutive rows along axis=1
        assert np.all(np.diff(exposed_time, axis=1) >= -1e-12)

    def test_known_exposure_time_column_zero(self, exposed_time):
        """For column 0, the exposure times (min) follow the fiber degrade times.

        With the synthetic degrade schedule and to_fortran_edge_index indexing:
          row 1 is exposed at 50/60 min, row 2 at 100/60 min, row 3 at 175/60 min.
        """
        expected_rows = np.array([0.0, 50 / 60, 100 / 60, 175 / 60])
        for sim in range(exposed_time.shape[0]):
            np.testing.assert_allclose(
                exposed_time[sim, :, 0], expected_rows, rtol=1e-6
            )


# ===========================================================================
# find_degradation_fronts
# ===========================================================================


class TestFindDegradationFronts:
    """Tests for :func:`find_degradation_fronts`."""

    def test_returns_list_of_length_n_sims(self, deg_fronts, stub_run):
        """Outer list has one element per simulation."""
        assert len(deg_fronts) == stub_run.macro_params.macro_simulations

    def test_each_element_is_list_of_length_cols(self, deg_fronts, stub_run):
        """Each simulation's front data is a list with one array per column."""
        cols = stub_run.macro_params.cols
        for sim_fronts in deg_fronts:
            assert len(sim_fronts) == cols

    def test_each_column_array_has_two_rows(self, deg_fronts):
        """Each column array has shape (2, n_events): row 0 = times, row 1 = y-positions."""
        for sim_fronts in deg_fronts:
            for col_front in sim_fronts:
                assert col_front.ndim == 2
                assert col_front.shape[0] == 2

    def test_front_times_are_non_negative(self, deg_fronts):
        """Front event times (row 0) must be non-negative."""
        for sim_fronts in deg_fronts:
            for col_front in sim_fronts:
                assert np.all(col_front[0] >= 0.0)

    def test_y_distances_are_positive(self, deg_fronts):
        """y-distance values (row 1) must be > 0 (row 0 is never in the front)."""
        for sim_fronts in deg_fronts:
            for col_front in sim_fronts:
                assert np.all(col_front[1] > 0.0)

    def test_y_distances_are_multiples_of_pore_size(self, deg_fronts, stub_run):
        """y-distances must be integer multiples of pore_size (µm)."""
        pore_um = stub_run.macro_params.pore_size.to("microns").magnitude
        for sim_fronts in deg_fronts:
            for col_front in sim_fronts:
                ratios = col_front[1] / pore_um
                np.testing.assert_allclose(ratios, np.round(ratios), atol=1e-9)

    def test_three_front_events_per_column(self, deg_fronts, stub_run):
        """With the synthetic schedule every column exposes 3 rows (rows 1, 2, 3)."""
        for sim_fronts in deg_fronts:
            for col_front in sim_fronts:
                assert col_front.shape == (2, 3)


# ===========================================================================
# mean_front_velocity
# ===========================================================================


class TestMeanFrontVelocity:
    """Tests for :func:`mean_front_velocity`."""

    def test_returns_tuple_of_two_floats(self, stub_run):
        """Return value is a 2-tuple of Python floats."""
        result = mean_front_velocity(stub_run)
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], float)
        assert isinstance(result[1], float)

    def test_mean_velocity_is_positive(self, stub_run):
        """Mean front velocity must be positive (front moves away from source)."""
        mean_v, _ = mean_front_velocity(stub_run)
        assert mean_v > 0.0

    def test_std_velocity_is_non_negative(self, stub_run):
        """Std of front velocity must be non-negative."""
        _, std_v = mean_front_velocity(stub_run)
        assert std_v >= 0.0


# ===========================================================================
# find_row_deg_fraction
# ===========================================================================


class TestFindRowDegFraction:
    """Tests for :func:`find_row_deg_fraction`."""

    def test_returns_list_of_length_n_sims(self, stub_run):
        """Return value is a list with one array per simulation."""
        result = find_row_deg_fraction(stub_run)
        assert len(result) == stub_run.macro_params.macro_simulations

    def test_each_element_has_correct_shape(self, stub_run):
        """Each array has shape (n_saves, fiber_rows - 1)."""
        result = find_row_deg_fraction(stub_run)
        n_saves = len(stub_run.data.macroscale_out[0].snapshot_time[:])
        fiber_rows = stub_run.macro_params.fiber_rows
        for arr in result:
            assert arr.shape == (n_saves, fiber_rows - 1)

    def test_values_in_zero_to_one(self, row_deg):
        """All per-row undegraded fractions lie in [0, 1]."""
        for arr in row_deg:
            assert np.all(arr >= -1e-12)
            assert np.all(arr <= 1.0 + 1e-12)

    def test_all_rows_intact_at_first_save(self, row_deg):
        """At save-point 0 no fibers have degraded, so all fractions are 1."""
        for arr in row_deg:
            np.testing.assert_allclose(arr[0], 1.0)

    def test_first_row_fully_degraded_at_frame_one(self, row_deg):
        """At save-point 1, row 0 of the fibrin region is fully degraded."""
        # Row 0 corresponds to fibers 0-7 which all degrade at t=50 < tsave[1]=75.
        for arr in row_deg:
            assert arr[1, 0] == pytest.approx(0.0)


# ===========================================================================
# find_front
# ===========================================================================


class TestFindFront:
    """Tests for :func:`find_front`."""

    def test_returns_list_of_length_n_sims(self, stub_run):
        """Return value is a list with one array per simulation."""
        result = find_front(stub_run, FRONT_THRESHOLD)
        assert len(result) == stub_run.macro_params.macro_simulations

    def test_each_element_has_shape_n_saves(self, stub_run):
        """Each array has shape (n_saves,)."""
        result = find_front(stub_run, FRONT_THRESHOLD)
        n_saves = len(stub_run.data.macroscale_out[0].snapshot_time[:])
        for arr in result:
            assert arr.shape == (n_saves,)

    def test_values_within_valid_range(self, stub_run):
        """Front indices lie in [0, fiber_rows - 1]."""
        fiber_rows = stub_run.macro_params.fiber_rows
        result = find_front(stub_run, FRONT_THRESHOLD)
        for arr in result:
            assert np.all(arr >= 0)
            assert np.all(arr <= fiber_rows - 1)

    def test_front_at_row_zero_when_all_intact(self, stub_run):
        """At save-point 0 all rows are intact, so the front is at row 0."""
        result = find_front(stub_run, FRONT_THRESHOLD)
        for arr in result:
            assert arr[0] == 0

    def test_front_at_last_row_when_no_row_meets_threshold(self, stub_run):
        """When no row meets the threshold the front defaults to fiber_rows - 1."""
        fiber_rows = stub_run.macro_params.fiber_rows
        # At save-point 4 all fibers are degraded; no row meets threshold=0.5
        result = find_front(stub_run, FRONT_THRESHOLD)
        for arr in result:
            assert arr[-1] == fiber_rows - 1


# ===========================================================================
# fiber_degradation_linear_extrapolation
# ===========================================================================


class TestFiberDegradationLinearExtrapolation:
    """Tests for :func:`fiber_degradation_linear_extrapolation`."""

    def test_returns_array_of_correct_shape(self, stub_run):
        """Returns array of shape (n_sims, max_t, fiber_rows-1, full_row)."""
        p = stub_run.macro_params
        result = fiber_degradation_linear_extrapolation(stub_run)
        n_sims = p.macro_simulations
        max_t = min(
            len(stub_run.data.macroscale_out[sim].snapshot_time[:])
            for sim in range(n_sims)
        )
        assert result.shape == (
            p.macro_simulations,
            max_t,
            p.fiber_rows - 1,
            p.full_row,
        )

    def test_values_in_zero_to_one(self, stub_run):
        """All remaining-fibrin fractions lie in [0, 1]."""
        result = fiber_degradation_linear_extrapolation(stub_run)
        assert np.all(result >= -1e-9)
        assert np.all(result <= 1.0 + 1e-9)

    def test_initial_fiber_amount_is_one(self, stub_run):
        """At save-point 0 every fiber is fully intact (remaining fraction = 1)."""
        result = fiber_degradation_linear_extrapolation(stub_run)
        np.testing.assert_allclose(result[:, 0, :, :], 1.0)

    def test_final_fiber_amount_is_zero(self, stub_run):
        """At the last save-point all fibers have degraded (remaining fraction = 0)."""
        result = fiber_degradation_linear_extrapolation(stub_run)
        np.testing.assert_allclose(result[:, -1, :, :], 0.0, atol=1e-9)


# ===========================================================================
# plot_degradation_percent
# ===========================================================================


class TestPlotDegradationPercent:
    """Tests for :func:`plot_degradation_percent`."""

    def test_returns_figure(self, stub_run):
        """Return value is a matplotlib Figure."""
        fig = plot_degradation_percent(stub_run, SLOPE_PAIRS, PERCENT_MARKERS)
        assert isinstance(fig, Figure)
        plt.close(fig)

    def test_figure_has_axes(self, stub_run):
        """The returned Figure contains at least one Axes."""
        fig = plot_degradation_percent(stub_run, SLOPE_PAIRS, PERCENT_MARKERS)
        assert len(fig.get_axes()) > 0
        plt.close(fig)

    def test_axes_has_lines(self, stub_run):
        """The axes contains one degradation curve per simulation."""
        n_sims = stub_run.macro_params.macro_simulations
        fig = plot_degradation_percent(stub_run, SLOPE_PAIRS, PERCENT_MARKERS)
        ax = fig.get_axes()[0]
        # One curve per sim plus one linear-fit segment per sim per slope pair
        n_expected_lines = n_sims + n_sims * len(SLOPE_PAIRS)
        assert len(ax.get_lines()) == n_expected_lines
        plt.close(fig)


# ===========================================================================
# plot_front_degradation
# ===========================================================================


class TestPlotFrontDegradation:
    """Tests for :func:`plot_front_degradation`."""

    def test_returns_figure(self, stub_run):
        """Return value is a matplotlib Figure."""
        fig = plot_front_degradation(stub_run)
        assert isinstance(fig, Figure)
        plt.close(fig)

    def test_figure_has_axes(self, stub_run):
        """The returned Figure contains at least one Axes."""
        fig = plot_front_degradation(stub_run)
        assert len(fig.get_axes()) > 0
        plt.close(fig)

    def test_axes_has_lines(self, stub_run):
        """The axes contains one line per simulation per column."""
        n_sims = stub_run.macro_params.macro_simulations
        cols = stub_run.macro_params.cols
        fig = plot_front_degradation(stub_run)
        ax = fig.get_axes()[0]
        assert len(ax.get_lines()) == n_sims * cols
        plt.close(fig)


# ===========================================================================
# Log file parsers
# ===========================================================================


class _LogRunStub:
    """Stub Run with os_path pointing at a temp directory."""

    def __init__(self, macro_params, os_path):
        self.macro_params = macro_params
        self.os_path = os_path
        self._cache = {}


def _write_log_file(directory, filename, content):
    path = os.path.join(directory, filename)
    with open(path, "w") as fh:
        fh.write(content)


LOG_CONTENT_TEMPLATE = """\
Some preamble text.
countmacrounbd= {macro}
countmicrounbd= {micro}
Processing time: {proc_time} sec
Total Binds: {binds}
Some trailing text.
"""


class TestGetUnbindAmounts:
    """Tests for :func:`get_unbind_amounts`."""

    @pytest.fixture
    def log_run(self, tmp_path, small_macro):
        """Stub run with fake log files in tmp_path/00 and tmp_path/01."""
        for sim, (macro_u, micro_u) in enumerate([(12, 3), (45, 6)]):
            sim_dir = tmp_path / f"{sim:02}"
            sim_dir.mkdir()
            content = LOG_CONTENT_TEMPLATE.format(
                macro=macro_u, micro=micro_u,
                proc_time="100.0", binds="999",
            )
            _write_log_file(sim_dir, f"macro_test_{sim:02}.txt", content)
        return _LogRunStub(small_macro, str(tmp_path))

    def test_returns_two_arrays(self, log_run):
        """Returns a 2-tuple of numpy arrays."""
        result = get_unbind_amounts(log_run, "_test")
        assert isinstance(result, tuple) and len(result) == 2

    def test_correct_macro_unbind_counts(self, log_run):
        """Correctly reads macroscale forced-unbind counts."""
        macro_unbinds, _ = get_unbind_amounts(log_run, "_test")
        np.testing.assert_array_equal(macro_unbinds, [12, 45])

    def test_correct_micro_unbind_counts(self, log_run):
        """Correctly reads microscale forced-unbind counts."""
        _, micro_unbinds = get_unbind_amounts(log_run, "_test")
        np.testing.assert_array_equal(micro_unbinds, [3, 6])


class TestGetProcessingTime:
    """Tests for :func:`get_processing_time`."""

    @pytest.fixture
    def log_run(self, tmp_path, small_macro):
        """Stub run with fake log files containing processing times."""
        for sim, proc in enumerate([123.45, 678.90]):
            sim_dir = tmp_path / f"{sim:02}"
            sim_dir.mkdir()
            content = LOG_CONTENT_TEMPLATE.format(
                macro=0, micro=0, proc_time=f"{proc:.2f}", binds="0",
            )
            _write_log_file(sim_dir, f"macro_test_{sim:02}.txt", content)
        return _LogRunStub(small_macro, str(tmp_path))

    def test_returns_float_array(self, log_run):
        """Returns a 1-D array of floats."""
        result = get_processing_time(log_run, "_test")
        assert result.dtype.kind == "f"

    def test_correct_processing_times(self, log_run):
        """Correctly reads processing times in seconds."""
        result = get_processing_time(log_run, "_test")
        np.testing.assert_allclose(result, [123.45, 678.90])

    def test_array_length_equals_n_sims(self, log_run, small_macro):
        """Array length equals macro_simulations."""
        result = get_processing_time(log_run, "_test")
        assert len(result) == small_macro.macro_simulations


class TestGetTotalBinds:
    """Tests for :func:`get_total_binds`."""

    @pytest.fixture
    def log_run(self, tmp_path, small_macro):
        """Stub run with fake log files containing total-bind counts."""
        for sim, binds in enumerate([1234, 5678]):
            sim_dir = tmp_path / f"{sim:02}"
            sim_dir.mkdir()
            content = LOG_CONTENT_TEMPLATE.format(
                macro=0, micro=0, proc_time="1.0", binds=binds,
            )
            _write_log_file(sim_dir, f"macro_test_{sim:02}.txt", content)
        return _LogRunStub(small_macro, str(tmp_path))

    def test_returns_integer_array(self, log_run):
        """Returns a 1-D array of integers."""
        result = get_total_binds(log_run, "_test")
        assert result.dtype.kind == "i" or result.dtype.kind == "u"

    def test_correct_bind_counts(self, log_run):
        """Correctly reads total bind counts."""
        result = get_total_binds(log_run, "_test")
        np.testing.assert_array_equal(result, [1234, 5678])

    def test_array_length_equals_n_sims(self, log_run, small_macro):
        """Array length equals macro_simulations."""
        result = get_total_binds(log_run, "_test")
        assert len(result) == small_macro.macro_simulations


# ===========================================================================
# _build_percent_markers
# ===========================================================================


class TestBuildPercentMarkers:
    """Tests for :func:`~lysis.analysis.degradation._build_percent_markers`."""

    def test_includes_bookends(self):
        from lysis.analysis.degradation import _build_percent_markers

        markers = _build_percent_markers([(20, 80)])
        assert 0.0 in markers
        assert 1.0 in markers

    def test_includes_endpoints(self):
        from lysis.analysis.degradation import _build_percent_markers

        markers = _build_percent_markers([(20, 80)])
        assert 0.20 in markers
        assert 0.80 in markers

    def test_sorted(self):
        from lysis.analysis.degradation import _build_percent_markers

        markers = _build_percent_markers([(20, 80), (50, 90)])
        assert markers == sorted(markers)

    def test_no_duplicates(self):
        from lysis.analysis.degradation import _build_percent_markers

        markers = _build_percent_markers([(20, 80), (20, 50), (50, 80)])
        assert len(markers) == len(set(markers))
