"""Analysis functions for fibrinolysis simulation Runs.

Provides functions for quantifying and visualising clot-degradation dynamics
from macroscale simulation output.  Functions that need fiber state use
:class:`~lysis.analysis.fiber_replay.FiberReplayCursor` internally and read
save-point times from ``run.data.macroscale_out[sim].snapshot_time``.

Typical workflow::

    from lysis.analysis.degradation import (
        find_degraded_fraction,
        find_degradation_marker_frames,
        find_degradation_marker_times,
        mean_degradation_rate,
        calculate_time_row_exposed,
        find_degradation_fronts,
        mean_front_velocity,
        plot_degradation_percent,
        plot_front_degradation,
    )

    percent_markers = [0.0, 0.25, 0.5, 0.75, 1.0]
    slope_pairs     = [(0.25, 0.75)]

    deg_fraction  = find_degraded_fraction(run)
    marker_frames = find_degradation_marker_frames(run, percent_markers)
    marker_times  = find_degradation_marker_times(run, percent_markers)
    rates         = degradation_rates(run, slope_pairs, percent_markers)
    deg_rate, offset, lag = mean_degradation_rate(run)
    fig = plot_degradation_percent(run, slope_pairs, percent_markers)
"""

__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = "GPLv3"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"

import os
import re
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure

if TYPE_CHECKING:
    from lysis.config.run import Run

from lysis.analysis.fiber_replay import FiberReplayCursor

__all__ = [
    "find_degraded_fraction",
    "find_degradation_marker_frames",
    "find_degradation_marker_times",
    "degradation_rates",
    "mean_degradation_rate",
    "calculate_time_row_exposed",
    "find_degradation_fronts",
    "per_sim_front_velocity",
    "mean_front_velocity",
    "find_row_deg_fraction",
    "find_front",
    "fiber_degradation_linear_extrapolation",
    "plot_degradation_percent",
    "plot_front_degradation",
    "get_unbind_amounts",
    "get_processing_time",
    "get_total_binds",
    "compute_run_statistics",
    "compute_degradation_marker_stats",
    "compute_degradation_rate_stats",
]


###############################################################################
# Cache keys
###############################################################################

_KEY_DEGRADED_FRACTION = "degraded_fraction"
_KEY_EXPOSED_TIME = "exposed_time"
_KEY_ROW_DEG = "row_deg"
_KEY_DEG_FRONTS = "deg_fronts"
_KEY_MEAN_DEG_RATE = "mean_deg_rate"
_KEY_MEAN_FRONT_VEL = "mean_front_vel"
_KEY_FRONT_VEL_PER_SIM = "front_vel_per_sim"
_KEY_FIBER_EXTRAP = "fiber_extrap"
_KEY_RUN_STATS = "run_stats"


###############################################################################
# Degradation fraction
###############################################################################


def find_degraded_fraction(
    run: "Run",
) -> list[np.ndarray]:
    """Calculate the fraction of fibrin fibers degraded at each save point.

    For each simulation, uses a :class:`FiberReplayCursor` to replay the
    event log and counts edges whose degrade time has passed at each save
    point, subtracts the fibrin-free (empty) edges, and normalises by the
    total number of fibrin fibers.

    :param run: Run object supplying macro grid parameters and data.
    :type run: Run
    :return: Per-simulation arrays of degraded fraction in ``[0, 1]``.
        ``degraded_fraction[sim]`` has shape ``(n_save,)``.
    :rtype: list[numpy.ndarray]
    """
    if _KEY_DEGRADED_FRACTION not in run._cache:
        n_sims = run.macro_params.macro_simulations
        empty_edges = run.macro_params.empty_edges
        total_fibers = run.macro_params.total_fibers

        degraded_fraction = []
        for sim in range(n_sims):
            tsave = run.data.macroscale_out[sim].snapshot_time[:]
            cursor = FiberReplayCursor(run, sim)
            run_frac = np.empty(tsave.shape[0], dtype=np.float64)
            for t_idx, t in enumerate(tsave):
                cursor.advance_to(t)
                run_frac[t_idx] = np.count_nonzero(cursor.state <= t)
            run_frac -= empty_edges
            degraded_fraction.append(run_frac / total_fibers)
        run._cache[_KEY_DEGRADED_FRACTION] = degraded_fraction
    return run._cache[_KEY_DEGRADED_FRACTION]


###############################################################################
# Degradation markers
###############################################################################


def find_degradation_marker_frames(
    run: "Run",
    percent_markers: list[float],
) -> np.ndarray:
    """Find the save-point frame index at which each degradation milestone is first reached.

    For each simulation and each milestone fraction, returns the index of the
    first save point at which the degraded fraction is greater than or equal to
    that milestone.  Calls :func:`find_degraded_fraction` internally.

    :param run: Run object supplying macro grid parameters and data.
    :type run: Run
    :param percent_markers: Degradation milestones to locate, e.g.
        ``[0.0, 0.25, 0.50, 0.75, 1.0]``.
    :type percent_markers: list[float]
    :return: Integer array of shape ``(n_sims, n_markers)``.
        ``marker_frames[sim, m]`` is the save-point index at which simulation
        ``sim`` first reaches milestone ``percent_markers[m]``.
    :rtype: numpy.ndarray
    """
    key = ("marker_frames", tuple(percent_markers))
    if key not in run._cache:
        degraded_fraction = find_degraded_fraction(run)
        n_sims = len(degraded_fraction)
        n_markers = len(percent_markers)
        marker_frames = np.empty((n_sims, n_markers), dtype=np.intp)
        for sim in range(n_sims):
            for m, threshold in enumerate(percent_markers):
                marker_frames[sim, m] = np.argmax(degraded_fraction[sim] >= threshold)
        run._cache[key] = marker_frames
    return run._cache[key]


def find_degradation_marker_times(
    run: "Run",
    percent_markers: list[float],
) -> np.ndarray:
    """Convert degradation milestones to elapsed time in minutes.

    Calls :func:`find_degradation_marker_frames` (which in turn calls
    :func:`find_degraded_fraction`) to locate the save-point frame for each
    milestone, then maps those frame indices to wall-clock times.

    :param run: Run object supplying data access.
    :type run: Run
    :param percent_markers: Degradation milestones to locate, e.g.
        ``[0.0, 0.25, 0.50, 0.75, 1.0]``.
    :type percent_markers: list[float]
    :return: Float array of shape ``(n_sims, n_markers)`` giving the elapsed
        time (minutes) at which each simulation reached each degradation
        milestone.
    :rtype: numpy.ndarray
    """
    key = ("marker_times", tuple(percent_markers))
    if key not in run._cache:
        marker_frames = find_degradation_marker_frames(run, percent_markers)
        n_sims, n_markers = marker_frames.shape
        marker_times = np.empty((n_sims, n_markers), dtype=np.float64)
        for sim in range(n_sims):
            tsave = run.data.macroscale_out[sim].snapshot_time[:]
            marker_times[sim] = tsave[marker_frames[sim]]
        run._cache[key] = marker_times / 60
    return run._cache[key]


###############################################################################
# Degradation rates
###############################################################################


def degradation_rates(
    run: "Run",
    slope_pairs: list[tuple[float, float]],
    percent_markers: list[float],
) -> np.ndarray:
    """Calculate degradation rates (fraction/min) between pairs of milestones.

    For each simulation and each pair of milestones, computes the finite-
    difference slope between the two marker frames and returns the value in
    units of fraction per minute.  Calls :func:`find_degraded_fraction` and
    :func:`find_degradation_marker_frames` internally.

    :param run: Run object supplying ``macro_simulations`` count and data.
    :type run: Run
    :param slope_pairs: List of ``(start_percent, end_percent)`` pairs defining
        the intervals over which to compute degradation rates.
        Example: ``[(0.25, 0.75)]`` computes one slope from 25 % to 75 %.
    :type slope_pairs: list[tuple[float, float]]
    :param percent_markers: The full list of milestone fractions used when
        computing marker frames.  Used to look up the frame indices
        corresponding to the endpoints of each slope pair.
    :type percent_markers: list[float]
    :return: Array of shape ``(n_sims, n_slopes)`` containing degradation rates
        in fraction/minute.  ``rates[sim, k]`` is the slope for simulation
        ``sim`` between the milestones defined by ``slope_pairs[k]``.
    :rtype: numpy.ndarray
    """
    key = ("rates", tuple(tuple(s) for s in slope_pairs), tuple(percent_markers))
    if key not in run._cache:
        degraded_fraction = find_degraded_fraction(run)
        marker_frames = find_degradation_marker_frames(run, percent_markers)
        n_sims = run.macro_params.macro_simulations
        rates = np.empty((n_sims, len(slope_pairs)), dtype=np.float64)
        for sim in range(n_sims):
            tsave = run.data.macroscale_out[sim].snapshot_time[:]
            for k, (start_pct, end_pct) in enumerate(slope_pairs):
                start_frame = marker_frames[sim, percent_markers.index(start_pct)]
                end_frame = marker_frames[sim, percent_markers.index(end_pct)]
                delta_frac = (
                    degraded_fraction[sim][end_frame] - degraded_fraction[sim][start_frame]
                )
                delta_time = tsave[end_frame] - tsave[start_frame]
                rates[sim, k] = delta_frac / delta_time * 60  # fraction/min
        run._cache[key] = rates
    return run._cache[key]


def mean_degradation_rate(
    run: "Run",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Estimate the degradation rate by fitting a line through the rapid-degradation phase.

    Identifies the rapid-degradation phase as the set of save points where the
    incremental change in degraded fraction exceeds half the maximum incremental
    change observed in that simulation.  Fits a degree-1 polynomial to those
    points.  Calls :func:`find_degraded_fraction` internally.

    :param run: Run object supplying ``macro_simulations`` count and data.
    :type run: Run
    :return: Tuple of three 1-D arrays of length ``n_sims``:

        - ``degradation_rate``: Slope of the linear fit (fraction/min).
        - ``offset``: y-intercept of the linear fit (fraction).
        - ``deg_start_time``: Time (min) of the first save point in the
          rapid-degradation phase.

    :rtype: tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]
    """
    if _KEY_MEAN_DEG_RATE not in run._cache:
        degraded_fraction = find_degraded_fraction(run)
        n_sims = run.macro_params.macro_simulations
        degradation_rate = np.empty(n_sims, dtype=np.float64)
        offset = np.empty(n_sims, dtype=np.float64)
        deg_start_time = np.empty(n_sims, dtype=np.float64)

        for r in range(n_sims):
            tsave = run.data.macroscale_out[r].snapshot_time[:]
            incremental = np.empty(degraded_fraction[r].shape[0], dtype=np.float64)
            incremental[0] = degraded_fraction[r][0]
            for t in range(1, degraded_fraction[r].shape[0]):
                incremental[t] = degraded_fraction[r][t] - degraded_fraction[r][t - 1]
            rapid_phase = incremental > incremental.max() / 2
            s = np.argmax(rapid_phase)
            b, m = np.polynomial.polynomial.polyfit(
                tsave[rapid_phase] / 60, degraded_fraction[r][rapid_phase], 1
            )
            degradation_rate[r] = m
            offset[r] = b
            deg_start_time[r] = tsave[s] / 60

        run._cache[_KEY_MEAN_DEG_RATE] = (degradation_rate, offset, deg_start_time)
    return run._cache[_KEY_MEAN_DEG_RATE]


###############################################################################
# Lysis front
###############################################################################


def calculate_time_row_exposed(
    run: "Run",
) -> np.ndarray:
    """Compute the time (min) at which each fiber row in each column first becomes exposed.

    A row is considered "exposed" once all fibers in the row immediately above it
    (closer to the tPA source) have degraded.  Row 0 is always exposed at time 0.
    For subsequent rows, the exposure time is the maximum of the exposure time of
    the preceding row and the final scheduled degradation time of the edge at
    position ``(i, j)``.

    Uses a :class:`FiberReplayCursor` advanced to the final save point to read
    fiber degrade times via 2-D ``(row, rank)`` indexing.

    :param run: Run object supplying grid parameters and data.
    :type run: Run
    :return: Array of shape ``(n_sims, rows - 1, cols)`` giving exposure time
        in minutes.  ``exposed_time[sim, i, j]`` is the time (min) at which
        row ``i`` in column ``j`` first became fully exposed in simulation
        ``sim``.
    :rtype: numpy.ndarray
    """
    if _KEY_EXPOSED_TIME not in run._cache:
        n_sims = run.macro_params.macro_simulations
        rows = run.macro_params.rows
        cols = run.macro_params.cols

        exposed_time = np.empty((n_sims, rows - 1, cols), dtype=np.float64)
        for sim in range(n_sims):
            tsave = run.data.macroscale_out[sim].snapshot_time[:]
            cursor = FiberReplayCursor(run, sim)
            cursor.advance_to(tsave[-1])
            state = cursor.state
            for j in range(cols):
                for i in range(rows - 1):
                    if i == 0:
                        exposed_time[sim, i, j] = 0
                    else:
                        exposed_time[sim, i, j] = max(
                            exposed_time[sim, i - 1, j], state[i, j]
                        )
        run._cache[_KEY_EXPOSED_TIME] = exposed_time / 60
    return run._cache[_KEY_EXPOSED_TIME]


def find_degradation_fronts(
    run: "Run",
) -> list[list[np.ndarray]]:
    """Track the spatial position of the lysis front over time in each column.

    For each column, collects ``(time, y-distance)`` pairs at which each
    successive fiber row first becomes exposed.  Rows that are never exposed
    before the end of the simulation are excluded.  y-distances are computed
    directly from the run's pore size.  Calls
    :func:`calculate_time_row_exposed` internally.

    :param run: Run object supplying grid parameters, pore size, and data.
    :type run: Run
    :return: Nested list ``deg_fronts[sim][col]`` where each element is an
        ``ndarray`` of shape ``(2, n_events)`` with row 0 giving event times
        (min) and row 1 giving y-distances (µm).
    :rtype: list[list[numpy.ndarray]]
    """
    if _KEY_DEG_FRONTS not in run._cache:
        exposed_time = calculate_time_row_exposed(run)
        n_sims = run.macro_params.macro_simulations
        rows = run.macro_params.rows
        cols = run.macro_params.cols
        grid_node_distance = (
            run.macro_params.pore_size + 2 * run.macro_params.micro_params.fiber_radius
        )
        y_distance = (np.arange(rows - 1) * grid_node_distance).to("microns").magnitude

        deg_fronts = []
        for sim in range(n_sims):
            tsave = run.data.macroscale_out[sim].snapshot_time[:]
            t_end_min = tsave[-1] / 60
            sim_fronts = []
            for j in range(cols):
                col_front = []
                for i in range(1, rows - 1):
                    if (
                        exposed_time[sim, i - 1, j]
                        < exposed_time[sim, i, j]
                        < t_end_min + 1
                    ):
                        col_front.append([exposed_time[sim, i, j], y_distance[i]])
                sim_fronts.append(np.array(col_front).T)
            deg_fronts.append(sim_fronts)
        run._cache[_KEY_DEG_FRONTS] = deg_fronts
    return run._cache[_KEY_DEG_FRONTS]


def per_sim_front_velocity(
    run: "Run",
) -> tuple[np.ndarray, np.ndarray]:
    """Compute per-simulation lysis-front velocity summaries (µm/min).

    For each simulation, fits a line to the ``(time, y-distance)`` data of
    every column and takes the slope as that column's front velocity.
    Returns the per-simulation mean and standard deviation across columns
    as two length-``n_sims`` arrays.  Calls :func:`find_degradation_fronts`
    internally.

    :param run: Run object supplying ``macro_simulations`` and ``cols`` counts.
    :type run: Run
    :return: Tuple ``(run_mean, run_std)`` of 1-D arrays of length
        ``n_sims``, each in µm/min.
    :rtype: tuple[numpy.ndarray, numpy.ndarray]
    """
    if _KEY_FRONT_VEL_PER_SIM not in run._cache:
        deg_fronts = find_degradation_fronts(run)
        n_sims = run.macro_params.macro_simulations
        cols = run.macro_params.cols

        run_mean = np.empty(n_sims, dtype=np.float64)
        run_std = np.empty(n_sims, dtype=np.float64)
        for sim in range(n_sims):
            front_velocity = np.empty(cols, dtype=np.float64)
            for j in range(cols):
                _b, m = np.polynomial.polynomial.polyfit(
                    deg_fronts[sim][j][0], deg_fronts[sim][j][1], 1
                )
                front_velocity[j] = m
            run_mean[sim] = np.mean(front_velocity)
            run_std[sim] = np.std(front_velocity)
        run._cache[_KEY_FRONT_VEL_PER_SIM] = (run_mean, run_std)
    return run._cache[_KEY_FRONT_VEL_PER_SIM]


def mean_front_velocity(
    run: "Run",
) -> tuple[float, float]:
    """Compute the mean and standard deviation of lysis-front velocity (µm/min).

    For each column in each simulation, fits a line to the ``(time, y-distance)``
    data and takes the slope as the column's front velocity.  Returns the grand
    mean and the mean standard deviation across all simulations.  Delegates to
    :func:`per_sim_front_velocity` for the per-simulation arrays.

    .. todo::

        Change to compute a single mean and std over all columns across all
        simulations rather than averaging per-simulation statistics.

    :param run: Run object supplying ``macro_simulations`` and ``cols`` counts.
    :type run: Run
    :return: Tuple ``(mean_velocity, std_velocity)`` in µm/min.
    :rtype: tuple[float, float]
    """
    if _KEY_MEAN_FRONT_VEL not in run._cache:
        run_mean, run_std = per_sim_front_velocity(run)
        run._cache[_KEY_MEAN_FRONT_VEL] = (
            float(np.mean(run_mean)), float(np.mean(run_std))
        )
    return run._cache[_KEY_MEAN_FRONT_VEL]


###############################################################################
# Per-row degradation and fiber-level interpolation
###############################################################################


def find_row_deg_fraction(
    run: "Run",
) -> list[np.ndarray]:
    """Compute the fraction of fibers NOT yet degraded in each row at each save point.

    Uses a :class:`FiberReplayCursor` to replay events at each save point and
    slices the fibrous region as ``cursor.state[empty_rows:rows-1, :]`` to
    count undegraded edges per row.

    :param run: Run object supplying grid parameters and data.
    :type run: Run
    :return: Per-simulation arrays of shape
        ``(n_save, fiber_rows - 1)`` containing the fraction of
        edges in each row that have NOT yet degraded at each save point.
    :rtype: list[numpy.ndarray]
    """
    if _KEY_ROW_DEG not in run._cache:
        n_sims = run.macro_params.macro_simulations
        empty_rows = run.macro_params.empty_rows
        rows = run.macro_params.rows
        full_row = run.macro_params.full_row

        row_deg = []
        for sim in range(n_sims):
            tsave = run.data.macroscale_out[sim].snapshot_time[:]
            cursor = FiberReplayCursor(run, sim)
            n_saves = tsave.shape[0]
            fiber_rows_minus_1 = rows - 1 - empty_rows
            sim_row_deg = np.empty((n_saves, fiber_rows_minus_1), dtype=np.float64)
            for t_idx, t in enumerate(tsave):
                cursor.advance_to(t)
                fibrous = cursor.state[empty_rows : rows - 1, :]
                sim_row_deg[t_idx] = np.count_nonzero(fibrous > t, axis=1) / full_row
            row_deg.append(sim_row_deg)
        run._cache[_KEY_ROW_DEG] = row_deg
    return run._cache[_KEY_ROW_DEG]


def find_front(
    run: "Run",
    front_threshold: float,
) -> list[np.ndarray]:
    """Locate the lysis-front row index from per-row undegraded fractions.

    For each simulation and save point, finds the first row whose undegraded
    fraction meets or exceeds ``front_threshold``.  If no row meets the
    threshold, the front is set to ``fiber_rows - 1`` (the last row).  Calls
    :func:`find_row_deg_fraction` internally.

    :param run: Run object supplying ``fiber_rows`` count.
    :type run: Run
    :param front_threshold: Minimum undegraded fraction for a row to be
        counted as part of the lysis front.
    :type front_threshold: float
    :return: Per-simulation integer arrays of shape ``(n_save,)`` giving the
        row index of the lysis front at each save point.
    :rtype: list[numpy.ndarray]
    """
    key = ("front", front_threshold)
    if key not in run._cache:
        row_deg = find_row_deg_fraction(run)
        fiber_rows = run.macro_params.fiber_rows
        fronts = []
        for sim_row_deg in row_deg:
            sim_fronts = np.argmax(sim_row_deg >= front_threshold, axis=1)
            sim_fronts[np.max(sim_row_deg >= front_threshold, axis=1) == 0] = fiber_rows - 1
            fronts.append(sim_fronts)
        run._cache[key] = fronts
    return run._cache[key]


def fiber_degradation_linear_extrapolation(
    run: "Run",
) -> np.ndarray:
    """Estimate fiber-level remaining-fibrin fractions by linear extrapolation.

    At each save point, uses the rate of change between successive save points
    to linearly extrapolate the fractional remaining fibrin for every edge.
    The result is stacked across simulations and clipped to the fibrin-
    containing region.

    Uses a :class:`FiberReplayCursor` to replay the event log at each save
    point, building 3-D ``(n_saves, rows, full_row)`` degrade-time snapshots
    with NaN padding for the last row.

    :param run: Run object supplying grid parameters and data.
    :type run: Run
    :return: Array of shape ``(n_sims, max_t, fiber_rows - 1, full_row)``
        containing the estimated remaining-fibrin fraction for each fiber at
        each save point.  ``max_t`` is the minimum number of save points across
        all simulations.
    :rtype: numpy.ndarray
    """
    if _KEY_FIBER_EXTRAP not in run._cache:
        n_sims = run.macro_params.macro_simulations
        empty_rows = run.macro_params.empty_rows
        rows = run.macro_params.rows
        fiber_rows = run.macro_params.fiber_rows
        full_row = run.macro_params.full_row

        f_deg_amt = []
        all_tsave = []
        for sim in range(n_sims):
            tsave = run.data.macroscale_out[sim].snapshot_time[:]
            all_tsave.append(tsave)
            n_saves = tsave.shape[0]
            cursor = FiberReplayCursor(run, sim)

            # Build snapshots of cursor state at each save point
            snapshots = np.empty((n_saves, rows, full_row), dtype=np.float64)
            for t_idx, t in enumerate(tsave):
                cursor.advance_to(t)
                snapshots[t_idx] = cursor.state

            # Build remaining-fibrin fraction array for the fibrous region
            # Shape: (n_saves, fiber_rows - 1, full_row)
            fibrous_rows = rows - 1 - empty_rows  # = fiber_rows - 1
            f_deg_amt_sim = np.ones(
                (n_saves, fibrous_rows, full_row), dtype=np.float64
            )
            f_deg_delta = np.zeros((fibrous_rows, full_row), dtype=np.float64)

            for i in range(1, n_saves):
                prev = snapshots[i - 1, empty_rows : rows - 1, :]
                curr = snapshots[i, empty_rows : rows - 1, :]
                changed = np.argwhere(curr < prev)
                if changed.size > 0:
                    r_idx, c_idx = changed[:, 0], changed[:, 1]
                    f_deg_delta[r_idx, c_idx] = f_deg_amt_sim[i - 1, r_idx, c_idx] / (
                        curr[r_idx, c_idx] - tsave[i - 1]
                    )
                f_deg_amt_sim[i] = np.maximum(
                    f_deg_amt_sim[i - 1]
                    - f_deg_delta * (tsave[i] - tsave[i - 1]),
                    0.0,
                )
            f_deg_amt.append(f_deg_amt_sim)

        max_t = min(len(t) for t in all_tsave)
        run._cache[_KEY_FIBER_EXTRAP] = np.stack(
            [f_deg_amt[sim][:max_t] for sim in range(n_sims)]
        )
    return run._cache[_KEY_FIBER_EXTRAP]


###############################################################################
# Plots
###############################################################################


def plot_degradation_percent(
    run: "Run",
    slope_pairs: list[tuple[float, float]],
    percent_markers: list[float],
) -> Figure:
    """Plot degradation curves with overlaid linear-fit segments for each simulation.

    For each simulation, draws the degraded fraction over time and overlays a
    linear-fit segment (in blue) for each slope pair.  Calls
    :func:`find_degraded_fraction`, :func:`find_degradation_marker_frames`, and
    :func:`degradation_rates` internally.

    :param run: Run object supplying ``macro_simulations`` count and data.
    :type run: Run
    :param slope_pairs: List of ``(start_percent, end_percent)`` pairs defining
        the intervals over which rates were computed.
    :type slope_pairs: list[tuple[float, float]]
    :param percent_markers: The full list of milestone fractions used when
        computing marker frames.
    :type percent_markers: list[float]
    :return: The matplotlib Figure containing the plot.
    :rtype: matplotlib.figure.Figure
    """
    degraded_fraction = find_degraded_fraction(run)
    marker_frames = find_degradation_marker_frames(run, percent_markers)
    rates = degradation_rates(run, slope_pairs, percent_markers)
    n_sims = run.macro_params.macro_simulations
    start_stop = [
        [percent_markers.index(mark) for mark in slope] for slope in slope_pairs
    ]

    all_tsave = [
        run.data.macroscale_out[sim].snapshot_time[:] for sim in range(n_sims)
    ]
    x_max = (max(all_tsave[sim][-1] for sim in range(n_sims)) // 60) + 1

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, x_max)
    ax.set_ylim(-0.1, 1.1)

    for sim in range(n_sims):
        tsave = all_tsave[sim]
        ax.plot(tsave / 60, degraded_fraction[sim])
        for k, (start_idx, end_idx) in enumerate(start_stop):
            start_frame = marker_frames[sim, start_idx]
            end_frame = marker_frames[sim, end_idx]
            x_seg = tsave[start_frame:end_frame] / 60
            y_seg = (
                rates[sim, k] * (x_seg - tsave[start_frame] / 60)
                + degraded_fraction[sim][start_frame]
            )
            ax.plot(x_seg, y_seg, color="b", alpha=0.5, zorder=0.1)

    return fig


def plot_front_degradation(
    run: "Run",
) -> Figure:
    """Plot the spatial lysis-front trajectories for all simulations.

    For each simulation and column, draws a line tracing the time (x-axis,
    minutes) at which the lysis front reached each successive y-position
    (y-axis, µm).  Calls :func:`find_degradation_fronts` internally.

    :param run: Run object supplying grid parameters, pore size, and data.
    :type run: Run
    :return: The matplotlib Figure containing the plot.
    :rtype: matplotlib.figure.Figure
    """
    deg_fronts = find_degradation_fronts(run)
    n_sims = run.macro_params.macro_simulations
    cols = run.macro_params.cols
    empty_rows = run.macro_params.empty_rows
    rows = run.macro_params.rows
    grid_node_distance = (
        run.macro_params.pore_size + 2 * run.macro_params.micro_params.fiber_radius
    )

    x_max = (
        max(
            run.data.macroscale_out[sim].snapshot_time[:][-1]
            for sim in range(n_sims)
        )
        // 60
    ) + 1

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_on()
    ax.set_xlim(0, x_max)
    ax.set_ylim(
        ((empty_rows - 1) * grid_node_distance).to("microns").magnitude,
        ((rows - 1) * grid_node_distance).to("microns").magnitude,
    )

    for sim in range(n_sims):
        for j in range(cols):
            ax.plot(deg_fronts[sim][j][0], deg_fronts[sim][j][1], linewidth=1)

    return fig


###############################################################################
# Log file parsers
###############################################################################


def get_unbind_amounts(
    run: "Run",
    file_code: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Read forced-unbind counts from macroscale log files.

    Parses ``macro{file_code}_{sim:02}.txt`` log files for each simulation
    and extracts the macroscale and microscale forced-unbind counts.

    :param run: Run object supplying ``os_path`` and ``macro_simulations``.
    :type run: Run
    :param file_code: File-code string used to construct log file names.
    :type file_code: str
    :return: Tuple ``(macro_unbinds, micro_unbinds)`` where each element is a
        1-D integer array of length ``n_sims``.
    :rtype: tuple[numpy.ndarray, numpy.ndarray]
    """
    key = ("unbind_amounts", file_code)
    if key not in run._cache:
        macro_pattern = re.compile(r"countmacrounbd=\s*(\d+)")
        micro_pattern = re.compile(r"countmicrounbd=\s*(\d+)")
        log_text = ""
        for sim in range(run.macro_params.macro_simulations):
            log_file = os.path.join(run.os_path, f"{sim:02}", f"macro{file_code}_{sim:02}.txt")
            with open(log_file) as fh:
                log_text += fh.read()
        macro_unbinds = np.array(re.findall(macro_pattern, log_text), dtype=int)
        micro_unbinds = np.array(re.findall(micro_pattern, log_text), dtype=int)
        run._cache[key] = (macro_unbinds, micro_unbinds)
    return run._cache[key]


def get_processing_time(
    run: "Run",
    file_code: str,
) -> np.ndarray:
    """Read simulation processing times (seconds) from macroscale log files.

    :param run: Run object supplying ``os_path`` and ``macro_simulations``.
    :type run: Run
    :param file_code: File-code string used to construct log file names.
    :type file_code: str
    :return: 1-D float array of processing times in seconds, one per simulation.
    :rtype: numpy.ndarray
    """
    key = ("processing_time", file_code)
    if key not in run._cache:
        pattern = re.compile(r"Processing time:\s*(\d+\.\d+)\s*sec")
        log_text = ""
        for sim in range(run.macro_params.macro_simulations):
            log_file = os.path.join(run.os_path, f"{sim:02}", f"macro{file_code}_{sim:02}.txt")
            with open(log_file) as fh:
                log_text += fh.read()
        run._cache[key] = np.array(re.findall(pattern, log_text), dtype=float)
    return run._cache[key]


def get_total_binds(
    run: "Run",
    file_code: str,
) -> np.ndarray:
    """Read total bind counts from macroscale log files.

    :param run: Run object supplying ``os_path`` and ``macro_simulations``.
    :type run: Run
    :param file_code: File-code string used to construct log file names.
    :type file_code: str
    :return: 1-D integer array of total bind counts, one per simulation.
    :rtype: numpy.ndarray
    """
    key = ("total_binds", file_code)
    if key not in run._cache:
        pattern = re.compile(r"Total Binds:\s*(\d+\.?\d*)\s*")
        log_text = ""
        for sim in range(run.macro_params.macro_simulations):
            log_file = os.path.join(run.os_path, f"{sim:02}", f"macro{file_code}_{sim:02}.txt")
            with open(log_file) as fh:
                log_text += fh.read()
        run._cache[key] = np.array(re.findall(pattern, log_text), dtype=int)
    return run._cache[key]


###############################################################################
# Summary statistics
###############################################################################


_DEFAULT_PERCENT_MARKERS = [0.0, 0.25, 0.5, 0.75, 1.0]


def compute_run_statistics(
    run: "Run",
    percent_markers: list[float] = None,
) -> pd.Series:
    """Compute key summary statistics for a Run as a labelled pandas Series.

    Collects six scalar statistics (each as a mean and standard deviation)
    across all macroscale simulations and returns them as a
    :class:`pandas.Series` with a two-level
    :class:`~pandas.MultiIndex` ``(metric_name, stat)`` where *stat* is either
    ``"Mean"`` or ``"Standard Deviation"``.  Multiple Series can be stacked
    into a :class:`~pandas.DataFrame` with one row per Run via
    ``pd.DataFrame([s1, s2, ...])``.

    The six metrics mirror the results block in the ``MASTER LOOP`` notebook
    cell (lines 73-96):

    - **Degradation rate (%/min)** — slope of the rapid-degradation linear
      fit, converted to percent per minute.  Calls
      :func:`mean_degradation_rate` internally.
    - **Lysis lag time (min)** — time of the first save point in the
      rapid-degradation phase.  Also from :func:`mean_degradation_rate`.
    - **Time to full clot degradation (min)** — time at which the last
      milestone in *percent_markers* (i.e. 100 %) is first reached.  Calls
      :func:`find_degradation_marker_times` internally.
    - **Percent of molecules that reached the back row** — fraction of tPA
      molecules whose first-passage time is positive, expressed as a
      percentage.  Averaged and std'd across simulations.
    - **First passage time (min)** — mean and std of all positive first-
      passage times (in minutes) pooled across all simulations.
    - **Front Velocity (microns/min)** — mean and std from
      :func:`mean_front_velocity`.

    :param run: Run object supplying grid parameters and data.
    :type run: Run
    :param percent_markers: Degradation milestones used to locate the
        full-degradation time.  Must include ``1.0`` as the last element.
        Defaults to ``[0.0, 0.25, 0.5, 0.75, 1.0]``.
    :type percent_markers: list[float], optional
    :return: Series indexed by
        ``(metric_name, "Mean" | "Standard Deviation")``.
    :rtype: pandas.Series
    """
    if percent_markers is None:
        percent_markers = _DEFAULT_PERCENT_MARKERS

    key = (_KEY_RUN_STATS, tuple(percent_markers))
    if key in run._cache:
        return run._cache[key]

    deg_rate, _offset, deg_start = mean_degradation_rate(run)
    marker_times = find_degradation_marker_times(run, percent_markers)
    front_mean, front_sd = mean_front_velocity(run)

    n_sims = run.macro_params.macro_simulations
    total_molecules = run.macro_params.total_molecules
    reach_pct = np.array(
        [
            np.count_nonzero(run.data.macroscale_out[sim].tpa_transit_time[:] > 0)
            / total_molecules
            * 100
            for sim in range(n_sims)
        ]
    )
    all_transit = np.concatenate(
        [run.data.macroscale_out[sim].tpa_transit_time[:] for sim in range(n_sims)]
    )
    fpt = all_transit[all_transit > 0] / 60

    index = pd.MultiIndex.from_tuples(
        [
            ("Degradation rate (%/min)", "Mean"),
            ("Degradation rate (%/min)", "Standard Deviation"),
            ("Lysis lag time (min)", "Mean"),
            ("Lysis lag time (min)", "Standard Deviation"),
            ("Time to full clot degradation (min)", "Mean"),
            ("Time to full clot degradation (min)", "Standard Deviation"),
            ("Percent of molecules that reached the back row", "Mean"),
            ("Percent of molecules that reached the back row", "Standard Deviation"),
            ("First passage time (min)", "Mean"),
            ("First passage time (min)", "Standard Deviation"),
            ("Front Velocity (microns/min)", "Mean"),
            ("Front Velocity (microns/min)", "Standard Deviation"),
        ]
    )
    values = [
        float(np.mean(deg_rate) * 100),
        float(np.std(deg_rate) * 100),
        float(np.mean(deg_start)),
        float(np.std(deg_start)),
        float(np.mean(marker_times[:, -1])),
        float(np.std(marker_times[:, -1])),
        float(np.mean(reach_pct)),
        float(np.std(reach_pct)),
        float(np.mean(fpt)),
        float(np.std(fpt)),
        front_mean,
        front_sd,
    ]
    result = pd.Series(values, index=index)
    run._cache[key] = result
    return result


###############################################################################
# Aggregated stat helpers (mean ± std across simulations)
###############################################################################


def compute_degradation_marker_stats(
    run: "Run",
    percent_markers: list[int],
) -> dict[int, tuple[float, float]]:
    """Compute mean and standard deviation of degradation-milestone times.

    Calls :func:`find_degradation_marker_times` and aggregates across
    simulations, returning results in minutes.

    :param run: Run object supplying macro grid parameters and data.
    :type run: Run
    :param percent_markers: Integer percentage milestones, e.g.
        ``[5, 20, 50, 80, 100]``.
    :type percent_markers: list[int]
    :return: Mapping from each integer milestone to
        ``(mean_min, std_min)``.
    :rtype: dict[int, tuple[float, float]]
    """
    pct_fractions = [m / 100 for m in percent_markers]
    times = find_degradation_marker_times(run, pct_fractions)
    return {
        m: (float(times[:, k].mean()), float(times[:, k].std()))
        for k, m in enumerate(percent_markers)
    }


def _build_percent_markers(intervals: list[tuple[int, int]]) -> list[float]:
    """Return a sorted list of fractional markers covering all interval endpoints.

    Always includes 0.0 and 1.0 as bookends.

    :param intervals: List of (start, end) integer percent pairs.
    :type intervals: list[tuple[int, int]]
    :return: Sorted list of fractions with no duplicates.
    :rtype: list[float]
    """
    endpoints = {0, 100}
    for start, end in intervals:
        endpoints.add(start)
        endpoints.add(end)
    return [p / 100 for p in sorted(endpoints)]


def compute_degradation_rate_stats(
    run: "Run",
    intervals: list[tuple[int, int]],
) -> dict[tuple[int, int], tuple[float, float]]:
    """Compute mean and standard deviation of degradation rates for each interval.

    Calls :func:`degradation_rates` and aggregates across simulations,
    converting results from fraction/min to %/min.

    :param run: Run object supplying macro grid parameters and data.
    :type run: Run
    :param intervals: List of ``(start, end)`` integer percentage pairs, e.g.
        ``[(20, 80), (20, 50), (50, 80)]``.
    :type intervals: list[tuple[int, int]]
    :return: Mapping from each ``(start, end)`` pair to
        ``(mean_pct_per_min, std_pct_per_min)``.
    :rtype: dict[tuple[int, int], tuple[float, float]]
    """
    slope_pairs = [(s / 100, e / 100) for s, e in intervals]
    pct_markers = _build_percent_markers(intervals)
    rates = degradation_rates(run, slope_pairs, pct_markers)
    rates_pct = rates * 100
    return {
        ivl: (float(rates_pct[:, k].mean()), float(rates_pct[:, k].std()))
        for k, ivl in enumerate(intervals)
    }
