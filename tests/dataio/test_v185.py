"""Tests for the v1.85.0 data specification and its converters.

v1.85.0 predates the one-directory-per-simulation convention: its
``macroscale_out`` collection is ``simulations_combined=True``, with every
simulation concatenated into a single top-level file per dataset.  Recovering
per-simulation arrays requires the ``Nsave`` vector, and three datasets need
special handling on the way to v1.90.0:

* ``f_deg_time`` -- v1.85.0's ``0.0`` sentinel is ambiguous (empty edge *or*
  unscheduled fibrin) and must be disambiguated by index;
* ``m_bind_t`` -- never written by v1.85.0, so it is reconstructed by
  differencing the ``m_bound`` snapshots;
* ``macro_log`` -- one combined log split on ``run number=`` markers.
"""

import numpy as np
import pytest

from lysis.config.constants import (
    CONST,
    EMPTY_EDGE_DEGRADE_TIME,
    UNSCHEDULED_DEGRADE_TIME,
    V185_UNSCHEDULED_DEGRADE_TIME,
)
from lysis.dataio.dataspec import CONVERSION_WARNINGS, dataspec, fortran_versions
from lysis.dataio.dataconvert import (
    combine_per_simulation,
    conversion_paths,
    convert_bind_events_to_bound,
    convert_f_deg_time_v185_to_v190,
    convert_f_deg_time_v190_to_v185,
    convert_m_bound_to_m_bind_t,
    data_converters,
    join_macro_log,
    params_converters,
    simulation_boundaries,
    split_combined_by_simulation,
    split_combined_by_snapshot,
    split_macro_log,
    _convert_params_v185_to_v190,
    _convert_params_v190_to_v185,
    _convert_params_v190_to_v195,
)

# Two simulations with different snapshot counts, so any off-by-one in the
# Nsave bookkeeping shows up as a shape mismatch rather than silently passing.
NSAVE = np.array([2, 3], dtype=np.int32)
SNAPSHOTS = int((NSAVE + 1).sum())  # 3 + 4 == 7
N_EDGES = 6
N_MOLECULES = 4
EMPTY_EDGES = 2


def _params(empty_edges=EMPTY_EDGES):
    """Build a minimal v1.85.0-style params dict."""
    return {
        "macro_params": {
            "last_empty_edge": empty_edges - 1,
            "total_edges": N_EDGES,
            "total_molecules": N_MOLECULES,
        }
    }


class TestSpec:
    """The v1.85.0 spec itself."""

    def test_registered(self):
        assert "v1.85.0" in dataspec
        assert "v1.85.0" in fortran_versions

    def test_macroscale_out_is_combined(self):
        assert dataspec["v1.85.0"]["macroscale_out"].simulations_combined is True
        # ...unlike every other Fortran version.
        assert dataspec["v1.90.0"]["macroscale_out"].simulations_combined is False

    def test_data_locations_have_no_simulation_directory(self):
        for name, spec in dataspec["v1.85.0"]["macroscale_out"].data.items():
            assert "{sim" not in spec.data_location, name
            assert "/" not in spec.data_location, name

    def test_nsave_is_a_vector(self):
        """Combined Nsave holds one count per simulation, not a scalar."""
        assert dataspec["v1.85.0"]["macroscale_out"].data["Nsave"].shape == (-1,)
        assert dataspec["v1.90.0"]["macroscale_out"].data["Nsave"].shape == ()

    def test_mfpt_declares_molecule_dimension(self):
        """mfpt is one row per simulation, so reshape needs the second axis."""
        assert dataspec["v1.85.0"]["macroscale_out"].data["mfpt"].shape == (
            -1,
            "macro_params.total_molecules",
        )

    def test_m_bind_t_removed_not_optional(self):
        """Must be absent so _convert_single_step reaches the reconstruction.

        An optional-but-missing dataset is skipped outright, which would
        silently drop the reconstructed binding history.
        """
        assert "m_bind_t" not in dataspec["v1.85.0"]["macroscale_out"].data
        assert "m_bind_t" in dataspec["v1.90.0"]["macroscale_out"].data

    def test_deg_is_optional_and_v185_only(self):
        assert dataspec["v1.85.0"]["macroscale_out"].data["deg"].optional
        assert "deg" not in dataspec["v1.90.0"]["macroscale_out"].data

    def test_other_collections_unchanged_from_v190(self):
        """microscale_out and macroscale_in are identical to v1.90.0."""

        def key(version, collection):
            return {
                name: (spec.data_location, spec.shape, spec.optional, spec.dtype)
                for name, spec in dataspec[version][collection].data.items()
            }

        for collection in ("microscale_out", "macroscale_in"):
            assert key("v1.85.0", collection) == key("v1.90.0", collection)
            assert (
                dataspec["v1.85.0"][collection].simulations_combined
                is dataspec["v1.90.0"][collection].simulations_combined
            )

    def test_conversion_warning_registered(self):
        """converted_from records only the FIRST hop, so v1.85.0 needs its own."""
        assert "v1.85.0" in CONVERSION_WARNINGS
        warning = CONVERSION_WARNINGS["v1.85.0"]
        assert "tpa_bind_events" in warning
        assert "f_deg_list" in warning


class TestConversionGraph:
    """Registration in the conversion registry and routing graph."""

    def test_both_directions_registered(self):
        """A one-way registration is silently dropped from the graph."""
        assert ("v1.85.0", "v1.90.0") in data_converters
        assert ("v1.90.0", "v1.85.0") in data_converters
        assert ("v1.85.0", "v1.90.0") in params_converters
        assert ("v1.90.0", "v1.85.0") in params_converters

    def test_path_to_hdf5(self):
        assert conversion_paths[("v1.85.0", "v2.0.0")] == [
            "v1.85.0",
            "v1.90.0",
            "v1.95.0",
            "v1.99.0",
            "v2.0.0",
        ]

    def test_path_from_hdf5(self):
        assert conversion_paths[("v2.0.0", "v1.85.0")] == [
            "v2.0.0",
            "v1.99.0",
            "v1.95.0",
            "v1.90.0",
            "v1.85.0",
        ]

    def test_v185_is_the_spanning_tree_root(self):
        """min(versions) picks v1.85.0 now; every pair must still route."""
        for source in ("v1.90.0", "v1.95.0", "v1.99.0", "v2.0.0"):
            for target in ("v1.90.0", "v1.95.0", "v1.99.0", "v2.0.0"):
                if source != target:
                    assert (source, target) in conversion_paths


class TestSplitAndCombine:
    """The generic split/join helpers."""

    def test_boundaries_account_for_the_initial_snapshot(self):
        """Simulation i occupies Nsave[i] + 1 rows, not Nsave[i]."""
        assert list(simulation_boundaries({"Nsave": NSAVE})) == [3]

    def test_split_by_snapshot(self):
        data = {"Nsave": NSAVE, "tsave": np.arange(SNAPSHOTS, dtype=np.float64)}
        parts = split_combined_by_snapshot(data, "tsave")
        assert [len(p) for p in parts] == [3, 4]
        assert np.array_equal(parts[0], [0, 1, 2])
        assert np.array_equal(parts[1], [3, 4, 5, 6])

    def test_split_by_snapshot_2d(self):
        data = {
            "Nsave": NSAVE,
            "m_loc": np.arange(SNAPSHOTS * N_MOLECULES).reshape(SNAPSHOTS, N_MOLECULES),
        }
        parts = split_combined_by_snapshot(data, "m_loc")
        assert [p.shape for p in parts] == [(3, N_MOLECULES), (4, N_MOLECULES)]

    def test_split_by_simulation(self):
        data = {"mfpt": np.arange(2 * N_MOLECULES).reshape(2, N_MOLECULES)}
        parts = split_combined_by_simulation(data, "mfpt")
        assert len(parts) == 2
        assert np.array_equal(parts[1], np.arange(N_MOLECULES, 2 * N_MOLECULES))

    def test_single_simulation(self):
        """A one-simulation run has no interior split points."""
        data = {"Nsave": np.array([2]), "tsave": np.arange(3.0)}
        parts = split_combined_by_snapshot(data, "tsave")
        assert len(parts) == 1
        assert np.array_equal(parts[0], [0, 1, 2])

    @pytest.mark.parametrize("stack", [False, True])
    def test_round_trip(self, stack):
        if stack:
            original = np.arange(2 * N_MOLECULES, dtype=np.float64).reshape(
                2, N_MOLECULES
            )
            parts = split_combined_by_simulation({"mfpt": original}, "mfpt")
        else:
            original = np.arange(SNAPSHOTS * N_EDGES, dtype=np.float64).reshape(
                SNAPSHOTS, N_EDGES
            )
            parts = split_combined_by_snapshot({"Nsave": NSAVE, "f": original}, "f")
        rebuilt = combine_per_simulation({"f": parts}, "f", stack=stack)
        assert np.array_equal(rebuilt, original)


class TestFDegTimeSentinel:
    """The index-restricted 0.0 <-> 9.9e100 remap.

    v1.85.0 wrote 0.0 for both empty edges and unscheduled fibrin; v1.90.0
    keeps 0.0 for empty edges only and marks unscheduled fibrin with 9.9e100.
    """

    def _combined(self):
        # Columns 0-1 are empty edges (always 0.0); 2-5 are fibrin.
        f = np.zeros((SNAPSHOTS, N_EDGES), dtype=np.float64)
        f[1:, 2] = 12.5  # one fibrin edge scheduled from snapshot 1 on
        f[2:, 3] = 30.0
        return {"Nsave": NSAVE, "f_deg_time": f, "params": _params()}

    def test_empty_edges_keep_their_zero(self):
        """A blanket value remap would destroy the empty-edge marker."""
        out = convert_f_deg_time_v185_to_v190(self._combined())
        for part in out:
            assert np.all(part[:, :EMPTY_EDGES] == EMPTY_EDGE_DEGRADE_TIME)

    def test_unscheduled_fibrin_gets_the_sentinel(self):
        out = convert_f_deg_time_v185_to_v190(self._combined())
        # Edges 4 and 5 are never scheduled in any snapshot.
        for part in out:
            assert np.all(part[:, 4:] == UNSCHEDULED_DEGRADE_TIME)

    def test_scheduled_values_survive(self):
        out = convert_f_deg_time_v185_to_v190(self._combined())
        assert out[0][1, 2] == 12.5
        assert out[0][2, 3] == 30.0
        # ...but the pre-schedule snapshot becomes the sentinel.
        assert out[0][0, 2] == UNSCHEDULED_DEGRADE_TIME

    def test_splits_by_simulation(self):
        out = convert_f_deg_time_v185_to_v190(self._combined())
        assert [p.shape for p in out] == [(3, N_EDGES), (4, N_EDGES)]

    def test_input_not_mutated(self):
        data = self._combined()
        before = data["f_deg_time"].copy()
        convert_f_deg_time_v185_to_v190(data)
        assert np.array_equal(data["f_deg_time"], before)

    def test_round_trip(self):
        data = self._combined()
        original = data["f_deg_time"].copy()
        per_sim = convert_f_deg_time_v185_to_v190(data)
        rebuilt = convert_f_deg_time_v190_to_v185({"f_deg_time": per_sim})
        assert np.array_equal(rebuilt, original)

    def test_no_empty_edges(self):
        """empty_edges == 0 (as in the v1.90.0 fixture) remaps everything."""
        data = self._combined()
        data["params"]["macro_params"]["last_empty_edge"] = -1
        out = convert_f_deg_time_v185_to_v190(data)
        assert np.all(out[0][0] == UNSCHEDULED_DEGRADE_TIME)

    def test_sentinels_are_distinct(self):
        """The whole remap is pointless if these ever collapse."""
        assert V185_UNSCHEDULED_DEGRADE_TIME == EMPTY_EDGE_DEGRADE_TIME
        assert UNSCHEDULED_DEGRADE_TIME != EMPTY_EDGE_DEGRADE_TIME


class TestMBindTReconstruction:
    """Rebuilding the binding event log from m_bound snapshots."""

    def _combined(self):
        # 7 snapshots x 4 molecules across 2 simulations (3 + 4 rows).
        bound = np.zeros((SNAPSHOTS, N_MOLECULES), dtype=np.int32)
        bound[1, 0] = 1  # sim 0: molecule 0 binds at snapshot 1...
        bound[2, 0] = 0  # ...and unbinds at snapshot 2
        bound[4:, 2] = 1  # sim 1: molecule 2 binds at its snapshot 1 and stays
        loc = np.arange(1, SNAPSHOTS * N_MOLECULES + 1, dtype=np.int32).reshape(
            SNAPSHOTS, N_MOLECULES
        )
        tsave = np.concatenate([np.arange(3.0), np.arange(4.0)])
        return {
            "Nsave": NSAVE,
            "m_bound": bound,
            "m_loc": loc,
            "tsave": tsave,
            "params": _params(),
        }

    def test_one_event_per_status_change(self):
        out = convert_m_bound_to_m_bind_t(self._combined())
        assert len(out) == 2
        assert len(out[0]) == 2  # bind then unbind
        assert len(out[1]) == 1  # binds once and stays bound: no repeat events

    def test_molecule_index_is_one_based(self):
        """v1.90.0 m_bind_t uses Fortran's 1-based molecule indexing."""
        out = convert_m_bound_to_m_bind_t(self._combined())
        assert out[0]["tPA Molecule Index"][0] == 1  # molecule 0
        assert out[1]["tPA Molecule Index"][0] == 3  # molecule 2

    def test_status_values_are_mol_status(self):
        out = convert_m_bound_to_m_bind_t(self._combined())
        assert out[0]["Molecule New Status"][0] == CONST.MOL_STATUS.BOUND
        assert out[0]["Molecule New Status"][1] == CONST.MOL_STATUS.UNBOUND

    def test_location_taken_from_matching_snapshot(self):
        data = self._combined()
        out = convert_m_bound_to_m_bind_t(data)
        # sim 0, molecule 0 binds at snapshot 1 -> m_loc[1, 0]
        assert out[0]["Grid Location Index"][0] == data["m_loc"][1, 0]

    def test_times_come_from_the_later_snapshot(self):
        out = convert_m_bound_to_m_bind_t(self._combined())
        assert out[0]["Simulation Time Elapsed"][0] == 1.0
        assert out[0]["Simulation Time Elapsed"][1] == 2.0

    def test_events_are_time_sorted(self):
        """replay_event_log_to_snapshot searchsorts on this column."""
        out = convert_m_bound_to_m_bind_t(self._combined())
        for events in out:
            times = events["Simulation Time Elapsed"]
            assert np.all(np.diff(times) >= 0)

    def test_molecules_bound_in_first_snapshot_are_seeded(self):
        """There is no earlier snapshot to differ from, so seed explicitly."""
        data = self._combined()
        data["m_bound"][0, 1] = 1  # bound at t=0, then unbound at snapshot 1
        out = convert_m_bound_to_m_bind_t(data)
        molecule = out[0][out[0]["tPA Molecule Index"] == 2]
        assert molecule["Simulation Time Elapsed"][0] == 0.0
        assert molecule["Molecule New Status"][0] == CONST.MOL_STATUS.BOUND

    def test_round_trip_through_replay(self):
        """Replaying the reconstructed log must reproduce m_bound exactly.

        This is the strongest available check: it uses the pre-existing
        inverse, :func:`convert_bind_events_to_bound`, rather than a
        reimplementation.
        """
        data = self._combined()
        events = convert_m_bound_to_m_bind_t(data)
        zero_based = []
        for log in events:
            log = log.copy()
            log["tPA Molecule Index"] -= 1  # replay expects 0-based indices
            zero_based.append(log)
        replayed = convert_bind_events_to_bound(
            {
                "m_bind_t": zero_based,
                "tsave": split_combined_by_snapshot(data, "tsave"),
                "params": data["params"],
            },
            input_dataset="m_bind_t",
            snapshot_dataset="tsave",
        )
        expected = split_combined_by_snapshot(data, "m_bound")
        for got, want in zip(replayed, expected):
            assert np.array_equal(got.astype(np.int32), want.astype(np.int32))

    def test_no_events(self):
        data = self._combined()
        data["m_bound"][:] = 0
        out = convert_m_bound_to_m_bind_t(data)
        assert all(len(events) == 0 for events in out)


class TestMacroLog:
    """Splitting the single combined log into per-simulation logs."""

    LOG = np.array(
        [
            " N=  93",
            " seed=  17109424",
            "  run number=           1",
            " sim 0 line",
            "  run number=           2",
            " sim 1 line a",
            " sim 1 line b",
        ]
    )

    def test_one_log_per_run_marker(self):
        out = split_macro_log({"macro_log": self.LOG})
        assert len(out) == 2

    def test_header_prepended_to_every_simulation(self):
        """Each per-simulation log must stay self-describing."""
        out = split_macro_log({"macro_log": self.LOG})
        for log in out:
            assert list(log[:2]) == [" N=  93", " seed=  17109424"]

    def test_body_goes_to_the_right_simulation(self):
        out = split_macro_log({"macro_log": self.LOG})
        assert " sim 0 line" in list(out[0])
        assert " sim 0 line" not in list(out[1])
        assert list(out[1][2:]) == [
            "  run number=           2",
            " sim 1 line a",
            " sim 1 line b",
        ]

    def test_no_markers_yields_one_log(self):
        """A truncated log degrades gracefully instead of raising."""
        out = split_macro_log({"macro_log": np.array([" N=  93", " nothing"])})
        assert len(out) == 1
        assert len(out[0]) == 2

    def test_round_trip(self):
        out = split_macro_log({"macro_log": self.LOG})
        assert list(join_macro_log({"macro_log": out})) == list(self.LOG)


class TestParams:
    """v1.85.0 legacy parameter names."""

    V185 = {
        "experiment_code": "2023-02-02-2200",
        "data_filenames": {"lysis_time": "lysismat.dat"},
        "micro_params": None,
        "macro_params": {
            "total_trials": 10,
            "seed": 17109424,
            "log_lvl": 30,
            "average_bind_time": 27.8,
            "last_empty_edge": 7783,
            "total_edges": 33545,
            "binding_rate": 0.1,
            "binding_sites": 427,
            "grid_node_distance": 1.0862,
            "microscale_runs": 50000,
            "input_data": ["unbinding_time"],
            "output_data": ["save_time"],
        },
    }

    def test_top_level_bookkeeping_dropped(self):
        """data_filenames is a dict and would be read as a params section."""
        out = _convert_params_v185_to_v190(self.V185)
        assert "experiment_code" not in out
        assert "data_filenames" not in out

    def test_macro_params_untouched(self):
        """v1.85.0 already uses the legacy spellings v1.90.0 expects.

        The ``total_trials``/``seed``/``log_lvl`` renames belong to the
        v1.90.0 -> v1.95.0 hop; duplicating them here would be redundant.
        """
        out = _convert_params_v185_to_v190(self.V185)["macro_params"]
        assert out == self.V185["macro_params"]

    def test_legacy_names_resolve_at_the_next_hop(self):
        v190 = _convert_params_v185_to_v190(self.V185)
        v195 = _convert_params_v190_to_v195(v190)["macro_params"]
        assert v195["macro_simulations"] == 10
        assert v195["macro_seed"] == 17109424
        assert v195["macro_log_lvl"] == 30

    def test_null_micro_params_passes_through(self):
        out = _convert_params_v185_to_v190(self.V185)
        assert out["micro_params"] is None

    def test_round_trip(self):
        out = _convert_params_v190_to_v185(_convert_params_v185_to_v190(self.V185))
        assert out["macro_params"] == self.V185["macro_params"]
