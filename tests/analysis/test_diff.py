"""Unit tests for :mod:`lysis.analysis.diff`.

Exercises the auto-detection + element-wise diff path that ``lysis diff``
uses: :func:`available_scales`, :func:`list_tables`, and
:func:`diff_runs` against lightweight stub DataStores.  Float-tolerance
behaviour is covered by the primitives in ``test_compare.py``; here we
check the dispatch wiring, scale-skip reporting, and log-table
exclusion.
"""

import numpy as np
import pytest

from lysis.analysis.diff import (
    available_scales,
    diff_runs,
    list_tables,
)


# ---------------------------------------------------------------------------
# Stub DataStore / Run helpers
# ---------------------------------------------------------------------------


class _StubDataset:
    """Wrap a numpy array for h5py-style slice access."""

    def __init__(self, data):
        self._data = np.asarray(data)

    def __getitem__(self, key):
        return self._data[key]


class _StubCollection:
    """Minimal stand-in for ``DataCollection`` (simulations_combined=True)."""

    def __init__(self, tables):
        self._tables = {k: _StubDataset(v) for k, v in tables.items()}

    @property
    def datasets(self):
        return list(self._tables.keys())

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        if name in self._tables:
            return self._tables[name]
        raise AttributeError(name)


class _StubSimView:
    """Minimal stand-in for ``SimulationView``."""

    def __init__(self, tables):
        self._tables = {k: _StubDataset(v) for k, v in tables.items()}

    @property
    def datasets(self):
        return list(self._tables.keys())

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        if name in self._tables:
            return self._tables[name]
        raise AttributeError(name)


class _StubMacroCollection:
    """Per-simulation collection indexed by integer simulation id."""

    def __init__(self, per_sim_tables):
        self._views = [_StubSimView(t) for t in per_sim_tables]

    def __getitem__(self, idx):
        return self._views[idx]


class _StubDataStore:
    """Minimal DataStore-like object exposing ``collections`` and dot-access."""

    def __init__(self, micro_tables=None, macro_per_sim=None):
        self._colls = {}
        if micro_tables is not None:
            self._colls["microscale_out"] = _StubCollection(micro_tables)
        if macro_per_sim is not None:
            self._colls["macroscale_out"] = _StubMacroCollection(macro_per_sim)

    @property
    def collections(self):
        return dict(self._colls)

    @property
    def microscale_out(self):
        return self._colls["microscale_out"]

    @property
    def macroscale_out(self):
        return self._colls["macroscale_out"]


class _StubParams:
    def __init__(self, macro_simulations):
        self.macro_simulations = macro_simulations


class _StubRunForData:
    def __init__(self, data, macro_simulations=0):
        self.data = data
        self.macro_params = _StubParams(macro_simulations)


def _make_run(micro_tables=None, macro_per_sim=None):
    macro_sims = 0 if macro_per_sim is None else len(macro_per_sim)
    return _StubRunForData(
        _StubDataStore(micro_tables=micro_tables, macro_per_sim=macro_per_sim),
        macro_simulations=macro_sims,
    )


# ---------------------------------------------------------------------------
# available_scales
# ---------------------------------------------------------------------------


class TestAvailableScales:
    def test_both_scales_present(self):
        run = _make_run(
            micro_tables={"x": np.array([1.0])},
            macro_per_sim=[{"y": np.array([2.0])}],
        )
        assert available_scales(run) == {"microscale_out", "macroscale_out"}

    def test_micro_only(self):
        run = _make_run(micro_tables={"x": np.array([1.0])})
        assert available_scales(run) == {"microscale_out"}

    def test_macro_only(self):
        run = _make_run(
            macro_per_sim=[{"y": np.array([2.0])}],
        )
        assert available_scales(run) == {"macroscale_out"}

    def test_neither_present(self):
        run = _make_run()
        assert available_scales(run) == set()


# ---------------------------------------------------------------------------
# list_tables
# ---------------------------------------------------------------------------


class TestListTables:
    def test_returns_sorted_labels_for_both_scales(self):
        run = _make_run(
            micro_tables={"fiber_degraded": np.array([True])},
            macro_per_sim=[{"snapshot_time": np.array([0.0, 1.0])}],
        )
        labels = list_tables(run)
        assert labels == sorted(labels)
        assert "microscale_out/fiber_degraded" in labels
        assert "macroscale_out[00]/snapshot_time" in labels

    def test_skips_log_tables(self):
        run = _make_run(
            micro_tables={
                "fiber_degraded": np.array([True]),
                "micro_log": np.array(["x"]),
            },
        )
        labels = list_tables(run)
        assert "microscale_out/fiber_degraded" in labels
        assert "microscale_out/micro_log" not in labels


# ---------------------------------------------------------------------------
# diff_runs
# ---------------------------------------------------------------------------


class TestDiffRunsIdenticalMicroOnly:
    def test_all_match_and_scales_compared_is_micro(self):
        run1 = _make_run(
            micro_tables={
                "fiber_degraded": np.array([True, False, True]),
                "tpa_leaving_time": np.array([1.0, 2.0, 3.0]),
            }
        )
        run2 = _make_run(
            micro_tables={
                "fiber_degraded": np.array([True, False, True]),
                "tpa_leaving_time": np.array([1.0, 2.0, 3.0]),
            }
        )
        result = diff_runs(run1, run2)
        assert set(result.keys()) == {
            "data_diff",
            "scales_compared",
            "scales_skipped_in_run1",
            "scales_skipped_in_run2",
        }
        assert result["scales_compared"] == ["microscale_out"]
        assert result["scales_skipped_in_run1"] == []
        assert result["scales_skipped_in_run2"] == []
        for entry in result["data_diff"].values():
            assert entry["status"] == "match"


class TestDiffRunsReportsDiff:
    def test_diff_location_reported(self):
        run1 = _make_run(
            micro_tables={"tpa_leaving_time": np.array([1.0, 2.0, 3.0])}
        )
        run2 = _make_run(
            micro_tables={"tpa_leaving_time": np.array([1.0, 2.0, 4.0])}
        )
        result = diff_runs(run1, run2)
        entry = result["data_diff"]["microscale_out/tpa_leaving_time"]
        assert entry["status"] == "diff"
        assert entry["location"] == (2,)
        assert entry["max_pct_diff"] > 0


class TestDiffRunsMissing:
    def test_extra_table_in_run2_reported_as_missing(self):
        run1 = _make_run(micro_tables={"only_in_1": np.array([1.0])})
        run2 = _make_run(micro_tables={"shared": np.array([1.0])})
        # run2 has a table that run1 doesn't
        run2.data.collections["microscale_out"]._tables["extra"] = _StubDataset(
            np.array([9.0])
        )
        result = diff_runs(run1, run2)["data_diff"]
        assert result["microscale_out/only_in_1"]["status"] == "missing"
        assert "run2" in result["microscale_out/only_in_1"]["detail"]
        assert result["microscale_out/extra"]["status"] == "missing"
        assert "run1" in result["microscale_out/extra"]["detail"]


class TestDiffRunsScaleMismatch:
    def test_run1_has_both_run2_micro_only_skips_macro(self):
        run1 = _make_run(
            micro_tables={"x": np.array([1.0])},
            macro_per_sim=[{"snapshot_time": np.array([0.0, 1.0])}],
        )
        run2 = _make_run(micro_tables={"x": np.array([1.0])})
        result = diff_runs(run1, run2)
        assert result["scales_compared"] == ["microscale_out"]
        assert result["scales_skipped_in_run2"] == ["macroscale_out"]
        assert result["scales_skipped_in_run1"] == []
        # Only micro labels present
        assert all(
            label.startswith("microscale_out/")
            for label in result["data_diff"].keys()
        )

    def test_run2_has_both_run1_micro_only_skips_macro(self):
        run1 = _make_run(micro_tables={"x": np.array([1.0])})
        run2 = _make_run(
            micro_tables={"x": np.array([1.0])},
            macro_per_sim=[{"snapshot_time": np.array([0.0, 1.0])}],
        )
        result = diff_runs(run1, run2)
        assert result["scales_compared"] == ["microscale_out"]
        assert result["scales_skipped_in_run1"] == ["macroscale_out"]
        assert result["scales_skipped_in_run2"] == []


class TestDiffRunsEmptyIntersection:
    def test_no_common_scales_raises(self):
        run1 = _make_run(micro_tables={"x": np.array([1.0])})
        run2 = _make_run(
            macro_per_sim=[{"snapshot_time": np.array([0.0, 1.0])}],
        )
        with pytest.raises(ValueError, match="No comparable data"):
            diff_runs(run1, run2)

    def test_neither_scale_present_raises(self):
        run1 = _make_run()
        run2 = _make_run()
        with pytest.raises(ValueError, match="No comparable data"):
            diff_runs(run1, run2)


class TestDiffRunsIncludesMacro:
    def test_macro_tables_with_sim_index_prefix(self):
        run1 = _make_run(
            micro_tables={"x": np.array([1.0])},
            macro_per_sim=[
                {"snapshot_time": np.array([0.0, 1.0])},
                {"snapshot_time": np.array([0.0, 1.0])},
            ],
        )
        run2 = _make_run(
            micro_tables={"x": np.array([1.0])},
            macro_per_sim=[
                {"snapshot_time": np.array([0.0, 1.0])},
                {"snapshot_time": np.array([0.0, 1.0])},
            ],
        )
        result = diff_runs(run1, run2)
        assert "macroscale_out[00]/snapshot_time" in result["data_diff"]
        assert "macroscale_out[01]/snapshot_time" in result["data_diff"]

    def test_macro_log_tables_skipped(self):
        run1 = _make_run(
            macro_per_sim=[{"macro_log": np.array(["x"]),
                           "snapshot_time": np.array([0.0])}],
        )
        run2 = _make_run(
            macro_per_sim=[{"macro_log": np.array(["y"]),
                           "snapshot_time": np.array([0.0])}],
        )
        labels = diff_runs(run1, run2)["data_diff"].keys()
        assert "macroscale_out[00]/snapshot_time" in labels
        assert all(not label.endswith("macro_log") for label in labels)
