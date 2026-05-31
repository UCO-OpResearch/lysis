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
    _extract_data_arrays,
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


def _wrap_table(value):
    """Wrap a table value, mirroring the real view layer.

    A ``None`` value models an *absent optional* dataset: it is still listed
    by ``.datasets`` but dot-access returns ``None`` (the issue #103
    sentinel) rather than an :class:`_StubDataset`.
    """
    return None if value is None else _StubDataset(value)


class _StubCollection:
    """Minimal stand-in for ``DataCollection`` (simulations_combined=True).

    A table value of ``None`` models an absent optional dataset: listed by
    ``.datasets`` but returned as ``None`` from dot-access.
    """

    def __init__(self, tables):
        self._tables = {k: _wrap_table(v) for k, v in tables.items()}

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
    """Minimal stand-in for ``SimulationView``.

    A table value of ``None`` models an absent optional dataset: listed by
    ``.datasets`` but returned as ``None`` from dot-access.
    """

    def __init__(self, tables):
        self._tables = {k: _wrap_table(v) for k, v in tables.items()}

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


# ---------------------------------------------------------------------------
# Absent optional datasets (issue #103)
# ---------------------------------------------------------------------------


class TestDiffRunsOptionalDatasets:
    """Absent optional datasets (dot-access returns ``None``) are skipped.

    Models the dispatcher-log scenario from issue #103: one run carries an
    optional dataset, the other lacks it.  ``lysis diff`` must complete
    rather than raising the ``component not found`` ``KeyError``.
    """

    def test_extract_skips_absent_optional_micro(self):
        """_extract_data_arrays drops a combined dataset whose access is None."""
        run = _make_run(
            micro_tables={
                "x": np.array([1.0, 2.0]),
                "micro_dispatcher_log": None,  # absent optional
            },
        )
        tables = _extract_data_arrays(run, {"microscale_out"})
        assert "microscale_out/x" in tables
        assert "microscale_out/micro_dispatcher_log" not in tables

    def test_extract_skips_absent_optional_macro(self):
        """_extract_data_arrays drops a per-sim dataset whose access is None."""
        run = _make_run(
            macro_per_sim=[
                {
                    "snapshot_time": np.array([0.0, 1.0]),
                    "macro_dispatcher_log": None,  # absent optional
                },
            ],
        )
        tables = _extract_data_arrays(run, {"macroscale_out"})
        assert "macroscale_out[00]/snapshot_time" in tables
        assert "macroscale_out[00]/macro_dispatcher_log" not in tables

    def test_diff_completes_when_optional_present_in_one_run(self):
        """A run with and a run without an optional dataset diff without error.

        The optional dataset is reported as present-in-one-run-only rather
        than raising, and shared required datasets still compare.
        """
        run1 = _make_run(
            micro_tables={
                "x": np.array([1.0, 2.0]),
                "micro_dispatcher_log": np.array(["log a"]),  # present
            },
        )
        run2 = _make_run(
            micro_tables={
                "x": np.array([1.0, 2.0]),
                "micro_dispatcher_log": None,  # absent optional
            },
        )
        result = diff_runs(run1, run2)
        diff = result["data_diff"]
        # Shared required dataset compared and matches.
        assert diff["microscale_out/x"]["status"] == "match"
        # Optional present only in run1 -> reported missing in run2.
        opt = diff["microscale_out/micro_dispatcher_log"]
        assert opt["status"] == "missing"
        assert opt["detail"] == "not in run2"

    def test_diff_completes_when_optional_absent_in_both(self):
        """Optional dataset absent from both runs is simply omitted."""
        run1 = _make_run(
            micro_tables={"x": np.array([1.0]), "micro_dispatcher_log": None},
        )
        run2 = _make_run(
            micro_tables={"x": np.array([1.0]), "micro_dispatcher_log": None},
        )
        result = diff_runs(run1, run2)
        labels = result["data_diff"].keys()
        assert "microscale_out/x" in labels
        assert "microscale_out/micro_dispatcher_log" not in labels
