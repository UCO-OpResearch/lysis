"""Tests for ``lysis.analysis.summary`` — DataFrame-building functions."""

import pandas as pd
import pytest

from lysis.analysis.summary import (
    MICRO_STATS_COLUMNS,
    _fmt_micro,
    compare_stats_table,
    deg_rate_table,
    deg_time_table,
    macro_stats_table,
    micro_stats_table,
    parameters_table,
)

# ---------------------------------------------------------------------------
# Shared test data
# ---------------------------------------------------------------------------

_MACRO_METRICS = [
    "Degradation rate (%/min)",
    "Lysis lag time (min)",
    "Time to full clot degradation (min)",
    "Percent of molecules that reached the back row",
    "First passage time (min)",
    "Front Velocity (microns/min)",
]

_DEFAULT_MARKERS = [5, 20, 50, 80, 100]
_DEFAULT_INTERVALS = [(20, 80), (20, 50), (50, 80)]


def _micro_stats(
    fibers=1234,
    lt_mean=42.0,
    lt_std=3.5,
    lt_med=41.5,
    tpa_mean=15.0,
    tpa_std=2.0,
    tpa_med=14.5,
):
    return {
        "fibers_degraded": fibers,
        "lysis_time_mean": lt_mean,
        "lysis_time_std": lt_std,
        "lysis_time_median": lt_med,
        "tpa_leaving_mean": tpa_mean,
        "tpa_leaving_std": tpa_std,
        "tpa_leaving_median": tpa_med,
    }


def _macro_stats(mean=1.234, std=0.056):
    index = pd.MultiIndex.from_tuples(
        [(m, s) for m in _MACRO_METRICS for s in ["Mean", "Standard Deviation"]]
    )
    values = [
        mean if s == "Mean" else std
        for m in _MACRO_METRICS
        for s in ["Mean", "Standard Deviation"]
    ]
    return pd.Series(values, index=index)


# ---------------------------------------------------------------------------
# _fmt_micro
# ---------------------------------------------------------------------------


class TestFmtMicro:
    """Moved from tests/cli/test_micro_stats.py::TestFmtMetric."""

    def test_fibers_degraded_formatted_with_comma(self):
        assert _fmt_micro("Fibers Degraded", _micro_stats(fibers=1234)) == "1,234"

    def test_mean_lysis_time_shows_pm(self):
        result = _fmt_micro("Mean Lysis Time (min)", _micro_stats(lt_mean=42.0, lt_std=3.5))
        assert "\u00b1" in result
        assert "42.000" in result
        assert "3.500" in result

    def test_median_lysis_time_no_pm(self):
        result = _fmt_micro("Median Lysis Time (min)", _micro_stats(lt_med=41.5))
        assert "\u00b1" not in result
        assert "41.500" in result

    def test_mean_tpa_leaving_shows_pm(self):
        result = _fmt_micro("Mean tPA Leaving Time (sec)", _micro_stats(tpa_mean=15.0, tpa_std=2.0))
        assert "\u00b1" in result
        assert "15.000" in result
        assert "2.000" in result

    def test_median_tpa_leaving_no_pm(self):
        result = _fmt_micro("Median tPA Leaving Time (sec)", _micro_stats(tpa_med=14.5))
        assert "\u00b1" not in result
        assert "14.500" in result

    def test_unknown_key_raises(self):
        with pytest.raises(KeyError):
            _fmt_micro("Unknown Metric", _micro_stats())


# ---------------------------------------------------------------------------
# micro_stats_table
# ---------------------------------------------------------------------------


class TestMicroStatsTable:
    def test_returns_dataframe(self):
        df = micro_stats_table({"run_A": _micro_stats()})
        assert isinstance(df, pd.DataFrame)

    def test_index_is_run_codes(self):
        df = micro_stats_table({"run_A": _micro_stats(), "run_B": _micro_stats()})
        assert list(df.index) == ["run_A", "run_B"]

    def test_columns_match_MICRO_STATS_COLUMNS(self):
        df = micro_stats_table({"run_A": _micro_stats()})
        assert list(df.columns) == MICRO_STATS_COLUMNS

    def test_fibers_degraded_formatted_with_comma(self):
        df = micro_stats_table({"run_A": _micro_stats(fibers=1234)})
        assert df.loc["run_A", "Fibers Degraded"] == "1,234"

    def test_mean_lysis_time_shows_pm(self):
        df = micro_stats_table({"run_A": _micro_stats(lt_mean=42.0, lt_std=3.5)})
        val = df.loc["run_A", "Mean Lysis Time (min)"]
        assert "42.000" in val and "\u00b1" in val

    def test_median_lysis_time_no_pm(self):
        df = micro_stats_table({"run_A": _micro_stats(lt_med=41.5)})
        val = df.loc["run_A", "Median Lysis Time (min)"]
        assert "41.500" in val and "\u00b1" not in val

    def test_tpa_leaving_time_in_seconds(self):
        df = micro_stats_table({"run_A": _micro_stats(tpa_mean=15.0, tpa_std=2.0)})
        val = df.loc["run_A", "Mean tPA Leaving Time (sec)"]
        assert "15.000" in val and "2.000" in val

    def test_multiple_runs_independent_values(self):
        df = micro_stats_table(
            {"run_X": _micro_stats(fibers=1000), "run_Y": _micro_stats(fibers=500)}
        )
        assert df.loc["run_X", "Fibers Degraded"] == "1,000"
        assert df.loc["run_Y", "Fibers Degraded"] == "500"

    def test_all_five_metrics_present(self):
        df = micro_stats_table({"run_A": _micro_stats()})
        assert len(df.columns) == 5


# ---------------------------------------------------------------------------
# macro_stats_table
# ---------------------------------------------------------------------------


class TestMacroStatsTable:
    def test_returns_dataframe(self):
        df = macro_stats_table({"run_A": _macro_stats()})
        assert isinstance(df, pd.DataFrame)

    def test_empty_dict_returns_empty_df(self):
        df = macro_stats_table({})
        assert df.empty

    def test_index_is_run_codes(self):
        df = macro_stats_table({"run_A": _macro_stats(), "run_B": _macro_stats()})
        assert list(df.index) == ["run_A", "run_B"]

    def test_columns_match_metrics(self):
        df = macro_stats_table({"run_A": _macro_stats()})
        assert list(df.columns) == _MACRO_METRICS

    def test_mean_pm_std_format(self):
        df = macro_stats_table({"run_A": _macro_stats(mean=1.234, std=0.056)})
        val = df.loc["run_A", _MACRO_METRICS[0]]
        assert "1.234" in val
        assert "\u00b1" in val
        assert "0.056" in val

    def test_three_decimal_places(self):
        df = macro_stats_table({"run_A": _macro_stats(mean=1.0, std=0.5)})
        val = df.loc["run_A", _MACRO_METRICS[0]]
        assert "1.000" in val
        assert "0.500" in val

    def test_multiple_runs_different_values(self):
        df = macro_stats_table(
            {"run_X": _macro_stats(1.0, 0.1), "run_Y": _macro_stats(2.0, 0.2)}
        )
        assert "1.000" in df.loc["run_X", _MACRO_METRICS[0]]
        assert "2.000" in df.loc["run_Y", _MACRO_METRICS[0]]

    def test_all_metrics_present(self):
        df = macro_stats_table({"run_A": _macro_stats()})
        assert len(df.columns) == len(_MACRO_METRICS)


# ---------------------------------------------------------------------------
# deg_rate_table
# ---------------------------------------------------------------------------


class TestDegRateTable:
    def _make_stats(self, mean=0.1234, std=0.0056):
        return {ivl: (mean, std) for ivl in _DEFAULT_INTERVALS}

    def test_returns_dataframe(self):
        df = deg_rate_table({"run_A": self._make_stats()}, _DEFAULT_INTERVALS)
        assert isinstance(df, pd.DataFrame)

    def test_index_is_run_codes(self):
        df = deg_rate_table(
            {"run_A": self._make_stats(), "run_B": self._make_stats()},
            _DEFAULT_INTERVALS,
        )
        assert list(df.index) == ["run_A", "run_B"]

    def test_columns_match_intervals(self):
        df = deg_rate_table({"run_A": self._make_stats()}, _DEFAULT_INTERVALS)
        assert list(df.columns) == ["20% to 80%", "20% to 50%", "50% to 80%"]

    def test_column_order_matches_intervals(self):
        intervals = [(50, 80), (20, 50)]
        stats = {ivl: (0.1, 0.01) for ivl in intervals}
        df = deg_rate_table({"run_A": stats}, intervals)
        assert list(df.columns) == ["50% to 80%", "20% to 50%"]

    def test_four_decimal_places(self):
        stats = {(20, 80): (0.1234, 0.0056)}
        df = deg_rate_table({"run_A": stats}, [(20, 80)])
        val = df.loc["run_A", "20% to 80%"]
        assert "0.1234" in val
        assert "0.0056" in val

    def test_mean_pm_std_format(self):
        df = deg_rate_table({"run_A": self._make_stats()}, _DEFAULT_INTERVALS)
        assert "\u00b1" in df.loc["run_A", "20% to 80%"]

    def test_multiple_runs(self):
        intervals = [(20, 80)]
        df = deg_rate_table(
            {"run_X": {(20, 80): (0.1, 0.01)}, "run_Y": {(20, 80): (0.2, 0.02)}},
            intervals,
        )
        assert "0.1000" in df.loc["run_X", "20% to 80%"]
        assert "0.2000" in df.loc["run_Y", "20% to 80%"]

    def test_no_units_in_column_names(self):
        df = deg_rate_table({"run_A": self._make_stats()}, _DEFAULT_INTERVALS)
        for col in df.columns:
            assert "/min" not in col and "%" in col


# ---------------------------------------------------------------------------
# deg_time_table
# ---------------------------------------------------------------------------


class TestDegTimeTable:
    def _make_stats(self, mean=42.0, std=3.5):
        return {m: (mean, std) for m in _DEFAULT_MARKERS}

    def test_returns_dataframe(self):
        df = deg_time_table({"run_A": self._make_stats()}, _DEFAULT_MARKERS)
        assert isinstance(df, pd.DataFrame)

    def test_index_is_run_codes(self):
        df = deg_time_table(
            {"run_A": self._make_stats(), "run_B": self._make_stats()},
            _DEFAULT_MARKERS,
        )
        assert list(df.index) == ["run_A", "run_B"]

    def test_columns_match_markers(self):
        df = deg_time_table({"run_A": self._make_stats()}, _DEFAULT_MARKERS)
        assert list(df.columns) == ["5%", "20%", "50%", "80%", "100%"]

    def test_column_order_matches_markers(self):
        markers = [100, 50, 20]
        stats = {m: (10.0, 1.0) for m in markers}
        df = deg_time_table({"run_A": stats}, markers)
        assert list(df.columns) == ["100%", "50%", "20%"]

    def test_two_decimal_places(self):
        stats = {50: (42.0, 3.5)}
        df = deg_time_table({"run_A": stats}, [50])
        val = df.loc["run_A", "50%"]
        assert "42.00" in val
        assert "3.50" in val

    def test_mean_pm_std_format(self):
        df = deg_time_table({"run_A": self._make_stats()}, _DEFAULT_MARKERS)
        assert "\u00b1" in df.loc["run_A", "50%"]

    def test_multiple_runs(self):
        markers = [50]
        df = deg_time_table(
            {"run_X": {50: (40.0, 2.0)}, "run_Y": {50: (50.0, 4.0)}},
            markers,
        )
        assert "40.00" in df.loc["run_X", "50%"]
        assert "50.00" in df.loc["run_Y", "50%"]

    def test_no_units_in_column_names(self):
        df = deg_time_table({"run_A": self._make_stats()}, _DEFAULT_MARKERS)
        for col in df.columns:
            assert "min" not in col


# ---------------------------------------------------------------------------
# parameters_table
# ---------------------------------------------------------------------------


class TestParametersTable:
    _PARAM_SPECS = [
        ("pore_size", "macro", "microns", None),
        ("cols", "macro", None, None),
        ("bind_rate_tPA", "micro", None, None),
    ]
    _NATURAL_UNITS = {"pore_size": "microns", "cols": None, "bind_rate_tPA": None}

    def _make_values(self, value="1.234", add_names=()):
        result = {attr_name: value for attr_name, *_ in self._PARAM_SPECS}
        result.update({name: value for name in add_names})
        return result

    def test_returns_dataframe(self):
        df = parameters_table(
            {"run_A": self._make_values()},
            self._PARAM_SPECS, [],
            self._NATURAL_UNITS, ["run_A"],
        )
        assert isinstance(df, pd.DataFrame)

    def test_has_multiindex(self):
        df = parameters_table(
            {"run_A": self._make_values()},
            self._PARAM_SPECS, [],
            self._NATURAL_UNITS, ["run_A"],
        )
        assert isinstance(df.index, pd.MultiIndex)

    def test_multiindex_names(self):
        df = parameters_table(
            {"run_A": self._make_values()},
            self._PARAM_SPECS, [],
            self._NATURAL_UNITS, ["run_A"],
        )
        assert df.index.names == ["section", "parameter"]

    def test_macro_params_in_macroscale_section(self):
        df = parameters_table(
            {"run_A": self._make_values()},
            self._PARAM_SPECS, [],
            self._NATURAL_UNITS, ["run_A"],
        )
        sections = df.index.get_level_values("section")
        assert "Macroscale" in sections.values

    def test_micro_params_in_microscale_section(self):
        df = parameters_table(
            {"run_A": self._make_values()},
            self._PARAM_SPECS, [],
            self._NATURAL_UNITS, ["run_A"],
        )
        sections = df.index.get_level_values("section")
        assert "Microscale" in sections.values

    def test_add_names_in_additional_section(self):
        df = parameters_table(
            {"run_A": self._make_values(add_names=["extra_param"])},
            self._PARAM_SPECS, ["extra_param"],
            self._NATURAL_UNITS, ["run_A"],
        )
        sections = df.index.get_level_values("section")
        assert "Additional" in sections.values

    def test_no_add_names_no_additional_section(self):
        df = parameters_table(
            {"run_A": self._make_values()},
            self._PARAM_SPECS, [],
            self._NATURAL_UNITS, ["run_A"],
        )
        sections = df.index.get_level_values("section")
        assert "Additional" not in sections.values

    def test_run_codes_as_columns(self):
        df = parameters_table(
            {"run_A": self._make_values(), "run_B": self._make_values("2.0")},
            self._PARAM_SPECS, [],
            self._NATURAL_UNITS, ["run_A", "run_B"],
        )
        assert list(df.columns) == ["run_A", "run_B"]

    def test_column_order_follows_run_codes_arg(self):
        df = parameters_table(
            {"run_A": self._make_values(), "run_B": self._make_values()},
            self._PARAM_SPECS, [],
            self._NATURAL_UNITS, ["run_B", "run_A"],
        )
        assert list(df.columns) == ["run_B", "run_A"]

    def test_param_label_includes_units(self):
        df = parameters_table(
            {"run_A": self._make_values()},
            self._PARAM_SPECS, [],
            self._NATURAL_UNITS, ["run_A"],
        )
        labels = df.index.get_level_values("parameter").tolist()
        assert "pore_size (microns)" in labels

    def test_param_label_no_units_when_none(self):
        df = parameters_table(
            {"run_A": self._make_values()},
            self._PARAM_SPECS, [],
            self._NATURAL_UNITS, ["run_A"],
        )
        labels = df.index.get_level_values("parameter").tolist()
        assert "cols" in labels

    def test_add_name_appears_in_additional(self):
        df = parameters_table(
            {"run_A": self._make_values(add_names=["extra_param"])},
            self._PARAM_SPECS, ["extra_param"],
            self._NATURAL_UNITS, ["run_A"],
        )
        labels = df.loc["Additional"].index.tolist()
        assert "extra_param" in labels

    def test_missing_value_gives_na(self):
        values = self._make_values()
        del values["cols"]
        df = parameters_table(
            {"run_A": values},
            self._PARAM_SPECS, [],
            self._NATURAL_UNITS, ["run_A"],
        )
        assert df.loc[("Macroscale", "cols"), "run_A"] == "N/A"


# ---------------------------------------------------------------------------
# compare_stats_table — data_diff rendering
# ---------------------------------------------------------------------------


class TestCompareStatsTableDataDiff:
    def test_match_renders_ok(self):
        results = {
            "run_A": {
                "data_diff": {
                    "microscale_out/fiber_degraded": {
                        "status": "match",
                        "max_pct_diff": 0.0,
                        "location": None,
                    }
                }
            }
        }
        df = compare_stats_table(results)
        assert df.loc["run_A", "microscale_out/fiber_degraded"] == "OK"

    def test_diff_renders_pct_and_location(self):
        results = {
            "run_A": {
                "data_diff": {
                    "microscale_out/tpa_leaving_time": {
                        "status": "diff",
                        "max_pct_diff": 1.5,
                        "location": (3, 7),
                    }
                }
            }
        }
        df = compare_stats_table(results)
        cell = df.loc["run_A", "microscale_out/tpa_leaving_time"]
        assert "+1.50%" in cell
        assert "(3, 7)" in cell

    def test_shape_mismatch_renders_detail(self):
        results = {
            "run_A": {
                "data_diff": {
                    "microscale_out/snapshot_time": {
                        "status": "shape_mismatch",
                        "max_pct_diff": None,
                        "location": None,
                        "detail": "(100,) vs (101,)",
                    }
                }
            }
        }
        df = compare_stats_table(results)
        cell = df.loc["run_A", "microscale_out/snapshot_time"]
        assert "shapes" in cell
        assert "(100,) vs (101,)" in cell

    def test_missing_renders_side(self):
        results = {
            "run_A": {
                "data_diff": {
                    "microscale_out/extra": {
                        "status": "missing",
                        "max_pct_diff": None,
                        "location": None,
                        "detail": "not in run2",
                    }
                }
            }
        }
        df = compare_stats_table(results)
        cell = df.loc["run_A", "microscale_out/extra"]
        assert "missing" in cell
        assert "not in run2" in cell

    def test_mixed_format_still_works(self):
        """Result dict with only data_diff (no ks/pct_diff) should still render."""
        results = {
            "run_A": {
                "data_diff": {
                    "microscale_out/x": {"status": "match", "max_pct_diff": 0.0, "location": None},
                }
            }
        }
        df = compare_stats_table(results)
        assert list(df.columns) == ["microscale_out/x"]
