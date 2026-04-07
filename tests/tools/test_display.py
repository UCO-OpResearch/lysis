"""Tests for ``lysis.tools.display`` — shared display utilities."""

from unittest.mock import MagicMock

import pandas as pd
import pytest

from lysis.tools.display import (
    emit_markdown,
    md_table,
    params_df_to_markdown,
    params_df_to_rich,
    stats_df_to_markdown,
    stats_df_to_rich,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _simple_stats_df():
    """A two-run, two-metric stats DataFrame with pre-formatted string values."""
    data = {
        "run_X": {"Metric A": "1.000 \u00b1 0.100", "Metric B": "2.000 \u00b1 0.200"},
        "run_Y": {"Metric A": "3.000 \u00b1 0.300", "Metric B": "4.000 \u00b1 0.400"},
    }
    return pd.DataFrame.from_dict(data, orient="index")


def _single_stats_df():
    """A one-run stats DataFrame (for single-file mode tests)."""
    data = {
        "run_A": {"Metric A": "1.000 \u00b1 0.100", "Metric B": "2.000 \u00b1 0.200"}
    }
    return pd.DataFrame.from_dict(data, orient="index")


def _simple_params_df():
    """A parameters DataFrame with Macroscale and Microscale sections."""
    index = pd.MultiIndex.from_tuples(
        [
            ("Macroscale", "pore_size (microns)"),
            ("Microscale", "bind_rate_tPA"),
        ],
        names=["section", "parameter"],
    )
    return pd.DataFrame(
        {"run_A": ["1.234", "0.001"], "run_B": ["2.345", "0.002"]},
        index=index,
    )


# ---------------------------------------------------------------------------
# md_table — consolidated from 5 duplicated TestMdTable* classes in CLI tests
# ---------------------------------------------------------------------------


class TestMdTable:
    def test_header_pipe_format(self):
        out = md_table(["Parameter", "Run1"], [["pore_size", "1.0135"]])
        assert out.splitlines()[0] == "| Parameter | Run1 |"

    def test_separator_row(self):
        out = md_table(["A", "B", "C"], [])
        assert out.splitlines()[1] == "| --- | --- | --- |"

    def test_data_row_present(self):
        out = md_table(["P", "V"], [["pore_size (microns)", "1.0135"]])
        assert "| pore_size (microns) | 1.0135 |" in out

    def test_multiple_rows(self):
        out = md_table(["X"], [["r1"], ["r2"], ["r3"]])
        assert out.count("| r") == 3

    def test_empty_rows_only_two_lines(self):
        out = md_table(["A", "B"], [])
        assert len(out.splitlines()) == 2

    def test_single_column(self):
        out = md_table(["Only"], [["value"]])
        assert out.splitlines()[0] == "| Only |"
        assert "| value |" in out

    def test_many_columns(self):
        headers = ["A", "B", "C", "D"]
        out = md_table(headers, [["1", "2", "3", "4"]])
        assert out.splitlines()[0] == "| A | B | C | D |"


# ---------------------------------------------------------------------------
# emit_markdown
# ---------------------------------------------------------------------------


class TestEmitMarkdown:
    def test_stdout_prints_content(self, capsys):
        emit_markdown("## Hello\n\nContent", "-", MagicMock())
        captured = capsys.readouterr()
        assert "## Hello" in captured.out
        assert "Content" in captured.out

    def test_file_is_created(self, tmp_path):
        out_file = tmp_path / "test.md"
        emit_markdown("# Test", str(out_file), MagicMock())
        assert out_file.exists()

    def test_file_has_correct_content(self, tmp_path):
        out_file = tmp_path / "test.md"
        emit_markdown("# Test Content", str(out_file), MagicMock())
        assert "# Test Content" in out_file.read_text()

    def test_file_always_ends_with_newline(self, tmp_path):
        out_file = tmp_path / "test.md"
        emit_markdown("no trailing newline", str(out_file), MagicMock())
        assert out_file.read_text().endswith("\n")

    def test_file_with_trailing_newline_not_doubled(self, tmp_path):
        out_file = tmp_path / "test.md"
        emit_markdown("content\n", str(out_file), MagicMock())
        content = out_file.read_text()
        assert not content.endswith("\n\n")

    def test_file_prints_confirmation(self, tmp_path):
        console = MagicMock()
        out_file = tmp_path / "test.md"
        emit_markdown("# Test", str(out_file), console)
        console.print.assert_called_once()
        assert "Markdown written to" in console.print.call_args[0][0]

    def test_bad_path_prints_error_not_raises(self):
        console = MagicMock()
        emit_markdown("# Test", "/no_such_dir/out.md", console)
        # Should not raise; should print an error message
        console.print.assert_called_once()
        assert "Error" in console.print.call_args[0][0]

    def test_stdout_does_not_call_console_print(self, capsys):
        console = MagicMock()
        emit_markdown("# Test", "-", console)
        console.print.assert_not_called()


# ---------------------------------------------------------------------------
# stats_df_to_markdown
# ---------------------------------------------------------------------------


class TestStatsDfToMarkdown:
    # --- single-file mode ---

    def test_single_file_starts_with_heading(self):
        df = _single_stats_df()
        md = stats_df_to_markdown(df, "Run", single_code="run_A")
        assert md.startswith("## run_A")

    def test_single_file_has_metric_value_header(self):
        df = _single_stats_df()
        md = stats_df_to_markdown(df, "Run", single_code="run_A")
        assert "| Metric | Value |" in md

    def test_single_file_two_column_table(self):
        df = _single_stats_df()
        md = stats_df_to_markdown(df, "Run", single_code="run_A")
        header_line = next(l for l in md.splitlines() if "Metric" in l)
        assert header_line.count("|") == 3  # | Metric | Value |

    def test_single_file_metric_names_as_rows(self):
        df = _single_stats_df()
        md = stats_df_to_markdown(df, "Run", single_code="run_A")
        assert "Metric A" in md
        assert "Metric B" in md

    def test_single_file_values_present(self):
        df = _single_stats_df()
        md = stats_df_to_markdown(df, "Run", single_code="run_A")
        assert "1.000" in md

    # --- directory mode ---

    def test_directory_run_codes_as_rows(self):
        df = _simple_stats_df()
        md = stats_df_to_markdown(df, "Run")
        assert "run_X" in md
        assert "run_Y" in md

    def test_directory_header_starts_with_index_header(self):
        df = _simple_stats_df()
        md = stats_df_to_markdown(df, "Run")
        assert md.splitlines()[0].startswith("| Run |")

    def test_directory_all_columns_in_header(self):
        df = _simple_stats_df()
        md = stats_df_to_markdown(df, "Run")
        header = md.splitlines()[0]
        assert "Metric A" in header
        assert "Metric B" in header

    def test_directory_values_present(self):
        df = _simple_stats_df()
        md = stats_df_to_markdown(df, "Run")
        assert "1.000" in md
        assert "3.000" in md

    def test_directory_no_heading_hash(self):
        df = _simple_stats_df()
        md = stats_df_to_markdown(df, "Run")
        assert not md.startswith("#")

    def test_custom_index_header(self):
        df = _simple_stats_df()
        md = stats_df_to_markdown(df, "Experiment")
        assert md.splitlines()[0].startswith("| Experiment |")


# ---------------------------------------------------------------------------
# params_df_to_markdown
# ---------------------------------------------------------------------------


class TestParamsDfToMarkdown:
    def test_macroscale_section_header(self):
        md = params_df_to_markdown(_simple_params_df())
        assert "### Macroscale Parameters" in md

    def test_microscale_section_header(self):
        md = params_df_to_markdown(_simple_params_df())
        assert "### Microscale Parameters" in md

    def test_macro_param_in_macro_section(self):
        md = params_df_to_markdown(_simple_params_df())
        macro_section = md.split("### Macroscale Parameters")[1].split("### Microscale")[0]
        assert "pore_size (microns)" in macro_section

    def test_micro_param_in_micro_section(self):
        md = params_df_to_markdown(_simple_params_df())
        micro_section = md.split("### Microscale Parameters")[1]
        assert "bind_rate_tPA" in micro_section

    def test_run_codes_appear_as_columns(self):
        md = params_df_to_markdown(_simple_params_df())
        assert "run_A" in md
        assert "run_B" in md

    def test_multiple_run_codes_in_same_header_row(self):
        md = params_df_to_markdown(_simple_params_df())
        header_lines = [l for l in md.splitlines() if "run_A" in l]
        assert any("run_B" in l for l in header_lines)

    def test_single_file_heading(self):
        md = params_df_to_markdown(_simple_params_df(), single_code="run_A")
        assert md.startswith("## run_A")

    def test_no_add_names_no_additional_section(self):
        md = params_df_to_markdown(_simple_params_df())
        assert "### Additional Parameters" not in md

    def test_additional_section_when_present(self):
        index = pd.MultiIndex.from_tuples(
            [
                ("Macroscale", "pore_size (microns)"),
                ("Additional", "extra_param"),
            ],
            names=["section", "parameter"],
        )
        df = pd.DataFrame({"run_A": ["1.234", "0.999"]}, index=index)
        md = params_df_to_markdown(df)
        assert "### Additional Parameters" in md
        assert "extra_param" in md

    def test_drop_removes_param(self):
        # Build a DataFrame missing pore_size (simulates --drop)
        index = pd.MultiIndex.from_tuples(
            [("Microscale", "bind_rate_tPA")],
            names=["section", "parameter"],
        )
        df = pd.DataFrame({"run_A": ["0.001"]}, index=index)
        md = params_df_to_markdown(df)
        assert "pore_size" not in md

    def test_single_file_mode_has_sections(self):
        md = params_df_to_markdown(_simple_params_df(), single_code="run_A")
        assert "### Macroscale Parameters" in md
        assert "### Microscale Parameters" in md


# ---------------------------------------------------------------------------
# stats_df_to_rich
# ---------------------------------------------------------------------------


class TestStatsDfToRich:
    def test_returns_rich_table(self):
        from rich.table import Table

        table = stats_df_to_rich(_simple_stats_df(), "Run")
        assert isinstance(table, Table)

    def test_column_count(self):
        # 1 index column + 2 metric columns
        table = stats_df_to_rich(_simple_stats_df(), "Run")
        assert len(table.columns) == 3

    def test_index_header_as_first_column(self):
        table = stats_df_to_rich(_simple_stats_df(), "MyRuns")
        assert table.columns[0].header == "MyRuns"

    def test_metric_columns_present(self):
        table = stats_df_to_rich(_simple_stats_df(), "Run")
        col_headers = [col.header for col in table.columns]
        assert "Metric A" in col_headers
        assert "Metric B" in col_headers

    def test_short_headers_override_column_names(self):
        short = {"Metric A": "A\n(short)", "Metric B": "B\n(short)"}
        table = stats_df_to_rich(_simple_stats_df(), "Run", short_headers=short)
        col_headers = [col.header for col in table.columns]
        assert "A\n(short)" in col_headers
        assert "Metric A" not in col_headers

    def test_row_count_matches_dataframe(self):
        table = stats_df_to_rich(_simple_stats_df(), "Run")
        assert table.row_count == 2

    def test_no_short_headers_uses_full_names(self):
        table = stats_df_to_rich(_simple_stats_df(), "Run", short_headers=None)
        col_headers = [col.header for col in table.columns]
        assert "Metric A" in col_headers


# ---------------------------------------------------------------------------
# params_df_to_rich
# ---------------------------------------------------------------------------


class TestParamsDfToRich:
    def test_returns_rich_table(self):
        from rich.table import Table

        table = params_df_to_rich(_simple_params_df())
        assert isinstance(table, Table)

    def test_first_column_is_parameter(self):
        table = params_df_to_rich(_simple_params_df())
        assert table.columns[0].header == "Parameter"

    def test_run_codes_as_additional_columns(self):
        table = params_df_to_rich(_simple_params_df())
        col_headers = [col.header for col in table.columns]
        assert "run_A" in col_headers
        assert "run_B" in col_headers

    def test_row_count_matches_parameter_count(self):
        table = params_df_to_rich(_simple_params_df())
        assert table.row_count == 2  # pore_size + bind_rate_tPA

    def test_column_count(self):
        # 1 Parameter column + 2 run columns
        table = params_df_to_rich(_simple_params_df())
        assert len(table.columns) == 3
