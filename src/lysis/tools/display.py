"""Common display utilities for CLI table output.

Provides shared helpers used by all table-producing CLI commands:

- :func:`md_table` — build a GitHub-flavored Markdown table string
- :func:`emit_markdown` — write Markdown to stdout or a file
- :func:`stats_df_to_rich` — render a runs-as-rows DataFrame as a Rich Table
- :func:`stats_df_to_markdown` — render a runs-as-rows DataFrame as Markdown
- :func:`params_df_to_rich` — render a parameters DataFrame as a Rich Table
- :func:`params_df_to_markdown` — render a parameters DataFrame as Markdown
- :func:`render_dataset_side_by_side` — page-through side-by-side diff of two arrays
"""

from __future__ import annotations

import pandas as pd


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------


def md_table(headers: list[str], rows: list[list[str]]) -> str:
    """Build a GitHub-flavored Markdown table string.

    :param headers: Column header strings.
    :type headers: list[str]
    :param rows: Data rows, each a list of cell strings.
    :type rows: list[list[str]]
    :return: GFM table as a string (no trailing newline).
    :rtype: str
    """
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join(lines)


def emit_markdown(md_text: str, markdown_out: str, console) -> None:
    """Write *md_text* to stdout or a file.

    :param md_text: Markdown content to write.
    :type md_text: str
    :param markdown_out: ``"-"`` to print to stdout, or a file path to write.
    :type markdown_out: str
    :param console: Rich Console used for status and error messages.
    """
    if markdown_out == "-":
        print(md_text)
    else:
        try:
            with open(markdown_out, "w", encoding="utf-8") as fh:
                fh.write(md_text)
                if not md_text.endswith("\n"):
                    fh.write("\n")
            console.print(f"[green]Markdown written to {markdown_out}[/green]")
        except OSError as e:
            console.print(f"[red]Error:[/red] Could not write to {markdown_out}: {e}")


# ---------------------------------------------------------------------------
# Stats DataFrame → Rich / Markdown
# (runs as rows, metrics as columns)
# ---------------------------------------------------------------------------


def stats_df_to_rich(
    df: pd.DataFrame,
    index_header: str,
    short_headers: dict[str, str] | None = None,
):
    """Render a stats DataFrame (runs as rows) as a Rich Table.

    :param df: DataFrame with run codes as index and formatted string values
        as cells.
    :type df: pandas.DataFrame
    :param index_header: Label for the index column (e.g. ``"Run"``).
    :type index_header: str
    :param short_headers: Optional mapping from full column name to an
        abbreviated display name.  Abbreviated names may contain ``\\n`` for
        Rich multi-line headers.  Defaults to using the full column names.
    :type short_headers: dict[str, str] or None
    :return: Configured Rich Table ready for ``console.print``.
    :rtype: rich.table.Table
    """
    from rich.table import Table

    table = Table(show_header=True, header_style="bold", show_lines=True)
    table.add_column(index_header, style="bold", no_wrap=True)
    for col in df.columns:
        display_name = (short_headers or {}).get(col, col)
        table.add_column(display_name, justify="right")

    for run_code, row in df.iterrows():
        table.add_row(str(run_code), *[str(v) for v in row])

    return table


def stats_df_to_markdown(
    df: pd.DataFrame,
    index_header: str,
    single_code: str | None = None,
) -> str:
    """Render a stats DataFrame (runs as rows) as a GFM Markdown table.

    In single-file mode (*single_code* given) produces a heading and a
    transposed ``Metric / Value`` table.  In directory mode produces a flat
    table with one row per Run.

    :param df: DataFrame with run codes as index and formatted string values.
    :type df: pandas.DataFrame
    :param index_header: Label for the run-code column in directory mode.
    :type index_header: str
    :param single_code: Run code for single-file mode, or ``None`` for
        directory mode.
    :type single_code: str or None
    :return: GFM Markdown table string.
    :rtype: str
    """
    if single_code is not None:
        row = df.loc[single_code]
        table_rows = [[col, str(val)] for col, val in row.items()]
        return "\n".join([
            f"## {single_code}",
            "",
            md_table(["Metric", "Value"], table_rows),
        ])
    else:
        headers = [index_header] + list(df.columns)
        rows = [[str(rc)] + [str(v) for v in row] for rc, row in df.iterrows()]
        return md_table(headers, rows)


# ---------------------------------------------------------------------------
# Parameters DataFrame → Rich / Markdown
# (parameters as MultiIndex rows, run codes as columns)
# ---------------------------------------------------------------------------


def params_df_to_rich(df: pd.DataFrame):
    """Render a parameters DataFrame as a Rich Table.

    The DataFrame must have a two-level :class:`pandas.MultiIndex` as its
    index, with level 0 being the section name (``"Macroscale"``,
    ``"Microscale"``, or ``"Additional"``) and level 1 being the parameter
    display label.  A section separator is inserted between each distinct
    section group.

    :param df: Parameters DataFrame as returned by
        :func:`~lysis.analysis.summary.parameters_table`.
    :type df: pandas.DataFrame
    :return: Configured Rich Table ready for ``console.print``.
    :rtype: rich.table.Table
    """
    from rich.table import Table

    table = Table(show_header=True, header_style="bold", show_lines=True)
    table.add_column("Parameter", style="bold", no_wrap=True)
    for col in df.columns:
        table.add_column(str(col), justify="right")

    current_section = None
    for (section, param_label), row in df.iterrows():
        if current_section is not None and section != current_section:
            table.add_section()
        current_section = section
        table.add_row(param_label, *[str(v) for v in row])

    return table


def params_df_to_markdown(
    df: pd.DataFrame,
    single_code: str | None = None,
) -> str:
    """Render a parameters DataFrame as GFM Markdown.

    Produces one ``### Section Parameters`` heading and table per section.

    :param df: Parameters DataFrame as returned by
        :func:`~lysis.analysis.summary.parameters_table`.
    :type df: pandas.DataFrame
    :param single_code: Run code to emit as a ``## heading`` in single-file
        mode, or ``None`` for directory mode.
    :type single_code: str or None
    :return: GFM Markdown string.
    :rtype: str
    """
    lines: list[str] = []
    if single_code is not None:
        lines += [f"## {single_code}", ""]

    run_codes = list(df.columns)

    for section in df.index.get_level_values(0).unique():
        section_df = df.loc[section]
        lines += [f"### {section} Parameters", ""]
        headers = ["Parameter"] + run_codes
        table_rows = [
            [str(param_label)] + [str(v) for v in row]
            for param_label, row in section_df.iterrows()
        ]
        lines.append(md_table(headers, table_rows))
        lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Side-by-side dataset diff
# ---------------------------------------------------------------------------


def _fmt_scalar(v) -> str:
    """Format a numpy scalar / Python value for a single-cell diff view.

    Floats use ``:.6g`` (with ``NaN`` for non-numeric); booleans render
    as ``True`` / ``False``; everything else falls back to ``str``.

    :param v: Scalar value.
    :return: Display string.
    :rtype: str
    """
    import numpy as np

    if isinstance(v, (np.bool_, bool)):
        return str(bool(v))
    if isinstance(v, (np.floating, float)):
        if np.isnan(v):
            return "NaN"
        return f"{v:.6g}"
    return str(v)


def render_dataset_side_by_side(
    arr1,
    arr2,
    label: str,
    console,
    file1_name: str = "File 1",
    file2_name: str = "File 2",
) -> None:
    """Render two arrays side-by-side in a Rich Table with differences highlighted.

    Differing rows get a red foreground; identical rows render in the
    default style.  For plain numeric/boolean arrays the table gets an
    additional ``% Diff`` column showing the symmetric percent
    difference.  For structured arrays (e.g. event logs) each struct
    field becomes a pair of columns ``<field> (1)`` / ``<field> (2)``,
    and diffs are highlighted per field.

    Shape mismatch is non-fatal: a warning is printed and only the
    first ``min(len1, len2)`` rows of the common leading dimension are
    shown.

    Output is piped through the Rich pager (``console.pager``); when
    *console* is a TTY this uses the system pager (typically ``less``),
    otherwise it prints directly.

    :param arr1: First array.
    :type arr1: numpy.ndarray
    :param arr2: Second array.
    :type arr2: numpy.ndarray
    :param label: Dataset label used as the table title.
    :type label: str
    :param console: Rich Console.
    :param file1_name: Column header suffix for *arr1*.
    :type file1_name: str
    :param file2_name: Column header suffix for *arr2*.
    :type file2_name: str
    """
    import numpy as np

    from rich.table import Table

    from lysis.analysis.compare import percent_difference

    shape1 = arr1.shape
    shape2 = arr2.shape
    if shape1 != shape2:
        console.print(
            f"[yellow]Note:[/yellow] shape mismatch {shape1} vs {shape2}; "
            f"showing min(len1, len2) leading-axis rows."
        )

    is_struct = arr1.dtype.names is not None

    if is_struct:
        fields = arr1.dtype.names
        table = Table(title=label, show_header=True, header_style="bold")
        table.add_column("Index", style="dim", no_wrap=True)
        for field in fields:
            table.add_column(f"{field}\n({file1_name})", justify="right")
            table.add_column(f"{field}\n({file2_name})", justify="right")

        n = min(len(arr1), len(arr2))
        for i in range(n):
            cells = [str(i)]
            for field in fields:
                v1 = arr1[field][i]
                v2 = arr2[field][i]
                s1 = _fmt_scalar(v1)
                s2 = _fmt_scalar(v2)
                if v1 == v2:
                    cells.extend([s1, s2])
                else:
                    cells.extend([f"[bold red]{s1}[/]", f"[bold red]{s2}[/]"])
            table.add_row(*cells)
    else:
        flat1 = arr1.ravel()
        flat2 = arr2.ravel()
        # Use arr1's shape for index rendering; if shapes differ, callers
        # were already warned above.
        shape = shape1
        multi_dim = len(shape) > 1

        table = Table(title=label, show_header=True, header_style="bold")
        table.add_column("Index", style="dim", no_wrap=True)
        table.add_column(file1_name, justify="right")
        table.add_column(file2_name, justify="right")
        table.add_column("% Diff", justify="right")

        n = min(flat1.size, flat2.size)
        for i in range(n):
            if multi_dim:
                idx_str = str(np.unravel_index(i, shape))
            else:
                idx_str = str(i)
            v1 = flat1[i]
            v2 = flat2[i]
            s1 = _fmt_scalar(v1)
            s2 = _fmt_scalar(v2)
            if v1 == v2:
                table.add_row(idx_str, s1, s2, "—")
            else:
                try:
                    pct = percent_difference(float(v1), float(v2))
                except (TypeError, ValueError):
                    pct = float("nan")
                pct_str = "NaN" if pct != pct else f"{pct:+.2f}%"
                table.add_row(
                    f"[yellow]{idx_str}[/]",
                    f"[bold red]{s1}[/]",
                    f"[bold red]{s2}[/]",
                    f"[bold red]{pct_str}[/]",
                )

    with console.pager(styles=True):
        console.print(table)
