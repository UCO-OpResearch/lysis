"""Regenerate the lysis internal-import network diagram.

Walks the Python files under ``src/lysis/``, extracts every internal import
(absolute or relative; including ``TYPE_CHECKING`` and ``try``/``except``
guarded imports), and writes ``docs/planning/lysis_internal_imports.dot``.
If the ``dot`` binary is on ``$PATH`` it also renders ``.svg`` and ``.png``
alongside the DOT source.

:Usage:
    python scripts/regen_import_graph.py
"""

from __future__ import annotations

import ast
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
PKG_ROOT = REPO_ROOT / "src" / "lysis"
PKG_NAME = "lysis"
OUT_DIR = REPO_ROOT / "docs" / "planning"
OUT_DOT = OUT_DIR / "lysis_internal_imports.dot"

# (cluster_key, label, fill, border). Order controls rendering order.
CLUSTERS = [
    ("top",       "lysis (top-level)", "#f5f5f5", "#aaaaaa"),
    ("config",    "lysis.config",      "#e8f0fe", "#4a6fa5"),
    ("dataio",    "lysis.dataio",      "#fff4e6", "#c98a4b"),
    ("geometry",  "lysis.geometry",    "#e8f5e9", "#4a8a4f"),
    ("execution", "lysis.execution",   "#fde8ea", "#b14a53"),
    ("analysis",  "lysis.analysis",    "#f3e8fd", "#7a4ba8"),
    ("tools",     "lysis.tools",       "#e8f7f8", "#3f8b90"),
    ("cli",       "lysis.cli",         "#fef6db", "#b69a2f"),
]
SUBPACKAGES = {key for key, *_ in CLUSTERS if key != "top"}


@dataclass(frozen=True)
class Edge:
    src: str
    dst: str
    kind: str          # "regular" | "type_checking" | "optional" | "parent"
    label: str = ""    # free-form label shown on the edge


def module_name_from_path(path: Path) -> str:
    rel = path.relative_to(PKG_ROOT).with_suffix("")
    parts = list(rel.parts)
    if parts and parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join((PKG_NAME, *parts)) if parts else PKG_NAME


def resolve_relative(package: str, level: int, name: str | None) -> str:
    """Mirror ``importlib._bootstrap._resolve_name`` for a known ``package``."""
    bits = package.rsplit(".", level - 1)
    if len(bits) < level:
        return ""  # attempted relative import beyond top-level
    base = bits[0]
    return f"{base}.{name}" if name else base


class ImportCollector(ast.NodeVisitor):
    def __init__(self, module_name: str, is_init: bool, modules: set[str]) -> None:
        self.module_name = module_name
        if is_init:
            self.package = module_name
        else:
            self.package = module_name.rsplit(".", 1)[0] if "." in module_name else ""
        self.modules = modules
        self.edges: list[Edge] = []
        self._tc_depth = 0
        self._try_stack: list[str] = []  # reason (may be "") per active try block

    @staticmethod
    def _is_type_checking(test: ast.expr) -> bool:
        if isinstance(test, ast.Name) and test.id == "TYPE_CHECKING":
            return True
        if isinstance(test, ast.Attribute) and test.attr == "TYPE_CHECKING":
            return True
        return False

    @staticmethod
    def _handles_import_error(handlers: list[ast.ExceptHandler]) -> bool:
        for h in handlers:
            t = h.type
            if t is None:
                return True
            names: list[ast.expr] = [t]
            if isinstance(t, ast.Tuple):
                names = list(t.elts)
            for n in names:
                if isinstance(n, ast.Name) and n.id in {
                    "ImportError", "ModuleNotFoundError",
                    "Exception", "BaseException",
                }:
                    return True
        return False

    @staticmethod
    def _first_external_import(body: list[ast.stmt]) -> str:
        for stmt in body:
            if isinstance(stmt, ast.Import):
                for a in stmt.names:
                    top = a.name.split(".", 1)[0]
                    if top != PKG_NAME:
                        return top
            elif isinstance(stmt, ast.ImportFrom) and stmt.level == 0 and stmt.module:
                top = stmt.module.split(".", 1)[0]
                if top != PKG_NAME:
                    return top
        return ""

    def visit_If(self, node: ast.If) -> None:
        if self._is_type_checking(node.test):
            self._tc_depth += 1
            for c in node.body:
                self.visit(c)
            self._tc_depth -= 1
            for c in node.orelse:
                self.visit(c)
        else:
            self.generic_visit(node)

    def visit_Try(self, node: ast.Try) -> None:
        if self._handles_import_error(node.handlers):
            reason = self._first_external_import(node.body)
            self._try_stack.append(reason)
            for c in node.body:
                self.visit(c)
            self._try_stack.pop()
            for h in node.handlers:
                for c in h.body:
                    self.visit(c)
            for c in node.orelse + node.finalbody:
                self.visit(c)
        else:
            self.generic_visit(node)

    def _classify(self, target: str) -> tuple[str, str]:
        if self._tc_depth:
            return "type_checking", "TYPE_CHECKING"
        if self._try_stack:
            reason = self._try_stack[-1]
            return "optional", f"optional ({reason})" if reason else "optional"
        parent = self.module_name.rsplit(".", 1)[0] if "." in self.module_name else ""
        if target == parent and target != self.module_name:
            return "parent", ""
        return "regular", ""

    def _emit(self, target: str) -> None:
        if target != PKG_NAME and not target.startswith(PKG_NAME + "."):
            return
        if target == self.module_name:
            return
        kind, label = self._classify(target)
        self.edges.append(Edge(self.module_name, target, kind, label))

    def visit_Import(self, node: ast.Import) -> None:
        for a in node.names:
            if a.name == PKG_NAME or a.name.startswith(PKG_NAME + "."):
                self._emit(a.name)

    def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
        if node.level:
            base = resolve_relative(self.package, node.level, node.module)
        else:
            base = node.module or ""
        if base != PKG_NAME and not base.startswith(PKG_NAME + "."):
            return
        targets: set[str] = set()
        for a in node.names:
            if a.name == "*":
                targets.add(base)
                continue
            candidate = f"{base}.{a.name}"
            if candidate in self.modules:
                targets.add(candidate)
            else:
                targets.add(base)
        for t in targets:
            self._emit(t)


def collect_edges() -> tuple[list[str], list[Edge]]:
    py_files = sorted(PKG_ROOT.rglob("*.py"))
    modules = {module_name_from_path(p) for p in py_files}
    edges: set[Edge] = set()
    for path in py_files:
        module_name = module_name_from_path(path)
        is_init = path.name == "__init__.py"
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        collector = ImportCollector(module_name, is_init, modules)
        collector.visit(tree)
        edges.update(collector.edges)
    return sorted(modules), sorted(edges, key=lambda e: (e.src, e.dst, e.kind))


def subpackage(module: str) -> str:
    if module == PKG_NAME:
        return "top"
    rest = module[len(PKG_NAME) + 1:]
    head = rest.split(".", 1)[0]
    return head if head in SUBPACKAGES else "top"


def short_label(module: str, cluster_key: str) -> str:
    if cluster_key == "top":
        return "__init__" if module == PKG_NAME else module[len(PKG_NAME) + 1:]
    prefix = f"{PKG_NAME}.{cluster_key}"
    return "__init__" if module == prefix else module[len(prefix) + 1:]


def edge_attrs(edge: Edge) -> str:
    if edge.kind == "type_checking":
        return 'style=dashed, label="TYPE_CHECKING"'
    if edge.kind == "optional":
        base = 'style=dashed, color="#b14a53"'
        return f'{base}, label="{edge.label}"' if edge.label else base
    if edge.kind == "parent":
        return 'style=dotted, color="#999999"'
    return 'style=solid'


def render_dot(modules: list[str], edges: list[Edge]) -> str:
    by_cluster: dict[str, list[str]] = {key: [] for key, *_ in CLUSTERS}
    for m in modules:
        by_cluster[subpackage(m)].append(m)

    def sort_key(m: str, cluster_key: str) -> tuple[int, str]:
        # __init__ first, then alphabetical
        return (0 if short_label(m, cluster_key) == "__init__" else 1, m)

    lines: list[str] = []
    lines.append(
        "// Internal import graph of the lysis Python package.\n"
        "// Auto-generated by scripts/regen_import_graph.py — do not edit by hand.\n"
        "//\n"
        "// Edge styles:\n"
        "//   solid       = regular runtime import\n"
        "//   dashed      = TYPE_CHECKING-only import\n"
        "//   dotted      = import of parent package (e.g. subcommand -> click group)\n"
        "//   dashed red  = conditional import guarded by try/except (optional dep)\n"
    )
    lines.append("digraph lysis_imports {")
    lines.append("    rankdir=LR;")
    lines.append("    newrank=true;")
    lines.append('    node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=10];')
    lines.append('    edge [fontname="Helvetica", fontsize=8, color="#555555"];')
    lines.append('    graph [fontname="Helvetica", fontsize=11, labeljust=l];')
    lines.append("")

    for cluster_key, cluster_label, fill, border in CLUSTERS:
        members = sorted(by_cluster[cluster_key], key=lambda m: sort_key(m, cluster_key))
        if not members:
            continue
        lines.append(f"    subgraph cluster_{cluster_key} {{")
        lines.append(f'        label="{cluster_label}";')
        lines.append(f'        style="rounded,filled"; fillcolor="{fill}"; color="{border}";')
        for m in members:
            lines.append(f'        "{m}" [label="{short_label(m, cluster_key)}", fillcolor="#ffffff"];')
        lines.append("    }")
        lines.append("")

    for e in edges:
        lines.append(f'    "{e.src}" -> "{e.dst}" [{edge_attrs(e)}];')

    lines.append("}")
    return "\n".join(lines) + "\n"


def main() -> int:
    modules, edges = collect_edges()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DOT.write_text(render_dot(modules, edges), encoding="utf-8")
    print(f"wrote {OUT_DOT.relative_to(REPO_ROOT)} "
          f"({len(modules)} modules, {len(edges)} edges)")

    dot = shutil.which("dot")
    if dot is None:
        print("note: `dot` not found on PATH — skipping SVG/PNG render", file=sys.stderr)
        return 0
    for fmt, extra in [("svg", []), ("png", ["-Gdpi=120"])]:
        out = OUT_DOT.with_suffix(f".{fmt}")
        subprocess.run(
            [dot, f"-T{fmt}", *extra, str(OUT_DOT), "-o", str(out)],
            check=True,
        )
        print(f"wrote {out.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
