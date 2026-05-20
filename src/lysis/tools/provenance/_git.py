"""Shared git helpers used by both provenance submodules.

Kept private so the package's public surface stays focused on the
two ``gather_*`` and ``verify_*`` entry points.
"""

import subprocess
from pathlib import Path


def _package_repo_root() -> Path | None:
    """Locate the lysis source repo root (the directory containing ``.git``)."""
    import lysis

    here = Path(lysis.__file__).resolve().parent
    for candidate in (here, *here.parents):
        if (candidate / ".git").exists():
            return candidate
    return None


def _git(args: list[str], repo_root: Path) -> str | None:
    """Run a read-only git command in the repo root. Return stripped stdout or ``None``."""
    try:
        result = subprocess.run(
            ["git", "-C", str(repo_root), *args],
            capture_output=True,
            text=True,
            check=False,
        )
    except (FileNotFoundError, OSError):
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip()
