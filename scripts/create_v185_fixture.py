#!/usr/bin/env python3
"""Create the committed test fixture for v1.85.0 Fortran data.

v1.85.0 predates the one-directory-per-simulation convention, so every
macroscale dataset is a single top-level file with all simulations
concatenated.  A v1.85.0 run also spans two directories: the microscale output
lives with the microscale run, while the macroscale output and the
micro-to-macro input files live with the macroscale run.

Reads from::

    data/2023-02-01-2200/   microscale_out
    data/2023-02-02-2200/   macroscale_in, macroscale_out

and writes a truncated copy (2 simulations x 3 snapshots) to
``tests/fixtures/fortran_v185_sample/``.

The source directories are opened READ-ONLY; nothing is ever written back to
``data/``.

Usage::

    uv run python scripts/create_v185_fixture.py

The script is idempotent: re-running it overwrites the existing fixture.
"""

import json
import pathlib
import re
import shutil

import numpy as np

# ─── Paths ────────────────────────────────────────────────────────────────────

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
MICRO_SRC_DIR = REPO_ROOT / "data" / "2023-02-01-2200"
MACRO_SRC_DIR = REPO_ROOT / "data" / "2023-02-02-2200"
FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "fortran_v185_sample"

# ─── File codes ───────────────────────────────────────────────────────────────
# Note the macroscale output uses a *different* code from the microscale output
# and the macroscale input files.

MICRO_FILE_CODE = "_PLG2_tPA01_Q2"
MACRO_IN_FILE_CODE = "_PLG2_tPA01_Q2"
MACRO_OUT_FILE_CODE = "_PLG2_tPA01_along_Q2"

# ─── Truncation ───────────────────────────────────────────────────────────────

N_SIMULATIONS = 2  # simulations kept
N_SNAPSHOTS = 3  # snapshot rows per simulation (t=0 plus 2 saves)
N_MICRO_RUNS = 2000  # microscale runs kept, to keep the fixture small

#: Marks the start of each simulation's block in the combined macro log.
RUN_MARKER = re.compile(r"^\s*run number\s*=")

# ─── Source parameters (read-only) ────────────────────────────────────────────

with open(MACRO_SRC_DIR / "params.json") as f:
    _raw_params = json.load(f)
_mp = _raw_params["macro_params"]

TOTAL_EDGES = _mp["total_edges"]
TOTAL_MOLECULES = _mp["total_molecules"]

print(f"Microscale source : {MICRO_SRC_DIR}")
print(f"Macroscale source : {MACRO_SRC_DIR}")
print(f"Fixture directory : {FIXTURE_DIR}")
print(f"total_edges={TOTAL_EDGES}, total_molecules={TOTAL_MOLECULES}")

FIXTURE_DIR.mkdir(parents=True, exist_ok=True)

# ─── params.json ──────────────────────────────────────────────────────────────
# Kept in genuine v1.85.0 shape -- legacy key spellings, the top-level
# bookkeeping written by the old Experiment class, and a null micro_params --
# so the fixture exercises the real v1.85.0 parameter path.

print("\n[1/5] Writing params.json …")

params_json = {
    "experiment_code": "fortran_v185_sample",
    "data_filenames": _raw_params.get("data_filenames", {}),
    "micro_params": None,
    "macro_params": dict(_mp, total_trials=N_SIMULATIONS, number_of_saves=N_SNAPSHOTS),
}
with open(FIXTURE_DIR / "params.json", "w") as f:
    json.dump(params_json, f, indent=2)
print(f"  Written: params.json (total_trials={N_SIMULATIONS})")

# ─── Microscale output ────────────────────────────────────────────────────────

print(f"\n[2/5] Truncating microscale files to {N_MICRO_RUNS} runs …")

MICRO_DTYPES = {
    "firstPLi": np.float64,
    "lasttPA": np.int32,
    "lyscomplete": np.uint32,
    "lysis": np.float64,
    "PLi": np.int32,
    "tPA_time": np.float64,
    "tPAPLiunbd": np.int32,
    "tPAunbind": np.int32,
}
for name, dtype in MICRO_DTYPES.items():
    src = MICRO_SRC_DIR / f"{name}{MICRO_FILE_CODE}.dat"
    if not src.exists():
        print(f"  MISSING: {src.name}")
        continue
    np.fromfile(src, dtype=dtype)[:N_MICRO_RUNS].tofile(
        FIXTURE_DIR / f"{name}{MICRO_FILE_CODE}.dat"
    )
    print(f"  Written: {name}{MICRO_FILE_CODE}.dat")

# The micro log is parsed for parameters, so keep "runs=" consistent with the
# truncated binaries.
micro_log = (MICRO_SRC_DIR / f"micro{MICRO_FILE_CODE}.txt").read_text().splitlines()
micro_log = [
    re.sub(r"(runs=)\s*\d+", rf"\g<1>{N_MICRO_RUNS:>12}", line) for line in micro_log
]
(FIXTURE_DIR / f"micro{MICRO_FILE_CODE}.txt").write_text("\n".join(micro_log) + "\n")
print(f"  Written: micro{MICRO_FILE_CODE}.txt (runs={N_MICRO_RUNS})")

# ─── Macroscale input ─────────────────────────────────────────────────────────
# Small text files with fixed bin counts; copied verbatim.

print("\n[3/5] Copying macroscale input files …")

for name in ("tPAleave", "tsectPA", "lysismat", "lenlysisvect"):
    src = MACRO_SRC_DIR / f"{name}{MACRO_IN_FILE_CODE}.dat"
    if src.exists():
        shutil.copy2(src, FIXTURE_DIR / f"{name}{MACRO_IN_FILE_CODE}.dat")
        print(f"  Copied: {name}{MACRO_IN_FILE_CODE}.dat")
    else:
        print(f"  MISSING: {src.name}")

# ─── Macroscale output ────────────────────────────────────────────────────────

print(
    f"\n[4/5] Truncating macroscale output to {N_SIMULATIONS} simulations "
    f"x {N_SNAPSHOTS} snapshots …"
)


def _src(name):
    return MACRO_SRC_DIR / f"{name}{MACRO_OUT_FILE_CODE}.dat"


def _dst(name):
    return FIXTURE_DIR / f"{name}{MACRO_OUT_FILE_CODE}.dat"


# Nsave drives every snapshot-indexed split, so it is truncated first and the
# rest of the datasets are sliced to match.
nsave_src = np.fromfile(_src("Nsave"), dtype=np.int32)
offsets = np.concatenate([[0], np.cumsum(nsave_src + 1)])
print(f"  Source Nsave: {nsave_src} ({offsets[-1]} snapshots total)")

nsave_out = np.full(N_SIMULATIONS, N_SNAPSHOTS - 1, dtype=np.int32)
nsave_out.tofile(_dst("Nsave"))
print(f"  Written: Nsave = {nsave_out}")

# Row indices to keep: the first N_SNAPSHOTS rows of each kept simulation.
keep_rows = np.concatenate(
    [offsets[sim] + np.arange(N_SNAPSHOTS) for sim in range(N_SIMULATIONS)]
)

SNAPSHOT_DATASETS = {
    "tsave": (np.float64, None),
    "f_deg_time": (np.float64, TOTAL_EDGES),
    "deg": (np.float64, TOTAL_EDGES),
    "m_loc": (np.int32, TOTAL_MOLECULES),
    "m_bound": (np.int32, TOTAL_MOLECULES),
}
for name, (dtype, width) in SNAPSHOT_DATASETS.items():
    src = _src(name)
    if not src.exists():
        print(f"  MISSING: {src.name}")
        continue
    array = np.fromfile(src, dtype=dtype)
    if width is not None:
        array = array.reshape(-1, width)
    array[keep_rows].tofile(_dst(name))
    print(f"  Written: {name} shape {array[keep_rows].shape}")

# mfpt holds one row per simulation, so it needs no Nsave bookkeeping.
mfpt = np.fromfile(_src("mfpt"), dtype=np.float64).reshape(-1, TOTAL_MOLECULES)
mfpt[:N_SIMULATIONS].tofile(_dst("mfpt"))
print(f"  Written: mfpt shape {mfpt[:N_SIMULATIONS].shape}")

# The combined macro log: keep the shared header plus the first N_SIMULATIONS
# "run number=" blocks.
log_lines = (MACRO_SRC_DIR / f"macro{MACRO_OUT_FILE_CODE}.txt").read_text().splitlines()
starts = [i for i, line in enumerate(log_lines) if RUN_MARKER.match(line)]
end = starts[N_SIMULATIONS] if len(starts) > N_SIMULATIONS else len(log_lines)
(FIXTURE_DIR / f"macro{MACRO_OUT_FILE_CODE}.txt").write_text(
    "\n".join(log_lines[:end]) + "\n"
)
print(
    f"  Written: macro log, {end} lines "
    f"({starts[0]} header + {N_SIMULATIONS} run blocks)"
)

# ─── README ───────────────────────────────────────────────────────────────────

print("\n[5/5] Writing README.md …")

tsave_out = np.fromfile(_dst("tsave"), dtype=np.float64)
(FIXTURE_DIR / "README.md").write_text(f"""v1.85.0 Fortran data fixture (truncated)
=========================================

Generated from:
  - microscale : data/2023-02-01-2200/
  - macroscale : data/2023-02-02-2200/
Script: scripts/create_v185_fixture.py

Micro file code     : {MICRO_FILE_CODE}
Macro-in file code  : {MACRO_IN_FILE_CODE}
Macro-out file code : {MACRO_OUT_FILE_CODE}
total_edges         : {TOTAL_EDGES}
total_molecules     : {TOTAL_MOLECULES}
Simulations         : {N_SIMULATIONS}
Snapshots per sim   : {N_SNAPSHOTS} (t=0 plus {N_SNAPSHOTS - 1} saves)
Microscale runs     : {N_MICRO_RUNS}
tsave               : {list(tsave_out)}

v1.85.0 stores every simulation in ONE top-level file per dataset, so
f_deg_time has shape ({N_SIMULATIONS * N_SNAPSHOTS}, {TOTAL_EDGES}) and
Nsave is a vector of length {N_SIMULATIONS}, not a scalar.

Note that this fixture deliberately keeps the v1.85.0 quirks that the
converters exist to handle:

  - m_bind_t is absent (v1.85.0 never wrote it)
  - f_deg_time uses 0.0 for BOTH empty edges and unscheduled fibrin
  - deg (a v1.85.0-only degradation-state array) is present
  - params.json uses the legacy key spellings and a null micro_params
""")
print("  Written: README.md")
print("\nDone.")
