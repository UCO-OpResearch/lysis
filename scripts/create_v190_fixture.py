#!/usr/bin/env python3
"""Create the committed test fixture for v1.90.0 Fortran macroscale data.

Reads from data/2026-02-28-1907/ and writes a truncated 3-snapshot copy to
tests/fixtures/fortran_v190_sample/.  Run this script once to regenerate
the fixture if the source data changes.

Usage::

    python scripts/create_v190_fixture.py

The script is idempotent: re-running it overwrites the existing fixture.
"""

import json
import pathlib
import shutil

import numpy as np

# ─── Paths ────────────────────────────────────────────────────────────────────

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC_DIR = REPO_ROOT / "data" / "2026-02-28-1907"
FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "fortran_v190_sample"

# ─── File codes ───────────────────────────────────────────────────────────────

MICRO_FILE_CODE = "_PLG2_tPA01_Q4"
MACRO_FILE_CODE = "_TK-L_307"

# ─── How many snapshots to include in the fixture ─────────────────────────────

N_SNAPSHOTS = 3  # rows in f_deg_time (t=0, t1, t2)

# ─── Derived constants (read from source data) ────────────────────────────────

# Read params.json to get total_edges and total_molecules
with open(SRC_DIR / "params.json") as f:
    _raw_params = json.load(f)
_mp = _raw_params["macro_params"]

TOTAL_EDGES = _mp["total_edges"]
TOTAL_MOLECULES = _mp["total_molecules"]
ROWS = _mp["rows"]
COLS = _mp["cols"]
EMPTY_ROWS = _mp["empty_rows"]
PORE_SIZE = _mp["pore_size"]       # cm
DIFFUSION_COEFF = _mp["diffusion_coeff"]  # cm^2/s
FORCED_UNBIND = _mp["forced_unbind"]
MOVING_PROB = _mp["moving_probability"]
MACRO_SIMULATIONS = _mp.get("total_trials") or _mp.get("macro_simulations", 10)
MACRO_SEED = _mp.get("seed") or _mp.get("macro_seed", 0)
SAVE_INTERVAL = _mp.get("save_interval", 100)
MACRO_VERSION = _mp.get("macro_version", "diffuse_into_and_along")
DUPLICATE_FORTRAN = _mp.get("duplicate_fortran", False)
PROCESSING_LIBRARY = _mp.get("processing_library", "numpy")

print(f"Source directory: {SRC_DIR}")
print(f"Fixture directory: {FIXTURE_DIR}")
print(f"total_edges={TOTAL_EDGES}, total_molecules={TOTAL_MOLECULES}, "
      f"rows={ROWS}, cols={COLS}")


# ─── Create output directory ──────────────────────────────────────────────────

FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
(FIXTURE_DIR / "00").mkdir(exist_ok=True)

print("\n[1/4] Writing params.json …")

# Write a params.json with fields in the format expected by MacroParameters.
params_json = {
    "macro_params": {
        "pore_size": PORE_SIZE,
        "diffusion_coeff": DIFFUSION_COEFF,
        "forced_unbind": FORCED_UNBIND,
        "cols": COLS,
        "rows": ROWS,
        "empty_rows": EMPTY_ROWS,
        "total_molecules": TOTAL_MOLECULES,
        # total_edges is a dependent parameter; include it so that parse_shape()
        # can resolve the f_deg_time shape before _validate_fortran_params runs.
        "total_edges": TOTAL_EDGES,
        "moving_probability": MOVING_PROB,
        "macro_simulations": MACRO_SIMULATIONS,
        "total_time": 0,
        "macro_seed": MACRO_SEED,
        "save_interval": SAVE_INTERVAL,
        "macro_version": MACRO_VERSION,
        "macro_log_lvl": 30,
        "duplicate_fortran": DUPLICATE_FORTRAN,
        "processing_library": PROCESSING_LIBRARY,
    }
}
with open(FIXTURE_DIR / "params.json", "w") as f:
    json.dump(params_json, f, indent=2)
print("  Written: params.json")

# ─── Copy microscale output files ─────────────────────────────────────────────

print("\n[2/4] Copying microscale files …")

micro_files = [
    f"micro{MICRO_FILE_CODE}.txt",
    f"firstPLi{MICRO_FILE_CODE}.dat",
    f"lasttPA{MICRO_FILE_CODE}.dat",
    f"lyscomplete{MICRO_FILE_CODE}.dat",
    f"lysis{MICRO_FILE_CODE}.dat",
    f"PLi{MICRO_FILE_CODE}.dat",
    f"tPA_time{MICRO_FILE_CODE}.dat",
    f"tPAPLiunbd{MICRO_FILE_CODE}.dat",
    f"tPAunbind{MICRO_FILE_CODE}.dat",
]
for fname in micro_files:
    src = SRC_DIR / fname
    dst = FIXTURE_DIR / fname
    if src.exists():
        shutil.copy2(src, dst)
        print(f"  Copied: {fname}")
    else:
        print(f"  MISSING: {fname}")

# ─── Read simulation 00 source data ──────────────────────────────────────────

print("\n[3/4] Reading and truncating simulation 00 …")

sim_src = SRC_DIR / "00"

# Nsave and tsave
nsave_src = np.fromfile(
    sim_src / f"Nsave{MACRO_FILE_CODE}_00.dat", dtype=np.int32
)
print(f"  Source Nsave: {nsave_src[0]}")

tsave_src = np.fromfile(
    sim_src / f"tsave{MACRO_FILE_CODE}_00.dat", dtype=np.float64
)
print(f"  Source tsave shape: {tsave_src.shape} (total snapshots)")

# f_deg_time (Nsave+1 rows × total_edges cols)
f_deg_time_src = np.fromfile(
    sim_src / f"f_deg_time{MACRO_FILE_CODE}_00.dat", dtype=np.float64
).reshape(nsave_src[0] + 1, TOTAL_EDGES)
print(f"  Source f_deg_time shape: {f_deg_time_src.shape}")

# m_loc and m_bound (Nsave+1 rows × total_molecules cols)
m_loc_src = np.fromfile(
    sim_src / f"m_loc{MACRO_FILE_CODE}_00.dat", dtype=np.int32
).reshape(nsave_src[0] + 1, TOTAL_MOLECULES)

m_bound_src = np.fromfile(
    sim_src / f"m_bound{MACRO_FILE_CODE}_00.dat", dtype=np.int32
).reshape(nsave_src[0] + 1, TOTAL_MOLECULES)

# ─── Write truncated simulation 00 ───────────────────────────────────────────

sim_dst = FIXTURE_DIR / "00"

# Nsave for truncated data (N_SNAPSHOTS - 1, since t=0 is included)
nsave_trunc = np.int32(N_SNAPSHOTS - 1)
nsave_trunc.tofile(sim_dst / f"Nsave{MACRO_FILE_CODE}_00.dat")
print(f"  Written: Nsave = {nsave_trunc} (3 snapshots total)")

# tsave (first N_SNAPSHOTS elements)
tsave_trunc = tsave_src[:N_SNAPSHOTS]
tsave_trunc.tofile(sim_dst / f"tsave{MACRO_FILE_CODE}_00.dat")
print(f"  Written: tsave shape {tsave_trunc.shape}: {tsave_trunc}")

# f_deg_time (first N_SNAPSHOTS rows)
f_deg_time_trunc = f_deg_time_src[:N_SNAPSHOTS]
f_deg_time_trunc.tofile(sim_dst / f"f_deg_time{MACRO_FILE_CODE}_00.dat")
print(f"  Written: f_deg_time shape {f_deg_time_trunc.shape}")

# m_loc (first N_SNAPSHOTS rows)
m_loc_trunc = m_loc_src[:N_SNAPSHOTS]
m_loc_trunc.tofile(sim_dst / f"m_loc{MACRO_FILE_CODE}_00.dat")
print(f"  Written: m_loc shape {m_loc_trunc.shape}")

# m_bound (first N_SNAPSHOTS rows)
m_bound_trunc = m_bound_src[:N_SNAPSHOTS]
m_bound_trunc.tofile(sim_dst / f"m_bound{MACRO_FILE_CODE}_00.dat")
print(f"  Written: m_bound shape {m_bound_trunc.shape}")

# Copy mfpt unchanged (1D, no snapshot dimension)
mfpt_src = np.fromfile(
    sim_src / f"mfpt{MACRO_FILE_CODE}_00.dat", dtype=np.float64
)
mfpt_src.tofile(sim_dst / f"mfpt{MACRO_FILE_CODE}_00.dat")
print(f"  Written: mfpt shape {mfpt_src.shape}")

# m_bind_t (text CSV; read all, write filtered to first N_SNAPSHOTS save time)
import csv
cutoff_time = float(tsave_trunc[-1])
m_bind_t_src = []
with open(sim_src / f"m_bind_t{MACRO_FILE_CODE}_00.dat", "r") as f:
    reader = csv.reader(f)
    for row in reader:
        if row and float(row[0].strip()) <= cutoff_time:
            m_bind_t_src.append(row)
with open(sim_dst / f"m_bind_t{MACRO_FILE_CODE}_00.dat", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(m_bind_t_src)
print(f"  Written: m_bind_t ({len(m_bind_t_src)} events ≤ t={cutoff_time:.2f}s)")

# macro_log: copy unchanged
macro_log_name = f"macro{MACRO_FILE_CODE}_00.txt"
if (sim_src / macro_log_name).exists():
    shutil.copy2(sim_src / macro_log_name, sim_dst / macro_log_name)
    print(f"  Copied: {macro_log_name}")

# ─── Summary ──────────────────────────────────────────────────────────────────

print("\n[4/4] Writing README …")

readme = f"""\
v1.90.0 Fortran data fixture (truncated)
=========================================

Generated from: data/2026-02-28-1907/
Script: scripts/create_v190_fixture.py

Micro file code : {MICRO_FILE_CODE}
Macro file code : {MACRO_FILE_CODE}
total_edges     : {TOTAL_EDGES}
total_molecules : {TOTAL_MOLECULES}
Snapshots       : {N_SNAPSHOTS} (t=0 plus {N_SNAPSHOTS - 1} saves)
tsave           : {list(tsave_trunc)}

Only simulation 00 is included.  The f_deg_time binary has shape
({N_SNAPSHOTS}, {TOTAL_EDGES}).  All other macroscale snapshot arrays
are also truncated to {N_SNAPSHOTS} rows.
"""
(FIXTURE_DIR / "README.md").write_text(readme)
print("  Written: README.md")

print(f"\nDone.  Fixture at: {FIXTURE_DIR}")
print("\nKey constants for tests/dataio/test_real_data.py:")
print(f"  V190_MICRO_FILE_CODE = \"{MICRO_FILE_CODE}\"")
print(f"  V190_MACRO_FILE_CODE = \"{MACRO_FILE_CODE}\"")
print(f"  V190_TOTAL_EDGES = {TOTAL_EDGES}")
print(f"  V190_TOTAL_MOLECULES = {TOTAL_MOLECULES}")
