"""Create truncated test fixture from real Fortran simulation data.

Reads data from data/2026-02-18-1723/ and writes a truncated subset to
tests/fixtures/fortran_sample/. The truncated data exercises all code paths
but is small enough to commit to git (~4.4 MB).

Truncation rules
----------------
- **Microscale**: All files copied as-is (full 50,000 values needed for
  generate_macroscale_in() validation).
- **Macroscale input**: All files copied as-is (aggregated data, not
  per-simulation).
- **Macroscale output**: Only simulation 00; truncated to 3 snapshots:
  - m_loc, m_bound: first 3 * total_molecules int32 values each
  - m_bind_t: all rows with timestamp < tsave[2]
  - f_deg_list: all rows with timestamp < tsave[2]
  - tsave: first 3 values
  - Nsave: scalar int32 = 2
  - mfpt: keep as-is (one value per molecule)
  - macro_log, params.json: copy as-is

Usage
-----
Run from the repository root::

    python scripts/create_test_fixture.py
"""

import shutil

import numpy as np
from pathlib import Path


SRC = Path("data/2026-02-18-1723")
DST = Path("tests/fixtures/fortran_sample")
MICRO_CODE = "_PLG2_tPA01_TB-xiii"
MACRO_CODE = "_TB-xiii__21_105"
TOTAL_MOLECULES = 21105
N_SNAPSHOTS = 3


def main():
    if not SRC.exists():
        raise FileNotFoundError(f"Source data not found: {SRC}")

    DST.mkdir(parents=True, exist_ok=True)
    (DST / "00").mkdir(exist_ok=True)

    # ── Microscale: copy everything as-is ──────────────────────────────
    # Log file
    shutil.copy2(SRC / f"micro{MICRO_CODE}.txt", DST / f"micro{MICRO_CODE}.txt")

    # Binary float64 datasets (full 50,000 values)
    for name in ["firstPLi", "lysis", "tPA_time"]:
        shutil.copy2(
            SRC / f"{name}{MICRO_CODE}.dat", DST / f"{name}{MICRO_CODE}.dat"
        )

    # Binary int32 datasets (full 50,000 values)
    for name in ["lasttPA", "lyscomplete", "PLi", "tPAPLiunbd", "tPAunbind"]:
        shutil.copy2(
            SRC / f"{name}{MICRO_CODE}.dat", DST / f"{name}{MICRO_CODE}.dat"
        )

    # ── Macroscale input: copy as-is ───────────────────────────────────
    for name in ["tPAleave", "tsectPA", "lysismat", "lenlysisvect", "neighbors"]:
        shutil.copy2(
            SRC / f"{name}{MICRO_CODE}.dat", DST / f"{name}{MICRO_CODE}.dat"
        )

    # ── Macroscale output: truncate to 1 sim, 3 snapshots ─────────────
    # params.json
    shutil.copy2(SRC / "params.json", DST / "params.json")

    # Macro log
    shutil.copy2(
        SRC / f"00/macro{MACRO_CODE}_00.txt",
        DST / f"00/macro{MACRO_CODE}_00.txt",
    )

    # tsave: first N_SNAPSHOTS values
    tsave = np.fromfile(SRC / f"00/tsave{MACRO_CODE}_00.dat", dtype=np.float64)
    tsave_trunc = tsave[:N_SNAPSHOTS]
    tsave_trunc.tofile(DST / f"00/tsave{MACRO_CODE}_00.dat")
    cutoff_time = tsave_trunc[-1]

    # Nsave: scalar int32 = N_SNAPSHOTS - 1 (number of save intervals)
    np.int32(N_SNAPSHOTS - 1).tofile(DST / f"00/Nsave{MACRO_CODE}_00.dat")

    # m_loc: first N_SNAPSHOTS * TOTAL_MOLECULES int32 values
    m_loc = np.fromfile(SRC / f"00/m_loc{MACRO_CODE}_00.dat", dtype=np.int32)
    m_loc[: N_SNAPSHOTS * TOTAL_MOLECULES].tofile(
        DST / f"00/m_loc{MACRO_CODE}_00.dat"
    )

    # m_bound: first N_SNAPSHOTS * TOTAL_MOLECULES int32 values
    m_bound = np.fromfile(SRC / f"00/m_bound{MACRO_CODE}_00.dat", dtype=np.int32)
    m_bound[: N_SNAPSHOTS * TOTAL_MOLECULES].tofile(
        DST / f"00/m_bound{MACRO_CODE}_00.dat"
    )

    # m_bind_t: all rows with timestamp < cutoff_time (CSV text)
    _filter_csv_by_time(
        SRC / f"00/m_bind_t{MACRO_CODE}_00.dat",
        DST / f"00/m_bind_t{MACRO_CODE}_00.dat",
        cutoff_time,
    )

    # f_deg_list: all rows with timestamp < cutoff_time (CSV text)
    _filter_csv_by_time(
        SRC / f"00/f_deg_list{MACRO_CODE}_00.dat",
        DST / f"00/f_deg_list{MACRO_CODE}_00.dat",
        cutoff_time,
    )

    # mfpt: keep as-is (one float64 per molecule, 21105 values)
    shutil.copy2(
        SRC / f"00/mfpt{MACRO_CODE}_00.dat",
        DST / f"00/mfpt{MACRO_CODE}_00.dat",
    )

    # ── Report ─────────────────────────────────────────────────────────
    total_bytes = sum(f.stat().st_size for f in DST.rglob("*") if f.is_file())
    print(f"Fixture created at {DST}")
    print(f"Total size: {total_bytes / 1024:.0f} KB ({total_bytes / 1024 / 1024:.1f} MB)")
    print(f"tsave cutoff: {cutoff_time}")


def _filter_csv_by_time(src_path, dst_path, cutoff_time):
    """Copy CSV rows where the first column (timestamp) is < cutoff_time."""
    kept = 0
    with open(src_path) as fin, open(dst_path, "w") as fout:
        for line in fin:
            timestamp = float(line.split(",", 1)[0])
            if timestamp < cutoff_time:
                fout.write(line)
                kept += 1
            else:
                break  # CSV is sorted by time
    print(f"  {src_path.name}: {kept} rows (t < {cutoff_time:.4f})")


if __name__ == "__main__":
    main()
