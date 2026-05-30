#!/usr/bin/env python
"""Backfill per-file SHA-256 sidecars for already-migrated HDF5 files (issue #61).

Companion one-off to ``migrate_provenance_attrs.py``.  The earliest migration
runs predated the per-file ``.sha256.txt`` sidecar, so this script reconstructs
that audit trail after the fact: for every backup the migration made, it
checksums the backup copy (== pre-edit) and the live file (== post-edit), writes
the sidecar next to the backup, and flags any dataset whose checksum differs
between the two (a post-hoc proof that the migration changed only attributes,
never data — expected: zero mismatches).

ARCHIVED: like its sibling this is a completed one-off, kept as a template.
See ``archive/README.rst``.  The SAFETY VALVE below blocks it from running
until the marked line is removed.
"""
import concurrent.futures as cf
import glob
import importlib.util
import os
import sys

# ===========================================================================
# SAFETY VALVE  (archived one-off — issue #61)
# ---------------------------------------------------------------------------
# This script already ran once.  While the line below is present it refuses to
# do anything (it writes .sha256.txt sidecars, so it is gated like its sibling).
# DELETE the single line below to actually run it when reusing as a template.
_SAFETY_VALVE = True
# ===========================================================================

# Reuse the migration script's checksum/manifest helpers from the sibling file.
_SIBLING = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "migrate_provenance_attrs.py")
_spec = importlib.util.spec_from_file_location("_migrate_sibling", _SIBLING)
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)

# (root, backup_dir) pairs to reconcile.  Adjust to match the actual backup
# locations used by the original migration run.
PAIRS = [
    (os.path.expanduser("~/git/UCO-OpResearch/lysis/data"),
     os.path.expanduser("~/prov_backup_local")),
    ("/shared/lysis-group", os.path.expanduser("~/prov_backup_shared")),
]


def worker(bdir, bf, live):
    backup_path = os.path.join(bdir, bf)
    pre = m._build_manifest(backup_path)        # backup == pre-edit
    if live is None:
        return ("ORPHAN", bf, [])
    post = m._build_manifest(live)              # live == post-edit
    m._write_sha_log(bdir, live, pre, post)     # sidecar name == bf + .sha256.txt
    diffs = [n for n, sig in pre["datasets"].items()
             if post["datasets"].get(n) != sig]
    if set(pre["datasets"]) != set(post["datasets"]):
        diffs.append("<dataset-set-changed>")
    return ("OK" if not diffs else "MISMATCH", bf, diffs)


def main():
    if globals().get("_SAFETY_VALVE", False):
        print(
            "SAFETY VALVE ENGAGED: this is an archived one-off; refusing to run.\n"
            "Delete the `_SAFETY_VALVE = True` line near the top to enable it "
            "(see archive/README.rst).",
            flush=True,
        )
        return 0

    tasks = []
    for root, bdir in PAIRS:
        if not os.path.isdir(bdir):
            continue
        livemap = {m._backup_basename(p): p
                   for p in glob.iglob(os.path.join(root, "**", "*.h5"), recursive=True)}
        for bf in os.listdir(bdir):
            if bf.endswith(".sha256.txt"):
                continue
            if os.path.exists(os.path.join(bdir, bf + ".sha256.txt")):
                continue  # already has a sidecar
            tasks.append((bdir, bf, livemap.get(bf)))

    print(f"backfilling {len(tasks)} sidecar(s) across {PAIRS}", flush=True)
    ok = orphan = 0
    mismatches = []
    with cf.ProcessPoolExecutor(max_workers=8) as ex:
        futs = {ex.submit(worker, *t): t[1] for t in tasks}
        done = 0
        for fut in cf.as_completed(futs):
            status, bf, diffs = fut.result()
            done += 1
            if status == "OK":
                ok += 1
            elif status == "ORPHAN":
                orphan += 1
                print(f"  ORPHAN (no live file): {bf}", flush=True)
            else:
                mismatches.append((bf, diffs))
                print(f"  !!! MISMATCH {bf}: {diffs}", flush=True)
            if done % 100 == 0:
                print(f"  ...{done}/{len(tasks)}", flush=True)

    print(f"\nDONE. sidecars OK={ok}  orphan={orphan}  mismatched={len(mismatches)}", flush=True)
    if mismatches:
        print("MISMATCHES (data differs between backup and live — investigate!):", flush=True)
        for bf, diffs in mismatches:
            print(f"  {bf}: {diffs}", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
