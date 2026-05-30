#!/usr/bin/env python
"""One-time, guard-railed migration of v2.0.0 HDF5 provenance attributes (#61).

Renames/retypes the provenance attributes on the ``micro_data`` / ``macro_data``
groups of existing v2.0.0 HDF5 files (see
:mod:`lysis.tools.provenance.migrate` for the exact mapping).  ``init_*`` attrs
and *all* datasets are left untouched.

Guardrails
----------
* **Dry-run by default.**  Nothing is written unless ``--apply`` is passed.
* **Idempotent.**  Files with no old attrs are classified ``already-migrated``
  and skipped; re-running after a successful pass is a no-op.
* **Scope-limited.**  Only ``dataspec_version == "v2.0.0"`` files that carry old
  provenance attrs *and* are writable by the current user are migrated.  Files
  owned by other users (not writable) are reported, never touched.
* **Per-file integrity verification.**  Before editing, a manifest is built of
  every root/group attribute *and a SHA-256 checksum of every dataset's bytes*.
  After the attr edits the file is closed, reopened read-only, and every dataset
  is re-checksummed.  If ANY dataset checksum, any unrelated attribute, or
  ``dataspec_version`` changed, the script **halts immediately**, prints the
  offending file, and exits non-zero — so you know exactly which file to restore
  from backup.  (Datasets are never accessed for writing; this is belt-and-braces.)
* **Backups.**  ``--backup-dir`` copies each file before editing.  Required under
  ``--apply`` unless ``--no-backup`` is given explicitly.

Usage
-----
    # Dry-run (default) over the standard roots, write a report:
    python scripts/migrate_provenance_attrs.py --report /tmp/prov_migration.txt

    # Apply, backing up every edited file first:
    python scripts/migrate_provenance_attrs.py --apply --backup-dir ~/prov_backup

    # Restrict to one root:
    python scripts/migrate_provenance_attrs.py ~/git/UCO-OpResearch/lysis/data --apply ...
"""

import argparse
import hashlib
import os
import pwd
import shutil
import sys
from datetime import datetime, timezone

import h5py
import numpy as np

# Import the shared pure transform from the installed package.
from lysis.tools.provenance.migrate import migrate_provenance_group, OLD_KEYS

DEFAULT_ROOTS = [
    os.path.expanduser("~/git/UCO-OpResearch/lysis/data"),
    "/shared/lysis-group",
]
SCALE_GROUPS = ("micro_data", "macro_data")
EXPECTED_VERSION = "v2.0.0"
_HASH_ROWS = 4096  # rows of the leading axis read per chunk when checksumming


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _decode(value):
    """Normalise an HDF5 attr value for stable comparison (decode bytes, unwrap np)."""
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return tuple(_decode(v) for v in value.tolist())
    return value


def _norm_attrs(attrs) -> dict:
    """Return a plain, comparison-stable dict from an h5py AttributeManager/dict."""
    return {k: _decode(v) for k, v in dict(attrs).items()}


def _hash_object_array(h, arr):
    """Stably fold an object/vlen array (e.g. h5py vlen strings) into *h*.

    ``ndarray.tobytes()`` on an ``object`` array hashes Python *pointers*
    (non-deterministic across reads), so iterate and feed each element's
    own bytes with a length-prefixed separator instead.
    """
    for item in arr.ravel():
        data = item if isinstance(item, bytes) else str(item).encode("utf-8")
        h.update(len(data).to_bytes(8, "little"))
        h.update(data)


def _dataset_checksum(ds) -> str:
    """SHA-256 of a dataset's contents, read in bounded chunks.

    Numeric/array datasets hash their raw bytes; object (vlen-string)
    datasets hash element-by-element so the digest is stable across reads.
    """
    h = hashlib.sha256()
    h.update(repr((ds.shape, str(ds.dtype))).encode())  # bind shape+dtype into the digest
    is_object = ds.dtype == object
    if ds.shape == ():  # scalar dataset
        if is_object:
            _hash_object_array(h, np.asarray(ds[()], dtype=object))
        else:
            h.update(np.asarray(ds[()]).tobytes())
        return h.hexdigest()
    if 0 in ds.shape:  # empty dataset — nothing to read
        return h.hexdigest()
    n = ds.shape[0]
    for i in range(0, n, _HASH_ROWS):
        chunk = ds[i : i + _HASH_ROWS]
        if is_object:
            _hash_object_array(h, np.asarray(chunk, dtype=object))
        else:
            h.update(np.ascontiguousarray(chunk).tobytes())
    return h.hexdigest()


def _build_manifest(path) -> dict:
    """Read-only manifest: root attrs, per-scale-group attrs, all dataset checksums."""
    manifest = {"root_attrs": {}, "group_attrs": {}, "datasets": {}}
    with h5py.File(path, "r") as f:
        manifest["root_attrs"] = _norm_attrs(f.attrs)
        for g in SCALE_GROUPS:
            if g in f:
                manifest["group_attrs"][g] = _norm_attrs(f[g].attrs)

        def _visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                manifest["datasets"][name] = (
                    obj.shape,
                    str(obj.dtype),
                    _dataset_checksum(obj),
                )

        f.visititems(_visit)
    return manifest


def _file_owner(path) -> str:
    try:
        return pwd.getpwuid(os.stat(path).st_uid).pw_name
    except (KeyError, OSError):
        return "?"


def _iter_h5(roots):
    for root in roots:
        for dirpath, _dirs, files in os.walk(root):
            for name in files:
                if name.endswith(".h5"):
                    yield os.path.join(dirpath, name)


# --------------------------------------------------------------------------- #
# Classification
# --------------------------------------------------------------------------- #
def classify(path) -> str:
    """Return one of the migration categories for *path* (read-only)."""
    try:
        with h5py.File(path, "r") as f:
            ver = _decode(f.attrs.get("dataspec_version"))
            if ver != EXPECTED_VERSION:
                return "not-v2.0.0"
            has_old = any(
                g in f and (OLD_KEYS & f[g].attrs.keys()) for g in SCALE_GROUPS
            )
    except Exception:  # noqa: BLE001 — any read failure is a category, not a crash
        return "read-error"
    if not has_old:
        return "already-migrated"
    if not os.access(path, os.W_OK):
        return "not-writable"
    return "will-migrate"


# --------------------------------------------------------------------------- #
# Migration + verification
# --------------------------------------------------------------------------- #
class IntegrityError(RuntimeError):
    """Raised when post-migration verification detects an unexpected change."""


def _expected_group_attrs(pre_attrs, sets, deletes) -> dict:
    expected = {k: v for k, v in pre_attrs.items() if k not in deletes}
    expected.update({k: _decode(v) for k, v in sets.items()})
    return expected


def migrate_file(path, *, backup_dir=None):
    """Migrate one file in place with full pre/post integrity verification.

    :raises IntegrityError: if any dataset checksum, unrelated attribute, or the
        dataspec_version differs after the edit.  The file's edits are NOT rolled
        back — recover from the backup and investigate.
    :return: A dict describing the per-group attribute changes applied.
    """
    pre = _build_manifest(path)

    if backup_dir is not None:
        os.makedirs(backup_dir, exist_ok=True)
        # Flatten the path into the backup name to avoid collisions across roots.
        flat = path.lstrip(os.sep).replace(os.sep, "__")
        shutil.copy2(path, os.path.join(backup_dir, flat))

    # Compute the planned changes per group and the expected post-state.
    plan = {}
    expected_group_attrs = {}
    for g, pre_attrs in pre["group_attrs"].items():
        sets, deletes = migrate_provenance_group(dict(pre_attrs))
        if sets or deletes:
            plan[g] = (sets, deletes)
            expected_group_attrs[g] = _expected_group_attrs(pre_attrs, sets, deletes)

    if not plan:  # nothing to do (shouldn't happen for a will-migrate file)
        return {}

    # Apply: attrs only — datasets are never opened for writing.
    with h5py.File(path, "a") as f:
        for g, (sets, deletes) in plan.items():
            group = f[g]
            for k, v in sets.items():
                group.attrs[k] = v
            for k in deletes:
                if k in group.attrs:
                    del group.attrs[k]

    # Verify: reopen read-only, rebuild manifest, compare exhaustively.
    post = _build_manifest(path)

    problems = []
    if post["root_attrs"] != pre["root_attrs"]:
        problems.append(f"root attrs changed: {pre['root_attrs']} -> {post['root_attrs']}")
    if set(post["datasets"]) != set(pre["datasets"]):
        problems.append("dataset set changed (added/removed datasets)")
    for name, sig in pre["datasets"].items():
        if post["datasets"].get(name) != sig:
            problems.append(f"dataset changed: {name}")
    for g, expected in expected_group_attrs.items():
        got = post["group_attrs"].get(g, {})
        if got != expected:
            problems.append(
                f"group '{g}' attrs mismatch:\n    expected={expected}\n    got     ={got}"
            )
    # Groups we did not plan to touch must be byte-for-byte identical.
    for g, pre_attrs in pre["group_attrs"].items():
        if g not in expected_group_attrs and post["group_attrs"].get(g) != pre_attrs:
            problems.append(f"untouched group '{g}' attrs changed unexpectedly")

    if problems:
        raise IntegrityError(
            f"Integrity check FAILED for:\n  {path}\n" + "\n".join(f"  - {p}" for p in problems)
        )

    return {g: plan[g] for g in plan}


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _format_plan(plan) -> str:
    lines = []
    for g, (sets, deletes) in plan.items():
        added = ", ".join(f"+{k}={v!r}" for k, v in sorted(sets.items()))
        removed = ", ".join(f"-{k}" for k in sorted(deletes))
        lines.append(f"      {g}: {added} | {removed}")
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "roots", nargs="*", default=None,
        help=f"Directories to scan (default: {DEFAULT_ROOTS}).",
    )
    parser.add_argument("--apply", action="store_true", help="Write changes (default: dry-run).")
    parser.add_argument("--backup-dir", default=None, help="Copy each file here before editing.")
    parser.add_argument(
        "--no-backup", action="store_true",
        help="Proceed under --apply without a backup dir (NOT recommended).",
    )
    parser.add_argument("--report", default=None, help="Write the full report to this path too.")
    args = parser.parse_args(argv)

    roots = args.roots if args.roots else DEFAULT_ROOTS
    if args.apply and args.backup_dir is None and not args.no_backup:
        parser.error(
            "--apply requires --backup-dir (or pass --no-backup to override). "
            "Backups are the recovery path if integrity verification fails."
        )

    buckets = {
        "will-migrate": [], "already-migrated": [], "not-v2.0.0": [],
        "not-writable": [], "read-error": [],
    }
    for path in _iter_h5(roots):
        buckets[classify(path)].append(path)

    out = []

    def emit(line=""):
        out.append(line)
        print(line)

    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    emit(f"# Provenance attr migration — {'APPLY' if args.apply else 'DRY-RUN'} — {stamp}")
    emit(f"# roots: {roots}")
    emit("")
    for cat, files in buckets.items():
        emit(f"{cat}: {len(files)}")
    emit("")

    if buckets["not-writable"]:
        emit("## not-writable (owned by another user — reported, NOT touched):")
        for p in sorted(buckets["not-writable"]):
            emit(f"  [{_file_owner(p)}] {p}")
        emit("")

    migrated, failed = 0, None
    emit("## will-migrate:")
    for path in sorted(buckets["will-migrate"]):
        if not args.apply:
            # Show the planned diff without writing.
            try:
                manifest = _build_manifest(path)
            except Exception as e:  # noqa: BLE001
                emit(f"  READ-ERROR {path}: {e}")
                continue
            plan = {}
            for g, pre_attrs in manifest["group_attrs"].items():
                sets, deletes = migrate_provenance_group(dict(pre_attrs))
                if sets or deletes:
                    plan[g] = (sets, deletes)
            emit(f"  WOULD MIGRATE {path}")
            emit(_format_plan(plan))
        else:
            try:
                plan = migrate_file(path, backup_dir=args.backup_dir)
            except IntegrityError as e:
                emit("")
                emit("!!! HALTING — integrity verification failed !!!")
                emit(str(e))
                failed = path
                break
            migrated += 1
            emit(f"  MIGRATED {path}")
            emit(_format_plan(plan))

    emit("")
    if failed:
        emit(f"ABORTED at {failed} after {migrated} successful file(s). Restore it from backup.")
    elif args.apply:
        emit(f"DONE. Migrated {migrated} file(s).")
    else:
        emit(f"DRY-RUN complete. {len(buckets['will-migrate'])} file(s) would be migrated.")

    if args.report:
        with open(args.report, "w") as fh:
            fh.write("\n".join(out) + "\n")
        print(f"\n(report written to {args.report})")

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
