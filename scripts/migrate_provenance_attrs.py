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
  provenance attrs *and* are writable are migrated.  By default a file must be
  OWNED by the current user; pass ``--include-other-owners`` to also migrate
  files merely made writable via group permissions (use only when the owner has
  agreed).  Everything skipped for ownership/permission reasons is reported.
* **Per-file integrity verification.**  Before editing, a manifest is built of
  every root/group attribute *and a SHA-256 checksum of every dataset's bytes*.
  After the attr edits the file is closed, reopened read-only, and every dataset
  is re-checksummed.  If ANY dataset checksum, any unrelated attribute, or
  ``dataspec_version`` changed, the script **halts immediately**, prints the
  offending file, and exits non-zero — so you know exactly which file to restore
  from backup.  (Datasets are never accessed for writing; this is belt-and-braces.)
* **Backups + SHA log.**  ``--backup-dir`` copies each file before editing and
  writes a ``<flattened>.sha256.txt`` sidecar listing every dataset's pre- and
  post-edit checksum.  A backup dir is required under ``--apply`` unless
  ``--no-backup`` is given explicitly.

Usage
-----
    # Dry-run (default) over the standard roots, write a report:
    python scripts/migrate_provenance_attrs.py --report /tmp/prov_migration.txt

    # Apply, backing up every edited file first, 8 workers:
    python scripts/migrate_provenance_attrs.py --apply --backup-dir ~/prov_backup -j 8

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


def _backup_basename(path) -> str:
    """Flattened backup name for *path* (collision-free across roots)."""
    return path.lstrip(os.sep).replace(os.sep, "__")


def _write_sha_log(backup_dir, path, pre, post) -> str:
    """Write a per-file SHA-256 audit log next to the backup.

    Two sections — every dataset's pre-edit checksum, then every dataset's
    post-edit checksum — so the lists line up for a quick diff without having
    to interleave per line.  Written *before* the integrity check raises, so
    the record survives even on a verification failure.

    :return: The path of the written ``.sha256.txt`` file.
    """
    log_path = os.path.join(backup_dir, _backup_basename(path) + ".sha256.txt")
    names = sorted(set(pre["datasets"]) | set(post["datasets"]))

    def _sha(manifest, name):
        entry = manifest["datasets"].get(name)
        return entry[2] if entry else "MISSING"

    lines = [
        f"# file: {path}",
        f"# datasets: {len(names)}    (column: dataset_name  sha256)",
        "",
        "[pre-edit]",
    ]
    lines += [f"{name}  {_sha(pre, name)}" for name in names]
    lines += ["", "[post-edit]"]
    lines += [f"{name}  {_sha(post, name)}" for name in names]
    with open(log_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return log_path


# --------------------------------------------------------------------------- #
# Classification
# --------------------------------------------------------------------------- #
def classify(path, *, include_other_owners=False) -> str:
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
    try:
        st = os.stat(path)
    except OSError:
        return "read-error"
    if not os.access(path, os.W_OK):
        return "not-writable"
    # By default only ever modify files the current user OWNS — never another
    # user's data, even when the shared group makes them writable.  The owner
    # opts in to group-writable files via --include-other-owners.
    if st.st_uid != os.getuid() and not include_other_owners:
        return "not-mine"
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
        shutil.copy2(path, os.path.join(backup_dir, _backup_basename(path)))

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

    # Persist the pre/post SHA audit log next to the backup BEFORE the integrity
    # check below can raise, so the record survives a verification failure.
    if backup_dir is not None:
        _write_sha_log(backup_dir, path, pre, post)

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
    parser.add_argument(
        "--include-other-owners", action="store_true",
        help="Also migrate files owned by another user but writable to you "
        "(e.g. group-writable shared files). Use only with the owner's consent.",
    )
    parser.add_argument("--report", default=None, help="Write the full report to this path too.")
    parser.add_argument(
        "--jobs", "-j", type=int, default=1,
        help="Parallel worker processes for --apply (default 1). Each worker "
        "migrates whole files independently — safe because no two workers ever "
        "touch the same file. SHA-256 verification is CPU-bound, so processes "
        "(not threads) give real speed-up. Has no effect on a dry-run.",
    )
    args = parser.parse_args(argv)

    roots = args.roots if args.roots else DEFAULT_ROOTS
    if args.apply and args.backup_dir is None and not args.no_backup:
        parser.error(
            "--apply requires --backup-dir (or pass --no-backup to override). "
            "Backups are the recovery path if integrity verification fails."
        )

    buckets = {
        "will-migrate": [], "already-migrated": [], "not-v2.0.0": [],
        "not-mine": [], "not-writable": [], "read-error": [],
    }
    for path in _iter_h5(roots):
        buckets[classify(path, include_other_owners=args.include_other_owners)].append(path)

    out = []

    def emit(line=""):
        out.append(line)
        print(line, flush=True)  # flush so `| tee` follows progress live

    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    emit(f"# Provenance attr migration — {'APPLY' if args.apply else 'DRY-RUN'} — {stamp}")
    emit(f"# roots: {roots}")
    if args.include_other_owners:
        emit("# include-other-owners: ON (migrating group-writable files you do not own)")
    emit("")
    for cat, files in buckets.items():
        emit(f"{cat}: {len(files)}")
    emit("")

    for cat in ("not-mine", "not-writable"):
        if buckets[cat]:
            emit(f"## {cat} (reported, NOT touched):")
            for p in sorted(buckets[cat]):
                emit(f"  [{_file_owner(p)}] {p}")
            emit("")

    migrated, failed = 0, None
    will = sorted(buckets["will-migrate"])
    emit("## will-migrate:")

    if not args.apply:
        # Dry-run: show the planned diff without writing — read group attrs
        # only (no dataset checksums; those are an --apply-time integrity check).
        for path in will:
            try:
                with h5py.File(path, "r") as f:
                    group_attrs = {
                        g: _norm_attrs(f[g].attrs) for g in SCALE_GROUPS if g in f
                    }
            except Exception as e:  # noqa: BLE001
                emit(f"  READ-ERROR {path}: {e}")
                continue
            plan = {}
            for g, pre_attrs in group_attrs.items():
                sets, deletes = migrate_provenance_group(dict(pre_attrs))
                if sets or deletes:
                    plan[g] = (sets, deletes)
            emit(f"  WOULD MIGRATE {path}")
            emit(_format_plan(plan))

    elif max(1, args.jobs) == 1:
        # Serial apply.
        for i, path in enumerate(will, 1):
            try:
                plan = migrate_file(path, backup_dir=args.backup_dir)
            except Exception as e:  # noqa: BLE001 — halt on ANY failure
                emit("")
                emit("!!! HALTING — verification failed !!!")
                emit(str(e))
                failed = path
                break
            migrated += 1
            emit(f"  MIGRATED [{i}/{len(will)}] {path}")
            emit(_format_plan(plan))

    else:
        # Parallel apply: one whole file per worker process.  No two workers
        # touch the same file, so there is no HDF5 race.  Halt on the first
        # failure — stop scheduling new work and name the offending file.
        import concurrent.futures as _cf  # noqa: PLC0415

        jobs = max(1, args.jobs)
        emit(f"  (running {len(will)} file(s) across {jobs} worker processes)")
        ex = _cf.ProcessPoolExecutor(max_workers=jobs)
        futs = {
            ex.submit(migrate_file, p, backup_dir=args.backup_dir): p for p in will
        }
        try:
            for fut in _cf.as_completed(futs):
                path = futs[fut]
                try:
                    plan = fut.result()
                except Exception as e:  # noqa: BLE001 — halt on ANY failure
                    emit("")
                    emit("!!! HALTING — verification failed !!!")
                    emit(str(e))
                    failed = path
                    break
                migrated += 1
                emit(f"  MIGRATED [{migrated}/{len(will)}] {path}")
                emit(_format_plan(plan))
        finally:
            # cancel_futures stops not-yet-started tasks; in-flight files (≤ jobs)
            # finish on their own.  Either way no new file is begun after a halt.
            ex.shutdown(wait=True, cancel_futures=bool(failed))

    emit("")
    if failed:
        emit(f"ABORTED at {failed} after {migrated} successful file(s). Restore it from backup.")
    elif args.apply:
        emit(f"DONE. Migrated {migrated} file(s).")
    else:
        emit(f"DRY-RUN complete. {len(will)} file(s) would be migrated.")

    if args.report:
        with open(args.report, "w") as fh:
            fh.write("\n".join(out) + "\n")
        print(f"\n(report written to {args.report})")

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
