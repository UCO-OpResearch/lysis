"""One-time migration of the v2.0.0 HDF5 provenance attribute schema.

Issue #61 renamed and retyped the provenance attributes stamped onto the
per-scale params groups (``micro_data`` / ``macro_data``) without bumping
``dataspec_version`` (the format is codified for v1.0.0).  Files written by
earlier v2.0.0 code therefore carry the *old* attribute names/types:

============================  =================================  ============
Old attr                      New attr                           Note
============================  =================================  ============
``execution_version``         ``pipeline_version``               rename
``execution_dirty`` (bool)    ``pipeline_dirty`` (3-state str)   bool→string
``execution_timestamp``       ``pipeline_timestamp``             rename
``execution_hostname``        ``pipeline_hostname``              rename
``execution_backend``         ``backend_type``                   rename + move
``binary_commit``             ``backend_commit``                 rename
``binary_dirty``              ``backend_dirty``                  rename
``binary_compiler``           ``backend_compiler``               rename
``binary_source``             ``backend_historical`` (= True)    drop prefix
``stale_binary_override``     ``stale_backend_override``         rename
============================  =================================  ============

``init_*`` attributes are deliberately left untouched.

This module exposes a single **pure** function, :func:`migrate_provenance_group`,
which operates on a plain ``dict`` snapshot of a group's ``attrs`` (no h5py, no
I/O) so it is trivially unit-testable and shared by the (now archived) one-off
driver ``archive/scripts/migrate_provenance_attrs.py``.  The pure transform is
kept in the package (not archived) so the logic stays unit-tested and reusable
as a reference for any future attribute migration.
"""

# Old attribute names (literals — the CONST constants no longer define them).
_OLD_EXECUTION_VERSION = "execution_version"
_OLD_EXECUTION_DIRTY = "execution_dirty"
_OLD_EXECUTION_TIMESTAMP = "execution_timestamp"
_OLD_EXECUTION_HOSTNAME = "execution_hostname"
_OLD_EXECUTION_BACKEND = "execution_backend"
_OLD_BINARY_COMMIT = "binary_commit"
_OLD_BINARY_DIRTY = "binary_dirty"
_OLD_BINARY_COMPILER = "binary_compiler"
_OLD_BINARY_SOURCE = "binary_source"
_OLD_STALE_BINARY_OVERRIDE = "stale_binary_override"

# New attribute names (kept in sync with lysis.config.constants.CONST).
_NEW_PIPELINE_VERSION = "pipeline_version"
_NEW_PIPELINE_DIRTY = "pipeline_dirty"
_NEW_PIPELINE_TIMESTAMP = "pipeline_timestamp"
_NEW_PIPELINE_HOSTNAME = "pipeline_hostname"
_NEW_BACKEND_TYPE = "backend_type"
_NEW_BACKEND_COMMIT = "backend_commit"
_NEW_BACKEND_DIRTY = "backend_dirty"
_NEW_BACKEND_COMPILER = "backend_compiler"
_NEW_BACKEND_HISTORICAL = "backend_historical"
_NEW_STALE_BACKEND_OVERRIDE = "stale_backend_override"

#: Direct one-to-one renames where the value is carried over verbatim.
_RENAMES = {
    _OLD_EXECUTION_VERSION: _NEW_PIPELINE_VERSION,
    _OLD_EXECUTION_TIMESTAMP: _NEW_PIPELINE_TIMESTAMP,
    _OLD_EXECUTION_HOSTNAME: _NEW_PIPELINE_HOSTNAME,
    _OLD_BINARY_COMMIT: _NEW_BACKEND_COMMIT,
    _OLD_BINARY_DIRTY: _NEW_BACKEND_DIRTY,
    _OLD_BINARY_COMPILER: _NEW_BACKEND_COMPILER,
    _OLD_STALE_BINARY_OVERRIDE: _NEW_STALE_BACKEND_OVERRIDE,
}

#: All old keys the migration consumes (and therefore deletes when present).
OLD_KEYS = frozenset(
    list(_RENAMES)
    + [
        _OLD_EXECUTION_DIRTY,
        _OLD_EXECUTION_BACKEND,
        _OLD_BINARY_SOURCE,
    ]
)

#: Keys whose presence means a group carries a Fortran backend stamp.
_BACKEND_STAMP_KEYS = (_OLD_BINARY_COMMIT, _OLD_BINARY_DIRTY, _OLD_BINARY_COMPILER)


def _coerce_str(value) -> str:
    """Return *value* as a ``str``, decoding ``bytes`` (HDF5 may hand back either)."""
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    return str(value)


def _dirty_bool_to_str(value) -> str:
    """Map an old bool ``execution_dirty`` to the 3-state string.

    Existing v2.0.0 files only ever stored ``False`` (the bool collapsed the
    ``"unknown"`` case into clean), so ``False`` → ``"clean"`` and ``True`` →
    ``"dirty"``.  If a file somehow already holds the 3-state string, it is
    passed through unchanged.
    """
    if isinstance(value, str):
        return value
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    return "dirty" if bool(value) else "clean"


def migrate_provenance_group(attrs: dict) -> "tuple[dict, set]":
    """Compute the provenance-attr migration for one scale group.

    :param attrs: A plain snapshot of an HDF5 group's ``attrs`` (e.g.
        ``dict(group.attrs)``).  Not mutated.
    :type attrs: dict
    :return: ``(sets, deletes)`` where *sets* maps new attribute names to the
        values to write, and *deletes* is the set of old attribute names to
        remove.  Both are empty when the group carries no old provenance attrs
        (so the call is idempotent — re-running on an already-migrated group is
        a no-op).
    :rtype: tuple[dict, set]
    """
    present_old = OLD_KEYS & attrs.keys()
    if not present_old:
        return {}, set()

    sets: dict = {}
    deletes: set = set(present_old)

    # 1. Verbatim renames.
    for old, new in _RENAMES.items():
        if old in attrs:
            sets[new] = attrs[old]

    # 2. execution_dirty (bool) -> pipeline_dirty (3-state string).
    if _OLD_EXECUTION_DIRTY in attrs:
        sets[_NEW_PIPELINE_DIRTY] = _dirty_bool_to_str(attrs[_OLD_EXECUTION_DIRTY])

    # 3. backend_type: prefer the per-scale execution_backend value; otherwise
    #    default to "fortran" when the group carries a Fortran backend stamp.
    if _OLD_EXECUTION_BACKEND in attrs:
        sets[_NEW_BACKEND_TYPE] = _coerce_str(attrs[_OLD_EXECUTION_BACKEND])
    elif any(k in attrs for k in _BACKEND_STAMP_KEYS):
        sets[_NEW_BACKEND_TYPE] = "fortran"

    # 4. binary_source ("historical:<sha>") -> backend_historical = True.
    #    The SHA already lives in backend_commit, so the prefix is dropped.
    if _OLD_BINARY_SOURCE in attrs:
        sets[_NEW_BACKEND_HISTORICAL] = True

    return sets, deletes
