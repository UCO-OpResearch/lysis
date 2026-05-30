"""Provenance metadata for lysis simulation outputs.

Two parallel concerns share this package:

- :mod:`~lysis.tools.provenance.execution` — Python-pipeline stamps (lysis
  version, dirty bit, timestamp, hostname) gathered at init time
  (:func:`gather_init_provenance`) and at HDF5 write time
  (:func:`gather_pipeline_provenance`).  Both functions scope their
  git checks to the ``src/lysis/`` subtree.
- :mod:`~lysis.tools.provenance.binary` — Fortran-binary staleness check
  run by :class:`~lysis.execution.fortran.FortranRunner` before launching
  a compiled simulation binary (:func:`verify_binary_matches_source`),
  plus :func:`gather_backend_provenance` which records the backend's
  embedded commit/dirty/compiler stamps unconditionally on every run.

The public API is re-exported here so callers can ``from
lysis.tools.provenance import gather_pipeline_provenance`` (or any
other public function) without reaching into a submodule.
"""

from .binary import (
    StaleBinaryError,
    allow_stale_from_env,
    gather_backend_provenance,
    gather_fortran_source_provenance,
    gather_historical_backend_provenance,
    query_binary_version,
    verify_binary_matches_source,
)
from .execution import (
    allow_commit_mismatch_from_env,
    allow_dirty_from_env,
    gather_pipeline_provenance,
    gather_init_provenance,
    mark_dirty_warning_emitted,
)

__all__ = [
    "StaleBinaryError",
    "allow_commit_mismatch_from_env",
    "allow_dirty_from_env",
    "allow_stale_from_env",
    "gather_backend_provenance",
    "gather_pipeline_provenance",
    "gather_fortran_source_provenance",
    "gather_historical_backend_provenance",
    "gather_init_provenance",
    "mark_dirty_warning_emitted",
    "query_binary_version",
    "verify_binary_matches_source",
]
