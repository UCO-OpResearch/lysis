"""Provenance metadata for lysis simulation outputs.

Two parallel concerns share this package:

- :mod:`~lysis.tools.provenance.execution` — Python-side stamps (lysis
  version, dirty bit, timestamp, hostname) gathered at HDF5 write time by
  :class:`~lysis.dataio.datastore.DataStore`.
- :mod:`~lysis.tools.provenance.binary` — Fortran-binary staleness check
  run by :class:`~lysis.execution.fortran.FortranRunner` before launching
  a compiled simulation binary, comparing the binary's embedded build
  stamp against the current ``src/fortran/`` source tree.

The public API is re-exported here so callers can ``from
lysis.tools.provenance import gather_execution_provenance`` (or
``verify_binary_matches_source``) without reaching into a submodule.
"""

from .binary import (
    StaleBinaryError,
    allow_stale_from_env,
    gather_fortran_source_provenance,
    query_binary_version,
    verify_binary_matches_source,
)
from .execution import gather_execution_provenance

__all__ = [
    "StaleBinaryError",
    "allow_stale_from_env",
    "gather_execution_provenance",
    "gather_fortran_source_provenance",
    "query_binary_version",
    "verify_binary_matches_source",
]
