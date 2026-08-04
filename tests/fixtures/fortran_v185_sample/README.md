v1.85.0 Fortran data fixture (truncated)
=========================================

Generated from:
  - microscale : data/2023-02-01-2200/
  - macroscale : data/2023-02-02-2200/
Script: scripts/create_v185_fixture.py

Micro file code     : _PLG2_tPA01_Q2
Macro-in file code  : _PLG2_tPA01_Q2
Macro-out file code : _PLG2_tPA01_along_Q2
total_edges         : 33545
total_molecules     : 43074
Simulations         : 2
Snapshots per sim   : 3 (t=0 plus 2 saves)
Microscale runs     : 2000
tsave               : [np.float64(0.0), np.float64(10.000303991916667), np.float64(20.00026558975), np.float64(0.0), np.float64(10.000303991916667), np.float64(20.00026558975)]

v1.85.0 stores every simulation in ONE top-level file per dataset, so
f_deg_time has shape (6, 33545) and
Nsave is a vector of length 2, not a scalar.

Note that this fixture deliberately keeps the v1.85.0 quirks that the
converters exist to handle:

  - m_bind_t is absent (v1.85.0 never wrote it)
  - f_deg_time uses 0.0 for BOTH empty edges and unscheduled fibrin
  - deg (a v1.85.0-only degradation-state array) is present
  - params.json uses the legacy key spellings and a null micro_params
