==========================================
Data Handling System Completion Plan
==========================================

Executive Summary
-----------------

The lysis project has a well-architected data handling system with clear separation of concerns across five main modules:

- **dataspec.py**: Complete specification definitions with ``DataSpec`` wrapper class,
  including v1.95.0, v1.99.0, and v2.0.0 specs (DONE)
- **fileops.py**: File I/O operations for all five storage formats, including Fortran
  log-file parsing (DONE)
- **dataconvert.py**: Format conversion logic with automatic multi-step routing via
  spanning tree (DONE)
- **datastore.py**: High-level Python API with read/write access, lazy macroscale_in
  generation (DONE)
- **paramcheck.py**: Strict parameter validation and Fortran log verification (DONE)

**Overall Status: ~97% Complete**

**Recent Updates (Feb 2026):**

- Implemented ``paramcheck`` module for strict parameter loading, Fortran log
  parsing, and parameter verification (Task #31)
- Added v1.95.0 data specification with copy/deepcopy support and ``replace()``
  method for ``DataSpec``
- Implemented v1.95.0 <-> v1.99.0 data converters and parameter conversion support
- Added automatic multi-step conversion routing via spanning tree — no longer
  limited to direct version pairs
- Created real data integration tests using truncated Fortran fixture data
  (68 tests in ``tests/data/test_real_data.py``)
- Integrated ``np_macroscale`` with ``DataStore`` for both macroscale_in reads
  and macroscale_out writes in v2.0.0 format
- Added lazy ``macroscale_in`` generation in ``DataStore``
- Implemented ``_read_file_parsed()`` for reading micro parameters from Fortran
  log files
- Updated Fortran spec to reflect that microparameters are parsed from log files
  (not stored as JSON)
- Added parameter validation against Fortran log files on every Fortran data load
- Enhanced CLI with multiple file_code support, parameter overrides, and
  separate aliases/overrides parameters
- Added unit tests for parameters module (53 tests) and paramcheck module
  (30 tests)
- Added ``np_macroscale`` test suite (40 tests)
- Removed filename references from constants and parameters (now in dataspec)
- Added null write function for datasets with ``None`` storage type
- Expanded test suite to 630 tests across 10 test files (from 291 across 6)


Architecture Overview
---------------------

.. code-block:: text

   +--------------------------------------------------+
   |   User Application Code                          |
   |   (Simulation scripts, analysis notebooks)       |
   +------------------+-------------------------------+
                      |
   +------------------v-------------------------------+
   |   DataStore API (datastore.py)                   |
   |   - Read/write: dot-access, collections, params  |
   |   - Per-simulation and combined views            |
   |   - HDF5 spec version validation on open         |
   |   - Initialization (create, add macroscale)      |
   |   - Lazy macroscale_in generation                |
   |   Status: COMPLETE                               |
   +------------------+-------------------------------+
                      |
   +------------------v-------------------------------+
   |   Parameter Validation (paramcheck.py)           |
   |   - Strict loading (missing-param detection)     |
   |   - Fortran log parsing & verification           |
   |   - Alias & override resolution                  |
   |   Status: COMPLETE                               |
   +------------------+-------------------------------+
                      |
   +------------------v-------------------------------+
   |   Data Conversion Layer (dataconvert.py)         |
   |   - Automatic multi-step routing (spanning tree) |
   |   - v1.95.0 <-> v1.99.0 <-> v2.0.0 conversion   |
   |   - Safe type casting with validation            |
   |   - Generic structured array converters          |
   |   - Event-log reconstruction (m_bound)           |
   |   Status: COMPLETE                               |
   +------------------+-------------------------------+
                      |
   +------------------v-------------------------------+
   |   Data I/O Layer (fileops.py)                    |
   |   - Format-agnostic read/write with routing      |
   |   - Collection-level multi-file handling         |
   |   - Fortran log-file parsing                     |
   |   - Parameter validation on Fortran data load    |
   |   - Format string lookup from constants          |
   |   Status: COMPLETE - All 5 formats supported     |
   +------------------+-------------------------------+
                      |
   +------------------v-------------------------------+
   |   Specification Layer (dataspec.py)              |
   |   - DataSpec wrapper with version tracking       |
   |   - v1.95.0, v1.99.0, v2.0.0 schema definitions |
   |   - Type/shape validation                        |
   |   - Hidden field auto-population                 |
   |   - Copy/deepcopy and replace() support          |
   |   Status: COMPLETE                               |
   +------------------+-------------------------------+
                      |
   +------------------v-------------------------------+
   |   Storage Backends                               |
   |   HDF5 (done)  Text (done)  Binary (done)       |
   |   JSON (done)  Fortran (done)  Parsed (done)     |
   +--------------------------------------------------+


Completed Work
--------------

The following items from the original plan have been completed.

1. HDF5 Attribute Writing Bug -- COMPLETED
++++++++++++++++++++++++++++++++++++++++++

:File: ``fileops.py``
:Task: #3

Fixed. The hard-coded ``[:6]`` slice has been replaced with proper path
construction.

2. Debug Print Statement -- COMPLETED
++++++++++++++++++++++++++++++++++++++

:File: ``dataconvert.py``
:Task: #4

Fixed. Debug print statement removed.

3. DataStore Read-Only Implementation -- COMPLETED
++++++++++++++++++++++++++++++++++++++++++++++++++

:File: ``datastore.py``
:Task: #1

The ``DataStore`` class provides read-only access to HDF5 simulation data:

- ``__init__()`` - Opens HDF5, validates spec version, detects collections,
  loads parameters, constructs ``DataCollection`` objects
- ``__getattr__()`` - Dot-access to collections (e.g. ``ds.microscale_out``)
- ``collections`` property - Lists available data collections
- ``micro_params`` / ``macro_params`` - Loaded parameter objects
- ``close()`` / context manager - Resource management
- ``SimulationView`` - Per-simulation dot-access to datasets
- ``DataCollection`` - Combined and per-simulation collection interface

Covered by 124 unit tests in ``tests/data/test_datastore.py``.

4. Comprehensive Unit Tests -- COMPLETED
+++++++++++++++++++++++++++++++++++++++++

:Task: #2

630 unit tests across 10 test files:

- ``tests/data/test_dataspec.py`` (126 tests) - Spec parsing, validation,
  DataSpec wrapper, hidden fields, v1.95.0 specs, copy/deepcopy, replace()
- ``tests/data/test_datastore.py`` (124 tests) - DataStore read/write API,
  initialization, collections, parameter loading, version validation,
  write-through, mode control, lazy macroscale_in generation
- ``tests/data/test_dataconvert.py`` (114 tests) - Conversion logic, safe
  type casting, round-trip integrity, multi-step routing, v1.95.0 converters
- ``tests/data/test_real_data.py`` (68 tests) - End-to-end integration tests
  with real truncated Fortran fixture data
- ``tests/config/test_parameters.py`` (53 tests) - Parameter parsing,
  validation, overlap checks
- ``tests/data/test_fileops.py`` (49 tests) - Read/write for all formats,
  Fortran log parsing, parameter validation on load
- ``tests/test_np_macroscale.py`` (40 tests) - Macroscale simulation with
  DataStore integration
- ``tests/config/test_paramcheck.py`` (30 tests) - Strict parameter loading,
  log parsing, verification
- ``tests/cli/test_convert.py`` (14 tests) - CLI convert command
- ``tests/cli/test_validate.py`` (12 tests) - CLI validate command

5. JSON Writer -- COMPLETED
++++++++++++++++++++++++++++

:File: ``fileops.py``
:Task: #5

``_write_file_json()`` is fully implemented. All five storage format writers
are now in the ``data_writers`` registry.

6. Data Validation in DataStore -- COMPLETED
++++++++++++++++++++++++++++++++++++++++++++

:File: ``datastore.py``
:Task: #7

Validation performed on DataStore open:

- Spec version validation (``dataspec_version`` root attribute)
- Dependency validation (macroscale_out requires microscale_out)
- Parameter loading and validation from HDF5 attributes

7. All Data Converters -- COMPLETED
++++++++++++++++++++++++++++++++++++

:Files: ``dataconvert.py``
:Tasks: #14, #15, #16, #17, #18

All 44 data converters are implemented in both directions:

- **Microscale**: All datasets, round-trip tested
- **Macroscale input**: All 5 datasets (bin_edge_proportions, bin_edge_tpa_leaving_time,
  binned_fiber_degrade_time, binned_fiber_degraded, edge_grid_neighbors)
- **Macroscale output**: All datasets including ``m_bound`` (event-log
  reconstruction via ``replay_event_log_to_snapshot``)
- **Macroscale input generation** (``generate_macroscale_in``) from microscale output

Generic converter infrastructure: ``convert_structured_grid_fields`` and
``convert_location_snapshot`` with field validation.

8. HDF5 Specification Versioning -- COMPLETED
++++++++++++++++++++++++++++++++++++++++++++++

:File: ``datastore.py``, ``fileops.py``
:Task: #11

- ``dataspec_version`` stored as root attribute in HDF5 files
- Version validated on every read and write
- ``COMPATIBLE_DATASPEC_VERSION`` module constant for tag-mismatch warnings

9. CLI Tools -- COMPLETED
+++++++++++++++++++++++++

:File: ``lysis/cli/``
:Task: #13

Click-based CLI with two commands:

- ``lysis convert`` - Data format conversion between spec versions, with
  support for multiple file codes, parameter overrides, and aliases
- ``lysis validate`` - Data validation against a specification

10. Macroscale Output Test Script -- COMPLETED
++++++++++++++++++++++++++++++++++++++++++++++

:Task: #19

Test script created for macroscale output conversion round-trip testing.

11. Package Reorganization -- COMPLETED
+++++++++++++++++++++++++++++++++++++++

:Tasks: #20, #21, #22, #25, #26, #27

- ``run.py`` moved from ``execution/`` to ``config/`` (ontology alignment)
- Circular import between ``geometry`` and ``execution`` resolved via
  ``TYPE_CHECKING`` guard in ``edge_grid.py``
- Backward-compatibility ``lysis/util/`` shim package removed
- All imports updated to canonical module paths
- Filename references removed from ``constants.py`` and ``parameters.py``
  (now contained in dataspec)

12. DataStore Write Access -- COMPLETED
+++++++++++++++++++++++++++++++++++++++

:File: ``datastore.py``
:Task: #30

The ``DataStore`` class now supports full read/write access:

- ``DataStore.__init__(mode=)`` accepts ``"r"`` (read-only, default),
  ``"a"`` (read/write), or ``"w"`` (create/truncate), passed directly to
  ``h5py.File``
- ``DataStore.create(run_code, path, micro_params)`` class method creates
  a new HDF5 file with microscale parameters and empty datasets, returned
  in ``"a"`` mode
- ``ds.initialize_macroscale(macro_params)`` adds macroscale parameters and
  empty per-simulation datasets in place; raises ``IOError`` if read-only
- Write-through via live ``h5py.Dataset`` objects: resize, slice-write,
  and incremental append in ``"a"`` mode
- ``DataStatus`` enum wired into ``DataStore`` (``status`` property,
  ``INITIALIZED`` state set on open)

Covered by 124 unit tests in ``tests/data/test_datastore.py``.

13. Parameter Validation and Verification -- COMPLETED
++++++++++++++++++++++++++++++++++++++++++++++++++++++

:File: ``paramcheck.py``, ``fileops.py``, ``parameters.py``
:Task: #31

The ``paramcheck`` module provides strict parameter loading and Fortran log
verification on top of ``parse_from_basedict()``:

- ``load_micro_params()`` / ``load_macro_params()`` — raise ``ValueError``
  when any independent parameter is missing from stored data, and verify that
  stored dependent parameters match values recalculated from independent ones
- ``parse_micro_log()`` / ``parse_macro_log()`` — parse Fortran simulator
  output files (``micro_*.txt``, ``macro_*.txt``) and return a dict of
  ``{python_name: value}``; unmapped numeric Fortran names raise ``ValueError``
- ``verify_micro_params()`` / ``verify_macro_params()`` — combine parsing with
  comparison against a parameter object, raising ``ValueError`` on mismatch
- ``_read_file_parsed()`` in ``fileops.py`` reads micro parameters from
  Fortran log files
- Fortran spec updated to reflect that microparameters are parsed from log
  files (not stored as JSON)
- Parameter validation against Fortran log files runs on every Fortran data
  load
- Separate aliases and overrides parameters for clean parameter resolution
- Fixed bug where empty parameters could overwrite others

Covered by 30 tests in ``tests/config/test_paramcheck.py`` and 53 tests in
``tests/config/test_parameters.py``.

14. v1.95.0 Data Specification and Converters -- COMPLETED
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

:Files: ``dataspec.py``, ``dataconvert.py``

Added full v1.95.0 support:

- v1.95.0 data specification definitions in ``dataspec.py``
- Copy/deepcopy support and ``replace()`` method for ``DataSpec``
- Bidirectional v1.95.0 <-> v1.99.0 data converters
- Parameter conversion between spec versions

15. Automatic Multi-Step Conversion Routing -- COMPLETED
++++++++++++++++++++++++++++++++++++++++++++++++++++++++

:File: ``dataconvert.py``

Conversion routing now uses a spanning tree to automatically discover
multi-step paths between any two spec versions. This replaces the previous
strategy of requiring all conversions to route through the v1.99.0 <-> v2.0.0
pair. Any combination of v1.95.0, v1.99.0, and v2.0.0 is supported
automatically.

16. Integration Tests with Real Data -- COMPLETED
+++++++++++++++++++++++++++++++++++++++++++++++++

:Task: #12

End-to-end integration tests using truncated real Fortran simulation data:

- 68 tests in ``tests/data/test_real_data.py``
- Test fixture creation script (``scripts/create_test_fixture.py``)
- Covers full v1.95.0 -> v1.99.0 -> v2.0.0 conversion pipelines
- Validates round-trip data integrity with real data

17. np_macroscale DataStore Integration -- COMPLETED
++++++++++++++++++++++++++++++++++++++++++++++++++++

:Files: ``np_macroscale.py``, ``datastore.py``, ``scripts/exec.py``

The macroscale simulation (``np_macroscale``) is now fully integrated with
the DataStore:

- Reads ``macroscale_in`` data directly from ``DataStore``
- Writes ``macroscale_out`` data to ``DataStore`` in v2.0.0 format
- Lazy ``macroscale_in`` generation in ``DataStore`` (generates from
  microscale output on demand)
- Updated exec script (``scripts/exec.py``) for DataStore workflow
- Event dtypes pulled from dataspec instead of being redefined locally

Covered by 40 tests in ``tests/test_np_macroscale.py``.


Remaining Work
--------------

18. Safe Float Type Conversion (MEDIUM)
+++++++++++++++++++++++++++++++++++++++

:File: ``dataconvert.py``
:Impact: Cannot convert float datasets between specifications
:Task: #6

Currently raises ``NotImplementedError`` for float-to-float dtype conversions
(e.g. float32 <-> float64). Note that float-to-int conversion *is* handled
by ``safe_np_int_conversion``.

19. Extract Hard-Coded Constants (LOW)
++++++++++++++++++++++++++++++++++++++

:File: ``dataconvert.py``
:Impact: Improves maintainability
:Task: #8

Hard-coded ``n_bins = 100`` should be configurable constant.

20. Optimize HDF5 Chunking (LOW)
+++++++++++++++++++++++++++++++++

:File: ``fileops.py``
:Impact: Better I/O performance for large datasets
:Task: #9

Currently auto-calculated by h5py, may not be optimal.

21. Complete Documentation (LOW)
++++++++++++++++++++++++++++++++

:Files: All
:Impact: Easier to use and maintain
:Task: #10

Many docstrings have "_description_" placeholders. Documentation for
``np_macroscale.py`` has been substantially updated. Sphinx API documentation
covers all subpackages. RST documentation in ``docs/source/usage/`` has known
syntax issues that need fixing.

22. Update Notebooks (LOW)
++++++++++++++++++++++++++

:Task: #29

Four Jupyter notebooks in ``notebooks/`` still reference the removed
``lysis.util`` package and need to be updated to use canonical module paths.

23. CI Pipeline (LOW)
+++++++++++++++++++++

:Task: #28

Set up continuous integration pipeline to execute the pytest unit test suite
automatically on commits and pull requests.

24. Integrate codeutil.py with DataStore (MEDIUM)
+++++++++++++++++++++++++++++++++++++++++++++++++

:File: ``execution/codeutil.py``
:Impact: Fortran subprocess wrappers cannot use DataStore directly

The Fortran subprocess execution wrappers in ``codeutil.py`` do not yet
reference ``DataStore``. The ``np_macroscale`` simulation is integrated
(item 17), but the Fortran execution helpers still use the older file-based
data flow. These need to be updated so that Fortran simulations can read
from and write to HDF5 via ``DataStore`` and ``dataconvert``.

25. Convert np_macroscale to Quantity (MEDIUM)
++++++++++++++++++++++++++++++++++++++++++++++

:File: ``np_macroscale.py``
:Impact: Dimensioned quantities not tracked at simulation level

``np_macroscale.py`` does not use ``pint.Quantity`` objects for dimensioned
parameters. Other modules (``parameters.py``, ``paramcheck.py``,
``edge_grid.py``, ``dataspec.py``) already use ``Quantity``. Converting
``np_macroscale`` would ensure unit consistency throughout the simulation.

26. forced_unbind Calculation (MEDIUM)
++++++++++++++++++++++++++++++++++++++

:File: ``np_macroscale.py``, ``dataconvert.py``
:Impact: Calculation logic needs review or implementation

The ``forced_unbind`` parameter is used in ``np_macroscale.py`` and computed
during data conversion in ``dataconvert.py``. The calculation task requires
further specification.

27. Update data_specification.rst to Include v1.95.0 (LOW)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

:File: ``docs/source/usage/data_specification.rst``
:Impact: Documentation does not reflect v1.95.0 spec support

The v1.95.0 data specification has been implemented in code (``dataspec.py``,
item 14), but the RST documentation has not yet been updated to describe the
v1.95.0 format.

28. Add Fortran Data Specifications for f_deg and Combined Macro (MEDIUM)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

:Files: ``dataspec.py``, ``dataconvert.py``
:Impact: Cannot import/convert f_deg or combined macro simulation data

Data specifications and converters are needed for the fiber degradation
(``f_deg``) and combined macroscale simulation Fortran formats.

29. MPI4Py for Parallel HDF5 Writes (LOW)
++++++++++++++++++++++++++++++++++++++++++

:Impact: No parallel I/O for multi-process simulations

Implement MPI-based parallel HDF5 writing via ``mpi4py`` to allow concurrent
writes from multiple simulation processes.


Implementation Priority
-----------------------

Phase 1: Core Functionality -- COMPLETED
+++++++++++++++++++++++++++++++++++++++++

All items delivered:

#. **Task #3**: Fix HDF5 attribute writing bug -- COMPLETED
#. **Task #4**: Remove debug print statement -- COMPLETED
#. **Task #1**: Implement DataStore read-only API -- COMPLETED
#. **Task #2**: Create comprehensive unit tests -- COMPLETED

**Deliverable**: Functional read-only data handling system with testing |checkmark|

Phase 2: Conversion & Validation -- COMPLETED
++++++++++++++++++++++++++++++++++++++++++++++

All items delivered:

5. **Task #14**: Complete convert_fiber_degrade_time -- COMPLETED
#. **Task #16**: Implement macroscale output data converters -- COMPLETED
#. **Task #17**: Implement generic macroscale output converters -- COMPLETED
#. **Task #5**: Implement JSON file writer -- COMPLETED
#. **Task #7**: Add data validation to DataStore -- COMPLETED
#. **Task #15**: Verify macroscale input data converters -- COMPLETED
#. **Task #18**: Implement m_bound conversion -- COMPLETED
#. **Task #19**: Create macroscale output test script -- COMPLETED

**Deliverable**: Complete format conversion and validation |checkmark|

Phase 3: Write Access & Polish -- COMPLETED
++++++++++++++++++++++++++++++++++++++++++++

All items delivered:

13. **Task #30**: Implement DataStore initialization and write access -- COMPLETED
#.  **Task #31**: Handle missing parameters in deserialization -- COMPLETED
#.  v1.95.0 data specification and converters -- COMPLETED
#.  Automatic multi-step conversion routing -- COMPLETED
#.  Integration tests with real Fortran data -- COMPLETED
#.  np_macroscale DataStore integration -- COMPLETED

**Deliverable**: Read-write data handling system with parameter validation |checkmark|

Phase 4: Simulation Integration (Current)
++++++++++++++++++++++++++++++++++++++++++

Active work:

19. Integrate ``codeutil.py`` with DataStore (item 24)
#.  Convert ``np_macroscale`` to ``Quantity`` (item 25)
#.  ``forced_unbind`` calculation (item 26)
#.  **Task #6**: Add safe float type conversion (item 18)
#.  **Task #8**: Extract hard-coded constants (item 19)
#.  Add Fortran data specs for f_deg and combined macro (item 28)

**Deliverable**: Full simulation pipeline using DataStore

Phase 5: Quality & Infrastructure
++++++++++++++++++++++++++++++++++

Polish and optimization:

25. **Task #9**: Optimize HDF5 chunk sizes (item 20)
#.  **Task #10**: Complete documentation and docstrings (item 21)
#.  **Task #29**: Update notebooks to remove lysis.util imports (item 22)
#.  **Task #28**: Set up CI pipeline (item 23)
#.  Update ``data_specification.rst`` for v1.95.0 (item 27)
#.  MPI4Py for parallel HDF5 writes (item 29)

**Deliverable**: Production-ready, well-documented system with CI


Current Capabilities (Working Now)
-----------------------------------

Ready for Use
+++++++++++++

#. Reading Fortran v1.95.0 and v1.99.0 formats (text, binary, JSON, parsed)
#. Reading HDF5 v2.0.0 format
#. Writing to HDF5 v2.0.0 format (with spec version metadata)
#. Writing to text, binary, and JSON v1.99.0 format
#. Reading HDF5 via ``DataStore`` with dot-access, per-simulation views,
   and parameter loading
#. HDF5 specification versioning (``dataspec_version`` attribute, validated
   on open/read/write)
#. ``DataSpec`` wrapper class with version tracking, hidden field propagation,
   copy/deepcopy support, and ``replace()`` method
#. Microscale data conversion v1.95.0 <-> v1.99.0 <-> v2.0.0 (all datasets,
   round-trip tested with real data)
#. Macroscale input generation from microscale output (round-trip tested)
#. Macroscale input data conversion v1.95.0 <-> v1.99.0 <-> v2.0.0 (all datasets)
#. Macroscale output data conversion v1.99.0 <-> v2.0.0 (all datasets including m_bound)
#. Automatic multi-step conversion routing via spanning tree
#. Safe type conversion for integers, booleans, strings/objects, and float-to-int
#. Strict parameter loading with missing-parameter detection
#. Fortran log-file parsing and parameter verification
#. Parameter loading and merging with alias/override resolution
#. Dynamic shape resolution from parameters
#. Multi-simulation file handling
#. Format string lookup for text file output (``NUMPY_SAVETXT_FORMATS``)
#. CLI commands: ``lysis convert``, ``lysis validate`` (with multiple file
   codes, parameter overrides, aliases)
#. DataStore initialization: ``create()`` for new files,
   ``initialize_macroscale()`` for adding macroscale in place
#. DataStore read/write mode (``"r"``, ``"a"``, ``"w"``) with read-only
   protection
#. Write-through to HDF5 datasets (resize, slice-write, incremental append)
#. Lazy ``macroscale_in`` generation in DataStore
#. ``np_macroscale`` reads from and writes to DataStore
#. 630 automated unit tests (124 for DataStore, 114 for dataconvert,
   68 for real-data integration)

Needs Work
++++++++++

#. Float-to-float type conversion (not implemented)
#. ``codeutil.py`` Fortran subprocess wrappers not integrated with DataStore
#. ``np_macroscale`` does not use ``Quantity`` for dimensioned parameters
#. ``forced_unbind`` calculation needs review
#. Fortran data specs for f_deg and combined macro simulations
#. Notebook imports (still reference removed ``lysis.util``)

Not Implemented
+++++++++++++++

#. HDF5 chunk optimization
#. Streaming I/O for large datasets
#. Parallel I/O (MPI support)
#. CI pipeline


Design Strengths
----------------

#. **Separation of Concerns**: Clear module boundaries across spec, I/O,
   conversion, validation, and store layers
#. **Format Agnostic**: Easy to add new storage formats; all five backends
   plus parsed format complete
#. **Version Management**: Tag system for forward compatibility; spec version
   stored in HDF5 files and validated on access; three spec versions supported
   (v1.95.0, v1.99.0, v2.0.0)
#. **Automatic Conversion Routing**: Spanning tree discovers multi-step paths
   between any pair of spec versions automatically
#. **Immutable Specs**: ``DataSpec`` wrapper prevents runtime changes to schemas
   while providing dict-like access, version tracking, and copy/replace support
#. **Parameter Validation**: Strict loading detects missing parameters;
   Fortran log parsing and verification ensure consistency
#. **Parameter Integration**: Dynamic shapes computed at I/O time
#. **Dispatcher Pattern**: Clean routing to format-specific functions
#. **Generic Converters**: Reusable ``convert_structured_grid_fields`` and
   ``convert_location_snapshot`` with field validation
#. **Comprehensive Testing**: 630 unit tests covering all modules with
   automated testing, including real-data integration tests
#. **Clean Package Structure**: Canonical imports with no shim layers;
   ``config/`` houses Run per project ontology

Design Weaknesses
-----------------

#. **Tight Path Coupling**: Format strings in spec definitions
#. **Hard-Coded Values**: Magic numbers should be constants (e.g. ``n_bins = 100``)
#. **No CI**: Tests must be executed manually


Risk Assessment
---------------

.. list-table::
   :header-rows: 1
   :widths: 30 12 12 46

   * - Risk
     - Likelihood
     - Impact
     - Mitigation
   * - DataStore write bugs due to untested write path
     - LOW
     - HIGH
     - MITIGATED: 124 tests cover initialization, write-through, and mode control (Task #30)
   * - Parameter mismatches between code and stored data
     - LOW
     - HIGH
     - MITIGATED: paramcheck module validates on every load (Task #31)
   * - Float conversion blocks workflows
     - LOW
     - MEDIUM
     - Implement safe conversion (Phase 4) - Task #6
   * - Poor chunking hurts performance at scale
     - LOW
     - LOW
     - Profile and optimize (Phase 5) - Task #9
   * - No CI means regressions go undetected
     - MEDIUM
     - MEDIUM
     - Set up CI pipeline (Phase 5) - Task #28
   * - Stale notebook code confuses users
     - LOW
     - LOW
     - Update notebook imports (Phase 5) - Task #29
   * - Missing f_deg / combined macro specs block workflows
     - MEDIUM
     - MEDIUM
     - Add Fortran data specs (Phase 4) - item 28


Success Criteria
----------------

Phase 1 Success -- ACHIEVED
++++++++++++++++++++++++++++

- [x] HDF5 attribute bug fixed and verified
- [x] No debug print statements in code
- [x] DataStore read-only methods implemented and tested
- [x] Unit test coverage across all modules
- [x] Can read full simulation datasets via DataStore
- [x] Round-trip conversion preserves data integrity (all data collections)

Phase 2 Success -- ACHIEVED
++++++++++++++++++++++++++++

- [x] Integer and boolean types convert correctly
- [x] String and object types convert correctly
- [x] JSON reading and writing both work
- [x] Data validation catches spec violations (version, dependency)
- [x] Fiber degrade conversion handles all cases
- [x] All macroscale input data converters implemented
- [x] Macroscale output data converters implemented (including m_bound)
- [x] m_bound conversion via event-log reconstruction
- [x] Full bidirectional conversion working for all datasets
- [x] Macroscale output conversion round-trip tested

Phase 3 Success -- ACHIEVED
++++++++++++++++++++++++++++

- [x] DataStore can create new HDF5 files with proper structure
- [x] DataStore can write simulation datasets
- [x] Missing parameters detected and handled during deserialization
- [x] Fortran log files parsed and verified against parameters
- [x] v1.95.0 data specification and converters implemented
- [x] Multi-step conversion routing works automatically
- [x] Integration tests cover real data conversion pipelines
- [x] np_macroscale reads from and writes to DataStore

Phase 4 Success
+++++++++++++++

- [ ] ``codeutil.py`` Fortran wrappers use DataStore
- [ ] ``np_macroscale`` uses ``Quantity`` for dimensioned parameters
- [ ] ``forced_unbind`` calculation reviewed/implemented
- [ ] Float-to-float type conversion works
- [ ] Constants extracted from code
- [ ] Fortran data specs for f_deg and combined macro simulations

Phase 5 Success
+++++++++++++++

- [ ] All docstrings complete with examples
- [ ] User guide documentation reviewed (including v1.95.0)
- [ ] HDF5 I/O performance benchmarked
- [ ] CI pipeline executes tests on every commit
- [ ] Notebooks updated to canonical imports
- [ ] Code review passes with no major issues


Recommended Next Steps
----------------------

#. **Immediate**:

   - Integrate ``codeutil.py`` Fortran subprocess wrappers with DataStore
     (item 24)
   - Add safe float type conversion (Task #6)

#. **Short-term**:

   - Convert ``np_macroscale`` to ``Quantity`` (item 25)
   - Review/implement ``forced_unbind`` calculation (item 26)
   - Extract hard-coded constants (Task #8)

#. **Medium-term**:

   - Add Fortran data specifications for f_deg and combined macro (item 28)
   - Update ``data_specification.rst`` for v1.95.0 (item 27)
   - Set up CI pipeline (Task #28)

#. **Longer-term**:

   - Complete documentation and docstrings (Task #10)
   - Update notebooks to remove ``lysis.util`` imports (Task #29)
   - Optimize HDF5 chunk sizes (Task #9)
   - MPI4Py for parallel HDF5 writes (item 29)


Conclusion
----------

The data handling system has reached a mature state with all core functionality
and the parameter validation system complete:

#. All data converters are complete and tested across three spec versions
   (v1.95.0, v1.99.0, v2.0.0), with automatic multi-step conversion routing
   via spanning tree
#. The ``DataStore`` provides a full read/write API with dot-access to
   collections, per-simulation views, parameter loading, spec version
   validation, initialization (``create()``, ``initialize_macroscale()``),
   write-through to HDF5 datasets, and lazy macroscale_in generation
#. The ``paramcheck`` module provides strict parameter loading with
   missing-parameter detection, Fortran log parsing, and parameter
   verification — resolving the silent-default-fallback problem
#. ``DataSpec`` wrapper class provides version-aware specification management
   with hidden field auto-population, copy/deepcopy, and replace() support
#. HDF5 specification versioning ensures version consistency on every
   read and write operation
#. All five storage backends (HDF5, text, binary, JSON, Fortran) plus
   parsed log-file format support both reading and writing
#. 630 unit tests cover all modules, including 68 real-data integration
   tests with truncated Fortran fixture data
#. CLI tools (``lysis convert``, ``lysis validate``) provide command-line
   access to conversion and validation with multi-file and override support
#. Clean package structure with ``run.py`` in ``config/`` per project ontology
   and no backward-compatibility shim layers
#. ``np_macroscale`` is fully integrated with DataStore for both reads and
   writes

With Phase 3 now complete, the focus shifts to **simulation integration**
(Phase 4): connecting the Fortran execution wrappers to DataStore, adding
``Quantity`` support to ``np_macroscale``, and adding data specs for remaining
Fortran simulation types. The remaining Phase 5 items are polish and
infrastructure work.

:Remaining Effort: ~3-4 weeks for Phases 4 and 5
:Next Milestone: codeutil.py DataStore integration and float conversion (Phase 4)
:Production Readiness: Phase 5 completion
