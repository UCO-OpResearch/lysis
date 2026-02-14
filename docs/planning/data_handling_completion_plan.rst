==========================================
Data Handling System Completion Plan
==========================================

Executive Summary
-----------------

The lysis project has a well-architected data handling system with clear separation of concerns across four main modules:

- **dataspec.py**: Complete specification definitions with ``DataSpec`` wrapper class (DONE)
- **fileops.py**: File I/O operations for all five storage formats (DONE)
- **dataconvert.py**: Format conversion logic (DONE - all 44 converters implemented)
- **datastore.py**: High-level Python API (Read-only access implemented; write access pending)

**Overall Status: ~90% Complete**

**Recent Updates:**

- Implemented ``DataStore`` read-only API with dot-access, collection navigation,
  lazy ``h5py.Dataset`` access, and parameter loading
- Implemented HDF5 specification versioning (``dataspec_version`` root attribute
  validated on every read/write)
- Added ``DataSpec`` wrapper class with ``version`` property and hidden field
  auto-population on child specs
- Implemented all 44 data converters (including ``m_bound`` event-log reconstruction)
- Implemented JSON file writer (``_write_file_json``)
- Created comprehensive unit test suite (226 tests across 6 test files)
- Created ``lysis`` CLI with ``convert`` and ``validate`` commands
- Moved ``run.py`` from ``execution/`` to ``config/`` per ontology (a Run is a
  configuration concept, not an execution concept)
- Removed backward-compatibility ``lysis/util/`` shim package; all imports now
  use canonical module paths
- Resolved circular import between ``geometry`` and ``execution`` packages
- Added safe string/object type conversion (``safe_np_string_conversion``)
- Improved ``safe_np_int_conversion`` to handle float arrays with integer values
- Added ``NUMPY_SAVETXT_FORMATS`` dictionary and ``get_savetxt_format()`` to constants
- Updated ``_write_file_text()`` to use format lookup from constants
- Established conversion routing strategy: all specs route through v1.99.0 <-> v2.0.0
- Created test scripts for microscale, macroscale input, and macroscale output
  conversion pipelines
- Updated data specification documentation


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
   |   - Read-only: dot-access, collections, params   |
   |   - Per-simulation and combined views            |
   |   - HDF5 spec version validation on open         |
   |   Status: READ-ONLY COMPLETE; write access TODO  |
   +------------------+-------------------------------+
                      |
   +------------------v-------------------------------+
   |   Data Conversion Layer (dataconvert.py)         |
   |   - Bidirectional v1.99.0 <-> v2.0.0 conversion |
   |   - Safe type casting with validation            |
   |   - Generic structured array converters          |
   |   - Event-log reconstruction (m_bound)           |
   |   Status: COMPLETE - All 44 converters done      |
   +------------------+-------------------------------+
                      |
   +------------------v-------------------------------+
   |   Data I/O Layer (fileops.py)                    |
   |   - Format-agnostic read/write with routing      |
   |   - Collection-level multi-file handling         |
   |   - Format string lookup from constants          |
   |   Status: COMPLETE - All 5 formats supported     |
   +------------------+-------------------------------+
                      |
   +------------------v-------------------------------+
   |   Specification Layer (dataspec.py)              |
   |   - DataSpec wrapper with version tracking       |
   |   - Schema definitions with versions & tags      |
   |   - Type/shape validation                        |
   |   - Hidden field auto-population                 |
   |   Status: COMPLETE                               |
   +------------------+-------------------------------+
                      |
   +------------------v-------------------------------+
   |   Storage Backends                               |
   |   HDF5 (done) Text (done) Binary (done)          |
   |   JSON (done)  Fortran (done)                    |
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

Covered by 62 unit tests in ``tests/data/test_datastore.py``.

4. Comprehensive Unit Tests -- COMPLETED
+++++++++++++++++++++++++++++++++++++++++

:Task: #2

226 unit tests across 6 test files:

- ``tests/data/test_dataspec.py`` (55 tests) - Spec parsing, validation,
  DataSpec wrapper, hidden fields
- ``tests/data/test_datastore.py`` (62 tests) - DataStore read-only API,
  collections, parameter loading, version validation
- ``tests/data/test_fileops.py`` (44 tests) - Read/write for all formats
- ``tests/data/test_dataconvert.py`` (39 tests) - Conversion logic, safe
  type casting, round-trip integrity
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

- ``lysis convert`` - Data format conversion between spec versions
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


Remaining Work
--------------

12. DataStore Write Access (HIGH)
+++++++++++++++++++++++++++++++++

:File: ``datastore.py``
:Impact: Cannot create new HDF5 files or write simulation data
:Task: #30

The ``DataStore`` currently supports read-only access. Remaining work:

- **Initialization**: Create new HDF5 files with proper structure, write
  ``dataspec_version`` attribute, create collection groups
- **Write access**: ``__setattr__()`` or explicit write methods to store
  datasets, update status tracking (``DataStatus`` enum exists but is unused)
- **Append**: Extend existing datasets (e.g. adding simulations)
- **Delete/overwrite**: Remove or replace datasets

This is the primary remaining gap before the system is fully usable for
simulation workflows.

13. Safe Float Type Conversion (MEDIUM)
+++++++++++++++++++++++++++++++++++++++

:File: ``dataconvert.py``
:Impact: Cannot convert float datasets between specifications
:Task: #6

Currently raises ``NotImplementedError`` for float-to-float dtype conversions
(e.g. float32 <-> float64). Note that float-to-int conversion *is* handled
by ``safe_np_int_conversion``.

14. Extract Hard-Coded Constants (LOW)
++++++++++++++++++++++++++++++++++++++

:File: ``dataconvert.py``
:Impact: Improves maintainability
:Task: #8

Hard-coded ``n_bins = 100`` should be configurable constant.

15. Optimize HDF5 Chunking (LOW)
+++++++++++++++++++++++++++++++++

:File: ``fileops.py``
:Impact: Better I/O performance for large datasets
:Task: #9

Currently auto-calculated by h5py, may not be optimal.

16. Complete Documentation (LOW)
++++++++++++++++++++++++++++++++

:Files: All
:Impact: Easier to use and maintain
:Task: #10

Many docstrings have "_description_" placeholders. Documentation for
``np_macroscale.py`` has been substantially updated. Sphinx API documentation
covers all subpackages. RST documentation in ``docs/source/usage/`` has known
syntax issues that need fixing.

17. Integration Tests (LOW)
+++++++++++++++++++++++++++

:Task: #12

End-to-end tests covering full v1.99.0 -> v2.0.0 -> v1.99.0 round trip with
multiple simulations and edge cases. The existing manual test scripts and 226
unit tests provide a starting point.

18. Update Notebooks (LOW)
++++++++++++++++++++++++++

:Task: #29

Four Jupyter notebooks in ``notebooks/`` still reference the removed
``lysis.util`` package and need to be updated to use canonical module paths.

19. CI Pipeline (LOW)
+++++++++++++++++++++

:Task: #28

Set up continuous integration pipeline to execute the pytest unit test suite
automatically on commits and pull requests.


Implementation Priority
-----------------------

Phase 1: Core Functionality -- COMPLETED
+++++++++++++++++++++++++++++++++++++++++

All items delivered:

#. **Task #3**: Fix HDF5 attribute writing bug -- COMPLETED
#. **Task #4**: Remove debug print statement -- COMPLETED
#. **Task #1**: Implement DataStore read-only API -- COMPLETED
#. **Task #2**: Create comprehensive unit tests (226 tests) -- COMPLETED

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

Phase 3: Write Access & Polish (Current)
+++++++++++++++++++++++++++++++++++++++++

Active work:

13. **Task #30**: Implement DataStore initialization and write access
#.  **Task #6**: Add safe float type conversion
#.  **Task #8**: Extract hard-coded constants
#.  **Task #29**: Update notebooks to remove lysis.util imports

**Deliverable**: Read-write data handling system

Phase 4: Quality & Infrastructure
++++++++++++++++++++++++++++++++++

Polish and optimization:

17. **Task #9**: Optimize HDF5 chunk sizes
#.  **Task #10**: Complete documentation and docstrings
#.  **Task #12**: Create integration tests
#.  **Task #28**: Set up CI pipeline

**Deliverable**: Production-ready, well-documented system with CI


Current Capabilities (Working Now)
-----------------------------------

Ready for Use
+++++++++++++

#. Reading Fortran v1.99.0 format (text, binary, JSON)
#. Reading HDF5 v2.0.0 format
#. Writing to HDF5 v2.0.0 format (with spec version metadata)
#. Writing to text, binary, and JSON v1.99.0 format
#. Reading HDF5 via ``DataStore`` with dot-access, per-simulation views,
   and parameter loading
#. HDF5 specification versioning (``dataspec_version`` attribute, validated
   on open/read/write)
#. ``DataSpec`` wrapper class with version tracking and hidden field propagation
#. Microscale data conversion v1.99.0 <-> v2.0.0 (all datasets, round-trip tested)
#. Macroscale input generation from microscale output (round-trip tested)
#. Macroscale input data conversion v1.99.0 <-> v2.0.0 (all datasets)
#. Macroscale output data conversion v1.99.0 <-> v2.0.0 (all datasets including m_bound)
#. Safe type conversion for integers, booleans, strings/objects, and float-to-int
#. Parameter loading and merging
#. Dynamic shape resolution from parameters
#. Multi-simulation file handling
#. Format string lookup for text file output (``NUMPY_SAVETXT_FORMATS``)
#. CLI commands: ``lysis convert``, ``lysis validate``
#. 226 automated unit tests

Needs Work
++++++++++

#. DataStore write access (initialization, dataset writing, append, delete)
#. Float-to-float type conversion (not implemented)
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
   conversion, and store layers
#. **Format Agnostic**: Easy to add new storage formats; all five backends complete
#. **Version Management**: Tag system for forward compatibility; spec version
   stored in HDF5 files and validated on access
#. **Immutable Specs**: ``DataSpec`` wrapper prevents runtime changes to schemas
   while providing dict-like access and version tracking
#. **Parameter Integration**: Dynamic shapes computed at I/O time
#. **Dispatcher Pattern**: Clean routing to format-specific functions
#. **Generic Converters**: Reusable ``convert_structured_grid_fields`` and
   ``convert_location_snapshot`` with field validation
#. **Conversion Routing**: All specs route through v1.99.0 <-> v2.0.0 pair,
   extensible to future versions
#. **Comprehensive Testing**: 226 unit tests covering all modules
#. **Clean Package Structure**: Canonical imports with no shim layers;
   ``config/`` houses Run per project ontology

Design Weaknesses
-----------------

#. **No Write Access**: ``DataStore`` is read-only; cannot create or modify
   HDF5 files through the high-level API
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
     - MEDIUM
     - HIGH
     - Implement write access with comprehensive tests (Phase 3) - Task #30
   * - Float conversion blocks workflows
     - LOW
     - MEDIUM
     - Implement safe conversion (Phase 3) - Task #6
   * - Poor chunking hurts performance at scale
     - LOW
     - LOW
     - Profile and optimize (Phase 4) - Task #9
   * - No CI means regressions go undetected
     - MEDIUM
     - MEDIUM
     - Set up CI pipeline (Phase 4) - Task #28
   * - Stale notebook code confuses users
     - LOW
     - LOW
     - Update notebook imports (Phase 3) - Task #29


Success Criteria
----------------

Phase 1 Success -- ACHIEVED
++++++++++++++++++++++++++++

- [x] HDF5 attribute bug fixed and verified
- [x] No debug print statements in code
- [x] DataStore read-only methods implemented and tested
- [x] Unit test coverage across all modules (226 tests)
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

Phase 3 Success
+++++++++++++++

- [ ] DataStore can create new HDF5 files with proper structure
- [ ] DataStore can write simulation datasets
- [ ] Float-to-float type conversion works
- [ ] Constants extracted from code
- [ ] Notebooks updated to canonical imports

Phase 4 Success
+++++++++++++++

- [ ] Integration tests cover full workflows
- [ ] All docstrings complete with examples
- [ ] User guide documentation reviewed
- [ ] HDF5 I/O performance benchmarked
- [ ] CI pipeline executes tests on every commit
- [ ] Code review passes with no major issues


Recommended Next Steps
----------------------

#. **Immediate**:

   - Implement DataStore initialization and write access (Task #30)
   - Add safe float type conversion (Task #6)

#. **Short-term**:

   - Extract hard-coded constants (Task #8)
   - Update notebooks to remove ``lysis.util`` imports (Task #29)

#. **Medium-term**:

   - Set up CI pipeline (Task #28)
   - Complete documentation and docstrings (Task #10)

#. **Longer-term**:

   - Create integration tests with sample data (Task #12)
   - Optimize HDF5 chunk sizes (Task #9)


Conclusion
----------

The data handling system has reached a mature state with all core functionality
implemented:

#. All 44 data converters are complete and tested, including ``m_bound``
   event-log reconstruction
#. The ``DataStore`` provides a clean read-only API with dot-access to
   collections, per-simulation views, parameter loading, and spec version
   validation
#. ``DataSpec`` wrapper class provides version-aware specification management
   with hidden field auto-population
#. HDF5 specification versioning ensures version consistency on every
   read and write operation
#. All five storage backends (HDF5, text, binary, JSON, Fortran) support
   both reading and writing
#. 226 unit tests cover all modules with automated testing
#. CLI tools (``lysis convert``, ``lysis validate``) provide command-line
   access to conversion and validation
#. Clean package structure with ``run.py`` in ``config/`` per project ontology
   and no backward-compatibility shim layers

The primary remaining gap is **DataStore write access** (Task #30), which will
enable creating new HDF5 files and writing simulation data through the
high-level API. Once implemented, the system will be fully functional for
end-to-end simulation workflows.

:Remaining Effort: ~3-5 weeks for all remaining phases
:Next Milestone: DataStore write access (Phase 3)
:Production Readiness: Phase 4 completion
