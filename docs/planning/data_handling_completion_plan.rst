==========================================
Data Handling System Completion Plan
==========================================

Executive Summary
-----------------

The lysis project has a well-architected data handling system with clear separation of concerns across four main modules:

- **dataspec.py**: Complete specification definitions (DONE)
- **fileops.py**: File I/O operations (DONE - JSON writer missing)
- **dataconvert.py**: Format conversion logic (90% complete - 43/44 converters implemented)
- **datastore.py**: High-level Python API (Interface only, no implementation)

**Overall Status: ~75% Complete**

**Recent Updates:**

- Fixed HDF5 attribute writing bug (used hard-coded slice)
- Removed debug print statement
- Implemented all microscale data converters with round-trip testing
- Implemented macroscale input generation pipeline (``generate_macroscale_in``)
- Implemented generic macroscale output converters (``convert_structured_grid_fields``,
  ``convert_location_snapshot``) using ``functools.partial`` in the converter registry
- Added safe string/object type conversion (``safe_np_string_conversion``)
- Improved ``safe_np_int_conversion`` to handle float arrays with integer values
- Added ``NUMPY_SAVETXT_FORMATS`` dictionary and ``get_savetxt_format()`` to constants
- Updated ``_write_file_text()`` to use format lookup from constants
- Improved ``read_data_collection()`` parameter handling
- Established conversion routing strategy: all specs route through v1.99.0 <-> v2.0.0
- Created test scripts for microscale and macroscale input conversion pipelines
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
   |   - High-level interface: get/set/append         |
   |   - Status tracking: LOADED, SAVED, FILLED       |
   |   Status: PLACEHOLDER - Methods not implemented  |
   +------------------+-------------------------------+
                      |
   +------------------v-------------------------------+
   |   Data Conversion Layer (dataconvert.py)         |
   |   - Bidirectional v1.99.0 <-> v2.0.0 conversion |
   |   - Safe type casting with validation            |
   |   - Generic structured array converters          |
   |   Status: 90% - Only m_bound stub remains        |
   +------------------+-------------------------------+
                      |
   +------------------v-------------------------------+
   |   Data I/O Layer (fileops.py)                    |
   |   - Format-agnostic read/write with routing      |
   |   - Collection-level multi-file handling         |
   |   - Format string lookup from constants          |
   |   Status: MOSTLY COMPLETE - JSON writer missing  |
   +------------------+-------------------------------+
                      |
   +------------------v-------------------------------+
   |   Specification Layer (dataspec.py)              |
   |   - Schema definitions with versions             |
   |   - Type/shape validation                        |
   |   Status: COMPLETE                               |
   +------------------+-------------------------------+
                      |
   +------------------v-------------------------------+
   |   Storage Backends                               |
   |   HDF5 (done) Text (done) Binary (done)          |
   |   JSON (read only)                               |
   +--------------------------------------------------+


Critical Gaps (Blocking System Use)
-----------------------------------

1. DataStore Implementation (CRITICAL)
++++++++++++++++++++++++++++++++++++++

:File: ``datastore.py``
:Impact: Users cannot use the high-level API
:Task: #1

All core methods are placeholder ``pass`` statements:

- ``__getattr__()`` - Cannot read data by name
- ``__setattr__()`` - Cannot write data by name
- ``status()`` - Cannot query data state
- ``delete()`` - Cannot remove datasets
- ``overwrite()`` - Cannot replace datasets
- ``append()`` - Cannot extend arrays

**Estimated Effort**: 2-3 days

2. Missing Tests (CRITICAL)
++++++++++++++++++++++++++++

:Impact: Cannot verify correctness or catch regressions
:Task: #2

No unit tests exist for:

- fileops.py (read/write functions)
- dataconvert.py (conversion logic)
- datastore.py (API methods)
- dataspec.py (shape parsing, validation)

Manual test scripts exist for microscale and macroscale input conversion pipelines
(see ``src/python/scripts/test_microscale_conversion.py`` and
``test_macroscale_in_conversion.py``), but these are not automated unit tests.

**Estimated Effort**: 3-4 days

3. HDF5 Attribute Writing Bug (HIGH) -- COMPLETED
++++++++++++++++++++++++++++++++++++++++++++++++++

:File: ``fileops.py``
:Task: #3

Fixed. The hard-coded ``[:6]`` slice has been replaced with proper path
construction.

4. Debug Print Statement (HIGH) -- COMPLETED
+++++++++++++++++++++++++++++++++++++++++++++

:File: ``dataconvert.py``
:Task: #4

Fixed. Debug print statement removed.


Important Gaps (Reduces Robustness)
-----------------------------------

5. JSON Writer Not Implemented (MEDIUM)
+++++++++++++++++++++++++++++++++++++++

:File: ``fileops.py``
:Impact: Cannot write parameters to JSON via standard interface
:Task: #5

``_write_file_json()`` is mapped to ``_not_implemented()`` in the ``data_writers``
registry.

**Estimated Effort**: 2-3 hours

6. Float Type Conversion Missing (MEDIUM)
++++++++++++++++++++++++++++++++++++++++++

:File: ``dataconvert.py``
:Impact: Cannot convert float datasets between specifications
:Task: #6

Currently raises ``NotImplementedError`` for float-to-float dtype conversions
(e.g. float32 <-> float64). Note that float-to-int conversion *is* now handled
by ``safe_np_int_conversion``.

**Estimated Effort**: 3-4 hours

7. Data Validation Missing (MEDIUM)
+++++++++++++++++++++++++++++++++++++

:File: ``datastore.py``
:Impact: Invalid data can be silently stored
:Task: #7

Need validation in:

- ``__init__()`` - Check existing HDF5 matches spec
- ``import_fortran_micro_data()`` - Validate imported data
- New ``validate()`` method for comprehensive checks

**Estimated Effort**: 1-2 days

8. Fiber Degrade Conversion (MEDIUM) -- COMPLETED
++++++++++++++++++++++++++++++++++++++++++++++++++

:File: ``dataconvert.py``
:Tasks: #14, #17

Implemented via the generic ``convert_structured_grid_fields()`` helper. This
function handles bidirectional conversion of structured arrays between 1D Fortran
grid indices and 2D row/rank coordinates. It also validates that all fields in both
dtypes are handled, raising ``NotImplementedError`` for unrecognized fields.

9. m_bound Conversion (MEDIUM)
++++++++++++++++++++++++++++++

:File: ``dataconvert.py``
:Impact: Cannot write m_bound.dat when converting v2.0.0 to v1.99.0
:Task: #18

The only remaining converter stub. The ``m_bound`` dataset (molecule binding status
at each snapshot) is not stored directly in v2.0.0 and would need to be
reconstructed from ``tpa_bind_events`` and ``snapshot_time`` by replaying the event
log.

**Estimated Effort**: 4-6 hours

10. Macroscale Input Data Converters (MEDIUM)
+++++++++++++++++++++++++++++++++++++++++++++

:File: ``dataconvert.py``
:Impact: Converters exist but are not tested as standalone conversions
:Task: #15

All 5 macroscale input datasets have implemented converters in both directions:

- bin_edge_proportions <-> tPAleave
- bin_edge_tpa_leaving_time <-> tsectPA
- binned_fiber_degrade_time <-> lysismat
- binned_fiber_degraded <-> lenlysisvect
- edge_grid_neighbors <-> neighbors

The ``generate_macroscale_in()`` function (which generates macroscale input from
microscale output) is also implemented and tested. The remaining work is to verify
that the standalone converters (``convert_data`` calls for individual datasets)
handle all edge cases.

**Estimated Effort**: 1 day (verification and edge cases)

11. Macroscale Output Data Converters (MEDIUM) -- MOSTLY COMPLETED
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

:File: ``dataconvert.py``
:Tasks: #16, #17

All macroscale output converters are implemented except ``m_bound`` (see item 9):

- snapshot_time <-> tsave/Nsave -- direct mapping
- tpa_bind_events <-> m_bind_t -- ``convert_structured_grid_fields``
- tpa_location_snapshot <-> m_loc -- ``convert_location_snapshot``
- tpa_transit_time <-> mfpt -- direct mapping
- fiber_degrade_time <-> f_deg_list -- ``convert_structured_grid_fields``
- macro_log -- direct mapping

A test script for these converters still needs to be written (Task #19).

12. Macroscale Output Test Script (MEDIUM)
++++++++++++++++++++++++++++++++++++++++++

:Task: #19

Create a test script (similar to ``test_microscale_conversion.py`` and
``test_macroscale_in_conversion.py``) that reads Fortran v1.99.0 macroscale output
files, converts to v2.0.0, converts back to v1.99.0, writes to disk, and compares
with originals.

**Estimated Effort**: 4-6 hours


Quality Improvements (Nice to Have)
-----------------------------------

13. Extract Hard-Coded Constants (LOW)
++++++++++++++++++++++++++++++++++++++

:Files: ``dataconvert.py``
:Impact: Improves maintainability
:Task: #8

Hard-coded ``n_bins = 100`` should be configurable constant.

**Estimated Effort**: 1-2 hours

14. Optimize HDF5 Chunking (LOW)
+++++++++++++++++++++++++++++++++

:File: ``fileops.py``
:Impact: Better I/O performance for large datasets
:Task: #9

Currently auto-calculated by h5py, may not be optimal.

**Estimated Effort**: 1-2 days (requires analysis)

15. Complete Documentation (LOW)
++++++++++++++++++++++++++++++++

:Files: All
:Impact: Easier to use and maintain
:Task: #10

Many docstrings have "_description_" placeholders. Documentation for
``np_macroscale.py`` has been substantially updated. RST documentation in
``docs/usage/`` has known syntax issues that need fixing.

**Estimated Effort**: 2-3 days


Long-Term Enhancements
----------------------

16. Specification Versioning in Files (FUTURE)
+++++++++++++++++++++++++++++++++++++++++++++++

:Task: #11

Store version metadata in HDF5 root attributes for:

- Version mismatch detection
- Automatic migration
- Audit trails

**Estimated Effort**: 1 day

17. Integration Tests (FUTURE)
++++++++++++++++++++++++++++++

:Task: #12

End-to-end tests covering:

- Full v1.99.0 -> v2.0.0 -> v1.99.0 round trip
- Multiple simulations
- Edge cases and error conditions

The existing manual test scripts provide a starting point for this work.

**Estimated Effort**: 2-3 days

18. CLI Tools (FUTURE)
++++++++++++++++++++++

:Task: #13

Command-line interface for:

- Data format conversion
- Validation and inspection
- Export to Fortran format

**Estimated Effort**: 2-3 days


Implementation Priority
-----------------------

Phase 1: Critical (2-3 weeks)
++++++++++++++++++++++++++++++

Must complete before system is usable:

#. **Task #3**: Fix HDF5 attribute writing bug (2-4 hours) -- COMPLETED
#. **Task #4**: Remove debug print statement (5 minutes) -- COMPLETED
#. **Task #1**: Implement DataStore core methods (2-3 days)
#. **Task #2**: Create comprehensive unit tests (3-4 days)

**Deliverable**: Functional data handling system with basic testing

Phase 2: Important (2-3 weeks)
+++++++++++++++++++++++++++++++

Improves robustness and completeness:

5. **Task #14**: Complete convert_fiber_degrade_time -- COMPLETED
#. **Task #16**: Implement macroscale output data converters -- COMPLETED (except m_bound)
#. **Task #17**: Implement generic macroscale output converters -- COMPLETED
#. **Task #5**: Implement JSON file writer (2-3 hours)
#. **Task #6**: Add safe float type conversion (3-4 hours)
#. **Task #7**: Add data validation to DataStore (1-2 days)
#. **Task #15**: Verify macroscale input data converters (1 day)
#. **Task #18**: Implement m_bound conversion (4-6 hours)
#. **Task #19**: Create macroscale output test script (4-6 hours)

**Deliverable**: Robust system handling all data types and complete format conversion

Phase 3: Quality (1 week)
++++++++++++++++++++++++++

Polish and optimization:

14. **Task #8**: Extract hard-coded constants (1-2 hours)
#.  **Task #9**: Optimize HDF5 chunk sizes (1-2 days)
#.  **Task #10**: Complete documentation and docstrings (2-3 days)

**Deliverable**: Well-documented, optimized system

Phase 4: Long-Term (2-3 weeks)
+++++++++++++++++++++++++++++++

Advanced features and tooling:

17. **Task #11**: Add specification versioning (1 day)
#.  **Task #12**: Create integration tests (2-3 days)
#.  **Task #13**: Create CLI for data operations (2-3 days)

**Deliverable**: Production-ready system with tools


Current Capabilities (Working Now)
-----------------------------------

Ready for Use
+++++++++++++

#. Reading Fortran v1.99.0 format (text, binary, JSON)
#. Reading HDF5 v2.0.0 format
#. Writing to HDF5 v2.0.0 format
#. Writing to text and binary v1.99.0 format
#. Microscale data conversion v1.99.0 <-> v2.0.0 (all datasets, round-trip tested)
#. Macroscale input generation from microscale output (round-trip tested)
#. Macroscale input data conversion v1.99.0 <-> v2.0.0 (all datasets)
#. Macroscale output data conversion v1.99.0 <-> v2.0.0 (all datasets except m_bound)
#. Safe type conversion for integers, booleans, and strings/objects
#. Parameter loading and merging
#. Dynamic shape resolution from parameters
#. Multi-simulation file handling
#. Format string lookup for text file output (``NUMPY_SAVETXT_FORMATS``)

Partial / Needs Work
++++++++++++++++++++

#. DataStore API (interface exists, no implementation)
#. m_bound converter (requires event log reconstruction)
#. Float-to-float type conversion (not implemented)
#. JSON writing (read only currently)
#. Data validation (no checks performed)
#. Automated unit tests (manual test scripts exist)

Not Implemented
+++++++++++++++

#. HDF5 chunk optimization
#. Specification versioning in files
#. Streaming I/O for large datasets
#. Parallel I/O (MPI support)
#. Command-line tools


Design Strengths
----------------

#. **Separation of Concerns**: Clear module boundaries
#. **Format Agnostic**: Easy to add new storage formats
#. **Version Management**: Tag system for forward compatibility
#. **Immutable Specs**: Prevents runtime changes to schemas
#. **Parameter Integration**: Dynamic shapes computed at I/O time
#. **Dispatcher Pattern**: Clean routing to format-specific functions
#. **Generic Converters**: Reusable ``convert_structured_grid_fields`` and
   ``convert_location_snapshot`` with field validation
#. **Conversion Routing**: All specs route through v1.99.0 <-> v2.0.0 pair,
   extensible to future versions

Design Weaknesses
-----------------

#. **Incomplete DataStore**: High-level API is placeholder only
#. **No Automated Tests**: Only manual test scripts exist
#. **No Validation Pipeline**: Can store invalid data
#. **Tight Path Coupling**: Format strings in spec definitions
#. **No Version Recording**: HDF5 files don't record spec version
#. **Hard-Coded Values**: Magic numbers should be constants (e.g. ``n_bins = 100``)


Risk Assessment
---------------

.. list-table::
   :header-rows: 1
   :widths: 30 12 12 46

   * - Risk
     - Likelihood
     - Impact
     - Mitigation
   * - DataStore bugs due to no implementation
     - HIGH
     - HIGH
     - Complete implementation + tests (Phase 1)
   * - m_bound conversion blocks full round-trip
     - MEDIUM
     - MEDIUM
     - Implement event log reconstruction (Phase 2) - Task #18
   * - No automated tests means regression bugs
     - HIGH
     - HIGH
     - Create comprehensive test suite (Phase 1) - Task #2
   * - Float conversion blocks workflows
     - LOW
     - MEDIUM
     - Implement safe conversion (Phase 2) - Task #6
   * - Poor chunking hurts performance
     - LOW
     - LOW
     - Profile and optimize (Phase 3) - Task #9
   * - Version mismatches cause confusion
     - MEDIUM
     - MEDIUM
     - Add version metadata (Phase 4) - Task #11


Success Criteria
----------------

Phase 1 Success
+++++++++++++++

- [x] HDF5 attribute bug fixed and verified
- [x] No debug print statements in code
- [ ] All DataStore methods implemented and tested
- [ ] Unit test coverage > 80% for all modules
- [ ] Can read/write full simulation datasets
- [x] Round-trip conversion preserves data integrity (microscale, macroscale input)

Phase 2 Success
+++++++++++++++

- [x] Integer and boolean types convert correctly
- [ ] Float-to-float type conversion works
- [x] String and object types convert correctly
- [ ] JSON reading and writing both work
- [ ] Data validation catches spec violations
- [x] Fiber degrade conversion handles all cases
- [x] All macroscale input data converters implemented
- [x] Macroscale output data converters implemented (except m_bound)
- [ ] m_bound conversion implemented
- [ ] Full bidirectional conversion working for all datasets
- [ ] Macroscale output conversion round-trip tested

Phase 3 Success
+++++++++++++++

- [ ] All docstrings complete with examples
- [ ] User guide documentation exists
- [ ] HDF5 I/O performance benchmarked
- [ ] Constants extracted from code
- [ ] Code review passes with no major issues

Phase 4 Success
+++++++++++++++

- [ ] Integration tests cover full workflows
- [ ] CLI tools work for common operations
- [ ] Version metadata in all HDF5 files
- [ ] Migration path documented for future versions
- [ ] Performance benchmarks meet targets


Recommended Next Steps
----------------------

#. **Immediate**:

   - Implement m_bound conversion (Task #18)
   - Create macroscale output test script (Task #19)

#. **Short-term**:

   - Start implementing DataStore methods (Task #1)
   - Begin unit test creation (Task #2)

#. **Medium-term**:

   - Complete DataStore implementation
   - Achieve basic test coverage
   - Implement JSON writer (Task #5)
   - Add float type conversion (Task #6)

#. **Longer-term**:

   - Verify macroscale input standalone converters (Task #15)
   - Add data validation (Task #7)
   - Continue with Phase 3 and Phase 4 tasks


Conclusion
----------

The data handling system has a solid architectural foundation with clear specifications
and working I/O for all major formats. Significant progress has been made since the
initial plan:

#. All microscale data converters are implemented and round-trip tested
#. Macroscale input generation pipeline is complete and tested
#. 43 of 44 data converters are implemented (only m_bound remains)
#. Generic converter infrastructure (``convert_structured_grid_fields``,
   ``convert_location_snapshot``) provides reusable, validated conversion for
   structured arrays with grid location fields
#. Safe type conversion covers integers, booleans, strings, objects, and
   float-to-int
#. Conversion routing strategy established for future spec version extensibility

The primary remaining gaps are:

#. Incomplete DataStore API implementation, which blocks high-level use
#. Missing automated unit tests (manual test scripts exist)
#. The m_bound converter (requires event log reconstruction)

Completing Phase 1 (2-3 weeks) will result in a functional system suitable for
production use with basic I/O. The remaining Phase 2 work is lighter than originally
estimated, as most converters are now complete.

:Estimated Total Effort: 6-10 weeks for all phases
:Minimum Viable System: 2-3 weeks (Phase 1 only)
:Full Conversion Support: 3-5 weeks (Phases 1-2)
