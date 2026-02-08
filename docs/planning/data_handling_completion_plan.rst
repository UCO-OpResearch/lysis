==========================================
Data Handling System Completion Plan
==========================================

Executive Summary
-----------------

The lysis project has a well-architected data handling system with clear separation of concerns across four main modules:

- **dataspec.py**: Complete specification definitions (DONE)
- **fileops.py**: File I/O operations (DONE - minor bugs fixed)
- **dataconvert.py**: Format conversion logic (40% complete - microscale done, macroscale missing)
- **datastore.py**: High-level Python API (Interface only, no implementation)

**Overall Status: 50-60% Complete**

**Recent Updates:**

- Fixed HDF5 attribute writing bug (used hard-coded slice)
- Removed debug print statement
- Added converter stubs for 20+ missing datasets with clear error messages

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
   |   Status: PARTIAL - Microscale done, Macroscale  |
   +------------------+-------------------------------+
                      |
   +------------------v-------------------------------+
   |   Data I/O Layer (fileops.py)                    |
   |   - Format-agnostic read/write with routing      |
   |   - Collection-level multi-file handling         |
   |   Status: MOSTLY COMPLETE - Minor bugs           |
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
   |   HDF5, Text, Binary (done)  JSON (read only)   |
   +--------------------------------------------------+


Critical Gaps (Blocking System Use)
-----------------------------------

1. DataStore Implementation (CRITICAL)
++++++++++++++++++++++++++++++++++++++

:File: ``datastore.py``
:Impact: Users cannot use the high-level API

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

No unit tests exist for:

- fileops.py (read/write functions)
- dataconvert.py (conversion logic)
- datastore.py (API methods)
- dataspec.py (shape parsing, validation)

**Estimated Effort**: 3-4 days

3. HDF5 Attribute Writing Bug (HIGH)
+++++++++++++++++++++++++++++++++++++

:File: ``fileops.py`` line 326
:Impact: Parameters may not be written correctly to HDF5

Suspicious code:

.. code-block:: python

   group_location = spec.data_location.format(...)[:6] + "params"

Hard-coded ``[:6]`` slice is fragile and likely incorrect.

**Estimated Effort**: 2-4 hours

4. Debug Print Statement (HIGH)
++++++++++++++++++++++++++++++++

:File: ``dataconvert.py`` line 164
:Impact: Unwanted console output in production

.. code-block:: python

   print(lysis_time[n_bins])  # Remove this

**Estimated Effort**: 5 minutes


Important Gaps (Reduces Robustness)
-----------------------------------

5. JSON Writer Not Implemented (MEDIUM)
+++++++++++++++++++++++++++++++++++++++

:File: ``fileops.py``
:Impact: Cannot write parameters to JSON via standard interface

``_write_file_json()`` is marked as ``_not_implemented()``.

**Estimated Effort**: 2-3 hours

6. Float Type Conversion Missing (MEDIUM)
++++++++++++++++++++++++++++++++++++++++++

:File: ``dataconvert.py`` line 275
:Impact: Cannot convert float datasets between specifications

Currently raises ``NotImplementedError`` for float dtypes.

**Estimated Effort**: 3-4 hours

7. Data Validation Missing (MEDIUM)
+++++++++++++++++++++++++++++++++++++

:File: ``datastore.py``
:Impact: Invalid data can be silently stored

Need validation in:

- ``__init__()`` - Check existing HDF5 matches spec
- ``import_fortran_micro_data()`` - Validate imported data
- New ``validate()`` method for comprehensive checks

**Estimated Effort**: 1-2 days

8. Incomplete Fiber Degrade Conversion (MEDIUM)
++++++++++++++++++++++++++++++++++++++++++++++++

:File: ``dataconvert.py`` lines 168-200
:Impact: May not handle all edge cases correctly

Function exists but marked as incomplete from notebook origin.

**Estimated Effort**: 4-6 hours

9. Missing Macroscale Input Data Converters (MEDIUM)
++++++++++++++++++++++++++++++++++++++++++++++++++++

:File: ``dataconvert.py``
:Impact: Cannot convert macroscale input data between specifications

Currently stubs that raise NotImplementedError for 5 datasets (bidirectional):

- bin_edge_proportions <-> tPAleave
- bin_edge_tpa_leaving_time <-> tsectPA
- binned_fiber_degrade_time <-> lysismat
- binned_fiber_degraded <-> lenlysisvect
- edge_grid_neighbors <-> neighbors

Some logic may already exist in ``generate_macroscale_in()`` and can be refactored.

**Estimated Effort**: 1-2 days

10. Missing Macroscale Output Data Converters (MEDIUM)
++++++++++++++++++++++++++++++++++++++++++++++++++++++

:File: ``dataconvert.py``
:Impact: Cannot convert macroscale simulation results between specifications

Currently stubs that raise NotImplementedError for 8 datasets (bidirectional):

- snapshot_time <-> tsave/Nsave
- tpa_bind_events <-> m_bind_t
- tpa_location_snapshot <-> m_loc/m_bound
- tpa_transit_time <-> mfpt
- fiber_degrade_time <-> f_deg_list (partial - see Task #14)

See TODO comment about cell 10 of H5-File-Builder.ipynb for tpa_bind_events.

**Estimated Effort**: 2-3 days


Quality Improvements (Nice to Have)
-----------------------------------

11. Extract Hard-Coded Constants (LOW)
++++++++++++++++++++++++++++++++++++++

:Files: ``dataconvert.py``
:Impact: Improves maintainability

Hard-coded ``n_bins = 100`` should be configurable constant.

**Estimated Effort**: 1-2 hours

12. Optimize HDF5 Chunking (LOW)
+++++++++++++++++++++++++++++++++

:File: ``fileops.py`` line 293
:Impact: Better I/O performance for large datasets

Currently auto-calculated by h5py, may not be optimal.

**Estimated Effort**: 1-2 days (requires analysis)

13. Complete Documentation (LOW)
++++++++++++++++++++++++++++++++

:Files: All
:Impact: Easier to use and maintain

Many docstrings have "_description_" placeholders.

**Estimated Effort**: 2-3 days


Long-Term Enhancements
----------------------

14. Specification Versioning in Files (FUTURE)
+++++++++++++++++++++++++++++++++++++++++++++++

Store version metadata in HDF5 root attributes for:

- Version mismatch detection
- Automatic migration
- Audit trails

**Estimated Effort**: 1 day

15. Integration Tests (FUTURE)
++++++++++++++++++++++++++++++

End-to-end tests covering:

- Full v1.99.0 -> v2.0.0 -> v1.99.0 round trip
- Multiple simulations
- Edge cases and error conditions

**Estimated Effort**: 2-3 days

16. CLI Tools (FUTURE)
++++++++++++++++++++++

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

#. **Task #1**: Implement DataStore core methods (2-3 days)
#. **Task #2**: Create comprehensive unit tests (3-4 days)
#. **Task #3**: Fix HDF5 attribute writing bug (2-4 hours) -- COMPLETED
#. **Task #4**: Remove debug print statement (5 minutes) -- COMPLETED

**Deliverable**: Functional data handling system with basic testing

Phase 2: Important (2-3 weeks)
+++++++++++++++++++++++++++++++

Improves robustness and completeness:

5. **Task #5**: Implement JSON file writer (2-3 hours)
#. **Task #6**: Add safe float type conversion (3-4 hours)
#. **Task #7**: Add data validation to DataStore (1-2 days)
#. **Task #8**: Complete convert_fiber_degrade_time (4-6 hours)
#. **Task #15**: Implement macroscale input data converters (1-2 days)
#. **Task #16**: Implement macroscale output data converters (2-3 days)

**Deliverable**: Robust system handling all data types and complete format conversion

Phase 3: Quality (1 week)
++++++++++++++++++++++++++

Polish and optimization:

11. **Task #11**: Extract hard-coded constants (1-2 hours)
#.  **Task #9**: Optimize HDF5 chunk sizes (1-2 days)
#.  **Task #10**: Complete documentation and docstrings (2-3 days)

**Deliverable**: Well-documented, optimized system

Phase 4: Long-Term (2-3 weeks)
+++++++++++++++++++++++++++++++

Advanced features and tooling:

14. **Task #14**: Add specification versioning (1 day)
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
#. Basic data conversion v1.99.0 <-> v2.0.0
#. Parameter loading and merging
#. Dynamic shape resolution from parameters
#. Multi-simulation file handling

Partial / Needs Work
++++++++++++++++++++

#. DataStore API (interface exists, no implementation)
#. Data converters (microscale done, macroscale input not done, macroscale output not done)
#. Float type conversion (not implemented)
#. JSON writing (read only currently)
#. Data validation (no checks performed)
#. Error handling (generic, not specific)

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

Design Weaknesses
-----------------

#. **Incomplete Implementation**: Many placeholder methods
#. **Missing Documentation**: Extensive "_description_" stubs
#. **No Validation Pipeline**: Can store invalid data
#. **Tight Path Coupling**: Format strings in spec definitions
#. **No Version Recording**: HDF5 files don't record spec version
#. **Hard-Coded Values**: Magic numbers should be constants


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
   * - Missing converters block workflows
     - HIGH
     - HIGH
     - Implement all converters (Phase 2) - Tasks #15, #16
   * - HDF5 attr bug corrupts parameters
     - LOW
     - HIGH
     - Fixed in Task #3
   * - Float conversion blocks workflows
     - LOW
     - MEDIUM
     - Implement safe conversion (Phase 2)
   * - No tests means regression bugs
     - HIGH
     - HIGH
     - Create comprehensive test suite (Phase 1)
   * - Poor chunking hurts performance
     - LOW
     - LOW
     - Profile and optimize (Phase 3)
   * - Version mismatches cause confusion
     - MEDIUM
     - MEDIUM
     - Add version metadata (Phase 4)


Success Criteria
----------------

Phase 1 Success
+++++++++++++++

- [x] HDF5 attribute bug fixed and verified
- [x] No debug print statements in code
- [ ] All DataStore methods implemented and tested
- [ ] Unit test coverage > 80% for all modules
- [ ] Can read/write full simulation datasets
- [ ] Round-trip conversion preserves data integrity

Phase 2 Success
+++++++++++++++

- [ ] All data types (int, float, bool) convert correctly
- [ ] JSON reading and writing both work
- [ ] Data validation catches spec violations
- [ ] Fiber degrade conversion handles all cases
- [ ] All macroscale input data converters implemented
- [ ] All macroscale output data converters implemented
- [ ] Full bidirectional conversion working for all datasets
- [ ] Error messages are specific and actionable

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

#. **Immediate (Today)**:

   - Remove debug print (Task #4) -- COMPLETED
   - Fix HDF5 attribute bug (Task #3) -- COMPLETED

#. **This Week**:

   - Start implementing DataStore methods (Task #1)
   - Begin unit test creation (Task #2)

#. **Next Week**:

   - Complete DataStore implementation
   - Achieve basic test coverage
   - Implement JSON writer (Task #5)

#. **Following Weeks**:

   - Complete all Phase 2 data converters (Tasks #8, #15, #16)
   - Continue with Phase 2 and Phase 3 tasks
   - Build out integration tests
   - Plan Phase 4 enhancements


Conclusion
----------

The data handling system has a solid architectural foundation with clear specifications
and working I/O for most formats. The primary gaps are:

#. Incomplete DataStore API implementation, which blocks high-level use
#. Missing data converters for macroscale data (20+ datasets)

**Recent Progress:**

- HDF5 attribute writing bug fixed
- Debug print statement removed
- All missing converter stubs added with clear error messages

Completing Phase 1 (2-3 weeks) will result in a functional system suitable for
production use with basic I/O. Phase 2 (2-3 weeks) adds complete bidirectional
conversion between all data formats, making the system fully functional for all
workflows.

The modular design makes it straightforward to complete each component independently,
and the existing test cases in notebooks can be converted to unit tests relatively
easily.

:Estimated Total Effort: 8-12 weeks for all phases
:Minimum Viable System: 2-3 weeks (Phase 1 only)
:Full Conversion Support: 4-6 weeks (Phases 1-2)
