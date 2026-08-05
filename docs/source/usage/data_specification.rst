=========================
Data Specifications
=========================

*Note: This document is written with NumPy indexing, which is zero-indexed.
This is especially noteworthy since the Fortran code is one-indexed.*

*Note: Grid Locations in this specification are stored in a two-dimension, 0-indexed system.
For more information, see the documentation of the* :class:`~lysis.geometry.edge_grid.EdgeGrid` *class.*

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
v2.0.0 (First HDF5-based specification)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Files
+++++++++++

``YYYY-MM-DD-hhmm.h5``
  HDF5 file containing all data for a run

Group Structure
--------------------

``micro_data``
  Contains microscale output datasets.

``macro_data``
  Contains a group for each macroscale simulation's datasets, named sim_00, sim_01, sim_02, ...

``log_files``
  Contains log files for all executions of code.  The binary's own code-output
  logs live at the top level (``micro_log``, ``macro_log__sim_XX``);
  scheduler-captured **worker** logs from Slurm-dispatched runs live in the
  ``dispatcher`` subgroup (optional — absent for direct execution).  The Slurm
  *master* log is left on disk in ``<run dir>/.slurm/`` and is not ingested.
  On a successful run the on-disk worker ``.out`` files are removed with the
  staging directory (the HDF5 becomes the record); on failure they are kept on
  disk, uningested.

Provenance attributes
---------------------

Each per-scale params group (``micro_data`` / ``macro_data``) carries
provenance attributes recording *when*, *where*, and *with what code* the
data was produced.  Three families are stamped, plus the root
``dataspec_version`` attribute.  All three ``*_dirty`` fields use the same
3-state string — ``"clean"``, ``"dirty"``, or ``"unknown"`` (git unavailable
or repo root not found) — so the unknown case is never collapsed into clean.

``init_*`` — the ``src/lysis/`` Python layer at **init** time (stamped by
``init-experiment`` / ``init-macroscale``):

  ``init_version``, ``init_dirty``, ``init_timestamp``, ``init_hostname``

``pipeline_*`` — the ``src/lysis/`` Python **pipeline** that orchestrated
init/import at HDF5-import time (stamped by ``run-micro`` / ``run-macro``).
This is the wrapper that ran the simulation, *not* the engine that executed
it:

  ``pipeline_version``, ``pipeline_dirty``, ``pipeline_timestamp``,
  ``pipeline_hostname``

``backend_*`` — the simulation **engine** that actually produced the output.
For Fortran runs these come from the binary's ``--version`` output:

  ``backend_type``
    ``"fortran"`` (or ``"python"`` for the future Python backend).
  ``backend_commit``, ``backend_dirty``, ``backend_compiler``
    The ``src/fortran/`` commit, dirty bit, and compiler string the binary
    was built from/with.
  ``backend_historical``
    ``True`` only when the binary was rebuilt from an older commit via
    ``run-{micro,macro} --fortran-commit <ref>``; absent otherwise.  The SHA
    itself is in ``backend_commit``.
  ``stale_backend_override``
    ``True`` only when a binary↔source staleness mismatch was explicitly
    overridden (``--allow-stale-binary`` / ``LYSIS_ALLOW_STALE_BINARY=1``);
    absent otherwise.

Microscale datasets
++++++++++++++++++++++++++++++++

``pli_first_time``
  Time to first plasmin

  :Data Type: 
    NumPy 64-bit float (``f8``)
  :Units:
    seconds
  :Dimensions: 
    (``micro_params.simulations``,)
  :Fortran Name:
    ``firstPLi``


``tpa_final_num``
  The number of tPA molecules in the fiber at the end of the simulation.

  :Data Type: 
    NumPy unsigned 8-bit integer (``u1``)
  :Units:
    None
  :Dimensions: 
    (``micro_params.simulations``,)
  :Fortran Name:
    ``lasttPA``

``fiber_degraded``
  Whether or not lysis is complete at the end of the simulation. That
  is, at least ``micro_params.snap_proportion`` of binding doublets were degraded.
  1 for yes, 0 for no.

  :Data Type: 
    NumPy boolean (``?``)
  :Units:
    None
  :Dimensions: 
    (``micro_params.simulations``,)
  :Fortran Name:
    ``lysis_complete``

``sim_final_time``
  The amount of time elapsed in each simulation. If lysis completed,
  this is the time at which that occurred.

  :Data Type: 
    NumPy 64-bit float (``f8``)
  :Units:
    seconds
  :Dimensions: 
    (``micro_params.simulations``,)
  :Fortran Name:
    ``lysis.dat``

``pli_generated_num``
  The number of Plasmin molecules generated in the fiber
  by the end of the simulation.

  :Data Type: 
    NumPy 16-bit unsigned integer (``u2``)
  :Units:
    None
  :Dimensions: 
    (``micro_params.simulations``,)
  :Fortran Name:
    ``PLi``

``tpa_leaving_time``
  The simulation time elapsed when the tPA molecule leaves the system
  or infinity if the simulation ends
  with the tPA molecule still bound. 
  
  :Data Type: 
    NumPy 64-bit float (``f8``)
  :Units:
    seconds
  :Dimensions: 
    (``micro_params.simulations``,)
  :Fortran Name:
    ``tPA_time``

``tpa_unbound_by_pli``
  Whether or not tPA was forced to unbind by plasmin-mediated
  degradation of fibrin.
  1 for yes, 0 for no.

  :Data Type: 
    NumPy boolean (``?``)
  :Units:
    None
  :Dimensions: 
    (``micro_params.simulations``,)
  :Fortran Name:
    ``tPAPLiunbd``

``tpa_unbound_kinetic``
  Whether or not tPA unbinds from the fiber by itself (kinetically).
  1 for yes, 0 for no.

  :Data Type: 
    NumPy boolean (``?``)
  :Units:
    None
  :Dimensions: 
    (``micro_params.simulations``,)
  :Fortran Name:
    ``tPAunbind``

Macro to Micro datasets
+++++++++++++++++++++++

:NOTE: These files should NOT be stored in this version of the data specification, 
       but should be generated when needed from the microscale data.

``bin_edge_proportions``
  The bin boundaries of the Cumulative Density Function for the tPA leaving time distribution.
  That is, ``bin_edge_proportions[i]`` is the fraction of simulations in bins 0 through ``i+1``.
  There are 100 bins evenly distributed along the interval :math:`[0, 1]`.

  :Data Type:
    NumPy 64-bit float (``f8``)
  :Units:
    None
  :Dimensions:
    (101,)
  :Fortran Name:
    ``tPAleave.dat``

``bin_edge_tpa_leaving_time``
  The times at the bin boundaries of the tPA leaving time distribution.
  That is, i% of all microscale simulations had their tPA leave in a time <= tsectPA[i].

  :Data Type:
    NumPy 64-bit float (``f8``)
  :Units:
    seconds
  :Dimensions:
    (101,)
  :Fortran Name:
    ``tsectPA.dat``

``binned_fiber_degrade_time``
  The fiber degrade time of all microscale simulations, binned by their tPA leaving time,
  sorted by their degradation time.
  That is, column ``i`` of this matrix contains the fiber degrade time of all simulations whose
  tPA leaving time falls in the bin i of the tPA leaving time distribution.
  Put another way, if a simulation has a tPA leaving time in the interval 
  [``bin_edge_tpa_leaving_time[i]``, ``bin_edge_tpa_leaving_time[i+1]``], 
  then that simulation's fiber degrade time will be found in column ``i``.
  In that column, the fiber degrade times are sorted in increasing order.

  All simulations where degradation did not occur are assigned a value of infinity (numpy.inf).

  :Data Type:
    NumPy 64-bit float (``f8``)
  :Units:
    seconds
  :Dimensions:
    (``micro_params.simulations`` // 100, 100)
  :Fortran Name:
    ``lysismat.dat``

``binned_fiber_degraded``
  The number of simulations in a tPA leaving time bin, where full lysis of the fiber occurs.
  That is, ``binned_fiber_degraded[i]`` is the 0-indexed location of the first ``numpy.inf``
  entry in column ``i`` of ``binned_fiber_degrade_time``.  Because each column is sorted in
  increasing order and ``numpy.inf`` marks "no lysis", every entry above that location is a
  genuine fiber degrade time, so the 0-indexed location *is* the count of degraded
  simulations in that bin.

  This is **not** the same convention as its Fortran counterpart ``lenlysisvect.dat``, which
  stores the equivalent location 1-indexed.  The two differ by exactly one::

      binned_fiber_degraded[i] == lenlysisvect[i] - 1

  :Data Type:
    NumPy 16-bit unsigned integer (``u2``)
  :Units:
    None
  :Dimensions:
    (100,)
  :Fortran Name:
    ``lenlysisvect.dat``

``edge_grid_neighbors``
  This file contains the fortran (1-D), 0-indexed location of the edges neighboring each edge in the
  edge grid.
  That is, edge_grid_neighbors[i, j] is the index of the jth neighbor of edge i
  
  *For more detail on how the edge grid co-ordinates are defined in both Fortran and Python,
  see the documentation of the* :class:`~lysis.geometry.edge_grid.EdgeGrid` *class.*

  :Data Type:
    NumPy 32-bit unsigned integer (``u4``)
  :Units:
    None
  :Dimensions:
    (``macro_params.total_edges``, 8)
  :Fortran Name:
    ``neighbors.dat``

Macroscale datasets
+++++++++++++++++++


``fiber_degrade_time``
  A list of updates to the degrade time of any fiber in the model.
  Each row represents one update to one fiber and contains 4 entries:

  #. The simulation time elapsed when the event occurred.
  #. The python location row index of the fiber on which the event occurred.
  #. The python location rank index of the fiber on which the event occurred.
  #. The new degrade time for the fiber. (the time at which the fiber will degrade if no further updates occur).

  :Data Type: 
    (NumPy float64 (``f8``), NumPy uint32 (``u4``), NumPy uint32 (``u4``), NumPy float64 (``f8``))
  :Units:
    (seconds, None, None, seconds)
  :Dimensions: 
    (number of binding events, 4)
  :Fortran Name:
    ``f_deg_list``

``tpa_bind_events``
  A list of tPA molecule bindings and their associated metrics.
  Each row represents one update to one tPA molecule and contains 5 entries:

  #. The simulation time elapsed when the event occurred.
  #. The index of the tPA molecule which was involved in the event.
  #. The new status of the tPA molecule (see below for details).
  #. The python location (fiber index) row at which the event occurred. 
  #. The python location (fiber index) rank at which the event occurred. 

  Valid molecule statuses are:

  0. Unbound. Not bound to any fibrin.
  1. Bound. Bound to an intact fiber in the lattice.
  2. Macro-unbound (unbinding by degradation). Still attached to a large fibrin
     fragment that has separated from the lattice because the fiber was fully degraded.
  3. Micro-unbound (forced unbinding). Still attached to a small fibrin fragment
     that has separated from the lattice due to plasmin-mediated degradation of the
     binding site.

  :Data Type: 
    (NumPy float64 (``f8``), NumPy uint64 (``u8``), NumPy uint8 (``u1``), NumPy uint32 (``u4``), NumPy uint32 (``u4``))
  :Units:
    (seconds, None, None, None, None)
  :Dimensions: 
    (number of binding events, 5) 
  :Fortran Name:
    ``m_bind_t``

``tpa_location_snapshot``
  An array, giving the location (fiber index) row and rank of each molecule at point when a save is 
  made. That is, ``tpa_location_snapshot[i, :, j]`` is the pair of location coordinates for tPA 
  molecule ``i`` when snapshot ``j`` is recorded.

  :Data Type: 
    NumPy uint32 (``u4``)
  :Units:
    None 
  :Dimensions: 
    (number of molecules, 2, number of snapshots)
  :Fortran Name:
    ``m_loc``

``tpa_transit_time``
  The simulation time elapsed when each tPA molecule reached the back row of the 
  fiber grid for the first time.

  :Data Type: 
    NumPy float64 (``f8``)
  :Units:
    seconds
  :Dimensions: 
    (number of molecules,)
  :Fortran Name:
    ``mfpt``

``snapshot_time``
  The simulation time elapsed at the point when each snapshot is recorded.

  :Data Type: 
    NumPy float64 (``f8``)
  :Units:
    seconds
  :Dimensions: 
    (number of snapshots,)
  :Fortran Name:
    ``tsave``

Log datasets
++++++++++++++++

``micro_log``
  The log from the microscale execution, stored one line per row.

  :Data Type:
    String (``S``)
  :Dimensions:
    (log events,)

``macro_log__sim_XX``
  The log from the macroscale execution of simulation XX, stored one line per row.

  :Data Type:
    String (``S``)
  :Dimensions:
    (log events,)

``dispatcher/micro_dispatcher_log``
  The scheduler-captured stdout/stderr of the microscale **Slurm worker**
  job(s) — the identifying header, ``set -x`` trace, and Fortran stdout —
  concatenated across array tasks into one combined log, stored one line per
  row.  Complements ``micro_log`` (the binary's own code-output log).  Present
  only for Slurm-dispatched runs; ``optional`` and absent for direct
  execution.

  :Data Type:
    String (``S``)
  :Dimensions:
    (log events,)

``dispatcher/macro_dispatcher_log__sim_XX``
  As above, for macroscale simulation XX's array task (one dataset per
  simulation, since a macro array task maps 1:1 to a simulation).  ``optional``
  and absent for direct execution.

  :Data Type:
    String (``S``)
  :Dimensions:
    (log events,)


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
v1.99.0 (Last Fortran-based specification)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*Note: Grid Locations here are stored in a single-dimension, 1-indexed system.
For more information, see the documentation of the* :class:`~lysis.geometry.edge_grid.EdgeGrid` *class.*

Folder root
++++++++++++++++
Purpose:
  Stores all files related to a run

Name:
  First 13 characters are the experiment code, YYYY-MM-DD-hh.
  Last two characters are a two-digit number giving the run's position
  in the experiment

Run files
++++++++++++++++

``README.rst``
  Brief explanation of the experiment and the source of the data

``params.json``
  The parameters of the run (scenario & mechanism combined)

``job.slurm``
  The Slurm script that executed the code

``job.slurm-XXXXXXXX.out``
  The output of the slurm script that executed
  the code. XXXXXXXX is the job number with higher numbers indicating
  later runs

Microscale files
++++++++++++++++

``micro.f90``
  Fortran source code for the microscale model.
  This is the code that was used in this run.

  :File Type: 
    Plain text

``micro.txt``
  Output log for microscale simulation.

  :Data Type: 
    Plain text

``firstPLi.dat``
  Time to first plasmin
  
  :File Type:
    Binary
  :Data Type: 
    double precision (``f8``)
  :Units:
    seconds
  :Dimensions: 
    (``simulations``,)
  :Fortran Name:
    ``firstPLi``


``lasttPA.dat``
  The number of tPA molecules in the fiber at the end of the simulation.

  :File Type:
    Binary
  :Data Type: 
    integer (``i4``)
  :Units:
    None
  :Dimensions: 
    (``simulations``,)
  :Fortran Name:
    ``ltPA = [tPA(count)]``

``lyscomplete.dat``
  Whether or not lysis is complete at the end of the simulation. That
  is, at least ``snap_proportion`` of binding doublets were degraded.
  1 for yes, 0 for no.

  :File Type:
    Binary
  :Data Type: 
    integer (``u4``)
  :Units:
    None
  :Dimensions: 
    (``simulations``,)
  :Fortran Name:
    ``lysiscomplete``

``lysis.dat``
  The amount of time elapsed in each simulation. If lysis completed,
  this is the time at which that occurred.

  :File Type:
    Binary
  :Data Type: 
    double precision (``f8``)
  :Units:
    seconds
  :Dimensions: 
    (``simulations``,)
  :Fortran Name:
    ``lysis_time = [tvals(count)]``

``PLi.dat``
  The number of Plasmin molecules generated in the fiber
  by the end of the simulation.

  :File Type:
    Binary
  :Data Type: 
    integer (``i4``)
  :Units:
    None
  :Dimensions: 
    (``simulations``,)
  :Fortran Name:
    ``Plasmin = [PLi(count)]``

``tPA_time.dat``
  The simulation time elapsed when the tPA molecule leaves the system
  or infinity if the simulation ends
  with the tPA molecule still bound. 
  
  :File Type:
    Binary
  :Data Type: 
    double precision (``f8``)
  :Units:
    seconds
  :Dimensions: 
    (``simulations``,)
  :Fortran Name:
    ``tPA_time = [tvals(loca)]``

``tPAPLiunbd.dat``
  Whether or not tPA was forced to unbind by plasmin-mediated
  degradation of fibrin.
  1 for yes, 0 for no.

  :File Type:
    Binary
  :Data Type: 
    integer (``i4``)
  :Units:
    None
  :Dimensions: 
    (``simulations``,)
  :Fortran Name:
    ``tPAPLiunbd = [reaction3]``

``tPAunbind.dat``
  Whether or not tPA unbinds from the fiber by itself (kinetically).
  1 for yes, 0 for no.

  :File Type:
    Binary
  :Data Type: 
    integer (``i4``)
  :Units:
    None
  :Dimensions: 
    (``simulations``,)
  :Fortran Name:
    ``tPAunbind = [reaction2]``

Macro to Micro files
++++++++++++++++++++

:NOTE: These files are needed to input the data from the microscale model into the Fortran 
       macroscale code. Historically these were generated with the Matlab ``micro_to_macro`` code.

*We continue to use 0-indexing for arrays here to be consistent with the rest of this document, 
even though these files are only used in Fortran which is 1-indexed.*

``tPAleave.dat``
  The bin boundaries of the Cumulative Density Function for the tPA leaving time distribution.
  There are 100 bins evenly distributed along the interval :math:`[0, 1]`.

  :File Type:
    Space-delimited Text
  :Data Type:
    double precision (``f8``)
  :Units:
    None
  :Dimensions:
    (101,)

``tsectPA.dat``
  The times at the bin boundaries of the tPA leaving time distribution.
  That is, i% of all microscale simulations had their tPA leave in a time <= tsectPA[i].

  :File Type:
    Space-delimited Text
  :Data Type:
    double precision (``f8``)
  :Units:
    seconds
  :Dimensions:
    (101,)

``lysismat.dat``
  The fiber degrade time of all microscale simulations, binned by their tPA leaving time,
  sorted by their degradation time.
  That is, column i of this matrix contains the fiber degrade time of all simulations whose
  tPA leaving time falls in the bin i of the tPA leaving time distribution.
  Put another way, if a simulation has a tPA leaving time in the interval 
  [``tsectPA[i]``, ``tsectPA[i+1]``], then that simulation's fiber degrade time will be found
  in column i.
  In that column, the fiber degrade times are sorted in increasing order.

  All simulations where degradation did not occur are assigned a value of 6,000 seconds.

  :File Type:
    Space-delimited Text
  :Data Type:
    double precision (``f8``)
  :Units:
    seconds
  :Dimensions:
    (``simulations`` // 100, 100)

``lenlysisvect.dat``
  The 1-indexed location of the first ``6000`` entry in each column of ``lysismat``.
  That is, ``lenlysisvect[i]`` is the row number, counting from one, of the first simulation
  in column ``i`` of ``lysismat`` for which full lysis of the fiber did *not* occur.

  Because each column of ``lysismat`` is sorted in increasing order and ``6000`` marks
  "no lysis", every entry above that row is a genuine fiber degrade time.  The number of
  simulations in bin ``i`` where full lysis occurred is therefore ``lenlysisvect[i] - 1``,
  **not** ``lenlysisvect[i]``.  See ``binned_fiber_degraded`` for the 0-indexed Python
  equivalent, which does hold the count directly.

  .. note::

     The original MATLAB pre-processing (``archive/matlab/micro_to_macro.m``) writes the
     magic value ``999`` for any bin whose column contains no ``6000`` entry at all — that
     is, a bin in which every simulation degraded.  ``999`` is neither a count nor a valid
     row index, and must be special-cased when reading MATLAB-generated files.  Files
     written by the Python pipeline never contain this value.

  :File Type:
    Space-delimited Text
  :Data Type:
    double precision (``f8``)
  :Units:
    None
  :Dimensions:
    (100,)

``neighbors.dat``
  This file contains the 1-indexed location of the fibers neighboring each fiber in the
  Fortran grid.
  This data was historically generated inside the Fortran code, 
  but it is now generated by the Python wrapper.

  *Note: The matrix is (``num``, 8) when generated in Python, 
  stored on disk as a flattened (1-D) array in row-major order, 
  and read into memory in Fortran as a (8, ``num``) array*

  :File Type:
    Line-delimited Text
  :Data Type:
    integer (``i4``)
  :Units:
    None
  :Dimensions:
    (``num`` * 8,)


Macroscale files
++++++++++++++++


Subfolders
----------
:Purpose:
  Stores all files related to an individual macroscale simulation

:Name: 
  A two-digit number giving the macroscale simulation's position
  in the run array

``macro.txt``
  Output log for the macroscale simulation.

  :Data Type: 
    Plain text

``f_deg_list.dat``
  A list of updates to the degrade time of any fiber in the model.
  Each row represents one update to one fiber and contains 3 entries:

  #. The simulation time elapsed when the event occurred.
  #. The fortran location index of the fiber on which the event occurred.
  #. The new degrade time for the fiber. (the time at which the fiber will degrade if no further 
     updates occur).

  :File Type:
    Comma-delimited Text
  :Data Type: 
    (double precision (``f8``), integer (``i4``), double precision (``f8``))
  :Units:
    (seconds, None, seconds)
  :Dimensions: 
    (number of binding events which result in an updated degrade time, 3)
  :Fortran Name:
    ``(t, V(1, j), t_degrade(V(1, j)))``

``m_bind_t.dat``
  A list of tPA molecule bindings and their associated metrics.
  Each row represents one update to one tPA molecule and contains 4 entries:

  #. The simulation time elapsed when the event occurred.
  #. The 1-based index of the tPA molecule which was involved in the event.
  #. The new status of the tPA molecule (see below for details).
  #. The fortran location (fiber index) at which the event occurred. 

  Valid molecule statuses are:

  0. Unbound. Not bound to any fibrin.
  1. Bound. Bound to an intact fiber in the lattice.
  2. Macro-unbound (unbinding by degradation). Still attached to a large fibrin
     fragment that has separated from the lattice because the fiber was fully degraded.
  3. Micro-unbound (forced unbinding). Still attached to a small fibrin fragment
     that has separated from the lattice due to plasmin-mediated degradation of the
     binding site.

  :File Type:
    Comma-delimited Text
  :Data Type: 
    (double precision (``f8``), integer (``i4``), integer (``i4``), integer (``i4``))
  :Units:
    (seconds, None, None, None)
  :Dimensions: 
    (number of binding and unbinding events, 4) 
  :Fortran Name:
    ``(t, j, [0-3], V(1, j))``

``m_loc.dat``
  An array, giving the fortran location (fiber index) of each molecule at point when a save is made.
  That is, ``m_loc[i, j]`` is the fortran location of tPA molecule ``j`` when 
  snapshot ``i`` is recorded.

  :File Type:
    Binary
  :Data Type: 
    integer (``i4``)
  :Units:
    None 
  :Dimensions: 
    (``cNsave``, ``M``)
  :Fortran Name:
    ``V(1, :)``

``m_bound.dat``
  An array, giving the bound/unbound status of each molecule at point when a save is made.
  That is, ``m_bound[i, j]`` is 1 if tPA molecule ``j`` is bound to a fiber when 
  snapshot ``i`` is recorded, 0 else.

  :File Type:
    Binary
  :Data Type: 
    integer (``i4``)
  :Units:
    None 
  :Dimensions: 
    (``cNsave``, ``M``)
  :Fortran Name:
    ``V(2, :)``

``mfpt.dat``
  The simulation time elapsed when each tPA molecule reached the back row of the 
  fiber grid for the first time.

  :File Type:
    Binary
  :Data Type: 
    double precision (``f8``)
  :Units:
    seconds
  :Dimensions: 
    (``M``,)
  :Fortran Name:
    ``mfpt``

``tsave.dat``
  The simulation time elapsed at the point when each snapshot is recorded.

  :File Type:
    Binary
  :Data Type: 
    double precision (``f8``)
  :Units:
    seconds
  :Dimensions: 
    (``cNsave``,)
  :Fortran Name:
    ``t``

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
v1.95.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Identical to v1.99.0 with the following exceptions:

Parameters
++++++++++++++++

- Microscale parameters are parsed from the microscale log file
  (``micro{file_code}.txt``) instead of ``params.json``.
- All parameters are stored as plain ``numpy.float64`` values
  instead of ``pint.Quantity`` objects with units.

``Nsave.dat``
  The number of snapshots recorded in the simulation.

  :File Type:
    Binary
  :Data Type:
    integer (``i4``)
  :Units:
    None
  :Dimensions:
    Scalar ``()``
  :Fortran Name:
    ``Nsavevect = [cNsave]``

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
v1.90.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Identical to v1.95.0 with the following exception in macroscale output:

Macroscale Output
+++++++++++++++++

- ``f_deg_list.dat`` **does not exist** in v1.90.0.
- ``f_deg_time.dat`` replaces it with a binary snapshot array recording the
  scheduled degradation time of every fiber at each snapshot.

``f_deg_time.dat``
  A binary snapshot array of fiber scheduled-degradation times.  Each row
  corresponds to one snapshot; each column corresponds to one fiber edge.
  Edges that have not yet been scheduled for degradation hold the Fortran
  sentinel value ``9.9e100``.

  :File Type:
    Binary
  :Data Type:
    double precision (``f8``)
  :Units:
    seconds
  :Dimensions:
    (``cNsave``, ``num``) — one row per snapshot, one column per fiber edge
    (``num = total_edges``)
  :Fortran Name:
    ``t_degrade``
  :Initial value:
    ``9.9e100`` (edges not yet scheduled for degradation)

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
v1.85.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Identical to v1.90.0 except that **every macroscale output simulation is
concatenated into a single top-level file per dataset**
(``simulations_combined=True``).  v1.85.0 predates the
one-directory-per-simulation convention, so there are no ``{sim:02}/``
subdirectories and no ``_{sim:02}`` filename suffixes.

``microscale_out`` and ``macroscale_in`` are unchanged from v1.90.0.

Run layout
++++++++++

A v1.85.0 run spans **two directories**, with different file codes:

============================  ==========================  =========================
Collection                    Directory                   File code
============================  ==========================  =========================
``microscale_out``            microscale run directory    e.g. ``_PLG2_tPA01_Q2``
``macroscale_in``             macroscale run directory    e.g. ``_PLG2_tPA01_Q2``
``macroscale_out``            macroscale run directory    e.g. ``_PLG2_tPA01_along_Q2``
============================  ==========================  =========================

Both :func:`~lysis.dataio.fileops.read_data_collection` and
:meth:`~lysis.dataio.datastore.DataStore.import_collection` accept a path and
file code per collection, so a v1.85.0 run is imported with one call per
directory.

Recovering per-simulation data
++++++++++++++++++++++++++++++

``Nsave.dat`` is the only record of where one simulation ends and the next
begins.  Simulation ``i`` occupies ``Nsave[i] + 1`` snapshot rows -- the extra
row is the initial state at ``t = 0``.  Datasets partition three different
ways:

- **Snapshot-indexed** (``tsave``, ``f_deg_time``, ``deg``, ``m_loc``,
  ``m_bound``): split at ``cumsum(Nsave + 1)``.
- **Simulation-indexed** (``mfpt``): exactly one row per simulation.
- **The combined log** (``macro{file_code}.txt``): split at lines matching
  ``run number=``.  The shared header preceding the first marker is prepended
  to every simulation's log so each stays self-describing.

All binary datasets are written one contiguous 1-D block at a time -- one per
snapshot, or one per simulation for ``mfpt`` -- so they are read in C order,
not Fortran column-major order.

``Nsave.dat``
  One snapshot count per simulation, rather than the single scalar used by
  v1.90.0 and later.

  :File Type:
    Binary
  :Data Type:
    integer (``i4``)
  :Units:
    None
  :Dimensions:
    (``macro_simulations``,)
  :Fortran Name:
    ``Nsavevect``

``mfpt.dat``
  First-passage times, one row per simulation.  Written once per run, after
  the time loop, so the file is run-major.

  :File Type:
    Binary
  :Data Type:
    double precision (``f8``)
  :Units:
    seconds
  :Dimensions:
    (``macro_simulations``, ``M``) — ``M = total_molecules``
  :Fortran Name:
    ``mfpt``

``f_deg_time.dat``
  As v1.90.0, **but with a different sentinel convention**.  v1.85.0
  initialises the whole ``t_degrade`` vector to zero and never distinguishes
  empty edges, so ``0.0`` is ambiguous:

  - for the first ``empty_edges`` columns it means "empty (ghost) edge";
  - for every other column it means "fibrin edge, no tPA has landed yet".

  v1.90.0 keeps ``0.0`` for empty edges only and marks unscheduled fibrin with
  ``9.9e100``.  Conversion therefore remaps by **index**, never by value:
  a blanket value remap would destroy the empty-edge marker.

  :Initial value:
    ``0.0`` (both empty edges and edges not yet scheduled for degradation)

``deg.dat``
  **v1.85.0 only** (optional).  The degradation *state* of each edge at each
  snapshot, as distinct from the ``f_deg_time`` *schedule*: ``0`` = intact,
  ``-t`` = degraded at time ``t``, ``-1`` = empty (ghost) edge.  v1.90.0
  removed the underlying Fortran ``degrade`` array in favour of ``t_degrade``,
  so this dataset is dropped on conversion.  It carries no independent
  information: it is derivable from ``f_deg_time``, ``tsave`` and
  ``empty_edges``.

  :File Type:
    Binary
  :Data Type:
    double precision (``f8``)
  :Units:
    seconds (negated), or the ``-1`` empty-edge marker
  :Dimensions:
    (``cNsave``, ``num``) — one row per snapshot, one column per fiber edge
  :Fortran Name:
    ``degrade``

``m_bind_t.dat``
  **Does not exist** in v1.85.0 — the Fortran ``open`` statement for it is
  commented out.  Conversion to v1.90.0 reconstructs the event log by
  differencing consecutive ``m_bound`` snapshots and taking each molecule's
  location from the matching ``m_loc`` snapshot.  This is lossy in two ways:
  event times are quantised to the snapshot interval, and v1.85.0 recorded
  only bound/unbound, so ``MICRO_UNBOUND`` and ``MACRO_UNBOUND`` are both
  reported as ``UNBOUND``.

Parameters
++++++++++

``params.json`` was written by the pre-package ``Experiment`` class.  It uses
the same legacy key spellings as v1.90.0 (``total_trials``, ``seed``,
``log_lvl``), carries a null ``micro_params``, and adds two top-level
bookkeeping keys — ``experiment_code`` and ``data_filenames`` — that are
dropped on conversion.  ``data_filenames`` in particular is a mapping, so it
would otherwise be mistaken for a parameter section.
