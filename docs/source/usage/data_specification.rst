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
The same family describes both backends; ``backend_type`` says which one:

  ``backend_type``
    ``"fortran"`` for the compiled binary, or ``"python"`` for the in-process
    NumPy backend (``lysis run-macro --backend python``).
  ``backend_commit``, ``backend_dirty``, ``backend_compiler``
    For the **Fortran** backend these come from the binary's ``--version``
    output: the ``src/fortran/`` commit, dirty bit, and compiler string the
    binary was built from/with.  For the **Python** backend the engine *is*
    the ``src/lysis/`` package, so ``backend_commit`` / ``backend_dirty``
    record that package's git state (mirroring ``pipeline_*``) and
    ``backend_compiler`` holds the interpreter and NumPy version that ran the
    model (e.g. ``"CPython 3.11.15; NumPy 2.4.4"``).
  ``backend_historical``
    Fortran-only.  ``True`` only when the binary was rebuilt from an older
    commit via ``run-{micro,macro} --fortran-commit <ref>``; absent
    otherwise.  The SHA itself is in ``backend_commit``.
  ``stale_backend_override``
    Fortran-only.  ``True`` only when a binary↔source staleness mismatch was
    explicitly overridden (``--allow-stale-binary`` /
    ``LYSIS_ALLOW_STALE_BINARY=1``); absent otherwise.

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
  That is, ``binned_fiber_degraded[i]`` is the 1-indexed location of the first ``6000`` entry in 
  column ``i`` of ``binned_fiber_degrade_time``.

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
  The number of simulations in a tPA leaving time bin, where full lysis of the fiber occurs.
  That is, ``lenlysisvect[i]`` is the 1-indexed location of the first ``6000`` entry in 
  column ``i`` of ``lysismat``.

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
