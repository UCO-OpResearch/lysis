=========================
Data Specifications
=========================

*Note: This document is written with NumPy indexing, which is zero-indexed.
This is especially noteworthy since the Fortran code is one-indexed.*

*Note: Grid Locations in this specification are stored in a two-dimension, 0-indexed system. 
For more information, see the documentation of the Python lysis.util.EdgeGrid class (link below)*

https://github.com/UCO-OpResearch/lysis/blob/c1e6b2a92758fb8620d6f5a2223976a1478c3231/src/python/lysis/util/edge_grid.py

^^^^^^^^^^^^^^^^^
v2.0.0 (First HDF5-based specification)
^^^^^^^^^^^^^^^^^

Files
+++++++++++

``YYYY-MM-DD-hhmm.h5``
  HDF5 file containing all data for a run

Group Structure
----------

``micro_group``
  Contains microscale output datasets.

``macro_group``
  Contains a group for each macroscale simulations datasets, numbered 00, 01, 02, ...

``log_files``
  Contains log files for all executions of code

Microscale datasets
++++++++++++++++

``pli_first_time``
  Time to first plasmin

  :Data Type: 
    NumPy 64-bit float (``f8``)
  :Units:
    seconds
  :Dimensions: 
    (``micro_params.simulations``,)
  :Fortran equivalent:
    ``firstPLi``


``tpa_final_num``
  The number of tPA molecules in the fiber at the end of the simulation.

  :Data Type: 
    NumPy unsigned 8-bit integer (``u1``)
  :Units:
    None
  :Dimensions: 
    (``micro_params.simulations``,)
  :Fortran equivalent:
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
  :Fortran equivalent:
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
  :Fortran equivalent:
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
  :Fortran equivalent:
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
  :Fortran equivalent:
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
  :Fortran equivalent:
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
  :Fortran equivalent:
    ``tPAunbind``

Macro to Micro datasets
++++++++++++++++++++

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
  :Fortran equivalent:
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
  :Fortran equivalent:
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
  :Fortran equivalent:
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
  :Fortran equivalent:
    ``lenlysisvect.dat``

  ``edge_grid_neighbors``
  This file contains the fortran (1-D), 0-indexed location of the edges neighboring each edge in the
  edge grid.
  That is, edge_grid_neighbors[i, j] is the index of the jth neighbor of edge i
  
  *For more detail on how the edge grid co-ordinates are defined in both Fortran and Python,
  see the documentation of the Python lysis.util.EdgeGrid class (link below)*

  https://github.com/UCO-OpResearch/lysis/blob/c1e6b2a92758fb8620d6f5a2223976a1478c3231/src/python/lysis/util/edge_grid.py

  :Data Type:
    NumPy 32-bit unsigned integer (``u4``)
  :Units:
    None
  :Dimensions:
    (``macro_params.total_edges``, 8)
  :Fortran equivalent:
    ``neighbors.dat``

Macroscale datasets
++++++++++++++++


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
  :Fortran equivalent:
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

  0. Kinetically Unbound.
  1. Bound to an undegraded fiber.
  2. Bound to a degraded fiber (macro-unbound).
  3. Bound to a fiber degradation product (micro-unbound).

  :Data Type: 
    (NumPy float64 (``f8``), NumPy uint64 (``u8``), NumPy uint8 (``u1``), NumPy uint32 (``u4``), NumPy uint32 (``u4``))
  :Units:
    (seconds, None, None, None, None)
  :Dimensions: 
    (number of binding events, 5) 
  :Fortran equivalent:
    ``m_bind_t``

``tpa_location_snapshot``
  An array, giving the location (fiber index) row and rank of each molecule at point when a save is made.
  That is, ``tpa_location_snapshot[i, :, j]`` is the pair of location coordinates for tPA molecule ``i`` when 
  snapshot ``j`` is recorded.

  :Data Type: 
    NumPy uint32 (``u4``)
  :Units:
    None 
  :Dimensions: 
    (number of molecules, 2, number of snapshots)
  :Fortran equivalent:
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
  :Fortran equivalent:
    ``mfpt``

``snapshot_time``
  The simulation time elapsed at the point when each snapshot is recorded.

  :Data Type: 
    NumPy float64 (``f8``)
  :Units:
    seconds
  :Dimensions: 
    (number of snapshots,)
  :Fortran equivalent:
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


^^^^^^^^^^^^^^^^^^
v1.99.0 (Last Fortran-based specification)
^^^^^^^^^^^^^^^^^^
*Note: Grid Locations here are stored in a single-dimension, 1-indexed system. 
For more information, see the documentation of the Python lysis.util.EdgeGrid class (link below)*

https://github.com/UCO-OpResearch/lysis/blob/c1e6b2a92758fb8620d6f5a2223976a1478c3231/src/python/lysis/util/edge_grid.py

Folder root
-----------
Purpose:
  Stores all files related to a run

Name:
  First 13 characters are the experiment code, YYYY-MM-DD-hh.
  Last two caracters are a two-digit number giving the run's position
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
  :Fortran Variable:
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
  :Fortran Variable:
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
  :Fortran Variable:
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
  :Fortran Variable:
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
  :Fortran Variable:
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
  :Fortran Variable:
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
  :Fortran Variable:
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
  :Fortran name:
    ``tPAunbind = [reaction2]``

Macro to Micro files
++++++++++++++++++++

:NOTE: These files are needed to input the data from the microscale model into the Fortran macroscale code.
Historically these were generated with the Matlab ``micro_to_macro`` code.

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
  in column i
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
    integer (``i4``)
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
  and read intp memory in Fortran as a (8, ``num``) array*

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
  A two-digit number giving the macroscale simulations's position
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
  #. The new degrade time for the fiber. (the time at which the fiber will degrade if no further updates occur).

  :File Type:
    Comma-delimited Text
  :Data Type: 
    (double precision (``f8``), integer (``i4``), double precision (``f8``))
  :Units:
    (seconds, None, seconds)
  :Dimensions: 
    (number of binding events which result in an updated degrade time, 3)
  :Fortran Variable:
    ``(t, V(1, j), t_degrade(V(1, j)))``

``m_bind_t.dat``
  A list of tPA molecule bindings and their associated metrics.
  Each row represents one update to one tPA molecule and contains 4 entries:

  #. The simulation time elapsed when the event occurred.
  #. The 1-based index of the tPA molecule which was involved in the event.
  #. The new status of the tPA molecule (see below for details).
  #. The fortran location (fiber index) at which the event occurred. 

  Valid molecule statuses are:

  0. Kinetically Unbound.
  1. Bound to an undegraded fiber.
  2. Bound to a degraded fiber (macro-unbound).
  3. Bound to a fiber degradation product (micro-unbound).

  :File Type:
    Comma-delimited Text
  :Data Type: 
    (double precision (``f8``), integer (``i4``), integer (``i4``), integer (``i4``))
  :Units:
    (seconds, None, None, None)
  :Dimensions: 
    (number of binding and unbinding events, 4) 
  :Fortran Variable:
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
  :Fortran Variable:
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
  :Fortran Variable:
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
  :Fortran Variable:
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
  :Fortran name:
    ``t``

``Nsave.dat``
  The number of snapshots is recorded in the simulation.

  :File Type:
    Binary
  :Data Type: 
    integer (``i4``)
  :Units:
    None
  :Dimensions: 
    (``simulations``,)
  :Fortran name:
    ``Nsavevect = [cNsave]``
