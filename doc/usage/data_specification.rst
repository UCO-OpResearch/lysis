=========================
Data Specification v2.0.0
=========================

*Note: This document is written with NumPy indexing, which is zero-indexed.
This is especially noteworthy since the Fortran code is one-indexed.*

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

  :Data Type: 
    Plain text

``micro.txt``
  Output log for microscale simulation.

  :Data Type: 
    Plain text

``pli_first_time``
  Time to first plasmin

  :Data Type: 
    NumPy float64 (``f8``)
  :Units:
    seconds
  :Dimensions: 
    (``micro_params.simulations``,)
  :Fortran name:
    ``firstPLi``


``tpa_final_num``
  The number of tPA molecules in the fiber at the end of the simulation.

  :Data Type: 
    NumPy unsigned int 1 byte (``u1``)
  :Units:
    None
  :Dimensions: 
    (``micro_params.simulations``,)
  :Fortran name:
    ``lasttPA``
``fiber_degraded``
  Whether or not lysis is complete at the end of the simulation. That
  is, at least :math:`\frac{2}{3}` of binding doublets were degraded.
  1 for yes, 0 for no.

  :Data Type: 
    NumPy boolean (``?``)
  :Units:
    None
  :Dimensions: 
    (``micro_params.simulations``,)
  :name:
    ``lysis_complete``

``sim_final_time``
  The amount of time elapsed in each simulation. If lysis completed,
  this is the time at which that occurred.

  :Data Type: 
    NumPy float64 (``f8``)
  :Units:
    seconds
  :Dimensions: 
    (``micro_params.simulations``,)
  :Fortran name:
    ``lysis.dat``

``pli_generated_num``
  The number of Plasmin molecules generated in the fiber
  by the end of the simulation.

  :Data Type: 
    NumPy integer (``u2``)
  :Units:
    None
  :Dimensions: 
    (``micro_params.simulations``,)
    :Fortran name:
    ``PLi``

``tpa_leaving_time``
  The simulation time elapsed when the tPA molecule leaves the system
  or infinity if the simulation ends
  with the tPA molecule still bound. 
  
  :Data Type: 
    NumPy float64 (``f8``)
  :Units:
    seconds
  :Dimensions: 
    (``micro_params.simulations``,)
  :Fortran name:
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
  :fortran name:
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
  :Fortran name:
    ``tPAunbind``


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

``fiber_degrade_time``
  A list of updates to the degrade time of any fiber in the model.
  Each row represents one update to one fiber and contains 3 entries:

  #. The simulation time elapsed when the event occurred.
  #. The location index of the fiber on which the event occurred.
  #. The new degrade time for the fiber. (the time at which the fiber will degrade if no further updates occur).

  :Data Type: 
    (NumPy float (``f8``), NumPy int32 (``u4``), NumPy float (``f8``))
  :Units:
    (seconds, None, seconds)
  :Dimensions: 
    (number of binding events, 3)
  :Fortran name:
    ``f_deg_list``

``tpa_bind_events``
  A list of tPA molecule bindings and their associated metrics.
  Each row represents one update to one tPA molecule and contains 4 entries:

  #. The simulation time elapsed when the event occurred.
  #. The index of the tPA molecule which was involved in the event.
  #. The new status of the tPA molecule (see below for details).
  #. The location (fiber index) at which the event occurred. 

  Valid molecule statuses are:

  0. Kinetically Unbound.
  1. Bound to an undegraded fiber.
  2. Bound to a degraded fiber (macro-unbound).
  3. Bound to a fiber degradation product (micro-unbound).

  :Data Type: 
    (NumPy float64 (``f8``), NumPy integer (``u8``), NumPy integer (``u1``), NumPy integer (``u4``))
  :Units:
    (seconds, None, None, None)
  :Dimensions: 
    (number of binding events, 4) 
  :Fortran name:
    ``m_bind_t``

``tpa_location_snapshot``
  An array, giving the location (fiber index) of each molecule at point when a save is made.
  That is, ``tpa_location_snapshot[i, j]`` is the location of tPA molecule ``i`` when 
  snapshot ``j`` is recorded.

  :Data Type: 
    NumPy integer (``i4``)
  :Units:
    None 
  :Dimensions: 
    (number of molecules, number of snapshots)
  :Fortran name:
    ``m_loc``

``tpa_transit_time``
  The simulation time elapsed when each tPA molecule reached the back row of the 
  fiber grid for the first time.

  :Data Type: 
    NumPy float (``f8``)
  :Units:
    seconds
  :Dimensions: 
    (number of molecules,)
  :Fortran name:
    ``mfpt``

``snapshot_time``
  The simulation time elapsed at the point when each snapshot is recorded.

  :Data Type: 
    NumPy float64 (``f8``)
  :Units:
    seconds
  :Dimensions: 
    (number of snapshots,)
  :Fortran name:
    ``tsave``


Macro to Micro files
++++++++++++++++++++

:NOTE: These files are needed to input the data from the microscale model into the Fortran macroscale code.
Historically these were generated with the Matlab ``micro_to_macro`` code.
These files should NOT be stored in this version of the data specification, 
but should be generated when needed from the microscale data.

*We continue to use 0-indexing here to be consistent with the rest of this document, 
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
    (``micro_params.simulations`` // 100, 100)

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
  but it is now generated by the Python wrapper

  *For more detail on how the edge grid co-ordinates are defined in both Fortran and Python,
  see the documentation of the Python lysis.util.EdgeGrid class (link below)*

  https://github.com/UCO-OpResearch/lysis/blob/c1e6b2a92758fb8620d6f5a2223976a1478c3231/src/python/lysis/util/edge_grid.py

  :File Type:
    Space-delimited Text
  :Data Type:
    integer (``i4``)
  :Units:
    None
  :Dimensions:
    (8, ``macro_params.total_edges``)
