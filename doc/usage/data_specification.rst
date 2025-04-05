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
