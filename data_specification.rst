=========================
Data Specification v1.9.0
=========================

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
  Output log for the run.

  :Data Type: 
    Plain text

``firstPLi.dat``
  Time to first plasmin

  :Data Type: 
    Binary (double)
  :Units:
    seconds
  :Dimensions: 
    1 x ``micro_params.simulations``

``lasttPA.dat``
  The number of tPA molecules in the fiber at the end of the simulation.

  :Data Type: 
    Binary (int32)
  :Units:
    None
  :Dimensions: 
    1 x ``micro_params.simulations``

``lyscomplete.dat``
  Whether or not lysis is complete at the end of the simulation.
  1 for yes, 0 for no.

  :Data Type: 
    Binary (int32)
  :Units:
    None
  :Dimensions: 
    1 x ``micro_params.simulations``

``lysis.dat``
  The amount of time elapsed in each simulation. If lysis completed,
  this is the time at which that occurred.

  :Data Type: 
    Binary (double)
  :Units:
    seconds
  :Dimensions: 
    1 x ``micro_params.simulations``

``PLi.dat``
  The number of Plasmin molecules generated in the fiber
  by the end of the simulation.

  :Data Type: 
    Binary (int32)
  :Units:
    None
  :Dimensions: 
    1 x ``micro_params.simulations``

``tPA_time.dat``
  The simulation time elapsed when the tPA molecule leaves the system.

  :Data Type: 
    Binary (double)
  :Units:
    seconds
  :Dimensions: 
    1 x ``micro_params.simulations``

``tPAPLiunbd.dat``
  Whether or not tPA was forced to unbind by plasmin-mediated
  degradation of fibrin.
  1 for yes, 0 for no.

  :Data Type: 
    Binary (int32)
  :Units:
    None
  :Dimensions: 
    1 x ``micro_params.simulations``

``tPAunbind.dat``
  Whether or not tPA unbinds from the fiber by itself.
  1 for yes, 0 for no.

  :Data Type: 
    Binary (int32)
  :Units:
    None
  :Dimensions: 
    1 x ``micro_params.simulations``


Subfolders
----------
:Purpose:
  Stores all files related to a macroscale simulation

:Name: 
  A two-digit number giving the macroscale simulations's position
  in the run array
