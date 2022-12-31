***********
Lysis
***********

Clot Lysis Simulation by Dr. Brittany Bannish of the University of Central Oklahoma

File Organization
===========

./src
^^^^^^^^^^^
Source code

./src/fortran
"""""""""""
The original microscale and macroscale models

./src/matlab
"""""""""""
Code used to pre- and post-process data for the Fortran model code

Also includes Jupyter Notebook versions of Matlab code. These will run using the Octave kernel for Jupyter.

./src/python
"""""""""""
An in-progress conversion of the model to Python

* Current features
    * Macroscale
        * Import/export of experiment parameters and data

* Future features
    * Microscale
        * Import/export of experiment parameters and data
        * Full conversion from Fortran
    * Macroscale
        * OpenMPI threading of the macroscale grid
        * CUDA threading of the macroscale grid


./src/cpp
"""""""""""
An incomplete conversion of the macroscale model to C++

./src/c
"""""""""""
An incomplete conversion of the macroscale to C.

This code uses MPI to multithread the macroscale grid

Also includes KISS random number generator used by all other implementations.

./data
^^^^^^^^^^^
The stored output of experimental runs.

Each experiment is in its own folder named by its experiment code (usually the date and time in YYYY-MM-DD-hhmm format).

Each experiment folder should contain a params.json file containing the experiment parameters.

./lib
^^^^^^^^^^^
Compiled shared libraries

* kiss.so - The KISS random number generator for use by Python code

References
===========

Bannish, Brittany E., James P. Keener, and Aaron L. Fogelson. "Modelling fibrinolysis: a 3D stochastic multiscale model."
*Mathematical medicine and biology: a journal of the IMA* 31.1 (2014): 17-44. https://doi.org/10.1093/imammb/dqs029