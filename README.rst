*****
Lysis
*****

Clot lysis simulation based on the work of Dr. Brittany Bannish at the
University of Central Oklahoma.

Lysis is a computational model for studying the breakdown (lysis) of blood
clots. It implements the multiscale stochastic model of fibrinolysis described
by Bannish, Keener, and Fogelson (see References_) and provides Python tooling
to initialise, execute, convert, validate, and analyse simulation data.

The model has two coupled scales:

* **Microscale** — a single-fiber stochastic model of tPA binding, unbinding,
  and fiber degradation.
* **Macroscale** — diffusion and lysis across a rectilinear edge grid that
  represents the clot, driven by the microscale results.

The original microscale and macroscale models are written in Fortran. The
``lysis`` Python package wraps and extends them with data I/O, format
conversion, and analysis tools, and is gradually reimplementing the models in
Python (with optional GPU acceleration via CuPy).


Documentation
=============

Full documentation — including the project ontology, the on-disk data
specification, and guides for initialising and executing simulations — is
published on Read the Docs at https://lysis.readthedocs.io and lives in
``docs/source/``.

To build the docs locally::

    uv sync --extra notebooks
    uv run sphinx-build -b html docs/source docs/build/html

Source repository: https://github.com/UCO-OpResearch/lysis


Installation
============

The project uses `uv <https://docs.astral.sh/uv/>`_ for dependency management
and Python 3.11. To create the environment and install the package (editable)
into ``.venv/``::

    uv sync

Optional dependency groups:

* ``--extra test`` — pytest and the test suite dependencies
* ``--extra notebooks`` — Jupyter, NetworkX, and notebook formatting tools
* ``--extra gpu`` — CuPy and NVTX for GPU-accelerated macroscale code

For example, to set up for running the tests::

    uv sync --extra test

Run any command in the environment with ``uv run`` (e.g. ``uv run lysis``,
``uv run pytest``).


Command-line interface
======================

Installing the package provides the ``lysis`` command. Run ``uv run lysis
--help`` or ``uv run lysis <command> --help`` for full details.

``init-experiment``
    Initialise an Experiment from a parameter CSV file.
``init-macroscale``
    Initialise macroscale structure in HDF5 files.
``run-micro``
    Execute the Fortran microscale simulation for a Run.
``run-macro``
    Execute the Fortran macroscale simulation for a Run.
``convert``
    Convert simulation data between specification formats.
``validate``
    Validate simulation data against a specification.
``parameters``
    Print parameter tables for one or more Runs.
``micro-stats``
    Print microscale statistics for one or more simulations.
``macro-stats``
    Print macroscale summary statistics.
``deg-rate``
    Print degradation-rate tables.
``deg-time``
    Print degradation-time tables.
``compare``
    Compare Runs across two folders or two ``.h5`` files.
``diff``
    Diff HDF5 data tables between two Runs or folders.
``rename``
    Rename an Experiment folder or a Run HDF5 file.


File Organization
=================

``./src/lysis``
    The main Python package for the simulation. Submodules include
    ``config`` (Scenario constants, parameters, Run container), ``dataio``
    (data specifications, HDF5 I/O, format conversion), ``geometry`` (edge
    grid and coordinate transforms), ``execution`` (Fortran subprocess
    wrappers), ``analysis`` (degradation analysis and plotting), ``cli`` (the
    ``lysis`` command), and ``tools`` (utilities, RNG, Slurm helpers).

``./src/fortran``
    The original microscale and macroscale Fortran models.

``./src/matlab``
    MATLAB code used to pre- and post-process data for the Fortran models,
    along with Jupyter Notebook versions that run under the Octave kernel.

``./src/cpp``
    An incomplete conversion of the macroscale model to C++.

``./src/c``
    An incomplete conversion of the macroscale model to C, using MPI to
    multithread the macroscale grid. Also includes the KISS random number
    generator used by the other implementations.

``./scripts``
    Standalone scripts for executing Fortran models, generating test
    fixtures, and other data-processing tasks.

``./tests``
    The pytest test suite.

``./notebooks``
    Jupyter notebooks for analysis and experiment work.

``./docs``
    Sphinx documentation source (``docs/source``) and planning notes.

``./bin``
    Compiled Fortran executables. *(Generated locally; not tracked in git.)*

``./lib``
    Compiled shared libraries (e.g. ``kiss.so``, the KISS random number
    generator for use by Python code). *(Generated locally; not tracked in
    git.)*

``./data``
    The stored output of experimental runs. *(Generated locally; not tracked
    in git.)* Each Experiment is a folder containing an ``experiment.json``
    file — recording the Experiment's name, description, and list of Runs —
    plus one HDF5 (``.h5``) file per Run. Each Run's ``.h5`` file holds that
    Run's parameters together with all of its microscale and macroscale
    simulation data.

    These HDF5 files can be inspected with the ``lysis`` command-line tools or
    with general-purpose HDF5 viewers such as
    `HDFView <https://www.hdfgroup.org/download-hdfview/>`_.


.. _References:

References
==========

Bannish, Brittany E., James P. Keener, and Aaron L. Fogelson. "Modelling
fibrinolysis: a 3D stochastic multiscale model." *Mathematical medicine and
biology: a journal of the IMA* 31.1 (2014): 17-44.
https://doi.org/10.1093/imammb/dqs029
