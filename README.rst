*****
Lysis
*****

Clot lysis simulation based on the work of Dr. Brittany Bannish at the
University of Central Oklahoma.

.. image:: https://github.com/UCO-OpResearch/lysis/actions/workflows/tests.yml/badge.svg
   :target: https://github.com/UCO-OpResearch/lysis/actions/workflows/tests.yml
   :alt: Tests

.. image:: https://readthedocs.org/projects/lysis/badge/?version=latest
   :target: https://lysis.readthedocs.io/en/latest/
   :alt: Documentation Status

.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0
   :alt: License: GPL v3

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.20406003.svg
   :target: https://doi.org/10.5281/zenodo.20406003
   :alt: DOI

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


Quick start
===========

A minimal end-to-end workflow, driven entirely through the ``lysis`` CLI::

    # Build an Experiment (one HDF5 file per Run) from a parameter CSV
    uv run lysis init-experiment params.csv ./data

    # Initialise the macroscale edge-grid structure in each Run's HDF5 file
    uv run lysis init-macroscale ./data/<experiment>

    # Execute the Fortran microscale, then macroscale, simulations
    uv run lysis run-micro ./data/<experiment>
    uv run lysis run-macro ./data/<experiment>

    # Inspect and compare results
    uv run lysis parameters ./data/<experiment>
    uv run lysis compare ./data/<experiment-a> ./data/<experiment-b>

See the `command-line interface`_ section below and the usage guide under
``docs/source/usage/`` for the full set of commands and options.


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

``./src/c``
    The KISS random number generator (``kiss.c`` / ``kiss.h``) --- the
    canonical source compiled into ``bin/kiss.o`` and ``lib/kiss.so`` and
    used by the Fortran and Python implementations. (The abandoned C/MPI
    macroscale port that used to live here now resides in ``archive/c/``.)

``./scripts``
    Standalone scripts for executing Fortran models, generating test
    fixtures, and other data-processing tasks.

``./tests``
    The pytest test suite.

``./notebooks``
    Jupyter notebooks for analysis and experiment work.

``./archive``
    Unmaintained, preserved-for-reference code that is no longer part of the
    build: the C/OpenMPI and C++ macroscale ports, the CuPy/GPU
    proof-of-concept, superseded Fortran macroscale variants, and the
    original MATLAB/Octave pre- and post-processing pipeline. See
    ``archive/README.rst``.

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


Releases and publications
=========================

Tagged releases are listed on the `GitHub Releases page
<https://github.com/UCO-OpResearch/lysis/releases>`__, where each release's
notes record the exact code state, the paper(s) generated with it, and the
associated Zenodo dataset where one exists.

* **v0.1.0** — initial Python wrapper around the Fortran macroscale model;
  state used for Risman et al. (2024) and Bannish et al. (2024) below.
* **v0.2.0** — protofibril packing-density paper; state used for Risman et al.
  (2025) below.
* **v0.3.0** — Python wrapper for the Fortran microscale (development release,
  no associated paper).


References
==========

The model is described in:

Bannish, Brittany E., James P. Keener, and Aaron L. Fogelson. "Modelling
fibrinolysis: a 3D stochastic multiscale model." *Mathematical medicine and
biology: a journal of the IMA* 31.1 (2014): 17-44.
https://doi.org/10.1093/imammb/dqs029

Papers that used this code (see the matching releases above):

Risman, R. A., Paynter, B., Percoco, V., Shroff, M., Bannish, B. E., &
Tutwiler, V. (2024). Internal fibrinolysis of fibrin clots is driven by pore
expansion. *Scientific Reports*, 14(1), 2623.
https://doi.org/10.1038/s41598-024-52844-4
(dataset: https://doi.org/10.5281/zenodo.8115180)

Bannish, B. E., Paynter, B., Risman, R. A., Shroff, M., & Tutwiler, V. (2024).
The effect of plasmin-mediated degradation on fibrinolysis and tissue
plasminogen activator diffusion. *Biophysical Journal*, 123(5), 610-621.
https://doi.org/10.1016/j.bpj.2024.02.002

Risman, R. A., Percoco, V., Paynter, B., Bannish, B. E., & Tutwiler, V. (2025).
Protofibril packing density of individual fibers alters fibrinolysis.
*Research and Practice in Thrombosis and Haemostasis*, 9(2), 102708.
https://doi.org/10.1016/j.rpth.2025.102708
(dataset: https://doi.org/10.5281/zenodo.15151618)


License
=======

Lysis is distributed under the GNU General Public License v3.0. See the
`LICENSE <LICENSE>`__ file for the full text.


How to cite
===========

If you use this software in your research, please cite the modelling paper
(Bannish, Keener, and Fogelson, 2014, above) together with the software itself
via its archived Zenodo record:

    https://doi.org/10.5281/zenodo.20406003

This DOI always resolves to the latest archived version; each tagged release
also has its own version-specific DOI on Zenodo. Releases and their associated
papers and datasets are listed on the `GitHub Releases page
<https://github.com/UCO-OpResearch/lysis/releases>`__.
