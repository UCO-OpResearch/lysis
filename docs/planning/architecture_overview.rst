=====================================
Lysis Project Architecture Overview
=====================================

:Audience: New team members (software development background assumed)
:Last Updated: February 2026

.. contents:: Table of Contents
   :depth: 2
   :local:


What This Project Does
======================

The lysis project simulates **fibrinolysis** --- the enzymatic breakdown of
blood clots.  Blood clots are meshes of fibrin fibres; the body dissolves
them by deploying the enzyme **tPA** (tissue-type plasminogen activator),
which converts inactive **plasminogen (PLG)** into **plasmin**, the enzyme
that actually cuts fibrin.  Too-slow lysis causes heart attacks and strokes;
too-fast lysis causes excessive bleeding.  The simulation helps researchers
understand how molecular-level chemistry produces the clot-scale "lysis
front" observed in experiments.

The model has two scales:

Microscale
    Simulates a single fibre cross-section at molecular resolution.  A
    square grid of binding locations tracks individual tPA, PLG, and plasmin
    molecules as they bind, activate, degrade fibrin, and crawl to
    neighbouring sites.  Uses the Gillespie algorithm for stochastic
    event-time selection.  Hundreds of independent microscale Simulations
    produce statistical distributions of tPA leaving times, plasmin counts,
    and fibre lysis times.

Macroscale
    Simulates the full 3-D fibrin clot as a rectilinear lattice (the
    "edge grid"), where each lattice edge represents one fibre.  tPA
    molecules diffuse across the lattice, bind to fibres, and cause
    degradation according to the microscale distributions.  The macroscale
    model consumes the microscale output; it cannot execute without it.

In short, the microscale captures *chemistry on a single fibre*, and the
macroscale captures *transport and geometry across the whole clot*.  Data
flows **one way**: microscale output feeds macroscale input.

For the formal definitions of *Run*, *Experiment*, *Scenario*, *Mechanism*,
and other project-specific terms, see the
`Ontology <../source/usage/ontology.rst>`_.


Why This Project Exists
=======================

The original simulation was written entirely in Fortran (microscale) and
C/C++ (macroscale).  It works, and it has produced published results.  But
it has several practical problems:

1. **Accessibility.**  Few new students read Fortran or C fluently.
   Modifying the code requires deep familiarity with both the science and the
   low-level implementation.

2. **Mechanism proliferation.**  To try a new hypothesis (a new *Mechanism*
   in our terminology), the current workflow is to **copy the entire
   simulation codebase** and edit the copy.  Over time this has produced
   multiple divergent forks that are difficult to compare or maintain.

3. **Data fragmentation.**  Different code versions produce data in
   different formats (text files, binary dumps, JSON, Fortran log files).
   Comparing results across Runs and Experiments requires manual data
   wrangling.

4. **Reproducibility.**  Published results must remain reproducible even as
   the code evolves.  Legacy data from earlier code versions must be
   importable into the current system.

The lysis Python package addresses all four problems.


Architectural Goals
===================

The project has four top-level goals, roughly in priority order.

Goal 1: One Codebase, Many Mechanisms
--------------------------------------

Users should be able to execute Simulations with different Scenarios (numeric
parameters) and different Mechanisms (logic/rule variations) **from a single
codebase**, without copying files.

The envisioned approach is a **strategy pattern** inspired by functional
programming.  The core simulation loop remains fixed; Mechanism-specific
behaviour is injected via **swappable functions** selected by configuration
flags.

A concrete example already exists in the macroscale Python code
(``np_macroscale.py``), where the ``duplicate_fortran`` flag controls 13
branch points --- RNG selection, neighbourhood calculation, molecule
placement, random-number generation order, and several others.  Today these
are inline ``if``/``else`` blocks.  The goal is to extract each branch into a
named function and let users compose a Mechanism by selecting which functions
to use:

.. code-block:: text

   Mechanism = {
       rng:                  "kiss"  |  "numpy_mt",
       neighborhood:         "fortran_order"  |  "python_order",
       molecule_placement:   "fortran_flat"  |  "python_grid",
       move_logic:           "fortran_reversed"  |  "python_standard",
       ...
   }

This would allow new hypotheses (e.g., a new unbinding rule) to be added as
a single function, registered in a table, and selected at execution time ---
without touching the rest of the simulation code.

This goal is **aspirational**: no pluggable-mechanism infrastructure exists
yet.  Building it is a key area where a software engineering perspective can
contribute.

Goal 2: Unified Data Format
-----------------------------

All simulation data --- past, present, and future --- should be stored in a
single, efficient format: **HDF5**.  A well-defined data specification
(versioned, validated) ensures that any analysis tool can read any dataset
regardless of which code version or Mechanism produced it.

This goal is **largely achieved**.  The ``lysis.dataio`` package provides:

- A versioned data specification system (``dataspec``) with three supported
  versions (v1.95.0, v1.99.0, v2.0.0).
- Automatic multi-step conversion routing between any pair of spec versions.
- Read/write access to HDF5 files through the ``DataStore`` API.
- Import support for all legacy formats (Fortran text, binary, JSON, and
  parsed log files).

See the `Data Handling Completion Plan <data_handling_completion_plan.rst>`_
for full details.

Goal 3: Standard Metrics and Visualizations
--------------------------------------------

Once data is in a unified format, standard analysis and visualizations can
be applied to any Run or Experiment for direct comparison.  This includes
degradation-front animations, molecule tracking plots, and summary
statistics.

Several Jupyter notebooks in ``notebooks/`` already perform analysis and
visualization (e.g., ``Degradation Animation.ipynb``,
``Degradation Front.ipynb``, ``F-Macro Compare.ipynb``).  However, these are
ad-hoc and have not been consolidated into a reusable analysis layer.
Standardizing and documenting these tools is future work.

Goal 4: Fortran-to-Python Migration
-------------------------------------

The primary development task is to **port simulation logic from Fortran/C to
Python** to make the code more accessible.  Current status:

:Microscale: **Fortran only.**  No Python implementation exists yet.
:Macroscale: **Python written, not yet validated.**  The Python macroscale
    (``np_macroscale.py``) has been implemented with a ``duplicate_fortran``
    mode that replicates the Fortran logic step-by-step for validation
    purposes.  Until validation is complete, the Fortran macroscale
    continues to be used in production.

Validation means confirming that the Python code produces **identical
results** to the Fortran code for the same inputs and RNG seeds.  This is
necessary because published results must remain reproducible, and because any
discrepancy could indicate a bug in either implementation.

During the migration, the Fortran code must continue to work.  The two
implementations coexist in the same repository:

- ``src/fortran/`` --- Fortran microscale source
- ``src/c/``, ``src/cpp/`` --- C/C++ macroscale source
- ``src/lysis/`` --- Python package (macroscale, data handling, tools)


Key Design Decisions
====================

Development Speed Over Runtime Performance
--------------------------------------------

This is scientific research code.  Any given version of the simulation is
only executed a few dozen times before the code is modified again.
**Development speed matters much more than runtime performance.**  Python was
chosen over Fortran for new code specifically because it is faster to write,
read, and modify.

That said, performance is not ignored: NumPy vectorization, CuPy GPU support
(optional), and efficient HDF5 I/O are used where they matter.

HPC with Local Testing
-----------------------

Production Simulations execute on **HPC clusters via SLURM**.  The codebase
includes SLURM utility functions (``lysis.tools.slurm``) and execution
scripts (``scripts/``).  However, the code must also execute on local
workstations for development and testing.

Validation Against Legacy Data
-------------------------------

Every change must be validated against existing results.  The data system
supports three specification versions precisely so that legacy Fortran data
can be imported, converted, and compared with new Python output.  The
``duplicate_fortran`` mode in the macroscale is another expression of this
principle: it exists solely to produce bit-identical output for validation.


Repository Layout
=================

.. code-block:: text

   lysis/
   +-- src/
   |   +-- lysis/              Python package (installed via pip install -e .)
   |   |   +-- config/         Scenario parameters, Mechanism parameters, Run
   |   |   +-- dataio/         Data specs, I/O, conversion, DataStore
   |   |   +-- geometry/       Edge grid, coordinate transforms
   |   |   +-- execution/      Fortran subprocess wrappers
   |   |   +-- tools/          Utilities (RNG, SLURM, misc)
   |   |   +-- cli/            Command-line interface (convert, validate)
   |   |   +-- np_macroscale.py   Python macroscale simulation (NumPy)
   |   |   +-- molecule.py        Molecule tracking
   |   +-- fortran/            Fortran microscale source
   |   +-- c/                  C macroscale source (original)
   |   +-- cpp/                C++ macroscale source (variant)
   +-- tests/                  pytest test suite (630+ tests)
   +-- notebooks/              Jupyter notebooks (analysis, visualization)
   +-- scripts/                Execution and utility scripts
   +-- archive/                Unmaintained experiments (stale GPU/CuPy port)
   +-- docs/
   |   +-- source/usage/       Sphinx documentation (RST)
   |   +-- planning/           Planning documents (not published)
   |   +-- papers/             Reference publications
   +-- data/                   Simulation data (git-ignored)


Current Progress
================

Data Handling System (~97% Complete)
-------------------------------------

The data handling system is the most mature part of the Python codebase.
It provides:

- **DataStore API** --- read/write access to HDF5 files with dot-notation,
  per-Simulation views, parameter loading, and spec version validation.
- **Three data specification versions** (v1.95.0, v1.99.0, v2.0.0) with
  automatic multi-step conversion routing.
- **Five storage format backends** (HDF5, text, binary, JSON, Fortran
  parsed) for importing legacy data.
- **Parameter validation** --- strict loading with missing-parameter
  detection, Fortran log parsing, and cross-verification.
- **CLI tools** --- ``lysis convert`` and ``lysis validate`` for command-line
  data operations.
- **630+ automated tests** covering all modules, including 68 integration
  tests with real Fortran fixture data.

For the complete status and remaining items, see the
`Data Handling Completion Plan <data_handling_completion_plan.rst>`_.

Macroscale Simulation (Written, Pending Validation)
----------------------------------------------------

The Python macroscale simulation (``np_macroscale.py``) is functionally
complete and integrated with the DataStore for both reads and writes.  It
includes a ``duplicate_fortran`` mode for step-by-step validation against
the Fortran implementation.  Validation is in progress.

Microscale Simulation (Fortran Only)
-------------------------------------

The microscale model exists only in Fortran.  Porting it to Python is a
long-term goal but is not currently in progress.

Analysis and Visualization (Ad-Hoc)
-------------------------------------

Jupyter notebooks exist for various analyses (degradation animations,
molecule tracking, macroscale comparison).  These are not yet organized into
a reusable analysis framework.


What Comes Next
===============

Phase 4: Simulation Integration (Current)
-------------------------------------------

- Integrate Fortran subprocess wrappers (``codeutil.py``) with the DataStore
  so that Fortran Simulations read/write HDF5.
- Validate the Python macroscale against Fortran output.
- Add dimensioned quantities (``pint.Quantity``) to the macroscale
  simulation for unit consistency.
- Add data specifications for remaining Fortran simulation types (``f_deg``,
  combined macro).

Phase 5: Quality and Infrastructure
-------------------------------------

- Set up a CI pipeline for automated testing.
- Complete documentation and docstrings.
- Update Jupyter notebooks to current package structure.
- Optimize HDF5 chunking for large datasets.

Future: Pluggable Mechanisms
-----------------------------

- Design and implement the strategy-pattern infrastructure for swappable
  Mechanism functions.
- Extract the ``duplicate_fortran`` branches into registered function
  alternatives.
- Build a Mechanism configuration system that lets users compose Mechanisms
  from function selections.
- Port the microscale model to Python.
- Standardize the analysis/visualization layer.


Getting Started
===============

**Fortran tools:**  Setup and execution guides for the Fortran microscale and
macroscale simulations are in ``docs/source/usage/`` (see
``fortran_microscale.rst`` and ``fortran_macroscale.rst``).

**Python package:**  Getting-started documentation for the Python package and
the HDF5 data system does not yet exist and needs to be created.  In the
meantime:

- Install the package in development mode: ``pip install -e .`` from the
  repository root.
- Execute the test suite: ``pytest`` from the repository root.
- Explore the ``notebooks/`` directory for working examples of data
  handling and analysis.

**Vocabulary:**  Read the `Ontology <../source/usage/ontology.rst>`_ early.
The project uses specific terms (*Run*, *Experiment*, *Scenario*,
*Mechanism*, *Simulation*) with precise meanings that differ from casual
usage.  In particular, avoid saying "run code" --- say "execute code" instead,
because *Run* has a specific technical meaning in this project.


Key Documents
=============

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Document
     - Purpose
   * - `Ontology <../source/usage/ontology.rst>`_
     - Definitions of all project-specific terms
   * - `Data Handling Completion Plan <data_handling_completion_plan.rst>`_
     - Detailed status and architecture of the data system
   * - `Data Specification <../source/usage/data_specification.rst>`_
     - Technical specification of the HDF5 data format
   * - `Git Workflow Guide <git_workflow_guide.rst>`_
     - Team git workflow and branching conventions
   * - Bannish et al. (2014)
     - Original paper describing the fibrinolysis model (in ``docs/papers/``)
