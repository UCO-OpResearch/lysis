============================
Experiment Initialisation
============================

Overview
--------

An **Experiment** is a collection of :term:`Runs <Run>` that are compared to
elucidate cause-and-effect relationships.  The ``lysis init-experiment``
command (and its Python API counterpart :meth:`~lysis.config.experiment.Experiment.from_csv`)
initialises a new Experiment from a CSV file:

* Each **row** in the CSV describes one Run's parameters.
* The system resolves any missing *dependent* parameters algebraically
  (so you can supply, for example, an off-rate instead of a dissociation
  constant).
* It validates that all provided values are mutually consistent.
* It creates the experiment folder, writes ``experiment.json``, and
  produces an HDF5 file for each Run containing the microscale parameters
  and empty dataset stubs.

.. code-block:: text

    {data_root}/
        {experiment_name}/
            experiment.json          ← name, description, run list
            {run_code_0}.h5          ← micro parameters + empty micro stubs
            {run_code_1}.h5
            ...

.. note::

   The HDF5 files created here contain only the **microscale** structure.
   Macroscale parameters and dataset stubs are added later by
   :meth:`~lysis.dataio.datastore.DataStore.initialize_macroscale`, which
   must be called **after** all microscale Simulations have completed.
   This is required because the ``forced_unbind`` parameter can only be
   computed from microscale output.


CSV Format
----------

Column headers
~~~~~~~~~~~~~~

Each column header must be a Python parameter name from
:class:`~lysis.config.parameters.MicroParameters` or
:class:`~lysis.config.parameters.MacroParameters`, or one of the special
metadata columns listed below.  Use **Python** names (``fiber_radius``,
``pore_size``), not Fortran names.  Unrecognised headers raise a
:exc:`ValueError`.

Cell values
~~~~~~~~~~~

Values may include a `Pint <https://pint.readthedocs.io/>`_-compatible
unit string::

    72.7 nm
    1.0135e-4 centimeter
    0.1 / micromolar / second

If a cell contains a bare number with no units, the canonical unit for
that parameter is assumed (see :meth:`~lysis.config.parameters.Parameters.units`).

Blank cells and omitted columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A **blank cell** means "not provided for this row".  An **omitted column**
is equivalent to every cell in that column being blank.  The parameter
resolver will attempt to derive the missing value algebraically from
whatever is provided.  If it cannot be derived, the parameter's default
value is used.

Special metadata columns
~~~~~~~~~~~~~~~~~~~~~~~~~

These columns are not passed to the parameter resolver:

``run_code``
    Override the auto-generated run code for that row.  If absent, a
    timestamp-based code is generated (e.g. ``2026-04-07-1422-00``).

``run_description``
    Prose note stored in ``experiment.json`` alongside the run code.


Parameter Resolution
--------------------

The :mod:`~lysis.config.param_resolver` module encodes the algebraic
relationships between parameters as `SymPy <https://www.sympy.org/>`_
equations.  This means you can supply **any consistent subset** of the
related parameters and the system will solve for the rest.

Kinetic rate constants example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The three tPA kinetic parameters are related by:

.. math::

    k_{\text{off}} = k_{\text{on}} \times K_d

You may provide any **two** of the three; the third is derived
automatically.  The following CSV fragments are all equivalent:

.. code-block:: text

    # Provide bind_rate_tPA + diss_const_tPA_woPLG → unbind_rate_tPA_woPLG solved
    bind_rate_tPA,diss_const_tPA_woPLG
    0.1 / micromolar / second,0.36 micromolar

    # Provide bind_rate_tPA + unbind_rate_tPA_woPLG → diss_const_tPA_woPLG solved
    bind_rate_tPA,unbind_rate_tPA_woPLG
    0.1 / micromolar / second,0.036 / second

    # Provide diss_const_tPA_woPLG + unbind_rate_tPA_woPLG → bind_rate_tPA solved
    diss_const_tPA_woPLG,unbind_rate_tPA_woPLG
    0.36 micromolar,0.036 / second

Grid geometry example
~~~~~~~~~~~~~~~~~~~~~

``full_row`` is defined as ``3 * cols - 1``.  You may provide either
``cols`` or ``full_row``; the other is solved automatically.

Conflict detection
~~~~~~~~~~~~~~~~~~

If you provide **all three** related values and they are inconsistent,
a :exc:`~lysis.config.param_resolver.ParameterConflict` is raised
(with a message identifying the conflicting parameter and both the
provided and expected values).

All-rows validation
~~~~~~~~~~~~~~~~~~~

Errors across **all** CSV rows are collected before raising, so
researchers see every problem at once rather than fixing them one by one.


CLI Usage
---------

.. code-block:: bash

    lysis init-experiment CSV_PATH DATA_ROOT [OPTIONS]

Arguments
~~~~~~~~~

``CSV_PATH``
    Path to the parameter CSV file.

``DATA_ROOT``
    Parent directory under which the experiment folder is created.
    Created automatically if it does not exist.

Options
~~~~~~~

.. option:: --name NAME

    Experiment name (used as the folder name under ``DATA_ROOT``).
    Defaults to the CSV filename stem.

.. option:: --description TEXT

    Optional prose description stored in ``experiment.json``.

.. option:: --dry-run

    Validate and resolve all parameters but do **not** create any files
    or folders.  Prints a summary of what would be created.

.. option:: --no-progress

    Suppress progress indicators.

Examples
~~~~~~~~

.. code-block:: bash

    # Basic usage — experiment folder named after the CSV
    lysis init-experiment runs.csv /data/experiments/

    # Custom name and description
    lysis init-experiment runs.csv /data/experiments/ \
        --name fiber-radius-sweep \
        --description "Effect of fiber radius on lysis time"

    # Validate without writing files
    lysis init-experiment runs.csv /data/experiments/ --dry-run


Python API Usage
----------------

.. code-block:: python

    from lysis.config.experiment import Experiment

    exp = Experiment.from_csv(
        csv_path="my_runs.csv",
        data_root="/data/experiments/",
        name="fiber-radius-sweep",
        description="Effect of fiber radius on lysis time",
    )
    for run in exp.runs:
        print(run.run_code, run.micro_params.fiber_radius)

:meth:`~lysis.config.experiment.Experiment.from_csv` returns the
:class:`~lysis.config.experiment.Experiment` object with all
:class:`~lysis.config.run.Run` objects populated (and, unless
``dry_run=True``, HDF5 files with microscale structure written to disk).

After all microscale Simulations complete, call
:meth:`~lysis.dataio.datastore.DataStore.initialize_macroscale` for each
Run to compute ``forced_unbind`` and write the macroscale structure:

.. code-block:: python

    from lysis.dataio.datastore import DataStore

    # After microscale Simulations are done:
    for run in exp.runs:
        with DataStore(run.run_code, exp.path, mode="a") as ds:
            ds.initialize_macroscale(run.macro_params)


Error Messages
--------------

:exc:`~lysis.config.param_resolver.ParameterConflict`
    One or more CSV rows contain mutually inconsistent parameter values.
    The message lists every conflicting parameter, the value you provided,
    and the value that was computed from the other parameters.

    **Fix**: remove or correct the over-constrained value.  Usually you
    only need to supply **two** of the three related kinetic parameters
    (the third will be derived).

:exc:`~lysis.config.param_resolver.UnderdeterminedParameters`
    The system cannot derive a value for one or more parameters from what
    was provided, and no default exists.

    **Fix**: add the missing parameter (or one from which it can be
    derived) to the CSV.

:exc:`FileExistsError`
    The experiment folder already exists.

    **Fix**: use ``--name`` to choose a different folder name, or delete
    the existing folder first.

:exc:`ValueError` (unrecognised column)
    A CSV header does not match any known parameter name or metadata column.

    **Fix**: check the spelling.  Parameter names are Python attribute
    names (e.g. ``fiber_radius``, ``pore_size``), not Fortran names.


Template CSV
------------

A template CSV covering every recognised parameter is provided at
``docs/usage/experiment_template.csv``.  It contains three example rows:

1. **all-defaults** — every independent parameter set to its default value.
2. **minimal** — only a handful of parameters provided; all others are
   derived from defaults.
3. **dep-koff** — kinetic off-rates (``unbind_rate_tPA_woPLG``,
   ``unbind_rate_PLG_intact``) specified instead of the corresponding
   dissociation constants; the resolver back-calculates the missing values.

.. literalinclude:: ../../../docs/usage/experiment_template.csv
   :language: text
   :caption: docs/usage/experiment_template.csv
