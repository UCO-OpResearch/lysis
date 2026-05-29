========================
Fiber Size Conventions
========================

*The size of a fiber is referred to in several different ways throughout this
project: sometimes as a radius, sometimes as a diameter; sometimes in microns,
sometimes in nanometres, and (in one place) in centimetres. This document
collects those conventions in one place so that any given reference can be
read unambiguously.*

The one rule
============

The stored parameter is **always a radius**, and its **canonical unit is
microns**. Everything else you will encounter is one of three things:

* an *input* convention (how a fiber size is specified),
* a *display* convention (how it is printed), or
* a *derived* value (something computed from the radius).

The canonical fiber is the ``Q2`` fiber:

.. math::

   72.7\ \text{nm diameter}
   \;\longrightarrow\;
   36.35\ \text{nm radius}
   \;\longrightarrow\;
   0.03635\ \mu\text{m radius}

The expression ``72.7/2/1000`` appears verbatim as the default ``radius`` in
the current Fortran source, and ``Q_("72.7/2 nanometers")`` is the default
:py:attr:`fiber_radius` in the Python :class:`~lysis.config.parameters.MicroParameters`.

.. warning::

   The fiber bundle sizes quoted in the Fortran header comments and in the
   literature (e.g. *"72.7 nm diameter fibers"*) are **diameters**. The
   ``radius`` variable holds **half** of that value. The ``FIBER_TYPES`` table
   in :mod:`lysis.config.constants` carries this same reminder:
   *"The fiber radius is half the fiber bundle diameter listed in the source
   code."*

Radius vs. diameter
===================

Used as a **radius** (the stored form):

* :py:attr:`MicroParameters.fiber_radius <lysis.config.parameters.MicroParameters>`
  and the Fortran ``radius`` variable (``micro_rates.f90``,
  ``macro_diffuse_into_and_along__{external,internal}.f90``).
* The microscale binding-site volume, which uses the cross-sectional area
  :math:`\pi r^2`: ``vol3 = 0.001*(radius)**2*pi`` in ``micro_rates.f90``.
* The Python concentration calculations ``protein_per_fiber`` and
  ``fibrin_conc_per_fiber``, which both use ``fiber_radius**2``.

Used as a **diameter** (:math:`2r`):

* The grid node spacing. In Python::

     grid_node_distance = pore_size + 2 * fiber_radius

  and the equivalent in the macroscale Fortran::

     dist = delx*10000 + radius*2

  Here ``delx*10000`` converts the pore size from centimetres to microns, and
  ``radius*2`` is the full fiber diameter. For the canonical fiber this gives
  :math:`1.0135 + 0.0727 = 1.0862\ \mu\text{m}`.

Units
=====

Three unit systems are in play. The value is the same physical quantity in
each; only the representation differs (Pint reconciles them automatically on
the Python side).

:microns:

   The **canonical** internal unit. The radius is stored here and grid spacing
   is computed here.

:nanometres:

   How fibers are **specified** and **displayed**. The ``FIBER_TYPES`` table
   lists ``Q0``--``Q4`` in nm, the CLI prints ``fiber_radius`` in nm, and tests
   provide values in nm.

:centimetres:

   Used **only** for the pore size (Fortran ``delx`` :math:`= 1.0135\times10^{-4}`
   cm). The Python :py:attr:`pore_size` default is written ``Q_("1.0135 um")``,
   and the conversion to centimetres happens when the value is passed to the
   Fortran code.

.. note::

   The pore size is **not** a fiber dimension. It is the gap *between* fibers.
   The center-to-center node spacing is the pore gap plus one full fiber
   diameter, which is why ``grid_node_distance = pore_size + 2 * fiber_radius``.

Standard fiber types
====================

The ``FIBER_TYPES`` table in :mod:`lysis.config.constants` maps each fiber code
to a **radius** (in nanometres) and a row width. The corresponding diameters
(twice the radius) are the values quoted in the literature and Fortran comments.

=========  ==============  ================  ===================
Code       Radius (nm)     Diameter (nm)     Nodes per micro row
=========  ==============  ================  ===================
Q0         23.0            46.0              4
Q1         28.7            57.4              5
Q2         36.35           72.7              7
Q3         40.65           81.3              8
TF-v       52.55           105.1             5
TF-vii     52.55           105.1             7
TF-x       52.55           105.1             10
TB-xi      61.5            123.0             11
TB-xiii    61.5            123.0             13
Q4         72.7            145.4             13
=========  ==============  ================  ===================

How to read any reference
=========================

* A ``radius`` / ``fiber_radius`` reference is a **radius in microns** (or in
  nm if it is a specification or display point).
* ``2 * fiber_radius``, ``radius*2``, or any ``dist`` / ``grid_node_distance``
  is a **diameter** added to the pore gap for node spacing.
* A ``**2`` on the radius is :math:`r^2` for an area or volume.
* ``delx`` / ``pore_size`` is never the fiber itself; it is the gap, and it is
  the one value carried in centimetres on the Fortran side.

How the Fortran currently uses radius and diameter
==================================================

Microscale (``micro_rates.f90``)
--------------------------------

The radius enters a single calculation -- the volume of one "location" (a
cleavage site), modelled as a 1 nm length of fiber::

   vol3 = 0.001*(radius)**2*pi      ! µm³ : 0.001 µm length × π r²

``vol3`` then feeds the local binding-site concentration inside the ``movetpa``
subroutine::

   bs0 = exposed/vol3/602.2          ! current binding-site concentration, µM
   k0  = ktPAon*bs0                  ! binding rate
   p_rebind = ...                    ! 3-D rebinding-escape probability

``p_rebind`` is the probability that a just-unbound tPA molecule rebinds to the
fiber before diffusing away in three dimensions. **Microscale rebinding is not
implemented**: if a random draw falls below ``p_rebind`` the code only prints a
warning (``'... must add rebinding to code'``) and continues, and with the
standard parameters ``p_rebind`` is ≈ 0. So ``radius`` currently has no material
effect on microscale output -- it parameterises a rebinding probability the
model assumes negligible. (A separate ``bs`` assignment near line 868 is also
computed from ``vol3`` but never read -- it is vestigial.)

Macroscale (``macro_diffuse_into_and_along__*.f90``)
----------------------------------------------------

The radius is combined with the pore size into the grid node spacing::

   dist = delx*10000 + radius*2      ! center-to-center node spacing, µm

but ``dist`` is **only echoed to stdout** -- the one line that would use it in a
calculation is commented out -- so it does not enter the dynamics. The timestep
is set by the **pore size** (``delx``), the move probability ``q``, and the
diffusion coefficient -- not by the radius::

   tstep = q*delx**2/(12*Diff)

The quantity that actually drives macroscale binding is the binding-site
concentration ``bs``, supplied as a separate ``--bs`` argument and used in the
binding-time draw ``bind(j) = t - log(r1)/(kon*bs) - tstep/2``. The macro's own
``--radius`` argument is therefore inert.

.. warning::

   In the Fortran code ``radius`` and ``bs`` are **independent inputs** -- each
   binary reads them as two unrelated numbers and never checks that they are
   consistent with each other. **In reality they are dependent**: the
   binding-site concentration is a function of the fiber radius
   (``bs`` :math:`\propto 1/\text{radius}^2`, via ``fibrin_conc_per_fiber`` /
   ``binding_sites``). If you invoke a Fortran binary directly with mismatched
   ``--radius`` and ``--bs`` it will run happily and produce physically
   inconsistent results. Always derive the two together -- by hand, or
   preferably by letting the Python pipeline compute ``binding_sites`` from
   ``fiber_radius`` (see below).

In the Fortran--Python pipeline
-------------------------------

The Python layer is what keeps the radius and the binding-site concentration
consistent:

#. ``fiber_radius`` is a ``MicroParameters`` field (canonical unit microns).
#. ``binding_sites`` is **derived** from it -- ``fibrin_conc_per_fiber`` (and
   hence ``binding_sites``) is computed as :math:`\propto 1/\text{fiber\_radius}^2`
   in ``lysis.config.param_resolver`` / ``parameters.py``. Changing the radius
   automatically updates the binding-site concentration.
#. **Microscale execution** (``FortranMicro``): the generic argument builder
   emits ``--radius <µm>`` whenever the run's radius differs from the class
   default. ``bs`` is *not* passed to the micro binary (it computes -- but does
   not use -- its own).
#. **Macroscale execution** (``FortranMacro``): ``_post_arguments`` passes
   **both** ``--radius`` (from ``micro_params.fiber_radius``) and ``--bs`` (from
   ``micro_params.binding_sites``) to the macro binary. The two are guaranteed
   consistent because both come from the same ``MicroParameters`` object.
#. The physical node spacing used for analysis (e.g. lysis-front position and
   velocity) is recomputed in Python as
   ``grid_node_distance = pore_size + 2 * fiber_radius`` -- this, not the Fortran
   ``dist``, is what downstream plots use.

Because the pipeline always recomputes ``binding_sites`` from ``fiber_radius``,
running through the ``lysis`` CLI (``run-micro`` / ``run-macro``) guarantees the
consistency that the Fortran itself does not enforce.

Historical note
===============

The two scales did not always agree on how fiber size was represented. The
microscale Fortran has stored an explicit ``radius`` (in microns) since 2017,
while the macroscale Fortran originally baked the diameter into a hardcoded
``dist = 1.0862`` micron grid spacing with no fiber-size variable at all. The
Python parameters adopted Pint units and the computed
``grid_node_distance = pore_size + 2 * fiber_radius`` in early 2024; the
macroscale Fortran was converted to take a ``radius`` input and compute
``dist = delx*10000 + radius*2`` in mid-2025, finally matching the Python
definition. (A 2024 fix also corrected experiment-setup code that had specified
fiber sizes in microns instead of nanometres -- a 1000x error.) The
:doc:`historical_fortran` archive variants predate the conversion and still
carry the old hardcoded ``dist``.

Values observed in the data archive
===================================

*This section records the result of a one-time audit (2026-05-28) of the local*
``data/`` *archive: every* ``.txt`` *log,* ``.log`` *file,* ``.json`` *parameter
file, and* ``.h5`` *attribute set was scanned for* ``radius``/``diameter``
*references (773 references across 250 files). It documents which fiber sizes
actually appear in stored runs and a few data-quality issues to be aware of.*

Every stored value is a **radius**. The only **diameter** appearances are
Fortran's grid-spacing arithmetic (``radius*2``) and its microscale diagnostic
print (``Fiber diameter: ... nm``). The distinct ``fiber_radius`` values found,
mapped to the standard fiber types:

==================  ============  ===========  ============
Stored value        = Diameter    Fiber type   Sources
==================  ============  ===========  ============
0.02875 micron      57.4 nm       Q1           json, h5
0.03635 micron      72.7 nm       Q2           json, h5, logs
0.05255 micron      105.1 nm      TF-v/vii/x   json, h5, logs
0.0615 micron       123.0 nm      TB-xi/xiii   json, h5, logs
0.0727 micron       145.4 nm      Q4           json, h5
72.7 nanometer      145.4 nm      Q4           h5 attr (one run)
==================  ============  ===========  ============

``fibrinogen_radius`` (``0.0012 micron``) and ``protofibril_radius``
(``0.0024 micron``) are constant across every run. The microscale Fortran logs
also print the size as a **diameter in nm** (``105.1 nm``, ``123.0 nm`` -- twice
the radius) alongside the node count.

.. note::

   File timestamps are **not** reliable creation dates. Several runs were
   re-converted years after they executed, so a file's ``mtime`` (and even its
   ``birth`` time) reflects the last conversion, not the original simulation.
   The run-code folder name (``YYYY-MM-DD-HHMM``) is the dependable nominal
   creation date; the file ``mtime`` is "last modified/converted".

Known data-quality issues
-------------------------

* **Runs** ``2024-09-02-1411`` **and** ``2024-09-02-1412`` -- the microscale
  log records ``Setting radius = 61.5``, i.e. the Fortran binary was invoked
  with a radius of **61.5 microns** (the nanometre magnitude dropped into a
  micron field -- 1000x too large), even though the re-converted
  ``params.json``/``.h5`` for those runs show the correct ``0.0615 micron``.
  These runs executed on 2024-09-02, the day *before* the units fix, so the
  buggy invocation survives in the original Fortran log. Treat their simulation
  output with caution.
* **Run** ``2024-09-19-1420`` -- ``params.json`` is empty/corrupt (JSON parse
  error) and its parameters cannot be read.
* **Serialisation variation** -- ``2026-02-28-1907.h5`` stores ``fiber_radius``
  in nanometres (``72.7 nanometer``) rather than microns like every other file.
  This is equivalent (= ``0.0727 micron``, Q4), just a different unit on the
  stored ``Quantity``.

.. seealso::

   :doc:`ontology` for the definitions of *Scenario*, *Edge grid*, and related
   terms, and :doc:`fortran_microscale` / :doc:`fortran_macroscale` for the
   full Fortran-to-Python parameter references.
