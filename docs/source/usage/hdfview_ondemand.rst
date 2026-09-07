======================================================
Browsing Lysis Data with HDFView on Open OnDemand
======================================================

This guide takes you from "I have a Buddy account" to "I can see what's
inside a Lysis HDF5 file" — its tree of groups and datasets, their shapes,
and the parameters and provenance stamped on them as attributes — without
writing any code. It covers launching HDFView through Buddy's Open OnDemand
web portal.

When To Use This
-------------------

HDFView is for **browsing and eyeballing** a file you don't already know the
shape of: what groups exist, what a dataset's shape and type are, what an
attribute says. It's the fastest way to answer "what's actually in this
file?" before you write a line of analysis code against it.

It is not the right tool for analysis itself. Once you know what a file
contains and want to compute something from it:

- For statistics and plots over one or more files, use
  :doc:`rstudio_ondemand` (R + ``hdf5r``) or Python + ``h5py`` directly.
- For quick summary numbers without writing any analysis code at all, the
  ``lysis`` command line's ``micro-stats`` / ``macro-stats`` commands may
  already do what you need — see :doc:`statistics_cli`.

Before You Start
-------------------

This guide needs no local software installation and no clone of the
``lysis`` repository — HDFView runs entirely in your browser, against
absolute paths under ``/shared/lysis-group/``. You will need:

- A Buddy HPC account. If you don't have one yet, see below.
- Membership in the ``lysis-group`` Unix group — this is what grants read
  access to ``/shared/lysis-group/``. Check whether you're already a member
  from a Buddy terminal (e.g. a terminal opened through OnDemand):

  .. code-block:: bash

      groups

  and look for ``lysis-group`` in the output.

  .. important::

     If it's missing, **do not email hpc@uco.edu directly to request
     access.** Access is authorised by **Dr. Bannish (Lysis Group PI)**, who
     contacts the UCO HPC group on your behalf. Ask Dr. Bannish first; once
     she has requested access for you, any follow-up can go directly between
     you and the HPC group.
- A web browser. Nothing else — no SSH key, no local install.

  If you *also* want the ``lysis`` command-line tool (for its
  ``micro-stats``, ``macro-stats``, or ``parameters`` commands), that
  additionally requires an SSH key and a repository clone — see
  :doc:`system_onboarding` for that setup. It is not needed for anything in
  this guide.

Launching the HDFView App
----------------------------

1. Go to the Open OnDemand portal, https://ondemand.hpc.uco.edu, and sign in
   with your Buddy account credentials.

2. From the top menu bar, open the **Interactive Apps** menu, then choose
   **HDFView**, listed under the **Data** category (alongside OpenRefine,
   ParaView, QGIS, and QualCoder). This app launches an HDFView interactive
   instance on Buddy — you interact with the HDFView GUI through a VNC
   session in your browser.

   .. note::

      This is a different location from RStudio Server, which is listed
      under **Mathematics** (see :doc:`rstudio_ondemand`).

3. On the launch form:

   - **HDFView Version** — leave this at ``default``. The dropdown only
     offers two entries, ``default`` and ``3.4.1-GCC-14.3.0-Java-21``, and
     ``default`` currently resolves to that same build. Worth knowing:
     3.4.1 is also the version of the ``HDFView-3.4.1-Windows.msi``
     installer the Lysis group has been distributing for local use, so what
     this guide describes matches what you'd see with a local install too.
   - **Additional modules** — leave this blank. HDFView already loads the
     HDF5 and HDF4 libraries it needs; this field is only for loading some
     *other*, unrelated LMod module into the same session, which this guide
     doesn't need. Don't load anything by hand yourself — the app takes care
     of the HDF5/HDF4 setup.
   - **Queue** — leave this at its default, ``general``, for normal
     browsing. If you know you're about to open a very large file (see
     :ref:`hdfview-troubleshooting`), consider **high-mem** instead, which
     gives the same 16 cores but 256 GB of RAM rather than 64 GB.
   - **Number of hours** — the form accepts 1-48 and defaults to **2**.

     .. important::

        This default is much shorter than RStudio's 12-hour default. Raise
        it now if you expect to be browsing for a while — you can't change
        it after launch, only cancel and relaunch.

   Click **Launch**.

4. Your session is submitted as a batch job and briefly enters a **queued**
   state before becoming **running** — on Buddy this is usually well under a
   minute. Track it on the **My Interactive Sessions** page.

5. Once **running**, the session's card shows the compute node (**Host**),
   start time (**Created at**), time left before the wall time expires
   (**Time Remaining**), a **Session ID**, and a **Cancel** button — the
   same fields RStudio's session card uses (see :doc:`rstudio_ondemand`).
   Where RStudio's card has a **Connect to RStudio Server** button, this
   card's button reads **Launch HDFView**; click it to open the HDFView
   desktop in a new browser tab.

   Your session ends automatically when its wall time runs out. A
   **Cancel** button on the session card ends it early.

Opening a Lysis File
-----------------------

.. warning::

   **Everything under** ``/shared/lysis-group/`` **is read-only, and
   HDFView is the one tool in this documentation set that can write.**
   Every other guide here (RStudio/``hdf5r``, Python/``h5py``) only reads
   Lysis data by construction. HDFView's own description is blunter: it
   "lets you view and **modify** datasets, attributes and images." If you
   open a shared file in its default editing mode and save, you can corrupt
   data that someone else's published results depend on — with no undo.
   Treat every step below as mandatory, not optional.

Set HDFView to read-only mode, every session, before opening anything
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. From the HDFView menu bar, open **Tools → Preferences...**.
2. The **General Settings** page opens by default (it's also the first
   entry in the list on the left, alongside **HDF Settings** and **View
   Modules Settings**).
3. Find the **Default File Access Mode** control and select **Read Only**
   (the other option is **Read/Write**).
4. Click **Apply and Close**.

Do this at the start of every HDFView session. HDFView's own documentation
describes the out-of-the-box behaviour for this setting inconsistently even
across its own chapters, so don't assume it's already set the way you left
it last time — check it explicitly.

As a second line of defense, HDFView also lets you choose the access mode
per file at open time: its file-open dialog offers **Open As Read-Only**
and **Open As Read/Write** as separate choices. Prefer **Open As Read-Only**
for anything under ``/shared/lysis-group/`` even if you've already set the
global preference above.

Navigating to the shared data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use **File → Open** (with read-only mode set as above) and browse to
``/shared/lysis-group/experiments/austin-runs/raw/``. For what's actually
in each subfolder — which runs, which parameters — see the README at
``/shared/lysis-group/experiments/austin-runs/README.rst`` (a plain-text
file, not part of this Sphinx build; open it in any text editor or HDFView's
own file browser).

.. note::

   Filesystem permissions under ``/shared/lysis-group/`` give **partial**,
   inconsistent coverage — some files there are already mode-protected
   against writes (e.g. ``-r--r--r--`` files inside a ``dr-xr-xr-x``
   directory, which blocks even renaming), so a write might simply fail
   with a permission error. Others are not locked down nearly as tightly —
   for example, the files under
   ``experiments/austin-runs/raw/austin-runs-corrected/`` are
   individually read-only (``-r--r-----``) but sit in a directory that
   ``lysis-group`` members can still write to, so renaming or deleting a
   file there is not blocked by permissions alone. **Don't treat filesystem
   permissions as your safety net.** They're inconsistent, and the whole
   point of the read-only steps above is to not need them.

Reading What You See
------------------------

HDFView's tree view shows the raw on-disk HDF5 structure, which uses
different names from this project's :doc:`ontology`. Mapping one to the
other:

- An HDF5 **group** is a :term:`Data Collection` — but the on-disk group
  names (``micro_data``, ``macro_data``, ``log_files``) are **not** the same
  as the Data Collection names used in :doc:`ontology` and
  :doc:`data_specification` (``microscale_out``, ``macroscale_in``,
  ``macroscale_out``). The former are literal HDF5 group names; the latter
  are the conceptual categories they implement.
- An HDF5 **dataset** is a :term:`Dataset` in the project sense — one table
  of data, one row per Simulation for microscale datasets.
- HDF5 **attributes**, shown in HDFView's metadata panel under the
  **Object Attribute Info** tab when an object is selected in the tree, hold
  two kinds of project-specific information:

  - **Scenario and Mechanism parameters** — attached to the ``micro_data``
    and ``macro_data`` groups (e.g. ``fiber_radius``, ``micro_seed``,
    ``macro_version``). See :doc:`data_specification` for what each one
    means.
  - **Provenance stamps** — recording when, where, and with what code the
    data was produced. The file-root attributes, visible by selecting the
    file itself (the top node in the tree) rather than any group inside it,
    always include ``dataspec_version`` (which format version the file is
    in) and, for a file converted from an older format, ``converted_from``.
    See :doc:`data_specification`'s "Provenance attributes" section for the
    rest (the ``init_*``, ``pipeline_*``, and ``backend_*`` families, on the
    ``micro_data`` / ``macro_data`` groups).

A Guided First Look
-----------------------

The following walks through one microscale file and one file with
macroscale structure. Both paths and everything reported about their
contents below were verified with ``h5py`` in read-only mode before writing
this guide.

Microscale file
~~~~~~~~~~~~~~~~~~

Open, read-only:

::

    /shared/lysis-group/experiments/austin-runs/raw/austin-runs-corrected/2131.h5

In the tree, expand it. You'll see two groups: ``log_files`` and
``micro_data``. Expand ``micro_data`` — it holds 8 datasets, each of length
50,000 (one entry per Simulation): ``fiber_degraded``, ``pli_first_time``,
``pli_generated_num``, ``sim_final_time``, ``tpa_final_num``,
``tpa_leaving_time``, ``tpa_unbound_by_pli``, ``tpa_unbound_kinetic``.

Click the file's root node and check the **Object Attribute Info** tab: you
should see ``dataspec_version = v2.0.0`` and ``converted_from = v1.99.0`` —
this file was converted from an earlier format version.

Click the ``micro_data`` group itself (not a dataset inside it) and check
its attributes: among others, ``fiber_radius``, ``micro_seed``,
``micro_simulations = 50000``, and ``backend_type = fortran`` — the Scenario
and Mechanism parameters and provenance stamps this Run used.

A file with macroscale structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

   The austin-runs data itself is entirely microscale (no ``macro_data``
   group at all in those files). The example below uses a file from a
   **different** Lysis dataset to show the macroscale shape.

Open, read-only:

::

    /shared/lysis-group/bpaynter/data/lysis-front-pre-lat/Q1.h5

Expand ``macro_data``: instead of datasets directly, you'll see ten
subgroups, ``sim_00`` through ``sim_09`` — one per macroscale Simulation.
Expand any one of them to find its datasets: ``fiber_degrade_time``,
``snapshot_time``, ``tpa_bind_events``, ``tpa_location_snapshot``, and
``tpa_transit_time``. Click the ``macro_data`` group itself to see its
attributes, including ``macro_version = diffuse_into_and_along`` (the
Mechanism used) and ``total_time = 0 second`` — the project's sentinel for
"run to completion" rather than a fixed simulated duration.

This file also has its own ``micro_data`` group (the microscale run that
fed it), so you can repeat the microscale steps above on it too.

Exporting
------------

For a quick look at a dataset's actual values beyond what the tree and
metadata panel show, double-click a dataset to open it in HDFView's table
viewer, then use **Import/Export Data → Export data to Text File** to save
it to a plain text file. The delimiter used is whatever's set in
**Tools → Preferences... → General Settings → Data Delimiter** (tab by
default; change it to comma there first if you want a true CSV).

This writes a *new* file to wherever you choose — it doesn't touch the
source file, so it's compatible with having opened the file read-only.
Anything beyond a one-off export like this — combining several files,
computing summary statistics, plotting — belongs in :doc:`rstudio_ondemand`
or Python/``h5py`` rather than HDFView.

.. _hdfview-troubleshooting:

Troubleshooting
-------------------

Session won't start / stays queued
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Buddy is a shared cluster; a queued session is waiting for resources. On
Buddy, interactive sessions typically start within seconds. If a session
stays queued for several minutes or longer, check its status on **My
Interactive Sessions**, and contact hpc@uco.edu if it persists.

The VNC screen is blank, very slow, or oddly scaled
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Try reconnecting from the **My Interactive Sessions** page — click
**Launch HDFView** again to reopen the VNC tab. If the window seems the
wrong size for your screen, resizing your browser window before
reconnecting can help.

"Permission denied" opening anything under ``/shared/lysis-group/``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You are not a member of the ``lysis-group`` Unix group. Run ``groups`` in a
terminal to confirm, then see "Before You Start" above.

A file won't open at all
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make sure you're pointed at an actual ``.h5`` file, not a directory or a
``.slurm`` log. If the open file dialog isn't showing the file you expect,
check its filter is set to show HDF5 files (or "All Files").

I opened a shared file read/write by mistake
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Close it **without saving** (don't use **File → Save**), then reopen it
using the read-only steps above. If you already saved a change, stop and
tell the Lysis maintainers which file and what you changed — don't try to
undo it yourself, since HDFView has no undo across a save.

A large file is very slow to expand in the tree
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some Lysis files are large — for example,
``/shared/lysis-group/bpaynter/data/lysis-front/TB-ix__1_582_867.h5`` is
about 7.7 GB. Expanding a large group or opening a large dataset's table
view can take a while. If it's consistently unworkable, relaunch with the
**high-mem** queue (see "Launching the HDFView App" above), or switch to
reading just the slice you need with :doc:`rstudio_ondemand` or ``h5py``
instead of browsing the whole thing visually.

See Also
-----------

- :doc:`rstudio_ondemand` — for actually computing statistics over Lysis
  data, once you know what's in the file.
- :doc:`statistics_cli` — command-line summary statistics, if you don't
  need to write any analysis code at all.
- :doc:`ontology` — the terms (Run, Scenario, Mechanism, Dataset, Data
  Collection, ...) used throughout this guide.
- :doc:`data_specification` — the full description of every dataset and
  attribute in a Lysis HDF5 file.
- ``/shared/lysis-group/experiments/austin-runs/README.rst`` — describes the
  austin-runs dataset itself (which runs, which parameters); not part of
  this Sphinx build.
