======================================================
Analysing Lysis Data in RStudio on Open OnDemand
======================================================

This guide takes you from "I have a Buddy account" to "I have a Lysis HDF5
dataset open as a data frame in RStudio and can start doing statistics." It
covers launching RStudio through Buddy's Open OnDemand web portal and reading
the project's HDF5 output files with the R package ``hdf5r``.


.. _before-you-start-group:

Before You Start
-----------------

This guide needs no local software installation and no clone of the
``lysis`` repository — everything in it runs against absolute paths under
``/shared/lysis-group/`` through your browser. You will need:

- A Buddy HPC account. If you don't have one yet, see below — it's arranged
  through the same process as ``lysis-group`` membership.

  If you *also* want to use the ``lysis`` command-line tool (for example its
  ``micro-stats``, ``macro-stats``, or ``parameters`` commands) rather than
  working purely through RStudio, that additionally requires an SSH key and
  a clone of the repository — see :doc:`system_onboarding` for that setup.
  It is not needed for anything in this guide.
- Membership in the ``lysis-group`` Unix group. That group membership is
  what grants read access to ``/shared/lysis-group/``, where the Lysis
  project's shared data lives. Check whether you're already a member by
  running, from a Buddy terminal (e.g. an SSH session, or a terminal opened
  through OnDemand):

  .. code-block:: bash

      groups
      # or
      id

  and looking for ``lysis-group`` in the output.

  .. important::

     If it's missing, **do not email hpc@uco.edu directly to request
     access.** Access is authorised by **Dr. Bannish (Lysis Group PI)**, who
     contacts the UCO HPC group at hpc@uco.edu on your behalf. Ask Dr.
     Bannish first; once she has requested access for you, any follow-up
     (account issues, questions) can go directly between you and the HPC
     group.
- A web browser. No local software installation is required — everything in
  this guide runs on Buddy through your browser.

Launching the RStudio App
--------------------------

1. Go to the Open OnDemand portal:

   https://ondemand.hpc.uco.edu

   Sign in with your Buddy account credentials. The dashboard's top menu bar
   has **Apps**, **Files**, **Jobs**, **Clusters**, **Interactive Apps**,
   **My Interactive Sessions**, and **All Apps** menus, plus **Help** and
   your logged-in username on the right.

2. From the top menu bar, open the **Interactive Apps** menu, then choose
   **RStudio Server**, listed under the **Mathematics** category.

3. On the RStudio Server launch form:

   - **RStudio/R Version** — choose exactly:

     ::

         RStudio 2024.12.1 (R 4.5.0) + tidyverse

     .. important::

        Pick this version specifically. It is the version the Lysis
        maintainers have confirmed carries the ``hdf5r`` R package, which
        this guide depends on. Other RStudio versions listed in the dropdown
        may not have ``hdf5r`` available — if you pick a different one,
        ``library(hdf5r)`` may fail (see :ref:`rstudio-troubleshooting`).

   - **Additional modules** — leave this blank. It's for loading extra
     Lmod modules beyond what the chosen RStudio version already provides,
     which you don't need for this guide.
   - **Queue** — leave this at its default, ``general``. That's a sensible
     choice for reading and summarising data, as opposed to running new
     simulations. If you need more memory than that gives you, see
     :ref:`rstudio-troubleshooting`.
   - **Number of hours** — the form accepts 1-48 and defaults to ``12``.
     The default is more than enough for a typical working session; lower
     it if you know you'll finish sooner.

   There is no separate field for CPU cores or account on this form — Buddy
   fixes those for the RStudio Server app (the resulting session runs with a
   fixed core count; see step 5).

   Click **Launch**.

4. Your session is submitted as a batch job and briefly enters a
   **queued** state before becoming **running**. On Buddy, interactive
   sessions typically start very quickly — often in well under a minute —
   so you may not see the queued state at all. You can watch its status on
   the **My Interactive Sessions** page (also reached from the top menu
   bar).

5. Once the session is **running**, its card on **My Interactive Sessions**
   shows the compute node it's running on (**Host**), when it started
   (**Created at**), how much time remains before its wall time expires
   (**Time Remaining**), and a **Connect to RStudio Server** button. Click
   it to open RStudio in a new browser tab.

   Your session ends automatically when its requested wall time runs out —
   save your work before then. A **Cancel** button on the session card ends
   it early.

   .. tip::

      When you're done, shut down RStudio itself first — the red power
      button in the top-right corner of the RStudio window — *before*
      clicking Cancel on the session card (or letting the wall time expire).
      This isn't required, but it avoids a harmless workspace-restore error
      message on your next launch; see :ref:`rstudio-troubleshooting`.

Orienting Yourself in RStudio
-------------------------------

Once connected, you have a normal RStudio Server IDE: Console, Environment,
Files, and Plots panes, running on a Buddy compute node.

Confirm your working directory from the R console:

.. code-block:: r

    getwd()

Reach the shared Lysis data either from the console:

.. code-block:: r

    list.files("/shared/lysis-group/")

or from the **Files** pane, using its "Go To Folder" option to navigate to
``/shared/lysis-group/``.

.. warning::

   **Everything under** ``/shared/lysis-group/`` **is read-only.** Never open
   a file there in a writable mode, and never write output files into that
   tree. Open HDF5 files with ``mode = "r"`` only, and save anything you
   produce (CSVs, plots, RDS files, ...) to your own home directory instead.

   .. code-block:: r

       # Correct — read-only open, output goes to your home directory
       h5f <- H5File$new("/shared/lysis-group/experiments/austin-runs/raw/austin-runs-corrected/some_run.h5",
                          mode = "r")
       write.csv(my_summary, "~/some_run_summary.csv")

       # INCORRECT — never do this against /shared/lysis-group/
       h5f <- H5File$new("/shared/lysis-group/experiments/austin-runs/raw/austin-runs-corrected/some_run.h5",
                          mode = "r+")   # write-capable mode on shared data

Reading Lysis HDF5 Files with hdf5r
--------------------------------------

``hdf5r`` is an R interface to the HDF5 C library. This section covers the
handful of calls you need to browse and read a Lysis output file. For what
the groups and datasets inside a Lysis HDF5 file actually *mean*, see
:doc:`data_specification` and :doc:`ontology` rather than this guide — this
section only covers the R mechanics.

Loading the Library
~~~~~~~~~~~~~~~~~~~~

.. code-block:: r

    library(hdf5r)
    packageVersion("hdf5r")

If ``library(hdf5r)`` fails with "there is no package called 'hdf5r'", you
most likely launched a different RStudio version — see
:ref:`rstudio-troubleshooting`.

Opening a File Read-Only
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: r

    path <- "/shared/lysis-group/experiments/austin-runs/raw/austin-runs-corrected/KdtPAnoplg_M4.h5"
    h5f <- H5File$new(path, mode = "r")

``mode = "r"`` is mandatory for anything under ``/shared/lysis-group/`` — see
the warning above.

Listing Contents
~~~~~~~~~~~~~~~~~

.. code-block:: r

    h5f$ls(recursive = TRUE)
    # or, just the names at the top level:
    names(h5f)

``ls(recursive = TRUE)`` returns a data frame with one row per group and
dataset, including each dataset's shape (``dataset.dims``) and HDF5 type
class.

Note that the group names you'll see here (``micro_data``, ``macro_data``,
``log_files``) are not the same as the Data Collection names used in
:doc:`data_specification` and :doc:`ontology` (``microscale_out``,
``macroscale_in``, ``macroscale_out``) — the latter are conceptual
categories, the former are the literal on-disk HDF5 group names.

Reading a Dataset
~~~~~~~~~~~~~~~~~~

Read a whole dataset with ``[]``:

.. code-block:: r

    fiber_degraded <- h5f[["micro_data/fiber_degraded"]][]

For a large dataset, read only the slice you need instead of pulling the
whole thing into memory — index like a normal R array:

.. code-block:: r

    first_five <- h5f[["micro_data/tpa_leaving_time"]][1:5]

This matters most for the macroscale collections, where some datasets carry
one row per snapshot across a whole simulation; see :doc:`data_specification`
for which datasets are large.

Macroscale datasets are also nested one level deeper than microscale ones:
each macroscale simulation gets its own subgroup, ``macro_data/sim_00``,
``macro_data/sim_01``, and so on. List the simulations present in a file
with ``names(h5f[["macro_data"]])``, then read a dataset from one of them
with, e.g., ``h5f[["macro_data/sim_00/snapshot_time"]][]``.

.. note::

   **The austin-runs data itself is entirely microscale** — every file
   under ``experiments/austin-runs/raw/austin-runs-corrected`` and
   ``experiments/austin-runs/raw/austin-old-data-imported`` (320 ``.h5``
   files, all checked) has only ``micro_data`` and ``log_files``, no ``macro_data``
   group at all. The macroscale shape above is real, but it was verified
   against a file from a **different** Lysis dataset:
   ``/shared/lysis-group/bpaynter/data/lysis-front-pre-lat/Q1.h5``, opened
   read-only. That file has ``macro_data/sim_00`` through
   ``macro_data/sim_09`` (10 simulations), each holding
   ``fiber_degrade_time``, ``snapshot_time``, ``tpa_bind_events``,
   ``tpa_location_snapshot``, and ``tpa_transit_time`` — matching
   :doc:`data_specification`'s macroscale dataset list. If you're working
   specifically with austin-runs files, you will only ever see
   ``micro_data``.

Reading Attributes
~~~~~~~~~~~~~~~~~~~~

Lysis HDF5 files carry provenance and scenario/mechanism parameters as HDF5
attributes, at the file root and on the ``micro_data`` / ``macro_data``
groups (see :doc:`data_specification`, "Provenance attributes").

.. code-block:: r

    # All attributes on an object, as a named list
    h5attributes(h5f)                    # root-level (file-wide) attributes
    h5attributes(h5f[["micro_data"]])    # scenario/mechanism parameters

    # A single named attribute
    h5attr(h5f[["micro_data"]], "snap_proportion")

    # Just the attribute names, without reading their values
    h5attr_names(h5f[["micro_data"]])

Chunking and Compression
~~~~~~~~~~~~~~~~~~~~~~~~~~

Lysis HDF5 datasets are stored chunked and gzip-compressed. You don't need
to do anything differently because of this — ``hdf5r`` decompresses
transparently on read — but it's why a dataset's on-disk size can be much
smaller than its in-memory size once loaded.

Closing the File
~~~~~~~~~~~~~~~~~~

.. code-block:: r

    h5f$close_all()

Always close a file when you're done with it. An HDF5 file left open by a
crashed or abandoned R session can produce a confusing "file already open"
error the next time you or someone else opens it — see
:ref:`rstudio-troubleshooting`.

Converting to a Data Frame
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most Lysis microscale datasets are simple 1-D arrays, one entry per
simulation, so they combine directly into a data frame:

.. code-block:: r

    df <- data.frame(
      fiber_degraded   = h5f[["micro_data/fiber_degraded"]][],
      tpa_leaving_time = h5f[["micro_data/tpa_leaving_time"]][]
    )
    summary(df)

Worked Example
----------------

The following is a complete, copy-pasteable script that opens one real file
from the austin-runs data, lists its structure, reads two datasets, reads a
parameter attribute, and computes a summary statistic.

.. code-block:: r

    library(hdf5r)

    path <- "/shared/lysis-group/experiments/austin-runs/raw/austin-runs-corrected/KdtPAnoplg_M4.h5"
    h5f <- H5File$new(path, mode = "r")

    # Structure: group/dataset names, shapes, and HDF5 type classes
    h5f$ls(recursive = TRUE)[, c("name", "obj_type", "dataset.dims")]

    # File-wide provenance attributes
    h5attributes(h5f)

    # A scenario parameter from the micro_data group
    h5attr(h5f[["micro_data"]], "snap_proportion")

    # Read two full microscale datasets into a data frame
    df <- data.frame(
      fiber_degraded   = h5f[["micro_data/fiber_degraded"]][],
      tpa_leaving_time = h5f[["micro_data/tpa_leaving_time"]][]
    )
    summary(df)

    h5f$close_all()

.. note::

   **This snippet was actually executed and verified**, using R 4.5.1 with
   ``hdf5r`` 1.3.12 against this exact file on a Buddy login node (via the
   ``R-bundle-CRAN/2025.10-foss-2025a`` module, not the OnDemand RStudio
   session itself, which this guide's author cannot access programmatically).
   The real output was:

   .. code-block:: text

       > h5f$ls(recursive = TRUE)[, c("name", "obj_type", "dataset.dims")]
                                   name    obj_type dataset.dims
       1                      log_files   H5I_GROUP         <NA>
       2            log_files/micro_log H5I_DATASET          100
       3                     micro_data   H5I_GROUP         <NA>
       4      micro_data/fiber_degraded H5I_DATASET        50000
       5      micro_data/pli_first_time H5I_DATASET        50000
       6   micro_data/pli_generated_num H5I_DATASET        50000
       7      micro_data/sim_final_time H5I_DATASET        50000
       8       micro_data/tpa_final_num H5I_DATASET        50000
       9    micro_data/tpa_leaving_time H5I_DATASET        50000
       10 micro_data/tpa_unbound_by_pli H5I_DATASET        50000
       11 micro_data/tpa_unbound_kinetic H5I_DATASET       50000

       > h5attributes(h5f)
       $dataspec_version
       [1] "v2.0.0"
       $converted_from
       [1] "v1.99.0"

       > h5attr(h5f[["micro_data"]], "snap_proportion")
       [1] 0.6666667

       > summary(df)
        fiber_degraded  tpa_leaving_time
        Mode :logical   Min.   :0.0000021
        FALSE:49820     1st Qu.:0.0896290
        TRUE :180       Median :0.2167437
                        Mean   :0.3158810
                        3rd Qu.:0.4368315
                        Max.   :3.9777182

   That is: out of the 50,000 microscale simulations recorded in this file,
   180 (0.36%) reached full fiber degradation before the simulation ended.

Where the Data Is
--------------------

The austin-runs data lives in one place:
``/shared/lysis-group/experiments/austin-runs/``. Its ``raw/`` folder holds
the three data sets — ``austin-runs-corrected`` (the authoritative data, and
the folder used in the worked example above), ``austin-old-data-imported``
(the historical baseline), and ``old_data_compiled`` (the raw Fortran
archive, kept for provenance only). **Start with**
``austin-runs-corrected``.

For what's actually in each of these — which runs, which parameters, how
they relate to each other — see the README at
``/shared/lysis-group/experiments/austin-runs/README.rst``, which describes
the dataset in detail; this guide only covers how to read the file format.

.. note::

   These folders moved on 2026-09-07, out of per-person directories
   (``wpumphrey/``, ``austin_segrest/``) and into ``experiments/``. The old
   paths still work — they are now symlinks to the new locations — so an
   older script or email will not break. Use the new paths in anything you
   write.

hdf5r Reference Pointers
---------------------------

- CRAN package page: https://cran.r-project.org/package=hdf5r
  (current version 1.3.12 at time of writing, matching the version installed
  on Buddy)
- Package vignette ("Introduction to the hdf5r package"):
  https://cran.r-project.org/web/packages/hdf5r/vignettes/hdf5r.html
- Reference manual: https://cran.r-project.org/web/packages/hdf5r/hdf5r.pdf
  (PDF) or https://cran.r-project.org/web/packages/hdf5r/refman/hdf5r.html
  (HTML)
- GitHub repository: https://github.com/hhoeflin/hdf5r
- The HDF Group's own explanation of the HDF5 data model (groups, datasets,
  attributes), for readers who haven't used HDF5 before:
  https://support.hdfgroup.org/documentation/hdf5/latest/_h5_d_m__u_g.html

.. _rstudio-troubleshooting:

Troubleshooting
------------------

Session won't start / stays queued
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Buddy is a shared cluster; a queued session is waiting for resources to
become free. On Buddy, interactive sessions typically start within seconds
(well under a minute), so a session that stays queued for several minutes or
longer is unusual. Check the **My Interactive Sessions** page for its
current status, and if it stays queued for an extended time, contact
hpc@uco.edu.

A workspace-restore error appears when RStudio starts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This happens if a previous session wasn't shut down cleanly — for example,
the session's wall time expired, or it was cancelled from the session card
while RStudio was still running, instead of being shut down from inside
RStudio first (see the tip in "Launching the RStudio App", step 5). It's
generally harmless and can be dismissed; going forward, you can avoid it by
shutting down RStudio itself (the red power button in the top-right corner
of the RStudio window) before ending the session.

``library(hdf5r)`` fails with "there is no package called 'hdf5r'"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You likely launched a different RStudio version. End this session and relaunch,
selecting exactly ``RStudio 2024.12.1 (R 4.5.0) + tidyverse`` on the launch
form.

"Permission denied" reading a file under ``/shared/lysis-group/``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You are not a member of the ``lysis-group`` Unix group. Run ``groups`` in a
terminal to confirm, then see :ref:`Before You Start <before-you-start-group>`.

"File already open" / HDF5 lock errors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This usually means a previous R session (yours, on another tab, or a
crashed one) still holds the file open. Make sure you call
``h5f$close_all()`` when you're done with a file, restart your R session
(Session > Restart R) if a handle seems stuck, and avoid opening the same
file from two R sessions at once.

Running out of memory reading a dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some macroscale datasets are much larger than the microscale datasets shown
in this guide's worked example (see :doc:`data_specification`). There's no
separate memory field on the launch form (see "Launching the RStudio App",
step 3) — try these first:

- Reading a subset instead of the whole dataset (see "Reading a Dataset"
  above).
- Removing large objects you no longer need with ``rm()``, then
  ``gc()``, before reading the next one.
- Processing a large dataset in chunks (read a slice, summarise it, discard
  it, move to the next slice) rather than holding the whole thing in memory
  at once.

If you genuinely need more RAM than that gets you, the **Queue** dropdown
is the mechanism for it: relaunch RStudio with **high-mem** selected
instead of ``general``. It provides nodes with the same 16 cores but
256 GB of RAM instead of the standard 64 GB.

See Also
-----------

- :doc:`ontology` — the terms (Run, Scenario, Mechanism, Dataset, Data
  Collection, ...) used throughout this guide and the rest of the Lysis
  documentation.
- :doc:`data_specification` — the full description of every dataset and
  attribute in a Lysis HDF5 file.
- :doc:`hdfview_ondemand` — a complementary guide for browsing Lysis HDF5
  files visually with HDFView, useful alongside this one for spot-checking a
  file's structure before writing R code against it.
