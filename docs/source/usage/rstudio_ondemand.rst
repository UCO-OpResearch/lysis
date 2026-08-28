======================================================
Analysing Lysis Data in RStudio on Open OnDemand
======================================================

This guide takes you from "I have a Buddy account" to "I have a Lysis HDF5
dataset open as a data frame in RStudio and can start doing statistics." It
covers launching RStudio through Buddy's Open OnDemand web portal and reading
the project's HDF5 output files with the R package ``hdf5r``.

.. note::

   This guide uses angle-bracketed placeholders like ``<_username_>`` for
   values you must fill in with your own — replace the whole token, brackets
   and underscores included, with the same value everywhere it appears.

   A token written ``<_TODO: ...__>`` instead marks a detail this guide does
   not have yet (it needs a maintainer to fill it in, usually from a
   screenshot of the live site). If you hit one, treat the surrounding steps
   as provisional and ask the Lysis maintainers for the current value.

.. _before-you-start-group:

Before You Start
-----------------

You will need:

- A Buddy HPC account. If you don't have one yet, or haven't confirmed you
  can SSH in and clone the repository, see :doc:`system_onboarding` first —
  this guide assumes that part is done.
- Membership in the ``lysis-group`` Unix group. That group membership is
  what grants read access to ``/shared/lysis-group/``, where the Lysis
  project's shared data lives. Check whether you're already a member by
  running, from a Buddy terminal (e.g. an SSH session, or a terminal opened
  through OnDemand):

  .. code-block:: bash

      groups
      # or
      id

  and looking for ``lysis-group`` in the output. If it's missing,
  ``<_TODO: who should a student contact to be added to lysis-group -- is it
  hpc@uco.edu (Buddy's general shared-directory request process), or the
  Lysis project maintainer directly?__>``
- A web browser. No local software installation is required — everything in
  this guide runs on Buddy through your browser.

Launching the RStudio App
--------------------------

1. Go to the Open OnDemand portal:

   https://ondemand.hpc.uco.edu

   Sign in with your Buddy account credentials.

   ``<_TODO: Screenshot 1 -- the OnDemand dashboard landing page, showing the
   top menu bar, so this step can describe exactly what a first-time user
   sees.__>``

2. From the top menu bar, open the **Interactive Apps** menu and find the
   RStudio tile.

   ``<_TODO: Screenshot 2 -- the expanded Interactive Apps dropdown, showing
   every app tile and its exact label, so this guide can name the RStudio
   tile precisely (there may be more than one RStudio entry).__>``

3. On the RStudio launch form, you will be asked to choose an RStudio
   **version**. Choose exactly:

   ::

       RStudio 2024.12.1 (R 4.5.0) + tidyverse

   .. important::

      Pick this version specifically. It is the version the Lysis
      maintainers have confirmed carries the ``hdf5r`` R package, which this
      guide depends on. Other RStudio versions listed in the form may not
      have ``hdf5r`` available — if you pick a different one, ``library(hdf5r)``
      may fail (see :ref:`rstudio-troubleshooting`).

   The rest of the form asks for session resources — typically a queue or
   partition, the number of hours to reserve, and the number of CPU cores.
   Sensible values for browsing and summarising Lysis HDF5 files (as opposed
   to running new simulations) are modest, but this guide does not yet have
   confirmed field names or recommended values for Buddy's form.

   ``<_TODO: Screenshot 3 -- the full RStudio launch form with every field
   visible (queue/partition, wall time, cores, account, and any other
   fields), so each can be documented by its real name along with a sensible
   default for a data-analysis (not simulation) workload.__>``

4. Click **Launch**. Your session is submitted as a batch job and starts in
   a **queued** state, moves to **starting**, and finally **running** once a
   compute node is assigned. You can watch this on the **My Interactive
   Sessions** page (also reached from the top menu bar).

   ``<_TODO: Screenshot 4 -- the My Interactive Sessions page showing one
   queued and one running session, so this guide can describe roughly how
   long a student should expect to wait and what each state looks like on
   Buddy specifically.__>``

5. Once the session is **running**, its session card shows a **Connect to
   RStudio Server** button (or similar) plus the time remaining before the
   session's wall time expires. Click it to open RStudio in a new browser
   tab.

   ``<_TODO: Screenshot 5 -- the RStudio session card itself, showing the
   Connect button, the time-remaining indicator, and the delete/cancel
   control, so this step can name them exactly.__>``

   Your session ends automatically when its requested wall time runs out —
   save your work before then. You can also end it early from the session
   card.

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
       h5f <- H5File$new("/shared/lysis-group/wpumphrey/austin-runs-corrected/some_run.h5",
                          mode = "r")
       write.csv(my_summary, "~/some_run_summary.csv")

       # INCORRECT — never do this against /shared/lysis-group/
       h5f <- H5File$new("/shared/lysis-group/wpumphrey/austin-runs-corrected/some_run.h5",
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

    path <- "/shared/lysis-group/wpumphrey/austin-runs-corrected/KdtPAnoplg_M4.h5"
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

    path <- "/shared/lysis-group/wpumphrey/austin-runs-corrected/KdtPAnoplg_M4.h5"
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

The austin-runs data lives under ``/shared/lysis-group/``, split across
three folders: ``austin_segrest/old_data_compiled``,
``austin_segrest/austin-old-data-imported``, and
``wpumphrey/austin-runs-corrected`` (the folder used in the worked example
above). For what's actually in each of these — which runs, which
parameters, how they relate to each other — see the README at
``/shared/lysis-group/austin-runs.rst``, which describes the dataset in
detail; this guide only covers how to read the file format.

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
become free. Check the **My Interactive Sessions** page for its current
status. ``<_TODO: is there a typical/expected queue time on Buddy worth
mentioning here, or a specific queue/partition known to start faster for
short interactive sessions?__>``

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
in this guide's worked example (see :doc:`data_specification`). If reading a
whole dataset with ``[]`` exhausts your session's memory, read a subset
instead (see "Reading a Dataset" above), or request a session with more
memory on the launch form.

See Also
-----------

- :doc:`ontology` — the terms (Run, Scenario, Mechanism, Dataset, Data
  Collection, ...) used throughout this guide and the rest of the Lysis
  documentation.
- :doc:`data_specification` — the full description of every dataset and
  attribute in a Lysis HDF5 file.
- ``<_TODO: link to the Buddy Desktop / HDFView guide once merged__>`` — a
  complementary guide for browsing Lysis HDF5 files visually with HDFView,
  useful alongside this one for spot-checking a file's structure before
  writing R code against it.
