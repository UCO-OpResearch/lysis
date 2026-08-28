System Onboarding
=================

.. note::

   Throughout this guide, anything written as ``<_like_this_>`` is a
   **placeholder** — you must replace the *whole* token, angle brackets and
   underscores included, with your own value.

   For example, if your HPC username is ``jsmith``, the command

   .. code-block:: bash

       ssh <_username_>@hpc.uco.edu

   becomes

   .. code-block:: bash

       ssh jsmith@hpc.uco.edu

   Leaving the angle brackets in place (running ``ssh <_username_>@hpc.uco.edu``
   literally) is the single most common mistake first-time readers make —
   if a command fails with something like "could not resolve hostname", check
   that you have not accidentally left a placeholder token in it.

Which Path Are You On?
-----------------------

This guide covers three different situations — skim the list below and
follow only the steps that apply to you. Steps and sections marked
**Development path only** can be skipped if you are on the "Analysis +
the ``lysis`` CLI" path.

- **Analysis only** — reading existing data with RStudio/hdf5r or
  HDFView, no local repository needed. Skip this page entirely and go to
  the RStudio or HDFView guide instead (forthcoming in this documentation
  set).
- **Analysis + the** ``lysis`` **CLI** — running commands such as
  micro-stats, macro-stats, or inspecting parameters, without executing
  new Simulations. Follow Steps 1, 2, 4, and 5 below, noting these
  differences:

  - In Step 4, clone over **HTTPS** instead of SSH — the repository is
    public, so no GitHub SSH key is needed. **Skip Step 3.**
  - In Step 5, a bare ``uv sync --frozen`` is enough — you do not need
    ``--extra test``.
  - **Skip Step 6** — these commands read data that already exists and
    never touch the Fortran binaries.

- **Development** — executing new Simulations, or changing the Python or
  Fortran source. Follow the whole page, in order.

This guide walks through setting up your local environment and HPC account
to run simulations on Buddy.

| Buddy is the UCO HPC cluster used for running simulations and computational workloads. Users connect remotely via SSH (Secure Shell — an encrypted remote-login protocol) and develop either directly on the cluster or through a remote IDE workflow.
| *Any IDE setup can be skipped by running VS Code or PyCharm on Buddy through* `Buddy OnDemand <https://ondemand.hpc.uco.edu>`_
|

Note that two separate SSH keys are used in this setup. An **SSH key** is a
matched pair of files (a private half you keep secret, and a public half you
hand out) that lets you log in without typing a password every time — you
generate the pair once, keep the private half where it was created, and
register the public half with whatever you want to connect *to*.

- Your local machine SSH key allows you to connect to Buddy
- Your Buddy SSH key allows Buddy to authenticate with GitHub

Prerequisites
--------------

Before starting, ensure you have:

- An HPC account on Buddy
- Git installed locally
- A GitHub account
- **Development path only.** Your IDE of choice (VS Code used here)
- Membership in the ``lysis-group`` Unix group, which grants read access
  to ``/shared/lysis-group/``. Check with:

  .. code-block:: bash

      groups | grep lysis-group

  If this prints nothing, you are not yet a member — contact your HPC
  administrator to be added before continuing.

Step 1: Configure Local SSH Config
----------------------------------

Configuring your SSH settings streamlines future remote connections to
Buddy. An SSH config file (usually ``~/.ssh/config`` on Linux/macOS, or
``C:\Users\{User}\.ssh\config`` on Windows) lets you type a short name like
``node-<_node_number_>`` instead of the full connection details every time
— you'll see where that node number comes from in `Connect to Buddy`_
below.

See the `VS Code Setup`_ section below for the IDE-specific configuration
instructions.

Step 2: Create Local SSH Key and Add It to Buddy
------------------------------------------------

1. Open a terminal **on your local machine** (not on Buddy).

2. Run:

   .. code-block:: bash

       ssh-keygen -t ed25519

   You should see prompts asking where to save the key and whether to set a
   passphrase.

3. Press Enter at each prompt to accept the default file location and skip
   the passphrase, unless you have a specific reason to set one.

4. Add the newly created public key to your Buddy account, using the
   instructions for your operating system below.

5. On first connection to Buddy, SSH may prompt you to confirm the host
   fingerprint. Type ``yes`` to continue.

Linux or MacOS
^^^^^^^^^^^^^^

Run:

.. code-block:: bash

    ssh-copy-id <_username_>@hpc.uco.edu

Replace ``<_username_>`` with your HPC username.

You should see a summary ending in something like
``Number of key(s) added: 1``, followed by a suggestion to test the login —
you can do that with the command in `Verify SSH Connection`_ below.

Windows
^^^^^^^

From your local PowerShell, run:

.. code-block:: powershell

    cat ~/.ssh/id_ed25519.pub | ssh <_username_>@hpc.uco.edu "mkdir -p ~/.ssh && chmod 700 ~/.ssh && cat >> ~/.ssh/authorized_keys && chmod 600 ~/.ssh/authorized_keys"

This command produces no output on success. You will be prompted for your
Buddy password once (this is the last time you should need it).

Step 3: Create Buddy SSH Key and Add It to GitHub
--------------------------------------------------

**Development path only.** If you're on the "Analysis + lysis CLI" path,
skip this step — Step 4 below shows you how to clone over HTTPS instead,
which needs no GitHub SSH key.

1. From the Buddy terminal, generate a new SSH key:

   .. code-block:: bash

       ssh-keygen -t ed25519

   As in Step 2, press Enter to accept the defaults.

2. Display the public key:

   .. code-block:: bash

       cat ~/.ssh/id_ed25519.pub

   You should see a single line of text starting with ``ssh-ed25519 AAAA...``.
   Copy this whole line.

3. While logged into GitHub in a web browser, navigate to:

   https://github.com/settings/keys

4. Click **New SSH Key**.
5. Enter a title of your choice.
6. Paste the copied public key into the **Key** field.
7. Save the key.

See the `VS Code Setup`_ section below for connecting your IDE to the
GitHub repository.

Step 4: Cloning the GitHub Repository
-------------------------------------

1. From the Buddy terminal, configure your Git identity (this labels the
   commits you make; it does not need to match your GitHub account exactly,
   but should be recognisable to your teammates):

   .. code-block:: bash

       git config --global user.name "<_full_name_>"
       git config --global user.email <_email_address_>

   These commands produce no output on success.

2. Clone the repository. Which method to use depends on your path from
   `Which Path Are You On?`_ above:

   **Analysis + lysis CLI path — clone over HTTPS** (no SSH key needed,
   the repository is public):

   .. code-block:: bash

       git clone https://github.com/UCO-OpResearch/lysis.git

   **Development path — clone over SSH**, using the Buddy → GitHub key
   from Step 3:

   - In a web browser, navigate to the GitHub repository page:
     https://github.com/UCO-OpResearch/lysis
   - Click the green **<> Code** dropdown (this ``<>`` is GitHub's own
     button label, not a placeholder for you to fill in).
   - Select **SSH**.
   - Copy the SSH repository URL shown. It will look like
     ``git@github.com:UCO-OpResearch/lysis.git``.
   - Back in the Buddy terminal, clone:

     .. code-block:: bash

         git clone <_ssh_repository_url_>

   Either way, you should see output starting with ``Cloning into
   'lysis'...``, followed by progress lines (``Receiving objects: ...``,
   ``Resolving deltas: ...``), and ending without an error.

   This creates a new directory named ``lysis`` in your current location.
   The rest of this guide refers to that directory as ``<_project_dir_>``
   (by default this is ``~/lysis``, if you cloned from your home
   directory).

Step 5: Configure ``.bashrc`` and Initialize ``uv``
---------------------------------------------------

1. **Development path only.** Add the Intel compiler module to your
   ``.bashrc`` (skip this on the "Analysis + lysis CLI" path — you won't
   be compiling anything):

   .. code-block:: bash

       echo 'module load intel-compilers/2023' >> ~/.bashrc

   Buddy uses a system called **LMod** to manage optional software
   ("modules"): software is not on your ``PATH`` until you ``module load``
   it. This line loads the Intel Fortran compiler automatically every time
   you open a new terminal on Buddy, since you will need it to build the
   Fortran simulation code in Step 6.

2. Reload the terminal afterward (close and reopen it, or run
   ``source ~/.bashrc``) so the change takes effect.

3. Install ``uv`` by running:

   .. code-block:: bash

       curl -LsSf https://astral.sh/uv/install.sh | sh

   ``uv`` is the Python package and environment manager this project uses
   (in place of Conda, which is no longer used here). Every Python command
   in this project should be prefixed with ``uv run``, e.g. ``uv run lysis
   --help``.

   You should see a message confirming ``uv`` was installed, and a note
   that it added itself to your shell's startup file. As with Step 5.2,
   reload your terminal (or ``source ~/.bashrc``) so the ``uv`` command
   becomes available.

4. Install the project dependencies:

   .. code-block:: bash

       cd <_project_dir_>
       uv sync --frozen --extra test

   If you're on the "Analysis + lysis CLI" path, a bare
   ``uv sync --frozen`` is enough — skip ``--extra test``.

   This creates an isolated Python environment (a **virtual environment**,
   stored in ``.venv/``) and installs the ``lysis`` package into it in
   **editable** mode — meaning changes you pull in with ``git pull`` take
   effect immediately, with no separate reinstall step. ``--frozen``
   installs exactly the versions recorded in ``uv.lock`` instead of
   re-resolving them (the project convention — never omit it), and
   ``--extra test`` additionally installs the packages needed to run the
   test suite.

   You should see a series of ``+ package==version`` lines as dependencies
   download, ending without an error. This can take a minute or two the
   first time.

Step 6: Build the Fortran Binaries
------------------------------------

**Development path only.** Skip this whole step on the "Analysis + lysis
CLI" path — those commands never touch the compiled Fortran binaries.

The Python package does not compile the Fortran simulation code for you.
You must build it yourself before executing a Simulation for the first
time, and again any time the Fortran source changes (for example, after a
``git pull`` that touches ``src/fortran/``).

1. Confirm the Intel compiler module is loaded. It should already be, if
   you completed Step 5.1 and opened a new terminal since:

   .. code-block:: bash

       module load intel-compilers/2023

2. From the project directory, build everything:

   .. code-block:: bash

       cd <_project_dir_>
       make

   You should see a sequence of compiler invocations, one per binary, for
   example:

   .. code-block:: text

       mkdir -p bin
       gcc -std=c99 -c ./src/c/kiss.c -o ./bin/kiss.o
       ifort -r8 -mcmodel medium -traceback -fpe0 -diag-disable=10448 ./bin/kiss.o ./src/fortran/version_stamp.f90 ./src/fortran/micro_rates.f90 -o ./bin/micro_rates
       ifort ... -o ./bin/macro_diffuse_into_and_along__internal
       ifort ... -o ./bin/macro_diffuse_into_and_along__external
       mkdir -p lib
       gcc -std=c99 -fPIC -shared -o ./lib/kiss.so ./src/c/kiss.c

   No output line should contain ``Error`` or ``Fatal Error``. ``make``
   creates the ``bin/`` and ``lib/`` directories itself — you do not need
   to create them first.

3. Confirm the binaries were created:

   .. code-block:: bash

       ls bin/ lib/

   You should see ``micro_rates``,
   ``macro_diffuse_into_and_along__internal``, and
   ``macro_diffuse_into_and_along__external`` under ``bin/``, and
   ``kiss.so`` under ``lib/``.

Keep the Intel compiler module loaded for any later command that executes a
Simulation (``lysis run-micro`` / ``lysis run-macro``), not just for the
build — the compiled binaries need the Intel runtime libraries to execute,
not only to compile.

.. note::

   ``lysis run-micro`` and ``lysis run-macro`` do **not** rebuild the
   Fortran binaries for you, and they do not call ``make`` automatically.
   Instead, each execution checks that the compiled binary you point it at
   (with ``--executable``) matches the current ``src/fortran/`` source
   *before* it starts a Simulation. If the binary is older — for example,
   after a ``git pull`` brought in Fortran changes since your last
   ``make`` — the command stops immediately with an error beginning:

   .. code-block:: text

       Fortran binary at <path> reports build <commit> (<clean/dirty>);
       src/fortran/ is at <commit> (<clean/dirty>).

   This means "your compiled binaries are older than the checked-out
   source." The fix is always to rebuild:

   .. code-block:: bash

       cd <_project_dir_>
       make

   The error message also mentions an override that lets the check be
   skipped; that override exists for advanced historical-reproduction
   workflows (see :doc:`data_specification`) and is not a substitute for
   rebuilding — do not use it to make this error go away. (See also
   `Fortran Binary Is Out Of Date`_ in Troubleshooting.)

Setup Verification
--------------------

Verify the following before continuing:

- You can SSH into Buddy without entering a password
- Group membership check passes: ``groups | grep lysis-group`` prints a
  line containing ``lysis-group``
- You can clone the GitHub repository from Buddy without error (via SSH
  or HTTPS, per your path from `Which Path Are You On?`_)
- ``uv sync --frozen`` (add ``--extra test`` on the Development path)
  completes successfully
- **Development path only.** ``make`` completes with no ``Error`` or
  ``Fatal Error`` lines, and ``bin/`` contains ``micro_rates``,
  ``macro_diffuse_into_and_along__internal``, and
  ``macro_diffuse_into_and_along__external``
- **Development path only.** The test suite runs:
  ``uv run --extra test pytest`` completes (this can take a while — some
  tests execute the compiled Fortran binaries). Do **not** run a bare
  ``uv run pytest`` — since ``pytest`` is an optional extra, a bare
  ``uv run pytest`` re-syncs your environment, drops it, and fails with:

  .. code-block:: text

      error: Failed to spawn: `pytest`
        Caused by: No such file or directory (os error 2)

- **Development path only.** IDE can open a remote SSH session

VS Code Setup
---------------

**Development path only.** Skip this whole section on the "Analysis +
lysis CLI" path.

SSH Configuration
^^^^^^^^^^^^^^^^^^

This configuration will need to be done once on each computer connected to
Buddy.

1. Open VS Code
2. Select the **Open a Remote Window** button in the bottom-left corner of the application
3. Choose **SSH**
4. Select **Configure SSH Hosts...**
5. Choose the SSH config file, typically:

   .. code-block:: text

       C:\Users\{User}\.ssh\config

Add the following configuration, replacing ``<_username_>`` with your HPC
username:

.. code-block:: text

    Host node-*
        ProxyJump hpc.uco.edu
        User <_username_>

    Host hpc.uco.edu
        User <_username_>

Connect to Buddy
^^^^^^^^^^^^^^^^^^

These steps will be done each time you connect to Buddy.

.. code-block:: powershell

    ssh <_username_>@hpc.uco.edu
    salloc --exclusive

- Leave this PowerShell open while using VS Code.
- ``salloc`` requests a compute node allocation from Buddy's job scheduler
  (Slurm). You should see a message like
  ``salloc: Nodes node-<_node_number_> are ready for job``. Note the node
  number allocated to you — you will need it in the next step, and it will
  be different each time you run ``salloc``.

In VS Code, connect to the node

1. Select **Open a Remote Window**
2. Select **Connect to Host**
3. Type ``node-<_node_number_>``, replacing ``<_node_number_>`` with the
   number noted above
4. If it asks for the platform, select **Linux**
5. If it shows *node-<_node_number_> has fingerprint...*, select **Continue**

Once it finishes connecting

1. Select **Open Folder** from the File dropdown.
2. Enter ``<_project_dir_>`` and click **OK**
3. Answer any trust requests that pop up


Troubleshooting
-----------------

SSH Permission Denied
^^^^^^^^^^^^^^^^^^^^^^^

If you receive an SSH permission error:

- Ensure your public key exists in ``~/.ssh/authorized_keys``
- Verify your SSH config uses the correct username
- Confirm your SSH agent is running

Verify SSH Connection
^^^^^^^^^^^^^^^^^^^^^^^

Test your connection with:

.. code-block:: bash

    ssh <_username_>@hpc.uco.edu

Verify GitHub Authentication
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Test GitHub SSH access with:

.. code-block:: bash

    ssh -T git@github.com

Permission Denied When Cloning the Repository
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If ``git clone <_ssh_repository_url_>`` (Step 4, Development path) fails
with a permission or "publickey" error, this almost always means the
**Buddy** SSH key from Step 3 was never added to GitHub, or was added
under the wrong GitHub account. This is a different key from the one in
Step 2 — see the note near the top of this page about the two separate
SSH keys. Check it using `Verify GitHub Authentication`_ above. (If you
cloned over HTTPS on the "Analysis + lysis CLI" path, this does not apply
to you — HTTPS cloning of this public repository needs no key at all.)

``uv: command not found``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``uv`` installer (Step 5.3) adds itself to your shell's startup file,
but that only takes effect in **new** terminal sessions. Close and reopen
your terminal, or run ``source ~/.bashrc``, then try again.

``module: command not found``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``module`` command is provided by Buddy's LMod system and only exists
**on Buddy**, not on your local machine. Make sure you are actually
connected to Buddy (your terminal prompt should show the cluster
hostname) before running ``module load`` commands.

If you are on Buddy and still see this error, your shell may not have
loaded LMod's initialization script; contact HPC support.

Fortran Binary Is Out Of Date
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If ``lysis run-micro`` or ``lysis run-macro`` (Development path) stops
immediately with an error starting ``Fortran binary at ... reports
build ...``, your compiled binaries are older than the checked-out
``src/fortran/`` source — most often because a ``git pull`` brought in
Fortran changes since your last ``make``. See the note in `Step 6: Build
the Fortran Binaries`_ above: the fix is always to rebuild with ``make``,
never the override flag mentioned in the error message.
