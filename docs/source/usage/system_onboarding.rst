System Onboarding
=================

This guide walks through setting up your local environment and HPC account
to run simulations on Buddy.

| Buddy is the UCO HPC cluster used for running simulations and computational workloads. Users connect remotely via SSH and develop either directly on the cluster or through a remote IDE workflow.
| *Any IDE setup can be skipped by running VS Code or PyCharm on Buddy through* `Buddy OnDemand <https://ondemand.hpc.uco.edu>`_
|

Note that two separate SSH keys are used in this setup:

- Your local machine SSH key allows you to connect to Buddy
- Your Buddy SSH key allows Buddy to authenticate with GitHub

Prerequisites
--------------

Before starting, ensure you have:

- An HPC account on Buddy
- Git installed locally
- A GitHub account
- Your IDE of choice (VS Code used here)

Step 1: Configure Local SSH Config
----------------------------------

Configuring your SSH settings streamlines future remote connections to Buddy.

See the IDE-specific configuration instructions below.

Step 2: Create Local SSH Key and Add It to Buddy
------------------------------------------------

Create a local SSH key by running the following command in your local terminal:

.. code-block:: bash

    ssh-keygen -t ed25519

Press Enter to accept the default file location when prompted.
On first connection, SSH may prompt you to confirm the host fingerprint.
Type ``yes`` to continue.

Add the newly created public key to your Buddy account.

Linux or MacOS
^^^^^^^^^^^^^^

Run:

.. code-block:: bash

    ssh-copy-id USERNAME@hpc.uco.edu

Replace ``USERNAME`` with your HPC username.

Windows
^^^^^^^

From your local PowerShell, run:

.. code-block:: powershell

    cat .ssh/id_ed25519.pub | ssh USERNAME@hpc.uco.edu "mkdir -p ~/.ssh && chmod 700 ~/.ssh && cat >> ~/.ssh/authorized_keys && chmod 600 ~/.ssh/authorized_keys"

Step 3: Create Buddy SSH Key and Add It to GitHub
--------------------------------------------------

From the Buddy terminal, generate a new SSH key:

.. code-block:: bash

    ssh-keygen -t ed25519

Then display the public key:

.. code-block:: bash

    cat id_ed25519.pub

Copy the output.

Next, while logged into GitHub:

1. Navigate to:

   https://github.com/settings/keys

2. Click **New SSH Key**
3. Enter a title of your choice
4. Paste the copied public key into the **Key** field
5. Save the key

See the IDE-specific documentation below for connecting your IDE
to the GitHub repository.

Step 4: Cloning the GitHub Repository
-------------------------------------

From the Buddy terminal, configure your Git User

.. code-block:: bash

    git config --global user.name "Your Name"
    git config --global user.email your.email@example.com

1. Navigate to the GitHub repository page
2. Click the green **<> Code** dropdown
3. Select **SSH**
4. Copy the SSH repository URL

Go back to the Buddy terminal

.. code-block:: bash

    git clone <paste SSH repository URL>

Step 5: Configure ``.bashrc`` and Initialize ``uv``
---------------------------------------------------

Add the Intel compiler module to your ``.bashrc``:

.. code-block:: bash

    echo 'module load intel-compilers/2023' >> ~/.bashrc

Reload the terminal afterward.

Install ``uv`` by running:

.. code-block:: bash

    curl -LsSf https://astral.sh/uv/install.sh | sh

Install the project dependencies:

.. code-block:: bash

    cd ~/lysis
    uv sync --extra test

Setup Verification
==================

Verify the following before continuing:

- You can SSH into Buddy without entering a password
- You can clone the GitHub repository from Buddy without error
- ``uv sync --extra test`` completes successfully
- IDE can open a remote SSH session

VS Code Setup
==============

SSH Configuration
------------------

1. Open VS Code
2. Select **Open a Remote Window**
3. Choose **SSH**
4. Select **Configure SSH Hosts...**
5. Choose the SSH config file, typically:

   .. code-block:: text

       C:\Users\{User}\.ssh\config

Add the following configuration, replacing ``USERNAME`` with your HPC username:

.. code-block:: text

    Host node-*
        ProxyJump hpc.uco.edu
        User USERNAME

    Host hpc.uco.edu
        User USERNAME


Troubleshooting
================

SSH Permission Denied
----------------------

If you receive an SSH permission error:

- Ensure your public key exists in ``~/.ssh/authorized_keys``
- Verify your SSH config uses the correct username
- Confirm your SSH agent is running

Verify SSH Connection
----------------------

Test your connection with:

.. code-block:: bash

    ssh USERNAME@hpc.uco.edu

Verify GitHub Authentication
-----------------------------

Test GitHub SSH access with:

.. code-block:: bash

    ssh -T git@github.com
