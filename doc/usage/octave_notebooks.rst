-----------------------
Lysis Model
-----------------------
Octave Notebook Usage
---------------------

First-time setup
+++++++++++++++++

#. Log into https://buddy.uco.edu.

#. Go to "Files", "Home Folder", and then navigate to your lysis code folder.

#. Click on the "Open in Terminal" button.

#. Run the command `git branch` and make sure you are on the `micro-wrapper` branch.

#. Type `source doc/usage/octave_setup.sh`

#. Once the script completes without error, type `exit` and close the tab.

#. Back on the Buddy interface, click on "Interactive Apps" and then on "Jupyter".

#. Make sure "Jupyter Session Type" says "Jupyter Lab".

#. Check the "View advanced options" box.

#. Under "Jupyter Version", select "Custom (User supplied Jupyter)".

#. In the "Commands" box, make sure the following lines are typed in. 
   (Be sure to scroll down to see all commands)

   - `module load Anaconda3`

   - `module load Octave`

   - `source ~/.bashrc`

   - `conda activate lysis`

#. Click the "Launch" button.

Running Octave Notebooks
++++++++++++++++++++++++

#. Log into https://buddy.uco.edu

#. Click on "Interactive Apps" and select "Jupyter".

#. All of the settings should remain from last time, 
   otherwise run through the "Initial Setup" steps again.

#. Click "Launch".

#. This will take you to a new page with a list of jobs.
   Your latest job will be at the top and will take a minute or so
   to start up. Once it starts, a "Connect to Jupyter" button will appear.

#. Click on the "Connect to Jupyter" button.

#. In the list of folders on the left, navigate to your `lysis` folder,
   then go to `src\matlab` and open the `micro_to_macro.ipynb` file.

#. Set the `run_code` that you want to analyze and then hit the ">>" 
   button at the top of the tab. This will restart the kernel and then
   run all code in the workbook. It will prompt for confirmation.

