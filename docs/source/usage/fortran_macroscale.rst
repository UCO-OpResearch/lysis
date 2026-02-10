-----------------------
Macroscale Model
-----------------------
Fortran Script Usage
---------------------

First-time setup
+++++++++++++++++

#. Open a web browser, go to https://ondemand.hpc.uco.edu, and log in.

#. Click on "Files" and then "Home Directory"

#. If you want your lysis code and data in a separate folder/directory,
   create that now using the "New Directory" button.

#. Click on your new folder (if you created one) and then 
   click "Open in Terminal".

#. If you've not yet reset your password, type "passwd" and change it
   
   *NOTE:* When you type your password, nothing will show on the screen.
   This is normal and for security reasons.

#. Copy and paste the following command, then hit enter.
   
   ``git clone https://github.com/UCO-OpResearch/lysis.git``
    
#. Type ``cd lysis`` and hit enter, 
   then run the command ``git checkout macro-wrapper``.

#. When the ``git`` command finishes without error,
   type ``exit`` and hit enter.

#. Return to the tab with your files on it

#. Click on the lysis folder that should be there now 
   (if not, click "Refresh")

#. Hit the "Copy Path" button in the top right.

   You will need this later, so paste it somewhere temporarily for now
   if you don't trust yourself to keep it.

#. Create two new folders called "bin" and "data"

#. At the very top of the page, click "Jobs", then "Job Composer".

#. Click the button that says, "+New Job", and then "From Default Template"

#. Click on the "Job Options" button, set the name of the job to 
   "Macroscale Template", and click "Save".

#. On the bottom-right, is a window titled "main_job.sh". 
   Click on the "Open Editor" button.

#. Copy and paste the text from the ``macro_fortran_run.sh`` file
   (found in the ``scripts/`` directory of the repository) into the editor window.

#. In line 22 (that starts ``LYSIS_ROOT=``), paste the path of the lysis
   folder that you copied earier. Make sure there are no spaces on that line.
   If your path has spaces, put double-quotes around the path.

#. Click the "Save" button in the top left and close the tab.

Running a job
+++++++++++++++++

#. Log into Buddy OnDemand (https://ondemand.hpc.uco.edu) and go to the Job Composer 
   ("Jobs" -> "Job Composer").

#. In the list of jobs, select your "Macroscale Template" job, 
   then click "+New Job" and "From Selected Job".

#. Click on "Job Options", set the name of the job to the Run Code,
   then click "Save".

#. In the "main_job.sh" window, click "Open Editor".

#. In line 23, (that starts ``MICRO_FILE_CODE=``) type the Run Code for the input microscale model 
   with no spaces.

#. In line 24, (that starts ``MACRO_RUN_CODE=``) type the Run Code for the macroscale model with no 
   spaces.

#. In lines 27-28, (that start ``N=`` and ``F``) type the dimensions (in nodes) of the edge grid 
   for the macroscale model with no spaces.

#. Do NOT enter the value of ``frac_forced`` manually in line 51. 
   This will be calculated by the micro_to_macro.py script and passed automatically.

#. Add any parameters that you want to be different from the defaults,
   between line 37 (that starts ``--frac_forced``) and the line that starts
   ``> data``. These MUST have the following format:

   - Start with ``--``, immediately followed by the name of the parameter 
     (see below).

   - Add a space after the name of the parameter, then type the value of
      the parameter without units.

   - Finally, place a backslash (``\``) at the end of each line.

   - There must *NOT* be any blank lines between parameters, 
      or between the parameters and the line that starts ``> data``.

#. Click "Save", then close the tab to return to the Job Composer.

#. Click "> Submit"

#. You can immediately start work on another job.

#. Once the status of the job changes to "Completed" or "Failed",
   check the ``macro########.out`` file in the "Job Details" window
   to make sure there are no errors.

#. You can find the output data of the macroscale code in the data folder
   in a folder named with the ``MACRO_RUN_CODE``. 
   There will be a separate numbered folder for each simulation and its data.

Parameters
+++++++++++++++++

Physical Parameters
#####################################

:radius:
   
   :Description: The radius of each fiber in the model.

   :Default Value: 72.7/2 nanometers

   :Units: microns

:delx:

   :Description: Pore size (distance between fibers/nodes)

   :Default Value: 1.0135 microns

   :Units: centimeters

:Diff:

   :Description: Diffusion coefficient

   :Default Value: 5.0e-7 cm^2/s

   :Units: square centimeters per second

:avgwait:

   :Description: This is the average time a tPA molecule stays bound to fibrin. 
      For now I'm using 27.8 to be 1/0.036, the value in the absence of PLG.

   :Default Value: 27.8 seconds

   :Units: seconds

:frac_forced:

   :Description: Fraction of times tPA was forced to unbind in microscale model.

   :Default Value: 0.0852

   :Units: None

:bs:

   :Description: Concentration of binding sites.

   :Default Value: 427

   :Units: micromolar

:kon:

   :Description: The binding rate of tPA, :math:`k^\text{on}_\text{tPA}`, to fibrin.

   :Default Value: 0.1 (micromolar*sec)^-1

   :Units: per micromolar per sec


Model Parameters
#####################################

:N:

   :Description: The number of lattice nodes in each (horizontal) row

   :Default Value: 93

   :Units: None

:F:

   :Description: The number of lattice nodes in each (vertical) column.

   :Default Value: 121

   :Units: None

:Ffree:

   :Description: The 1st node in vertical direction containing fibers.
      So if Ffree = 10, then rows 0-9 have no fibers, there's one more 
      row of fiber-free planar vertical edges, and then the row with index 
      'Ffree' (e.g. 11th) is a full row of fibers.

   :Default Value: 29

   :Units: None

:M:

   :Description: The total number of tPA molecules

   :Default Value: 43074

   :Units: None

:q:

   :Description: The probability of moving.

   :Default Value: 0.2

   :Units: None


Experimental Parameters
#####################################

:simulations:

   :Description: The number of independent trials run in the macroscale model.

   :Default Value: 10
   
   :Units: None

:tf:

   :Description: Total running time for model.

   :Default Value: 20 minutes
   
   :Units: seconds

:nummicro:

   :Description: The number of independent trials run in the microscale model.

   :Default Value: 500
   
   :Units: hundreds of simulations

:seed:

   :Description: Seed for the random number generator

   :Default Value: 0 (randomly drawn)
   
   :Units: None

:save_interval:

   :Description: How often to record data from the model.

   :Default Value: 10 seconds
   
   :Units: seconds
