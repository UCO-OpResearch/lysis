-----------------------
Microscale Model
-----------------------
Fortran Script Usage
---------------------

First-time setup
+++++++++++++++++

#. Open a web browser, go to buddy.uco.edu, and log in.

#. Click on "Files" and then "Home Directory"

#. If you want your lysis code and data in a separate folder/directory,
   create that now using the "New Directory" button.

#. Click on your new folder (if you created one) and then 
   click "Open in Terminal".

#. If you've not yet reset your password, type "passwd" and change it
   
   *NOTE:* When you type your password, nothing will show on the screen.
   This is normal and for security reasons.

#. Copy and paste the following command, then hit enter.
   
   `git clone https://github.com/UCO-OpResearch/lysis.git`
    

#. When the `git clone` command finishes without error, 
   type `exit` and hit enter.

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
   "Microscale Template", and click "Save".

#. On the bottom-right, is a window titled "main_job.sh". 
   Click on the "Open Editor" button.

#. Copy and paste the text from the `micro_fortran_run.sh` file (link below)
   into the editor window.

   https://github.com/UCO-OpResearch/lysis/blob/micro-wrapper/doc/usage/micro_fortran_run.sh

#. In line 18 (that starts `LYSIS_ROOT=`), paste the path of the lysis
   folder that you copied earier. Make sure there are no spaces that line.
   If your path has spaces, put double-quotes around it.

#. Click the "Save" button in the top right and close the tab.

Running a job
+++++++++++++++++

#. Log into Buddy (https://buddy.uco.edu) and go to the Job Composer 
   ("Jobs" -> "Job Composer").

#. In the list of jobs, select your "Microscale Template" job, 
   then click "+New Job" and "From Selected Job".

#. Click on "Job Options", set the name of the job to the Run Code,
   then click "Save".

#. In the "main_job.sh" window, click "Open Editor".

#. In line 19, (that starts `RUN_CODE=`) type the Run Code with no spaces.

#. Add any parameters that you want to be different from the defaults,
   between line 30 (that starts `--outFileCode`) and the line that starts
   `> data`. These MUST have the following format:

   - Start with `--`, immediately followed by the name of the parameter 
     (see below).

   - Add a space after the name of the parameter, then type the value of
      the parameter without units.

   - Finally, place a `\\` at the end of each line.

   - There must *NOT* be any blank lines between parameters, 
      or between the parameters and the line that starts `> data`.

#. Click "Save", then close the tab to return to the Job Composer.

#. Click "> Submit"

#. You can immediately start work on another job.

#. Once the status of the job changes to "Completed" or "Failed",
   check the `micro_rates_########.out` file in the "Job Details" window
   to make sure there are no errors.

#. You can find the output data of the microscale code in the data folder
   in a folder named with the `RUN_CODE`.

Parameters
+++++++++++++++++

Physical Parameters
#####################################

:radius:
   
   :Description: The radius of each fiber in the model.

   :Default Value: 72.7/2 nanometers

   :Units: microns

:KdtPAyesplg:
   
   :Description: The dissociation constant of tPA, :math:`k^D_\text{tPA}`, to fibrin 
      in the presence of PLG.

   :Default Value: 0.02 micromolar

   :Units: micromolar

:KdtPAnoplg:

   :Description: The dissociation constant of tPA, :math:`k^D_\text{tPA}`, to fibrin
      in the absence of PLG.

   :Default Value: 0.36 micromolar

   :Units: micromolar


:KdPLGintact:

   :Description: The dissociation constant of PLG, :math:`k^D_\text{PLG}`, to intact fibrin.

   :Default Value: 38 micromolar

   :Units: micromolar

:KdPLGnicked:

   :Description: The dissociation constant of PLG, :math:`k^D_\text{PLG}`, to nicked fibrin.

   :Default Value: 2.2 micromolar

   :Units: micromolar

:ktPAon:

   :Description: The binding rate of tPA, :math:`k^\text{on}_\text{tPA}`, to fibrin.
   
   :Default Value: 0.1 (micromolar*sec)^-1

   :Units: (micromolar*sec)^-1

:kPLGon:

   :Description: The binding rate of PLG, :math:`k^\text{on}_\text{PLG}`, to fibrin.

   :Default Value: 0.1 (micromolar*sec)^-1
   
   :Units: (micromolar*sec)^-1

:freeplg:

   :Description: The concentration of free plasminogen.

   :Default Value: 2 micromolar
   
   :Units: micromolar

:kdeg:

   :Description: The plasmin-mediated rate of fibrin degradation.

   :Default Value: 5 sec^-1
   
   :Units: sec^-1


:kplioff:

   :Description: The unbinding rate of PLi, :math:`k^\text{off}_\text{PLi}`, 
      from fibrin.

   :Default Value: 57.6 sec^-1
   
   :Units: sec^-1

:kapcat:

   :Description: The catalytic rate constant, :math:`k_\text{cat}^\text{ap}`, 
      for activation of PLG into PLI.

   :Default Value: 0.1 sec^-1
   
   :Units: sec^-1

:kncat:

   :Description: The catalytic rate constant, :math:`k_\text{cat}^\text{n}`, 
      for the PLi-mediated rate of exposure of new binding sites.

   :Default Value: 5 sec^-1
   
   :Units: sec^-1



Model Parameters
#####################################

:nodes:

   :Description: The number of protofibrils in one row of the lattice inside one
      fiber.

   :Default Value: 7
   
   :Units: None

:snap_proportion:

   :Description: The proportion of doublets that need to be degraded before the
      fiber snaps.

   :Default Value: 0.6666666666667
   
   :Units: None

Experimental Parameters
#####################################

:runs:

   :Description: The number of independent trials run in the microscale model.

   :Default Value: 50_000
   
   :Units: None

:seed:

   :Description: Seed for the random number generator

   :Default Value: 0 (randomly drawn)
   
   :Units: None
