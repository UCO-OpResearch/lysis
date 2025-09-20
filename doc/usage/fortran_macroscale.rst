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
   
   `git clone https://github.com/UCO-OpResearch/lysis.git`
    
#. Type `cd lysis` and hit enter, 
   then run the command `git checkout macro-wrapper`.

#. When the `git` command finishes without error,
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
   "Macroscale Template", and click "Save".

#. On the bottom-right, is a window titled "main_job.sh". 
   Click on the "Open Editor" button.

#. Copy and paste the text from the `macro_fortran_run.sh` file (link below)
   into the editor window.

   https://github.com/UCO-OpResearch/lysis/blob/macro-wrapper/doc/usage/macro_fortran_run.sh

#. In line 18 (that starts `LYSIS_ROOT=`), paste the path of the lysis
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

#. In line 22, (that starts `MICRO_RUN_CODE=`) type the Run Code for the microscale model with no spaces.

#. In line 23, (that starts `MACRO_RUN_CODE=`) type the Run Code for the macroscale model with no spaces.

#. Add any parameters that you want to be different from the defaults,
   between line 37 (that starts `--outFileCode`) and the line that starts
   `> data`. These MUST have the following format:

   - Start with `--`, immediately followed by the name of the parameter 
     (see below).

   - Add a space after the name of the parameter, then type the value of
      the parameter without units.

   - Finally, place a backslash (``\``) at the end of each line.

   - There must *NOT* be any blank lines between parameters, 
      or between the parameters and the line that starts `> data`.

#. Click "Save", then close the tab to return to the Job Composer.

#. Click "> Submit"

#. You can immediately start work on another job.

#. Once the status of the job changes to "Completed" or "Failed",
   check the `macro_########.out` file in the "Job Details" window
   to make sure there are no errors.

#. You can find the output data of the macroscale code in the data folder
   in a folder named with the `MACRO_RUN_CODE`. 
   There will be a separate numbered folder for each simulation and its data.

Parameters
+++++++++++++++++

Physical Parameters
#####################################

:radius:
   
   :Description: The radius of each fiber in the model.

   :Default Value: 72.7/2 nanometers

   :Units: microns




Model Parameters
#####################################

:microscale_nodes:

   :Description: The number of protofibrils in one row of the lattice inside one
      fiber.

   :Default Value: 7
   
   :Units: None



Experimental Parameters
#####################################

:simulations:

   :Description: The number of independent trials run in the microscale model.

   :Default Value: 50_000
   
   :Units: None

:seed:

   :Description: Seed for the random number generator

   :Default Value: 0 (randomly drawn)
   
   :Units: None


case ('runCode')
case ('inFileCode')
case ('outFileCode')
case ('N')
case ('F')
case ('Ffree')
case ('simulations')
case ('M')
case ('tf')
case ('nummicro')
case ('kon')
case ('frac_forced')
case ('avgwait')
case ('q')
case ('delx')
case ('Diff')
case ('bs')
case ('radius')
case ('seed')
case ('save_interval')
