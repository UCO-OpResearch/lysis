-----------------------
Microscale Model
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
   then run the command ``git checkout micro-wrapper``.

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
   "Microscale Template", and click "Save".

#. On the bottom-right, is a window titled "main_job.sh". 
   Click on the "Open Editor" button.

#. Copy and paste the text from the ``micro_fortran_run.sh`` file
   (found in the ``scripts/`` directory of the repository) into the editor window.

#. In line 18 (that starts ``LYSIS_ROOT=``), paste the path of the lysis
   folder that you copied earier. Make sure there are no spaces on that line.
   If your path has spaces, put double-quotes around the path.

#. Click the "Save" button in the top left and close the tab.

Running a job
+++++++++++++++++

#. Log into Buddy OnDemand (https://ondemand.hpc.uco.edu) and go to the Job Composer 
   ("Jobs" -> "Job Composer").

#. In the list of jobs, select your "Microscale Template" job, 
   then click "+New Job" and "From Selected Job".

#. Click on "Job Options", set the name of the job to the Run Code,
   then click "Save".

#. In the "main_job.sh" window, click "Open Editor".

#. In line 19, (that starts ``RUN_CODE=``) type the Run Code with no spaces.

#. Add any parameters that you want to be different from the defaults,
   between line 30 (that starts ``--outFileCode``) and the line that starts
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
   check the ``micro_rates_########.out`` file in the "Job Details" window
   to make sure there are no errors.

#. You can find the output data of the microscale code in the data folder
   in a folder named with the ``RUN_CODE``.

Parameters
+++++++++++++++++

Physical Parameters
#####################################

:radius:

   :Description: The radius of each fiber in the model.

   :Default Value: 72.7/2 nanometers

   :Units: microns

   :Python Name: ``fiber_radius``

:KdtPAyesplg:

   :Description: The dissociation constant of tPA, :math:`k^D_\text{tPA}`, to fibrin
      in the presence of PLG.

   :Default Value: 0.02 micromolar

   :Units: micromolar

   :Python Name: ``diss_const_tPA_wPLG``

:KdtPAnoplg:

   :Description: The dissociation constant of tPA, :math:`k^D_\text{tPA}`, to fibrin
      in the absence of PLG.

   :Default Value: 0.36 micromolar

   :Units: micromolar

   :Python Name: ``diss_const_tPA_woPLG``


:KdPLGintact:

   :Description: The dissociation constant of PLG, :math:`k^D_\text{PLG}`, to intact fibrin.

   :Default Value: 38 micromolar

   :Units: micromolar

   :Python Name: ``diss_const_PLG_intact``

:KdPLGnicked:

   :Description: The dissociation constant of PLG, :math:`k^D_\text{PLG}`, to nicked fibrin.

   :Default Value: 2.2 micromolar

   :Units: micromolar

   :Python Name: ``diss_const_PLG_nicked``

:ktPAon:

   :Description: The binding rate of tPA, :math:`k^\text{on}_\text{tPA}`, to fibrin.

   :Default Value: 0.1 (micromolar*sec)^-1

   :Units: (micromolar*sec)^-1

   :Python Name: ``bind_rate_tPA``

:kplgon:

   :Description: The binding rate of PLG, :math:`k^\text{on}_\text{PLG}`, to fibrin.

   :Default Value: 0.1 (micromolar*sec)^-1

   :Units: (micromolar*sec)^-1

   :Python Name: ``bind_rate_PLG``

:freeplg:

   :Description: The concentration of free plasminogen.

   :Default Value: 2 micromolar

   :Units: micromolar

   :Python Name: ``conc_free_PLG``

:kdeg:

   :Description: The plasmin-mediated rate of fibrin degradation.

   :Default Value: 5 sec^-1

   :Units: sec^-1

   :Python Name: ``deg_rate_fibrin``


:kplioff:

   :Description: The unbinding rate of PLi, :math:`k^\text{off}_\text{PLi}`,
      from fibrin.

   :Default Value: 57.6 sec^-1

   :Units: sec^-1

   :Python Name: ``unbind_rate_PLi``

:kapcat:

   :Description: The catalytic rate constant, :math:`k_\text{cat}^\text{ap}`,
      for activation of PLG into PLI.

   :Default Value: 0.1 sec^-1

   :Units: sec^-1

   :Python Name: ``activation_rate_PLG``

:kncat:

   :Description: The catalytic rate constant, :math:`k_\text{cat}^\text{n}`,
      for the PLi-mediated rate of exposure of new binding sites.

   :Default Value: 5 sec^-1

   :Units: sec^-1

   :Python Name: ``exposure_rate_binding_site``



Model Parameters
#####################################

:nodes:

   :Description: The number of protofibrils in one row of the lattice inside one
      fiber.

   :Default Value: 7

   :Units: None

   :Python Name: ``nodes_in_micro_row``

:snap_proportion:

   :Description: The proportion of doublets that need to be degraded before the
      fiber snaps.

   :Default Value: 0.6666666666667

   :Units: None

   :Python Name: ``snap_proportion``

Experimental Parameters
#####################################

:simulations:

   :Description: The number of independent trials run in the microscale model.

   :Default Value: 50,000

   :Units: None

   :Python Name: ``micro_simulations``

:seed:

   :Description: Seed for the random number generator. Stored Python-side as
      :class:`numpy.uint32`; cast internally to a signed ``INTEGER*4`` for the
      Fortran CLI and reinterpreted as ``uint_least32_t`` by the C KISS RNG,
      so a ``np.uint32`` value round-trips bit-exactly.

   :Default Value: 0 (randomly drawn)

   :Units: None

   :Python Name: ``micro_seed`` (Fortran ``INTEGER*4`` bits reinterpreted as ``np.uint32``)
