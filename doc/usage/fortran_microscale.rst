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

#. Copy and paste the text from the `micro_fortran_run.sh` file 
   into the editor window

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

   - Finally, place a `\` at the end of each line.

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
   
   :Description: The dissociation constant of tPA, :math:`k^D_\\text{tPA}`, to fibrin 
      in the presence of PLG.

   :Default Value: 0.02 micromolar

   :Units: micromolar

:KdtPAnoplg:

diss_const_tPA_woPLG: Quantity = Q_("")
   :Description: The dissociation constant of tPA, :math:`k^D_\\text{tPA}`, to fibrin
      in the absence of PLG.

   :Default Value: 0.36 micromolar

   :Units: micromolar


:KdPLGintact:
diss_const_PLG_intact: Quantity = Q_("38 micromolar")
"""The dissociation constant of PLG, :math:`k^D_\\text{PLG}`, to intact fibrin.

:Units: micromolar
:Fortran: """

diss_const_PLG_nicked: Quantity = Q_("2.2 micromolar")
"""The dissociation constant of PLG, :math:`k^D_\\text{PLG}`, to nicked fibrin.

:Units: micromolar
:Fortran: KdPLGnicked"""

bind_rate_tPA: Quantity = Q_("0.1 (micromolar*sec)^-1")
"""The binding rate of tPA, :math:`k^\\text{on}_\\text{tPA}`, to fibrin.

:Units: (micromolar*sec)^-1
:Fortran: ktPAon"""

bind_rate_PLG: Quantity = Q_("0.1 (micromolar*sec)^-1")
"""The binding rate of PLG, :math:`k^\\text{on}_\\text{PLG}`, to fibrin.

:Units: (micromolar*sec)^-1
:Fortran: kPLGon"""

conc_free_PLG: Quantity = Q_("2 micromolar")
"""The concentration of free plasminogen.

:Units: micromolar
:Fortran: freeplg"""

deg_rate_fibrin: Quantity = Q_("5 sec^-1")
"""The plasmin-mediated rate of fibrin degradation.

:Units: sec^-1
:Fortran: kdeg"""

unbind_rate_PLG_intact: Quantity = field(init=False)
"""The unbinding rate of PLG, :math:`k^\\text{off}_\\text{PLG}`, 
from intact fibrin.

:Units: sec^-1
:Fortran: kplgoff"""

unbind_rate_PLG_nicked: Quantity = field(init=False)
"""The unbinding rate of PLG, :math:`k^\\text{off}_\\text{PLG}`, 
from nicked fibrin.

:Units: sec^-1
:Fortran: kplgoffnick"""

unbind_rate_PLi: Quantity = Q_("57.6 sec^-1")
"""The unbinding rate of PLi, :math:`k^\\text{off}_\\text{PLi}`, 
from fibrin.

:Units: sec^-1
:Fortran: kplioff"""

unbind_rate_tPA_wPLG: Quantity = field(init=False)
"""The unbinding rate of tPA, :math:`k^\\text{off}_\\text{tPA}`, 
from fibrin in the presence of PLG.

:Units: sec^-1
:Fortran: kaoff12"""

unbind_rate_tPA_woPLG: Quantity = field(init=False)
"""The unbinding rate of tPA, :math:`k^\\text{off}_\\text{tPA}`, 
from fibrin in the absence of PLG.

:Units: sec^-1
:Fortran: kaoff10"""

activation_rate_PLG: Quantity = Q_("0.1 sec^-1")
"""The catalytic rate constant, :math:`k_\\text{cat}^\\text{ap}`, 
for activation of PLG into PLI.

:Units: sec^-1
:Fortran: kapcat"""

exposure_rate_binding_site: Quantity = Q_("5 sec^-1")
"""The catalytic rate constant, :math:`k_\\text{cat}^\\text{n}`, 
for the PLi-mediated rate of exposure of new binding sites.

:Units: sec^-1
:Fortran: kncat"""

protein_per_fiber: Quantity = field(init=False)
"""The fraction of protein in each fiber (by volume?)

:Units: %
:Fortran: None"""

fibrin_conc_per_fiber: Quantity = field(init=False)
"""The concentration of fibrin in each fiber

:Units: micromolar
:Fortran: None"""

binding_sites: Quantity = field(init=False)  # int = 427
"""Concentration of binding sites.

:Units: micromolar
:Fortran: bs"""

#####################################
# Model Parameters
#####################################

nodes_in_row: int = 7
"""The number of protofibrils in one row of the lattice inside one
fiber.

:Units: None
:Fortran: nodes"""

snap_proportion: float = 2.0/3.0
"""The proportion of doublets that need to be degraded before the
fiber snaps.

:Units: None
:Fortran: snap_proportion"""

#####################################
# Experimental Parameters
#####################################

simulations: int = 50_000
"""The number of independent trials run in the microscale model.

:Units: None
:Fortran: runs"""

seed: int = 0
"""Seed for the random number generator

:Units: None
:Fortran: seed"""