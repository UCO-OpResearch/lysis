from datetime import datetime
import os


class Experiment(object):
    def __init__(self, data_root, experiment_code=-1):
        if not os.path.isdir(data_root):
            raise RuntimeError('Data folder not found.')
        self.os_data_root = data_root
        if experiment_code == -1:
            now = datetime.now()
            self.experiment_code = str.join('', [str(now.year),
                                                 '_',
                                                 str(now.month),
                                                 '_',
                                                 str(now.day),
                                                 '_',
                                                 str(now.hour),
                                                 str(now.day)
                                                 ])
        else:
            self.experiment_code = experiment_code

        self.os_path = os.path.join(data_root, str(self.experiment_code))
        self.os_param_file = os.path.join(self.os_path, 'params.json')
        if os.path.isfile(self.os_param_file):
            raise RuntimeError('Experiment already exists. Use load_experiment() instead.')

        self.micro_params = None
        self.macro_params = None

    def __str__(self):
        values = self.to_dict()
        output = ''
        nl = os.linesep
        key_len = 0
        for k in values.keys():
            key_len = max(len(k), key_len)
        key_len += 1
        values.pop('micro_params', None)
        values.pop('macro_params', None)
        for k, v in values.items():
            output += f'{k:<{key_len}}: {v}' + nl
        tab = ' '*(key_len+2)
        for param in ['micro_params', 'macro_params']:
            output += f'{param:<{key_len}}: '
            if self.__dict__[param] is None:
                output += 'None' + nl
            else:
                param_string = str(self.__dict__[param])
                output += tab.join(param_string.splitlines(True))
        return output

    def initialize_macro_param(self, macro_ver, params=None):
        my_params = MacroParameters(macro_ver)

        #####################################
        # Physical Parameters
        #####################################
        # The tPA binding rate. Units of inverse (micromolar*sec)
        # (kon)
        my_params.binding_rate = 1.0e-2
        # Pore size (distance between fibers/nodes), measured in centimeters
        # (delx)
        my_params.pore_size = 1.0135e-4
        # Diffusion coefficient, measured in cm^2/s
        # (Diff)
        my_params.diffusion_coeff = 5.0e-7
        # Concentration of binding sites in micromolar
        # (bs)
        my_params.binding_sites = 4.27e+2
        # Distance from the start of one fiber to the next, in microns
        # because distance between nodes is 1.0135 micron
        # and diameter of 1 fiber is 0.0727 micron
        # (dist)
        my_params.grid_symmetry_distance = 1.0862

        #####################################
        # Model Parameters
        #####################################
        # The number of lattice nodes in each (horizontal) row
        # (N)
        my_params.rows = 93
        # The number of lattice nodes in each (vertical) column
        # (F)
        my_params.cols = 121
        # 1st node in vertical direction containing fibers.
        # So if firstFiberRow = 10, then rows 1-9 have no fibers,
        # there's one more row of fiber-free planar vertical edges,
        # and then the row starting with the firstFiberRow-th (e.g. 10th) vertical node
        # is a full row of fibers
        # (Ffree-1)
        my_params.first_fiber_row = 29 - 1
        # The total number of tPA molecules:
        #      43074 is Colin's [tPA]=0.6 nM
        #      86148 is Colin's [tPA]=1.2 nM
        # (M)
        my_params.total_molecules = 43074
        # The probability of moving.
        # Make sure it is small enough that we've converged.
        # (q)
        my_params.moving_probability = 0.2

        #####################################
        # Experimental Parameters
        #####################################
        # The number of independent trials to be run
        # (stats)
        my_params.total_trials = 10
        # Total running time for model in seconds
        # (tf)
        my_params.total_time = 10 * 60

        # Seed for the random number generator
        # (seed)
        my_params.seed = -2137354075
        # State for the random number generator
        # (state)
        my_params.state = [129281, 362436069, 123456789, my_params.seed]
        #####################################
        # Data Parameters
        #####################################
        # Data file names
        my_params.unbinding_time_filename = "tsectPA"
        my_params.lysis_time_filename = "lysismat"
        my_params.total_lyses_filename = "lenlysisvect"
        # Data size parameters
        my_params.lysis_blocks = 100
        my_params.unbinds_per_block = 500
        my_params.max_lyses_per_block = 283

        if isinstance(params, dict):
            my_params.__dict__.update(params)

        self.macro_params = my_params

    def to_dict(self):
        output = self.__dict__.copy()
        for k in list(output.keys()):
            if k[:2] == 'os':
                output.pop(k, None)
        output["macro_params"] = self.macro_params.to_dict()
        return output

#class MicroParameters(object):



class MacroParameters(object):
    def __init__(self, macro_ver):
        self.macro_version = macro_ver
        self.__rows = 0
        self.__cols = 0
        self.__first_fiber_row = 0
        self.__moving_probability = 1
        self.__pore_size = 1
        self.__diffusion_coeff = 1
        self.__total_time = 0
        self.__calculate_dependent_parameters()

    def __str__(self):
        values = self.to_dict()
        output = ''
        nl = os.linesep
        key_len = 0
        for k in values.keys():
            key_len = max(len(k), key_len)
        key_len += 1
        for k, v in values.items():
            output += f'{k:<{key_len}}: {v}' + nl
        return output

    @property
    def rows(self):
        return self.__rows

    @rows.setter
    def rows(self, rows):
        self.__rows = rows
        self.__calculate_dependent_parameters()

    @property
    def cols(self):
        return self.__cols

    @cols.setter
    def cols(self, cols):
        self.__cols = cols
        self.__calculate_dependent_parameters()

    @property
    def first_fiber_row(self):
        return self.__first_fiber_row

    @first_fiber_row.setter
    def first_fiber_row(self, first_fiber_row):
        self.__first_fiber_row = first_fiber_row
        self.__calculate_dependent_parameters()

    @property
    def moving_probability(self):
        return self.__moving_probability

    @moving_probability.setter
    def moving_probability(self, moving_probability):
        self.__moving_probability = moving_probability
        self.__calculate_dependent_parameters()

    @property
    def pore_size(self):
        return self.__pore_size

    @pore_size.setter
    def pore_size(self, pore_size):
        self.__pore_size = pore_size
        self.__calculate_dependent_parameters()

    @property
    def diffusion_coeff(self):
        return self.__diffusion_coeff

    @diffusion_coeff.setter
    def diffusion_coeff(self, diffusion_coeff):
        self.__diffusion_coeff = diffusion_coeff
        self.__calculate_dependent_parameters()

    @property
    def total_time(self):
        return self.__total_time

    @total_time.setter
    def total_time(self, total_time):
        self.__total_time = total_time
        self.__calculate_dependent_parameters()

    def __calculate_dependent_parameters(self):
        # Fibers in a full row of nodes
        self.full_row = 3 * self.__rows - 1
        # Number of all x- and z-fibers in a row
        self.xz_row = 2 * self.__rows - 1
        # The total number of fibers in the model
        # (num)
        self.total_fibers = ((2 * self.__rows - 1) * self.__cols
                             + self.__rows * (self.__cols - 1))
        # The last edge number without fibrin
        # (enoFB-1)
        self.last_ghost_fiber = ((3 * self.__cols - 1)
                                 * self.__first_fiber_row - 1)
        # The length of one timestep, in seconds
        # (tstep)
        self.time_step = (self.__moving_probability
                          * self.__pore_size ** 2
                          / (12 * self.__diffusion_coeff))
        # The total number of timesteps
        self.total_time_steps = self.__total_time / self.time_step

    def to_dict(self):
        output = {}
        for k, v in self.__dict__.items():
            if '__' in k:
                key = k.split('__')[1]
            else:
                key = k
            output[key] = v
        return output



