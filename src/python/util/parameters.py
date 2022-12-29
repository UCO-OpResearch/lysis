from datetime import datetime
import os
import json


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

    def initialize_macro_param(self, params=None):
        self.macro_params = MacroParameters(params)

    def to_dict(self):
        output = {'experiment_code': self.experiment_code,
                  'micro_params': None,
                  'macro_params': None
                  }
        if self.macro_params is not None:
            output["macro_params"] = self.macro_params.to_dict()
        if self.micro_params is not None:
            output["micro_params"] = self.micro_params.to_dict()
        return output
    
    def to_file(self):
        with open(self.os_param_file, 'w') as file:
            json.dump(self.to_dict(), file)

    def read_file(self):
        if os.path.isfile(self.os_param_file):
            with open(self.os_param_file, 'r') as file:
                params = json.load(file)
                macro_params = params.pop('macro_params', None)
                if macro_params is not None:
                    self.macro_params = MacroParameters(macro_params)
                else:
                    self.macro_params = None
        else:
            raise RuntimeError('Experiment parameter file not found.')


# class MicroParameters(object):


class MacroParameters(object):
    defaults = {
        'macro_version':           'diffuse_into_and_along',
        #####################################
        # Physical Parameters
        #####################################
        # The tPA binding rate. Units of inverse (micromolar*sec)
        # (kon)
        'binding_rate':            1.0e-2,
        # Pore size (distance between fibers/nodes), measured in centimeters
        # (delx)
        'pore_size':               1.0135e-4,
        # Diffusion coefficient, measured in cm^2/s
        # (Diff)
        'diffusion_coeff':         5.0e-7,
        # Concentration of binding sites in micromolar
        # (bs)
        'binding_sites':           4.27e+2,
        # Distance from the start of one fiber to the next, in microns
        # because distance between nodes is 1.0135 micron
        # and diameter of 1 fiber is 0.0727 micron
        # (dist)
        'grid_symmetry_distance':  1.0862,

        #####################################
        # Model Parameters
        #####################################
        # The number of lattice nodes in each (horizontal) row
        # (N)
        'rows':                    93,
        # The number of lattice nodes in each (vertical) column
        # (F)
        'cols':                    121,
        # 1st node in vertical direction containing fibers.
        # So if firstFiberRow = 10, then rows 1-9 have no fibers,
        # there's one more row of fiber-free planar vertical edges,
        # and then the row starting with the firstFiberRow-th (e.g. 10th) vertical node
        # is a full row of fibers
        # (Ffree-1)
        'first_fiber_row':         29 - 1,
        # The total number of tPA molecules:
        #      43074 is Colin's [tPA]=0.6 nM
        #      86148 is Colin's [tPA]=1.2 nM
        # (M)
        'total_molecules':         43074,
        # The probability of moving.
        # Make sure it is small enough that we've converged.
        # (q)
        'moving_probability':      0.2,

        #####################################
        # Experimental Parameters
        #####################################
        # The number of independent trials to be run
        # (stats)
        'total_trials':            10,
        # Total running time for model in seconds
        # (tf)
        'total_time':              10 * 60,

        # Seed for the random number generator
        # (seed)
        'seed':                    -2137354075,
        # State for the random number generator
        # (state)
        'state':                   [129281, 362436069, 123456789, None],
        #####################################
        # Data Parameters
        #####################################
        # Data file names
        'unbinding_time_filename': "tsectPA.dat",
        'lysis_time_filename':     "lysismat.dat",
        'total_lyses_filename':    "lenlysisvect.dat",
        # Data size parameters
        'lysis_blocks':            100,
        'unbinds_per_block':       500,
        'max_lyses_per_block':     283,
    }

    independent_parameters = defaults.keys()

    dependent_parameters = [
        # Fibers in a full row of nodes
        'full_row',
        # Number of all x- and z-fibers in a row
        'xz_row',
        # The total number of fibers in the model
        # (num)
        'total_fibers',
        # The last edge number without fibrin
        # (enoFB-1)
        'last_ghost_fiber',
        # The length of one timestep, in seconds
        # (tstep)
        'time_step',
        # The total number of timesteps
        'total_time_steps',
        ]

    def __init__(self, params=None):
        self.params = {}
        self.set_default_parameters()

        if params is not None:
            for k, v in params.items():
                if k in self.independent_parameters:
                    self[k] = v

        self.__calculate_dependent_parameters()

    def __getitem__(self, item):
        return self.params[item]

    # In the future, this should be rewritten so that values cannot be changed after initialization
    def __setitem__(self, key, value):
        if key in self.independent_parameters:
            self.params[key] = value
            self.__calculate_dependent_parameters()
        elif key in self.dependent_parameters:
            raise RuntimeError('Cannot set the value of a dependent parameter directly.')
        else:
            raise KeyError('Parameter not found.')

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

    def set_default_parameters(self):
        self.params = self.defaults.copy()

    def __calculate_dependent_parameters(self):
        self.params['full_row'] = 3 * self['rows'] - 1
        self.params['xz_row'] = 2 * self['rows'] - 1
        self.params['total_fibers'] = ((2 * self['rows'] - 1) * self['cols']
                             + self['rows'] * (self['cols'] - 1))
        self.params['last_ghost_fiber'] = ((3 * self['cols'] - 1)
                                 * self['first_fiber_row'] - 1)
        self.params['time_step'] = (self['moving_probability']
                          * self['pore_size'] ** 2
                          / (12 * self['diffusion_coeff']))
        self.params['total_time_steps'] = self['total_time'] / self['time_step']
        if self['state'][3] is None:
            self['state'][3] = self['seed']

    def to_dict(self):
        return self.params

    
