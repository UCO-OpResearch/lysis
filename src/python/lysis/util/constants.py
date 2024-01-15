from enum import Enum, IntEnum, unique, Flag

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2023, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__version__ = "0.2"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"

default_filenames = {
    ## Microscale Out
    "lysis_complete_time": "lysis",  # Fortran: lysis_time
    "tPA_leaving_time": "tPA_time",  # Fortran: tPA_time
    "PLi_generated": "PLi",  # Fortran: Plasmin
    "lysis_completed": "lyscomplete",  # Fortran: lysiscomplete
    "tPA_kinetic_unbound": "tPAunbind",  # Fortran: tPAunbind
    "tPA_forced_unbound": "tPAPLiunbind",  # Fortran: tPAPLiunbd
    "tPA_still_bound": "lasttPA",  # Fortran: ltPA
    "first_PLi": "firstPLi",  # Fortran: firstPLi
    ## Macroscale In
    "unbinding_time_dist": "tsectPA",  # Fortran: tsec1
    # 'leaving_time': "tPAleave",  # Fortran: CDFtPA
    "lysis_time_dist": "lysismat",  # Fortran: lysismat
    "total_lyses": "lenlysisvect",  # Fortran: lenlysismat
    ## Macroscale Out
    "degradation_state": "deg",  # Fortran: degnext
    "molecule_location": "m_loc",
    "molecule_state": "m_bound",
    "save_time": "tsave",  # Fortran: tsave
}


class Const:
    def __init__(self):
        self.BOUND = BoundaryDirection
        self.BOUND_COND = BoundaryCondition
        self.DIR = FiberDirection
        self.NEIGHBORHOOD = Neighbors()


class Neighbors:
    def __init__(self):
        self.X = ((-1, -1, 0, 0, 0, 0, 0, 0), (-2, 1, -2, 1, -1, -1, 2, 2))
        self.Y = ((0, 0, 1, 1, 0, 1, 0, 1), (1, 1, 1, 1, -1, -1, 2, 2))
        self.Z = ((-1, -1, 0, 0, 0, 0, 0, 0), (-1, -1, -1, -1, -2, -2, 1, 1))
        self.TOP_REFL = (0, 0, -1, -1, 0, 0, 0, 0)
        self.BOTTOM_REFL = (1, 1, 0, 0, 0, 0, 0, 0)
        self.LEFT_REFL = (0, 0, 0, 0, 3, 3, 0, 0)
        self.RIGHT_REFL = (0, 0, 0, 0, 0, 0, -3, -3)


@unique
class RandomDraw(IntEnum):
    BINDING_TIME_WHEN_UNBINDING = 0
    BINDING_TIME_WHEN_MOVING = 1
    MICRO_UNBIND = 2
    MOVE = 3
    UNBINDING_TIME = 4
    LYSIS_TIME = 5
    CONFLICT_RESOLUTION = 6
    RESTRICTED_MOVE = 7


@unique
class BoundaryDirection(IntEnum):
    TOP = 0
    BOTTOM = 1
    LEFT = 2
    RIGHT = 3
    FRONT = 4
    BACK = 5


@unique
class BoundaryCondition(Enum):
    REFLECTING = 0
    PERIODIC = 1
    CONTINUING = 2


@unique
class FiberDirection(Enum):
    UP = 1
    DOWN = -1
    LEFT = 2
    RIGHT = -2
    OUT = 3
    IN = -3


@unique
class ExpComponent(Flag):
    NONE = 0
    MICRO = 1
    MICRO_POSTPROCESSING = 2
    MACRO = 4
    MACRO_POSTPROCESSING = 8
    ALL = MICRO | MICRO_POSTPROCESSING | MACRO | MACRO_POSTPROCESSING
