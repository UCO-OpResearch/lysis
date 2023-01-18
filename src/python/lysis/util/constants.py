from enum import Enum, IntEnum, unique

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2022, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"

default_filenames = {
    "unbinding_time": "tsectPA.dat",  # Fortran: tsec1
    # 'leaving_time': "tPAleave.dat",  # Fortran: CDFtPA
    "lysis_time": "lysismat.dat",  # Fortran: lysismat
    "total_lyses": "lenlysisvect.dat",  # Fortran: lenlysismat
    "degradation_state": "deg.dat",  # Fortran: degnext
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
