from enum import Enum, unique


class Const(dict):
    """Encapsulates a dictionary so that its elements can be accessed as properties.

    e.g.,
        >>> my_const = Const()
        >>> my_const['test'] = 'hat'
        >>> my_const.test
        'hat'
    """

    def __init__(self):
        """Creates a blank Const object."""
        # Call the super-constructor
        super(Const, self).__init__()
        # This is where the wizardry happens.
        # We define the classes parameters as the dictionary itself
        self.__dict__ = self


@unique
class BoundaryDirection(Enum):
    TOP = 0
    BOTTOM = 1
    LEFT = 2
    RIGHT = 3


@unique
class BoundaryConditions(Enum):
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


def get_constants() -> Const:
    CONST = Const()
    CONST.BOUND = BoundaryDirection
    CONST.BOUND_COND = BoundaryConditions
    CONST.DIR = FiberDirection
    return CONST

