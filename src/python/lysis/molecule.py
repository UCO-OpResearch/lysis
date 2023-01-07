from dataclasses import dataclass, field

__author__ = "Brittany Bannish and Bradley Paynter"
__copyright__ = "Copyright 2022, Brittany Bannish"
__credits__ = ["Brittany Bannish", "Bradley Paynter"]
__license__ = ""
__version__ = "0.1"
__maintainer__ = "Bradley Paynter"
__email__ = "bpaynter@uco.edu"
__status__ = "Development"


@dataclass
class Molecule:
    index: int = field(default=None, compare=False)
    location_i: int = field(default=None, compare=False)
    location_j: int = field(default=None, compare=False)
    bound: bool = field(default=False, compare=False)
    leaving_time: float = field(default=float('inf'), compare=False)
    waiting_time: float = field(default=0, compare=False)
    binding_time: float = field(default=float('inf'), compare=True)
    unbound_by_degradation: float = field(default=0, compare=False)
    time_to_reach_back_row: float = field(default=float('inf'), compare=False)

