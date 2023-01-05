from dataclasses import dataclass

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
    index: int = None
    location_i: int = None
    location_j: int = None
    bound: bool = False
    leaving_time: float = float('inf')
    waiting_time: float = 0
    binding_time: float = float('inf')
    unbound_by_degradation: float = 0
    time_to_reach_back_row: float = float('inf')

