from dataclasses import asdict, dataclass, field

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
    location_i: int = None
    location_j: int = None
    bound: bool = False
    leaving_time: float = float('inf')
    waiting_time: float = 0
    binding_time: float = float('inf')
    unbound_by_degradation: bool = False

