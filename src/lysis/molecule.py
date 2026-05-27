from dataclasses import dataclass, field


@dataclass
class Molecule:
    index: int = field(default=None, compare=False)
    location: int = field(default=None, compare=False)
    bound: bool = field(default=False, compare=False)
    leaving_time: float = field(default=float("inf"), compare=False)
    waiting_time: float = field(default=0, compare=False)
    binding_time: float = field(default=float("inf"), compare=True)
    unbound_by_degradation: float = field(default=0, compare=False)
    time_to_reach_back_row: float = field(default=float("inf"), compare=False)
