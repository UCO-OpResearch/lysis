"""Abstract base classes for every strategy category.

STATUS: placeholder.

Each ABC's abstract methods are the real "contract" that keeps this codebase
coupled-on-purpose: a new subclass that doesn't implement one raises at
instantiation time, not three timesteps into a run.

Target shape:

    class MoveStrategy(abc.ABC):
        @abc.abstractmethod
        def should_move(self, total_molecules, bound): ...
        @abc.abstractmethod
        def pick_neighbor(self, free_to_move): ...

    class BindStrategy(abc.ABC):
        ...

    class UnbindStrategy(abc.ABC):
        ...
        
"""
