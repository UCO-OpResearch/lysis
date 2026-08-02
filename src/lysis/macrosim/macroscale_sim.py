"""The thin orchestrator that macroscale.py's MacroscaleSim will shrink into.

STATUS: placeholder until strategies/services exist.

Target shape, once filled in:

    class MacroscaleSim:
        def __init__(self, run, draws, move_strategy, bind_strategy, unbind_strategy, recorder):
            # every dependency is passed in, already constructed by factory.py.
            # MacroscaleSim itself never branches on duplicate_fortran.
            ...

        def step(self, ...):
            # order: bind -> unbind -> move -> record
            ...

Keeping this thin is the point: once everything below is proven, this file
should contain almost no logic of its own, just sequencing.
"""
