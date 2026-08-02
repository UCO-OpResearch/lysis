"""Composition root -- the ONLY module allowed to branch on duplicate_fortran.

STATUS: placeholder until strategies/services exist.

This is the enforced coupling point: every strategy category and service 
must get exactly one entry here per mode. The structural
test in test_factory_registry.py will walk each ABC's
registered subclasses and fail if one isn't wired up here -- so a
half-finished new strategy can't silently slip through.

Target shape, once strategies/services exist:

    def make_draw_source(duplicate_fortran: bool, seed: int): ...
    def make_move_strategy(duplicate_fortran: bool, draws, moving_probability): ...
    def make_bind_strategy(duplicate_fortran: bool, draws, ...): ...
    def make_unbind_strategy(duplicate_fortran: bool, draws, ...): ...
    def make_recorder(run): ...

    def build_simulation(run) -> "MacroscaleSim":
        # the single place that reads run.macro_params.duplicate_fortran
        # and assembles a fully-wired MacroscaleSim.
        ...
"""
