"""Modular replacement for the monolithic macroscale.py.

Package layout (see macrosim/README.md for the contributor guide):

    macrosim/
        macroscale_sim.py       -- MacroscaleSim: owns simulation state and
                                    the per-timestep loop, delegating every
                                    mode-dependent decision to its injected
                                    strategies.
        factory.py               -- composition root; the only place
                                    duplicate_fortran is read, via a
                                    native/fortran registry per category.
        strategies/
            base.py              -- MoveStrategy / BindStrategy /
                                    UnbindStrategy / ConflictResolutionStrategy
                                    ABCs
            move.py               -- NativeMoveStrategy / FortranMoveStrategy
            bind.py               -- DefaultBindStrategy (single
                                    implementation -- bind() has no actual
                                    Native/Fortran divergence)
            unbind.py             -- NativeUnbindStrategy / FortranUnbindStrategy
            conflict_resolution.py -- Native/FortranConflictResolutionStrategy
        services/
            random_draw.py        -- RandomDrawSource + Native/FortranDrawSource
            binding_time.py        -- BindingTimeFactory (batch binding-time
                                    generation, shared by both modes)
"""