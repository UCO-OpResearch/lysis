"""Modular replacement for the monolithic macroscale.py.

Package layout (see simulation/README.md for the full contributor guide):

    macrosim/
        macroscale_sim.py   -- thin orchestrator
        factory.py          -- composition root; the ONLY place duplicate_fortran
                                is read to decide which concrete classes get built
        strategies/
            base.py          -- MoveStrategy / BindStrategy / UnbindStrategy ABCs
            move.py          -- NativeMoveStrategy / FortranMoveStrategy
            bind.py          -- Native/Fortran bind strategies
            unbind.py        -- UnbindDegradation / UnbindTime variants
        services/
            random_draw.py   -- RandomDrawSource protocol + Native/Fortran impls
            recording.py     -- Recorder service

Nothing here is wired up yet -- this file exists so the target shape is
reviewable before any logic moves out of macroscale.py.
"""
