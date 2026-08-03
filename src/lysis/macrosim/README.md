# macrosim package

Modular replacement for the monolithic `macroscale.py`.

## Adding a new strategy variant (e.g. a third Move mode)

1. Subclass the relevant ABC in `strategies/base.py` (e.g. `MoveStrategy`).
   Implement every abstract method -- Python will refuse to instantiate
   your class otherwise.
2. Add an entry to the matching registry dict in `factory.py`
   (e.g. `_MOVE_STRATEGIES`). This is the only file that selects between
   variants; nothing else should branch on mode.
3. If the new mode needs a different selector than `duplicate_fortran`
   (e.g. a mode string instead of a bool), update the one line in
   `build_simulation` that computes `mode` -- the registries and `make_*`
   functions don't need to change.

There's currently no automated check that a new subclass actually gets
registered in step 2. If it matters later, a small test walking each ABC's
`__subclasses__()` against its registry dict would close it.

## Why some categories have only one concrete class

`BindStrategy` has a single implementation, `DefaultBindStrategy`, with no
Native/Fortran split. Unlike Move, Unbind, and ConflictResolution, bind()'s
random draws use an identical formula in both modes -- the only difference
is where the number comes from, which `RandomDrawSource` already handles.

## Why services are separate from strategies

`services/` (`random_draw.py`, `binding_time.py`) holds dependencies that
strategies use but that aren't themselves polymorphic per mode -- one
`RandomDrawSource` instance and one `BindingTimeFactory` instance are shared
across every strategy for a given simulation, injected via constructor
rather than each strategy building its own.

`BindingTimeFactory` in particular is shared unconditionally by both modes
(see its docstring) -- it's not part of the Native/Fortran split at all,
just a batching optimization on top of a formula both modes agree on.