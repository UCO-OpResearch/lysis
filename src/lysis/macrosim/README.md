# microsim package

Modular replacement for the monolithic `macroscale.py`. This file is the fast
reference for "how do I add a new variant without missing a step."

## Adding a new strategy variant (e.g. a third Move mode)

1. Subclass the relevant ABC in `strategies/base.py` (e.g. `MoveStrategy`).
   Implement every abstract method, otherwise Python will refuse to instantiate
   your class.
2. Add a construction branch for it in `factory.py`'s `make_*` function.
3. Run `test_factory_registry.py` It will fail if your new subclass isn't 
   reachable from the factory, or if the factory references something that
   isn't a proper subclass.
4. Add/update the regression comparison in the test suite if the new mode
   changes simulation output.

If you only do step 1, the registry test is what catches it -- that test
existing and being run in CI is what keeps this pattern from silently
rotting into half-finished subclasses.

## Why services are separate from strategies

`services/` (RNG, recording) holds dependencies that strategies use but
that aren't themselves polymorphic per-mode in the same way -- there's one
`RandomDrawSource` instance per simulation, selected once at construction,
and every strategy takes it as a constructor argument rather than
constructing its own.
