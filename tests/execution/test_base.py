"""Unit tests for :mod:`lysis.execution.base` — SimulationRunner."""

import pytest

from lysis.execution.base import SimulationRunner


class TestSimulationRunner:
    """Tests for :class:`SimulationRunner`."""

    def test_cannot_instantiate_directly(self):
        """SimulationRunner is abstract and must not be instantiatable."""
        with pytest.raises(TypeError):
            SimulationRunner()

    def test_subclass_without_execute_is_abstract(self):
        """A subclass that omits execute() must also be uninstantiatable."""
        class Incomplete(SimulationRunner):
            pass

        with pytest.raises(TypeError):
            Incomplete()

    def test_concrete_subclass_instantiates(self):
        """A subclass that implements execute() can be instantiated."""
        class Concrete(SimulationRunner):
            def execute(self) -> None:
                pass

        obj = Concrete()
        assert isinstance(obj, SimulationRunner)

    def test_execute_is_abstract(self):
        """execute must be listed in __abstractmethods__."""
        assert "execute" in SimulationRunner.__abstractmethods__
