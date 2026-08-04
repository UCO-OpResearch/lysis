"""Unit tests for :mod:`lysis.tools.kiss`.

Cover the ctypes wrapper around the Marsaglia KISS C library:

* ``_kiss_library_path`` — resolves ``lib/kiss.so`` at the repository root, and
  a missing library raises a helpful ``OSError``.
* ``getstate`` / ``setstate`` / ``seed`` — state round-trip, exact rewind, the
  seed occupying only the fourth state word, and 32-bit truncation.
* ``random`` — scalar and array forms, ``float64`` output (the NumPy 2
  regression: ``np.float_`` was removed), array-equals-scalar-loop equivalence.
* ``integers`` — one- and two-argument bounds, array form, determinism.
* ``kiss32`` / ``mscw`` — 32-bit range; ``mscw`` is independent of the state.
* Algorithm conformance — the C library's output is checked, bit for bit,
  against a pure-Python reference implementation of KISS transcribed from
  ``src/c/kiss.c``.

The C library keeps its state in file-scope globals, so every generator shares
one stream.  Tests that care about specific values therefore pin the state with
``setstate`` rather than relying on construction order.
"""

from pathlib import Path

import numpy as np
import pytest

from lysis.tools.kiss import KissRandomGenerator, _kiss_library_path

pytestmark = pytest.mark.skipif(
    not _kiss_library_path().is_file(),
    reason=f"{_kiss_library_path()} not built; run 'make shared'",
)

_MASK32 = 0xFFFFFFFF

#: The C library's compiled-in default state (``DEFAULT_C``, ``DEFAULT_JSR``,
#: ``DEFAULT_X`` from ``src/c/kiss.h``) with an arbitrary seed word.
_KNOWN_STATE = (129281, 362436069, 123456789, 20260804)


# ---------------------------------------------------------------------------
#  Pure-Python reference implementation (transcribed from src/c/kiss.c)
# ---------------------------------------------------------------------------


class _ReferenceKiss:
    """Independent re-implementation of the KISS generator in ``src/c/kiss.c``.

    Used to verify that the compiled library really runs the documented
    algorithm, rather than pinning golden values captured from the library
    itself (which would pass even if the library were wrong).
    """

    #: ``1 / (2**64 - 1)``, the scaling applied by ``urcw1_``.
    _SCALE = 1.0 / 18446744073709551615.0

    def __init__(self, state):
        self.c, self.jsr, self.x, self.y = (int(w) for w in state)

    def kiss32(self) -> int:
        """Return the next 32-bit draw, mirroring ``kiss32_``."""
        # Congruential generator
        self.y = (69069 * self.y + 12345) & _MASK32

        # 3-shift-register (xorshift) generator
        jsr = self.jsr
        jsr ^= (jsr << 13) & _MASK32
        jsr ^= jsr >> 17
        jsr ^= (jsr << 5) & _MASK32
        self.jsr = jsr & _MASK32

        # Multiply-with-carry generator (64-bit intermediate)
        t = 333333314 * self.x + self.c
        self.c = t >> 32
        self.x = (t + self.c) & _MASK32
        if self.x < self.c:
            self.x = (self.x + 1) & _MASK32
            self.c = (self.c + 1) & _MASK32
        self.x = (~self.x + 1) & _MASK32

        return (self.x + self.y + self.jsr) & _MASK32

    def random(self) -> float:
        """Return the next U(0,1) value, mirroring ``urcw1_``.

        Statement order matches the C source so the double-precision rounding
        is identical.
        """
        result = 4294967296.0 * float(self.kiss32())
        result += float(self.kiss32())
        result *= self._SCALE
        return min(result, 1.0)


@pytest.fixture
def kiss():
    """A generator pinned to :data:`_KNOWN_STATE`.

    The underlying C state is global, so pinning it here makes each test
    independent of the order the suite happens to run in.
    """
    generator = KissRandomGenerator()
    generator.setstate(_KNOWN_STATE)
    return generator


class TestLibraryLocation:
    """The shared library is found relative to the repository root."""

    def test_path_points_at_repo_lib(self):
        """``lib/kiss.so`` sits three parents above ``tools/kiss.py``."""
        path = _kiss_library_path()
        assert path.name == "kiss.so"
        assert path.parent.name == "lib"
        # The repo root is the parent of both lib/ and src/lysis/
        assert (path.parent.parent / "src" / "lysis" / "tools" / "kiss.py").is_file()

    def test_path_is_absolute(self):
        assert _kiss_library_path().is_absolute()

    def test_missing_library_raises_helpful_oserror(self, monkeypatch, tmp_path):
        """A missing kiss.so names the path and the build command."""
        missing = tmp_path / "lib" / "kiss.so"
        monkeypatch.setattr(
            "lysis.tools.kiss._kiss_library_path", lambda: Path(missing)
        )
        with pytest.raises(OSError, match="make shared"):
            KissRandomGenerator(1)


class TestState:
    """getstate/setstate/seed manipulate the four-word KISS state."""

    def test_getstate_returns_four_uint32_words(self, kiss):
        state = kiss.getstate()
        assert len(state) == 4
        assert all(isinstance(word, int) for word in state)
        assert all(0 <= word <= _MASK32 for word in state)

    def test_setstate_round_trip(self, kiss):
        kiss.setstate(_KNOWN_STATE)
        assert kiss.getstate() == _KNOWN_STATE

    def test_setstate_rewinds_the_sequence(self, kiss):
        """Restoring a saved state reproduces the sequence exactly."""
        state = kiss.getstate()
        first = kiss.random(10)
        kiss.setstate(state)
        second = kiss.random(10)
        assert np.array_equal(first, second)

    def test_setstate_truncates_to_32_bits(self, kiss):
        """Words wider than 32 bits are truncated, as documented."""
        kiss.setstate((0, 0, 0, 2**32 + 7))
        assert kiss.getstate()[3] == 7

    def test_seed_sets_only_the_fourth_word(self, kiss):
        head = kiss.getstate()[:3]
        kiss.seed(4242)
        assert kiss.getstate() == (*head, 4242)

    def test_seed_accepts_numpy_integers(self, kiss):
        kiss.seed(np.int64(4242))
        assert kiss.getstate()[3] == 4242

    def test_constructor_seed_lands_in_state(self):
        assert KissRandomGenerator(31337).getstate()[3] == 31337

    def test_constructor_without_seed_uses_clock(self, kiss):
        """A seedless generator takes its seed word from ``mscw``."""
        head = kiss.getstate()[:3]
        seedless = KissRandomGenerator()
        state = seedless.getstate()
        # The clock-based seed replaced only the fourth word
        assert state[:3] == head
        assert 0 <= state[3] <= _MASK32

    def test_state_is_shared_across_instances(self, kiss):
        """All instances share the C library's file-scope state (documented)."""
        other = KissRandomGenerator(1234)
        assert kiss.getstate() == other.getstate()
        other.random()
        assert kiss.getstate() == other.getstate()


class TestRandom:
    """random() returns U(0,1) scalars and float64 arrays."""

    def test_scalar_is_a_float_in_the_unit_interval(self, kiss):
        value = kiss.random()
        assert isinstance(value, float)
        assert 0.0 <= value <= 1.0

    def test_array_shape_and_dtype(self, kiss):
        """The array path is the NumPy 2 regression: ``np.float_`` was removed."""
        out = kiss.random(7)
        assert isinstance(out, np.ndarray)
        assert out.shape == (7,)
        assert out.dtype == np.float64

    def test_array_values_in_the_unit_interval(self, kiss):
        out = kiss.random(1000)
        assert np.all(out >= 0.0)
        assert np.all(out <= 1.0)

    def test_array_accepts_numpy_integer_size(self, kiss):
        """Sizes from e.g. ``np.count_nonzero`` are numpy ints, not ints."""
        out = kiss.random(np.int64(5))
        assert out.shape == (5,)

    def test_zero_size_returns_empty_array(self, kiss):
        out = kiss.random(0)
        assert out.shape == (0,)
        assert out.dtype == np.float64

    def test_array_matches_successive_scalar_draws(self, kiss):
        """The vectorised C loop yields the same stream as repeated urcw1."""
        state = kiss.getstate()
        vector = kiss.random(16)
        kiss.setstate(state)
        scalars = np.array([kiss.random() for _ in range(16)])
        assert np.array_equal(vector, scalars)

    def test_draws_are_distinct(self, kiss):
        """A large sample has no repeats (a stuck generator would repeat)."""
        out = kiss.random(1000)
        assert len(np.unique(out)) == 1000

    def test_sample_mean_is_near_one_half(self, kiss):
        """Crude uniformity check over a large sample."""
        out = kiss.random(50000)
        assert out.mean() == pytest.approx(0.5, abs=0.01)


class TestIntegers:
    """integers() mimics numpy.random.Generator.integers()."""

    def test_single_argument_uses_zero_as_lower_bound(self, kiss):
        for _ in range(100):
            value = kiss.integers(10)
            assert isinstance(value, int)
            assert 0 <= value < 10

    def test_two_argument_bounds(self, kiss):
        for _ in range(100):
            value = kiss.integers(5, 15)
            assert 5 <= value < 15

    def test_array_form(self, kiss):
        out = kiss.integers(0, 8, size=200)
        assert isinstance(out, np.ndarray)
        assert out.shape == (200,)
        assert np.issubdtype(out.dtype, np.integer)
        assert np.all(out >= 0)
        assert np.all(out < 8)

    def test_degenerate_range_is_constant(self, kiss):
        assert np.all(kiss.integers(3, 4, size=20) == 3)

    def test_deterministic_given_state(self, kiss):
        state = kiss.getstate()
        first = kiss.integers(0, 1000, size=50)
        kiss.setstate(state)
        assert np.array_equal(first, kiss.integers(0, 1000, size=50))

    def test_covers_the_range(self, kiss):
        """Every value in a small range appears in a large sample."""
        out = kiss.integers(8, size=2000)
        assert set(np.unique(out)) == set(range(8))


class TestRawGenerators:
    """The C functions exposed directly: kiss32 and mscw."""

    def test_kiss32_is_a_32_bit_value(self, kiss):
        for _ in range(100):
            value = kiss.kiss32()
            assert 0 <= value <= _MASK32

    def test_kiss32_advances_the_state(self, kiss):
        before = kiss.getstate()
        kiss.kiss32()
        assert kiss.getstate() != before

    def test_mscw_is_a_32_bit_value(self, kiss):
        assert 0 <= kiss.mscw() <= _MASK32

    def test_mscw_does_not_disturb_the_state(self, kiss):
        """The clock-based seed helper is independent of the KISS stream."""
        before = kiss.getstate()
        kiss.mscw()
        assert kiss.getstate() == before


class TestAlgorithmConformance:
    """The compiled library matches an independent Python implementation.

    These are the tests that would catch a miscompiled or substituted
    ``kiss.so``, and that pin the exact stream the Fortran model reproduces.
    """

    def test_kiss32_sequence(self, kiss):
        reference = _ReferenceKiss(_KNOWN_STATE)
        assert [kiss.kiss32() for _ in range(100)] == [
            reference.kiss32() for _ in range(100)
        ]

    def test_state_evolution(self, kiss):
        reference = _ReferenceKiss(_KNOWN_STATE)
        for _ in range(20):
            kiss.kiss32()
            reference.kiss32()
            assert kiss.getstate() == (
                reference.c,
                reference.jsr,
                reference.x,
                reference.y,
            )

    def test_random_scalars_bit_for_bit(self, kiss):
        reference = _ReferenceKiss(_KNOWN_STATE)
        for _ in range(100):
            assert kiss.random() == reference.random()

    def test_random_consumes_two_draws_per_value(self, kiss):
        """``urcw1_`` builds each double from two 32-bit draws."""
        state = kiss.getstate()
        value = kiss.random()
        kiss.setstate(state)
        high, low = kiss.kiss32(), kiss.kiss32()
        expected = (4294967296.0 * float(high) + float(low)) / 18446744073709551615.0
        assert value == pytest.approx(expected, rel=1e-15)

    def test_random_array_bit_for_bit(self, kiss):
        reference = _ReferenceKiss(_KNOWN_STATE)
        expected = np.array([reference.random() for _ in range(64)])
        assert np.array_equal(kiss.random(64), expected)

    @pytest.mark.parametrize("seed", [0, 1, 42, 123456789, _MASK32])
    def test_matches_reference_from_default_state(self, kiss, seed):
        """Conformance holds across the seed range, from the default state."""
        state = (129281, 362436069, 123456789, seed)
        kiss.setstate(state)
        reference = _ReferenceKiss(state)
        assert np.array_equal(
            kiss.random(32), np.array([reference.random() for _ in range(32)])
        )
