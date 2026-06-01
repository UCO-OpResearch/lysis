"""Unit tests for :mod:`lysis.tools.seedcodec`.

Cover the public API and the Base58 helpers:

* ``encode_seed`` / ``parse_seed`` — round-trip across the 32-bit boundary, the
  decimal-vs-``base58:`` representation policy, legacy-int and signed-negative input,
  case-insensitive sigil, and malformed-input rejection.
* ``_b58encode_int`` / ``_b58decode_int`` — alphabet correctness, zero, round-trip.
* ``random_entropy`` — wide, distinct draws.
"""

import numpy as np
import pytest

from lysis.tools.seedcodec import (
    BASE58_PREFIX,
    _b58decode_int,
    _b58encode_int,
    encode_seed,
    parse_seed,
    random_entropy,
)

_UINT32_MAX = 0xFFFFFFFF
_WIDE = 0x0123456789ABCDEF0123456789ABCDEF  # 128-bit


class TestRoundTrip:
    """parse_seed(encode_seed(n)) == n across the representation boundary."""

    @pytest.mark.parametrize(
        "n", [0, 1, 12345, _UINT32_MAX, _UINT32_MAX + 1, _WIDE, 2**128 - 1]
    )
    def test_round_trip(self, n):
        assert parse_seed(encode_seed(n)) == n

    def test_small_values_stay_decimal(self):
        """Values that fit in 32 bits encode as a bare decimal string."""
        assert encode_seed(0) == "0"
        assert encode_seed(12345) == "12345"
        assert encode_seed(_UINT32_MAX) == str(_UINT32_MAX)

    def test_wide_values_use_base58(self):
        """Values beyond 32 bits encode with the base58: sigil."""
        encoded = encode_seed(_UINT32_MAX + 1)
        assert encoded.startswith(BASE58_PREFIX)
        assert parse_seed(encoded) == _UINT32_MAX + 1


class TestParseSeed:
    """parse_seed accepts ints, decimals, and base58: strings; rejects junk."""

    def test_bare_decimal(self):
        assert parse_seed("12345") == 12345

    def test_legacy_int(self):
        assert parse_seed(99) == 99
        assert parse_seed(np.uint32(99)) == 99

    def test_negative_int_folds_to_uint32(self):
        """A signed-negative legacy seed maps to its uint32 bit pattern."""
        assert parse_seed(-1) == _UINT32_MAX
        assert parse_seed(-559038737) == 0xDEADBEEF
        assert parse_seed("-559038737") == 0xDEADBEEF

    def test_base58_sigil_case_insensitive(self):
        payload = encode_seed(_WIDE)[len(BASE58_PREFIX):]
        assert parse_seed("base58:" + payload) == _WIDE
        assert parse_seed("BASE58:" + payload) == _WIDE

    def test_whitespace_trimmed(self):
        assert parse_seed("  4242  ") == 4242

    @pytest.mark.parametrize("bad", ["", "   ", "not-a-number", "12x34", "base58:"])
    def test_malformed_raises(self, bad):
        with pytest.raises(ValueError):
            parse_seed(bad)

    def test_base58_invalid_char_raises(self):
        """0, O, I, l are not in the Base58 alphabet."""
        with pytest.raises(ValueError):
            parse_seed("base58:0OIl")

    def test_wrong_type_raises(self):
        with pytest.raises(TypeError):
            parse_seed(1.5)


class TestBase58Helpers:
    """Direct checks on the integer Base58 codec."""

    def test_zero_is_first_char(self):
        assert _b58encode_int(0) == "1"
        assert _b58decode_int("1") == 0

    def test_negative_encode_raises(self):
        with pytest.raises(ValueError):
            _b58encode_int(-1)

    def test_empty_decode_raises(self):
        with pytest.raises(ValueError):
            _b58decode_int("")

    @pytest.mark.parametrize("n", [0, 1, 57, 58, 59, 12345, _WIDE])
    def test_int_round_trip(self, n):
        assert _b58decode_int(_b58encode_int(n)) == n


class TestRandomEntropy:
    """random_entropy draws wide, distinct values."""

    def test_distinct_and_wide(self):
        draws = {random_entropy() for _ in range(8)}
        assert len(draws) == 8
        assert all(d >= 0 for d in draws)
        assert max(draws) > _UINT32_MAX  # SeedSequence entropy is ~128-bit
