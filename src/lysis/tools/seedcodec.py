"""Encode and decode RNG seed entropy for the ``micro_seed`` / ``macro_seed`` fields.

The seed fields persist the *entropy* handed to :class:`numpy.random.SeedSequence`.
A single field carries two interchangeable textual forms, disambiguated on read so
the stored representation never needs a dataspec version bump:

- **Bare decimal** (e.g. ``"12345"``) — a plain integer entropy.  A legacy
  ``np.uint32`` seed is a bare decimal, so old files and hand-written CSV cells
  reproduce with no conversion.
- **``base58:``-prefixed** (e.g. ``"base58:3xK9q"``) — a full-width entropy drawn
  from the OS, encoded compactly so a human can copy or retype it.

Both forms decode to one non-negative Python ``int`` fed to ``SeedSequence``, so the
Fortran and Python backends stay bit-identical (they share this single parse).

Encoding policy (:func:`encode_seed`): values that fit in 32 bits stay bare decimal
(legacy-friendly); wider entropy is emitted as ``base58:``.

Base58 uses the **Bitcoin alphabet** (omitting the visually ambiguous ``0``, ``O``,
``I``, ``l``) and operates on the integer directly (repeated division by 58), *not*
on a byte string — so there is no leading-zero-byte convention to worry about; the
integer ``0`` encodes to the first alphabet character (``"1"``).

Negative integers are interpreted as 32-bit two's-complement and folded to their
``uint32`` bit pattern (``x & 0xFFFFFFFF``), matching the ``|uint32`` Fortran tag and
preserving reproduction of legacy signed-seed values.
"""

__all__ = ["parse_seed", "encode_seed", "random_entropy", "BASE58_PREFIX"]

import numpy as np

#: Sigil marking a Base58-encoded entropy payload.  ``:`` is not a Base58 character,
#: is CSV- and shell-safe, and reads like a URI scheme.
BASE58_PREFIX = "base58:"

#: Bitcoin Base58 alphabet (no ``0``, ``O``, ``I``, ``l``).
_ALPHABET = "123456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz"
_ALPHABET_INDEX = {ch: i for i, ch in enumerate(_ALPHABET)}

#: 32-bit mask used to fold signed-negative legacy seeds to their uint32 bit pattern.
_UINT32_MASK = 0xFFFFFFFF


def _normalize_int(value: int) -> int:
    """Fold a signed integer to its ``uint32`` bit pattern; pass non-negatives through.

    :param value: An integer seed/entropy value.
    :type value: int
    :return: ``value & 0xFFFFFFFF`` when ``value`` is negative, else ``value``.
    :rtype: int
    """
    return value & _UINT32_MASK if value < 0 else value


def _b58encode_int(value: int) -> str:
    """Encode a non-negative integer as Base58 (Bitcoin alphabet).

    :param value: A non-negative integer.
    :type value: int
    :return: The Base58 representation; ``0`` encodes to ``"1"``.
    :rtype: str
    :raises ValueError: If ``value`` is negative.
    """
    if value < 0:
        raise ValueError(f"cannot Base58-encode a negative integer: {value}")
    if value == 0:
        return _ALPHABET[0]
    digits = []
    while value:
        value, rem = divmod(value, 58)
        digits.append(_ALPHABET[rem])
    return "".join(reversed(digits))


def _b58decode_int(payload: str) -> int:
    """Decode a Base58 (Bitcoin alphabet) string to a non-negative integer.

    :param payload: The Base58 characters (without the ``base58:`` prefix).
    :type payload: str
    :return: The decoded non-negative integer.
    :rtype: int
    :raises ValueError: If ``payload`` is empty or contains a non-Base58 character.
    """
    if not payload:
        raise ValueError("empty Base58 payload")
    result = 0
    for ch in payload:
        try:
            result = result * 58 + _ALPHABET_INDEX[ch]
        except KeyError:
            raise ValueError(
                f"invalid Base58 character {ch!r} in seed payload {payload!r}"
            ) from None
    return result


def parse_seed(token) -> int:
    """Parse a stored/typed seed into the non-negative integer entropy.

    Accepts an integer (legacy numeric HDF5 attr) or a string in one of two forms:
    a bare decimal, or a ``base58:``-prefixed Base58 payload.  The sigil match is
    case-insensitive; the payload is case-sensitive.  Negative integers (legacy
    signed seeds) are folded to their ``uint32`` bit pattern.

    :param token: An ``int``/:class:`numpy.integer`, or a ``str`` seed.
    :type token: int or numpy.integer or str
    :return: The non-negative integer entropy to feed to ``SeedSequence``.
    :rtype: int
    :raises ValueError: If ``token`` is an empty/blank string or malformed.
    :raises TypeError: If ``token`` is neither an integer nor a string.
    """
    if isinstance(token, (int, np.integer)):
        return _normalize_int(int(token))
    if isinstance(token, str):
        text = token.strip()
        if not text:
            raise ValueError("blank seed string has no entropy to parse")
        if text[: len(BASE58_PREFIX)].lower() == BASE58_PREFIX:
            return _b58decode_int(text[len(BASE58_PREFIX):])
        try:
            return _normalize_int(int(text, 10))
        except ValueError:
            raise ValueError(f"malformed seed token: {token!r}") from None
    raise TypeError(f"seed must be an int or str, got {type(token).__name__}")


def encode_seed(value) -> str:
    """Encode an integer entropy into its canonical stored string form.

    Values that fit in 32 bits are kept as a bare decimal (legacy-friendly); wider
    entropy is emitted in ``base58:`` form.  ``parse_seed(encode_seed(n)) == n`` for
    every non-negative ``n`` (and for negative ``n`` after the ``uint32`` fold).

    :param value: The integer entropy to encode.
    :type value: int or numpy.integer
    :return: The canonical seed string.
    :rtype: str
    """
    n = _normalize_int(int(value))
    if n <= _UINT32_MASK:
        return str(n)
    return BASE58_PREFIX + _b58encode_int(n)


def random_entropy() -> int:
    """Draw fresh OS entropy via :class:`numpy.random.SeedSequence`.

    :return: A non-negative integer drawn from a freshly seeded ``SeedSequence``.
    :rtype: int
    """
    return int(np.random.SeedSequence().entropy)
