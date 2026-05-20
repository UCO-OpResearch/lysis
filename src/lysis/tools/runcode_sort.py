"""Intelligent sorting for run codes containing roman numerals or python-style integers.

Run codes such as ``TF-vii__447_798`` embed two kinds of numeric tokens that sort
incorrectly when treated as plain strings:

- **Roman numerals** (e.g. ``v``, ``vii``, ``ix``, ``xi``, ``xiii``) — sorted by
  numeric value rather than alphabetically.
- **Python-style integers** (e.g. ``1_582_867``, ``9_951``) — digits optionally
  separated by underscores as thousands separators, converted to ``int`` before
  comparison.

:func:`smart_sort` tokenises a collection of codes, detects positions where
*every* code carries the same kind of numeric token, and builds a composite sort
key that sorts those positions numerically while leaving purely textual positions
in lexicographic order.  If the codes cannot be tokenised consistently (different
token counts) the function falls back to a plain lexicographic sort.
"""

__author__ = "Bradley Paynter"
__copyright__ = "Copyright 2026, Bradley Paynter"
__version__ = "0.1"
__status__ = "Development"

__all__ = ["smart_sort"]

import re

# ---------------------------------------------------------------------------
# Tokenisation
# ---------------------------------------------------------------------------

# Delimiter: __ (double underscore) or - (hyphen).  __ is tried first so that a
# double underscore is consumed as one delimiter rather than two hyphens.
_DELIM_RE = re.compile(r"(?:__|-)")


def _tokenize(code: str) -> list[str]:
    """Split *code* on ``__`` or ``-``, preserving single-underscore digit groups."""
    return _DELIM_RE.split(code)


# ---------------------------------------------------------------------------
# Roman numeral helpers
# ---------------------------------------------------------------------------

# Valid roman numeral pattern (case-insensitive).  fullmatch is used so that
# plain strings containing only Roman characters (e.g. "d", "c") are still
# accepted only when they form a structurally valid numeral.
_ROMAN_RE = re.compile(
    r"^M{0,4}(?:CM|CD|D?C{0,3})(?:XC|XL|L?X{0,3})(?:IX|IV|V?I{0,3})$",
    re.IGNORECASE,
)
_ROMAN_MAP = {"I": 1, "V": 5, "X": 10, "L": 50, "C": 100, "D": 500, "M": 1000}


def _is_roman(token: str) -> bool:
    """Return *True* if *token* is a non-empty, structurally valid roman numeral."""
    return bool(token) and _ROMAN_RE.fullmatch(token) is not None


def _roman_to_int(s: str) -> int:
    """Convert a roman numeral string to an integer (case-insensitive)."""
    result, prev = 0, 0
    for ch in reversed(s.upper()):
        val = _ROMAN_MAP[ch]
        result += val if val >= prev else -val
        prev = val
    return result


# ---------------------------------------------------------------------------
# Python-style integer helpers
# ---------------------------------------------------------------------------

# Matches digit strings optionally separated by single underscores used as
# thousands separators (e.g. 1_582_867, 9_951).  Must start and end with a digit.
_PYINT_RE = re.compile(r"^\d+(_\d+)*$")


def _is_pyint(token: str) -> bool:
    """Return *True* if *token* is a python-style integer literal."""
    return _PYINT_RE.fullmatch(token) is not None


def _pyint_value(token: str) -> int:
    """Convert a python-style integer string (with ``_`` separators) to ``int``."""
    return int(token.replace("_", ""))


# ---------------------------------------------------------------------------
# Position classification
# ---------------------------------------------------------------------------

_ROMAN = "roman"
_PYINT = "pyint"
_STR = "str"


def _classify_position(tokens: list[str]) -> str:
    """Return the detected type for all tokens at a single code position.

    :param tokens: one token per code, all at the same positional index.
    :return: ``'roman'``, ``'pyint'``, or ``'str'``.
    """
    if all(_is_roman(t) for t in tokens):
        return _ROMAN
    if all(_is_pyint(t) for t in tokens):
        return _PYINT
    return _STR


# ---------------------------------------------------------------------------
# Sort key construction
# ---------------------------------------------------------------------------

def _make_key(tokenized: list[str], position_types: list[str]) -> tuple:
    """Build a sort-key tuple for one tokenised code.

    Each element is an ``int`` for numeric positions or a ``str`` for textual
    positions.  Because every code at a given position has the same Python type,
    cross-type comparisons never occur within the same tuple slot.

    :param tokenized: list of token strings for a single code.
    :param position_types: parallel list of position type labels.
    :return: tuple suitable for use as a ``sorted()`` key.
    """
    key = []
    for token, ptype in zip(tokenized, position_types):
        if ptype == _ROMAN:
            key.append(_roman_to_int(token))
        elif ptype == _PYINT:
            key.append(_pyint_value(token))
        else:
            key.append(token)
    return tuple(key)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def smart_sort(codes: list[str]) -> list[str]:
    """Sort *codes* intelligently, treating embedded numerics as numbers.

    The algorithm:

    1. Tokenise every code by splitting on ``__`` and ``-``.
    2. If the codes have inconsistent token counts, fall back to
       lexicographic ``sorted()``.
    3. For each token position, classify it as ``'roman'``, ``'pyint'``,
       or ``'str'`` based on whether *all* codes carry that type of token there.
    4. Build a composite sort key that converts numeric positions to ``int``
       while keeping textual positions as ``str``.
    5. Return the codes sorted by that key.

    :param codes: run code strings to sort.
    :type codes: list[str]
    :return: new list with codes in intelligent sort order.
    :rtype: list[str]
    """
    if not codes:
        return []

    tokenized = [_tokenize(c) for c in codes]
    lengths = {len(t) for t in tokenized}

    if len(lengths) > 1:
        # Inconsistent structure — fall back to plain lexicographic sort
        return sorted(codes)

    n_pos = next(iter(lengths))
    position_types = [
        _classify_position([t[i] for t in tokenized]) for i in range(n_pos)
    ]

    pairs = [
        (_make_key(tok, position_types), code)
        for tok, code in zip(tokenized, codes)
    ]
    pairs.sort(key=lambda x: x[0])
    return [code for _, code in pairs]
