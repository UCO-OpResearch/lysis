"""Unit tests for :mod:`lysis.tools.runcode_sort`.

Tests cover every public and private symbol in the module:

* ``_tokenize``             — delimiter splitting, double-underscore priority,
                              single-underscore preservation inside numbers
* ``_is_roman``             — valid roman numerals, invalid strings, edge cases
* ``_roman_to_int``         — additive and subtractive notation, case-insensitivity
* ``_is_pyint``             — valid python-style integers, reject malformed tokens
* ``_pyint_value``          — correct integer value after stripping underscores
* ``_classify_position``    — unanimous roman, pyint, or str; mixed falls back
* ``smart_sort``            — primary public API: roman ordering, integer ordering,
                              compound codes, lexicographic fallback, edge cases
"""

import pytest

from lysis.tools.runcode_sort import (
    _classify_position,
    _is_pyint,
    _is_roman,
    _pyint_value,
    _roman_to_int,
    _tokenize,
    smart_sort,
)

# ---------------------------------------------------------------------------
# _tokenize
# ---------------------------------------------------------------------------


class TestTokenize:
    """_tokenize splits run codes on ``__`` and ``-``, preserving single ``_``."""

    def test_hyphen_delimiter(self):
        """A single hyphen splits into two tokens."""
        assert _tokenize("TB-xi") == ["TB", "xi"]

    def test_double_underscore_delimiter(self):
        """A double underscore splits into two tokens."""
        assert _tokenize("TB__xi") == ["TB", "xi"]

    def test_compound_code(self):
        """Real run code produces three tokens."""
        assert _tokenize("TF-vii__447_798") == ["TF", "vii", "447_798"]

    def test_single_underscore_preserved_in_number(self):
        """Single underscores inside a number token are not consumed as delimiters."""
        assert _tokenize("TB-xi__1_582_867") == ["TB", "xi", "1_582_867"]

    def test_double_underscore_consumed_before_hyphens(self):
        """``__`` is matched as one delimiter, not as two separate tokens."""
        tokens = _tokenize("A-B__C-D")
        assert tokens == ["A", "B", "C", "D"]

    def test_no_delimiter(self):
        """A code with no recognised delimiter returns a single token."""
        assert _tokenize("ABCDEF") == ["ABCDEF"]

    def test_multiple_hyphens(self):
        """Multiple hyphens each act as individual delimiters."""
        assert _tokenize("a-b-c") == ["a", "b", "c"]


# ---------------------------------------------------------------------------
# _is_roman
# ---------------------------------------------------------------------------


class TestIsRoman:
    """_is_roman accepts structurally valid roman numerals only."""

    # --- valid ---
    @pytest.mark.parametrize("token", ["i", "v", "x", "l", "c", "d", "m"])
    def test_single_valid_letter(self, token):
        """Each single valid roman letter is recognised (case-insensitive)."""
        assert _is_roman(token)
        assert _is_roman(token.upper())

    @pytest.mark.parametrize(
        "token,expected_int",
        [
            ("iv", 4),
            ("ix", 9),
            ("xi", 11),
            ("vii", 7),
            ("xiii", 13),
            ("xl", 40),
            ("xc", 90),
            ("cd", 400),
            ("cm", 900),
        ],
    )
    def test_valid_compound_numerals(self, token, expected_int):
        """Common compound roman numerals are accepted."""
        assert _is_roman(token)

    def test_case_insensitive_upper(self):
        """Uppercase roman numeral strings are valid."""
        assert _is_roman("XI")
        assert _is_roman("VII")

    def test_case_insensitive_mixed(self):
        """Mixed-case roman numeral strings are valid."""
        assert _is_roman("Xi")
        assert _is_roman("iX")

    # --- invalid ---
    def test_empty_string_rejected(self):
        """Empty string is not a roman numeral."""
        assert not _is_roman("")

    def test_plain_string_rejected(self):
        """A regular word is not a roman numeral."""
        assert not _is_roman("hello")

    def test_string_with_non_roman_letters_rejected(self):
        """A string containing letters outside the roman set is rejected."""
        assert not _is_roman("TB")
        assert not _is_roman("TF")

    def test_digits_rejected(self):
        """A digit string is not a roman numeral."""
        assert not _is_roman("42")

    def test_underscore_rejected(self):
        """A token with underscore is not a roman numeral."""
        assert not _is_roman("1_582_867")

    def test_structurally_invalid_rejected(self):
        """Strings that use only roman letters but violate structure are rejected."""
        # 'iiii' is not valid standard notation
        assert not _is_roman("iiii")
        # 'vv' is invalid (only one V allowed per numeral)
        assert not _is_roman("vv")


# ---------------------------------------------------------------------------
# _roman_to_int
# ---------------------------------------------------------------------------


class TestRomanToInt:
    """_roman_to_int converts roman numerals to their integer values."""

    @pytest.mark.parametrize(
        "token,value",
        [
            ("i", 1),
            ("v", 5),
            ("x", 10),
            ("l", 50),
            ("c", 100),
            ("d", 500),
            ("m", 1000),
        ],
    )
    def test_single_letters(self, token, value):
        """Each single roman letter maps to its standard value."""
        assert _roman_to_int(token) == value

    @pytest.mark.parametrize(
        "token,value",
        [
            ("iv", 4),
            ("ix", 9),
            ("xl", 40),
            ("xc", 90),
            ("cd", 400),
            ("cm", 900),
        ],
    )
    def test_subtractive_pairs(self, token, value):
        """Subtractive notation is handled correctly."""
        assert _roman_to_int(token) == value

    @pytest.mark.parametrize(
        "token,value",
        [
            ("vii", 7),
            ("xi", 11),
            ("xiii", 13),
            ("xiv", 14),
            ("xlii", 42),
            ("mcmxcix", 1999),
        ],
    )
    def test_compound_numerals(self, token, value):
        """Compound roman numerals convert to the correct integer."""
        assert _roman_to_int(token) == value

    def test_case_insensitive(self):
        """_roman_to_int is case-insensitive."""
        assert _roman_to_int("XI") == _roman_to_int("xi") == 11

    def test_lysis_front_roman_ordering(self):
        """The roman numerals used in lysis-front run codes sort correctly."""
        tokens = ["xi", "ix", "xiii", "v", "vii", "x"]
        values = [_roman_to_int(t) for t in tokens]
        assert sorted(tokens, key=_roman_to_int) == ["v", "vii", "ix", "x", "xi", "xiii"]
        assert sorted(values) == [5, 7, 9, 10, 11, 13]


# ---------------------------------------------------------------------------
# _is_pyint
# ---------------------------------------------------------------------------


class TestIsPyint:
    """_is_pyint recognises python-style integer literals."""

    # --- valid ---
    @pytest.mark.parametrize(
        "token",
        ["0", "42", "1000", "9_951", "316_573", "949_720", "1_582_867"],
    )
    def test_valid_tokens(self, token):
        """Bare digits and underscore-separated digit groups are accepted."""
        assert _is_pyint(token)

    def test_large_number(self):
        """Large numbers with multiple underscore groups are accepted."""
        assert _is_pyint("1_000_000_000")

    # --- invalid ---
    def test_leading_underscore_rejected(self):
        """Tokens with a leading underscore are not python integers."""
        assert not _is_pyint("_123")

    def test_trailing_underscore_rejected(self):
        """Tokens with a trailing underscore are not python integers."""
        assert not _is_pyint("123_")

    def test_double_underscore_rejected(self):
        """Double underscores within the token are not python integers."""
        assert not _is_pyint("1__000")

    def test_roman_numeral_rejected(self):
        """Roman numerals are not python integers."""
        assert not _is_pyint("xi")

    def test_plain_string_rejected(self):
        """Plain alphabetic strings are not python integers."""
        assert not _is_pyint("hello")

    def test_empty_string_rejected(self):
        """Empty string is not a python integer."""
        assert not _is_pyint("")

    def test_float_rejected(self):
        """Floating-point strings are not python integers."""
        assert not _is_pyint("3.14")


# ---------------------------------------------------------------------------
# _pyint_value
# ---------------------------------------------------------------------------


class TestPyintValue:
    """_pyint_value converts python-style integer strings to int."""

    def test_plain_digits(self):
        """Bare digit strings convert to the expected integer."""
        assert _pyint_value("42") == 42

    def test_single_underscore(self):
        """Single underscore is stripped before conversion."""
        assert _pyint_value("9_951") == 9951

    def test_multiple_underscores(self):
        """Multiple underscores are all stripped before conversion."""
        assert _pyint_value("1_582_867") == 1582867

    def test_three_groups(self):
        """Three-group numbers convert correctly."""
        assert _pyint_value("316_573") == 316573

    @pytest.mark.parametrize(
        "token,expected",
        [
            ("9_951", 9951),
            ("21_105", 21105),
            ("149_266", 149266),
            ("316_573", 316573),
            ("447_798", 447798),
            ("746_330", 746330),
            ("949_720", 949720),
            ("1_582_867", 1582867),
        ],
    )
    def test_lysis_front_values(self, token, expected):
        """All integer tokens from lysis-front run codes convert correctly."""
        assert _pyint_value(token) == expected


# ---------------------------------------------------------------------------
# _classify_position
# ---------------------------------------------------------------------------


class TestClassifyPosition:
    """_classify_position identifies the dominant type across all tokens at one position."""

    def test_all_roman(self):
        """All roman numeral tokens → 'roman'."""
        assert _classify_position(["v", "vii", "ix", "x", "xi", "xiii"]) == "roman"

    def test_all_pyint(self):
        """All python-style integer tokens → 'pyint'."""
        assert _classify_position(["9_951", "316_573", "949_720", "1_582_867"]) == "pyint"

    def test_all_plain_strings(self):
        """All non-numeric tokens → 'str'."""
        assert _classify_position(["TB", "TF"]) == "str"

    def test_mixed_roman_and_string_falls_back(self):
        """A mix of roman and plain strings → 'str'."""
        assert _classify_position(["xi", "hello"]) == "str"

    def test_mixed_pyint_and_string_falls_back(self):
        """A mix of integers and plain strings → 'str'."""
        assert _classify_position(["42", "hello"]) == "str"

    def test_mixed_roman_and_pyint_falls_back(self):
        """A mix of roman and integer tokens → 'str' (not a pure numeric type)."""
        assert _classify_position(["xi", "42"]) == "str"

    def test_single_roman_token(self):
        """A list with a single roman numeral → 'roman'."""
        assert _classify_position(["ix"]) == "roman"

    def test_single_pyint_token(self):
        """A list with a single python integer → 'pyint'."""
        assert _classify_position(["1_000"]) == "pyint"


# ---------------------------------------------------------------------------
# smart_sort — edge cases
# ---------------------------------------------------------------------------


class TestSmartSortEdgeCases:
    """smart_sort handles degenerate inputs gracefully."""

    def test_empty_list(self):
        """Empty input returns an empty list."""
        assert smart_sort([]) == []

    def test_single_element(self):
        """A one-element list is returned unchanged."""
        assert smart_sort(["TB-xi__1_582_867"]) == ["TB-xi__1_582_867"]

    def test_already_sorted(self):
        """A correctly-ordered list is returned in the same order."""
        codes = ["TB-ix__21_105", "TB-ix__316_573", "TB-ix__949_720"]
        assert smart_sort(codes) == codes

    def test_inconsistent_token_counts_fallback(self):
        """Codes with different token counts fall back to plain lexicographic sort."""
        codes = ["z-long__extra__token", "a", "m-short"]
        assert smart_sort(codes) == sorted(codes)

    def test_all_plain_strings(self):
        """Codes with no numeric tokens are sorted lexicographically."""
        codes = ["zebra", "apple", "mango"]
        assert smart_sort(codes) == ["apple", "mango", "zebra"]

    def test_returns_new_list(self):
        """smart_sort does not modify the input list."""
        codes = ["b", "a", "c"]
        original = codes[:]
        smart_sort(codes)
        assert codes == original


# ---------------------------------------------------------------------------
# smart_sort — roman numeral ordering
# ---------------------------------------------------------------------------


class TestSmartSortRoman:
    """smart_sort orders roman numeral tokens by numeric value, not alphabetically."""

    def test_simple_roman_codes(self):
        """Codes with only a roman numeral are sorted by value."""
        codes = ["xiii", "v", "xi", "vii", "ix", "x"]
        assert smart_sort(codes) == ["v", "vii", "ix", "x", "xi", "xiii"]

    def test_roman_vs_alpha_differs(self):
        """Smart order differs from alphabetical for roman numerals."""
        codes = ["xiii", "v", "xi", "ix", "x", "vii"]
        assert smart_sort(codes) != sorted(codes)

    def test_prefixed_roman(self):
        """Roman numeral in the second token position is sorted numerically."""
        codes = ["A-xiii", "A-v", "A-xi", "A-ix"]
        assert smart_sort(codes) == ["A-v", "A-ix", "A-xi", "A-xiii"]

    def test_roman_uppercase(self):
        """Uppercase roman numeral tokens sort by value."""
        codes = ["XIII", "V", "IX", "XI"]
        assert smart_sort(codes) == ["V", "IX", "XI", "XIII"]


# ---------------------------------------------------------------------------
# smart_sort — python-style integer ordering
# ---------------------------------------------------------------------------


class TestSmartSortPyint:
    """smart_sort orders python-style integer tokens by numeric value."""

    def test_simple_pyint_codes(self):
        """Integer-only codes sort numerically, not lexicographically."""
        codes = ["1_582_867", "9_951", "316_573", "949_720", "21_105", "149_266"]
        assert smart_sort(codes) == [
            "9_951", "21_105", "149_266", "316_573", "949_720", "1_582_867"
        ]

    def test_pyint_vs_alpha_differs(self):
        """Numeric order differs from lexicographic for underscore-separated numbers."""
        codes = ["1_582_867", "9_951", "316_573"]
        # Lexicographic: "1_582_867" < "316_573" < "9_951"
        assert sorted(codes) == ["1_582_867", "316_573", "9_951"]
        # Smart: 9951 < 316573 < 1582867
        assert smart_sort(codes) == ["9_951", "316_573", "1_582_867"]

    def test_plain_integers_without_underscores(self):
        """Plain digit tokens (no underscores) also sort numerically."""
        codes = ["100", "9", "42", "1000"]
        assert smart_sort(codes) == ["9", "42", "100", "1000"]


# ---------------------------------------------------------------------------
# smart_sort — compound codes (lysis-front run codes)
# ---------------------------------------------------------------------------


class TestSmartSortCompound:
    """smart_sort correctly handles the two-level numeric embedding in lysis-front codes."""

    # The 24 codes from data/lysis-front/ in the order they appear in the notebook.
    LYSIS_FRONT_CODES = [
        "TB-xi__1_582_867", "TB-xi__949_720", "TB-xi__316_573", "TB-xi__21_105",
        "TF-x__746_330", "TF-x__447_798", "TF-x__149_266", "TF-x__9_951",
        "TF-vii__746_330", "TF-vii__447_798", "TF-vii__149_266", "TF-vii__9_951",
        "TF-v__746_330", "TF-v__447_798", "TF-v__149_266", "TF-v__9_951",
        "TB-ix__1_582_867", "TB-ix__949_720", "TB-ix__316_573", "TB-ix__21_105",
        "TB-xiii__1_582_867", "TB-xiii__949_720", "TB-xiii__316_573", "TB-xiii__21_105",
    ]

    EXPECTED_ORDER = [
        # TB group (alpha before TF), roman ix(9) < xi(11) < xiii(13),
        # integer ascending within each roman group
        "TB-ix__21_105", "TB-ix__316_573", "TB-ix__949_720", "TB-ix__1_582_867",
        "TB-xi__21_105", "TB-xi__316_573", "TB-xi__949_720", "TB-xi__1_582_867",
        "TB-xiii__21_105", "TB-xiii__316_573", "TB-xiii__949_720", "TB-xiii__1_582_867",
        # TF group, roman v(5) < vii(7) < x(10), integer ascending
        "TF-v__9_951", "TF-v__149_266", "TF-v__447_798", "TF-v__746_330",
        "TF-vii__9_951", "TF-vii__149_266", "TF-vii__447_798", "TF-vii__746_330",
        "TF-x__9_951", "TF-x__149_266", "TF-x__447_798", "TF-x__746_330",
    ]

    def test_full_lysis_front_order(self):
        """All 24 lysis-front codes sort into the expected intelligent order."""
        import random
        shuffled = self.LYSIS_FRONT_CODES[:]
        random.shuffle(shuffled)
        assert smart_sort(shuffled) == self.EXPECTED_ORDER

    def test_roman_position_sorted_numerically(self):
        """Within a fixed prefix, roman numeral position determines group order."""
        codes = ["TB-xi__1_000", "TB-ix__1_000", "TB-xiii__1_000"]
        result = smart_sort(codes)
        assert result == ["TB-ix__1_000", "TB-xi__1_000", "TB-xiii__1_000"]

    def test_integer_position_sorted_numerically(self):
        """Within a fixed prefix and roman, integer position sorts numerically."""
        codes = ["TB-xi__1_582_867", "TB-xi__9_951", "TB-xi__316_573"]
        result = smart_sort(codes)
        assert result == ["TB-xi__9_951", "TB-xi__316_573", "TB-xi__1_582_867"]

    def test_roman_primary_over_integer(self):
        """Roman numeral position is primary; a lower roman value comes first
        even when its integer is larger."""
        # ix (9) with large int vs xi (11) with small int → ix row first
        codes = ["TB-xi__9_951", "TB-ix__1_582_867"]
        assert smart_sort(codes) == ["TB-ix__1_582_867", "TB-xi__9_951"]

    def test_tb_before_tf_lexicographic(self):
        """The first token (TB vs TF) is sorted lexicographically (TB < TF)."""
        codes = ["TF-v__9_951", "TB-xiii__1_582_867"]
        assert smart_sort(codes) == ["TB-xiii__1_582_867", "TF-v__9_951"]

    def test_deterministic_on_repeated_calls(self):
        """smart_sort returns the same order on repeated calls."""
        codes = self.LYSIS_FRONT_CODES[:]
        assert smart_sort(codes) == smart_sort(codes)
