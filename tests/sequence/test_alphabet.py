import re

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

from seedmatch.errors import InputError
from seedmatch.sequence import reverse_complement
from seedmatch.sequence.alphabet import MASKS, decode_dna, encode_dna

IUPAC = "ACGTMRWSYKVHDBN"
COMPLEMENT_PAIRS = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "M": "K",
    "K": "M",
    "R": "Y",
    "Y": "R",
    "V": "B",
    "B": "V",
    "H": "D",
    "D": "H",
    "W": "W",
    "S": "S",
    "N": "N",
}


def test_masks_hold_the_four_bases_and_the_union_of_each_ambiguity_code():
    a, c, g, t = 1, 2, 4, 8
    assert dict(MASKS) == {
        "A": a,
        "C": c,
        "G": g,
        "T": t,
        "M": a | c,
        "R": a | g,
        "W": a | t,
        "S": c | g,
        "Y": c | t,
        "K": g | t,
        "V": a | c | g,
        "H": a | c | t,
        "D": a | g | t,
        "B": c | g | t,
        "N": 15,
    }


def test_masks_cannot_be_changed():
    with pytest.raises(TypeError):
        MASKS["A"] = 0


def test_reverse_complement_of_wobble_seed_matches_biostrings():
    assert reverse_complement("YYAYAGAG") == "CTCTRTRR"


@given(st.text(alphabet=IUPAC))
def test_reverse_complement_is_an_involution(s):
    assert reverse_complement(reverse_complement(s)) == s


@pytest.mark.parametrize(("letter", "expected"), COMPLEMENT_PAIRS.items())
def test_reverse_complement_of_each_letter_swaps_it_with_its_pair(letter, expected):
    assert reverse_complement(letter) == expected


@pytest.mark.parametrize("letter", IUPAC)
def test_complement_of_a_letter_covers_the_complement_of_each_base_of_the_letter(letter):
    bases = [b for b in "ACGT" if MASKS[b] & MASKS[letter]]
    expected = sum(MASKS[COMPLEMENT_PAIRS[b]] for b in bases)
    assert MASKS[reverse_complement(letter)] == expected


def test_reverse_complement_changes_lowercase_to_uppercase():
    assert reverse_complement("acgtn") == "NACGT"


def test_reverse_complement_of_the_empty_string_is_empty():
    assert reverse_complement("") == ""


def test_reverse_complement_rejects_a_letter_that_is_not_dna():
    with pytest.raises(InputError, match=r"'U' at position 3"):
        reverse_complement("ACUG")


def test_lowercase_bases_encode_as_uppercase():
    assert np.array_equal(encode_dna("acgtn"), encode_dna("ACGTN"))


def test_encode_gives_the_mask_of_each_letter_as_uint8():
    codes = encode_dna(IUPAC)
    assert codes.dtype == np.uint8
    assert codes.tolist() == [MASKS[letter] for letter in IUPAC]


def test_encode_of_the_empty_string_is_an_empty_uint8_array():
    codes = encode_dna("")
    assert codes.dtype == np.uint8
    assert codes.shape == (0,)


def test_encoded_array_can_be_changed_by_the_caller():
    codes = encode_dna("ACGT")
    codes[0] = 15
    assert codes.tolist() == [15, 2, 4, 8]


def test_uracil_in_dna_is_rejected_with_its_position():
    with pytest.raises(InputError, match=r"'U' at position 3"):
        encode_dna("ACUG")


@pytest.mark.parametrize(
    ("text", "letter", "position"),
    [
        ("acgu", "'u'", 4),
        ("AC-G", "'-'", 3),
        ("A.CG", "'.'", 2),
        ("AC GT", "' '", 3),
        ("ACG\n", r"'\\n'", 4),
        ("ACG\tT", r"'\\t'", 4),
        ("XACG", "'X'", 1),
        ("ACGé", "'é'", 4),
        ("ACG\U0001f600T", "'\U0001f600'", 4),
        ("ACé\U0001f600U", "'é'", 3),
    ],
)
def test_encode_rejects_each_byte_that_is_not_an_iupac_dna_letter(text, letter, position):
    with pytest.raises(InputError, match=rf"{letter} at position {position}\b"):
        encode_dna(text)


@pytest.mark.parametrize("letter", ["-", "+", "."])
def test_d017_gap_mask_and_unknown_letters_are_an_error(letter):
    """D-017: Biostrings keeps these letters in a DNAString, but the port stops with InputError."""
    text = f"AC{letter}GT"
    expected = rf"invalid DNA letter '{re.escape(letter)}' at position 3"
    with pytest.raises(InputError, match=expected):
        encode_dna(text)
    with pytest.raises(InputError, match=expected):
        reverse_complement(text)


def test_error_message_names_the_expected_letters():
    with pytest.raises(InputError, match=IUPAC):
        encode_dna("U")


@given(st.text(alphabet=IUPAC + IUPAC.lower()))
def test_decode_gives_back_the_uppercase_text(s):
    assert decode_dna(encode_dna(s)) == s.upper()


def test_decode_of_an_empty_array_is_the_empty_string():
    assert decode_dna(np.zeros(0, dtype=np.uint8)) == ""


@pytest.mark.parametrize("value", [0, 16, 255])
def test_decode_rejects_a_value_that_is_not_a_mask(value):
    codes = np.array([1, 2, value], dtype=np.uint8)
    with pytest.raises(InputError, match=rf"mask {value} at position 3"):
        decode_dna(codes)
