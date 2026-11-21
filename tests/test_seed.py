import dataclasses
import re

import pytest
from helpers import GUIDE, golden_rows

import seedmatch
from seedmatch import CustomSeed, Guide, Seed, SeedSpec, Target
from seedmatch.errors import InputError
from seedmatch.seed import seed_window

MIN_GUIDE_LENGTH = {
    Seed.MER7M8: 8,
    Seed.MER8: 8,
    Seed.MER6: 8,
    Seed.MER7A1: 8,
    Seed.FULL: 8,
    Seed.MER18: 19,
    Seed.MER15: 16,
}
SEED_NAMES = ["mer7m8", "mer8", "mer6", "mer7A1", "Full", "18mer", "15mer"]


def row_id(row):
    window = f"-{row.start}-{row.stop}" if row.start else ""
    return f"{row.guide[:8]}-{row.seed}{window}-wobble{row.wobble}"


@pytest.mark.parametrize("row", golden_rows("seed/targets.tsv"), ids=row_id)
def test_target_matches_r(row):
    seed = CustomSeed("custom", int(row.start), int(row.stop)) if row.start else Seed(row.seed)
    t = Guide(row.guide).target(seed, wobble=row.wobble == "TRUE")
    assert (t.sequence, t.seed_rna) == (row.target, row.seed_rna)


def test_seed_values_are_the_r_seed_names():
    assert [s.value for s in Seed] == SEED_NAMES


@pytest.mark.parametrize("value", ["bogus", "MER8", "full", 8, None])
def test_seed_rejects_a_value_that_is_not_a_seed_name(value):
    message = (
        f"seed {value!r} is not the name of a canonical seed; "
        f"expected one of {', '.join(SEED_NAMES)}"
    )
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        Seed(value)


@pytest.mark.parametrize("seed", list(Seed))
def test_d072_custom_seed_rejects_the_name_of_a_canonical_seed(seed):
    """D-072: R uses the canonical seed and ignores the window, but the port stops."""
    message = (
        f"custom seed name {seed.value!r} is the name of a canonical seed; "
        f"expected another name, or Seed.{seed.name}"
    )
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        CustomSeed(seed.value, 2, 5)


def test_custom_seed_accepts_a_canonical_name_in_another_case():
    """R compares the seed name with the canonical names with the case."""
    target = Guide(GUIDE).target(CustomSeed("MER8", 2, 8))
    assert (target.name, target.sequence) == ("MER8", "CTCTATA")


def test_top_level_package_exports_the_seed_types():
    assert (seedmatch.Guide, seedmatch.Seed, seedmatch.CustomSeed, seedmatch.Target) == (
        Guide,
        Seed,
        CustomSeed,
        Target,
    )
    assert SeedSpec == Seed | CustomSeed


def test_target_name_is_the_seed_value_or_the_custom_name():
    guide = Guide(GUIDE)
    assert guide.target(Seed.MER7A1).name == "mer7A1"
    assert guide.target(CustomSeed("g2-g9", 2, 9)).name == "g2-g9"


def test_guide_stores_uppercase_letters():
    assert Guide("uuauagagcaagaacacuguuuu").sequence == GUIDE
    assert Guide("UuAuAgAgCaAgAaCaCuGuUuU") == Guide(GUIDE)


def test_guide_is_immutable():
    guide = Guide(GUIDE)
    with pytest.raises(dataclasses.FrozenInstanceError):
        guide.sequence = "AAAAAAAA"


@pytest.mark.parametrize("thymine", ["T", "t"])
def test_guide_rejects_thymine_and_names_uracil(thymine):
    with pytest.raises(InputError, match=f"'{thymine}' at position 3.*use 'U'"):
        Guide(f"UU{thymine}UAGAGCAAGAACACUGUUUU")


def test_guide_rejects_thymine_at_position_one():
    with pytest.raises(InputError, match=r"'T' at position 1.*use 'U'"):
        Guide("TTAUAGAGCAAGAACACUGUUUU")


@pytest.mark.parametrize(
    ("text", "letter", "position"),
    [
        ("UUAUAG GCAAGAACACUGUUUU", " ", 7),
        ("UUAUAGAGCAAGAACACUGUUU-", "-", 23),
        ("UUßUAGAGCAAGAACACUGUUUU", "ß", 3),
    ],
)
def test_guide_rejects_a_letter_that_is_not_a_u_c_g_a_and_names_it(text, letter, position):
    pattern = f"{re.escape(repr(letter))} at position {position}.*expected.*A, C, G, or U"
    with pytest.raises(InputError, match=pattern):
        Guide(text)


@pytest.mark.parametrize("seed", list(Seed))
def test_canonical_seed_needs_its_minimum_guide_length(seed):
    minimum = MIN_GUIDE_LENGTH[seed]
    guide = GUIDE[:minimum]
    assert Guide(guide).target(seed).sequence
    message = f"{seed.value} needs a guide of at least {minimum} nt, got {minimum - 1}"
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        Guide(guide[:-1]).target(seed)


def test_canonical_seed_needs_eight_nucleotides():
    with pytest.raises(InputError, match="at least 8"):
        Guide("UUAUAGA").target(Seed.MER7M8)


def test_eighteen_mer_message_names_the_seed_and_both_lengths():
    with pytest.raises(InputError, match=r"^18mer needs a guide of at least 19 nt, got 18$"):
        Guide("UUAUAGAGCAAGAACACU").target(Seed.MER18)


def test_custom_window_must_lie_inside_the_guide():
    with pytest.raises(InputError, match=r"stop 24 .* length 23"):
        Guide("UUAUAGAGCAAGAACACUGUUUU").target(CustomSeed("x", 2, 24))


@pytest.mark.parametrize("letter", list("MRWSYKVHDBNmrwsykvhdbn"))
def test_d019_guide_rejects_each_ambiguity_letter(letter):
    """D-019: R accepts an IUPAC ambiguity letter in a guide; the port raises InputError."""
    pattern = f"{re.escape(repr(letter))} at position 3.*expected A, C, G, or U"
    with pytest.raises(InputError, match=pattern):
        Guide(f"UU{letter}UAGAGCAAGAACACUGUUUU")


@pytest.mark.parametrize("start", [0, -1])
def test_custom_window_start_is_one_based(start):
    with pytest.raises(InputError, match=rf"start {start} .*1"):
        Guide(GUIDE).target(CustomSeed("x", start, 8))


@pytest.mark.parametrize(("start", "stop"), [(8, 2), (24, 23)])
def test_d020_custom_window_start_must_not_pass_stop(start, stop):
    """D-020: R returns an empty target for start after stop; the port raises InputError."""
    with pytest.raises(InputError, match=rf"start {start} .*stop {stop}"):
        Guide(GUIDE).target(CustomSeed("x", start, stop))


def test_custom_window_may_hold_one_nucleotide_and_may_span_the_guide():
    guide = Guide(GUIDE)
    assert guide.target(CustomSeed("one", 1, 1)).sequence == "A"
    assert guide.target(CustomSeed("all", 1, 23)).seed_rna == GUIDE


def test_wobble_leaves_the_fixed_a_of_mer8_and_mer7a1():
    guide = Guide("UUUUUUUUUUUUUUUUUUUUUUU")
    assert guide.target(Seed.MER8, wobble=True).sequence == "RRRRRRRA"
    assert guide.target(Seed.MER7A1, wobble=True).sequence == "RRRRRRA"


def test_d018_lowercase_guide_gets_the_wobble_target_of_the_uppercase_guide():
    """D-018: R ignores the wobble for a lowercase guide; the port applies it."""
    lower = Guide("uuauagagcaagaacacuguuuu").target(Seed.MER7M8, wobble=True)
    assert lower.sequence == "CTCTRTR"
    assert lower == Guide(GUIDE).target(Seed.MER7M8, wobble=True)


@pytest.mark.parametrize("seed", ["mer7m8", "Mer7m8"])
def test_target_rejects_a_seed_name_as_text(seed):
    message = f"seed {seed!r} has type str; expected a Seed or a CustomSeed"
    with pytest.raises(InputError, match=re.escape(message)):
        Guide(GUIDE).target(seed)


def test_target_rejects_a_wobble_that_is_not_a_bool():
    with pytest.raises(InputError, match=re.escape("wobble 'no' has type str; expected a bool")):
        Guide(GUIDE).target(Seed.MER7M8, wobble="no")


@pytest.mark.parametrize(
    ("start", "stop", "message"),
    [
        (2.5, 8, "start 2.5 has type float; expected an integer"),
        (True, 8, "start True has type bool; expected an integer"),
        (2, "8", "stop '8' has type str; expected an integer"),
    ],
)
def test_custom_seed_rejects_a_window_position_that_is_not_an_integer(start, stop, message):
    with pytest.raises(InputError, match=re.escape(message)):
        CustomSeed("x", start, stop)


def test_seed_window_gives_g1_to_the_a_of_mer8_and_mer7a1():
    assert seed_window(Seed.MER8, 23) == (1, 8)
    assert seed_window(Seed.MER7A1, 23) == (1, 7)


@pytest.mark.parametrize(
    ("seed", "length", "window"),
    [
        (Seed.MER7M8, 23, (2, 8)),
        (Seed.MER6, 23, (2, 7)),
        (Seed.FULL, 23, (2, 23)),
        (Seed.FULL, 30, (2, 30)),
        (Seed.MER18, 23, (2, 19)),
        (Seed.MER15, 23, (2, 16)),
        (CustomSeed("custom", 3, 12), 23, (3, 12)),
    ],
)
def test_seed_window_gives_the_guide_positions_of_the_target(seed, length, window):
    assert seed_window(seed, length) == window
