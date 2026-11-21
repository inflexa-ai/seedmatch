import math
from itertools import product

import numpy as np
import pandas as pd
import pytest
from helpers import DATA, GUIDE, exactly, golden_note, golden_rows, sirna_seqs
from hypothesis import given, settings
from hypothesis import strategies as st

import seedmatch
from seedmatch import CustomSeed, Guide, Seed, SequenceSet, reverse_complement
from seedmatch.errors import InputError
from seedmatch.noncanonical import CompensatoryOptions, compensatory_sites
from seedmatch.search import MAX_PATTERN_LENGTH

GOLDEN_FILE = "noncanonical/compensatory.tsv"
GOLDEN = golden_rows(GOLDEN_FILE)
CASES = pd.read_csv(DATA / "compensatory_cases.tsv", sep="\t", dtype=str, keep_default_na=False)
COLUMNS = [
    "seq_id",
    "gene_id",
    "seed_start",
    "seed_end",
    "seed_mismatches",
    "region_start",
    "region_end",
    "region_length",
    "guide_start",
    "guide_end",
    "region_target",
    "distance",
    "gc_content",
]
TEXT_COLUMNS = ["seq_id", "gene_id", "region_target"]
INTEGER_COLUMNS = [
    "seed_start",
    "seed_end",
    "seed_mismatches",
    "region_start",
    "region_end",
    "region_length",
    "guide_start",
    "guide_end",
    "distance",
]
ROW_COLUMNS = [
    "seed_start",
    "seed_end",
    "seed_mismatches",
    "region_start",
    "region_end",
    "guide_start",
    "guide_end",
    "region_target",
    "distance",
]
SORT_COLUMNS = ["seed_start", "region_start", "region_end", "guide_start"]
ORACLE_COLUMNS = [
    *INTEGER_COLUMNS[:8],
    "region_target",
    "distance",
    "gc_content",
]

# The port arguments of each R golden case of the test data, from the plan.
R_ARGUMENTS = {
    "vignette": {"seed": Seed.MER7M8, "options": {"first_position": 8}},
    "seed_mm0": {"seed": Seed.MER7M8, "options": {"first_position": 8, "max_seed_mismatch": 0}},
    "min8_max10": {
        "seed": Seed.MER7M8,
        "options": {"first_position": 8, "min_length": 8, "max_length": 10},
    },
    "limit40": {"seed": Seed.MER7M8, "options": {"first_position": 8, "max_distance": 38}},
    "region_mm1": {
        "seed": Seed.MER7M8,
        "options": {"first_position": 8, "min_length": 9, "max_region_mismatch": 1},
    },
    "mer6": {"seed": Seed.MER6, "options": {}},
}


@pytest.mark.parametrize("name", ["CompensatoryOptions", "compensatory_sites"])
def test_the_package_exports_the_compensatory_search(name):
    assert name in seedmatch.__all__
    assert getattr(seedmatch, name) is getattr(seedmatch.noncanonical, name)


def case_inputs(case, **changes):
    """Return the guide, the set, and the options of a case of `compensatory_cases.tsv`.

    The options have the arguments of the R call, with `first_position=8` as R does. A change
    replaces one of them. R `limit = L` is `max_distance = L - 2`.
    """
    part = CASES[CASES["case"] == case]
    first = part.iloc[0]
    sequences = SequenceSet.from_records(zip(part["seq_id"], part["sequence"], strict=True))
    arguments = {
        "min_length": int(first.min_length),
        "max_length": None if first.max_length == "NA" else int(first.max_length),
        "first_position": 8,
        "max_distance": int(first.limit) - 2,
        "max_seed_mismatch": int(first.seed_mismatch),
    }
    return Guide(first.guide), sequences, CompensatoryOptions(**(arguments | changes))


def sites_of(case, **changes):
    guide, sequences, options = case_inputs(case, **changes)
    return compensatory_sites(guide, sequences, options=options)


def rows_of(frame):
    return list(frame[ROW_COLUMNS].itertuples(index=False, name=None))


def golden_of(case):
    return [row for row in GOLDEN if row.case == case]


def r_rows(case):
    """Return the golden rows of a case with the port names, in the port order.

    The R columns `Seed_End` and `Seed_Start` hold the start and the end of the seed match, and
    `Compensatory_End` and `Compensatory_Start` hold the start and the end of the region.
    """
    index = {name: number for number, name in enumerate(sirna_seqs().ids)}
    rows = [
        {
            "index": index[row.TX_ID],
            "seq_id": row.TX_ID,
            "seed_start": int(row.Seed_End),
            "seed_end": int(row.Seed_Start),
            "seed_mismatches": int(row.mismatch_in_seed),
            "region_start": int(row.Compensatory_End),
            "region_end": int(row.Compensatory_Start),
            "region_length": int(row.CompensatoryLength),
            "guide_start": int(row.relative_position_guide_seq_start),
            "guide_end": int(row.relative_position_guide_seq_end),
            "region_target": row.Compensatory_Region_in_Transcript,
            "distance": int(row.Distance),
            "gc_content": float(row.GC_Content),
        }
        for row in golden_of(case)
    ]
    return sorted(rows, key=lambda row: (row["index"], *(row[name] for name in SORT_COLUMNS)))


@pytest.mark.parametrize("case", R_ARGUMENTS)
def test_compensatory_sites_match_r_D004(case):
    """R returns 0 rows for each input (note `original_vignette`), and it stops at line 99.

    The golden file comes from a patched copy of the function. The patch skips the early return and
    keeps the rest of the code. The note `line_99` gives the error of the original code after that
    return, and the golden rows are the R rows of each case.
    """
    original = golden_note(GOLDEN_FILE, r"^# original_vignette: .* gives (.*)$")
    assert original == "0 rows and the message: No compensatory regions found."
    line_99 = golden_note(GOLDEN_FILE, r"^# line_99: .* stops with: (.*)$")
    assert line_99.startswith("unable to find an inherited method for function 'reverseComplement'")
    arguments = R_ARGUMENTS[case]
    got = compensatory_sites(
        Guide(GUIDE),
        sirna_seqs(),
        seed=arguments["seed"],
        options=CompensatoryOptions(**arguments["options"]),
    )
    expected = r_rows(case)
    assert len(got) == len(expected) > 0
    for row, reference in zip(got.to_dict("records"), expected, strict=True):
        assert row["gene_id"] == reference["seq_id"]
        for name in ["seq_id", *INTEGER_COLUMNS, "region_target"]:
            assert row[name] == reference[name], (name, row)
        assert math.isclose(row["gc_content"], reference["gc_content"], rel_tol=1e-9)


def test_the_rows_come_in_the_order_of_the_set():
    guide, options = Guide(GUIDE), CompensatoryOptions(first_position=8)
    got = compensatory_sites(guide, sirna_seqs(), options=options)
    index = {name: number for number, name in enumerate(sirna_seqs().ids)}
    keys = [
        (index[row.seq_id], *(getattr(row, name) for name in SORT_COLUMNS))
        for row in got.itertuples()
    ]
    assert len(keys) == 104
    assert keys == sorted(keys)
    assert isinstance(got.index, pd.RangeIndex)
    site = "GGGGGTTCTTGCGGGGGGGGGGCTCTATAGGGGG"
    reversed_ids = SequenceSet.from_records([("b", site), ("a", site)])
    assert compensatory_sites(guide, reversed_ids)["seq_id"].tolist() == ["b", "a"]


def test_a_set_without_a_site_gives_an_empty_frame():
    """The patched R stops with `undefined columns selected`, and the original R returns 0 rows."""
    stops = golden_note(GOLDEN_FILE, r"^# no_site: .* stops with: (.*)$")
    assert stops == "undefined columns selected"
    got = sites_of("no_site")
    assert got.empty
    assert list(got.columns) == COLUMNS
    assert got[INTEGER_COLUMNS].dtypes.tolist() == [np.int64] * len(INTEGER_COLUMNS)
    assert got[TEXT_COLUMNS].dtypes.tolist() == [object] * 3
    assert got["gc_content"].dtype == np.float64
    assert isinstance(got.index, pd.RangeIndex)


def test_a_seed_match_without_a_region_gives_an_empty_frame():
    got = compensatory_sites(Guide(GUIDE), SequenceSet.from_records([("a", "GGGGCTCTATAGGGG")]))
    assert got.empty
    assert list(got.columns) == COLUMNS
    assert got[INTEGER_COLUMNS].dtypes.tolist() == [np.int64] * len(INTEGER_COLUMNS)
    assert got[TEXT_COLUMNS].dtypes.tolist() == [object] * 3
    assert got["gc_content"].dtype == np.float64


@pytest.mark.parametrize(
    "records",
    [[], [("e", ""), ("s", "CTCTAT"), ("t", "G")]],
    ids=["empty set", "empty and short sequences"],
)
def test_an_empty_set_an_empty_sequence_and_a_short_sequence_give_an_empty_frame(records):
    got = compensatory_sites(Guide(GUIDE), SequenceSet.from_records(records))
    assert got.empty
    assert list(got.columns) == COLUMNS
    assert got[INTEGER_COLUMNS].dtypes.tolist() == [np.int64] * len(INTEGER_COLUMNS)
    assert got[TEXT_COLUMNS].dtypes.tolist() == [object] * 3
    assert got["gc_content"].dtype == np.float64
    assert isinstance(got.index, pd.RangeIndex)


def test_the_columns_the_dtypes_and_the_gene_ids_of_a_site():
    site = "GGGGGTTCTTGCGGGGGGGGGGCTCTATAGGGGG"
    sequences = SequenceSet.from_records([("t1", site)], gene_ids=["gene1"])
    got = compensatory_sites(Guide(GUIDE), sequences, options=CompensatoryOptions(first_position=8))
    assert list(got.columns) == COLUMNS
    assert got[INTEGER_COLUMNS].dtypes.tolist() == [np.int64] * len(INTEGER_COLUMNS)
    assert got[TEXT_COLUMNS].dtypes.tolist() == [object] * 3
    assert got["gc_content"].dtype == np.float64
    assert got.to_dict("records") == [
        {
            "seq_id": "t1",
            "gene_id": "gene1",
            "seed_start": 23,
            "seed_end": 29,
            "seed_mismatches": 0,
            "region_start": 5,
            "region_end": 12,
            "region_length": 8,
            "guide_start": 8,
            "guide_end": 15,
            "region_target": "GTTCTTGC",
            "distance": 10,
            "gc_content": 50.0,
        }
    ]


def test_a_seed_match_past_the_end_counts_the_outside_position():
    (reference,) = golden_of("seed_past_end")
    assert (reference.Seed_End, reference.Seed_Start, reference.mismatch_in_seed) == (
        "23",
        "29",
        "1",
    )
    assert rows_of(sites_of("seed_past_end")) == [(23, 29, 1, 5, 12, 8, 15, "GTTCTTGC", 10)]


def test_each_sequence_is_searched_for_a_repeated_id_D052():
    """R selects the sequences by name, thus it gives 4 rows of the first sequence and none of the
    second (note `repeated_name_r`)."""
    first_site = {("23", "5")}
    assert {
        (row.Seed_End, row.Compensatory_End) for row in golden_of("repeated_name_r")
    } == first_site
    assert len(golden_of("repeated_name_r")) == 4
    got = sites_of("repeated_name_r")
    assert got["seq_id"].tolist() == ["g", "g"]
    assert rows_of(got) == [
        (23, 29, 0, 5, 12, 8, 15, "GTTCTTGC", 10),
        (58, 64, 0, 30, 37, 8, 15, "GTTCTTGC", 20),
    ]


def test_a_region_of_one_seed_match_does_not_drop_a_region_of_another_D053():
    """R keeps the longest region of each transcript and region end, and it ignores the seed match.

    The golden has one row, at `Seed_End` 30. The seed match at 13 has a shorter region with the
    same start, 6. The filter by region end keeps both rows, and the second filter, by region
    start, drops the row of the seed match at 13.
    """
    (reference,) = golden_of("two_seed_matches_r")
    assert reference.Seed_End == "30"
    assert rows_of(sites_of("two_seed_matches_r")) == [
        (13, 19, 0, 6, 11, 10, 15, "GTTCTT", 1),
        (30, 36, 0, 6, 13, 8, 15, "GTTCTTGC", 16),
    ]


def test_a_region_inside_a_longer_region_of_its_seed_match_is_dropped():
    got = compensatory_sites(Guide(GUIDE), sirna_seqs())
    row = got[got["seq_id"] == "ENSRNOG00000016275"]
    assert rows_of(row) == [(73, 79, 0, 59, 71, 10, 22, "AAACAGTGTTCTT", 1)]
    assert row["region_length"].tolist() == [13]


def test_each_row_gives_the_guide_positions_of_its_region_D054():
    """R joins each match of the region in the guide with a comma, and it keeps one row for each
    region match. The golden has the texts `8,12` and `13,17`, and `2,8` and `7,13`.
    """
    repeated = golden_of("repeated_region_r")
    assert [row.relative_position_guide_seq_start for row in repeated] == ["8,12"] * 2
    assert [row.relative_position_guide_seq_end for row in repeated] == ["13,17"] * 2
    assert rows_of(sites_of("repeated_region_r")) == [
        (24, 30, 0, 7, 12, 8, 13, "GCTTGC", 11),
        (24, 30, 0, 7, 12, 12, 17, "GCTTGC", 11),
    ]
    (at_g2,) = golden_of("region_at_g2_r")
    assert (at_g2.relative_position_guide_seq_start, at_g2.relative_position_guide_seq_end) == (
        "2,8",
        "7,13",
    )
    assert rows_of(sites_of("region_at_g2_r")) == [(24, 30, 0, 7, 12, 8, 13, "TCTATA", 11)]


def test_seed_mismatches_is_the_number_of_mismatches_D054():
    (reference,) = golden_of("seed_mm2_r")
    assert reference.mismatch_in_seed == "1"
    got = sites_of("seed_mm2_r")
    assert got["seed_mismatches"].tolist() == [2]
    assert rows_of(got) == [(23, 29, 2, 5, 12, 8, 15, "GTTCTTGC", 10)]


def test_the_region_starts_after_the_seed_D056():
    """R always starts the region at g8, and 10 of the 104 rows of the vignette start there."""
    vignette = golden_of("vignette")
    assert len(vignette) == 104
    assert sum(row.relative_position_guide_seq_start == "8" for row in vignette) == 10
    after_seed = compensatory_sites(Guide(GUIDE), sirna_seqs())
    assert len(after_seed) == 96
    assert after_seed["gene_id"].nunique() == 51
    assert after_seed["guide_start"].min() == 9
    from_g8 = compensatory_sites(
        Guide(GUIDE), sirna_seqs(), options=CompensatoryOptions(first_position=8)
    )
    assert len(from_g8) == 104
    assert rows_of(sites_of("repeated_region_r", first_position=None)) == [
        (24, 30, 0, 7, 12, 12, 17, "GCTTGC", 11)
    ]


def test_the_default_region_start_follows_the_seed_window_D056():
    mer6 = compensatory_sites(Guide(GUIDE), sirna_seqs(), seed=Seed.MER6)
    from_g8 = compensatory_sites(
        Guide(GUIDE),
        sirna_seqs(),
        seed=Seed.MER6,
        options=CompensatoryOptions(first_position=8),
    )
    pd.testing.assert_frame_equal(mer6, from_g8)
    assert mer6["guide_start"].min() == 8


def test_compensatory_options_reject_lengths_that_select_no_region_D055():
    """R gives 3 rows of length 6 for `min_length` 6 and `max_length` 3 (case `min_above_max_r`).

    A `min_length` of 17 searches the region g8 to g23 three times (note `min_17`). R with
    `limit = 2` gives no row.
    """
    assert [row.CompensatoryLength for row in golden_of("min_above_max_r")] == ["6"] * 3
    note = golden_note(GOLDEN_FILE, r"^# min_17: .* are (.*)$")
    assert note == ",".join(["GCAAGAACACUGUUUU"] * 3)
    message = "max_length 3 is less than min_length 6; expected None or an integer of 6 or more"
    with pytest.raises(InputError, match=exactly(message)):
        CompensatoryOptions(min_length=6, max_length=3)
    with pytest.raises(InputError, match=exactly(message)):
        case_inputs("min_above_max_r")
    guide, sequences = Guide(GUIDE), sirna_seqs()
    message = "min_length 17 is more than the 15 positions from g9 to g23; expected at most 15"
    with pytest.raises(InputError, match=exactly(message)):
        compensatory_sites(guide, sequences, options=CompensatoryOptions(min_length=17))
    message = "max_length 20 is more than the 15 positions from g9 to g23; expected at most 15"
    with pytest.raises(InputError, match=exactly(message)):
        compensatory_sites(guide, sequences, options=CompensatoryOptions(max_length=20))
    message = "max_distance 0 is less than 1; expected 1 or more"
    with pytest.raises(InputError, match=exactly(message)):
        CompensatoryOptions(max_distance=0)


def test_the_message_names_one_position_in_the_singular_D055():
    one = "min_length 2 is more than the 1 position from g23 to g23; expected at most 1"
    with pytest.raises(InputError, match=exactly(one)):
        compensatory_sites(
            Guide(GUIDE),
            sirna_seqs(),
            options=CompensatoryOptions(min_length=2, first_position=23),
        )
    two = "max_length 3 is more than the 2 positions from g22 to g23; expected at most 2"
    with pytest.raises(InputError, match=exactly(two)):
        compensatory_sites(
            Guide(GUIDE),
            sirna_seqs(),
            options=CompensatoryOptions(min_length=1, max_length=3, first_position=22),
        )


def test_a_window_above_the_pattern_limit_raises_before_a_search_D055(monkeypatch):
    def search(*_args, **_kwargs):
        message = "the function searched before it checked the longest region"
        raise AssertionError(message)

    monkeypatch.setattr("seedmatch.noncanonical.compensatory.match_ranges", search)
    guide = Guide("A" * (MAX_PATTERN_LENGTH + 9))
    message = (
        f"the longest region has {MAX_PATTERN_LENGTH + 1} letters; "
        f"expected at most {MAX_PATTERN_LENGTH}, the pattern limit of the search"
    )
    with pytest.raises(InputError, match=exactly(message)):
        compensatory_sites(guide, sirna_seqs())
    options = CompensatoryOptions(max_length=MAX_PATTERN_LENGTH + 1)
    with pytest.raises(InputError, match=exactly(message)):
        compensatory_sites(guide, sirna_seqs(), options=options)


def test_a_window_at_the_pattern_limit_is_searched():
    guide = Guide("A" * (MAX_PATTERN_LENGTH + 8))
    options = CompensatoryOptions(min_length=MAX_PATTERN_LENGTH, max_length=MAX_PATTERN_LENGTH)
    assert compensatory_sites(guide, sirna_seqs(), options=options).empty


def test_a_short_guide_has_no_room_for_the_region_D055():
    message = "min_length 6 is more than the 5 positions from g9 to g13; expected at most 5"
    with pytest.raises(InputError, match=exactly(message)):
        compensatory_sites(Guide("UUAUAGAGCAAGA"), sirna_seqs())


def test_a_region_can_fill_each_position_after_the_start():
    long = SequenceSet.from_records([("a", "GG" + "AAAACAGTGTTCTTG" + "G" * 10 + "CTCTATA" + "GG")])
    options = CompensatoryOptions(min_length=15, max_length=15)
    got = compensatory_sites(Guide(GUIDE), long, options=options)
    assert rows_of(got) == [(28, 34, 0, 3, 17, 9, 23, "AAAACAGTGTTCTTG", 10)]
    short = SequenceSet.from_records([("a", "GG" + "TCTTG" + "G" * 10 + "CTCTATA")])
    options = CompensatoryOptions(min_length=5, max_length=5)
    got = compensatory_sites(Guide("UUAUAGAGCAAGA"), short, options=options)
    assert rows_of(got) == [(18, 24, 0, 3, 7, 9, 13, "TCTTG", 10)]


def test_the_default_region_starts_after_a_custom_seed_D056():
    sequences = SequenceSet.from_records([("a", "TGTTCT" + "G" * 10 + "TGCTCTATA")])
    got = compensatory_sites(Guide(GUIDE), sequences, seed=CustomSeed("g2_g10", 2, 10))
    assert rows_of(got) == [(17, 25, 0, 1, 6, 11, 16, "TGTTCT", 10)]


def test_a_seed_that_ends_at_the_guide_end_leaves_no_region():
    message = (
        "the region starts at g24, but the guide ends at g23; "
        "expected a seed or a first_position that leaves a region"
    )
    with pytest.raises(InputError, match=exactly(message)):
        compensatory_sites(Guide(GUIDE), sirna_seqs(), seed=Seed.FULL)


def test_compensatory_options_reject_a_first_position_below_1():
    message = "first_position 0 is less than 1; expected 1 or more"
    with pytest.raises(InputError, match=exactly(message)):
        CompensatoryOptions(first_position=0)


def test_compensatory_sites_reject_a_first_position_after_the_guide():
    message = (
        "the region starts at g24, but the guide ends at g23; "
        "expected a seed or a first_position that leaves a region"
    )
    with pytest.raises(InputError, match=exactly(message)):
        compensatory_sites(
            Guide(GUIDE), sirna_seqs(), options=CompensatoryOptions(first_position=24)
        )


def test_the_last_position_of_the_guide_is_a_first_position_of_a_region_of_one_base():
    options = CompensatoryOptions(min_length=1, first_position=23)
    sequences = SequenceSet.from_records([("a", "GGGGG" + "A" + "G" * 18 + "CTCTATA" + "GGGGG")])
    got = compensatory_sites(Guide(GUIDE), sequences, options=options)
    assert rows_of(got) == [(25, 31, 0, 6, 6, 23, 23, "A", 18)]


@pytest.mark.parametrize(
    ("field", "value", "kind"),
    [
        ("min_length", True, "bool"),
        ("min_length", "6", "str"),
        ("min_length", 6.0, "float"),
        ("max_length", "9", "str"),
        ("first_position", 8.5, "float"),
        ("max_distance", None, "NoneType"),
        ("max_seed_mismatch", False, "bool"),
        ("max_region_mismatch", "0", "str"),
    ],
)
def test_compensatory_options_reject_a_value_that_is_not_an_integer(field, value, kind):
    message = f"{field} {value!r} has type {kind}; expected an integer"
    with pytest.raises(InputError, match=exactly(message)):
        CompensatoryOptions(**{field: value})


def test_compensatory_options_store_plain_ints():
    options = CompensatoryOptions(
        min_length=np.int64(7),
        max_length=np.int8(9),
        first_position=np.int32(8),
        max_distance=np.int64(40),
        max_seed_mismatch=np.int8(2),
        max_region_mismatch=np.int8(1),
    )
    for name in (
        "min_length",
        "max_length",
        "first_position",
        "max_distance",
        "max_seed_mismatch",
        "max_region_mismatch",
    ):
        assert type(getattr(options, name)) is int
    assert (options.min_length, options.max_length, options.first_position) == (7, 9, 8)
    assert CompensatoryOptions().max_length is None
    assert CompensatoryOptions().first_position is None


@pytest.mark.parametrize("field", ["min_length", "first_position", "max_distance"])
def test_compensatory_options_reject_a_value_below_1(field):
    message = f"{field} -3 is less than 1; expected 1 or more"
    with pytest.raises(InputError, match=exactly(message)):
        CompensatoryOptions(**{field: -3})


@pytest.mark.parametrize("field", ["max_seed_mismatch", "max_region_mismatch"])
def test_compensatory_options_reject_a_negative_mismatch(field):
    message = f"{field} -1 is negative; expected 0 or more"
    with pytest.raises(InputError, match=exactly(message)):
        CompensatoryOptions(**{field: -1})


def test_a_max_distance_above_the_buffer_does_not_overflow():
    huge = compensatory_sites(
        Guide(GUIDE), sirna_seqs(), options=CompensatoryOptions(max_distance=2**63 - 1)
    )
    long = compensatory_sites(
        Guide(GUIDE), sirna_seqs(), options=CompensatoryOptions(max_distance=10**8)
    )
    pd.testing.assert_frame_equal(huge, long)
    assert huge["distance"].max() > 98


def oracle_matches(pattern, text, limit):
    """Return (start, end, mismatches) of each match of the pattern with at most `limit` mismatches.

    A match of a pattern of length m starts at each position from 2 - m to the length of the text.
    A position outside the text, and a letter that is not equal, is a mismatch. The letter N of a
    text is a mismatch too.
    """
    found = []
    for start in range(2 - len(pattern), len(text) + 1):
        mismatches = sum(
            not 1 <= start + offset <= len(text) or text[start + offset - 1] != letter
            for offset, letter in enumerate(pattern)
        )
        if mismatches <= limit:
            found.append((start, start + len(pattern) - 1, mismatches))
    return found


def oracle_pairs(guide_text, text, options):
    """Return a dict for each pair of a seed match and a region match of one text."""
    seed_target = Guide(guide_text).target(Seed.MER7M8).sequence
    seeds = oracle_matches(seed_target, text, options.max_seed_mismatch)
    first = 9
    longest = len(guide_text) - first + 1 if options.max_length is None else options.max_length
    pairs = []
    for size in range(options.min_length, longest + 1):
        for a in range(first, len(guide_text) - size + 2):
            pattern = reverse_complement(guide_text[a - 1 : a - 1 + size].replace("U", "T"))
            regions = oracle_matches(pattern, text, options.max_region_mismatch)
            for (region_start, region_end, _), (seed_start, seed_end, mismatches) in product(
                regions, seeds
            ):
                distance = seed_start - region_end - 1
                if 1 <= distance <= options.max_distance:
                    percent = (pattern.count("G") + pattern.count("C")) / size * 100
                    row = [seed_start, seed_end, mismatches, region_start, region_end, size]
                    row += [a, a + size - 1, pattern, distance, percent]
                    pairs.append(dict(zip(ORACLE_COLUMNS, row, strict=True)))
    return pairs


def oracle_rows(guide_text, records, options):
    """Return the rows of the rules, one pair at a time, with the port column names and `index`.

    A pair is dropped if a pair of the same seed match has a region that contains its region and is
    longer.
    """
    rows = []
    for index, (_name, text) in enumerate(records):
        pairs = oracle_pairs(guide_text, text, options)
        rows += [
            {"index": index, **pair}
            for pair in pairs
            if not any(
                other["seed_start"] == pair["seed_start"]
                and other["region_start"] <= pair["region_start"]
                and pair["region_end"] <= other["region_end"]
                and other["region_length"] > pair["region_length"]
                for other in pairs
            )
        ]
    return sorted(rows, key=lambda row: (row["index"], *(row[name] for name in SORT_COLUMNS)))


@st.composite
def oracle_cases(draw):
    guide = "".join(draw(st.lists(st.sampled_from("ACGU"), min_size=10, max_size=16)))
    seed_target = Guide(guide).target(Seed.MER7M8).sequence
    windows = [
        reverse_complement(guide[first - 1 : last].replace("U", "T"))
        for first in range(9, len(guide) + 1)
        for last in range(first + 1, len(guide) + 1)
    ]
    pieces = st.one_of(
        st.text(alphabet="ACGTN", max_size=6),
        st.just(seed_target),
        st.builds(
            lambda text, at: text[:at] + "N" + text[at + 1 :],
            st.just(seed_target),
            st.integers(0, 6),
        ),
        st.builds(
            lambda text, at: text[:at] + "A" + text[at + 1 :],
            st.just(seed_target),
            st.integers(0, 6),
        ),
        st.sampled_from(windows),
    )
    records = []
    for number in range(draw(st.integers(1, 3))):
        text = "".join(draw(st.lists(pieces, max_size=7)))[:50]
        records.append((draw(st.sampled_from(["a", "b"])) if number else "a", text))
    min_length = draw(st.integers(2, 4))
    max_length = draw(st.one_of(st.none(), st.integers(min_length, min_length + 3)))
    options = CompensatoryOptions(
        min_length=min_length,
        max_length=max_length,
        max_distance=draw(st.integers(1, 20)),
        max_seed_mismatch=draw(st.integers(0, 2)),
        max_region_mismatch=draw(st.integers(0, 1)),
    )
    return guide, records, options


@settings(deadline=None)
@given(oracle_cases())
def test_compensatory_sites_equal_an_oracle(case):
    guide, records, options = case
    sequences = SequenceSet.from_records(records)
    room = len(guide) - 9 + 1
    if options.min_length > room or (options.max_length or 0) > room:
        with pytest.raises(InputError, match="positions from g9 to"):
            compensatory_sites(Guide(guide), sequences, options=options, workers=1)
        return
    got = compensatory_sites(Guide(guide), sequences, options=options, workers=1)
    expected = oracle_rows(guide, records, options)
    assert len(got) == len(expected)
    for row, reference in zip(got.to_dict("records"), expected, strict=True):
        assert row["seq_id"] == records[reference["index"]][0]
        for name in [*INTEGER_COLUMNS, "region_target"]:
            assert row[name] == reference[name], (name, row, reference)
        assert math.isclose(row["gc_content"], reference["gc_content"], rel_tol=1e-9)


def test_a_window_stops_at_the_end_of_its_sequence():
    """The seed of `b` follows the region of `a` in the buffer of the set, at the distance 1.

    A window that ran into `b` would pair them. The region g9 to g15 ends at the position 11 of `a`.
    """
    sequences = SequenceSet.from_records([("a", "GGGGGTTCTTGC"), ("b", "CTCTATAGGGG")])
    assert compensatory_sites(Guide(GUIDE), sequences).empty
    alone = SequenceSet.from_records([("a", "GGGGGTTCTTGCCTCTATAGGGG")])
    assert rows_of(compensatory_sites(Guide(GUIDE), alone)) == [
        (13, 19, 0, 5, 11, 9, 15, "GTTCTTG", 1)
    ]


def test_a_seed_match_before_the_start_pairs_with_no_region():
    """The seed match of `b` starts at 0, thus its key lies in the last base of `a`."""
    sequences = SequenceSet.from_records([("a", "GGGGGTTCTTGCGGGGG"), ("b", "TCTATAGGGG")])
    assert compensatory_sites(Guide(GUIDE), sequences).empty


def test_a_max_seed_mismatch_above_the_seed_length_adds_no_match():
    site = "GGGGGTTCTTGCGGGGGGGGGGCTCTATAGGGGG"
    sequences = SequenceSet.from_records([("a", site), ("b", site[::-1])])
    seven = compensatory_sites(
        Guide(GUIDE), sequences, options=CompensatoryOptions(max_seed_mismatch=7)
    )
    huge = compensatory_sites(
        Guide(GUIDE), sequences, options=CompensatoryOptions(max_seed_mismatch=2**63 - 1)
    )
    pd.testing.assert_frame_equal(seven, huge)
    assert seven["seed_mismatches"].value_counts().sort_index().to_dict() == {
        0: 1,
        3: 2,
        6: 3,
        7: 21,
    }


def test_a_seed_match_at_the_last_base_pairs_with_a_region_before_it():
    """The seed match at 14 starts at the last base of the sequence, with 6 outside positions."""
    sequences = SequenceSet.from_records([("a", "GGGGGTTCTTGGGC")])
    got = compensatory_sites(
        Guide(GUIDE), sequences, options=CompensatoryOptions(max_seed_mismatch=6)
    )
    assert rows_of(got) == [
        (12, 18, 6, 5, 10, 10, 15, "GTTCTT", 1),
        (14, 20, 6, 5, 11, 9, 15, "GTTCTTG", 2),
    ]


def test_an_n_in_the_seed_site_is_a_mismatch():
    sequences = SequenceSet.from_records([("a", "GGGGGTTCTTGCGGGGGGGGGGCTCNATAGGGGG")])
    got = compensatory_sites(Guide(GUIDE), sequences)
    assert got["seed_mismatches"].tolist() == [1]
    assert len(got) == 1
    none = compensatory_sites(
        Guide(GUIDE), sequences, options=CompensatoryOptions(max_seed_mismatch=0)
    )
    assert none.empty


def test_an_n_in_the_region_is_a_mismatch():
    sequences = SequenceSet.from_records([("a", "GGGGGTTCNTGCGGGGGGGGGGCTCTATAGGGGG")])
    assert compensatory_sites(Guide(GUIDE), sequences).empty
    got = compensatory_sites(
        Guide(GUIDE), sequences, options=CompensatoryOptions(max_region_mismatch=1)
    )
    assert rows_of(got) == [
        (23, 29, 0, 3, 8, 12, 17, "GTGTTC", 14),
        (23, 29, 0, 5, 11, 9, 15, "GTTCTTG", 11),
    ]


def test_a_region_past_the_start_counts_the_outside_position():
    sequences = SequenceSet.from_records([("a", "TTCTTGCGGGGGGGGGGCTCTATAGGGGG")])
    options = CompensatoryOptions(min_length=7, max_length=7, max_region_mismatch=1)
    got = compensatory_sites(Guide(GUIDE), sequences, options=options)
    assert rows_of(got) == [(18, 24, 0, 0, 6, 9, 15, "GTTCTTG", 11)]


def test_compensatory_sites_take_the_options_in_position_like_wobble_bulge_matches():
    options = CompensatoryOptions(min_length=8)
    by_position = compensatory_sites(Guide(GUIDE), sirna_seqs(), options)
    pd.testing.assert_frame_equal(
        by_position, compensatory_sites(Guide(GUIDE), sirna_seqs(), options=options)
    )
    assert len(by_position) < len(compensatory_sites(Guide(GUIDE), sirna_seqs()))
    with pytest.raises(TypeError, match="positional"):
        compensatory_sites(Guide(GUIDE), sirna_seqs(), options, Seed.MER6)


def test_compensatory_sites_do_not_depend_on_workers():
    one = compensatory_sites(Guide(GUIDE), sirna_seqs(), workers=1)
    four = compensatory_sites(Guide(GUIDE), sirna_seqs(), workers=4)
    pd.testing.assert_frame_equal(one, four)
    assert len(one) == 96


def test_compensatory_sites_reject_a_workers_below_1():
    message = "workers 0 is less than 1; expected 1 or more"
    with pytest.raises(InputError, match=exactly(message)):
        compensatory_sites(Guide(GUIDE), sirna_seqs(), workers=0)


def test_compensatory_sites_reject_a_guide_of_another_type():
    message = "guide has type str; expected a Guide"
    with pytest.raises(InputError, match=exactly(message)):
        compensatory_sites(GUIDE, sirna_seqs())


def test_compensatory_sites_reject_sequences_of_another_type():
    message = "sequences has type list; expected a SequenceSet"
    with pytest.raises(InputError, match=exactly(message)):
        compensatory_sites(Guide(GUIDE), [("a", "ACGT")])


def test_compensatory_sites_reject_options_of_another_type():
    message = "options has type dict; expected a CompensatoryOptions"
    with pytest.raises(InputError, match=exactly(message)):
        compensatory_sites(Guide(GUIDE), sirna_seqs(), options={"min_length": 6})


def test_compensatory_sites_reject_a_seed_of_another_type():
    message = "seed 'mer7m8' has type str; expected a Seed or a CustomSeed"
    with pytest.raises(InputError, match=exactly(message)):
        compensatory_sites(Guide(GUIDE), sirna_seqs(), seed="mer7m8")
