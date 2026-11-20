import numpy as np
import pandas as pd
import pytest
from helpers import GUIDE, exactly, golden_note, golden_rows, sirna_de, sirna_seqs
from hypothesis import given, settings
from hypothesis import strategies as st

import seedmatch
from seedmatch import Guide, SearchOptions, Seed, SequenceSet, find_matches
from seedmatch.errors import InputError
from seedmatch.noncanonical import (
    WobbleBulgeOptions,
    annotate_wobble_bulge_counts,
    wobble_bulge_matches,
    wobble_bulge_patterns,
)
from seedmatch.sequence.alphabet import MASKS

GOLDEN_FILE = "noncanonical/wobble_bulge_r.tsv"
CASES = golden_rows(GOLDEN_FILE)
SEQUENCE_OF = {row.case: row.sequence for row in CASES}
ROW_COLUMNS = ["start", "end", "bulge_after", "bulge_size"]
MATCH_COLUMNS = ["seq_id", "gene_id", "start", "end", "width", "bulge_after", "bulge_size"]
INTEGER_COLUMNS = ["start", "end", "width", "bulge_after", "bulge_size"]

# A row is (start, end, bulge_after, bulge_size). The rows come from a regular-expression oracle of
# the rules, and the hypothesis test below recomputes them.
EXPECTED = {
    "perfect": [(9, 15, 0, 0)],
    "wobble_g4": [(9, 15, 0, 0)],
    "mismatch_t4": [],
    "bulge_after_g5": [(9, 16, 5, 1)],
    "bulge_after_g4": [],
    "bulge_4_after_g7": [(9, 19, 7, 4)],
    "bulge_5_after_g7": [],
    "one_start_two_patterns": [(9, 15, 0, 0)],
    "five_wobbles": [(5, 11, 0, 0)],
}
DEFAULT_PATTERNS = [
    ("CTCTRTR", 0, 0),
    ("CTCNTRTR", 5, 1),
    ("CTNCTRTR", 6, 1),
    ("CNTCTRTR", 7, 1),
    ("CTCNNTRTR", 5, 2),
    ("CTNNCTRTR", 6, 2),
    ("CNNTCTRTR", 7, 2),
    ("CTCNNNTRTR", 5, 3),
    ("CTNNNCTRTR", 6, 3),
    ("CNNNTCTRTR", 7, 3),
    ("CTCNNNNTRTR", 5, 4),
    ("CTNNNNCTRTR", 6, 4),
    ("CNNNNTCTRTR", 7, 4),
]


@pytest.mark.parametrize(
    "name",
    [
        "WobbleBulgeOptions",
        "annotate_wobble_bulge_counts",
        "wobble_bulge_matches",
        "wobble_bulge_patterns",
    ],
)
def test_the_package_exports_the_wobble_and_bulge_search(name):
    assert name in seedmatch.__all__
    assert getattr(seedmatch, name) is getattr(seedmatch.noncanonical, name)


def case_row(case):
    (row,) = (row for row in CASES if row.case == case)
    return row


def matches_of(case, options=WobbleBulgeOptions()):
    row = case_row(case)
    return wobble_bulge_matches(
        Guide(row.guide), SequenceSet.from_records([(case, row.sequence)]), options
    )


def rows_of(frame):
    return list(frame[ROW_COLUMNS].itertuples(index=False, name=None))


def test_the_expected_rows_cover_the_golden_cases():
    assert [row.case for row in CASES] == list(EXPECTED)


def test_wobble_bulge_patterns_of_the_default_guide():
    patterns = wobble_bulge_patterns(Guide(GUIDE))
    assert list(patterns.columns) == ["pattern", "bulge_after", "bulge_size"]
    assert list(patterns.itertuples(index=False, name=None)) == DEFAULT_PATTERNS
    assert patterns[["bulge_after", "bulge_size"]].dtypes.tolist() == [np.int64, np.int64]
    assert isinstance(patterns.index, pd.RangeIndex)


@pytest.mark.parametrize(("max_bulge", "bulge_free"), [(0, 3), (1, 0), (2, 5), (2, 6), (3, 7)])
def test_wobble_bulge_patterns_obey_the_options(max_bulge, bulge_free):
    patterns = wobble_bulge_patterns(Guide(GUIDE), WobbleBulgeOptions(max_bulge, bulge_free))
    assert len(patterns) == 1 + max_bulge * max(0, 6 - bulge_free)
    assert patterns["pattern"].iloc[0] == "CTCTRTR"
    assert patterns[["bulge_after", "bulge_size"]].dtypes.tolist() == [np.int64, np.int64]
    assert isinstance(patterns.index, pd.RangeIndex)
    if max_bulge > 0:
        bulged = patterns.iloc[1:]
        assert set(bulged["bulge_after"]) == set(range(bulge_free + 2, 8))
    assert patterns["pattern"].str.len().tolist() == (7 + patterns["bulge_size"]).tolist()


def test_wobble_bulge_options_store_plain_ints():
    options = WobbleBulgeOptions(max_bulge=np.int64(2), bulge_free=np.int8(1))
    assert (type(options.max_bulge), type(options.bulge_free)) == (int, int)
    assert (options.max_bulge, options.bulge_free) == (2, 1)


@pytest.mark.parametrize("row", CASES, ids=lambda row: row.case)
def test_wobble_bulge_matches_of_the_cases_D005(row):
    """R stops for each call (note `pipeline_error`), thus the rows come from the port rules."""
    assert rows_of(matches_of(row.case)) == EXPECTED[row.case]


def test_each_r_call_stops_D005():
    """The notes `pipeline_error` and `pipeline_error_bulges_lim_1` give the error of R."""
    error = "object 'bulge_allowance' not found"
    assert golden_note(GOLDEN_FILE, r"^# pipeline_error: .* stops with: (.*)$") == error
    assert (
        golden_note(GOLDEN_FILE, r"^# pipeline_error_bulges_lim_1: .* stops with: (.*)$") == error
    )


def test_one_site_counts_once_D005():
    row = case_row("perfect")
    got = annotate_wobble_bulge_counts(
        pd.DataFrame({"gene_id": ["perfect"]}),
        SequenceSet.from_records([("perfect", row.sequence)]),
        Guide(GUIDE),
    )
    assert got["overall_matches"].tolist() == [1]
    assert (row.r_rule_bulges_lim_1, row.r_rule_bulges_lim_4) == ("12", "24")


def test_a_mismatch_is_not_a_bulge_D005():
    assert matches_of("mismatch_t4").empty
    row = case_row("mismatch_t4")
    assert (row.r_rule_bulges_lim_1, row.r_rule_bulges_lim_4) == ("4", "12")


def test_an_extra_target_base_is_a_bulge_D005():
    for case in ("bulge_after_g5", "bulge_4_after_g7"):
        assert len(matches_of(case)) == 1
        row = case_row(case)
        assert (row.r_rule_bulges_lim_1, row.r_rule_bulges_lim_4) == ("0", "0")


def test_no_bulge_is_next_to_g2_to_g4_D005():
    """R drops an N after t1, but it keeps an N after t3 or t5 (note `check_bulges_nt`).

    The N after t1 is between g8 and g7. The N after t3 is between g6 and g5, and the N after t5
    is between g4 and g3. The port keeps the gaps after t1 and t3, and drops the gap after t5,
    which is next to g4.
    """
    assert golden_note(GOLDEN_FILE, r"^# check_bulges_nt: .* keeps (.*)$") == "CUCUANUA,CUCNUAUA"
    assert matches_of("bulge_after_g4").empty
    got = matches_of("bulge_after_g4", WobbleBulgeOptions(bulge_free=2))
    assert rows_of(got) == [(9, 16, 4, 1)]


def test_a_bulge_above_max_bulge_is_not_a_site():
    assert matches_of("bulge_5_after_g7").empty
    got = matches_of("bulge_5_after_g7", WobbleBulgeOptions(max_bulge=5))
    assert rows_of(got) == [(9, 20, 7, 5)]


def test_each_set_of_wobbles_is_in_one_pattern_D005():
    """R stops for fewer than 4 U, and no seed of R has R at 5 positions.

    The notes `wobble_seeds_0_u` and `wobble_seeds_2_u` give the errors of R for a guide with no U
    and with 2 U. The note `wobble_seeds_5_u_in_seed` gives the seeds of R for 5 U in the seed.
    """
    no_u = wobble_bulge_patterns(Guide("CAGAGAGCAAGAACACAGACAAC"))
    assert no_u["pattern"].iloc[0] == "GCTCTCT"
    two_u = wobble_bulge_patterns(Guide("AUAGAGAGCAAGAACACAGAUAA"))
    assert two_u["pattern"].iloc[0] == "CTCTCTR"
    assert rows_of(matches_of("five_wobbles")) == [(5, 11, 0, 0)]
    assert (
        golden_note(GOLDEN_FILE, r"^# wobble_seeds_0_u: .* stops with: (.*)$") == "bad value of n"
    )
    assert (
        golden_note(GOLDEN_FILE, r"^# wobble_seeds_2_u: .* stops with: (.*)$")
        == "v is either non-atomic or too short"
    )
    assert "R" * 5 not in golden_note(GOLDEN_FILE, r"^# wobble_seeds_5_u_in_seed: .* gives (.*)$")


def test_each_distinct_start_counts_once_D005():
    got = wobble_bulge_matches(Guide(GUIDE), sirna_seqs())
    ttr = got[got["seq_id"] == "ENSRNOG00000016275"]
    assert rows_of(ttr) == [(69, 79, 6, 4), (73, 79, 0, 0)]
    assert len(matches_of("one_start_two_patterns")) == 1


def test_max_bulge_0_gives_the_matches_of_the_wobble_target():
    got = wobble_bulge_matches(Guide(GUIDE), sirna_seqs(), WobbleBulgeOptions(max_bulge=0))
    target = Guide(GUIDE).target(Seed.MER7M8, wobble=True)
    unfixed = SearchOptions(fixed_pattern=False, fixed_subject=False)
    expected = find_matches(target, sirna_seqs(), unfixed)
    assert len(expected) == 91
    pd.testing.assert_frame_equal(got[expected.columns], expected)
    assert got[["bulge_after", "bulge_size"]].eq(0).all(axis=None)


def test_each_canonical_match_is_a_row_without_a_bulge():
    canonical = find_matches(Guide(GUIDE).target(Seed.MER7M8), sirna_seqs())
    assert len(canonical) == 20
    got = wobble_bulge_matches(Guide(GUIDE), sirna_seqs())
    rows = set(zip(got["seq_id"], got["start"], got["end"], got["bulge_size"], strict=True))
    for row in canonical.itertuples():
        assert (row.seq_id, row.start, row.end, 0) in rows


def test_the_fixture_gives_the_counts_of_the_oracle():
    got = wobble_bulge_matches(Guide(GUIDE), sirna_seqs())
    assert len(got) == 757
    assert got["seq_id"].nunique() == 145
    assert got["bulge_size"].value_counts().sort_index().tolist() == [91, 189, 162, 159, 156]


def oracle_patterns(guide, max_bulge, bulge_free):
    """Return (pattern, bulge_after, bulge_size) in the order of the rules.

    Target position k pairs with the guide position g(9 - k), and a guide U pairs with A or G.
    """
    pairs = {"A": "T", "C": "G", "G": "C", "U": "R"}
    target = "".join(pairs[guide[8 - k]] for k in range(1, 8))
    patterns = [(target, 0, 0)]
    for size in range(1, max_bulge + 1):
        for after in range(bulge_free + 2, 8):
            cut = 8 - after
            patterns.append((target[:cut] + "N" * size + target[cut:], after, size))
    return patterns


def oracle_rows(guide, records, max_bulge, bulge_free):
    """Return (seq_id, start, end, bulge_after, bulge_size) for each distinct start of a match."""
    patterns = oracle_patterns(guide, max_bulge, bulge_free)
    rows = []
    for seq_id, text in records:
        for start in range(len(text)):
            for pattern, after, size in patterns:
                end = start + len(pattern)
                window = text[start:end]
                if len(window) == len(pattern) and all(
                    MASKS[letter] & MASKS[found]
                    for letter, found in zip(pattern, window, strict=True)
                ):
                    rows.append((seq_id, start + 1, end, after, size))
                    break
    return rows


@st.composite
def oracle_cases(draw):
    guide = "".join(draw(st.lists(st.sampled_from("ACGU"), min_size=8, max_size=23)))
    target = Guide(guide).target(Seed.MER7M8, wobble=True).sequence
    records = []
    for number in range(draw(st.integers(0, 4))):
        text = draw(st.text(alphabet="ACGTNRY", max_size=40))
        if draw(st.booleans()):
            site = "".join(draw(st.sampled_from("AG")) if c == "R" else c for c in target)
            cut = draw(st.integers(0, len(site)))
            site = site[:cut] + draw(st.text(alphabet="ACGT", max_size=4)) + site[cut:]
            at = draw(st.integers(0, len(text)))
            text = text[:at] + site + text[at:]
        records.append((f"s{number}", text))
    return guide, records, draw(st.integers(0, 3)), draw(st.integers(0, 7))


@settings(deadline=None)
@given(oracle_cases())
def test_wobble_bulge_matches_equal_an_oracle(case):
    guide, records, max_bulge, bulge_free = case
    got = wobble_bulge_matches(
        Guide(guide),
        SequenceSet.from_records(records),
        WobbleBulgeOptions(max_bulge, bulge_free),
        workers=1,
    )
    columns = ["seq_id", "start", "end", "bulge_after", "bulge_size"]
    assert list(got[columns].itertuples(index=False, name=None)) == oracle_rows(
        guide, records, max_bulge, bulge_free
    )


def test_an_n_run_matches_each_pattern():
    got = wobble_bulge_matches(Guide(GUIDE), SequenceSet.from_records([("n", "N" * 11)]))
    assert got["start"].tolist() == [1, 2, 3, 4, 5]
    assert got["bulge_size"].tolist() == [0] * 5


@pytest.mark.parametrize(
    "records", [[], [("e", ""), ("s", "CTCTAT")]], ids=["empty set", "short sequences"]
)
def test_empty_and_short_sequences_give_an_empty_frame(records):
    got = wobble_bulge_matches(Guide(GUIDE), SequenceSet.from_records(records))
    assert got.empty
    assert list(got.columns) == MATCH_COLUMNS
    assert got[INTEGER_COLUMNS].dtypes.tolist() == [np.int64] * 5
    assert got[["seq_id", "gene_id"]].dtypes.tolist() == [object, object]
    assert isinstance(got.index, pd.RangeIndex)


def test_the_matches_have_the_columns_the_dtypes_and_the_gene_ids():
    seqs = SequenceSet.from_records(
        [("t1", SEQUENCE_OF["perfect"]), ("t2", SEQUENCE_OF["bulge_after_g5"])],
        gene_ids=["g1", "g2"],
    )
    got = wobble_bulge_matches(Guide(GUIDE), seqs)
    assert list(got.columns) == MATCH_COLUMNS
    assert got[INTEGER_COLUMNS].dtypes.tolist() == [np.int64] * 5
    assert got.to_numpy().tolist() == [["t1", "g1", 9, 15, 7, 0, 0], ["t2", "g2", 9, 16, 8, 5, 1]]
    assert isinstance(got.index, pd.RangeIndex)


def test_a_repeated_sequence_id_keeps_the_rows_and_the_counts_of_each_sequence():
    site = SEQUENCE_OF["perfect"]
    records = [("x", site), ("x", site + "CTCTATAGG")]
    seqs = SequenceSet.from_records(records, gene_ids=["g1", "g2"])
    got = wobble_bulge_matches(Guide(GUIDE), seqs)
    assert got[["seq_id", "gene_id", "start", "end"]].to_numpy().tolist() == [
        ["x", "g1", 9, 15],
        ["x", "g2", 9, 15],
        ["x", "g2", 24, 30],
    ]
    de = pd.DataFrame({"gene_id": ["g1", "g2"]})
    assert annotate_wobble_bulge_counts(de, seqs, Guide(GUIDE))["overall_matches"].tolist() == [
        1,
        2,
    ]
    one_gene = SequenceSet.from_records(records, gene_ids=["g", "g"])
    got = annotate_wobble_bulge_counts(pd.DataFrame({"gene_id": ["g"]}), one_gene, Guide(GUIDE))
    assert got["overall_matches"].tolist() == [2]


def test_annotate_wobble_bulge_counts_give_the_gene_maximum_D003():
    seqs = SequenceSet.from_records(
        [("a", SEQUENCE_OF["perfect"]), ("b", "GGCTCTATAGGCTCTATAGG")], gene_ids=["g", "g"]
    )
    got = annotate_wobble_bulge_counts(pd.DataFrame({"gene_id": ["g"]}), seqs, Guide(GUIDE))
    assert got["overall_matches"].tolist() == [2]


def test_annotate_wobble_bulge_counts_equal_the_rows_of_each_gene():
    matches = wobble_bulge_matches(Guide(GUIDE), sirna_seqs())
    rows = matches.groupby("gene_id").size()
    got = annotate_wobble_bulge_counts(sirna_de(), sirna_seqs(), Guide(GUIDE))
    assert got["overall_matches"].dtype == np.int64
    assert got["overall_matches"].tolist() == [rows.get(gene, 0) for gene in got["gene_id"]]
    assert got["overall_matches"].sum() == len(matches)


def test_annotate_wobble_bulge_counts_name_the_column_and_keep_the_shared_genes_rule():
    seqs = SequenceSet.from_records([("perfect", SEQUENCE_OF["perfect"])])
    de = pd.DataFrame({"gene_id": ["perfect", "absent"], "sites": [9, 9]})
    shared = annotate_wobble_bulge_counts(de, seqs, Guide(GUIDE), column="sites")
    assert shared.to_dict("list") == {"gene_id": ["perfect"], "sites": [1]}
    every = annotate_wobble_bulge_counts(de, seqs, Guide(GUIDE), column="sites", shared_genes=False)
    assert every["sites"].tolist() == [1, 0]


def test_annotate_wobble_bulge_counts_check_the_de_table_before_the_searches(monkeypatch):
    def search(*_args, **_kwargs):
        message = "the function searched before it checked the DE table"
        raise AssertionError(message)

    monkeypatch.setattr("seedmatch.noncanonical.wobble.match_ranges", search)
    seqs = SequenceSet.from_records([("perfect", SEQUENCE_OF["perfect"])])
    message = "DE table has no 'gene_id' column; found ['x']"
    with pytest.raises(InputError, match=exactly(message)):
        annotate_wobble_bulge_counts(pd.DataFrame({"x": [1]}), seqs, Guide(GUIDE))
    message = "column 'gene_id' would overwrite the gene IDs; expected another column name"
    with pytest.raises(InputError, match=exactly(message)):
        annotate_wobble_bulge_counts(
            pd.DataFrame({"gene_id": ["perfect"]}), seqs, Guide(GUIDE), column="gene_id"
        )


def test_wobble_bulge_matches_do_not_depend_on_workers():
    one = wobble_bulge_matches(Guide(GUIDE), sirna_seqs(), workers=1)
    four = wobble_bulge_matches(Guide(GUIDE), sirna_seqs(), workers=4)
    pd.testing.assert_frame_equal(one, four)


def test_wobble_bulge_matches_reject_a_workers_below_1():
    message = "workers 0 is less than 1; expected an int of 1 or more"
    with pytest.raises(InputError, match=exactly(message)):
        wobble_bulge_matches(Guide(GUIDE), sirna_seqs(), workers=0)


def test_wobble_bulge_options_reject_a_negative_max_bulge():
    message = "max_bulge -1 is negative; expected 0 or more"
    with pytest.raises(InputError, match=exactly(message)):
        WobbleBulgeOptions(max_bulge=-1)


@pytest.mark.parametrize(
    ("max_bulge", "letters"), [(19994, 20001), (2**63 - 1, 9223372036854775814)]
)
def test_wobble_bulge_options_reject_a_max_bulge_above_the_pattern_limit(max_bulge, letters):
    message = f"max_bulge {max_bulge} makes a pattern of {letters} letters; expected at most 19993"
    with pytest.raises(InputError, match=exactly(message)):
        WobbleBulgeOptions(max_bulge=max_bulge)
    assert WobbleBulgeOptions(max_bulge=19993).max_bulge == 19993


@pytest.mark.parametrize("bulge_free", [8, -1])
def test_wobble_bulge_options_reject_bulge_free_above_7(bulge_free):
    message = f"bulge_free {bulge_free} is out of range; expected an integer from 0 to 7"
    with pytest.raises(InputError, match=exactly(message)):
        WobbleBulgeOptions(bulge_free=bulge_free)


@pytest.mark.parametrize("field", ["max_bulge", "bulge_free"])
def test_wobble_bulge_options_reject_a_bool(field):
    message = f"{field} True has type bool; expected an integer"
    with pytest.raises(InputError, match=exactly(message)):
        WobbleBulgeOptions(**{field: True})


def test_wobble_bulge_matches_reject_a_guide_of_another_type():
    message = "guide has type str; expected a Guide"
    with pytest.raises(InputError, match=exactly(message)):
        wobble_bulge_matches("UUAUAGAG", sirna_seqs())


def test_wobble_bulge_matches_reject_sequences_of_another_type():
    message = "sequences has type list; expected a SequenceSet"
    with pytest.raises(InputError, match=exactly(message)):
        wobble_bulge_matches(Guide(GUIDE), [("a", "ACGT")])


def test_wobble_bulge_matches_reject_options_of_another_type():
    message = "options has type SearchOptions; expected a WobbleBulgeOptions"
    with pytest.raises(InputError, match=exactly(message)):
        wobble_bulge_matches(Guide(GUIDE), sirna_seqs(), SearchOptions())


def test_the_other_functions_reject_a_guide_and_options_of_another_type():
    de = pd.DataFrame({"gene_id": ["a"]})
    with pytest.raises(InputError, match=exactly("guide has type str; expected a Guide")):
        wobble_bulge_patterns("UUAUAGAG")
    with pytest.raises(InputError, match=exactly("guide has type str; expected a Guide")):
        annotate_wobble_bulge_counts(de, sirna_seqs(), "UUAUAGAG")
    message = "options has type dict; expected a WobbleBulgeOptions"
    with pytest.raises(InputError, match=exactly(message)):
        wobble_bulge_patterns(Guide(GUIDE), {})
    with pytest.raises(InputError, match=exactly(message)):
        annotate_wobble_bulge_counts(de, sirna_seqs(), Guide(GUIDE), options={})
    message = "sequences has type str; expected a SequenceSet"
    with pytest.raises(InputError, match=exactly(message)):
        annotate_wobble_bulge_counts(de, "ACGT", Guide(GUIDE))
