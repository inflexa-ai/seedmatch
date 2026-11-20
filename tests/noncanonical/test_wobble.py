import inspect

import numpy as np
import pandas as pd
import pytest
from helpers import GUIDE, exactly, golden_note, golden_rows, sirna_de, sirna_seqs
from hypothesis import given, settings
from hypothesis import strategies as st

import seedmatch
from seedmatch import Guide, SearchOptions, Seed, SequenceSet, e_value_null, e_values, find_matches
from seedmatch import _pool as pool_module
from seedmatch import evalue as evalue_module
from seedmatch.errors import InputError
from seedmatch.evalue import _shuffled_guides
from seedmatch.noncanonical import (
    WobbleBulgeOptions,
    annotate_wobble_bulge_counts,
    wobble_bulge_matches,
    wobble_bulge_null,
    wobble_bulge_patterns,
)
from seedmatch.search import MAX_PATTERN_LENGTH
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
        "wobble_bulge_null",
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


def test_annotate_wobble_bulge_counts_replace_a_column_with_the_name_in_place_D015():
    de = pd.DataFrame(
        {"gene_id": ["perfect", "absent"], "overall_matches": [99, 99], "other": ["x", "y"]}
    )
    seqs = SequenceSet.from_records([("perfect", SEQUENCE_OF["perfect"])])
    got = annotate_wobble_bulge_counts(de, seqs, Guide(GUIDE), shared_genes=False)
    assert got.columns.tolist() == ["gene_id", "overall_matches", "other"]
    assert got["overall_matches"].tolist() == [1, 0]
    assert de["overall_matches"].tolist() == [99, 99]


@pytest.mark.parametrize("shared_genes", [True, False])
@pytest.mark.parametrize(("column", "name"), [(None, "overall_matches"), ("sites", "sites")])
def test_annotate_wobble_bulge_counts_reject_a_repeated_gene_id_or_count_column_D051(
    shared_genes, column, name
):
    seqs = SequenceSet.from_records([("perfect", SEQUENCE_OF["perfect"])])
    arguments = {"column": column} if column else {}
    de = pd.DataFrame([["perfect", 5, 6]], columns=["gene_id", name, name])
    message = f"DE table has 2 columns '{name}'; expected one column for each name"
    with pytest.raises(InputError, match=exactly(message)):
        annotate_wobble_bulge_counts(de, seqs, Guide(GUIDE), shared_genes=shared_genes, **arguments)
    de = pd.DataFrame([["perfect", "x"]], columns=["gene_id", "gene_id"])
    message = "DE table has 2 columns 'gene_id'; expected one column for each name"
    with pytest.raises(InputError, match=exactly(message)):
        annotate_wobble_bulge_counts(de, seqs, Guide(GUIDE), shared_genes=shared_genes, **arguments)


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


@pytest.mark.parametrize("above", [1, 2**63 - 1 - MAX_PATTERN_LENGTH])
def test_wobble_bulge_options_reject_a_max_bulge_above_the_pattern_limit(above):
    largest = MAX_PATTERN_LENGTH - 7
    max_bulge = largest + above
    message = (
        f"max_bulge {max_bulge} makes a pattern of {7 + max_bulge} letters; "
        f"expected at most {largest}"
    )
    with pytest.raises(InputError, match=exactly(message)):
        WobbleBulgeOptions(max_bulge=max_bulge)
    assert WobbleBulgeOptions(max_bulge=largest).max_bulge == largest


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


def multi_sequence_set():
    """Return 21 random sequences that share 7 genes, so that a gene has sequences apart."""
    letters = np.random.default_rng(5).choice(list("ACGT"), size=(21, 1500))
    records = [(f"t{i}", "".join(row)) for i, row in enumerate(letters)]
    return SequenceSet.from_records(records, gene_ids=[f"g{(i * 5) % 7}" for i in range(21)])


def oracle_shuffles(guide, count, seed):
    """Draw the shuffles with the rule of D-008, as a reader of the docstring does."""
    generator = np.random.default_rng(seed)
    letters = guide.sequence
    return [
        Guide("".join(letters[i] for i in generator.permutation(len(letters))))
        for _ in range(count)
    ]


def test_wobble_bulge_null_counts_each_shuffle_like_annotate_wobble_bulge_counts_D008():
    sequences = multi_sequence_set()
    guide = Guide(GUIDE)
    options = WobbleBulgeOptions(max_bulge=2, bulge_free=2)
    null = wobble_bulge_null(guide, sequences, options, rng=7, n_random=5)
    genes = list(dict.fromkeys(sequences.gene_ids))
    table = pd.DataFrame({"gene_id": genes})
    expected = {
        j: annotate_wobble_bulge_counts(
            table, sequences, shuffle, options=options, shared_genes=False
        )["overall_matches"].to_numpy()
        for j, shuffle in enumerate(oracle_shuffles(guide, 5, 7))
    }
    assert null.index.tolist() == genes
    assert null.index.name == "gene_id"
    assert null.columns.tolist() == [0, 1, 2, 3, 4]
    assert (null.dtypes == np.int32).all()
    assert null.to_numpy().max() > 0
    for j, counts in expected.items():
        assert null[j].tolist() == counts.tolist()


def test_wobble_bulge_null_uses_the_default_options_of_the_wobble_search():
    sequences = multi_sequence_set()
    guide = Guide(GUIDE)
    null = wobble_bulge_null(guide, sequences, rng=3, n_random=3)
    genes = list(dict.fromkeys(sequences.gene_ids))
    for j, shuffle in enumerate(oracle_shuffles(guide, 3, 3)):
        counts = annotate_wobble_bulge_counts(
            pd.DataFrame({"gene_id": genes}), sequences, shuffle, shared_genes=False
        )["overall_matches"]
        assert null[j].tolist() == counts.tolist()
    default = inspect.signature(wobble_bulge_null).parameters
    assert default["options"].default == WobbleBulgeOptions()
    assert default["n_random"].default == 1000
    assert default["workers"].default is None


def test_wobble_bulge_null_takes_the_gene_maximum_of_the_sequences_D003():
    sequences = multi_sequence_set()
    guide = Guide(GUIDE)
    null = wobble_bulge_null(guide, sequences, rng=2, n_random=4)
    for j, shuffle in enumerate(oracle_shuffles(guide, 4, 2)):
        rows = wobble_bulge_matches(shuffle, sequences).groupby("seq_id").size()
        by_sequence = pd.Series(
            [rows.get(name, 0) for name in sequences.ids], index=sequences.gene_ids
        )
        assert null[j].to_dict() == by_sequence.groupby(level=0).max().to_dict()


def test_wobble_bulge_null_draws_the_shuffles_of_e_value_null_D008(monkeypatch):
    calls = []
    draw = _shuffled_guides

    def recording_draw(guide, count, generator):
        drawn = []
        calls.append(drawn)
        for shuffle in draw(guide, count, generator):
            drawn.append(shuffle)
            yield shuffle

    monkeypatch.setattr(evalue_module, "_shuffled_guides", recording_draw)
    guide = Guide(GUIDE)
    sequences = sirna_seqs()
    e_value_null(guide, sequences, rng=11, n_random=6)
    wobble_bulge_null(guide, sequences, rng=11, n_random=6)
    assert len(calls) == 2
    assert calls[0] == calls[1] == oracle_shuffles(guide, 6, 11)


def test_wobble_bulge_null_is_the_same_for_the_same_seed():
    sequences = multi_sequence_set()
    guide = Guide(GUIDE)
    first = wobble_bulge_null(guide, sequences, rng=7, n_random=6)
    again = wobble_bulge_null(guide, sequences, rng=7, n_random=6)
    from_generator = wobble_bulge_null(guide, sequences, rng=np.random.default_rng(7), n_random=6)
    other = wobble_bulge_null(guide, sequences, rng=8, n_random=6)
    pd.testing.assert_frame_equal(again, first)
    pd.testing.assert_frame_equal(from_generator, first)
    assert not other.equals(first)


def test_wobble_bulge_null_does_not_depend_on_workers():
    sequences = multi_sequence_set()
    guide = Guide(GUIDE)
    one = wobble_bulge_null(guide, sequences, rng=7, n_random=30, workers=1)
    four = wobble_bulge_null(guide, sequences, rng=7, n_random=30, workers=4)
    pd.testing.assert_frame_equal(four, one)


def test_wobble_bulge_null_runs_the_patterns_of_a_shuffle_in_the_pool_of_the_shuffles(monkeypatch):
    """The pool runs over the shuffles, and the 13 patterns of a shuffle run one after the other."""
    pools = []
    executor = pool_module.ThreadPoolExecutor

    def counting_pool(*args, **kwargs):
        pools.append(kwargs.get("max_workers", args[0] if args else None))
        return executor(*args, **kwargs)

    monkeypatch.setattr(pool_module, "ThreadPoolExecutor", counting_pool)
    wobble_bulge_null(Guide(GUIDE), sirna_seqs(), rng=1, n_random=5, workers=3)
    assert pools == [3]


def test_wobble_bulge_null_of_a_set_without_sequences_has_no_row():
    empty = SequenceSet.from_records([])
    null = wobble_bulge_null(Guide(GUIDE), empty, rng=1, n_random=3)
    assert null.shape == (0, 3)
    assert null.index.name == "gene_id"
    assert (null.dtypes == np.int32).all()


def test_wobble_bulge_null_gives_e_values_from_0_to_1_for_the_wobble_count():
    sequences = sirna_seqs()
    guide = Guide(GUIDE)
    counted = annotate_wobble_bulge_counts(sirna_de(), sequences, guide)
    null = wobble_bulge_null(guide, sequences, rng=0, n_random=20)
    scored = e_values(counted, null, column="overall_matches")
    assert len(scored) == len(counted)
    assert scored["e_value"].between(0, 1).all()
    assert len(set(scored["e_value"])) > 1


def test_wobble_bulge_null_checks_the_types_with_the_messages_of_the_wobble_search():
    with pytest.raises(InputError, match=exactly("guide has type str; expected a Guide")):
        wobble_bulge_null("UUAUAGAG", sirna_seqs(), rng=1)
    message = "sequences has type list; expected a SequenceSet"
    with pytest.raises(InputError, match=exactly(message)):
        wobble_bulge_null(Guide(GUIDE), [("a", "ACGT")], rng=1)
    message = "options has type SearchOptions; expected a WobbleBulgeOptions"
    with pytest.raises(InputError, match=exactly(message)):
        wobble_bulge_null(Guide(GUIDE), sirna_seqs(), SearchOptions(), rng=1)


def test_wobble_bulge_null_checks_rng_n_random_and_workers_with_the_messages_of_e_value_null():
    guide, sequences = Guide(GUIDE), sirna_seqs()
    with pytest.raises(InputError, match=r"^rng 1.5 has type float; expected a numpy\.random"):
        wobble_bulge_null(guide, sequences, rng=1.5)
    with pytest.raises(InputError, match=exactly("n_random 0 is less than 1; expected 1 or more")):
        wobble_bulge_null(guide, sequences, rng=1, n_random=0)
    message = "n_random 1.5 has type float; expected an integer"
    with pytest.raises(InputError, match=exactly(message)):
        wobble_bulge_null(guide, sequences, rng=1, n_random=1.5)
    message = "workers 0 is less than 1; expected an int of 1 or more"
    with pytest.raises(InputError, match=exactly(message)):
        wobble_bulge_null(guide, sequences, rng=1, workers=0)


def test_wobble_bulge_null_raises_before_any_draw_for_a_guide_that_has_no_wobble_target():
    generator = np.random.default_rng(1)
    before = generator.bit_generator.state
    with pytest.raises(InputError, match=r"^mer7m8 needs a guide of at least 8 nt, got 7$"):
        wobble_bulge_null(Guide("UUAUAGA"), sirna_seqs(), rng=generator)
    assert generator.bit_generator.state == before
