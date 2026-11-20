from dataclasses import FrozenInstanceError

import numpy as np
import pytest
from helpers import GUIDE, sirna_seqs
from hypothesis import given, settings
from hypothesis import strategies as st
from strategies import subjects_for

from seedmatch import Guide, Seed, Target
from seedmatch.errors import InputError
from seedmatch.search import SearchOptions, count_matches, find_matches, match_ranges
from seedmatch.search.compat import compat_table
from seedmatch.sequence import SequenceSet
from seedmatch.sequence.alphabet import MASKS


def s(**named):
    return SequenceSet.from_records(named.items())


# R probes with Biostrings 2.78.0 give the expected values of the tests down to
# `test_overlapping_matches_count_separately`, and of
# `test_indels_with_max_mismatch_at_the_pattern_length_gives_the_mismatch_result` (the R values of
# `EMPTY_INSIDE` in `test_kernels.py`). The other tests follow the search rules or check the shape
# of the result.


def test_mismatch_match_can_extend_past_both_ends():
    got = find_matches("ACGT", s(a="CGTAAAAAAAAACGT", b="TTTTTTACG"), SearchOptions(max_mismatch=1))
    assert got[["seq_id", "start", "end"]].values.tolist() == [
        ["a", 0, 3],
        ["a", 12, 15],
        ["b", 7, 10],
    ]


def test_indel_match_stays_inside_the_sequence():
    got = find_matches("ACGT", s(b="TTTTTTACG"), SearchOptions(max_mismatch=1, indels=True))
    assert got[["start", "end"]].values.tolist() == [[7, 9]]


def test_ranges_obey_the_indels_flag_D001():
    seqs = s(a="TTTTACGACGTTTT")
    assert find_matches("ACGTACG", seqs, SearchOptions(max_mismatch=1)).empty
    with_indels = find_matches("ACGTACG", seqs, SearchOptions(max_mismatch=1, indels=True))
    assert with_indels[["start", "end"]].values.tolist() == [[5, 10]]


def test_ranges_obey_min_mismatch_D001():
    assert find_matches(
        "ACGAACG", s(a="TTTTACGAACGTTTT"), SearchOptions(max_mismatch=1, min_mismatch=1)
    ).empty


def test_n_in_subject_matches_each_base_when_subject_is_not_fixed():
    seqs = s(a="ACNTGG")
    assert count_matches("ACGT", seqs).tolist() == [0]
    unfixed = SearchOptions(fixed_pattern=False, fixed_subject=False)
    assert count_matches("ACGT", seqs, unfixed).tolist() == [1]


@pytest.mark.parametrize("fixed_pattern", [True, False])
def test_a_run_of_n_in_the_subject_matches_each_base_when_subject_is_not_fixed(fixed_pattern):
    seqs = s(a="ANNNNT")
    unfixed = SearchOptions(fixed_pattern=fixed_pattern, fixed_subject=False)
    assert count_matches("ACGT", seqs, unfixed).tolist() == [3]
    assert find_matches("ACGT", seqs, unfixed)[["start", "end"]].values.tolist() == [
        [1, 4],
        [2, 5],
        [3, 6],
    ]
    fixed = SearchOptions(fixed_pattern=fixed_pattern, fixed_subject=True)
    assert count_matches("ACGT", seqs, fixed).tolist() == [0]


def test_overlapping_matches_count_separately():
    assert count_matches("AA", s(a="AAAA")).tolist() == [3]


def test_indels_with_no_mismatch_is_the_exact_search():
    seqs = s(a="ACGTT", b="ACTT")
    options = SearchOptions(indels=True)
    assert count_matches("ACGT", seqs, options).tolist() == [1, 0]
    assert find_matches("ACGT", seqs, options)[["seq_id", "start", "end"]].values.tolist() == [
        ["a", 1, 4]
    ]


def test_indels_with_max_mismatch_at_the_pattern_length_gives_the_mismatch_result():
    seqs = s(a="GG", b="", c="CAC")
    options = SearchOptions(max_mismatch=2, indels=True)
    assert count_matches("AC", seqs, options).tolist() == [3, 1, 4]
    assert find_matches("AC", seqs, options)[["seq_id", "start", "end"]].values.tolist() == [
        ["a", 0, 1],
        ["a", 1, 2],
        ["a", 2, 3],
        ["b", 0, 1],
        ["c", 0, 1],
        ["c", 1, 2],
        ["c", 2, 3],
        ["c", 3, 4],
    ]


def test_counts_have_one_int64_entry_for_each_sequence():
    got = count_matches("AC", s(a="ACAC", b="GG", c=""))
    assert got.dtype == np.int64
    assert got.tolist() == [2, 0, 0]


def test_find_matches_gives_the_ids_the_gene_ids_and_the_width():
    seqs = SequenceSet.from_records(
        [("t1", "ACAC"), ("t2", "GG"), ("t3", "AC")], gene_ids=["g1", "g1", "g2"]
    )
    got = find_matches("AC", seqs)
    assert list(got.columns) == ["seq_id", "gene_id", "start", "end", "width"]
    assert got.values.tolist() == [
        ["t1", "g1", 1, 2, 2],
        ["t1", "g1", 3, 4, 2],
        ["t3", "g2", 1, 2, 2],
    ]
    assert got[["start", "end", "width"]].dtypes.tolist() == [np.int64] * 3


def test_width_of_an_indel_match_is_its_range_width():
    got = find_matches("ACGT", s(b="TTTTTTACG"), SearchOptions(max_mismatch=1, indels=True))
    assert got["width"].tolist() == [3]


@pytest.mark.parametrize("seqs", [s(a="TTTT", b=""), SequenceSet.from_records([])])
def test_no_match_gives_an_empty_frame_with_the_columns_and_dtypes(seqs):
    got = find_matches("ACGT", seqs)
    assert got.empty
    assert list(got.columns) == ["seq_id", "gene_id", "start", "end", "width"]
    assert got[["start", "end", "width"]].dtypes.tolist() == [np.int64] * 3
    assert count_matches("ACGT", seqs).tolist() == [0] * len(seqs)


def test_a_target_is_searched_by_its_sequence():
    target = Target(name="t", sequence="AC", seed_rna="GU")
    assert count_matches(target, s(a="ACAC")).tolist() == [2]
    assert find_matches(target, s(a="ACAC"))["start"].tolist() == [1, 3]


def test_lowercase_pattern_and_subject_match_as_uppercase():
    assert count_matches("acg", s(a="ACGacg")).tolist() == [2]


def test_indels_with_min_mismatch_is_rejected():
    with pytest.raises(InputError, match="min_mismatch must be 0 when indels"):
        SearchOptions(max_mismatch=1, min_mismatch=1, indels=True)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"max_mismatch": -1}, r"max_mismatch -1.*0 or more"),
        ({"min_mismatch": -2}, r"min_mismatch -2.*0 or more"),
        ({"max_mismatch": 1, "min_mismatch": 2}, r"min_mismatch 2.*max_mismatch 1"),
    ],
)
def test_options_reject_a_limit_that_is_not_valid(kwargs, message):
    with pytest.raises(InputError, match=message):
        SearchOptions(**kwargs)


@pytest.mark.parametrize("field", ["max_mismatch", "min_mismatch"])
@pytest.mark.parametrize(
    "bad",
    [1.5, float("nan"), float("inf"), True, "1", None, np.array(1.5), np.array([2]), np.True_],
)
def test_options_reject_a_limit_that_is_not_an_integer(field, bad):
    with pytest.raises(InputError, match=rf"^{field} .*expected an integer"):
        SearchOptions(**{field: bad})


@pytest.mark.parametrize("field", ["max_mismatch", "min_mismatch"])
@pytest.mark.parametrize("big", [2**63, np.uint64(2**64 - 1), 10**30])
def test_options_reject_a_limit_that_the_kernels_cannot_hold(field, big):
    with pytest.raises(InputError, match=rf"^{field} {big} is too large; expected at most 2\^63"):
        SearchOptions(**{field: big})


@pytest.mark.parametrize("field", ["indels", "fixed_pattern", "fixed_subject"])
@pytest.mark.parametrize("bad", ["no", 1, 0, None])
def test_options_reject_a_flag_that_is_not_a_bool(field, bad):
    with pytest.raises(InputError, match=rf"^{field} .*expected a bool"):
        SearchOptions(**{field: bad})


def test_numpy_integers_give_the_result_of_plain_integers_and_are_stored_as_int():
    seqs = s(a="ACGTAC", b="TTACGTTT", c="AAGTAA")
    plain = SearchOptions(max_mismatch=2)
    options = SearchOptions(max_mismatch=np.int32(2), min_mismatch=np.int8(0))
    assert type(options.max_mismatch) is int
    assert type(options.min_mismatch) is int
    assert options == plain
    assert (
        count_matches("ACGT", seqs, options).tolist() == count_matches("ACGT", seqs, plain).tolist()
    )
    assert count_matches("ACGT", seqs, plain).tolist() != count_matches("ACGT", seqs).tolist()
    with_indels = SearchOptions(max_mismatch=np.int64(1), indels=True)
    assert type(with_indels.max_mismatch) is int
    assert (
        count_matches("ACGT", seqs, with_indels).tolist()
        == count_matches("ACGT", seqs, SearchOptions(max_mismatch=1, indels=True)).tolist()
    )


def test_numpy_bools_are_stored_as_bool():
    options = SearchOptions(indels=np.True_, fixed_pattern=np.False_, fixed_subject=np.True_)
    assert (options.indels, options.fixed_pattern, options.fixed_subject) == (True, False, True)
    assert all(
        type(flag) is bool
        for flag in (options.indels, options.fixed_pattern, options.fixed_subject)
    )


def test_options_are_frozen():
    with pytest.raises(FrozenInstanceError):
        SearchOptions().max_mismatch = 1


@pytest.mark.parametrize("pattern", ["", Target(name="t", sequence="", seed_rna="")])
def test_empty_pattern_is_rejected(pattern):
    with pytest.raises(InputError, match="empty pattern"):
        count_matches(pattern, s(a="ACGT"))
    with pytest.raises(InputError, match="empty pattern"):
        find_matches(pattern, s(a="ACGT"))


def test_pattern_with_a_letter_that_is_not_dna_is_rejected():
    with pytest.raises(InputError, match=r"'U' at position 3"):
        count_matches("ACU", s(a="ACGT"))


PATTERN = "ACGTAC"
API_SETTINGS = settings(deadline=None, max_examples=200)

SEQUENCE_SETS = st.lists(subjects_for(PATTERN, 60), max_size=6).map(
    lambda texts: SequenceSet.from_records((f"s{i}", text) for i, text in enumerate(texts))
)


@st.composite
def search_options(draw):
    indels = draw(st.booleans())
    max_mismatch = draw(st.integers(0, 4))
    min_mismatch = 0 if indels else draw(st.integers(0, min(1, max_mismatch)))
    return SearchOptions(
        max_mismatch=max_mismatch,
        min_mismatch=min_mismatch,
        indels=indels,
        fixed_pattern=draw(st.booleans()),
        fixed_subject=draw(st.booleans()),
    )


def distance(pattern, text, start, end, options):
    """The edit distance of the pattern and the range [start, end] of the text, 1-based.

    Without indels, it counts the pattern positions that are outside the text or incompatible.
    With indels, it is the Levenshtein distance of the pattern and the part of the range that is
    inside the text.
    """
    table = compat_table(fixed_pattern=options.fixed_pattern, fixed_subject=options.fixed_subject)

    def compatible(p, t):
        return bool(table[MASKS[p], MASKS[t]])

    if not options.indels:
        return sum(
            not (0 <= start - 1 + a < len(text) and compatible(letter, text[start - 1 + a]))
            for a, letter in enumerate(pattern)
        )
    window = text[max(start, 1) - 1 : end]
    previous = list(range(len(window) + 1))
    for letter in pattern:
        current = [previous[0] + 1]
        for j, base in enumerate(window, start=1):
            substitution = previous[j - 1] + (not compatible(letter, base))
            current.append(min(substitution, previous[j] + 1, current[j - 1] + 1))
        previous = current
    return previous[-1]


def test_distance_helper_on_known_cases():
    plain = SearchOptions(max_mismatch=2)
    assert distance("ACGT", "ACGT", 1, 4, plain) == 0
    assert distance("ACGT", "CGTA", 0, 3, plain) == 1
    assert distance("ACGT", "ACCT", 1, 4, plain) == 1
    assert distance("ACGT", "AC", 1, 4, plain) == 2
    indels = SearchOptions(max_mismatch=2, indels=True)
    assert distance("ACGT", "ACT", 1, 3, indels) == 1
    assert distance("ACGT", "ACCGT", 1, 5, indels) == 1
    assert distance("ACGT", "TTTT", 1, 4, indels) == 3


@API_SETTINGS
@given(SEQUENCE_SETS, search_options())
def test_count_equals_number_of_ranges(seqs, options):
    counts = count_matches(PATTERN, seqs, options)
    ranges = find_matches(PATTERN, seqs, options)
    per_sequence = ranges.groupby("seq_id", sort=False).size().reindex(list(seqs.ids), fill_value=0)
    assert counts.tolist() == per_sequence.tolist()
    order = [seqs.ids.index(seq_id) for seq_id in ranges["seq_id"]]
    assert order == sorted(order)


@API_SETTINGS
@given(SEQUENCE_SETS, search_options())
def test_each_range_is_within_the_edit_limit(seqs, options):
    for row in find_matches(PATTERN, seqs, options).itertuples():
        text = seqs.sequence(seqs.ids.index(row.seq_id))
        assert row.width == row.end - row.start + 1
        found = distance(PATTERN, text, row.start, row.end, options)
        assert options.min_mismatch <= found <= options.max_mismatch


def test_the_pattern_length_limit_of_biostrings():
    # Biostrings 2.78.0 stops above 20000 letters in each search mode, and counts 10001 here.
    subject = s(x="A" * 30000)
    assert count_matches("A" * 20000, subject).tolist() == [10001]
    with pytest.raises(InputError, match="pattern has 20001 letters; expected at most 20000"):
        count_matches("A" * 20001, subject)


@pytest.mark.parametrize(
    "options",
    [
        SearchOptions(),
        SearchOptions(max_mismatch=1),
        SearchOptions(max_mismatch=1, indels=True),
    ],
)
def test_match_ranges_give_the_rows_of_find_matches(options):
    seqs = sirna_seqs()
    target = Guide(GUIDE).target(Seed.MER7M8)
    index, starts, ends = match_ranges(target, seqs, options)
    frame = find_matches(target, seqs, options)
    assert len(frame) > 0
    assert np.array(seqs.ids, dtype=object)[index].tolist() == frame["seq_id"].tolist()
    assert starts.tolist() == frame["start"].tolist()
    assert ends.tolist() == frame["end"].tolist()
    assert [array.dtype for array in (index, starts, ends)] == [np.int64] * 3


def test_match_ranges_name_the_sequence_of_a_repeated_id():
    seqs = SequenceSet.from_records([("x", "GG"), ("x", "AC"), ("y", "ACAC")])
    index, starts, ends = match_ranges("AC", seqs)
    assert index.tolist() == [1, 2, 2]
    assert starts.tolist() == [1, 1, 3]
    assert ends.tolist() == [2, 2, 4]
    assert find_matches("AC", seqs)["seq_id"].tolist() == ["x", "y", "y"]
