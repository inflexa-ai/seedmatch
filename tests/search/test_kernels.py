import re
from itertools import product

import corpus
import numpy as np
import oracle
import pytest
from corpus import Range
from hypothesis import assume, given, settings
from hypothesis import strategies as st
from strategies import dna, subjects_for

from seedmatch.search._kernels import hamming_counts, hamming_ranges, indel_counts, indel_ranges
from seedmatch.search.compat import compat_table
from seedmatch.sequence import SequenceSet
from seedmatch.sequence.alphabet import encode_dna

type Fixed = tuple[bool, bool]
type Hit = tuple[int, int, int]

FIXED: Fixed = (True, True)


def encoded(pattern, subjects, fixed):
    seqs = SequenceSet.from_records((f"s{i}", subject) for i, subject in enumerate(subjects))
    compat = compat_table(fixed_pattern=fixed[0], fixed_subject=fixed[1])
    return seqs.codes, seqs.offsets, encode_dna(pattern), compat


def hits(index, starts, ends) -> list[Hit]:
    assert index.dtype == starts.dtype == ends.dtype == np.int64
    return list(zip(index.tolist(), starts.tolist(), ends.tolist(), strict=True))


def hamming_hits(pattern, subjects, max_mm, min_mm, fixed=FIXED) -> list[Hit]:
    index, starts = hamming_ranges(*encoded(pattern, subjects, fixed), max_mm, min_mm)
    return hits(index, starts, starts + len(pattern) - 1)


def indel_hits(pattern, subjects, max_mm, fixed=FIXED) -> list[Hit]:
    return hits(*indel_ranges(*encoded(pattern, subjects, fixed), max_mm))


def hamming_count_list(pattern, subjects, max_mm, min_mm, fixed=FIXED) -> list[int]:
    counts = hamming_counts(*encoded(pattern, subjects, fixed), max_mm, min_mm)
    assert counts.dtype == np.int64
    return counts.tolist()


def indel_count_list(pattern, subjects, max_mm, fixed=FIXED) -> list[int]:
    counts = indel_counts(*encoded(pattern, subjects, fixed), max_mm)
    assert counts.dtype == np.int64
    return counts.tolist()


def kernel_hamming(pattern, subject, max_mm, min_mm, fixed) -> list[Range]:
    found = hamming_hits(pattern, [subject], max_mm, min_mm, fixed)
    return [(start, end) for _, start, end in found]


def kernel_indels(pattern, subject, max_mm, fixed) -> list[Range]:
    return [(start, end) for _, start, end in indel_hits(pattern, [subject], max_mm, fixed)]


PATTERNS = dna("ACGTN" * 4 + "MRY", min_size=1, max_size=12)
FIXED_PAIRS = st.tuples(st.booleans(), st.booleans())
CASES = PATTERNS.flatmap(
    lambda pattern: st.tuples(
        st.just(pattern), subjects_for(pattern, 60), st.integers(0, 3), FIXED_PAIRS
    )
)
SETS = PATTERNS.flatmap(
    lambda pattern: st.tuples(st.just(pattern), st.lists(subjects_for(pattern, 30), max_size=6))
)

KERNEL_SETTINGS = settings(max_examples=300)


@KERNEL_SETTINGS
@given(CASES, st.integers(0, 1))
def test_hamming_kernel_equals_oracle(case, min_mm):
    pattern, subject, max_mm, fixed = case
    assume(min_mm <= max_mm)
    expected = oracle.hamming(pattern, subject, max_mm, min_mm, *fixed)
    assert kernel_hamming(pattern, subject, max_mm, min_mm, fixed) == expected
    assert hamming_count_list(pattern, [subject], max_mm, min_mm, fixed) == [len(expected)]


@KERNEL_SETTINGS
@given(CASES)
def test_indel_kernel_equals_oracle(case):
    pattern, subject, max_mm, fixed = case
    expected = oracle.indels(pattern, subject, max_mm, *fixed)
    assert kernel_indels(pattern, subject, max_mm, fixed) == expected
    assert indel_count_list(pattern, [subject], max_mm, fixed) == [len(expected)]


@KERNEL_SETTINGS
@given(SETS, st.integers(0, 3), st.integers(0, 1), FIXED_PAIRS)
def test_ranges_of_a_set_equal_the_oracle_for_each_sequence(sets, max_mm, min_mm, fixed):
    pattern, subjects = sets
    assume(min_mm <= max_mm)
    assert hamming_hits(pattern, subjects, max_mm, min_mm, fixed) == [
        (i, *found)
        for i, subject in enumerate(subjects)
        for found in oracle.hamming(pattern, subject, max_mm, min_mm, *fixed)
    ]
    assert indel_hits(pattern, subjects, max_mm, fixed) == [
        (i, *found)
        for i, subject in enumerate(subjects)
        for found in oracle.indels(pattern, subject, max_mm, *fixed)
    ]


@KERNEL_SETTINGS
@given(SETS, st.integers(0, 3), st.integers(0, 1), FIXED_PAIRS)
def test_counts_equal_number_of_ranges_for_each_sequence(sets, max_mm, min_mm, fixed):
    pattern, subjects = sets
    assume(min_mm <= max_mm)
    args = encoded(pattern, subjects, fixed)
    index, _ = hamming_ranges(*args, max_mm, min_mm)
    expected = np.bincount(index, minlength=len(subjects))
    np.testing.assert_array_equal(hamming_counts(*args, max_mm, min_mm), expected)
    index, _, _ = indel_ranges(*args, max_mm)
    expected = np.bincount(index, minlength=len(subjects))
    np.testing.assert_array_equal(indel_counts(*args, max_mm), expected)


def test_empty_set_gives_empty_arrays():
    assert hamming_count_list("AC", [], 2, 0) == []
    assert hamming_hits("AC", [], 2, 0) == []
    # With max_mm 2, the indel kernel takes the mismatch path; with max_mm 1, the indel path.
    for max_mm in (1, 2):
        assert indel_count_list("AC", [], max_mm) == []
        assert indel_hits("AC", [], max_mm) == []


# The expected values of the tests below come from R probes with Biostrings 2.78.0.
EMPTY_INSIDE = ["GG", "", "CAC"]
EMPTY_INSIDE_MISMATCH_HITS = [
    (0, 0, 1),
    (0, 1, 2),
    (0, 2, 3),
    (1, 0, 1),
    (2, 0, 1),
    (2, 1, 2),
    (2, 2, 3),
    (2, 3, 4),
]


def test_empty_sequence_inside_a_set_keeps_its_mismatch_match():
    assert hamming_hits("AC", EMPTY_INSIDE, 2, 0) == EMPTY_INSIDE_MISMATCH_HITS
    assert hamming_count_list("AC", EMPTY_INSIDE, 2, 0) == [3, 1, 4]


def test_indel_kernel_gives_the_mismatch_result_if_max_mm_reaches_the_pattern_length():
    assert indel_hits("AC", EMPTY_INSIDE, 2) == EMPTY_INSIDE_MISMATCH_HITS
    assert indel_count_list("AC", EMPTY_INSIDE, 2) == [3, 1, 4]


def test_empty_sequence_inside_a_set_has_no_indel_match():
    assert indel_hits("AC", EMPTY_INSIDE, 1) == [(2, 1, 1), (2, 2, 3)]
    assert indel_count_list("AC", EMPTY_INSIDE, 1) == [0, 0, 2]


def test_mismatch_match_can_extend_past_both_ends_of_a_short_sequence():
    assert hamming_hits("ACGTACG", ["GTA"], 4, 0) == [(0, -1, 5)]
    assert hamming_count_list("ACGTACG", ["GTA"], 4, 0) == [1]
    assert indel_hits("ACG", ["G"], 3) == [(0, -1, 1), (0, 0, 2), (0, 1, 3)]
    assert indel_count_list("ACG", ["G"], 3) == [3]


def test_indel_match_of_a_short_sequence_stays_inside_it():
    assert indel_hits("ACGTACG", ["GTA", "", "G"], 4) == [(0, 1, 3)]
    assert indel_count_list("ACGTACG", ["GTA", "", "G"], 4) == [1, 0, 0]


def corpus_kernel_mismatches(name: str) -> list[tuple[object, ...]]:
    """Run the kernels once for each case and option set, with the max_mm of the option set."""
    r_counts = corpus.counts(name)
    compat = {
        fixed: compat_table(fixed_pattern=fixed[0], fixed_subject=fixed[1])
        for fixed in product((True, False), repeat=2)
    }
    arrays: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    mismatches = []
    for case, opts, expected in corpus.groups(name):
        if case.case_id not in arrays:
            seqs = SequenceSet.from_records([(case.case_id, case.subject)])
            arrays[case.case_id] = (seqs.codes, seqs.offsets, encode_dna(case.pattern))
        args = (*arrays[case.case_id], compat[opts.fixed_pattern, opts.fixed_subject])
        if opts.indels:
            counts = indel_counts(*args, opts.max_mm)
            _, starts, ends = indel_ranges(*args, opts.max_mm)
        else:
            counts = hamming_counts(*args, opts.max_mm, opts.min_mm)
            _, starts = hamming_ranges(*args, opts.max_mm, opts.min_mm)
            ends = starts + len(case.pattern) - 1
        found = list(zip(starts.tolist(), ends.tolist(), strict=True))
        r_count = r_counts[case.case_id, opts]
        if found != expected or counts.tolist() != [r_count]:
            report = (case.case_id, case.pattern, case.subject, opts, expected, found)
            mismatches.append((*report, r_count, counts.tolist()))
    return mismatches


@pytest.mark.parametrize("name", [corpus.RANDOM, corpus.SEEDS])
def test_kernels_equal_biostrings_on_corpus(name):
    assert corpus_kernel_mismatches(name) == []


LONG_LENGTHS = (63, 64, 65, 80)
ALL_FIXED = list(product((True, False), repeat=2))


def long_case(m, size):
    """Return a pattern of m letters, and a subject with an exact copy and a copy with 2 edits."""
    rng = np.random.default_rng(m)

    def letters(alphabet, count):
        return "".join(rng.choice(list(alphabet), size=count))

    pattern = letters("ACGTRYN", m)
    edited = list(pattern)
    for edit in rng.choice(["substitute", "insert", "delete"], size=2):
        at = int(rng.integers(len(edited)))
        if edit == "substitute":
            edited[at] = letters("ACGT", 1)
        elif edit == "insert":
            edited.insert(at, letters("ACGT", 1))
        else:
            del edited[at]
    filler = size - m - len(edited)
    first, second = sorted(int(cut) for cut in rng.integers(0, filler + 1, size=2))
    subject = (
        letters("ACGT", first)
        + pattern
        + letters("ACGT", second - first)
        + "".join(edited)
        + letters("ACGT", filler - second)
    )
    return pattern, subject


@pytest.mark.parametrize("m", LONG_LENGTHS)
@pytest.mark.parametrize("max_mm", [0, 1, 2, 3])
@pytest.mark.parametrize("fixed", ALL_FIXED)
def test_long_patterns_equal_the_oracle(m, max_mm, fixed):
    pattern, subject = long_case(m, 4 * m)
    expected = oracle.hamming(pattern, subject, max_mm, 0, *fixed)
    assert expected
    assert kernel_hamming(pattern, subject, max_mm, 0, fixed) == expected
    assert kernel_indels(pattern, subject, max_mm, fixed) == oracle.indels(
        pattern, subject, max_mm, *fixed
    )


@pytest.mark.parametrize("max_mm", [62, 63])
@pytest.mark.parametrize("fixed", ALL_FIXED)
def test_pattern_of_64_letters_with_max_mm_near_its_length_equals_the_oracle(max_mm, fixed):
    pattern, subject = long_case(64, 160)
    expected = oracle.indels(pattern, subject, max_mm, *fixed)
    assert expected
    assert kernel_indels(pattern, subject, max_mm, fixed) == expected


KERNEL_CALLS = [
    lambda c, o, p, t: hamming_counts(c, o, p, t, 1, 0),
    lambda c, o, p, t: hamming_ranges(c, o, p, t, 1, 0),
    lambda c, o, p, t: indel_counts(c, o, p, t, 1),
    lambda c, o, p, t: indel_ranges(c, o, p, t, 1),
]
COMPAT = compat_table(fixed_pattern=True, fixed_subject=True)


@pytest.mark.parametrize("call", KERNEL_CALLS)
def test_kernels_reject_offsets_outside_codes(call):
    message = "offsets do not fit codes of length 2; expected the offsets of a SequenceSet"
    with pytest.raises(ValueError, match="^" + re.escape(message) + "$"):
        call(np.array([1, 2], np.uint8), np.array([0, 3], np.int64), encode_dna("AC"), COMPAT)


@pytest.mark.parametrize("call", KERNEL_CALLS)
def test_kernels_read_a_code_above_15_as_its_low_four_bits(call):
    offsets = np.array([0, 3], np.int64)
    high = call(np.array([200, 1, 8], np.uint8), offsets, encode_dna("TA"), COMPAT)
    low = call(np.array([8, 1, 8], np.uint8), offsets, encode_dna("TA"), COMPAT)
    for found, expected in zip(np.atleast_2d(high), np.atleast_2d(low), strict=True):
        np.testing.assert_array_equal(found, expected)


INT64_MAX = 2**63 - 1


def test_limits_near_the_int64_maximum():
    assert hamming_count_list("AC", EMPTY_INSIDE, INT64_MAX, 0) == [3, 1, 4]
    assert hamming_count_list("AC", EMPTY_INSIDE, INT64_MAX, INT64_MAX) == [0, 0, 0]
    assert indel_hits("AC", EMPTY_INSIDE, INT64_MAX) == EMPTY_INSIDE_MISMATCH_HITS
