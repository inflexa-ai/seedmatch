from itertools import product

import corpus
import numpy as np
import oracle
import pytest
from corpus import Range
from hypothesis import assume, given, settings
from hypothesis import strategies as st
from numba import config
from strategies import dna, subjects_for

from seedmatch.search.compat import compat_table
from seedmatch.search.kernels import hamming_counts, hamming_ranges, indel_counts, indel_ranges
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

# The first example compiles the kernels, which takes longer than the default deadline.
KERNEL_SETTINGS = settings(deadline=None, max_examples=300)


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


# Without the JIT (the coverage run), the kernels run as Python code, and the complete corpora take
# about a minute. Each 13th group still compares the Python code with R.
CORPUS_STRIDE = 13 if config.DISABLE_JIT else 1


def corpus_kernel_mismatches(name: str) -> list[tuple[object, ...]]:
    """Run the kernels once for each case and option set, with the max_mm of the option set."""
    r_counts = corpus.counts(name)
    compat = {
        fixed: compat_table(fixed_pattern=fixed[0], fixed_subject=fixed[1])
        for fixed in product((True, False), repeat=2)
    }
    arrays: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    mismatches = []
    for case, opts, expected in corpus.groups(name)[::CORPUS_STRIDE]:
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
