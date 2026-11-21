from itertools import product
from typing import NamedTuple

import corpus
import numpy as np
import oracle
import pytest
from corpus import Group, Options, Range
from helpers import exactly, golden_rows

from seedmatch.search.compat import compat_table
from seedmatch.sequence.alphabet import MASKS

FIXED_PAIRS = list(product((True, False), repeat=2))
MAX_MM = 4


def option_sets(max_mms: range, indels_values: tuple[bool, ...]) -> list[Options]:
    return [
        Options(max_mm, min_mm, indels, fixed_pattern, fixed_subject)
        for max_mm, indels, (fixed_pattern, fixed_subject) in product(
            max_mms, indels_values, FIXED_PAIRS
        )
        for min_mm in ((0,) if indels else (0, 1))
        if min_mm <= max_mm
    ]


OPTION_SETS = {
    corpus.RANDOM: option_sets(range(MAX_MM + 1), (False, True)),
    corpus.SEEDS: option_sets(range(1, MAX_MM + 1), (True,)),
}


def mismatch_report(group: Group, found: list[Range]) -> tuple[object, ...]:
    case, opts, expected = group
    return case.case_id, case.pattern, case.subject, opts, expected, found


def test_compat_table_matches_biostrings():
    expected = {pair: np.zeros((16, 16), dtype=np.bool_) for pair in FIXED_PAIRS}
    seen = set()
    for row in golden_rows("search/iupac_compat.tsv"):
        pair = (row.fixed_pattern == "TRUE", row.fixed_subject == "TRUE")
        expected[pair][MASKS[row.pattern], MASKS[row.subject]] = row.match == "TRUE"
        seen.add((pair, row.pattern, row.subject))
    assert len(seen) == len(FIXED_PAIRS) * len(MASKS) ** 2
    for pair in FIXED_PAIRS:
        table = compat_table(fixed_pattern=pair[0], fixed_subject=pair[1])
        assert table.dtype == np.bool_
        np.testing.assert_array_equal(table, expected[pair], err_msg=str(pair))


@pytest.mark.parametrize("name", list(OPTION_SETS))
def test_corpus_has_each_case_with_each_valid_option_set(name):
    assert set(corpus.counts(name)) == set(product(corpus.cases(name), OPTION_SETS[name]))


@pytest.mark.parametrize("name", list(OPTION_SETS))
def test_biostrings_counts_equal_its_range_counts(name):
    counts, ranges = corpus.counts(name), corpus.ranges(name)
    assert set(ranges) <= set(counts)
    assert [key for key, count in counts.items() if count != len(ranges.get(key, []))] == []


def test_hamming_oracle_matches_biostrings_on_corpus():
    mismatches = []
    for group in corpus.groups(corpus.RANDOM):
        case, opts, expected = group
        if opts.indels:
            continue
        found = oracle.hamming(
            case.pattern,
            case.subject,
            opts.max_mm,
            opts.min_mm,
            opts.fixed_pattern,
            opts.fixed_subject,
        )
        if found != expected:
            mismatches.append(mismatch_report(group, found))
    assert mismatches == []


@pytest.mark.parametrize("name", list(OPTION_SETS))
def test_indel_oracle_matches_biostrings_on_corpus(name):
    mismatches = []
    for group in corpus.groups(name):
        case, opts, expected = group
        if not opts.indels:
            continue
        fixed = (opts.fixed_pattern, opts.fixed_subject)
        found = oracle.indels(case.pattern, case.subject, opts.max_mm, *fixed)
        if found != expected:
            mismatches.append(mismatch_report(group, found))
    assert mismatches == []


def test_biostrings_ranges_hold_each_documented_best_local_match():
    """Biostrings does not obey its "best local match" definition, but it reports each such match.

    The R ranges hold the documented matches in the same order, for 0 < max_mm < len(pattern).
    """
    tables: dict[tuple[str, bool, bool], oracle.EditTable] = {}
    mismatches = []
    for group in corpus.groups(corpus.RANDOM):
        case, opts, expected = group
        if not opts.indels or not 0 < opts.max_mm < len(case.pattern):
            continue
        fixed = (opts.fixed_pattern, opts.fixed_subject)
        key = (case.case_id, *fixed)
        if key not in tables:
            tables[key] = oracle.edit_distances(case.pattern, case.subject, MAX_MM, *fixed)
        documented = oracle.documented_best_local_matches(tables[key], opts.max_mm)
        if [r for r in expected if r in documented] != documented:
            mismatches.append(mismatch_report(group, documented))
    assert mismatches == []


class Probe(NamedTuple):
    search: str
    pattern: str
    subject: str
    max_mm: int
    min_mm: int
    fixed: tuple[bool, bool]
    r: list[Range]


FIXED = (True, True)
AMBIGUOUS = (False, False)

# Each value comes from an R probe with Biostrings 2.78.0.
PROBES = [
    Probe("hamming", "AA", "AAAA", 0, 0, FIXED, [(1, 2), (2, 3), (3, 4)]),
    Probe("hamming", "ACGT", "CGTAAAAAAAAACGT", 1, 0, FIXED, [(0, 3), (12, 15)]),
    Probe("hamming", "ACGT", "TTTTTTACG", 1, 0, FIXED, [(7, 10)]),
    Probe("hamming", "ACGTACG", "TTTTACGACGTTTT", 1, 0, FIXED, []),
    Probe("hamming", "ACGT", "ACNTGG", 0, 0, AMBIGUOUS, [(1, 4)]),
    Probe("hamming", "ACGT", "ACNTGG", 0, 0, FIXED, []),
    Probe("hamming", "AC", "ACAC", 2, 2, FIXED, [(0, 1), (2, 3), (4, 5)]),
    Probe("hamming", "AC", "GG", 3, 0, FIXED, [(0, 1), (1, 2), (2, 3)]),
    Probe("hamming", "AC", "", 2, 0, FIXED, [(0, 1)]),
    Probe("hamming", "ACGT", "", 4, 0, FIXED, [(-2, 1), (-1, 2), (0, 3)]),
    Probe("hamming", "A", "", 1, 0, FIXED, []),
    Probe("indels", "ACGT", "TTTTTTACG", 1, 0, FIXED, [(7, 9)]),
    Probe("indels", "ACGTACG", "TTTTACGACGTTTT", 1, 0, FIXED, [(5, 10)]),
    Probe("indels", "ACGT", "TTACNTTACGGT", 1, 0, AMBIGUOUS, [(3, 6), (8, 10)]),
    Probe("indels", "ACGT", "TTACGTTTACCGTTAGT", 1, 0, FIXED, [(3, 6), (11, 13), (15, 17)]),
    Probe("indels", "ACG", "GAGCG", 2, 0, FIXED, [(1, 1), (2, 3), (4, 5)]),
    Probe("indels", "ACGT", "ACGTACG", 0, 0, FIXED, [(1, 4)]),
    Probe("indels", "ACGT", "ACCGT", 1, 0, FIXED, [(3, 5)]),
    Probe("indels", "CTT", "CCTT", 1, 0, FIXED, [(1, 3), (2, 4)]),
    Probe("indels", "A", "CAC", 1, 0, FIXED, [(1, 1), (2, 2), (3, 3)]),
    Probe("indels", "AC", "GAGCG", 2, 0, FIXED, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]),
]


@pytest.mark.parametrize(
    "probe", PROBES, ids=lambda p: f"{p.search}-{p.pattern}-{p.subject}-{p.max_mm}"
)
def test_oracle_gives_the_ranges_of_r_probes(probe):
    if probe.search == "hamming":
        found = oracle.hamming(
            probe.pattern, probe.subject, probe.max_mm, probe.min_mm, *probe.fixed
        )
    else:
        found = oracle.indels(probe.pattern, probe.subject, probe.max_mm, *probe.fixed)
    assert found == probe.r


FIXED_KEYWORDS = {"fixed_pattern": True, "fixed_subject": True}


@pytest.mark.parametrize(
    "search",
    [
        lambda: oracle.hamming("", "ACGT", 1, 0, **FIXED_KEYWORDS),
        lambda: oracle.indels("", "ACGT", 1, **FIXED_KEYWORDS),
        lambda: oracle.indel_candidates("", "ACGT", 1, **FIXED_KEYWORDS),
        lambda: oracle.edit_distances("", "ACGT", 1, **FIXED_KEYWORDS),
    ],
    ids=["hamming", "indels", "indel_candidates", "edit_distances"],
)
def test_oracle_rejects_an_empty_pattern(search):
    with pytest.raises(ValueError, match=exactly("the pattern is empty")):
        search()


def test_documented_best_local_matches_rejects_a_limit_above_the_table_limit():
    table = oracle.edit_distances("ACGT", "TTACGT", 1, fixed_pattern=True, fixed_subject=True)
    with pytest.raises(ValueError, match="max_mm 2 is greater than the table limit 1"):
        oracle.documented_best_local_matches(table, 2)
