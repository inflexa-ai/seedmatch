"""The R search corpora in `tests/golden/search/`, grouped by case and by option set.

`RANDOM` holds short random cases with many edge cases. `SEEDS` holds patterns with the lengths
of the seeds and longer subjects.
"""

from collections import defaultdict
from dataclasses import dataclass
from functools import cache

from helpers import golden_rows

RANDOM = "corpus"
SEEDS = "seed_corpus"

type Range = tuple[int, int]


@dataclass(frozen=True)
class Options:
    max_mm: int
    min_mm: int
    indels: bool
    fixed_pattern: bool
    fixed_subject: bool

    @classmethod
    def of(cls, row) -> "Options":
        return cls(
            int(row.max_mm),
            int(row.min_mm),
            row.indels == "TRUE",
            row.fixed_pattern == "TRUE",
            row.fixed_subject == "TRUE",
        )


@dataclass(frozen=True)
class Case:
    case_id: str
    pattern: str
    subject: str


type Group = tuple[Case, Options, list[Range]]


@cache
def cases(corpus: str) -> dict[str, Case]:
    rows = golden_rows(f"search/{corpus}_cases.tsv")
    return {row.case_id: Case(row.case_id, row.pattern, row.subject) for row in rows}


@cache
def counts(corpus: str) -> dict[tuple[str, Options], int]:
    rows = golden_rows(f"search/{corpus}_counts.tsv")
    return {(row.case_id, Options.of(row)): int(row.count) for row in rows}


@cache
def ranges(corpus: str) -> dict[tuple[str, Options], list[Range]]:
    found: dict[tuple[str, Options], list[Range]] = defaultdict(list)
    for row in golden_rows(f"search/{corpus}_matches.tsv"):
        found[row.case_id, Options.of(row)].append((int(row.start), int(row.end)))
    return dict(found)


def groups(corpus: str) -> list[Group]:
    """Each case with each option set that R ran, and the R ranges in R order."""
    all_cases, all_ranges = cases(corpus), ranges(corpus)
    return [
        (all_cases[case_id], opts, all_ranges.get((case_id, opts), []))
        for case_id, opts in counts(corpus)
    ]
