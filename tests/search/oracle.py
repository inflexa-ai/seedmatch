"""A slow reference version of the Biostrings search, written from the search rules.

The arguments of `hamming` and `indels` mirror `matchPattern(pattern, subject, max.mismatch,
min.mismatch, fixed)`, and the kernel tests pass them in this order.

The module also holds the "best local match" definition of the Biostrings documentation. Biostrings
does not obey it, so only a test of the R corpus uses it.
"""
# ruff: noqa: FBT001, PLR0913, PLR0917

from dataclasses import dataclass

from seedmatch.search.compat import compat_table
from seedmatch.sequence.alphabet import MASKS

type Range = tuple[int, int]


def _compatible(
    pattern: str, subject: str, fixed_pattern: bool, fixed_subject: bool
) -> list[list[bool]]:
    """`[a][j]` is true if letter `a` of the pattern matches letter `j` of the subject."""
    table = compat_table(fixed_pattern=fixed_pattern, fixed_subject=fixed_subject).tolist()
    return [[table[MASKS[p]][MASKS[s]] for s in subject] for p in pattern]


def hamming(
    pattern: str,
    subject: str,
    max_mm: int,
    min_mm: int,
    fixed_pattern: bool,
    fixed_subject: bool,
) -> list[Range]:
    """Return the mismatch matches of search rule 2."""
    m, n = len(pattern), len(subject)
    ok = _compatible(pattern, subject, fixed_pattern, fixed_subject)
    found = []
    for start in range(2 - m, n + 1):
        mismatches = sum(not (0 <= start - 1 + a < n and ok[a][start - 1 + a]) for a in range(m))
        if min_mm <= mismatches <= max_mm:
            found.append((start, start + m - 1))
    return found


@dataclass(frozen=True)
class Candidate:
    """The candidate of one start of the indel search. The coordinates are 1-based and inclusive."""

    start: int
    end: int
    cost: int


def _right_alignment(
    ok: list[list[bool]], first: int, start: int, budget: int
) -> tuple[int, int] | None:
    """Align `pattern[first:]` to a subject substring that starts at the 0-based `start`.

    Return the least edit distance and the shortest width with it, or None if the least edit
    distance is greater than `budget`.
    """
    rest = range(first, len(ok))
    n = len(ok[0])
    # A width greater than len(rest) + budget has an edit distance greater than budget.
    widths = min(n - start, len(rest) + budget)
    row = list(range(widths + 1))
    for a in rest:
        previous, row = row, [row[0] + 1]
        for w in range(1, widths + 1):
            substitution = previous[w - 1] + (0 if ok[a][start + w - 1] else 1)
            row.append(min(substitution, previous[w] + 1, row[w - 1] + 1))
        # The least value of a row never decreases in the next row.
        if min(row) > budget:
            return None
    cost = min(row)
    return cost, row.index(cost)


def indel_candidates(
    pattern: str, subject: str, max_mm: int, fixed_pattern: bool, fixed_subject: bool
) -> list[Candidate]:
    """Return the candidate of each start of search rule 4 that costs max_mm or less.

    The candidate of a start does not depend on max_mm. Thus, the candidates for one limit also
    give the candidates for each lower limit.
    """
    m, n = len(pattern), len(subject)
    ok = _compatible(pattern, subject, fixed_pattern, fixed_subject)
    candidates = []
    for j in range(n):
        first = next((a for a in range(m) if ok[a][j]), None)
        if first is None or first > max_mm:
            continue
        right = _right_alignment(ok, first + 1, j + 1, max_mm - first)
        if right is not None:
            cost, width = right
            candidates.append(Candidate(j + 1, j + 1 + width, first + cost))
    return candidates


def scan(candidates: list[Candidate], max_mm: int) -> list[Range]:
    """Report the candidates that cost max_mm or less, with the scan of search rule 4."""
    found = []
    held: Candidate | None = None
    for candidate in candidates:
        if candidate.cost > max_mm:
            continue
        if held is not None and candidate.end <= held.end:
            if candidate.cost <= held.cost:
                held = candidate
            continue
        if held is not None:
            found.append((held.start, held.end))
        held = candidate
    if held is not None:
        found.append((held.start, held.end))
    return found


def indels(
    pattern: str, subject: str, max_mm: int, fixed_pattern: bool, fixed_subject: bool
) -> list[Range]:
    """Return the matches of the Biostrings indel search (search rules 4 and 5)."""
    if max_mm == 0 or max_mm >= len(pattern):
        return hamming(pattern, subject, max_mm, 0, fixed_pattern, fixed_subject)
    candidates = indel_candidates(pattern, subject, max_mm, fixed_pattern, fixed_subject)
    return scan(candidates, max_mm)


@dataclass(frozen=True)
class EditTable:
    """`rows[i][w]` is the edit distance of the pattern and `subject[i : i + w]`.

    Each row holds the widths from 0 to `pattern_length + max_mm`, if the subject is long enough.
    """

    pattern_length: int
    max_mm: int
    rows: list[list[int]]


def edit_distances(
    pattern: str, subject: str, max_mm: int, fixed_pattern: bool, fixed_subject: bool
) -> EditTable:
    """Return the edit distances, with the compatibility table for the substitutions."""
    m, n = len(pattern), len(subject)
    ok = _compatible(pattern, subject, fixed_pattern, fixed_subject)
    rows = []
    for i in range(n):
        # After letter j, column[a] is the edit distance of pattern[:a] and subject[i : j + 1].
        column = list(range(m + 1))
        row = [m]
        for j in range(i, min(n, i + m + max_mm)):
            next_column = [column[0] + 1]
            for a in range(1, m + 1):
                substitution = column[a - 1] + (0 if ok[a - 1][j] else 1)
                next_column.append(min(substitution, column[a] + 1, next_column[a - 1] + 1))
            column = next_column
            row.append(column[m])
        rows.append(row)
    return EditTable(m, max_mm, rows)


def documented_best_local_matches(table: EditTable, max_mm: int) -> list[Range]:
    """Return the "best local matches" of the Biostrings documentation of `matchPattern`.

    A substring S' of the subject is a match if each of these conditions is true:

    - (a) nedit(P, S') <= max_mm
    - (b) each shorter substring of S' has an edit distance greater than nedit(P, S')
    - (c) each longer substring of the subject that contains S' has an edit distance that is not
      less than nedit(P, S')
    """
    if max_mm > table.max_mm:
        message = f"max_mm {max_mm} is greater than the table limit {table.max_mm}"
        raise ValueError(message)
    rows, n = table.rows, len(table.rows)
    # A substring of width w has an edit distance of at least |m - w|. Thus, a substring with a
    # width outside [low, high] can be neither a match nor a substring that breaks (b) or (c).
    low = max(1, table.pattern_length - max_mm)
    high = table.pattern_length + max_mm

    def shortest(i: int, w: int, e: int) -> bool:
        return all(
            rows[inner][width] > e
            for width in range(w - 1, low - 1, -1)
            for inner in range(i, i + w - width + 1)
        )

    def closest(i: int, w: int, e: int) -> bool:
        return all(
            rows[outer][width] >= e
            for width in range(w + 1, high + 1)
            for outer in range(max(0, i + w - width), min(i, n - width) + 1)
        )

    found = []
    for i in range(n):
        for w in range(low, min(high, n - i) + 1):
            e = rows[i][w]
            if e <= max_mm and shortest(i, w, e) and closest(i, w, e):
                found.append((i + 1, i + w))
    return found
