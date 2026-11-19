"""Numba kernels of the Biostrings search, for the search rules 2 and 4 of `docs/design.md`.

A kernel takes the arrays of a `SequenceSet` (`codes` and `offsets`), the pattern as masks, and a
table of `compat_table`. It trusts its caller: the pattern has one letter or more, and
`0 <= min_mm <= max_mm`. The coordinates are 1-based and inclusive. The matches of each sequence
come in the order of Biostrings, and the sequences come in the order of the set.

A ranges kernel counts the matches of each sequence first, and then fills arrays of the exact size.
"""

import numpy as np
from numba import njit
from numpy.typing import NDArray


@njit(cache=True, nogil=True)
def _letter_table(pattern: NDArray[np.uint8], compat: NDArray[np.bool_]) -> NDArray[np.bool_]:
    """Return the table `[a, s]`: true if letter `a` of the pattern matches the subject mask `s`."""
    masks = compat.shape[1]
    table = np.empty((pattern.size, masks), dtype=np.bool_)
    for a in range(pattern.size):
        for s in range(masks):
            table[a, s] = compat[pattern[a], s]
    return table


@njit(cache=True, nogil=True)
def _sequence_index(counts: NDArray[np.int64]) -> NDArray[np.int64]:
    """Return the sequence index of each match, for `counts[i]` matches in sequence i."""
    return np.repeat(np.arange(counts.size), counts)


@njit(cache=True, nogil=True)
def _mismatches(
    codes: NDArray[np.uint8], lo: int, n: int, ok: NDArray[np.bool_], start: int, max_mm: int
) -> int:
    """Count the mismatches of the pattern at the 0-based `start` of `codes[lo : lo + n]`.

    Each pattern letter outside the sequence is a mismatch. The count stops after it passes max_mm.
    """
    m = len(ok)
    first = max(0, -start)
    stop = min(m, n - start)
    count = first + m - stop
    a = first
    while count <= max_mm and a < stop:
        if not ok[a, codes[lo + start + a]]:
            count += 1
        a += 1
    return count


@njit(cache=True, nogil=True)
def _hamming_sequence(
    codes: NDArray[np.uint8],
    lo: int,
    hi: int,
    ok: NDArray[np.bool_],
    max_mm: int,
    min_mm: int,
    starts: NDArray[np.int64] | None,
    at: int,
) -> int:
    """Find the matches of search rule 2 in `codes[lo:hi]`, and return their number.

    If `starts` is not None, put the 1-based start of match t in `starts[at + t]`.
    """
    found = 0
    for start in range(1 - len(ok), hi - lo):
        count = _mismatches(codes, lo, hi - lo, ok, start, max_mm)
        if min_mm <= count <= max_mm:
            if starts is not None:
                starts[at + found] = start + 1
            found += 1
    return found


@njit(cache=True, nogil=True)
def hamming_counts(
    codes: NDArray[np.uint8],
    offsets: NDArray[np.int64],
    pattern: NDArray[np.uint8],
    compat: NDArray[np.bool_],
    max_mm: int,
    min_mm: int,
) -> NDArray[np.int64]:
    """Return the number of matches of search rule 2 in each sequence."""
    ok = _letter_table(pattern, compat)
    counts = np.zeros(offsets.size - 1, dtype=np.int64)
    for i in range(counts.size):
        counts[i] = _hamming_sequence(
            codes, offsets[i], offsets[i + 1], ok, max_mm, min_mm, None, 0
        )
    return counts


@njit(cache=True, nogil=True)
def hamming_ranges(
    codes: NDArray[np.uint8],
    offsets: NDArray[np.int64],
    pattern: NDArray[np.uint8],
    compat: NDArray[np.bool_],
    max_mm: int,
    min_mm: int,
) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    """Return the sequence index and the 1-based start of each match of search rule 2.

    The end of a match is its start + len(pattern) - 1.
    """
    counts = hamming_counts(codes, offsets, pattern, compat, max_mm, min_mm)
    index = _sequence_index(counts)
    starts = np.empty(index.size, dtype=np.int64)
    ok = _letter_table(pattern, compat)
    at = 0
    for i in range(counts.size):
        at += _hamming_sequence(codes, offsets[i], offsets[i + 1], ok, max_mm, min_mm, starts, at)
    return index, starts


@njit(cache=True, nogil=True)
def _first_letters(ok: NDArray[np.bool_]) -> NDArray[np.int64]:
    """Return the index of the first pattern letter that matches each subject mask.

    If no letter matches a mask, its value is the pattern length.
    """
    m, masks = ok.shape
    first = np.full(masks, m, dtype=np.int64)
    for s in range(masks):
        for a in range(m):
            if ok[a, s]:
                first[s] = a
                break
    return first


@njit(cache=True, nogil=True)
def _right_alignment(
    codes: NDArray[np.uint8],
    begin: int,
    hi: int,
    ok: NDArray[np.bool_],
    letter: int,
    budget: int,
    row: NDArray[np.int64],
) -> tuple[int, int]:
    """Align the pattern letters from `letter` to `codes[begin : begin + w]`, for each width w.

    Return the least edit distance and the shortest width with it. If the least distance is
    greater than `budget`, return budget + 1 and the width 0.
    """
    rest = len(ok) - letter
    widths = min(hi - begin, rest + budget)
    cap = budget + 1
    # After k letters, the distance at width w is |k - w| or more. Thus only the band of widths
    # k - budget to k + budget can have a distance of budget or less, and the DP computes only that
    # band. Each value saturates at cap. A value of budget or less comes only from values of budget
    # or less, so the saturation keeps it exact.
    low, high = 0, min(widths, budget)
    for w in range(high + 1):
        row[w] = w
    best = 0
    for k in range(1, rest + 1):
        low, high = max(0, k - budget), min(widths, k + budget)
        if low > high:
            return cap, 0
        if high == k + budget:
            # Row k - 1 did not compute this width.
            row[high] = cap
        if low == 0:
            diagonal = row[0]
            row[0] = k
            left = best = k
        else:
            diagonal = row[low - 1]
            left = best = cap
        for w in range(max(low, 1), high + 1):
            up = row[w]
            value = diagonal if ok[letter + k - 1, codes[begin + w - 1]] else diagonal + 1
            value = min(value, up + 1, left + 1, cap)
            diagonal = up
            row[w] = left = value
            best = min(best, value)
        if best > budget:
            return cap, 0
    width = low
    while row[width] != best:
        width += 1
    return best, width


@njit(cache=True, nogil=True)
def _put(
    starts: NDArray[np.int64] | None,
    ends: NDArray[np.int64] | None,
    t: int,
    start: int,
    end: int,
) -> None:
    """Put a match at index t, if the arrays are not None."""
    if starts is not None and ends is not None:
        starts[t] = start
        ends[t] = end


@njit(cache=True, nogil=True)
def _indel_sequence(
    codes: NDArray[np.uint8],
    lo: int,
    hi: int,
    ok: NDArray[np.bool_],
    first: NDArray[np.int64],
    max_mm: int,
    row: NDArray[np.int64],
    starts: NDArray[np.int64] | None,
    ends: NDArray[np.int64] | None,
    at: int,
) -> int:
    """Find the matches of search rule 4 in `codes[lo:hi]`, and return their number.

    If `starts` and `ends` are not None, put match t at index at + t.
    """
    found = 0
    held = False
    held_start = held_end = held_cost = 0
    for j in range(hi - lo):
        before = first[codes[lo + j]]
        if before > max_mm:
            continue
        distance, width = _right_alignment(
            codes, lo + j + 1, hi, ok, before + 1, max_mm - before, row
        )
        cost = before + distance
        if cost > max_mm:
            continue
        end = j + 1 + width
        if held and end <= held_end:
            if cost <= held_cost:
                held_start, held_end, held_cost = j + 1, end, cost
            continue
        if held:
            _put(starts, ends, at + found, held_start, held_end)
            found += 1
        held, held_start, held_end, held_cost = True, j + 1, end, cost
    if held:
        _put(starts, ends, at + found, held_start, held_end)
        found += 1
    return found


@njit(cache=True, nogil=True)
def _indel_tables(
    pattern: NDArray[np.uint8], compat: NDArray[np.bool_], max_mm: int
) -> tuple[NDArray[np.bool_], NDArray[np.int64], NDArray[np.int64]]:
    """Return the letter table, the first letters, and a DP row for the indel search."""
    ok = _letter_table(pattern, compat)
    # A candidate aligns len(pattern) - 1 letters or fewer, with a budget of max_mm or less. Thus
    # its widths are 0 to len(pattern) - 1 + max_mm.
    row = np.empty(pattern.size + max_mm, dtype=np.int64)
    return ok, _first_letters(ok), row


@njit(cache=True, nogil=True)
def indel_counts(
    codes: NDArray[np.uint8],
    offsets: NDArray[np.int64],
    pattern: NDArray[np.uint8],
    compat: NDArray[np.bool_],
    max_mm: int,
) -> NDArray[np.int64]:
    """Return the number of matches of the indel search in each sequence (search rules 4 and 5)."""
    if max_mm == 0 or max_mm >= pattern.size:
        return hamming_counts(codes, offsets, pattern, compat, max_mm, 0)
    ok, first, row = _indel_tables(pattern, compat, max_mm)
    counts = np.zeros(offsets.size - 1, dtype=np.int64)
    for i in range(counts.size):
        counts[i] = _indel_sequence(
            codes, offsets[i], offsets[i + 1], ok, first, max_mm, row, None, None, 0
        )
    return counts


@njit(cache=True, nogil=True)
def indel_ranges(
    codes: NDArray[np.uint8],
    offsets: NDArray[np.int64],
    pattern: NDArray[np.uint8],
    compat: NDArray[np.bool_],
    max_mm: int,
) -> tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.int64]]:
    """Return the sequence index, the start, and the end of each match of the indel search."""
    if max_mm == 0 or max_mm >= pattern.size:
        index, starts = hamming_ranges(codes, offsets, pattern, compat, max_mm, 0)
        return index, starts, starts + (pattern.size - 1)
    counts = indel_counts(codes, offsets, pattern, compat, max_mm)
    index = _sequence_index(counts)
    starts = np.empty(index.size, dtype=np.int64)
    ends = np.empty(index.size, dtype=np.int64)
    ok, first, row = _indel_tables(pattern, compat, max_mm)
    at = 0
    for i in range(counts.size):
        at += _indel_sequence(
            codes, offsets[i], offsets[i + 1], ok, first, max_mm, row, starts, ends, at
        )
    return index, starts, ends
