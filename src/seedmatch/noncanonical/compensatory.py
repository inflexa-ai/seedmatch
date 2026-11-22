"""The 3' compensatory search: `get_compensatory_sequences` of R, with the rules of `D-004`."""

from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from typing import NamedTuple

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from seedmatch._checks import (
    as_limit,
    check_at_least,
    check_not_negative,
    check_type,
    concatenate_and_clear,
    shown,
)
from seedmatch._frames import with_text_dtype
from seedmatch._pool import map_in_pool, pool_width
from seedmatch.errors import InputError
from seedmatch.search import MAX_PATTERN_LENGTH, SearchOptions, match_ranges
from seedmatch.seed import Guide, Seed, SeedSpec, Target, seed_window
from seedmatch.sequence import SequenceSet, reverse_complement

_Ints = NDArray[np.int64]


@dataclass(frozen=True)
class CompensatoryOptions:
    """The options of the 3' compensatory search: the arguments of `get_compensatory_sequences`.

    A region is a window of the guide from g_a to g_b. It has a at least `first_position`, b at
    most the guide length, and a length from `min_length` to `max_length`. `max_length=None` is the
    R default `NA`: each position of the region. `first_position=None` is the guide position after
    the seed window, g9 for `Seed.MER7M8` (`D-056`). R always uses g8.

    `max_distance` is the largest number of bases between a region and its seed match. R `limit` L
    gives `max_distance` L - 2, thus the default 98 is the R default 100. `max_seed_mismatch` and
    `max_region_mismatch` are `mismatch_` and `mismatch_compensatory_region_` of R.

    `min_length`, `first_position`, and `max_distance` are at least 1, and `max_length` is at least
    `min_length`. The guide gives the largest value of a length (`D-055`).
    """

    min_length: int = 6
    max_length: int | None = None
    first_position: int | None = None
    max_distance: int = 98
    max_seed_mismatch: int = 1
    max_region_mismatch: int = 0

    def __post_init__(self) -> None:
        """Store the limits as plain ints, and raise `InputError` for a value that is not valid."""
        min_length = as_limit("min_length", self.min_length)
        max_length = None if self.max_length is None else as_limit("max_length", self.max_length)
        first = (
            None if self.first_position is None else as_limit("first_position", self.first_position)
        )
        max_distance = as_limit("max_distance", self.max_distance)
        seed_mismatch = as_limit("max_seed_mismatch", self.max_seed_mismatch)
        region_mismatch = as_limit("max_region_mismatch", self.max_region_mismatch)
        check_at_least("min_length", min_length, 1)
        if max_length is not None and max_length < min_length:
            message = (
                f"max_length {shown(max_length)} is less than min_length {shown(min_length)}; "
                f"expected None or an integer of {shown(min_length)} or more"
            )
            raise InputError(message)
        if first is not None:
            check_at_least("first_position", first, 1)
        check_at_least("max_distance", max_distance, 1)
        check_not_negative("max_seed_mismatch", seed_mismatch)
        check_not_negative("max_region_mismatch", region_mismatch)
        # The dataclass is frozen, so the constructor is the one writer of the plain ints.
        object.__setattr__(self, "min_length", min_length)
        object.__setattr__(self, "max_length", max_length)
        object.__setattr__(self, "first_position", first)
        object.__setattr__(self, "max_distance", max_distance)
        object.__setattr__(self, "max_seed_mismatch", seed_mismatch)
        object.__setattr__(self, "max_region_mismatch", region_mismatch)


class _SeedMatches(NamedTuple):
    """The seed matches that can have a region before them, sorted by `key`.

    The key of a match is the position of its start in the buffer of the set.
    """

    key: _Ints
    sequence: _Ints
    start: _Ints
    mismatches: _Ints


class _Pairs(NamedTuple):
    """The pairs of a seed match and a region match. `seed` is the position in `_SeedMatches`."""

    seed: _Ints
    begin: _Ints
    end: _Ints
    guide_start: _Ints
    guide_end: _Ints


def _first_position(guide: Guide, seed: SeedSpec, options: CompensatoryOptions) -> int:
    """Return the first guide position of a region, and check that a region fits.

    The guide must have room for the shortest region, and the longest region must not be longer
    than a pattern that a search accepts. The function checks both before the first search.
    """
    length = len(guide.sequence)
    first = options.first_position
    if first is None:
        first = seed_window(seed, length)[1] + 1
    if first > length:
        message = (
            f"the region starts at g{first}, but the guide ends at g{length}; "
            "expected a seed or a first_position that leaves a region"
        )
        raise InputError(message)
    positions = length - first + 1
    name, value = (
        ("min_length", options.min_length)
        if options.max_length is None
        else ("max_length", options.max_length)
    )
    if value > positions:
        noun = "position" if positions == 1 else "positions"
        message = (
            f"{name} {value} is more than the {positions} {noun} from g{first} to g{length}; "
            f"expected at most {positions}"
        )
        raise InputError(message)
    longest = positions if options.max_length is None else options.max_length
    if longest > MAX_PATTERN_LENGTH:
        message = (
            f"the longest region has {longest} letters; "
            f"expected at most {MAX_PATTERN_LENGTH}, the pattern limit of the search"
        )
        raise InputError(message)
    return first


def _region_target(guide: Guide, start: int, stop: int) -> str:
    """Return the DNA text that pairs with the window g_`start` to g_`stop` of the guide."""
    return reverse_complement(guide.sequence[start - 1 : stop].replace("U", "T"))


def _windows(
    guide: Guide, first: int, options: CompensatoryOptions
) -> Iterator[tuple[int, int, str]]:
    """Yield (guide start, guide end, region target) of each window, by length and then by start."""
    length = len(guide.sequence)
    longest = length - first + 1 if options.max_length is None else options.max_length
    for size in range(options.min_length, longest + 1):
        for start in range(first, length - size + 2):
            yield start, start + size - 1, _region_target(guide, start, start + size - 1)


def _seed_matches(
    target: Target, sequences: SequenceSet, max_mismatch: int, width: int
) -> _SeedMatches:
    """Return the seed matches with 0 to `max_mismatch` mismatches, with the mismatch count.

    A match with a start below 1 has no region before it, thus the function drops it. Its key
    would lie in an earlier sequence. A match has at most one mismatch for each letter of the
    target, thus a larger `max_mismatch` adds no search.
    """

    def search(mismatches: int) -> _SeedMatches:
        options = SearchOptions(max_mismatch=mismatches, min_mismatch=mismatches)
        index, start, _end = match_ranges(target, sequences, options)
        inside = start >= 1
        index, start = index[inside], start[inside]
        return _SeedMatches(
            sequences.offsets[index] + start - 1,
            index,
            start,
            np.full(start.size, mismatches, dtype=np.int64),
        )

    found = list(map_in_pool(search, range(min(max_mismatch, len(target.sequence)) + 1), width))
    key, index, start, mismatches = (
        concatenate_and_clear([match[field] for match in found]) for field in range(4)
    )
    order = np.argsort(key)
    return _SeedMatches(key[order], index[order], start[order], mismatches[order])


def _pair_regions(
    seeds: _SeedMatches, offsets: _Ints, reach: int, index: _Ints, end: _Ints
) -> tuple[_Ints, _Ints]:
    """Return the region and the seed of each pair of the region matches and the seed matches.

    A region match in sequence i with the end e pairs with the seed matches whose key is from
    `offsets[i] + e + 1` to `offsets[i] + e + reach`, and before the end of the sequence. Thus the
    distance is from 1 to `reach`, and a window does not run into the next sequence.
    """
    base = offsets[index] + end
    low = np.searchsorted(seeds.key, base + 1, side="left")
    top = np.minimum(base + reach, offsets[index + 1] - 1)
    high = np.searchsorted(seeds.key, top, side="right")
    counts = np.maximum(high - low, 0)
    region = np.repeat(np.arange(counts.size, dtype=np.int64), counts)
    seed = np.arange(region.size, dtype=np.int64) + np.repeat(
        low - (np.cumsum(counts) - counts), counts
    )
    return region, seed


def _find_pairs(
    sequences: SequenceSet,
    seeds: _SeedMatches,
    windows: Iterable[tuple[int, int, str]],
    options: CompensatoryOptions,
    width: int,
) -> _Pairs:
    """Return the pairs of the seed matches and the region matches of each window."""
    fields: tuple[list[_Ints], ...] = tuple([] for _ in _Pairs._fields)
    if seeds.key.size:
        # A distance cannot exceed the buffer, and the clip keeps `offsets + end + reach` in int64.
        reach = min(options.max_distance, int(sequences.offsets[-1]))
        region_options = SearchOptions(max_mismatch=options.max_region_mismatch)

        def search(window: tuple[int, int, str]) -> tuple[int, int, _Ints, _Ints, _Ints]:
            start, stop, pattern = window
            index, begin, end = match_ranges(pattern, sequences, region_options)
            return start, stop, index, begin, end

        for start, stop, index, begin, end in map_in_pool(search, windows, width):
            region, seed = _pair_regions(seeds, sequences.offsets, reach, index, end)
            if seed.size:
                parts = (
                    seed,
                    begin[region],
                    end[region],
                    np.full(seed.size, start, dtype=np.int64),
                    np.full(seed.size, stop, dtype=np.int64),
                )
                for field, part in zip(fields, parts, strict=True):
                    field.append(part)
    return _Pairs(
        *(
            concatenate_and_clear(field) if field else np.empty(0, dtype=np.int64)
            for field in fields
        )
    )


def _drop_contained(pairs: _Pairs) -> _Ints:
    """Return the positions of the pairs to keep.

    A pair is dropped if another pair of the same seed match has a region that contains its region
    and is longer (`D-053`). Pairs with the same region stay.
    """
    if not pairs.seed.size:
        return np.empty(0, dtype=np.int64)
    order = np.lexsort((-pairs.end, pairs.begin, pairs.seed))
    seed, begin, end = pairs.seed[order], pairs.begin[order], pairs.end[order]
    new_seed = np.ones(seed.size, dtype=np.bool_)
    new_seed[1:] = seed[1:] != seed[:-1]
    new_region = new_seed.copy()
    new_region[1:] |= (begin[1:] != begin[:-1]) | (end[1:] != end[:-1])
    first_of_region = np.maximum.accumulate(
        np.where(new_region, np.arange(seed.size, dtype=np.int64), 0)
    )
    # A region end is 1 or more (search rule 2), thus the levels of two seed matches do not overlap.
    level = (np.cumsum(new_seed) - 1) * (int(end.max()) + 1) + end
    before = np.concatenate((np.full(1, -1, dtype=np.int64), np.maximum.accumulate(level)[:-1]))
    return order[before[first_of_region] < level]


def _frame(
    guide: Guide, sequences: SequenceSet, target: Target, seeds: _SeedMatches, pairs: _Pairs
) -> pd.DataFrame:
    keep = _drop_contained(pairs)
    seed, begin, end = pairs.seed[keep], pairs.begin[keep], pairs.end[keep]
    start, stop = pairs.guide_start[keep], pairs.guide_end[keep]
    order = np.lexsort((start, end, begin, seed))
    seed, begin, end, start, stop = seed[order], begin[order], end[order], start[order], stop[order]
    # Each window has one code, thus a target and a GC content are made once for each window.
    shift = len(guide.sequence) + 1
    codes, inverse = np.unique(start * shift + stop, return_inverse=True)
    first, last = np.divmod(codes, shift)
    targets = [
        _region_target(guide, a, b) for a, b in zip(first.tolist(), last.tolist(), strict=True)
    ]
    percent = [(text.count("G") + text.count("C")) / len(text) * 100 for text in targets]
    seed_start = seeds.start[seed]
    index = seeds.sequence[seed]
    sites = pd.DataFrame(
        {
            "seq_id": np.array(sequences.ids, dtype=object)[index],
            "gene_id": np.array(sequences.gene_ids, dtype=object)[index],
            "seed_start": seed_start,
            "seed_end": seed_start + len(target.sequence) - 1,
            "seed_mismatches": seeds.mismatches[seed],
            "region_start": begin,
            "region_end": end,
            "region_length": end - begin + 1,
            "guide_start": start,
            "guide_end": stop,
            "region_target": np.array(targets, dtype=object)[inverse],
            "distance": seed_start - end - 1,
            "gc_content": np.array(percent, dtype=np.float64)[inverse],
        }
    )
    return with_text_dtype(sites, ("seq_id", "gene_id", "region_target"))


def compensatory_sites(
    guide: Guide,
    sequences: SequenceSet,
    options: CompensatoryOptions = CompensatoryOptions(),
    *,
    seed: SeedSpec = Seed.MER7M8,
    workers: int | None = None,
) -> pd.DataFrame:
    """Return the sites of the 3' compensatory search (`get_compensatory_sequences` of R).

    A site is a pair of a seed match and a region match of one sequence. The seed matches are the
    matches of `guide.target(seed)` with at most `options.max_seed_mismatch` mismatches, and a
    position outside the sequence is a mismatch. The regions are the windows of the guide from
    g_`first_position` on, with a length from `options.min_length` to `options.max_length`. A
    region is searched as the reverse complement of its window, with at most
    `options.max_region_mismatch` mismatches. The distance of a pair is the number of bases from
    the end of the region to the start of the seed match, from 1 to `options.max_distance`.
    A pair is dropped if a pair of the same seed match has a longer region that contains its
    region (`D-053`). Each sequence is searched, also if its ID occurs more than once (`D-052`).
    The search uses `fixed = TRUE` of R, thus a subject `N` is a mismatch.

    The function raises `InputError` before the first search if the guide has no room for a region,
    because `first_position` is after the guide end, or `min_length` or `max_length` is more than
    the number of positions (`D-055`). It also raises `InputError` if the longest region has more
    letters than a search accepts, 20000 (`D-055`).

    The columns are `seq_id`, `gene_id`, `seed_start`, `seed_end`, `seed_mismatches`,
    `region_start`, `region_end`, `region_length`, `guide_start`, `guide_end`, `region_target`,
    `distance`, and `gc_content`. The text columns have the text dtype of pandas, `gc_content` is
    `float64`, and the other columns are `int64`. The coordinates are 1-based and inclusive.
    `region_target` is the reverse complement of the DNA text of the window, the pattern of the
    search, and `gc_content` is its percentage of `G` and `C`. With `max_region_mismatch` above 0,
    the text in the sequence can differ from `region_target`. The columns `guide_start` and
    `guide_end` give the window of the row (`D-054`). The rows come in the order of the set, and
    then of `seed_start`, `region_start`, `region_end`, and `guide_start`, with a `RangeIndex`.
    `workers` is the width of the thread pool, and `None` is the number of CPUs. The result does not
    depend on `workers`.

    The function does `options.max_seed_mismatch` + 1 searches for the seed matches, and then one
    search for each window. A guide of N bases has W windows, the sum of N - `first_position` - l +
    2 over the lengths l. The default guide of 23 bases gives W = 55 and 57 searches, and
    `first_position=8` gives W = 66. A set without a seed match skips the window searches. W grows
    with the square of N, and `max_length` bounds it. The pairing takes O((S + R) log S + P log P)
    time for S seed matches, R region matches, and P pairs. A window has at most one match for each
    end, thus P is at most S x W x `options.max_distance`, which is 5,390 pairs for each seed match
    with the defaults. The memory is O(S + P) and the results that wait in the pool, because the
    function holds each pair until it has dropped the contained regions.

    Measured on the 18,614 canonical 3' UTRs of `rnor7` (17.9 million bases), on an Apple M2 Max
    with 12 CPUs that a coverage run shared: the defaults give 8,273 rows in 4,290 genes. One
    search takes 0.09 to 0.18 s, and the call takes 1.4 s with `workers=None` and 10 s with
    `workers=1`. The peak memory above the start is 12 to 17 MB, and the result holds 2.1 MB.

    P grows if `min_length` is short, if `max_region_mismatch` or `max_seed_mismatch` is high, or if
    the sequence has low complexity. The peak is from 150 to 390 bytes for each pair in these
    cases, and a large result adds to it. A guide of 23 `U` on a sequence of 2,000 `A` gives 10.4
    million pairs and 189,000 rows, and the peak is 1.6 GB. The default guide on the `rnor7` UTRs
    with `min_length=3` and `max_region_mismatch=1` gives 7.9 million pairs and 1.9 million rows,
    and the peak is 2.2 GB.
    """
    check_type("guide", guide, Guide)
    check_type("sequences", sequences, SequenceSet)
    check_type("options", options, CompensatoryOptions)
    width = pool_width(workers)
    target = guide.target(seed)
    first = _first_position(guide, seed, options)
    seeds = _seed_matches(target, sequences, options.max_seed_mismatch, width)
    windows = _windows(guide, first, options)
    pairs = _find_pairs(sequences, seeds, windows, options, width)
    return _frame(guide, sequences, target, seeds, pairs)
