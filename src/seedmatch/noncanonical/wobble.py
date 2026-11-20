"""The wobble and bulge search: `run_wobble_bulges_pipeline` of R, with the rules of `D-005`."""

from dataclasses import dataclass
from typing import Final

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from seedmatch._checks import (
    as_flag,
    as_limit,
    as_text,
    check_not_negative,
    check_type,
    concatenate_and_clear,
    shown,
)
from seedmatch._pool import map_in_pool, pool_width
from seedmatch.errors import InputError
from seedmatch.report.annotate import add_gene_counts, check_count_column
from seedmatch.search import SearchOptions, match_ranges
from seedmatch.search.api import _MAX_PATTERN_LENGTH
from seedmatch.seed import Guide, Seed
from seedmatch.sequence import SequenceSet

# The target of mer7m8 has 7 letters. Target position k pairs with the guide position g(9 - k), so
# t1 pairs with g8 and t7 pairs with g2.
_SEED_END: Final = 8
_TARGET_LENGTH: Final = _SEED_END - 1
_MAX_BULGE_FREE: Final = _SEED_END - 1
# R searches with fixed = FALSE (R/run_wobble_bulges_pipeline.R:39, :42), so N and R are wildcards.
_SEARCH_OPTIONS: Final = SearchOptions(fixed_pattern=False, fixed_subject=False)


@dataclass(frozen=True)
class WobbleBulgeOptions:
    """The options of the wobble and bulge search: `bulges_lim` and `nt` of R, with its defaults.

    `max_bulge` is the largest number of extra target bases in one bulge. The value 0 gives the
    search without a bulge. `bulge_free` is the number of seed positions from g2 that no bulge is
    next to, from 0 to 7. With the default 3, a bulge can lie only between g5 and g6, g6 and g7, or
    g7 and g8.

    The longest pattern has 7 + `max_bulge` letters, and a search accepts at most 20000 letters.
    Thus `max_bulge` is at most 19993, and a larger value raises `InputError`.
    """

    max_bulge: int = 4
    bulge_free: int = 3

    def __post_init__(self) -> None:
        """Store the limits as plain ints, and raise `InputError` for a value that is not valid."""
        max_bulge = as_limit("max_bulge", self.max_bulge)
        bulge_free = as_limit("bulge_free", self.bulge_free)
        check_not_negative("max_bulge", max_bulge)
        longest = _TARGET_LENGTH + max_bulge
        if longest > _MAX_PATTERN_LENGTH:
            message = (
                f"max_bulge {shown(max_bulge)} makes a pattern of {shown(longest)} letters; "
                f"expected at most {_MAX_PATTERN_LENGTH - _TARGET_LENGTH}"
            )
            raise InputError(message)
        if not 0 <= bulge_free <= _MAX_BULGE_FREE:
            message = (
                f"bulge_free {shown(bulge_free)} is out of range; "
                f"expected an integer from 0 to {_MAX_BULGE_FREE}"
            )
            raise InputError(message)
        # The dataclass is frozen, so the constructor is the one writer of the plain ints.
        object.__setattr__(self, "max_bulge", max_bulge)
        object.__setattr__(self, "bulge_free", bulge_free)


def _shapes(options: WobbleBulgeOptions) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    """Return `bulge_after` and `bulge_size` of each pattern, in the order of the patterns.

    The unbulged pattern comes first. The bulged patterns follow by size, and then by gap. A gap
    after g_i has i from `bulge_free` + 2 to 7. Thus no bulge is next to g2 to g(`bulge_free` + 1).
    """
    gaps = np.arange(options.bulge_free + 2, _SEED_END, dtype=np.int64)
    sizes = np.arange(1, options.max_bulge + 1, dtype=np.int64)
    after = np.concatenate((np.zeros(1, dtype=np.int64), np.tile(gaps, sizes.size)))
    size = np.concatenate((np.zeros(1, dtype=np.int64), np.repeat(sizes, gaps.size)))
    return after, size


def _pattern(target: str, after: int, size: int) -> str:
    """Return the target with `size` extra bases after g_`after`. Size 0 gives the target."""
    cut = _SEED_END - after
    return target[:cut] + "N" * size + target[cut:]


def wobble_bulge_patterns(
    guide: Guide, options: WobbleBulgeOptions = WobbleBulgeOptions()
) -> pd.DataFrame:
    """Return the patterns of the wobble and bulge search, in the order of the search.

    The first pattern is the target of `Seed.MER7M8` with each guide U as `R`, thus one pattern
    covers each set of G:U pairs (`D-005`). A bulge is a run of 1 to `options.max_bulge` extra
    target bases between two seed positions. Its pattern has the letter `N` for each extra base.
    The bulged patterns follow, in the order of the bulge size, and then of the position.

    The columns are `pattern` (text), `bulge_after`, and `bulge_size` (`int64`). `bulge_after` is
    the guide position g_i that the bulge follows, so that the bulge lies between g_i and g_(i+1).
    The first pattern has `bulge_after` 0 and `bulge_size` 0. The frame has a `RangeIndex` and
    1 + `max_bulge` x max(0, 6 - `bulge_free`) rows. Each row holds a pattern of up to 7 +
    `max_bulge` letters.

    The function builds all patterns at once. The time and the memory are the sum of the pattern
    lengths, O(`max_bulge`^2 x (6 - `bulge_free`)). For `max_bulge` 19993 and `bulge_free` 0, the
    patterns hold about 1.2 billion letters. The defaults give 13 short patterns.
    """
    check_type("guide", guide, Guide)
    check_type("options", options, WobbleBulgeOptions)
    target = guide.target(Seed.MER7M8, wobble=True).sequence
    after, size = _shapes(options)
    return pd.DataFrame(
        {
            "pattern": [
                _pattern(target, a, s) for a, s in zip(after.tolist(), size.tolist(), strict=True)
            ],
            "bulge_after": after,
            "bulge_size": size,
        }
    )


def _distinct_starts(
    guide: Guide, sequences: SequenceSet, options: WobbleBulgeOptions, workers: int | None
) -> tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.int64]]:
    """Return the sequence index, the start, and the pattern number of each distinct start.

    The rows have the order of the set and then of the start. A start that more than one pattern
    matches takes the first pattern of the order of `wobble_bulge_patterns`.
    """
    check_type("guide", guide, Guide)
    check_type("sequences", sequences, SequenceSet)
    check_type("options", options, WobbleBulgeOptions)
    width = pool_width(workers)
    target = guide.target(Seed.MER7M8, wobble=True).sequence
    after, size = _shapes(options)

    def search(shape: tuple[int, int]) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
        index, starts, _ends = match_ranges(_pattern(target, *shape), sequences, _SEARCH_OPTIONS)
        return index, starts

    indexes: list[NDArray[np.int64]] = []
    starts: list[NDArray[np.int64]] = []
    numbers: list[NDArray[np.int64]] = []
    shapes = zip(after.tolist(), size.tolist(), strict=True)
    for position, (found, begin) in enumerate(map_in_pool(search, shapes, width)):
        indexes.append(found)
        starts.append(begin)
        numbers.append(np.full(found.size, position, dtype=np.int64))
    index = concatenate_and_clear(indexes)
    start = concatenate_and_clear(starts)
    number = concatenate_and_clear(numbers)
    order = np.lexsort((number, start, index))
    index = index[order]
    start = start[order]
    number = number[order]
    first = np.ones(index.size, dtype=np.bool_)
    first[1:] = (index[1:] != index[:-1]) | (start[1:] != start[:-1])
    return index[first], start[first], number[first]


def wobble_bulge_matches(
    guide: Guide,
    sequences: SequenceSet,
    options: WobbleBulgeOptions = WobbleBulgeOptions(),
    *,
    workers: int | None = None,
) -> pd.DataFrame:
    """Return the sites of the wobble and bulge search (`run_wobble_bulges_pipeline` of R).

    The patterns come from `wobble_bulge_patterns`. Each pattern is searched with `fixed = FALSE`
    as R does, so a subject `N` matches each pattern letter. A sequence has a row for each distinct
    start of a match. The row takes the first pattern in pattern order that matches at this start,
    and its `end` and `width` come from this pattern. A bulge pattern and the unbulged pattern that
    start at different positions and end at the same position thus give two rows (`D-005`).

    The columns are `seq_id`, `gene_id`, `start`, `end`, `width`, `bulge_after`, and `bulge_size`.
    The integer columns are `int64`, and the coordinates are 1-based and inclusive. The rows come
    in the order of the set, and then of `start`, with a `RangeIndex`. `workers` is the width of the
    thread pool, and `None` is the number of CPUs. The result does not depend on `workers`.

    The time is 1 + `max_bulge` x max(0, 6 - `bulge_free`) searches of the set on `workers`
    threads, 13 for the defaults. The memory is O(the hits of all patterns): the function keeps the
    sequence index, the start, and the pattern number of each hit until it has sorted them. Each hit
    takes about 70 bytes. The worst case is a sequence of N, which each pattern matches at each
    position, thus there is one hit for each pattern and base. With the defaults, 4 million N take
    3.6 GB.

    Measured on the 18,614 canonical 3' UTRs of `rnor7` (17.9 million bases), on an Apple M2 Max
    with 12 CPUs: the defaults give 66,374 rows. One search takes 0.16 s, and the call takes 0.40 s
    with `workers=None` and 2.1 s with `workers=1`. The peak memory above the start is about 20 MB,
    and the result holds 11 MB.
    """
    index, start, number = _distinct_starts(guide, sequences, options, workers)
    after, size = _shapes(options)
    seq_ids = np.array(sequences.ids, dtype=object)
    gene_ids = np.array(sequences.gene_ids, dtype=object)
    width = _TARGET_LENGTH + size[number]
    return pd.DataFrame(
        {
            "seq_id": seq_ids[index],
            "gene_id": gene_ids[index],
            "start": start,
            "end": start + width - 1,
            "width": width,
            "bulge_after": after[number],
            "bulge_size": size[number],
        }
    )


# Each keyword-only parameter is an independent option of the search, so no parameter object fits.
def annotate_wobble_bulge_counts(  # noqa: PLR0913
    de: pd.DataFrame,
    sequences: SequenceSet,
    guide: Guide,
    *,
    column: str = "overall_matches",
    options: WobbleBulgeOptions = WobbleBulgeOptions(),
    shared_genes: bool = True,
    workers: int | None = None,
) -> pd.DataFrame:
    """Return a copy of the DE table with the number of wobble and bulge sites of each gene.

    The count of a sequence is its number of rows in `wobble_bulge_matches`. The count of a gene is
    the maximum over its sequences (`D-003`), and 0 for a gene without a sequence. The column is
    named `column`, and its default is the column of R. A column with that name is replaced in
    place (`D-015`). With `shared_genes`, the result keeps only the rows whose gene has a sequence,
    in the order of the input, and has a new index. A repeated `gene_id` column, or a repeated
    column with the name of the count column, raises `InputError` (`D-051`). The function checks
    the table before it searches. R gives one column for each wobble seed, and the port gives only
    the total (`D-005`).

    The time and the memory are those of `wobble_bulge_matches`, and O(sequences + genes) more: the
    function keeps one count for each sequence, and it does not build the frame of the sites.
    """
    shared_genes = as_flag("shared_genes", shared_genes)
    name = as_text("column", column)
    table = check_count_column(de, name)
    index, _start, _number = _distinct_starts(guide, sequences, options, workers)
    counts = np.bincount(index, minlength=len(sequences)).astype(np.int64, copy=False)
    return add_gene_counts(table, sequences, counts, name=name, shared_genes=shared_genes)
