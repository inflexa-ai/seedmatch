"""`full_search` and `full_search_ranges` of R: one guide, many seeds and edit distances.

The searches of one call are independent, so a thread pool runs them in parallel. The numba
kernels release the GIL. The output does not depend on the width of the pool.
"""

import os
from collections.abc import Callable, Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Final

import numpy as np
import pandas as pd

from seedmatch._checks import as_flag, as_limit, check_not_negative, shown
from seedmatch.errors import InputError
from seedmatch.search import SearchOptions, count_matches, find_matches
from seedmatch.seed import Guide, Seed, Target
from seedmatch.sequence import SequenceSet

DEFAULT_SEEDS: Final[tuple[Seed, ...]] = (
    Seed.FULL,
    Seed.MER18,
    Seed.MER15,
    Seed.MER8,
    Seed.MER7M8,
    Seed.MER7A1,
    Seed.MER6,
)
NONSEED_SEEDS: Final = frozenset({Seed.FULL, Seed.MER18, Seed.MER15})

_NO_MATCHES: Final = "No.matches"
# The order and the labels of `.get_best_match` (R/full_search.R:98-121). The first column with a
# count of 1 or more gives the label.
_BEST_MATCH_RANKING: Final = (
    *(
        (f"{seed.value}.mm{distance}", f"{seed.value}-D{distance}")
        for seed in (Seed.FULL, Seed.MER18, Seed.MER15)
        for distance in range(5)
    ),
    (f"{Seed.MER8.value}.mm0", "8mer"),
    (f"{Seed.MER7M8.value}.mm0", "7mer-m8"),
    (f"{Seed.MER7A1.value}.mm0", "7mer-A1"),
    (f"{Seed.MER6.value}.mm0", "6mer"),
)

_RANGES_DTYPES: Final = {
    "seq_id": "object",
    "gene_id": "object",
    "start": "int64",
    "end": "int64",
    "width": "int64",
    "query": "object",
    "target_seq": "object",
    "edit_distance": "int64",
    "indels": "bool",
    "group": "object",
}


def _items[T](name: str, values: tuple[T, ...] | list[T]) -> tuple[T, ...]:
    """Return the values as a tuple of one item or more.

    The order of the values gives the order of the output, so a set is not valid.
    """
    if not isinstance(values, tuple | list):
        message = (
            f"{name} {shown(values)} has type {type(values).__name__}; expected a tuple or a list"
        )
        raise InputError(message)
    if not values:
        message = f"{name} is empty; expected one value or more"
        raise InputError(message)
    return tuple(values)


def _check_seeds(seeds: tuple[Seed, ...] | list[Seed]) -> tuple[Seed, ...]:
    checked: list[Seed] = []
    for seed in _items("seeds", seeds):
        if not isinstance(seed, Seed):
            message = f"seeds holds {shown(seed)} of type {type(seed).__name__}; expected a Seed"
            raise InputError(message)
        if seed in checked:
            message = f"seeds repeats {seed.value!r}; expected each seed once"
            raise InputError(message)
        checked.append(seed)
    return tuple(checked)


def _check_distance(name: str, distance: int) -> int:
    limit = as_limit(name, distance)
    check_not_negative(name, limit)
    return limit


def _check_distances(name: str, distances: tuple[int, ...] | list[int]) -> tuple[int, ...]:
    checked: list[int] = []
    for distance in _items(name, distances):
        limit = _check_distance(name, distance)
        if limit in checked:
            message = f"{name} repeats {limit}; expected each distance once"
            raise InputError(message)
        checked.append(limit)
    return tuple(checked)


def _check_indel_options(flags: tuple[bool, ...] | list[bool]) -> tuple[bool, ...]:
    checked: list[bool] = []
    for flag in _items("indel_options", flags):
        value = as_flag("indel_options", flag)
        if value in checked:
            message = f"indel_options repeats {value}; expected each value once"
            raise InputError(message)
        checked.append(value)
    return tuple(checked)


def _pool_width(workers: int | None) -> int:
    if workers is None:
        # os.cpu_count() is None if the number of CPUs is unknown.
        return os.cpu_count() or 1
    width = as_limit("workers", workers)
    if width < 1:
        message = f"workers {width} is less than 1; expected an int of 1 or more"
        raise InputError(message)
    return width


@dataclass(frozen=True)
class FullSearchOptions:
    """The options of `full_search`, with the names of R in snake case.

    The defaults copy `R/full_search.R:26-35`. `nonseed_edit_distances` applies to the seeds in
    `NONSEED_SEEDS`, and `seed_edit_distances` applies to the other seeds. The distances are the
    `max_mismatch` of the searches. `fixed_pattern` and `fixed_subject` are the entries `pattern`
    and `subject` of the `fixed` argument of R. Each field is checked when the object is made.
    """

    seeds: tuple[Seed, ...] = DEFAULT_SEEDS
    nonseed_edit_distances: tuple[int, ...] = (0, 1, 2, 3, 4)
    seed_edit_distances: tuple[int, ...] = (0, 1)
    indels: bool = False
    wobble: bool = False
    fixed_pattern: bool = False
    fixed_subject: bool = False

    def __post_init__(self) -> None:
        """Make sure that the options are valid, and store plain tuples, ints, and bools."""
        # The dataclass is frozen, so the constructor is the one writer of the normalized values.
        set_field = object.__setattr__
        set_field(self, "seeds", _check_seeds(self.seeds))
        set_field(
            self,
            "nonseed_edit_distances",
            _check_distances("nonseed_edit_distances", self.nonseed_edit_distances),
        )
        set_field(
            self,
            "seed_edit_distances",
            _check_distances("seed_edit_distances", self.seed_edit_distances),
        )
        for name in ("indels", "wobble", "fixed_pattern", "fixed_subject"):
            set_field(self, name, as_flag(name, getattr(self, name)))


@dataclass(frozen=True)
class FullSearchRangesOptions:
    """The options of `full_search_ranges`, with the names of R in snake case.

    The defaults copy `R/full_search.R:148-155`. A seed outside `NONSEED_SEEDS` has no search with
    an edit distance above `max_seed_edit_distance`. Each field is checked when the object is made.
    """

    indel_options: tuple[bool, ...] = (True, False)
    edit_distances: tuple[int, ...] = (0, 1, 2, 3, 4)
    seeds: tuple[Seed, ...] = DEFAULT_SEEDS
    fixed_pattern: bool = False
    fixed_subject: bool = False
    max_seed_edit_distance: int = 1

    def __post_init__(self) -> None:
        """Make sure that the options are valid, and store plain tuples, ints, and bools."""
        # The dataclass is frozen, so the constructor is the one writer of the normalized values.
        set_field = object.__setattr__
        set_field(self, "indel_options", _check_indel_options(self.indel_options))
        set_field(self, "edit_distances", _check_distances("edit_distances", self.edit_distances))
        set_field(self, "seeds", _check_seeds(self.seeds))
        set_field(self, "fixed_pattern", as_flag("fixed_pattern", self.fixed_pattern))
        set_field(self, "fixed_subject", as_flag("fixed_subject", self.fixed_subject))
        set_field(
            self,
            "max_seed_edit_distance",
            _check_distance("max_seed_edit_distance", self.max_seed_edit_distance),
        )


@dataclass(frozen=True)
class _Search:
    name: str
    target: Target
    options: SearchOptions


def _run[T](work: Callable[[_Search], T], searches: Sequence[_Search], width: int) -> list[T]:
    """Return the result of each search, in the order of the searches."""
    with ThreadPoolExecutor(max_workers=width) as pool:
        return list(pool.map(work, searches))


def best_match(counts: pd.DataFrame) -> pd.Series:
    """Return the best match of each row of a table of `full_search` counts.

    The rank of a column and the label come from `.get_best_match` (`R/full_search.R:96-125`). The
    first column in rank order with a count of 1 or more gives the label of the row, and a row
    without such a count gets `No.matches`. A missing count is no match, as in R. A column that the
    table lacks is skipped (`D-012`).
    """
    labels = np.full(len(counts), _NO_MATCHES, dtype=object)
    present = [(column, label) for column, label in _BEST_MATCH_RANKING if column in counts.columns]
    # The loop goes from the lowest rank to the highest, so the highest rank writes last.
    for column, label in reversed(present):
        labels[(counts[column] >= 1).fillna(value=False).to_numpy(dtype=bool)] = label
    return pd.Series(labels, index=counts.index, name="best_match")


def full_search(
    guide: Guide,
    sequences: SequenceSet,
    options: FullSearchOptions = FullSearchOptions(),
    group: str = "siRNA",
    *,
    workers: int | None = None,
) -> pd.DataFrame:
    """Count the matches of each seed at each edit distance in each sequence (`full_search` of R).

    The columns are `seq_id`, `gene_id`, one count column `<seed>.mm<distance>` for each seed and
    distance in the order of the options, `group`, and `best_match`. A guide that is too short for
    a seed raises `InputError` before any search. `workers` is the width of the thread pool, and
    `None` is the number of CPUs. The result does not depend on `workers`.
    """
    width = _pool_width(workers)
    searches: list[_Search] = []
    for seed in options.seeds:
        target = guide.target(seed, wobble=options.wobble)
        distances = (
            options.nonseed_edit_distances if seed in NONSEED_SEEDS else options.seed_edit_distances
        )
        searches.extend(
            _Search(
                f"{seed.value}.mm{distance}",
                target,
                SearchOptions(
                    max_mismatch=distance,
                    indels=options.indels,
                    fixed_pattern=options.fixed_pattern,
                    fixed_subject=options.fixed_subject,
                ),
            )
            for distance in distances
        )
    counts = _run(lambda s: count_matches(s.target, sequences, s.options), searches, width)
    result = pd.DataFrame(
        {
            # A list of no text would make a float64 column and break the merge on gene_id.
            "seq_id": np.array(sequences.ids, dtype=object),
            "gene_id": np.array(sequences.gene_ids, dtype=object),
            **{search.name: found for search, found in zip(searches, counts, strict=True)},
        }
    )
    result["group"] = group
    result["best_match"] = best_match(result)
    return result


def full_search_ranges(
    guide: Guide,
    sequences: SequenceSet,
    options: FullSearchRangesOptions = FullSearchRangesOptions(),
    group: str = "siRNA",
    *,
    workers: int | None = None,
) -> pd.DataFrame:
    """Find the ranges of the matches of each search of a grid (`full_search_ranges` of R).

    The columns are `seq_id`, `gene_id`, `start`, `end`, `width`, `query`, `target_seq`,
    `edit_distance`, `indels`, and `group`. The searches follow the `expand.grid` order of R: the
    `indels` flag changes fastest, then the seed, then the edit distance. A search without a match
    has no row. A guide that is too short for a seed raises `InputError` before any search.
    `workers` is the width of the thread pool, and `None` is the number of CPUs. The result does
    not depend on `workers`.
    """
    width = _pool_width(workers)
    targets = {seed: guide.target(seed) for seed in options.seeds}
    searches = [
        _Search(
            f"{seed.value}.mm{distance}.indel{'TRUE' if indels else 'FALSE'}",
            targets[seed],
            SearchOptions(
                max_mismatch=distance,
                indels=indels,
                fixed_pattern=options.fixed_pattern,
                fixed_subject=options.fixed_subject,
            ),
        )
        for distance in options.edit_distances
        for seed in options.seeds
        for indels in options.indel_options
        if seed in NONSEED_SEEDS or distance <= options.max_seed_edit_distance
    ]
    found = _run(lambda s: find_matches(s.target, sequences, s.options), searches, width)
    pieces = []
    for search, matches in zip(searches, found, strict=True):
        if matches.empty:
            continue
        matches["query"] = search.name
        matches["target_seq"] = search.target.sequence
        matches["edit_distance"] = search.options.max_mismatch
        matches["indels"] = search.options.indels
        pieces.append(matches)
    if pieces:
        result = pd.concat(pieces, ignore_index=True)
    else:
        result = pd.DataFrame(
            {name: pd.Series(dtype=kind) for name, kind in _RANGES_DTYPES.items()}
        )
    result["group"] = group
    return result
