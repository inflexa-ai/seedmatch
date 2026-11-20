"""The e-value of a seed count against a null of shuffled guides (`generate_e_value` of R)."""

from collections.abc import Iterator
from typing import Final

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from seedmatch._checks import (
    as_data_frame,
    as_generator,
    as_limit,
    as_text,
    check_at_least,
    check_columns,
    check_gene_ids_are_text,
    check_numeric_columns,
    column_values,
    is_text,
    number_text,
    shown,
)
from seedmatch._pool import map_in_pool, pool_width
from seedmatch.errors import InputError
from seedmatch.search import SearchOptions, count_matches
from seedmatch.seed import Guide, Seed, SeedSpec, Target
from seedmatch.sequence import SequenceSet

_BLOCK_ROWS: Final = 4096


def _shuffled_guides(guide: Guide, count: int, generator: np.random.Generator) -> Iterator[Guide]:
    letters = np.frombuffer(guide.sequence.encode("ascii"), dtype=np.uint8)
    for _ in range(count):
        yield Guide(letters[generator.permutation(letters.size)].tobytes().decode("ascii"))


def _result_table(shuffles: int, genes: int) -> NDArray[np.int32]:
    return np.empty((shuffles, genes), dtype=np.int32)


def _is_count(values: NDArray[np.float64]) -> NDArray[np.bool_]:
    return np.isfinite(values) & (values >= 0)


# Each keyword-only parameter is an independent option of the null, so no parameter object fits.
def e_value_null(  # noqa: PLR0913
    guide: Guide,
    sequences: SequenceSet,
    *,
    rng: np.random.Generator | int,
    n_random: int = 1000,
    seed: SeedSpec = Seed.MER7M8,
    options: SearchOptions = SearchOptions(indels=True),
    workers: int | None = None,
) -> pd.DataFrame:
    """Return the seed counts of shuffled guides for each gene: the null of `e_values` (`D-008`).

    Shuffle j, for j = 0 to `n_random` - 1, has the letters of the guide in the order
    `rng.permutation(len(guide.sequence))`. Thus a shuffle keeps the length and the letters of the
    guide. The function draws the shuffles in the order j = 0, 1, and so on, in the thread of the
    caller, thus the result does not depend on `workers`. The count of a gene is the maximum of
    `count_matches` over the sequences of the gene (`D-003`), for the target of the shuffle and
    `options`. The default `options` are the default of `annotate_seed_counts`. A `seed` that does
    not fit the guide raises `InputError` before any draw.

    The result has the gene IDs of `sequences` as its index, named `gene_id`, in the order of the
    first sequence of each gene. It has one `int32` column for each shuffle, with the labels 0 to
    `n_random` - 1. `workers` is the width of the thread pool, and `None` is the number of CPUs.

    The time is `n_random` searches of the set on `workers` threads. The memory is the result,
    genes x `n_random` x 4 bytes, and O(`workers` x (genes + sequences)) more. For example, 20,000
    genes and 1000 shuffles take 80 MB for the result. The pool holds at most 2 `workers` + 1
    shuffles. Each of them has a vector of one count for each gene, and each search that runs holds
    one count for each sequence. The function draws a shuffle only when the pool has a free place,
    and it writes each vector into the result at once. It allocates the result before the first
    draw. Thus a table that the system cannot allocate raises `MemoryError` before any search, and
    `rng` does not move. A system that overcommits memory can accept a table that is larger than
    the memory.

    Measured on an Apple M2 Max with 12 CPUs: the default 1000 shuffles on the 18,614 canonical
    3' UTRs of `rnor7` (17.9 million bases) take 9 to 13 s with `workers=None`. One search takes
    0.09 s, thus the same null on one thread takes about 90 s. The result holds 71 MB, and the peak
    memory above the start is 75 MB.
    """
    generator = as_generator(rng)
    width = pool_width(workers)
    shuffles = as_limit("n_random", n_random)
    check_at_least("n_random", shuffles, 1)
    if not isinstance(options, SearchOptions):
        message = f"options has type {type(options).__name__}; expected a SearchOptions"
        raise InputError(message)
    guide.target(seed)

    codes, genes = pd.factorize(np.array(sequences.gene_ids, dtype=object))
    order = np.argsort(codes, kind="stable")
    sizes = np.bincount(codes, minlength=len(genes))
    starts = np.cumsum(sizes) - sizes

    def count_genes(target: Target) -> NDArray[np.int32]:
        counts = count_matches(target, sequences, options)
        # Each gene has a sequence, thus no group of reduceat is empty.
        return np.maximum.reduceat(counts[order], starts).astype(np.int32)

    # Each shuffle fills one row, so the table is the transpose of the frame. copy=False keeps the
    # frame on the memory of the table also under copy-on-write, the default of pandas 3.
    # The table exists before the first draw, thus a table that the system cannot allocate fails
    # before the generator of the caller moves.
    table = _result_table(shuffles, len(genes))
    targets = (shuffle.target(seed) for shuffle in _shuffled_guides(guide, shuffles, generator))
    for j, found in enumerate(map_in_pool(count_genes, targets, width)):
        table[j] = found
    return pd.DataFrame(
        table.T, index=pd.Index(genes, name="gene_id"), columns=pd.RangeIndex(shuffles), copy=False
    )


def _check_null(null: object) -> pd.DataFrame:
    frame = as_data_frame(
        null, name="null", expected="a pandas.DataFrame with one row for each gene"
    )
    index = frame.index
    if not is_text(index.to_series()):
        message = f"null index has dtype {index.dtype}; expected text gene IDs"
        raise InputError(message)
    if index.hasnans:
        message = "null index has a missing gene ID; expected text gene IDs"
        raise InputError(message)
    if not index.is_unique:
        gene = index[int(np.argmax(index.duplicated()))]
        message = (
            f"null has {int((index == gene).sum())} rows for gene {shown(gene)}; "
            "expected one row for each gene"
        )
        raise InputError(message)
    if len(frame) == 0:
        message = "null has no row; expected one row for each gene"
        raise InputError(message)
    if frame.shape[1] == 0:
        message = "null has no column; expected one column for each random guide"
        raise InputError(message)
    if not frame.columns.is_unique:
        label = frame.columns.tolist()[int(np.argmax(frame.columns.duplicated()))]
        message = f"null repeats the column label {shown(label)}; expected each label once"
        raise InputError(message)
    check_numeric_columns(frame, list(frame.columns))
    for start in range(0, len(frame), _BLOCK_ROWS):
        block = frame.iloc[start : start + _BLOCK_ROWS].to_numpy(dtype=np.float64, na_value=np.nan)
        valid = _is_count(block)
        if not valid.all():
            row, found = np.argwhere(~valid)[0]
            message = (
                f"null holds {number_text(float(block[row, found]))} "
                f"for gene {shown(frame.index[start + row])}; expected counts of 0 or more"
            )
            raise InputError(message)
    return frame


def e_values(de: pd.DataFrame, null: pd.DataFrame, *, column: str) -> pd.DataFrame:
    """Return the rows of `de` with the e-value of their count (`generate_e_value` of R).

    `de` has the columns `gene_id` (text) and `column` (numeric). `null` has one row for each gene,
    with the gene ID as the index, and one column for each random guide, for example the result of
    `e_value_null`. Genes match by `gene_id` (`D-008`). The e-value of a row is the number of
    columns of `null` whose count for the gene is at least the count of the row, divided by the
    number of columns. A count of 0 gives the e-value 1. Each row has its own e-value, also for a
    gene with two rows (`D-047`).

    The result is a copy of the rows of `de` whose gene is in the index of `null`, in the order of
    `de`, with a new index and the `float64` column `e_value`. A column `e_value` that is present
    gets the new values. A `de` without such a row raises `InputError`. A count that is NaN,
    infinite, or negative, in `de` or in `null`, raises `InputError`, and so does a gene with two
    rows in `null` (`D-047`). A count that is not an integer is valid. If `gene_id` or the name in
    `column` labels more than one column of `de`, the function raises `InputError`.

    The time is O(R x C) for R rows of `de` and C columns of `null`. The function compares blocks
    of 4096 rows, thus the work memory is O(4096 x C), and it makes no copy of `null`.
    """
    table = as_data_frame(de)
    name = as_text("column", column)
    check_columns(table, ("gene_id", name))
    check_gene_ids_are_text(table)
    check_numeric_columns(table, (name,))
    observed = column_values(table, name, _is_count, "a count of 0 or more")
    frame = _check_null(null)

    positions = frame.index.get_indexer(pd.Index(table["gene_id"].to_numpy(dtype=object)))
    found = positions >= 0
    if not found.any():
        message = (
            "de has no gene of the null; expected the gene IDs of the null in the column gene_id, "
            f"for example {shown(frame.index[0])}"
        )
        raise InputError(message)
    result = table.loc[found].reset_index(drop=True)
    rows, counts = positions[found], observed[found]
    fractions = np.empty(len(rows), dtype=np.float64)
    for start in range(0, len(rows), _BLOCK_ROWS):
        stop = start + _BLOCK_ROWS
        block = frame.iloc[rows[start:stop]].to_numpy(dtype=np.float64, na_value=np.nan)
        at_least = (block >= counts[start:stop, np.newaxis]).sum(axis=1)
        fractions[start:stop] = at_least / block.shape[1]
    result["e_value"] = fractions
    return result
