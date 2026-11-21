"""The k-mer counts of a set of sequences, as `Biostrings::oligonucleotideFrequency` gives them."""

import itertools
from collections.abc import Iterator
from dataclasses import dataclass
from typing import Final, NoReturn

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.sparse import coo_array, csr_array

from seedmatch._checks import (
    as_limit,
    check_at_least,
    check_type,
    concatenate_and_clear,
    shown,
)
from seedmatch.errors import InputError
from seedmatch.sequence.alphabet import MASKS
from seedmatch.sequence.sequence_set import SequenceSet

MAX_K: Final = 12
# A bin edge is at most 3.5 times the highest count, and it must fit in an int64.
_MAX_COUNT: Final = 2**61
_BLOCK_STARTS: Final = 2**20
_OTHER: Final = 4
_MASK_COUNT: Final = 16


def _digit_table() -> NDArray[np.uint8]:
    """Map each mask to the base-4 digit of A, C, G, or T. Each other mask maps to 4."""
    table = np.full(_MASK_COUNT, _OTHER, dtype=np.uint8)
    for digit, letter in enumerate("ACGT"):
        table[MASKS[letter]] = digit
    table.flags.writeable = False
    return table


_DIGITS: Final = _digit_table()


@dataclass(frozen=True, eq=False)
class KmerCounts:
    """The k-mer counts of a `SequenceSet`. Row i of `matrix` is sequence i, also for a repeated ID.

    Column j is the k-mer whose base-4 digits give j, with A = 0, C = 1, G = 2, and T = 3. The first
    letter is the most significant digit, thus the columns come in the order of Biostrings: `AA`,
    `AC`, `AG`, `AT`, `CA`, and so on. `matrix` is a `scipy.sparse.csr_array` of shape
    (len(ids), 4^k) and dtype `int64`, with sorted indices, no duplicate entry, and no stored 0.
    """

    k: int
    step: int
    ids: tuple[str, ...]
    gene_ids: tuple[str, ...]
    matrix: csr_array

    def kmers(self) -> tuple[str, ...]:
        """Return the label of each column, in column order.

        The tuple has 4^k strings: 16,384 for k = 7 and 16,777,216 (about 1 GB) for k = 12.
        """
        return tuple("".join(letters) for letters in itertools.product("ACGT", repeat=self.k))

    def to_frame(self) -> pd.DataFrame:
        """Return the layout of `sequence_kmer_counts` of R with `simplify.as = "matrix"`.

        The columns are `seq_id`, `gene_id`, and then one `int64` column for each k-mer, in column
        order. The frame has one row for each sequence and a `RangeIndex`.

        The frame is dense. It holds n · 4^k · 8 bytes, for n sequences: 3.0 GB for 23,000 sequences
        and k = 7. Use `matrix` for the counts, or `kmer_gene_counts`, as the small forms.
        """
        frame = pd.DataFrame(self.matrix.toarray(), columns=pd.Index(self.kmers()), copy=False)
        frame.insert(0, "gene_id", list(self.gene_ids))
        frame.insert(0, "seq_id", list(self.ids))
        return frame


def kmer_counts(sequences: SequenceSet, *, k: int = 7, step: int = 1) -> KmerCounts:
    """Count the k-mers of each sequence, as `sequence_kmer_counts` of R does.

    The defaults are the R defaults `width = 7` and `step = 1`. The windows of a sequence of length
    L start at the positions 1, 1 + step, 1 + 2 · step, and so on, and the last start is at most
    L - k + 1. A sequence shorter than k has no window. A window counts only if each of its letters
    is A, C, G, or T, in either case. A window with another IUPAC letter, for example N, does not
    count.

    `k` is from 1 to `MAX_K` (`D-043`). Each sequence has its own row, also if its ID occurs again
    (`D-044`).

    The function reads the window starts in blocks of 2^20 positions, and a block can split a
    sequence. The work memory of a block is about 68 bytes for each start, thus about 70 MB, and it
    does not grow with the set. The time is O(L · k) for the total length L, and O(W log W) for the
    sorts of the W windows of a block. The matrix holds at most one entry for each window, and an
    entry takes 12 bytes. The peak memory is about 2.3 times the size of the matrix, because the
    join of the blocks and the conversion to CSR each hold a copy.

    Measured on the 18,614 canonical 3' UTRs of rnor7 (17.9 million bases) with k = 7, on an Apple
    M2 Max: the call takes 1.2 s. The matrix has 15.7 million entries and 179 MB. The peak memory
    of the call is 421 MB.
    """
    check_type("sequences", sequences, SequenceSet)
    width = as_limit("k", k)
    if not 1 <= width <= MAX_K:
        message = f"k {shown(width)} is out of range; expected an integer from 1 to {MAX_K}"
        raise InputError(message)
    stride = as_limit("step", step)
    check_at_least("step", stride, 1)
    return KmerCounts(
        k=width,
        step=stride,
        ids=sequences.ids,
        gene_ids=sequences.gene_ids,
        matrix=_count_matrix(sequences, width, stride),
    )


def kmer_gene_counts(counts: KmerCounts) -> pd.Series:
    """Return the number of genes that have each k-mer, as the histogram of `kmer_hist` counts it.

    The result is an `int64` Series named `genes`. Its index is named `kmer` and holds
    `counts.kmers()`. Each k-mer has a value, also a k-mer that no gene has. The value of a k-mer
    is the number of distinct gene IDs with at least one sequence that has the k-mer. R counts the
    sequences (`D-042`).

    `kmer_hist` of R returns `integer(0)` and draws the histogram (`D-010`). R gives the count of
    the target as `gene.counts[target.seq]`. The port gives it as `result[target.sequence]`, for
    example `kmer_gene_counts(kmer_counts(sequences))[target.sequence]`.

    The work is the sparse product of a gene indicator matrix with `counts.matrix`. Each stored
    count is 1 or more, thus the structure of the product is the presence matrix of the genes. The
    memory is O(the entries of the matrix + 4^k). On the set of `kmer_counts` above, the call takes
    0.11 s, and its peak memory is 360 MB, because the product has an entry for each pair of a gene
    and a k-mer. For k = 12 the call takes about 3 s, and its peak memory is near 2 GB, because of
    the 16,777,216 labels.
    """
    check_type("counts", counts, KmerCounts)
    matrix = counts.matrix
    sequence_count = matrix.shape[0]
    codes, genes = pd.factorize(pd.Index(counts.gene_ids))
    indicator = csr_array(
        (np.ones(sequence_count, dtype=np.int64), (codes, np.arange(sequence_count))),
        shape=(len(genes), sequence_count),
    )
    # No stored count is negative. Thus a stored entry of the product is a gene with a sequence
    # that has the k-mer, and the product has the structure of the presence matrix of the genes.
    with_kmer = np.bincount((indicator @ matrix).indices, minlength=matrix.shape[1])
    return pd.Series(
        with_kmer.astype(np.int64), index=pd.Index(counts.kmers(), name="kmer"), name="genes"
    )


def kmer_histogram(gene_counts: pd.Series, *, max_bins: int = 50) -> pd.DataFrame:
    """Return the bins of the histogram of `kmer_hist`: how many k-mers have each gene count.

    `gene_counts` is the result of `kmer_gene_counts`, or another Series of integer counts with
    one value for each k-mer. A missing value, a negative count, and a count above 2^61 raise
    `InputError`, and so does a Series without a value. R draws `hist(gene.counts, breaks = 50)`,
    and its `pretty` breaks give bins of a fractional width, or more than 50 bins, for other counts
    (`D-066`). This function gives bins of integer counts. Let lo and hi be the lowest and the
    highest count. The width w of a bin is the first of 1, 2, 5, 10, 20, 50, 100, and so on, with
    `hi // w - lo // w + 1 <= max_bins`. Bin b holds the counts from `b * w` to `b * w + w - 1`,
    for b from `lo // w` to `hi // w`, thus each bin holds w counts, also the first.

    The result has the `int64` columns `low`, `high`, and `kmers`, one row for each bin from the
    bin of lo to the bin of hi, and a `RangeIndex`. A bin without a k-mer has a row with 0. `low`
    and `high` are the lowest and the highest count of the bin, and `kmers` is the number of
    k-mers with a count in the bin. `max_bins` is at least 1. The bins of R differ from these bins
    (`D-066`). For the test data with `max_bins=50`, the first bin of R holds the counts 0 to 2
    and has 3394 k-mers, and the first bin of this function holds 0 and 1 and has 2142.

    The time is O(n) for n k-mers, and the memory is O(n + `max_bins`), because one `np.bincount`
    counts the bins. For the 16,384 counts of k = 7 on the test data, the call takes 0.1 ms. For
    the 4^12 counts of k = 12, the call takes 0.05 seconds, and its peak memory above the input is
    128 MiB, one temporary array of n values. Machine: Apple M2 Max, Python 3.12.
    """
    check_type("gene_counts", gene_counts, pd.Series)
    bins = as_limit("max_bins", max_bins)
    check_at_least("max_bins", bins, 1)
    if not pd.api.types.is_integer_dtype(gene_counts.dtype):
        message = (
            f"gene_counts has dtype {gene_counts.dtype}; "
            "expected an integer dtype, as kmer_gene_counts gives"
        )
        raise InputError(message)
    if gene_counts.hasnans:
        label = gene_counts.index[int(np.argmax(gene_counts.isna().to_numpy()))]
        message = f"gene_counts has a missing value for k-mer {shown(label)}; expected a count"
        raise InputError(message)
    if gene_counts.empty:
        message = "gene_counts has no value; expected at least one k-mer"
        raise InputError(message)
    lowest, highest = int(gene_counts.min()), int(gene_counts.max())
    if lowest < 0:
        _raise_bad_count(gene_counts, gene_counts < 0, "0 or more")
    if highest > _MAX_COUNT:
        _raise_bad_count(gene_counts, gene_counts > _MAX_COUNT, "at most 2^61")
    width = next(w for w in _bin_widths() if highest // w - lowest // w + 1 <= bins)
    first = lowest // width
    counts = np.bincount(gene_counts.to_numpy(dtype=np.int64) // width - first)
    low = np.arange(first, first + counts.size, dtype=np.int64) * width
    return pd.DataFrame({"low": low, "high": low + (width - 1), "kmers": counts.astype(np.int64)})


def _bin_widths() -> Iterator[int]:
    """Yield the widths 1, 2, 5, 10, 20, 50, 100, and so on."""
    decade = 1
    while True:
        yield decade
        yield 2 * decade
        yield 5 * decade
        decade *= 10


def _raise_bad_count(gene_counts: pd.Series, offends: pd.Series, expected: str) -> NoReturn:
    """Raise `InputError` for the first count that `offends` marks. `expected` ends the message."""
    position = int(np.argmax(offends.to_numpy(dtype=bool)))
    message = (
        f"gene_counts holds {gene_counts.iloc[position]} "
        f"for k-mer {shown(gene_counts.index[position])}; expected {expected}"
    )
    raise InputError(message)


def _count_matrix(sequences: SequenceSet, k: int, step: int) -> csr_array:
    # Scipy keeps the dtype of the index arrays, so int32 gives the 12 bytes of an entry, not 16.
    # A column is below 4^k <= 2^24, and a row is below the number of sequences.
    index = np.int32 if len(sequences) <= np.iinfo(np.int32).max else np.int64
    # Each list starts with an empty array, because `np.concatenate` rejects a list of no array,
    # and a set with no window has no block.
    rows: list[NDArray[np.int32 | np.int64]] = [np.empty(0, dtype=index)]
    columns: list[NDArray[np.int32 | np.int64]] = [np.empty(0, dtype=index)]
    counts: list[NDArray[np.int64]] = [np.empty(0, dtype=np.int64)]
    last_start = len(sequences.codes) - k + 1
    for first in range(0, last_start, _BLOCK_STARTS):
        keys, block_counts = _count_block(
            sequences, first, min(first + _BLOCK_STARTS, last_start), k, step
        )
        block_rows, block_columns = np.divmod(keys, 4**k)
        rows.append(block_rows.astype(index))
        columns.append(block_columns.astype(index))
        counts.append(block_counts)
    # `concatenate_and_clear` empties its list, so the parts of one array are free before the next
    # array joins.
    row, column, count = (concatenate_and_clear(parts) for parts in (rows, columns, counts))
    # A sequence that two blocks split has an entry in each block. The COO form sums them.
    return coo_array((count, (row, column)), shape=(len(sequences), 4**k), dtype=np.int64).tocsr()


def _count_block(
    sequences: SequenceSet, first: int, stop: int, k: int, step: int
) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    """Count the windows that start at the positions first to stop - 1 of the buffer.

    Return the sorted keys (sequence · 4^k + window code) and the count of each key. The caller
    gives only starts with a whole window in the buffer.
    """
    offsets = sequences.offsets
    size = stop - first
    digits = _DIGITS[sequences.codes[first : stop + k - 1]]
    window = np.zeros(size, dtype=np.int64)
    for position in range(k):
        window *= 4
        window += digits[position : position + size]
    # A window has no other letter if the cumulative count is equal at both of its ends.
    others = np.zeros(size + k, dtype=np.int32)
    np.cumsum(digits == _OTHER, dtype=np.int32, out=others[1:])
    counted = others[k : k + size] == others[:size]
    # The first offset above a start is the end of the sequence of the start, also if an empty
    # sequence has the same offset.
    local = np.arange(first, stop, dtype=np.int64)
    sequence = np.searchsorted(offsets, local, side="right")
    counted &= local + k <= offsets[sequence]
    sequence -= 1
    local -= offsets[sequence]
    counted &= local % step == 0
    keys = sequence[counted]
    keys *= 4**k
    keys += window[counted]
    return np.unique(keys, return_counts=True)
