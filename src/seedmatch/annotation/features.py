"""The ranges of a feature of each transcript, and the selection of the transcripts of a gene."""

from dataclasses import dataclass
from itertools import compress
from typing import assert_never

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from seedmatch._checks import check_offsets, check_vector
from seedmatch.annotation.ensdb import ExonRows
from seedmatch.annotation.filters import FeatureType, PerGene
from seedmatch.errors import InputError

type _Ranges = tuple[NDArray[np.int64], NDArray[np.int64]]


@dataclass(frozen=True, eq=False)
class Features:
    """The ranges of a feature for each transcript or gene, with 1-based inclusive coordinates.

    Feature i has the ranges `starts[k]` to `ends[k]` for k in `offsets[i]:offsets[i + 1]`, in
    transcript order. Thus the ranges descend on the minus strand (-1). A feature can have no range,
    for example a transcript without a 3' UTR. The arrays are read-only.
    """

    ids: tuple[str, ...]
    gene_ids: tuple[str, ...]
    contigs: tuple[str, ...]
    strands: NDArray[np.int8]
    starts: NDArray[np.int64]
    ends: NDArray[np.int64]
    offsets: NDArray[np.int64]

    def __post_init__(self) -> None:
        """Make sure that the parts fit together, then make the arrays read-only."""
        for name, array, dtype in (
            ("strands", self.strands, np.int8),
            ("starts", self.starts, np.int64),
            ("ends", self.ends, np.int64),
            ("offsets", self.offsets, np.int64),
        ):
            check_vector(name, array, dtype)
        n = len(self.ids)
        lengths = {
            "gene_ids": len(self.gene_ids),
            "contigs": len(self.contigs),
            "strands": len(self.strands),
        }
        for name, length in lengths.items():
            if length != n:
                message = f"{name} has {length} entries; expected one for each of the {n} features"
                raise InputError(message)
        ranges = len(self.starts)
        if len(self.ends) != ranges:
            message = (
                f"ends has {len(self.ends)} entries; expected one for each of the {ranges} ranges"
            )
            raise InputError(message)
        check_offsets(self.offsets, n, "features", ranges, "the number of ranges")
        if (i := _first((self.strands != 1) & (self.strands != -1))) is not None:
            message = f"strands[{i}] is {self.strands[i]}; expected 1 or -1"
            raise InputError(message)
        if (i := _first(self.starts < 1)) is not None:
            message = f"starts[{i}] is {self.starts[i]}; expected 1 or more"
            raise InputError(message)
        if (i := _first(self.starts > self.ends)) is not None:
            message = (
                f"range {i} has start {self.starts[i]} after end {self.ends[i]}; "
                "expected a start that is not after the end"
            )
            raise InputError(message)
        for array in (self.strands, self.starts, self.ends, self.offsets):
            array.flags.writeable = False

    def __len__(self) -> int:
        """Return the number of features."""
        return len(self.ids)

    def widths(self) -> NDArray[np.int64]:
        """Return the sum of the range widths of each feature, thus the length of its sequence."""
        return np.diff(running_sum(self.ends - self.starts + 1)[self.offsets])

    def at_least(self, min_width: int) -> "Features":
        """Return the features whose sequence has at least `min_width` bases, in the same order."""
        keep = self.widths() >= min_width
        counts = np.diff(self.offsets)
        offsets = np.zeros(int(keep.sum()) + 1, dtype=np.int64)
        np.cumsum(counts[keep], out=offsets[1:])
        in_kept = np.repeat(keep, counts)
        selectors = keep.tolist()
        return Features(
            ids=tuple(compress(self.ids, selectors)),
            gene_ids=tuple(compress(self.gene_ids, selectors)),
            contigs=tuple(compress(self.contigs, selectors)),
            strands=self.strands[keep],
            starts=self.starts[in_kept],
            ends=self.ends[in_kept],
            offsets=offsets,
        )

    def to_frame(self) -> pd.DataFrame:
        """Return one row for each range, with the rank of the range in its feature from 1."""
        feature = np.repeat(np.arange(len(self)), np.diff(self.offsets))
        rank = np.arange(len(self.starts)) - self.offsets[feature] + 1
        return pd.DataFrame(
            {
                "id": np.array(self.ids, dtype=object)[feature],
                "gene_id": np.array(self.gene_ids, dtype=object)[feature],
                "contig": np.array(self.contigs, dtype=object)[feature],
                "strand": np.where(self.strands == 1, "+", "-").astype(object)[feature],
                "rank": rank,
                "start": self.starts,
                "end": self.ends,
            }
        )


def derive_features(rows: ExonRows, feature: FeatureType, per_gene: PerGene) -> Features:
    """Cut the feature from each exon of the rows, then select the transcripts of each gene.

    The 3' UTR, the 5' UTR, and the CDS have each transcript with a CDS. The exons have each
    transcript. `EACH_TRANSCRIPT` gives the transcripts in the order of the rows. `LONGEST` and
    `UNION` give the genes in the order of the gene ID.
    """
    starts, ends = _cut(rows, feature)
    bounds = _transcript_bounds(rows.tx_row)
    candidates = np.arange(len(bounds) - 1)
    if feature is not FeatureType.EXONS:
        candidates = candidates[rows.cds_start[bounds[:-1]] > 0]
    in_feature = starts <= ends
    match per_gene:
        case PerGene.EACH_TRANSCRIPT:
            chosen = candidates
        case PerGene.LONGEST:
            row_widths = np.where(in_feature, ends - starts + 1, 0)
            chosen = _longest(rows, bounds, candidates, np.diff(running_sum(row_widths)[bounds]))
        case PerGene.UNION:
            return _union(rows, bounds, candidates, (starts, ends), in_feature)
        case _:
            assert_never(per_gene)
    first = bounds[chosen]
    taken = _runs(first, bounds[chosen + 1])
    taken = taken[in_feature[taken]]
    counts = np.diff(running_sum(in_feature.astype(np.int64))[bounds])
    return Features(
        ids=tuple(rows.tx_id[first].tolist()),
        gene_ids=tuple(rows.gene_id[first].tolist()),
        contigs=tuple(rows.seq_name[first].tolist()),
        strands=rows.strand[first],
        starts=starts[taken],
        ends=ends[taken],
        offsets=running_sum(counts[chosen]),
    )


def _cut(rows: ExonRows, feature: FeatureType) -> _Ranges:
    """Return the part of each exon in the feature. A part with start > end is empty."""
    start, end = rows.exon_start, rows.exon_end
    match feature:
        case FeatureType.EXONS:
            return start, end
        case FeatureType.CDS:
            return np.maximum(start, rows.cds_start), np.minimum(end, rows.cds_end)
        case FeatureType.UTR3 | FeatureType.UTR5:
            # On the plus strand, the 3' UTR is above the CDS and the 5' UTR is below it. On the
            # minus strand, the two change places.
            above = rows.strand == (1 if feature is FeatureType.UTR3 else -1)
            return (
                np.where(above, np.maximum(start, rows.cds_end + 1), start),
                np.where(above, end, np.minimum(end, rows.cds_start - 1)),
            )
        case _:
            assert_never(feature)


def _longest(
    rows: ExonRows,
    bounds: NDArray[np.int64],
    candidates: NDArray[np.int64],
    widths: NDArray[np.int64],
) -> NDArray[np.int64]:
    """Return the transcript with the maximum width of each gene, in the order of the gene ID.

    A tie goes to the first transcript in the order of `Reference.transcripts`. The transcripts of
    a gene share its contig, thus this order is the order of `tx_seq_start` and the `tx` row.
    """
    first = bounds[candidates]
    _, gene = np.unique(rows.gene_id[first], return_inverse=True)
    order = np.lexsort((rows.tx_row[first], rows.tx_seq_start[first], -widths[candidates], gene))
    return candidates[order[_run_starts(gene[order])]]


def _union(
    rows: ExonRows,
    bounds: NDArray[np.int64],
    candidates: NDArray[np.int64],
    cut: _Ranges,
    in_feature: NDArray[np.bool_],
) -> Features:
    """Merge the ranges of the transcripts of each gene that overlap or touch.

    The genes come in the order of the gene ID, and the ranges of a gene in transcript order.
    """
    first = bounds[candidates]
    gene_ids, gene_first, gene = np.unique(
        rows.gene_id[first], return_index=True, return_inverse=True
    )
    gene_of_tx = np.full(len(bounds) - 1, -1)
    gene_of_tx[candidates] = gene
    row_gene = np.repeat(gene_of_tx, np.diff(bounds))
    taken = np.flatnonzero(in_feature & (row_gene >= 0))
    merged_gene, merged_starts, merged_ends = _merge(row_gene[taken], cut[0][taken], cut[1][taken])
    counts = np.bincount(merged_gene, minlength=len(gene_ids))
    offsets = running_sum(counts)
    gene_rows = first[gene_first]
    strands = rows.strand[gene_rows]
    position = np.arange(len(merged_starts))
    rank = position - offsets[merged_gene]
    place = np.where(strands[merged_gene] == 1, position, offsets[merged_gene + 1] - 1 - rank)
    starts = np.empty_like(merged_starts)
    ends = np.empty_like(merged_ends)
    starts[place] = merged_starts
    ends[place] = merged_ends
    ids = tuple(gene_ids.tolist())
    return Features(
        ids=ids,
        gene_ids=ids,
        contigs=tuple(rows.seq_name[gene_rows].tolist()),
        strands=strands,
        starts=starts,
        ends=ends,
        offsets=offsets,
    )


def _merge(
    gene: NDArray[np.int64], starts: NDArray[np.int64], ends: NDArray[np.int64]
) -> tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.int64]]:
    """Merge the ranges of each gene that overlap or touch, in the order of the gene and the start.

    Return the gene, the start, and the end of each merged range.
    """
    if not starts.size:
        return gene, starts, ends
    order = np.lexsort((starts, gene))
    gene, starts, ends = gene[order], starts[order], ends[order]
    # Each gene gets its own band of values, above the bands of the genes before it. Thus the
    # running maximum of the shifted ends is the running maximum of the ends of the same gene.
    band = gene * (int(ends.max()) + 2)
    reach = np.maximum.accumulate(ends + band) - band
    new = np.ones(len(starts), dtype=np.bool_)
    new[1:] = (gene[1:] != gene[:-1]) | (starts[1:] > reach[:-1] + 1)
    first = np.flatnonzero(new)
    return gene[first], starts[first], np.maximum.reduceat(ends, first)


def _first(mask: NDArray[np.bool_]) -> int | None:
    """Return the position of the first true value, or `None` if no value is true."""
    found = np.flatnonzero(mask)
    return int(found[0]) if found.size else None


def _transcript_bounds(tx_row: NDArray[np.int64]) -> NDArray[np.int64]:
    """Return the first row of each transcript, then the number of rows."""
    return np.append(_run_starts(tx_row), len(tx_row))


def _run_starts(values: NDArray[np.int64]) -> NDArray[np.int64]:
    """Return the position of each value that differs from the value before it."""
    if not values.size:
        return np.zeros(0, dtype=np.int64)
    return np.flatnonzero(np.append(True, values[1:] != values[:-1]))


def _runs(starts: NDArray[np.int64], stops: NDArray[np.int64]) -> NDArray[np.int64]:
    """Return the positions `starts[k]` to `stops[k] - 1` for each k, one run after the other."""
    lengths = stops - starts
    run_starts = running_sum(lengths)[:-1]
    return np.arange(int(lengths.sum())) - np.repeat(run_starts - starts, lengths)


def running_sum(values: NDArray[np.int64]) -> NDArray[np.int64]:
    """Return 0, then the sum of the first value, of the first two values, and so on."""
    total = np.zeros(len(values) + 1, dtype=np.int64)
    np.cumsum(values, out=total[1:])
    return total
