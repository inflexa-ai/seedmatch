"""The sequence of each feature, cut from a genome."""

import os
from typing import Final

import numpy as np
from numpy.typing import NDArray

from seedmatch.annotation.features import Features, running_sum
from seedmatch.annotation.genome import Genome, listed_names, outside_contig_error
from seedmatch.errors import InputError
from seedmatch.sequence.alphabet import reverse_complement_codes
from seedmatch.sequence.sequence_set import SequenceSet

_RANGES_PER_CHUNK: Final = 1 << 16


def extract_sequences(features: Features, genome: Genome) -> SequenceSet:
    """Return the sequence of each feature, from 5' to 3', as a `SequenceSet`.

    The sequence joins the ranges of the feature in list order. On the minus strand, each range is
    reverse complemented. This is `GenomicFeatures::extractTranscriptSeqs`. The IDs of the set are
    the IDs and the gene IDs of the features. A feature without a range has an empty sequence.
    `Genome.resolve` maps the contig of each feature that has a range. A contig that the genome
    lacks, and a range that ends after its contig, raise `InputError` before a read. The ranges of
    all features are read in the contig order of the genome, then by start, so a reader that keeps
    its last block reads each block once. The cost is linear in the total width of the ranges, plus
    a sort of the ranges.
    """
    contigs = _resolve_contigs(features, genome)
    names = list(genome.contigs)
    rank = {name: i for i, name in enumerate(names)}
    ranks_of_contigs = {contig: rank[found] for contig, found in contigs.items()}
    counts = np.diff(features.offsets)
    # A feature without a range has no resolved contig, and `np.repeat` drops its rank.
    feature_ranks = [ranks_of_contigs.get(contig, 0) for contig in features.contigs]
    ranks = np.repeat(np.array(feature_ranks, dtype=np.int64), counts)
    _check_ends(features, genome, ranks)
    strands = np.repeat(features.strands, counts)
    bounds = running_sum(features.ends - features.starts + 1)
    codes = np.empty(int(bounds[-1]), dtype=np.uint8)
    order = np.lexsort((features.starts, ranks))
    for first in range(0, len(order), _RANGES_PER_CHUNK):
        chunk = order[first : first + _RANGES_PER_CHUNK]
        rows = zip(
            ranks[chunk].tolist(),
            features.starts[chunk].tolist(),
            features.ends[chunk].tolist(),
            strands[chunk].tolist(),
            bounds[chunk].tolist(),
            strict=True,
        )
        for contig_rank, start, end, strand, position in rows:
            part = genome.codes(names[contig_rank], start, end)
            if strand == -1:
                part = reverse_complement_codes(part)
            codes[position : position + len(part)] = part
    offsets = bounds[features.offsets]
    return SequenceSet(ids=features.ids, gene_ids=features.gene_ids, codes=codes, offsets=offsets)


def _check_ends(features: Features, genome: Genome, ranks: NDArray[np.int64]) -> None:
    """Raise `InputError` for the first range that ends after its contig.

    `ranks` gives the index of the contig of each range in `genome.contigs`.
    """
    names = list(genome.contigs)
    lengths = np.array(list(genome.contigs.values()), dtype=np.int64)[ranks]
    outside = np.flatnonzero(features.ends > lengths)
    if len(outside):
        k = int(outside[0])
        feature = int(np.searchsorted(features.offsets, k, side="right")) - 1
        place = (
            f"contig {names[ranks[k]]!r} of {features.ids[feature]!r} "
            f"in the genome {os.fspath(genome.path)!r}"
        )
        raise outside_contig_error(
            place, int(features.starts[k]), int(features.ends[k]), int(lengths[k])
        )


def _resolve_contigs(features: Features, genome: Genome) -> dict[str, str]:
    """Map the contig of each feature with a range to its contig in the genome."""
    resolved: dict[str, str] = {}
    for i in np.flatnonzero(np.diff(features.offsets) > 0).tolist():
        contig = features.contigs[i]
        if contig in resolved:
            continue
        found = genome.resolve(contig)
        if found is None:
            message = (
                f"contig {contig!r} of {features.ids[i]!r} is not in the genome "
                f"{os.fspath(genome.path)!r}; expected one of {listed_names(genome.contigs)}"
            )
            raise InputError(message)
        resolved[contig] = found
    return resolved
