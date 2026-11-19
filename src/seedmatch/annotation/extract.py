"""The sequence of each feature, cut from a genome."""

import os
from itertools import pairwise

import numpy as np

from seedmatch.annotation.features import Features
from seedmatch.annotation.genome import Genome, listed_names
from seedmatch.errors import InputError
from seedmatch.sequence.alphabet import reverse_complement_codes
from seedmatch.sequence.sequence_set import SequenceSet


def extract_sequences(features: Features, genome: Genome) -> SequenceSet:
    """Return the sequence of each feature, from 5' to 3', as a `SequenceSet`.

    The sequence joins the ranges of the feature in list order. On the minus strand, each range is
    reverse complemented. This is `GenomicFeatures::extractTranscriptSeqs`. The IDs of the set are
    the IDs and the gene IDs of the features. A feature without a range has an empty sequence.
    `Genome.resolve` maps the contig of each feature that has a range. A contig that the genome
    lacks raises `InputError`. The cost is linear in the total width of the ranges.
    """
    contigs = _resolve_contigs(features, genome)
    codes = np.empty(int((features.ends - features.starts + 1).sum()), dtype=np.uint8)
    offsets = np.zeros(len(features) + 1, dtype=np.int64)
    bounds = features.offsets.tolist()
    position = 0
    for i, (contig, strand, (first, last)) in enumerate(
        zip(features.contigs, features.strands.tolist(), pairwise(bounds), strict=True)
    ):
        starts = features.starts[first:last].tolist()
        ends = features.ends[first:last].tolist()
        for start, end in zip(starts, ends, strict=True):
            part = genome.codes(contigs[contig], start, end)
            if strand == -1:
                part = reverse_complement_codes(part)
            codes[position : position + len(part)] = part
            position += len(part)
        offsets[i + 1] = position
    return SequenceSet(ids=features.ids, gene_ids=features.gene_ids, codes=codes, offsets=offsets)


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
