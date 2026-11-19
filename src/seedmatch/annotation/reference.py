"""A reference: one annotation and one genome that agree."""

import os
from contextlib import ExitStack
from typing import Self

import pandas as pd

from seedmatch._checks import shown
from seedmatch.annotation.ensdb import AnnotationDb
from seedmatch.annotation.features import Features, derive_features
from seedmatch.annotation.filters import FeatureType, PerGene, TranscriptFilter
from seedmatch.annotation.genome import Genome, listed_names
from seedmatch.annotation.genome_file import open_genome
from seedmatch.errors import InputError


class Reference:
    """An annotation and a genome. The object holds both open, and `close` releases both."""

    def __init__(self, annotation: AnnotationDb, genome: Genome) -> None:
        """Join the annotation and the genome, and make sure that they are of the same build.

        The annotation must have a gene. Some contig of the annotation must be in the genome, with
        the names that `Genome.resolve` maps. A contig that has a length in the annotation must
        have the same length in the genome. Otherwise `InputError` is raised, and the caller keeps
        and closes the two objects.
        """
        _check_agreement(annotation, genome)
        self.annotation = annotation
        self.genome = genome

    @classmethod
    def from_paths(cls, *, ensdb: str | os.PathLike[str], genome: str | os.PathLike[str]) -> Self:
        """Open the EnsDb file and the genome file (2bit or FASTA).

        If a file does not open, or if the two do not agree, the function closes what it opened
        and raises the error.
        """
        with ExitStack() as stack:
            annotation = AnnotationDb.open(ensdb)
            stack.callback(annotation.close)
            opened = open_genome(genome)
            stack.callback(opened.close)
            reference = cls(annotation, opened)
            stack.pop_all()
        return reference

    def transcripts(self, transcript_filter: TranscriptFilter = TranscriptFilter()) -> pd.DataFrame:
        """Return the transcripts that pass the filter, as `AnnotationDb.transcripts` does.

        The default filter is the default of `load_annotations`.
        """
        return self.annotation.transcripts(transcript_filter)

    def features(
        self,
        feature: FeatureType,
        transcript_filter: TranscriptFilter = TranscriptFilter(),
        *,
        per_gene: PerGene = PerGene.EACH_TRANSCRIPT,
    ) -> Features:
        """Return the ranges of the feature of the transcripts that pass the filter.

        `per_gene` selects the transcripts of each gene. The default filter is the default of
        `load_annotations`.
        """
        for name, value, kind in (
            ("feature", feature, FeatureType),
            ("per_gene", per_gene, PerGene),
        ):
            if not isinstance(value, kind):
                message = (
                    f"{name} {shown(value)} has type {type(value).__name__}; "
                    f"expected a {kind.__name__}"
                )
                raise InputError(message)
        rows = self.annotation.exon_rows(
            transcript_filter, coding_only=feature is not FeatureType.EXONS
        )
        return derive_features(rows, feature, per_gene)

    def close(self) -> None:
        """Close the annotation and the genome."""
        try:
            self.annotation.close()
        finally:
            self.genome.close()

    def __enter__(self) -> Self:
        """Return the reference, so that `with` closes it."""
        return self

    def __exit__(self, *_exc_info: object) -> None:
        """Close the reference."""
        self.close()


def _check_agreement(annotation: AnnotationDb, genome: Genome) -> None:
    names = annotation.contig_names()
    if not names:
        message = (
            f"the annotation {annotation.source!r} has no genes; "
            "expected an annotation with at least one gene"
        )
        raise InputError(message)
    if all(genome.resolve(name) is None for name in names):
        message = (
            f"no contig of the annotation {annotation.source!r} is in the genome "
            f"{os.fspath(genome.path)!r}; annotation contigs: {listed_names(names)}; "
            f"genome contigs: {listed_names(genome.contigs)}"
        )
        raise InputError(message)
    for name, length in annotation.contig_lengths().items():
        contig = genome.resolve(name)
        if contig is not None and genome.contigs[contig] != length:
            shown = f"{name!r}" if contig == name else f"{name!r} (genome contig {contig!r})"
            message = (
                f"contig {shown}: the genome {os.fspath(genome.path)!r} has length "
                f"{genome.contigs[contig]}, but the annotation {annotation.source!r} has length "
                f"{length}; expected the same genome build"
            )
            raise InputError(message)
