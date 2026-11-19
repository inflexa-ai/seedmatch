"""A reference: one annotation and one genome that agree."""

import os
from collections.abc import Callable
from contextlib import ExitStack
from pathlib import Path
from typing import Final, Self

import pandas as pd

from seedmatch._checks import as_limit, check_not_negative, shown
from seedmatch._io import local_path
from seedmatch.annotation.ensdb import AnnotationDb
from seedmatch.annotation.extract import extract_sequences
from seedmatch.annotation.features import Features, derive_features
from seedmatch.annotation.filters import FeatureType, PerGene, TranscriptFilter
from seedmatch.annotation.genome import Genome, listed_names
from seedmatch.annotation.genome_file import open_genome
from seedmatch.annotation.gtf import load_gtf
from seedmatch.errors import InputError
from seedmatch.sequence.sequence_set import SequenceSet

_ENSDB_SUFFIXES: Final = (".sqlite",)
_GENOME_SUFFIXES: Final = (".2bit", ".fa.gz")


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
        return cls._open(AnnotationDb.open, ensdb, genome)

    @classmethod
    def from_gtf(cls, gtf: str | os.PathLike[str], genome: str | os.PathLike[str]) -> Self:
        """Read an Ensembl or GENCODE GTF file and open the genome file (2bit or FASTA).

        The GTF has no contig lengths, so the check of `Reference` covers only the contig names.
        The errors and the closing are those of `from_paths`. `load_gtf` describes what the reader
        takes from the GTF.
        """
        return cls._open(load_gtf, gtf, genome)

    @classmethod
    def from_directory(cls, path: str | os.PathLike[str]) -> Self:
        """Open the EnsDb and the genome of a directory, as `from_paths` does.

        The directory has one file that ends with `.sqlite`, and one genome file that ends with
        `.2bit` or `.fa.gz`. The index files of a FASTA file are not candidates, and the reader of
        the FASTA file raises the error for a missing index. A directory that is absent raises
        `FileNotFoundError`. A file that is not a directory, and an absent or a repeated candidate,
        raise `InputError`.
        """
        directory = local_path(path)
        found = os.fspath(directory)
        if not directory.exists():
            message = (
                f"{found!r}: no such directory; expected a directory with an EnsDb and a genome"
            )
            raise FileNotFoundError(message)
        if not directory.is_dir():
            message = f"{found!r}: not a directory; expected a directory with an EnsDb and a genome"
            raise InputError(message)
        names = sorted(entry.name for entry in directory.iterdir() if entry.is_file())
        ensdb = _only_candidate(directory, names, _ENSDB_SUFFIXES, "EnsDb")
        genome = _only_candidate(directory, names, _GENOME_SUFFIXES, "genome")
        return cls.from_paths(ensdb=ensdb, genome=genome)

    @classmethod
    def _open(
        cls,
        read_annotation: Callable[[str | os.PathLike[str]], AnnotationDb],
        annotation: str | os.PathLike[str],
        genome: str | os.PathLike[str],
    ) -> Self:
        with ExitStack() as stack:
            opened_annotation = read_annotation(annotation)
            stack.callback(opened_annotation.close)
            opened_genome = open_genome(genome)
            stack.callback(opened_genome.close)
            reference = cls(opened_annotation, opened_genome)
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

    def sequences(
        self,
        feature: FeatureType,
        transcript_filter: TranscriptFilter = TranscriptFilter(),
        *,
        per_gene: PerGene = PerGene.EACH_TRANSCRIPT,
        min_width: int = 8,
    ) -> SequenceSet:
        """Return the sequence of the feature of each transcript or gene that passes the filter.

        The ID of a sequence is its transcript ID, or its gene ID with `PerGene.UNION`. The gene ID
        of a sequence is the gene of the transcript. A sequence on the minus strand is the reverse
        complement of the genome. `min_width` removes each sequence that is shorter than the width,
        also an empty one. The default is the default of `load_annotations`. `InputError` is raised
        for a contig that the genome lacks.
        """
        width = as_limit("min_width", min_width)
        check_not_negative("min_width", width)
        features = self.features(feature, transcript_filter, per_gene=per_gene)
        return extract_sequences(features.at_least(width), self.genome)

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


def _only_candidate(
    directory: Path, names: list[str], suffixes: tuple[str, ...], kind: str
) -> Path:
    """Return the one file of the directory whose name ends with a suffix, from the sorted names."""
    candidates = [name for name in names if name.endswith(suffixes)]
    expected = f"expected one file that ends with {' or '.join(suffixes)}"
    if not candidates:
        present = f"these files: {listed_names(names)}" if names else "no files"
        message = (
            f"{os.fspath(directory)!r}: no {kind} file; {expected}; the directory has {present}"
        )
        raise InputError(message)
    if len(candidates) > 1:
        message = (
            f"{os.fspath(directory)!r}: more than one {kind} file ({', '.join(candidates)}); "
            f"{expected}"
        )
        raise InputError(message)
    return directory / candidates[0]
