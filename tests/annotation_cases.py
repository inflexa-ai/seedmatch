import shutil
import sqlite3
from typing import NamedTuple

from helpers import DATA

from seedmatch.annotation.filters import FeatureType, PerGene, TranscriptFilter
from seedmatch.annotation.reference import Reference

ANNOTATION = DATA / "annotation"
ENSDB = ANNOTATION / "synthetic.sqlite"


class Case(NamedTuple):
    feature: FeatureType
    transcript_filter: TranscriptFilter
    per_gene: PerGene = PerGene.EACH_TRANSCRIPT
    min_width: int = 8


ALL = TranscriptFilter(canonical=False)
NOT_CODING = TranscriptFilter(canonical=False, protein_coding=False)

# The configs of tools/golden/R/annotation.R, without the configs whose names end with _r. The names
# of the goldens are tests/golden/annotation/<name>_features.tsv, and so on.
CASES: dict[str, Case] = {
    "utr3": Case(FeatureType.UTR3, TranscriptFilter()),
    "utr5": Case(FeatureType.UTR5, TranscriptFilter()),
    "exons": Case(FeatureType.EXONS, TranscriptFilter()),
    "cds": Case(FeatureType.CDS, TranscriptFilter()),
    "utr3_all": Case(FeatureType.UTR3, ALL),
    "utr3_min0": Case(FeatureType.UTR3, ALL, min_width=0),
    "exons_all": Case(FeatureType.EXONS, ALL),
    "exons_everything": Case(
        FeatureType.EXONS,
        TranscriptFilter(canonical=False, protein_coding=False, standard_chromosomes=False),
    ),
    "utr3_longest": Case(FeatureType.UTR3, ALL, PerGene.LONGEST),
    "utr3_longest_noncoding": Case(FeatureType.UTR3, NOT_CODING, PerGene.LONGEST),
    "utr3_union": Case(FeatureType.UTR3, ALL, PerGene.UNION),
    "utr3_union_min0": Case(FeatureType.UTR3, ALL, PerGene.UNION, min_width=0),
    "tsl1": Case(FeatureType.UTR3, TranscriptFilter(canonical=False, support_level=1)),
    "gene_ids": Case(
        FeatureType.UTR3,
        TranscriptFilter(canonical=False, gene_ids=("ENSSYNG00000000001", "ENSSYNG00000000002")),
    ),
    "tx_ids": Case(
        FeatureType.UTR3,
        TranscriptFilter(canonical=False, tx_ids=("ENSSYNT00000000002", "ENSSYNT00000000004")),
    ),
    "symbols": Case(FeatureType.UTR3, TranscriptFilter(canonical=False, symbols=("Syn1", "Syn5"))),
    "entrez": Case(
        FeatureType.UTR3, TranscriptFilter(canonical=False, entrez_ids=(1006, 1007, 1001))
    ),
    "nmd": Case(
        FeatureType.EXONS,
        TranscriptFilter(
            canonical=False, protein_coding=False, tx_biotypes=("nonsense_mediated_decay",)
        ),
    ),
    "lnc": Case(
        FeatureType.EXONS,
        TranscriptFilter(canonical=False, protein_coding=False, gene_biotypes=("lncRNA",)),
    ),
}


def fixture_reference() -> Reference:
    return Reference.from_paths(ensdb=ENSDB, genome=ANNOTATION / "synthetic.2bit")


def copy_ensdb(directory, *statements):
    """Copy the fixture EnsDb, run the statements on the copy, and return its path."""
    path = directory / "copy.sqlite"
    shutil.copy(ENSDB, path)
    with sqlite3.connect(path) as connection:
        for statement in statements:
            connection.execute(statement)
    connection.close()
    return path
