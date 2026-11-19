import hashlib
import os
from pathlib import Path

import numpy as np
import pytest
from annotation_cases import ALL, CASES
from helpers import golden

from seedmatch.annotation.filters import FeatureType, PerGene
from seedmatch.annotation.reference import Reference

pytestmark = pytest.mark.refs

ENSDB = "EnsDb.Rnorvegicus.v109.sqlite"
TWOBIT = "Rattus_norvegicus.mRatBN7.2.dna_sm.toplevel.2bit"
GTF = "Rattus_norvegicus.mRatBN7.2.109.gtf.gz"
RNOR7_FILES = {
    ENSDB: (179650560, "2453c59e1054157798df5a4114650911d046b5f7795af91a6a2171614ba0289b"),
    TWOBIT: (661989722, "73eacccd66a8843cced837d0cb0a272b42b7214b4a7c62a3c48963509ea914f3"),
    GTF: (18029537, "ff7d5ecd69ed7df5e641ec90383e5de9e10e9d82316a0570c71fa80a98e8edca"),
}
RNOR7_CONFIGS = ("utr3", "utr5", "exons", "cds", "utr3_longest", "utr3_union")


@pytest.fixture(scope="module")
def rnor7_dir():
    refs = os.environ.get("SEEDMATCH_TEST_REFS")
    if not refs:
        pytest.skip("SEEDMATCH_TEST_REFS is not set")
    return Path(refs) / "rnor7"


def present(directory, name):
    path = directory / name
    if not path.is_file():
        pytest.skip(f"{path} is absent")
    return path


@pytest.fixture(scope="module")
def rnor7(rnor7_dir):
    present(rnor7_dir, ENSDB)
    present(rnor7_dir, TWOBIT)
    with Reference.from_directory(rnor7_dir) as ref:
        yield ref


@pytest.fixture(scope="module")
def rnor7_gtf(rnor7_dir):
    with Reference.from_gtf(present(rnor7_dir, GTF), present(rnor7_dir, TWOBIT)) as ref:
        yield ref


def digest16(sequences):
    return [
        hashlib.sha256(sequences.sequence(i).encode("ascii")).hexdigest()[:16]
        for i in range(len(sequences))
    ]


def assert_same_sequences(got, expected):
    assert len(expected) > 0
    assert got.ids == expected.ids
    assert got.gene_ids == expected.gene_ids
    np.testing.assert_array_equal(got.offsets, expected.offsets)
    np.testing.assert_array_equal(got.codes, expected.codes)


@pytest.mark.parametrize("name", RNOR7_FILES)
def test_rnor7_files_are_the_pinned_files(rnor7_dir, name):
    path = present(rnor7_dir, name)
    size, sha256 = RNOR7_FILES[name]
    assert path.stat().st_size == size
    with path.open("rb") as file:
        assert hashlib.file_digest(file, "sha256").hexdigest() == sha256


@pytest.mark.parametrize("config", RNOR7_CONFIGS)
def test_rnor7_sequences_match_r(rnor7, config):
    case = CASES[config]
    got = rnor7.sequences(
        case.feature, case.transcript_filter, per_gene=case.per_gene, min_width=case.min_width
    )
    expected = golden(f"annotation/rnor7/{config}.tsv.gz")
    # For a union, R puts the gene ID in name_tx and NA in name_gene (D-002), so the gene IDs come
    # from name_tx.
    gene_ids = "name_tx" if case.per_gene is PerGene.UNION else "name_gene"
    assert list(got.ids) == expected["name_tx"].tolist()
    assert list(got.gene_ids) == expected[gene_ids].tolist()
    assert got.lengths().tolist() == expected["width"].astype(int).tolist()
    assert digest16(got) == expected["digest"].tolist()


@pytest.mark.parametrize("feature", list(FeatureType))
def test_rnor7_gtf_gives_the_sequences_of_the_ensdb(rnor7, rnor7_gtf, feature):
    assert_same_sequences(rnor7_gtf.sequences(feature, ALL), rnor7.sequences(feature, ALL))


def test_rnor7_gtf_union_gives_the_sequences_of_the_ensdb(rnor7, rnor7_gtf):
    assert_same_sequences(
        rnor7_gtf.sequences(FeatureType.UTR3, ALL, per_gene=PerGene.UNION),
        rnor7.sequences(FeatureType.UTR3, ALL, per_gene=PerGene.UNION),
    )


# The row order breaks a tie of the width and the start, and the GTF has another row order than the
# EnsDb. Thus a tie can select another transcript of the gene, and the test compares the widths.
def test_rnor7_gtf_longest_gives_the_widths_of_the_ensdb(rnor7, rnor7_gtf):
    expected = rnor7.sequences(FeatureType.UTR3, ALL, per_gene=PerGene.LONGEST)
    got = rnor7_gtf.sequences(FeatureType.UTR3, ALL, per_gene=PerGene.LONGEST)
    assert len(expected) > 0
    assert got.gene_ids == expected.gene_ids
    np.testing.assert_array_equal(got.lengths(), expected.lengths())
