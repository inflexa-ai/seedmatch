import hashlib
import os
import sqlite3
from pathlib import Path
from typing import NamedTuple

import numpy as np
import pandas as pd
import pytest
from annotation_cases import ALL, CASES
from helpers import golden

from seedmatch.annotation.filters import FeatureType, PerGene
from seedmatch.annotation.reference import Reference
from seedmatch.resources import Group, lookup

pytestmark = pytest.mark.refs


class Build(NamedTuple):
    gtf: str
    gtf_size: int
    gtf_sha256: str
    configs: tuple[str, ...]


RNOR7_CONFIGS = ("utr3", "utr5", "exons", "cds", "utr3_longest", "utr3_union")
# The registry has no GTF file, so this file pins the GTF of each build.
BUILDS = {
    "rnor7": Build(
        "Rattus_norvegicus.mRatBN7.2.109.gtf.gz",
        18029537,
        "ff7d5ecd69ed7df5e641ec90383e5de9e10e9d82316a0570c71fa80a98e8edca",
        RNOR7_CONFIGS,
    ),
    "hg38": Build(
        "Homo_sapiens.GRCh38.109.gtf.gz",
        54258835,
        "1ce4b88d9283f70220c318ebb3e654a9b5815316b9ba4c9b93564f29f10fab29",
        (*RNOR7_CONFIGS, "entrez", "tsl1"),
    ),
}


def pinned_files(build_id):
    twobit, ensdb = lookup(Group.REFERENCE, build_id).resources
    build = BUILDS[build_id]
    return {
        ensdb.name: (ensdb.size, ensdb.sha256),
        twobit.name: (twobit.size, twobit.sha256),
        build.gtf: (build.gtf_size, build.gtf_sha256),
    }


FILES = [(build_id, name) for build_id in BUILDS for name in pinned_files(build_id)]
CONFIGS = [(build_id, config) for build_id, build in BUILDS.items() for config in build.configs]


@pytest.fixture(scope="module")
def refs_dir():
    refs = os.environ.get("SEEDMATCH_TEST_REFS")
    if not refs:
        pytest.skip("SEEDMATCH_TEST_REFS is not set")
    return Path(refs)


# The GTF reference of hg38 needs about 1 GB to read. A module fixture with params groups its
# tests by build, so pytest closes the references of one build before it opens the next build.
@pytest.fixture(scope="module", params=list(BUILDS))
def build(request):
    return request.param


def present(directory, name):
    path = directory / name
    if not path.is_file():
        pytest.skip(f"{path} is absent")
    return path


def open_ensdb(refs_dir, build_id):
    directory = refs_dir / build_id
    twobit, ensdb = lookup(Group.REFERENCE, build_id).resources
    present(directory, ensdb.name)
    present(directory, twobit.name)
    return Reference.from_directory(directory)


@pytest.fixture(scope="module")
def ensdb_ref(refs_dir, build):
    with open_ensdb(refs_dir, build) as ref:
        yield ref


@pytest.fixture(scope="module")
def gtf_ref(refs_dir, build):
    directory = refs_dir / build
    twobit, _ = lookup(Group.REFERENCE, build).resources
    gtf = present(directory, BUILDS[build].gtf)
    with Reference.from_gtf(gtf, present(directory, twobit.name)) as ref:
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


@pytest.mark.parametrize(("build_id", "name"), FILES)
def test_files_are_the_pinned_files(refs_dir, build_id, name):
    path = present(refs_dir / build_id, name)
    size, sha256 = pinned_files(build_id)[name]
    assert path.stat().st_size == size
    with path.open("rb") as file:
        assert hashlib.file_digest(file, "sha256").hexdigest() == sha256


# The builds have different configs. pytest sets up an indirect fixture again for each pair of a
# build and a config, so this test opens the EnsDb itself.
@pytest.mark.parametrize(("build_id", "config"), CONFIGS)
def test_sequences_match_r(refs_dir, build_id, config):
    case = CASES[config]
    with open_ensdb(refs_dir, build_id) as ref:
        got = ref.sequences(
            case.feature, case.transcript_filter, per_gene=case.per_gene, min_width=case.min_width
        )
    expected = golden(f"annotation/{build_id}/{config}.tsv.gz")
    # For a union, R puts the gene ID in name_tx and NA in name_gene (D-002), so the gene IDs come
    # from name_tx.
    gene_ids = "name_tx" if case.per_gene is PerGene.UNION else "name_gene"
    assert list(got.ids) == expected["name_tx"].tolist()
    assert list(got.gene_ids) == expected[gene_ids].tolist()
    assert got.lengths().tolist() == expected["width"].astype(int).tolist()
    assert digest16(got) == expected["digest"].tolist()


@pytest.mark.parametrize("feature", list(FeatureType))
def test_gtf_gives_the_sequences_of_the_ensdb(ensdb_ref, gtf_ref, feature):
    assert_same_sequences(gtf_ref.sequences(feature, ALL), ensdb_ref.sequences(feature, ALL))


def test_gtf_union_gives_the_sequences_of_the_ensdb(ensdb_ref, gtf_ref):
    assert_same_sequences(
        gtf_ref.sequences(FeatureType.UTR3, ALL, per_gene=PerGene.UNION),
        ensdb_ref.sequences(FeatureType.UTR3, ALL, per_gene=PerGene.UNION),
    )


# The row order breaks a tie of the width and the start, and the GTF has another row order than the
# EnsDb. Thus a tie can select another transcript of the gene, and the test compares the widths.
def test_gtf_longest_gives_the_widths_of_the_ensdb(ensdb_ref, gtf_ref):
    expected = ensdb_ref.sequences(FeatureType.UTR3, ALL, per_gene=PerGene.LONGEST)
    got = gtf_ref.sequences(FeatureType.UTR3, ALL, per_gene=PerGene.LONGEST)
    assert len(expected) > 0
    assert got.gene_ids == expected.gene_ids
    np.testing.assert_array_equal(got.lengths(), expected.lengths())


# The row order breaks a tie of the contig and the start (D-030), and the GTF has another row order
# than the EnsDb. Thus the test compares the transcripts in the order of tx_id.
def test_gtf_gives_the_transcripts_of_the_ensdb(ensdb_ref, gtf_ref):
    expected = ensdb_ref.transcripts(ALL).sort_values("tx_id", ignore_index=True)
    got = gtf_ref.transcripts(ALL).sort_values("tx_id", ignore_index=True)
    assert len(expected) > 0
    pd.testing.assert_frame_equal(got, expected)


# The hg38 EnsDb also has the proteins of the patch scaffolds (CHR_*) and of the LRG transcripts,
# which the primary-assembly GTF lacks.
GTF_PROTEINS = {"rnor7": 45936, "hg38": 111317}


def test_gtf_gives_the_protein_genes_of_the_ensdb(refs_dir, build, ensdb_ref, gtf_ref):
    _, ensdb = lookup(Group.REFERENCE, build).resources
    uri = f"{(refs_dir / build / ensdb.name).resolve().as_uri()}?mode=ro&immutable=1"
    with sqlite3.connect(uri, uri=True) as connection:
        ids = [protein_id for (protein_id,) in connection.execute("SELECT protein_id FROM protein")]
    connection.close()
    expected = ensdb_ref.annotation.protein_genes(ids)
    got = gtf_ref.annotation.protein_genes(ids)
    assert len(got) == GTF_PROTEINS[build]
    assert got.items() <= expected.items()
