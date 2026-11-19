import re
import shutil

import numpy as np
import pytest
from annotation_cases import ALL, ANNOTATION, CASES, ENSDB, copy_ensdb, fixture_reference
from helpers import golden

import seedmatch
import seedmatch.annotation
from seedmatch.annotation.extract import extract_sequences
from seedmatch.annotation.fasta_genome import FastaGenome
from seedmatch.annotation.features import Features
from seedmatch.annotation.filters import FeatureType, PerGene, TranscriptFilter
from seedmatch.annotation.reference import Reference
from seedmatch.annotation.twobit import TwoBitGenome
from seedmatch.errors import InputError
from seedmatch.sequence import reverse_complement

FA_GZ = ("synthetic.fa.gz", "synthetic.fa.gz.fai", "synthetic.fa.gz.gzi")


def directory_with(directory, *names):
    """Copy the fixture files of the names into the directory, and return the directory."""
    for name in names:
        shutil.copy(ANNOTATION / name, directory / name)
    return directory


@pytest.mark.parametrize("config", CASES)
def test_sequences_match_r(config):
    case = CASES[config]
    with fixture_reference() as ref:
        got = ref.sequences(
            case.feature, case.transcript_filter, per_gene=case.per_gene, min_width=case.min_width
        )
    expected = golden(f"annotation/{config}_sequences.tsv")
    # R gives NA to each name_gene of a union (D-002), so the gene IDs are the names of the genes.
    gene_names = "name_tx" if case.per_gene is PerGene.UNION else "name_gene"
    assert list(got.ids) == expected["name_tx"].tolist()
    assert list(got.gene_ids) == expected[gene_names].tolist()
    assert got.lengths().tolist() == expected["width"].astype(int).tolist()
    assert [got.sequence(i) for i in range(len(got))] == expected["sequence"].tolist()


def test_a_transcript_without_the_feature_has_an_empty_sequence_when_min_width_is_0():
    expected = golden("annotation/utr3_min0_sequences.tsv")
    empty = expected.loc[expected["width"] == "0", "name_tx"].tolist()
    assert empty
    with fixture_reference() as ref:
        got = ref.sequences(FeatureType.UTR3, ALL, min_width=0)
    assert [i for i, length in zip(got.ids, got.lengths(), strict=True) if length == 0] == empty


def test_min_width_keeps_a_sequence_of_exactly_that_width():
    with fixture_reference() as ref:
        every = ref.sequences(FeatureType.UTR3, ALL, min_width=0)
        widest = int(every.lengths().max())
        got = ref.sequences(FeatureType.UTR3, ALL, min_width=widest)
    kept = [i for i, length in zip(every.ids, every.lengths(), strict=True) if length == widest]
    assert 0 < len(kept) < len(every)
    assert list(got.ids) == kept


def test_the_defaults_of_sequences_are_the_defaults_of_load_annotations():
    expected = golden("annotation/utr3_sequences.tsv")
    with fixture_reference() as ref:
        got = ref.sequences(FeatureType.UTR3)
    assert list(got.ids) == expected["name_tx"].tolist()
    assert [got.sequence(i) for i in range(len(got))] == expected["sequence"].tolist()


def test_extract_joins_the_ranges_in_list_order_and_complements_each_range_on_the_minus_strand():
    features = Features(
        ids=("plus", "minus", "empty"),
        gene_ids=("g1", "g2", "g3"),
        contigs=("1", "X", "1"),
        strands=np.array([1, -1, 1], dtype=np.int8),
        starts=np.array([1101, 1001, 200, 100], dtype=np.int64),
        ends=np.array([1105, 1012, 209, 104], dtype=np.int64),
        offsets=np.array([0, 2, 4, 4], dtype=np.int64),
    )
    with TwoBitGenome(ANNOTATION / "synthetic.2bit") as genome:
        got = extract_sequences(features, genome)
        first_one, second_one = genome.sequence("1", 1101, 1105), genome.sequence("1", 1001, 1012)
        first_x, second_x = genome.sequence("X", 200, 209), genome.sequence("X", 100, 104)
    assert got.ids == ("plus", "minus", "empty")
    assert got.gene_ids == ("g1", "g2", "g3")
    assert got.sequence(0) == first_one + second_one
    assert got.sequence(1) == reverse_complement(first_x) + reverse_complement(second_x)
    assert got.sequence(2) == ""
    assert got.lengths().tolist() == [17, 15, 0]


def test_a_filter_without_a_transcript_gives_an_empty_set():
    with fixture_reference() as ref:
        got = ref.sequences(FeatureType.UTR3, TranscriptFilter(tx_ids=("absent",)))
    assert len(got) == 0
    assert len(got.codes) == 0


def test_union_sequences_keep_the_gene_id_D002():
    r_names = golden("annotation/utr3_union_r_sequences.tsv")
    assert set(r_names["name_gene"]) == {"NA"}
    case = CASES["utr3_union"]
    with fixture_reference() as ref:
        got = ref.sequences(case.feature, case.transcript_filter, per_gene=case.per_gene)
    assert got.gene_ids == got.ids
    assert list(got.gene_ids) == r_names["name_tx"].tolist()
    assert all(gene_id.startswith("ENSSYNG") for gene_id in got.gene_ids)


def test_sequence_ids_are_the_transcript_ids_D009():
    expected = golden("annotation/utr3_all_sequences.tsv")
    with fixture_reference() as ref:
        got = ref.sequences(FeatureType.UTR3, ALL)
        transcripts = ref.transcripts(ALL).set_index("tx_id")
    assert list(got.ids) == expected["name_tx"].tolist()
    assert set(got.ids) <= set(transcripts.index)
    assert list(got.gene_ids) == transcripts.loc[list(got.ids), "gene_id"].tolist()
    assert list(got.gene_ids) == expected["name_gene"].tolist()


def test_chr_names_of_the_annotation_find_the_contigs_of_the_genome(tmp_path):
    path = copy_ensdb(
        tmp_path, "UPDATE gene SET seq_name = 'chr' || seq_name WHERE seq_name != 'SCAF1.1'"
    )
    expected = golden("annotation/utr3_all_sequences.tsv")
    with Reference.from_paths(ensdb=path, genome=ANNOTATION / "synthetic.2bit") as ref:
        assert "chr1" in ref.annotation.contig_names()
        got = ref.sequences(FeatureType.UTR3, ALL)
    assert list(got.ids) == expected["name_tx"].tolist()
    assert [got.sequence(i) for i in range(len(got))] == expected["sequence"].tolist()


def test_contig_absent_from_the_genome_raises(tmp_path):
    path = copy_ensdb(
        tmp_path, "UPDATE gene SET seq_name = 'Y' WHERE gene_id = 'ENSSYNG00000000004'"
    )
    genome = ANNOTATION / "synthetic.2bit"
    message = (
        f"contig 'Y' of 'ENSSYNT00000000006' is not in the genome {str(genome)!r}; "
        "expected one of 1, SCAF1.1, X"
    )
    with (
        Reference.from_paths(ensdb=path, genome=genome) as ref,
        pytest.raises(InputError, match=re.escape(message)),
    ):
        ref.sequences(FeatureType.UTR3)


@pytest.mark.parametrize(
    ("min_width", "message"),
    [
        (-1, "min_width -1 is negative; expected 0 or more"),
        (1.5, "min_width 1.5 has type float; expected an integer"),
        (True, "min_width True has type bool; expected an integer"),
    ],
)
def test_min_width_rejects_a_value_that_is_not_a_width(min_width, message):
    with fixture_reference() as ref, pytest.raises(InputError, match=re.escape(message)):
        ref.sequences(FeatureType.UTR3, min_width=min_width)


def test_sequences_reject_a_feature_of_another_type():
    message = "feature '3UTR' has type str; expected a FeatureType"
    with fixture_reference() as ref, pytest.raises(InputError, match=re.escape(message)):
        ref.sequences("3UTR")


@pytest.mark.parametrize("genome", [("synthetic.2bit",), FA_GZ], ids=["2bit", "fa.gz"])
def test_from_directory_finds_the_ensdb_and_the_genome(tmp_path, genome):
    directory_with(tmp_path, "synthetic.sqlite", *genome)
    expected = golden("annotation/utr3_sequences.tsv")
    with Reference.from_directory(tmp_path) as ref:
        got = ref.sequences(FeatureType.UTR3)
    assert list(got.ids) == expected["name_tx"].tolist()
    assert [got.sequence(i) for i in range(len(got))] == expected["sequence"].tolist()


def test_from_directory_takes_a_str_and_gives_the_reader_of_the_file(tmp_path):
    directory_with(tmp_path, "synthetic.sqlite", "synthetic.2bit")
    with Reference.from_directory(str(tmp_path)) as ref:
        assert isinstance(ref.genome, TwoBitGenome)
    directory_with(tmp_path, *FA_GZ)
    (tmp_path / "synthetic.2bit").unlink()
    with Reference.from_directory(tmp_path) as ref:
        assert isinstance(ref.genome, FastaGenome)


def test_from_directory_ignores_the_files_that_are_not_candidates(tmp_path):
    directory_with(tmp_path, "synthetic.sqlite", "synthetic.fa", "synthetic.fa.fai", *FA_GZ)
    directory_with(tmp_path, "synthetic.gtf", "synthetic.gencode.gtf.gz")
    (tmp_path / "notes.sqlite.txt").write_text("not a candidate")
    (tmp_path / "folder.2bit").mkdir()
    with Reference.from_directory(tmp_path) as ref:
        assert isinstance(ref.genome, FastaGenome)


def test_from_directory_with_two_genomes_raises():
    with pytest.raises(InputError, match="more than one genome file") as raised:
        Reference.from_directory(ANNOTATION)
    message = str(raised.value)
    assert message.startswith(f"{str(ANNOTATION)!r}: ")
    assert "synthetic.2bit" in message
    assert "synthetic.fa.gz" in message
    assert ".fai" not in message
    assert message.endswith("expected one file that ends with .2bit or .fa.gz")


def test_from_directory_with_two_ensdbs_raises(tmp_path):
    directory_with(tmp_path, "synthetic.sqlite", "synthetic.2bit")
    shutil.copy(ENSDB, tmp_path / "another.sqlite")
    message = (
        f"{str(tmp_path)!r}: more than one EnsDb file (another.sqlite, synthetic.sqlite); "
        "expected one file that ends with .sqlite"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        Reference.from_directory(tmp_path)


def test_from_directory_without_an_ensdb_lists_the_files_of_the_directory(tmp_path):
    directory_with(tmp_path, "synthetic.2bit", "synthetic.gtf")
    message = (
        f"{str(tmp_path)!r}: no EnsDb file; expected one file that ends with .sqlite; "
        "the directory has these files: synthetic.2bit, synthetic.gtf"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        Reference.from_directory(tmp_path)


def test_from_directory_without_a_genome_lists_the_files_of_the_directory(tmp_path):
    directory_with(tmp_path, "synthetic.sqlite", "synthetic.fa", "synthetic.fa.fai")
    message = (
        f"{str(tmp_path)!r}: no genome file; expected one file that ends with .2bit or .fa.gz; "
        "the directory has these files: synthetic.fa, synthetic.fa.fai, synthetic.sqlite"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        Reference.from_directory(tmp_path)


def test_from_directory_without_a_file_of_either_kind_in_an_empty_directory_raises(tmp_path):
    message = (
        f"{str(tmp_path)!r}: no EnsDb file; expected one file that ends with .sqlite; "
        "the directory has no files"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        Reference.from_directory(tmp_path)


def test_from_directory_lists_only_the_first_files_of_a_large_directory(tmp_path):
    for number in range(1, 13):
        (tmp_path / f"file{number:02d}.txt").write_text("not a candidate")
    (tmp_path / "subdir").mkdir()
    listed = ", ".join(f"file{number:02d}.txt" for number in range(1, 11))
    message = f"the directory has these files: {listed}, ..."
    with pytest.raises(InputError, match=re.escape(message)) as raised:
        Reference.from_directory(tmp_path)
    assert "subdir" not in str(raised.value)


def test_from_directory_with_a_fa_gz_without_its_indexes_gives_the_error_of_the_reader(tmp_path):
    directory_with(tmp_path, "synthetic.sqlite", "synthetic.fa.gz")
    message = f"{str(tmp_path / 'synthetic.fa.gz.fai')!r}: the index file does not exist"
    with pytest.raises(FileNotFoundError, match=re.escape(message)):
        Reference.from_directory(tmp_path)


def test_from_directory_with_a_missing_directory_raises(tmp_path):
    message = f"{str(tmp_path / 'absent')!r}: no such directory"
    with pytest.raises(FileNotFoundError, match=re.escape(message)):
        Reference.from_directory(tmp_path / "absent")


def test_from_directory_with_a_file_raises():
    message = f"{str(ENSDB)!r}: not a directory; expected a directory with an EnsDb and a genome"
    with pytest.raises(InputError, match=re.escape(message)):
        Reference.from_directory(ENSDB)


@pytest.mark.parametrize(
    "name",
    ["Reference", "TranscriptFilter", "FeatureType", "PerGene", "Features"],
)
def test_the_annotation_api_is_in_the_top_level_package(name):
    assert name in seedmatch.__all__
    assert getattr(seedmatch, name) is getattr(seedmatch.annotation, name)


@pytest.mark.parametrize("name", ["TwoBitGenome", "FastaGenome", "Genome", "open_genome"])
def test_the_genome_api_is_in_the_annotation_package(name):
    assert name in seedmatch.annotation.__all__
    assert hasattr(seedmatch.annotation, name)
