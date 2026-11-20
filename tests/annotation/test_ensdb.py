import json
import os
import re
import shutil
import sqlite3
import stat

import pandas as pd
import pytest
from annotation_cases import ANNOTATION, CASES, ENSDB, copy_ensdb, fixture_reference
from helpers import GOLDEN, golden

from seedmatch.annotation.ensdb import AnnotationDb
from seedmatch.annotation.fasta_genome import FastaGenome
from seedmatch.annotation.filters import STANDARD_CHROMOSOMES, TranscriptFilter
from seedmatch.annotation.genome import Genome
from seedmatch.annotation.reference import Reference
from seedmatch.errors import InputError

LENGTHS = {"1": 150002, "X": 5003, "SCAF1.1": 2001}
R_NAMES = {
    "tx_id": "tx_id",
    "gene_id": "gene_id",
    "seqnames": "contig",
    "start": "start",
    "end": "end",
    "strand": "strand",
    "tx_biotype": "tx_biotype",
    "tx_cds_seq_start": "cds_start",
    "tx_cds_seq_end": "cds_end",
    "tx_support_level": "support_level",
    "tx_is_canonical": "canonical",
}
INTEGERS = {"start", "end", "cds_start", "cds_end", "support_level"}
# T10 is a transcript of the gene whose row comes first in the gene table. T02 and T01 start at the
# same base, and their gene comes second.
D030_TIE = ["ENSSYNT00000000010", "ENSSYNT00000000002", "ENSSYNT00000000001"]
EVERYTHING_FLAGS = {"canonical": False, "protein_coding": False, "standard_chromosomes": False}
EVERYTHING = TranscriptFilter(**EVERYTHING_FLAGS)


class StubGenome(Genome):
    """A genome that has contigs and no bases. It counts the calls of `close`."""

    def __init__(self, lengths, path="stub.2bit"):
        super().__init__(path)
        self._lengths = lengths
        self.closed = 0

    @property
    def contigs(self):
        return self._lengths

    def _codes(self, contig, start, end):
        raise NotImplementedError

    def close(self):
        self.closed += 1


def write_genome(directory, lengths):
    """Write a plain FASTA file with one line for each contig, and its .fai. Return the path."""
    data = bytearray()
    fai = []
    for name, length in lengths.items():
        data += f">{name}\n".encode()
        fai.append(f"{name}\t{length}\t{len(data)}\t{length}\t{length + 1}\n")
        data += b"A" * length + b"\n"
    path = directory / "genome.fa"
    path.write_bytes(bytes(data))
    (directory / "genome.fa.fai").write_text("".join(fai))
    return path


def expected_value(column, text):
    if text == "NA":
        return None
    if column in INTEGERS:
        return int(text)
    if column == "canonical":
        return {"0": False, "1": True}[text]
    return text


def port_values(series):
    return [None if pd.isna(value) else value for value in series]


def in_port_order(expected):
    """Return the golden rows with the rows of the D-030 tie in the order of the port."""
    rows = list(range(len(expected)))
    tied = [row for row in rows if expected.loc[row, "tx_id"] in D030_TIE]
    in_order = sorted(tied, key=lambda row: D030_TIE.index(expected.loc[row, "tx_id"]))
    for slot, row in zip(tied, in_order, strict=True):
        rows[slot] = row
    return expected.iloc[rows].reset_index(drop=True)


def tie_order(tx_ids):
    return [tx_id for tx_id in tx_ids if tx_id in D030_TIE]


ALL_TRANSCRIPTS = in_port_order(golden("annotation/transcripts_exons_everything.tsv"))[
    "tx_id"
].tolist()


@pytest.mark.parametrize("config", CASES)
def test_transcripts_match_r(config):
    expected = golden(f"annotation/transcripts_{config}.tsv")
    if config == "exons_everything":
        expected = in_port_order(expected)  # the tie T10, T02, T01 (D-030)
    with fixture_reference() as ref:
        got = ref.transcripts(CASES[config].transcript_filter)
    assert len(got) == len(expected)
    for r_name, name in R_NAMES.items():
        wanted = [expected_value(name, text) for text in expected[r_name]]
        assert port_values(got[name]) == wanted, name


def test_cases_cover_each_golden_config_that_does_not_end_with_r():
    names = {path.stem for path in (GOLDEN / "annotation").glob("transcripts_*.tsv")}
    configs = {name.removeprefix("transcripts_") for name in names}
    assert set(CASES) == {config for config in configs if not config.endswith("_r")}


def test_reference_default_filter_is_the_default_of_load_annotations():
    expected = golden("annotation/transcripts_utr3.tsv")
    with fixture_reference() as ref:
        got = ref.transcripts()
    assert got["tx_id"].tolist() == expected["tx_id"].tolist()
    assert len(got) < len(ALL_TRANSCRIPTS)


def test_tie_order_does_not_depend_on_the_filter_D030():
    r_all = golden("annotation/transcripts_exons_all.tsv")["tx_id"]
    r_everything = golden("annotation/transcripts_exons_everything.tsv")["tx_id"]
    assert tie_order(r_all) == D030_TIE
    assert tie_order(r_everything) == ["ENSSYNT00000000002", "ENSSYNT00000000001", D030_TIE[0]]
    with fixture_reference() as ref:
        for name in ("exons_all", "exons_everything"):
            got = ref.transcripts(CASES[name].transcript_filter)
            assert tie_order(got["tx_id"]) == D030_TIE, name


def test_transcripts_have_the_columns_and_dtypes_of_the_port():
    with fixture_reference() as ref:
        got = ref.transcripts(EVERYTHING)
    assert {name: str(dtype) for name, dtype in got.dtypes.items()} == {
        "tx_id": "object",
        "gene_id": "object",
        "gene_name": "object",
        "contig": "object",
        "start": "int64",
        "end": "int64",
        "strand": "object",
        "tx_biotype": "object",
        "gene_biotype": "object",
        "cds_start": "Int64",
        "cds_end": "Int64",
        "support_level": "Int64",
        "canonical": "boolean",
    }
    assert got["tx_id"].tolist() == ALL_TRANSCRIPTS
    assert list(got.index) == list(range(len(got)))


def test_transcripts_give_the_names_and_the_biotypes_of_the_gene():
    with fixture_reference() as ref:
        got = ref.transcripts(EVERYTHING).set_index("tx_id")
    # The gene ENSSYNG00000000006 has the name '' in the EnsDb, as a gene without a name in rnor7.
    assert got["gene_name"].to_dict() == {
        "ENSSYNT00000000001": "Syn1",
        "ENSSYNT00000000002": "Syn1",
        "ENSSYNT00000000003": "Syn2",
        "ENSSYNT00000000004": "Syn2",
        "ENSSYNT00000000005": "Syn3",
        "ENSSYNT00000000006": "Syn4",
        "ENSSYNT00000000007": "Syn5",
        "ENSSYNT00000000008": "Syn5",
        "ENSSYNT00000000009": "",
        "ENSSYNT00000000010": "Syn7",
        "ENSSYNT00000000011": "Syn7",
        "ENSSYNT00000000012": "Syn8",
    }
    assert got.loc["ENSSYNT00000000005", "gene_biotype"] == "lncRNA"
    assert got.loc["ENSSYNT00000000011", "gene_biotype"] == "protein_coding"
    assert got.loc["ENSSYNT00000000011", "tx_biotype"] == "nonsense_mediated_decay"


def test_strand_is_the_strand_of_the_gene():
    with fixture_reference() as ref:
        got = ref.transcripts(EVERYTHING).set_index("tx_id")["strand"]
    assert got["ENSSYNT00000000003"] == "-"
    assert got["ENSSYNT00000000011"] == "-"
    assert got["ENSSYNT00000000001"] == "+"


def test_strand_that_is_not_plus_or_minus_raises(tmp_path):
    path = copy_ensdb(
        tmp_path, "UPDATE gene SET seq_strand = 0 WHERE gene_id = 'ENSSYNG00000000008'"
    )
    db = AnnotationDb.open(path)
    message = f"{os.fspath(path)!r}: transcript 'ENSSYNT00000000012' has strand 0; expected 1 or -1"
    with pytest.raises(InputError, match=re.escape(message)):
        db.transcripts(EVERYTHING)
    db.close()


def test_a_filter_without_a_match_gives_an_empty_frame_with_the_columns():
    with fixture_reference() as ref:
        got = ref.transcripts(TranscriptFilter(gene_ids=("none",)))
        everything = ref.transcripts(EVERYTHING)
    assert len(got) == 0
    assert list(got.columns) == list(everything.columns)
    assert got.dtypes.equals(everything.dtypes)


@pytest.mark.parametrize(
    ("transcript_filter", "expected"),
    [
        (
            TranscriptFilter(gene_ids=("ENSSYNG00000000007",), canonical=False),
            ["ENSSYNT00000000010"],
        ),
        (
            TranscriptFilter(
                gene_ids=("ENSSYNG00000000007",), canonical=False, protein_coding=False
            ),
            ["ENSSYNT00000000010", "ENSSYNT00000000011"],
        ),
        (
            TranscriptFilter(
                gene_ids=("ENSSYNG00000000007",),
                tx_biotypes=("nonsense_mediated_decay",),
                canonical=False,
                protein_coding=False,
            ),
            ["ENSSYNT00000000011"],
        ),
        (
            TranscriptFilter(
                gene_ids=("ENSSYNG00000000007",),
                canonical=False,
                protein_coding=False,
                support_level=2,
            ),
            ["ENSSYNT00000000010"],
        ),
        (
            TranscriptFilter(symbols=("Syn5", "Syn1"), entrez_ids=(1005,), canonical=False),
            ["ENSSYNT00000000007", "ENSSYNT00000000008"],
        ),
        (
            TranscriptFilter(
                gene_biotypes=("lncRNA", "protein_coding"),
                tx_biotypes=("lncRNA",),
                canonical=False,
                protein_coding=False,
            ),
            ["ENSSYNT00000000005"],
        ),
        (
            TranscriptFilter(
                standard_chromosomes=True, tx_ids=("ENSSYNT00000000006",), canonical=False
            ),
            [],
        ),
        (
            TranscriptFilter(
                standard_chromosomes=False, tx_ids=("ENSSYNT00000000006",), canonical=False
            ),
            ["ENSSYNT00000000006"],
        ),
    ],
    ids=[
        "gene",
        "gene without protein_coding",
        "gene+tx_biotype",
        "gene+support",
        "symbols+entrez",
        "biotypes",
        "std",
        "no std",
    ],
)
def test_fields_combine_with_and(transcript_filter, expected):
    with fixture_reference() as ref:
        got = ref.transcripts(transcript_filter)
    assert got["tx_id"].tolist() == expected


def test_require_tx_id_removes_the_transcript_without_an_id(tmp_path):
    path = copy_ensdb(tmp_path, "UPDATE tx SET tx_id = '' WHERE tx_id = 'ENSSYNT00000000012'")
    db = AnnotationDb.open(path)
    with_id = db.transcripts(TranscriptFilter(canonical=False))["tx_id"].tolist()
    without_check = db.transcripts(TranscriptFilter(canonical=False, require_tx_id=False))
    db.close()
    assert "" not in with_id
    assert len(with_id) == len(without_check) - 1
    assert "" in without_check["tx_id"].tolist()


def test_a_filter_without_a_condition_selects_each_transcript():
    off = {"standard_chromosomes": False, "require_tx_id": False, "protein_coding": False}
    with fixture_reference() as ref:
        got = ref.transcripts(TranscriptFilter(canonical=False, **off))
    assert got["tx_id"].tolist() == ALL_TRANSCRIPTS


def test_where_gives_one_condition_for_each_field_and_the_values_as_parameters():
    off = {"standard_chromosomes": False, "require_tx_id": False, "protein_coding": False}
    transcript_filter = TranscriptFilter(
        support_level=2, gene_ids=("G1", "G2"), entrez_ids=(7,), symbols=("S",), **off
    )
    db = AnnotationDb.open(ENSDB)
    condition, parameters = db.where(transcript_filter)
    db.close()
    assert condition == (
        "(tx.tx_is_canonical = 1) AND (tx.tx_support_level = ?)"
        " AND (tx.gene_id IN (SELECT value FROM json_each(?)))"
        " AND (gene.gene_name IN (SELECT value FROM json_each(?)))"
        " AND (tx.gene_id IN (SELECT gene_id FROM entrezgene"
        " WHERE entrezid IN (SELECT value FROM json_each(?))))"
    )
    assert parameters == [2, '["G1", "G2"]', '["S"]', "[7]"]


def test_where_of_the_defaults_filters_on_the_standard_chromosomes_as_one_parameter():
    db = AnnotationDb.open(ENSDB)
    condition, parameters = db.where(TranscriptFilter())
    db.close()
    assert condition == (
        "(gene.seq_name IN (SELECT value FROM json_each(?))) AND (tx.tx_id != '')"
        " AND (gene.gene_biotype = 'protein_coding' AND tx.tx_biotype = 'protein_coding')"
        " AND (tx.tx_is_canonical = 1)"
    )
    assert len(parameters) == 1
    assert set(json.loads(parameters[0])) == STANDARD_CHROMOSOMES


def test_filter_with_more_ids_than_the_parameter_limit_keeps_the_correct_transcripts():
    fixture_genes = [f"ENSSYNG{number:011d}" for number in range(1, 9)]
    others = [f"ENSSYNG9{number:010d}" for number in range(39992)]
    ids = [*others[:20000], *fixture_genes, *others[20000:]]
    assert len(ids) == 40000
    assert len(ids) > 32766  # the default limit of the bound parameters of SQLite
    with fixture_reference() as ref:
        got = ref.transcripts(TranscriptFilter(gene_ids=ids, **EVERYTHING_FLAGS))
    assert got["tx_id"].tolist() == ALL_TRANSCRIPTS


def test_filter_with_many_entrez_ids_keeps_the_correct_transcripts():
    ids = [*range(2000, 42000), 1006, 1007]
    with fixture_reference() as ref:
        got = ref.transcripts(TranscriptFilter(entrez_ids=ids, **EVERYTHING_FLAGS))
    assert got["tx_id"].tolist() == ["ENSSYNT00000000009"]


def test_a_value_with_a_quote_is_a_value_and_not_sql():
    with fixture_reference() as ref:
        got = ref.transcripts(TranscriptFilter(symbols=("x') OR ('1' = '1",), **EVERYTHING_FLAGS))
    assert len(got) == 0


def test_contig_lengths_and_names():
    db = AnnotationDb.open(ENSDB)
    lengths = db.contig_lengths()
    names = db.contig_names()
    db.close()
    assert lengths == LENGTHS
    assert names == ["1", "SCAF1.1", "X"]


def test_contig_lengths_without_the_chromosome_table_is_empty(tmp_path):
    path = copy_ensdb(tmp_path, "DROP TABLE chromosome")
    db = AnnotationDb.open(path)
    lengths = db.contig_lengths()
    db.close()
    assert lengths == {}


def test_a_file_without_the_optional_columns_gives_missing_values(tmp_path):
    path = copy_ensdb(
        tmp_path,
        "ALTER TABLE tx DROP COLUMN tx_is_canonical",
        "ALTER TABLE tx DROP COLUMN tx_support_level",
    )
    db = AnnotationDb.open(path)
    got = db.transcripts(TranscriptFilter(canonical=False))
    db.close()
    assert len(got) > 0
    assert got["support_level"].isna().all()
    assert got["canonical"].isna().all()
    assert str(got["support_level"].dtype) == "Int64"
    assert str(got["canonical"].dtype) == "boolean"


@pytest.mark.parametrize(
    ("column", "transcript_filter"),
    [
        ("tx_is_canonical", TranscriptFilter()),
        ("tx_support_level", TranscriptFilter(canonical=False, support_level=1)),
    ],
)
def test_a_filter_that_needs_an_absent_column_raises(tmp_path, column, transcript_filter):
    path = copy_ensdb(tmp_path, f"ALTER TABLE tx DROP COLUMN {column}")
    db = AnnotationDb.open(path)
    message = f"{os.fspath(path)!r}: the table 'tx' has no column {column!r}"
    with pytest.raises(InputError, match=re.escape(message)):
        db.transcripts(transcript_filter)
    db.close()


def test_canonical_on_a_file_without_the_column_raises(tmp_path):
    path = copy_ensdb(tmp_path, "ALTER TABLE tx DROP COLUMN tx_is_canonical")
    db = AnnotationDb.open(path)
    with pytest.raises(InputError, match=r"tx_is_canonical.*canonical=False"):
        db.where(TranscriptFilter(canonical=True))
    db.close()


def test_entrez_on_a_file_without_the_table_raises(tmp_path):
    path = copy_ensdb(tmp_path, "DROP TABLE entrezgene")
    db = AnnotationDb.open(path)
    message = f"{os.fspath(path)!r}: the file has no table 'entrezgene'"
    with pytest.raises(InputError, match=re.escape(message)):
        db.transcripts(TranscriptFilter(canonical=False, entrez_ids=(1006,)))
    db.close()


def test_file_that_is_not_sqlite_raises(tmp_path):
    path = tmp_path / "notes.sqlite"
    path.write_text("this is not an SQLite file. " * 100)
    message = f"{os.fspath(path)!r}: SQLite cannot read the file (file is not a database)"
    with pytest.raises(InputError, match=re.escape(message)):
        AnnotationDb.open(path)


def test_a_directory_is_not_an_ensdb(tmp_path):
    # SQLite words the reason by platform, so the test matches the part that the port writes.
    message = f"{os.fspath(tmp_path)!r}: SQLite cannot read the file ("
    with pytest.raises(InputError, match=re.escape(message)):
        AnnotationDb.open(tmp_path)


def test_missing_ensdb_raises_file_not_found(tmp_path):
    # With mode=ro, SQLite says only "unable to open database file".
    path = tmp_path / "absent.sqlite"
    with pytest.raises(FileNotFoundError, match=r"absent\.sqlite"):
        AnnotationDb.open(path)
    assert not path.exists()


def test_ensdb_url_raises():
    with pytest.raises(InputError, match="the path is a URL"):
        AnnotationDb.open("https://example.org/EnsDb.sqlite")


@pytest.mark.parametrize("table", ["gene", "tx", "tx2exon", "exon"])
def test_sqlite_file_without_a_required_table_raises(tmp_path, table):
    path = copy_ensdb(tmp_path, f"DROP TABLE {table}")
    message = f"{os.fspath(path)!r}: the file has no table {table!r}; expected the EnsDb tables"
    with pytest.raises(InputError, match=re.escape(message)):
        AnnotationDb.open(path)


REQUIRED_COLUMNS = {
    "gene": ["gene_id", "gene_name", "gene_biotype", "seq_name", "seq_strand"],
    "tx": [
        "tx_id",
        "tx_biotype",
        "tx_seq_start",
        "tx_seq_end",
        "tx_cds_seq_start",
        "tx_cds_seq_end",
        "gene_id",
    ],
    "tx2exon": ["tx_id", "exon_id", "exon_idx"],
    "exon": ["exon_id", "exon_seq_start", "exon_seq_end"],
}
OPTIONAL_COLUMNS = {"chromosome": ["seq_name", "seq_length"], "entrezgene": ["gene_id", "entrezid"]}
EACH_FIELD = TranscriptFilter(
    support_level=1,
    gene_ids=("ENSSYNG00000000001",),
    tx_ids=("ENSSYNT00000000001",),
    symbols=("Syn1",),
    tx_biotypes=("protein_coding",),
    gene_biotypes=("protein_coding",),
    entrez_ids=(1001,),
)


def test_the_queries_read_only_the_required_columns_and_the_optional_columns():
    connection = sqlite3.connect(f"{ENSDB.resolve().as_uri()}?mode=ro", uri=True)
    db = AnnotationDb(connection, "fixture")
    read = set()

    def record(action, table, column, *_):
        if action == sqlite3.SQLITE_READ:
            read.add((table, column))
        return sqlite3.SQLITE_OK

    connection.set_authorizer(record)
    db.contig_names()
    db.contig_lengths()
    for transcript_filter in (EACH_FIELD, EVERYTHING):
        db.transcripts(transcript_filter)
        db.exon_rows(transcript_filter, coding_only=True)
    db.close()
    # Each table of an EnsDb has a rowid, and the optional tx columns have their own check.
    optional = {("tx", "tx_is_canonical"), ("tx", "tx_support_level")}
    tables = {**REQUIRED_COLUMNS, **OPTIONAL_COLUMNS}
    assert {
        (table, column)
        for table, column in read - optional
        if table in tables and column != "ROWID"
    } == {(table, column) for table, columns in tables.items() for column in columns}


@pytest.mark.parametrize(
    ("table", "column"),
    [(table, column) for table, columns in OPTIONAL_COLUMNS.items() for column in columns],
)
def test_an_optional_table_without_a_column_that_the_queries_read_raises_at_open(
    tmp_path, table, column
):
    path = copy_ensdb(tmp_path, f"ALTER TABLE {table} RENAME COLUMN {column} TO renamed")
    message = (
        f"{os.fspath(path)!r}: the table {table!r} has no column {column!r}; "
        f"expected the EnsDb columns {', '.join(OPTIONAL_COLUMNS[table])}"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        AnnotationDb.open(path)


@pytest.mark.parametrize(
    ("table", "column"),
    [(table, column) for table, columns in REQUIRED_COLUMNS.items() for column in columns],
)
def test_a_required_table_without_a_column_that_the_queries_read_raises_at_open(
    tmp_path, table, column
):
    path = copy_ensdb(tmp_path, f"ALTER TABLE {table} RENAME COLUMN {column} TO renamed")
    message = (
        f"{os.fspath(path)!r}: the table {table!r} has no column {column!r}; "
        f"expected the EnsDb columns {', '.join(REQUIRED_COLUMNS[table])}"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        AnnotationDb.open(path)


@pytest.mark.parametrize(
    ("change", "problem"),
    [
        ("tx_seq_start = NULL", "tx.tx_seq_start None; expected an integer from 1 to 2^62 - 1"),
        ("tx_seq_start = 0", "tx.tx_seq_start 0; expected an integer from 1 to 2^62 - 1"),
        ("tx_seq_end = 'x'", "tx.tx_seq_end 'x'; expected an integer from 1 to 2^62 - 1"),
        ("tx_seq_end = 1.5", "tx.tx_seq_end 1.5; expected an integer from 1 to 2^62 - 1"),
        (
            f"tx_seq_end = {2**62}",
            f"tx.tx_seq_end {2**62}; expected an integer from 1 to 2^62 - 1",
        ),
        (
            "tx_cds_seq_start = 'x'",
            "tx.tx_cds_seq_start 'x'; expected NULL or an integer from 0 to 2^62 - 1",
        ),
        (
            f"tx_cds_seq_end = {2**62}",
            f"tx.tx_cds_seq_end {2**62}; expected NULL or an integer from 0 to 2^62 - 1",
        ),
        ("tx_support_level = 'x'", "tx.tx_support_level 'x'; expected NULL or an integer"),
        ("tx_support_level = 1.5", "tx.tx_support_level 1.5; expected NULL or an integer"),
        (
            "tx_seq_start = 81601",
            "start 81601 and end 81600; expected integers with start <= end",
        ),
    ],
)
def test_transcripts_with_a_number_that_is_not_valid_raise(tmp_path, change, problem):
    # The statement has a change that this test fixes, so no input reaches it.
    path = copy_ensdb(tmp_path, f"UPDATE tx SET {change} WHERE tx_id = 'ENSSYNT00000000012'")  # noqa: S608
    message = f"{os.fspath(path)!r}: transcript 'ENSSYNT00000000012' has {problem}"
    db = AnnotationDb.open(path)
    with pytest.raises(InputError, match=re.escape(message) + "$"):
        db.transcripts(EVERYTHING)
    db.close()


@pytest.mark.parametrize("level", [-3, 0, 99, 2**62])
def test_transcripts_keep_a_support_level_outside_the_range_of_the_filter(tmp_path, level):
    # The statement has a level that this test fixes, so no input reaches it.
    change = f"UPDATE tx SET tx_support_level = {level} WHERE tx_id = 'ENSSYNT00000000012'"  # noqa: S608
    path = copy_ensdb(tmp_path, change)
    db = AnnotationDb.open(path)
    frame = db.transcripts(EVERYTHING).set_index("tx_id")
    db.close()
    assert frame.loc["ENSSYNT00000000012", "support_level"] == level


def test_sqlite_file_without_the_tx_table_raises(tmp_path):
    path = tmp_path / "other.sqlite"
    with sqlite3.connect(path) as connection:
        connection.execute("CREATE TABLE gene (gene_id TEXT)")
    connection.close()
    with pytest.raises(InputError, match="table 'tx'"):
        AnnotationDb.open(path)


def test_ensdb_with_a_space_and_a_hash_in_its_name_opens(tmp_path):
    path = tmp_path / "an ensdb #1?.sqlite"
    shutil.copy(ENSDB, path)
    db = AnnotationDb.open(path)
    count = len(db.transcripts(EVERYTHING))
    db.close()
    assert count == len(ALL_TRANSCRIPTS)


def test_reference_opens_from_a_read_only_directory(tmp_path):
    directory = tmp_path / "refs"
    directory.mkdir()
    for name in ("synthetic.sqlite", "synthetic.2bit"):
        shutil.copy(ANNOTATION / name, directory / name)
        (directory / name).chmod(stat.S_IRUSR)
    before = sorted(path.name for path in directory.iterdir())
    directory.chmod(stat.S_IRUSR | stat.S_IXUSR)
    try:
        with Reference.from_paths(
            ensdb=directory / "synthetic.sqlite", genome=directory / "synthetic.2bit"
        ) as ref:
            got = ref.transcripts(EVERYTHING)
            ref.genome.sequence("1", 1, 8)
        after = sorted(path.name for path in directory.iterdir())
    finally:
        directory.chmod(stat.S_IRWXU)
    assert got["tx_id"].tolist() == ALL_TRANSCRIPTS
    assert after == before


def test_reference_from_paths_gives_the_annotation_and_the_genome():
    with fixture_reference() as ref:
        assert isinstance(ref.annotation, AnnotationDb)
        assert dict(ref.genome.contigs) == LENGTHS


def test_with_statement_closes_the_annotation_and_the_genome():
    genome = StubGenome(LENGTHS)
    ref = Reference(AnnotationDb.open(ENSDB), genome)
    with ref:
        pass
    assert genome.closed == 1
    with pytest.raises(sqlite3.ProgrammingError, match="closed"):
        ref.transcripts()


def test_close_closes_the_genome_even_if_the_annotation_fails(monkeypatch):
    genome = StubGenome(LENGTHS)
    ref = Reference(AnnotationDb.open(ENSDB), genome)
    annotation_close = AnnotationDb.close

    def failing(self):
        annotation_close(self)
        message = "cannot close"
        raise RuntimeError(message)

    monkeypatch.setattr(AnnotationDb, "close", failing)
    with pytest.raises(RuntimeError, match="cannot close"):
        ref.close()
    assert genome.closed == 1


def test_genome_with_another_length_raises(tmp_path):
    genome = write_genome(tmp_path, {**LENGTHS, "X": 4000})
    message = r"contig 'X': the genome .* has length 4000, but the annotation .* has length 5003"
    with pytest.raises(InputError, match=message):
        Reference.from_paths(ensdb=ENSDB, genome=genome)


def test_genome_with_chr_names_and_the_same_lengths_agrees():
    genome = StubGenome({"chr1": 150002, "chrX": 5003})
    with Reference(AnnotationDb.open(ENSDB), genome) as ref:
        assert ref.genome is genome


def test_genome_with_chr_names_and_another_length_names_both_contigs():
    genome = StubGenome({"chr1": 150002, "chrX": 4000})
    message = r"contig 'X' \(genome contig 'chrX'\): .* has length 4000, but .* has length 5003"
    with pytest.raises(InputError, match=message):
        Reference(AnnotationDb.open(ENSDB), genome)


def test_a_contig_of_the_annotation_that_the_genome_lacks_does_not_raise():
    with Reference(AnnotationDb.open(ENSDB), StubGenome({"1": 150002})) as ref:
        assert ref.transcripts(EVERYTHING)["tx_id"].tolist() == ALL_TRANSCRIPTS


def test_annotation_without_lengths_agrees_with_any_genome(tmp_path):
    path = copy_ensdb(tmp_path, "DROP TABLE chromosome")
    with Reference(AnnotationDb.open(path), StubGenome({"1": 7})) as ref:
        assert ref.annotation.contig_lengths() == {}


def test_annotation_without_genes_raises(tmp_path):
    path = copy_ensdb(tmp_path, "DELETE FROM gene")
    db = AnnotationDb.open(path)
    message = f"the annotation {os.fspath(path)!r} has no genes; expected an annotation with at"
    with pytest.raises(InputError, match=re.escape(message)):
        Reference(db, StubGenome(LENGTHS))
    db.close()


def test_genome_without_an_agreeing_contig_raises(tmp_path):
    genome = write_genome(tmp_path, {"2": 10, "other": 5})
    message = (
        r"no contig of the annotation .* is in the genome .*; "
        r"annotation contigs: 1, SCAF1.1, X; genome contigs: 2, other"
    )
    with pytest.raises(InputError, match=message):
        Reference.from_paths(ensdb=ENSDB, genome=genome)


def test_genome_without_an_agreeing_contig_lists_some_names_of_each(tmp_path):
    lengths = {f"contig{number:02d}": 5 for number in range(1, 13)}
    genome = write_genome(tmp_path, lengths)
    with pytest.raises(InputError) as raised:
        Reference.from_paths(ensdb=ENSDB, genome=genome)
    listed = ", ".join(sorted(lengths)[:10])
    assert f"genome contigs: {listed}, ..." in str(raised.value)


@pytest.mark.parametrize("failure", ["length", "format"])
def test_from_paths_closes_what_it_opened_when_it_fails(tmp_path, monkeypatch, failure):
    if failure == "length":
        genome = write_genome(tmp_path, {**LENGTHS, "X": 4000})
    else:
        genome = tmp_path / "notes.txt"
        genome.write_text("not a genome file")
    closed = []
    annotation_close = AnnotationDb.close
    genome_close = FastaGenome.close

    def record_annotation(self):
        closed.append("annotation")
        annotation_close(self)

    def record_genome(self):
        closed.append("genome")
        genome_close(self)

    monkeypatch.setattr(AnnotationDb, "close", record_annotation)
    monkeypatch.setattr(FastaGenome, "close", record_genome)
    messages = {
        "length": "contig 'X': the genome .* has length 4000",
        "format": r"notes\.txt': not a 2bit or FASTA file",
    }
    with pytest.raises(InputError, match=messages[failure]):
        Reference.from_paths(ensdb=ENSDB, genome=genome)
    assert closed == (["genome", "annotation"] if failure == "length" else ["annotation"])


def test_from_paths_closes_the_annotation_when_the_genome_is_missing(tmp_path, monkeypatch):
    closed = []
    annotation_close = AnnotationDb.close

    def record_annotation(self):
        closed.append("annotation")
        annotation_close(self)

    monkeypatch.setattr(AnnotationDb, "close", record_annotation)
    with pytest.raises(FileNotFoundError):
        Reference.from_paths(ensdb=ENSDB, genome=tmp_path / "absent.2bit")
    assert closed == ["annotation"]


def test_a_result_of_many_chunks_keeps_each_row_in_order(tmp_path):
    copies = 120000
    # The statement has a number that this test fixes, so no input reaches it.
    path = copy_ensdb(
        tmp_path,
        f"""WITH RECURSIVE n(i) AS (SELECT 1 UNION ALL SELECT i + 1 FROM n WHERE i < {copies})
        INSERT INTO tx SELECT 'COPY' || i, tx_biotype, tx_seq_start, tx_seq_end, tx_cds_seq_start,
            tx_cds_seq_end, gene_id, tx_support_level, tx_id_version, gc_content,
            tx_external_name, tx_is_canonical FROM tx, n WHERE tx.tx_id = 'ENSSYNT00000000012'""",  # noqa: S608
    )
    db = AnnotationDb.open(path)
    got = db.transcripts(EVERYTHING)
    db.close()
    # The copies have the gene, the start, and so the place of ENSSYNT00000000012. They follow it
    # in the order of their rows.
    at = ALL_TRANSCRIPTS.index("ENSSYNT00000000012") + 1
    wanted = [
        *ALL_TRANSCRIPTS[:at],
        *(f"COPY{number}" for number in range(1, copies + 1)),
        *ALL_TRANSCRIPTS[at:],
    ]
    assert got["tx_id"].tolist() == wanted
    assert str(got["support_level"].dtype) == "Int64"
    assert str(got["canonical"].dtype) == "boolean"
    assert list(got.index) == list(range(len(wanted)))
