import re
import sqlite3

import pandas as pd
import pytest
from annotation_cases import ANNOTATION, CASES, fixture_reference
from helpers import damaged_compressed_files

from seedmatch.annotation import gtf as gtf_module
from seedmatch.annotation.filters import TranscriptFilter
from seedmatch.annotation.gtf import load_gtf
from seedmatch.annotation.reference import Reference
from seedmatch.errors import InputError

GTF_CASES = {k: v for k, v in CASES.items() if v.transcript_filter.entrez_ids is None}
PAIRS = [("synthetic.gtf", "synthetic.2bit"), ("synthetic.gencode.gtf.gz", "synthetic.fa.gz")]
EVERYTHING = TranscriptFilter(canonical=False, protein_coding=False, standard_chromosomes=False)
# T10 is a transcript of the first gene row. T02 and T01 start at the same base, and their gene
# comes second.
D030_TIE = ["ENSSYNT00000000010", "ENSSYNT00000000002", "ENSSYNT00000000001"]
BIOTYPES = 'gene_biotype "protein_coding"; transcript_biotype "protein_coding";'
TX1 = f'gene_id "G1"; transcript_id "T1"; gene_name "Name"; {BIOTYPES}'
EXON1 = 'gene_id "G1"; transcript_id "T1"; exon_id "E1";'


def line(feature, start, end, attributes, *, strand="+"):
    return "\t".join(["1", "test", feature, str(start), str(end), ".", strand, ".", attributes])


def write_gtf(directory, *lines):
    path = directory / "test.gtf"
    path.write_text("\n".join(lines) + "\n")
    return path


def transcripts_of(directory, *lines):
    annotation = load_gtf(write_gtf(directory, *lines))
    try:
        return annotation.transcripts(EVERYTHING)
    finally:
        annotation.close()


def exons_of(directory, *lines):
    annotation = load_gtf(write_gtf(directory, *lines))
    try:
        return annotation.exon_rows(EVERYTHING, coding_only=False)
    finally:
        annotation.close()


def by_tx(frame):
    return frame.drop(columns="contig").sort_values("tx_id").reset_index(drop=True)


def tie_order(tx_ids):
    return [tx_id for tx_id in tx_ids if tx_id in D030_TIE]


@pytest.mark.parametrize(("gtf", "genome"), PAIRS)
@pytest.mark.parametrize("config", GTF_CASES)
def test_gtf_gives_the_sequences_of_the_ensdb(gtf, genome, config):
    case = GTF_CASES[config]
    options = {"per_gene": case.per_gene, "min_width": case.min_width}
    with (
        fixture_reference() as expected_ref,
        Reference.from_gtf(ANNOTATION / gtf, ANNOTATION / genome) as ref,
    ):
        expected = expected_ref.sequences(case.feature, case.transcript_filter, **options)
        got = ref.sequences(case.feature, case.transcript_filter, **options)
    assert len(expected) > 0
    assert got.ids == expected.ids
    assert got.gene_ids == expected.gene_ids
    assert [got.sequence(i) for i in range(len(got))] == [
        expected.sequence(i) for i in range(len(expected))
    ]


def test_ensembl_gtf_gives_the_transcripts_of_the_ensdb():
    with (
        fixture_reference() as expected_ref,
        Reference.from_gtf(ANNOTATION / "synthetic.gtf", ANNOTATION / "synthetic.2bit") as ref,
    ):
        expected = expected_ref.transcripts(EVERYTHING)
        got = ref.transcripts(EVERYTHING)
    assert len(got) == 12
    pd.testing.assert_frame_equal(got, expected)


def test_gencode_gtf_gives_the_transcripts_of_the_ensdb_with_the_contig_names_of_the_gtf():
    gtf = ANNOTATION / "synthetic.gencode.gtf.gz"
    with (
        fixture_reference() as expected_ref,
        Reference.from_gtf(gtf, ANNOTATION / "synthetic.fa.gz") as ref,
    ):
        expected = expected_ref.transcripts(EVERYTHING)
        got = ref.transcripts(EVERYTHING)
    assert len(got) == 12
    pd.testing.assert_frame_equal(by_tx(got), by_tx(expected))
    contigs = dict(zip(got["tx_id"], got["contig"], strict=True))
    assert contigs["ENSSYNT00000000001"] == "chr1"
    assert contigs["ENSSYNT00000000007"] == "chrX"
    assert contigs["ENSSYNT00000000006"] == "SCAF1.1"


@pytest.mark.parametrize(("gtf", "genome"), PAIRS)
def test_gtf_gives_the_tie_of_the_ensdb_for_each_filter_D030(gtf, genome):
    with Reference.from_gtf(ANNOTATION / gtf, ANNOTATION / genome) as ref:
        for name in ("exons_all", "exons_everything"):
            got = ref.transcripts(CASES[name].transcript_filter)
            assert tie_order(got["tx_id"]) == D030_TIE, name


def test_support_level_is_the_leading_integer():
    with Reference.from_gtf(ANNOTATION / "synthetic.gtf", ANNOTATION / "synthetic.2bit") as ref:
        got = ref.transcripts(EVERYTHING).set_index("tx_id")["support_level"]
    assert got["ENSSYNT00000000001"] == 1
    assert got["ENSSYNT00000000010"] == 2
    assert pd.isna(got["ENSSYNT00000000002"])


def test_entrez_filter_on_a_gtf_raises():
    gtf = ANNOTATION / "synthetic.gtf"
    message = (
        f"{str(gtf)!r}: the file has no table 'entrezgene'; "
        "expected a file with this table, or entrez_ids=None"
    )
    with (
        Reference.from_gtf(gtf, ANNOTATION / "synthetic.2bit") as ref,
        pytest.raises(InputError, match=re.escape(message)),
    ):
        ref.transcripts(TranscriptFilter(entrez_ids=(1001,)))


def test_gtf_without_an_agreeing_contig_raises(tmp_path):
    text = (ANNOTATION / "synthetic.gtf").read_text()
    renamed = [
        row if row.startswith("#") else "chrZ" + row[row.index("\t") :] for row in text.splitlines()
    ]
    gtf = write_gtf(tmp_path, *renamed)
    genome = ANNOTATION / "synthetic.2bit"
    message = (
        f"no contig of the annotation {str(gtf)!r} is in the genome {str(genome)!r}; "
        "annotation contigs: chrZ; genome contigs: 1, SCAF1.1, X"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        Reference.from_gtf(gtf, genome)


def test_gtf_without_a_transcript_raises_as_a_reference_without_genes(tmp_path):
    gtf = write_gtf(tmp_path, "# only a comment")
    message = f"the annotation {str(gtf)!r} has no genes; expected an annotation with at least one"
    with pytest.raises(InputError, match=re.escape(message)):
        Reference.from_gtf(gtf, ANNOTATION / "synthetic.2bit")


def test_from_gtf_takes_a_str_and_refuses_a_url():
    got = Reference.from_gtf(str(ANNOTATION / "synthetic.gtf"), str(ANNOTATION / "synthetic.2bit"))
    got.close()
    message = "'https://example.org/a.gtf': the path is a URL; expected a local path"
    with pytest.raises(InputError, match=re.escape(message)):
        Reference.from_gtf("https://example.org/a.gtf", ANNOTATION / "synthetic.2bit")


def test_a_missing_gtf_raises(tmp_path):
    path = tmp_path / "absent.gtf"
    message = f"{str(path)!r}: no such file; expected a GTF file"
    with pytest.raises(FileNotFoundError, match=re.escape(message)):
        load_gtf(path)


def test_a_gtf_that_is_not_utf8_raises(tmp_path):
    path = tmp_path / "bad.gtf"
    path.write_bytes(b"1\ttest\ttranscript\t\xff\n")
    with pytest.raises(InputError, match=re.escape(f"{str(path)!r}: the file is not valid UTF-8")):
        load_gtf(path)


def test_line_without_9_fields_raises_with_the_line_number(tmp_path):
    gtf = write_gtf(
        tmp_path,
        "# header",
        line("transcript", 100, 500, TX1),
        "1\ttest\texon\t100\t500\t.\t+\t.",
    )
    message = f"{str(gtf)!r}: line 3: the line has 8 tab-separated fields; expected 9"
    with pytest.raises(InputError, match=re.escape(message)):
        load_gtf(gtf)


def test_exon_without_a_transcript_line_raises(tmp_path):
    lines = [
        row
        for row in (ANNOTATION / "synthetic.gtf").read_text().splitlines()
        if not (
            row.split("\t")[2:3] == ["transcript"] and 'transcript_id "ENSSYNT00000000005"' in row
        )
    ]
    gtf = write_gtf(tmp_path, *lines)
    number = next(
        i for i, row in enumerate(lines, 1) if 'transcript_id "ENSSYNT00000000005"' in row
    )
    message = (
        f"{str(gtf)!r}: line {number}: the exon of transcript 'ENSSYNT00000000005' has no "
        "transcript line; expected a line with the feature transcript"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        load_gtf(gtf)


def test_a_repeated_transcript_line_raises(tmp_path):
    gtf = write_gtf(tmp_path, line("transcript", 100, 500, TX1), line("transcript", 100, 500, TX1))
    message = (
        f"{str(gtf)!r}: line 2: transcript 'T1' occurs twice; "
        "expected one transcript line for each transcript_id"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        load_gtf(gtf)


def test_an_exon_id_with_other_coordinates_raises(tmp_path):
    gtf = write_gtf(
        tmp_path,
        line("transcript", 100, 500, TX1),
        line("exon", 100, 500, EXON1),
        line("transcript", 100, 900, TX1.replace("T1", "T2")),
        line("exon", 100, 900, EXON1.replace("T1", "T2")),
    )
    message = (
        f"{str(gtf)!r}: exon 'E1' has the coordinates 100-500 at line 2 and 100-900 at line 4; "
        "expected the same coordinates for each line of an exon_id"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        load_gtf(gtf)


def test_transcripts_of_a_gene_share_an_exon_that_is_stored_once(tmp_path):
    second = TX1.replace("T1", "T2")
    rows = exons_of(
        tmp_path,
        line("transcript", 100, 900, TX1),
        line("exon", 100, 200, EXON1),
        line("exon", 800, 900, EXON1.replace("E1", "E2")),
        line("transcript", 100, 200, second),
        line("exon", 100, 200, EXON1.replace("T1", "T2")),
    )
    assert rows.tx_id.tolist() == ["T1", "T1", "T2"]
    assert rows.exon_start.tolist() == [100, 800, 100]
    assert rows.exon_end.tolist() == [200, 900, 200]


def test_exon_order_follows_the_strand_and_not_the_line_order(tmp_path):
    minus = TX1.replace("T1", "T2").replace("G1", "G2")
    rows = exons_of(
        tmp_path,
        line("transcript", 100, 900, TX1),
        line("exon", 800, 900, EXON1.replace("E1", "E2")),
        line("exon", 100, 200, EXON1),
        line("transcript", 100, 900, minus, strand="-"),
        line("exon", 100, 200, 'gene_id "G2"; transcript_id "T2"; exon_id "E3";'),
        line("exon", 800, 900, 'gene_id "G2"; transcript_id "T2"; exon_id "E4";'),
    )
    assert rows.tx_id.tolist() == ["T1", "T1", "T2", "T2"]
    assert rows.exon_start.tolist() == [100, 800, 800, 100]


def many_exons(count):
    lines = [line("transcript", 1, 10 * count, TX1, strand="-")]
    lines += [
        line("exon", 10 * i + 1, 10 * i + 5, f'gene_id "G1"; transcript_id "T1"; exon_id "E{i}";')
        for i in range(count)
    ]
    return lines


def test_a_transcript_with_more_exons_than_a_batch_keeps_each_exon_in_order(tmp_path):
    count = 10_005
    rows = exons_of(tmp_path, *many_exons(count))
    assert rows.exon_start.tolist() == [10 * i + 1 for i in reversed(range(count))]


def test_the_reader_writes_the_exon_lines_in_batches_of_10000(tmp_path, monkeypatch):
    batches = []
    reader_class = gtf_module._Reader  # noqa: SLF001 - the private reader does the batch writes
    flush = reader_class._flush  # noqa: SLF001 - the test wraps the batch write of the reader

    def counting_flush(reader):
        batches.append(len(reader._exons))  # noqa: SLF001 - the staged rows are the batch size
        flush(reader)

    monkeypatch.setattr(reader_class, "_flush", counting_flush)
    exons_of(tmp_path, *many_exons(25_000))
    assert batches == [10_000, 10_000, 5_000]


def test_a_comment_and_a_blank_line_are_skipped_in_the_middle_of_the_file(tmp_path):
    got = transcripts_of(
        tmp_path,
        "#!genome-build TEST",
        line("transcript", 100, 500, TX1),
        "# a comment",
        "",
        line("exon", 100, 500, EXON1),
    )
    assert got["tx_id"].tolist() == ["T1"]
    assert got["start"].tolist() == [100]
    assert got["end"].tolist() == [500]


def test_a_line_end_with_a_carriage_return_is_not_part_of_the_attributes(tmp_path):
    path = tmp_path / "crlf.gtf"
    path.write_bytes((line("transcript", 100, 500, TX1) + "\r\n").encode())
    annotation = load_gtf(path)
    got = annotation.transcripts(EVERYTHING)
    annotation.close()
    assert got["gene_biotype"].tolist() == ["protein_coding"]
    assert got["tx_biotype"].tolist() == ["protein_coding"]


def test_a_quoted_value_can_hold_a_semicolon(tmp_path):
    attributes = 'gene_id "G1"; transcript_id "T1"; gene_name "Syn; 1"; gene_biotype "lncRNA";'
    got = transcripts_of(tmp_path, line("transcript", 100, 500, attributes))
    assert got["gene_name"].tolist() == ["Syn; 1"]
    assert got["gene_biotype"].tolist() == ["lncRNA"]


def test_a_gene_without_a_name_has_an_empty_name(tmp_path):
    attributes = 'gene_id "G1"; transcript_id "T1"; gene_biotype "lncRNA";'
    got = transcripts_of(tmp_path, line("transcript", 100, 500, attributes))
    assert got["gene_name"].tolist() == [""]


@pytest.mark.parametrize(
    "tags", ['tag "Ensembl_canonical"; tag "basic";', 'tag "basic"; tag "Ensembl_canonical";']
)
def test_the_canonical_tag_can_be_any_of_a_repeated_tag(tmp_path, tags):
    other = TX1.replace("T1", "T2")
    got = transcripts_of(
        tmp_path,
        line("transcript", 100, 500, f"{TX1} {tags}"),
        line("transcript", 100, 500, f'{other} tag "basic";'),
    )
    assert got["canonical"].tolist() == [True, False]


def test_the_values_of_the_attributes_can_be_unquoted(tmp_path):
    attributes = f"gene_id G1; transcript_id T1; transcript_support_level 3; {BIOTYPES}"
    got = transcripts_of(tmp_path, line("transcript", 100, 500, attributes))
    assert got["tx_id"].tolist() == ["T1"]
    assert got["gene_id"].tolist() == ["G1"]
    assert got["support_level"].tolist() == [3]


def test_a_transcript_without_an_exon_has_no_exon_row_and_no_cds(tmp_path):
    got = transcripts_of(tmp_path, line("transcript", 100, 500, TX1))
    assert got["tx_id"].tolist() == ["T1"]
    assert got["cds_start"].isna().tolist() == [True]
    rows = exons_of(tmp_path, line("transcript", 100, 500, TX1))
    assert len(rows.tx_id) == 0


def test_a_cds_without_a_stop_codon_gives_the_range_of_the_cds_lines(tmp_path):
    cds = 'gene_id "G1"; transcript_id "T1";'
    got = transcripts_of(
        tmp_path,
        line("transcript", 100, 900, TX1),
        line("CDS", 300, 350, cds),
        line("CDS", 120, 200, cds),
    )
    assert got["cds_start"].tolist() == [120]
    assert got["cds_end"].tolist() == [350]


def test_the_stop_codon_extends_the_cds_and_a_start_codon_does_not(tmp_path):
    attributes = 'gene_id "G1"; transcript_id "T1";'
    got = transcripts_of(
        tmp_path,
        line("transcript", 100, 900, TX1),
        line("start_codon", 100, 102, attributes),
        line("CDS", 120, 200, attributes),
        line("stop_codon", 201, 203, attributes),
    )
    assert got["cds_start"].tolist() == [120]
    assert got["cds_end"].tolist() == [203]


def test_an_exon_line_before_its_transcript_line_is_read(tmp_path):
    rows = exons_of(
        tmp_path,
        line("exon", 100, 500, EXON1),
        line("transcript", 100, 500, TX1),
    )
    assert rows.tx_id.tolist() == ["T1"]
    assert rows.exon_start.tolist() == [100]


def test_the_lines_of_other_features_are_ignored(tmp_path):
    got = transcripts_of(
        tmp_path,
        line("gene", 100, 500, 'gene_id "G9";'),
        line("transcript", 100, 500, TX1),
        line("five_prime_utr", 100, 200, "gene_id; no pairs"),
        line("UTR", 100, 200, ""),
    )
    assert got["tx_id"].tolist() == ["T1"]


@pytest.mark.parametrize(
    ("row", "expected"),
    [
        (line("transcript", 100, 500, 'gene_id "G1"; gene_biotype "x";'), "has no transcript_id"),
        (line("transcript", 100, 500, 'transcript_id "T1";'), "has no gene_id"),
        (line("transcript", 100, 500, 'gene_id "G1"; transcript_id "";'), "has no transcript_id"),
        (line("transcript", "ten", 500, TX1), "start 'ten' and end '500'"),
        (line("transcript", 500, 100, TX1), "start '500' and end '100'"),
        (line("transcript", 0, 100, TX1), "start '0' and end '100'"),
        (line("transcript", 100, 500, TX1, strand="."), "strand '.'; expected + or -"),
        (line("exon", 100, 500, 'transcript_id "T1";'), "has no exon_id"),
        (line("exon", 100, 500, 'exon_id "E1";'), "has no transcript_id"),
        (line("CDS", 100, 500, 'gene_id "G1";'), "has no transcript_id"),
    ],
)
def test_a_line_with_a_bad_field_or_a_missing_attribute_raises_with_its_line_number(
    tmp_path, row, expected
):
    gtf = write_gtf(tmp_path, "# comment", row)
    with pytest.raises(InputError, match=re.escape(f"{str(gtf)!r}: line 2: ")) as raised:
        load_gtf(gtf)
    assert expected in str(raised.value)


@pytest.mark.parametrize("text", ["NA", "NA (assigned to previous version 3)"])
def test_a_support_level_that_starts_with_na_is_missing(tmp_path, text):
    attributes = f'{TX1} transcript_support_level "{text}";'
    got = transcripts_of(tmp_path, line("transcript", 100, 500, attributes))
    assert got["support_level"].isna().tolist() == [True]


@pytest.mark.parametrize("text", ["high", "NAN", ""])
def test_a_support_level_without_an_integer_or_na_raises(tmp_path, text):
    attributes = f'{TX1} transcript_support_level "{text}";'
    gtf = write_gtf(tmp_path, line("transcript", 100, 500, attributes))
    message = (
        f"{str(gtf)!r}: line 1: transcript_support_level {text!r}; "
        "expected a text that starts with NA or with an integer"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        load_gtf(gtf)


def test_a_gtf_without_the_canonical_tag_has_no_canonical_column(tmp_path):
    gtf = write_gtf(tmp_path, line("transcript", 100, 500, TX1))
    annotation = load_gtf(gtf)
    message = (
        f"{str(gtf)!r}: the table 'tx' has no column 'tx_is_canonical'; "
        "expected a file with this column, or canonical=False"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        annotation.transcripts(TranscriptFilter())
    got = annotation.transcripts(TranscriptFilter(canonical=False))
    annotation.close()
    assert got["canonical"].isna().tolist() == [True]


def test_a_gtf_without_a_support_level_has_no_support_level_column(tmp_path):
    gtf = write_gtf(tmp_path, line("transcript", 100, 500, TX1))
    annotation = load_gtf(gtf)
    message = (
        f"{str(gtf)!r}: the table 'tx' has no column 'tx_support_level'; "
        "expected a file with this column, or support_level=None"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        annotation.transcripts(TranscriptFilter(canonical=False, support_level=1))
    annotation.close()


def test_features_of_a_gtf_with_the_exons_of_a_minus_strand_gene_come_in_transcript_order(tmp_path):
    attributes = 'gene_id "G1"; transcript_id "T1";'
    gtf = write_gtf(
        tmp_path,
        line("transcript", 100, 900, TX1, strand="-"),
        line("exon", 100, 200, EXON1),
        line("exon", 800, 900, EXON1.replace("E1", "E2")),
        line("CDS", 150, 200, attributes),
        line("stop_codon", 147, 149, attributes),
    )
    annotation = load_gtf(gtf)
    rows = annotation.exon_rows(EVERYTHING, coding_only=True)
    annotation.close()
    assert rows.exon_start.tolist() == [800, 100]
    assert rows.cds_start.tolist() == [147, 147]
    assert rows.cds_end.tolist() == [200, 200]


@pytest.mark.parametrize(
    ("other", "contig", "strand"),
    [
        (
            line("transcript", 100, 500, TX1.replace("T1", "T2")).replace("1\ttest", "2\ttest", 1),
            "'2'",
            "'+'",
        ),
        (line("transcript", 100, 500, TX1.replace("T1", "T2"), strand="-"), "'1'", "'-'"),
    ],
    ids=["contig", "strand"],
)
def test_a_transcript_on_another_contig_or_strand_than_its_gene_raises(
    tmp_path, other, contig, strand
):
    gtf = write_gtf(tmp_path, line("transcript", 100, 500, TX1), other)
    message = (
        f"{str(gtf)!r}: line 2: transcript 'T2' is on contig {contig} and strand {strand}, "
        "but its gene 'G1' is on contig '1' and strand '+', from its first transcript line; "
        "expected one contig and one strand for each gene"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        load_gtf(gtf)


@pytest.mark.parametrize("name", damaged_compressed_files())
def test_a_damaged_compressed_gtf_raises_with_the_name_of_the_file(tmp_path, name):
    path = tmp_path / "damaged.gtf"
    path.write_bytes(damaged_compressed_files()[name][0])
    prefix = f"{str(path)!r}: the compressed data is damaged or truncated ("
    with pytest.raises(InputError, match=re.escape(prefix)):
        load_gtf(path)


def test_a_directory_as_a_gtf_raises(tmp_path):
    message = f"{str(tmp_path)!r}: the path is a directory; expected a file"
    with pytest.raises(InputError, match=re.escape(message)):
        load_gtf(tmp_path)


@pytest.mark.parametrize("fault", ["repeated transcript", "exon without a transcript"])
def test_load_gtf_closes_its_connection_when_it_raises(tmp_path, monkeypatch, fault):
    connections = []
    connect = sqlite3.connect

    def recording_connect(*args, **kwargs):
        connections.append(connect(*args, **kwargs))
        return connections[-1]

    monkeypatch.setattr(sqlite3, "connect", recording_connect)
    faults = {
        "repeated transcript": [line("transcript", 1, 9, TX1), line("transcript", 1, 9, TX1)],
        "exon without a transcript": [line("exon", 1, 9, EXON1)],
    }
    with pytest.raises(InputError):
        load_gtf(write_gtf(tmp_path, *faults[fault]))
    assert len(connections) == 1
    with pytest.raises(sqlite3.ProgrammingError, match="closed database"):
        connections[0].execute("SELECT 1")
