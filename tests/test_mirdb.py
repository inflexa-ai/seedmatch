import gzip
import os
import re

import pandas as pd
import pytest
from helpers import (
    DATA,
    NOT_UTF8,
    OBJECT_TEXT_ONLY,
    TEXT_DTYPE,
    assert_frame_matches_golden,
    golden_note,
)

import seedmatch
from seedmatch.errors import InputError
from seedmatch.mirdb import read_mirdb

SAMPLE = DATA / "mirdb_sample.txt"
LINE_FORMAT = "the miRDB ID, the RefSeq ID, and the score"
R_NAMES = {
    "mirdb_id": "miRDB.ID",
    "target_refseq_id": "target.REFSEQ.ID",
    "target_score": "target.score",
    "mirna_name": "miRNA.name",
}


def write_lines(tmp_path, lines, name="x.txt", ending="\n"):
    path = tmp_path / name
    # surrogateescape writes a "\udcff" in a line as the byte 0xFF, which is not UTF-8.
    path.write_bytes(
        "".join(f"{line}{ending}" for line in lines).encode("utf-8", "surrogateescape")
    )
    return path


def failure(path, text):
    return f"^{re.escape(f'{os.fspath(path)!r}: {text}')}$"


def test_read_mirdb_matches_r():
    got = read_mirdb(SAMPLE)
    assert list(got.columns) == [
        "mirdb_id",
        "target_refseq_id",
        "target_score",
        "species",
        "mirna_name",
    ]
    assert [str(dtype) for dtype in got.dtypes] == [
        TEXT_DTYPE,
        TEXT_DTYPE,
        "float64",
        TEXT_DTYPE,
        TEXT_DTYPE,
    ]
    assert isinstance(got.index, pd.RangeIndex)
    assert_frame_matches_golden(got.rename(columns=R_NAMES), "analyses/mirdb.tsv")


def test_read_mirdb_reads_gzip(tmp_path):
    path = tmp_path / "miRDB_v6.0_prediction_result.txt.gz"
    path.write_bytes(gzip.compress(SAMPLE.read_bytes()))
    pd.testing.assert_frame_equal(read_mirdb(path), read_mirdb(SAMPLE))


def test_read_mirdb_accepts_a_str_path():
    pd.testing.assert_frame_equal(read_mirdb(os.fspath(SAMPLE)), read_mirdb(SAMPLE))


def test_read_mirdb_keeps_only_the_species():
    got = read_mirdb(SAMPLE, species="rno")
    assert got["mirdb_id"].tolist() == ["rno-miR-192-5p", "rno-miR-192-5p", "rno-let-7a-2-3p"]
    assert got["target_refseq_id"].tolist() == ["NM_000301", "XM_000302", "NM_000303"]
    full = read_mirdb(SAMPLE)
    pd.testing.assert_frame_equal(got, full[full["species"] == "rno"].reset_index(drop=True))


@pytest.mark.parametrize("species", ["rat", "RNO", "rno-miR-192-5p"])
def test_read_mirdb_rejects_an_unknown_species(species):
    codes = "cfa, gga, hsa, lin4, mmu, rno"
    message = f"no row has the species {species!r}; expected one of {codes}"
    with pytest.raises(InputError, match=failure(SAMPLE, message)):
        read_mirdb(SAMPLE, species=species)


def test_read_mirdb_rejects_a_species_that_is_not_text():
    message = "species 5 has type int; expected a str or None"
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        read_mirdb(SAMPLE, species=5)


@pytest.mark.parametrize(
    "score",
    [
        "abc",
        "",
        "NA",
        "NaN",
        "Inf",
        "-inf",
        "1e999",
        "1_0",
        "6 0",
        "   ",
        "\u00a05",
        "\u0665",
        "\uff15\uff10",
    ],
)
def test_read_mirdb_rejects_a_score_that_is_not_a_number_D049(tmp_path, score):
    path = write_lines(tmp_path, ["hsa-a\tNM_1\t50", f"hsa-b\tNM_2\t{score}"])
    message = f"line 2: score {score!r} is not a finite number; expected a decimal number"
    with pytest.raises(InputError, match=failure(path, message)):
        read_mirdb(path)


def test_read_mirdb_checks_the_score_of_the_rows_of_the_other_species(tmp_path):
    path = write_lines(tmp_path, ["hsa-a\tNM_1\tabc", "rno-b\tNM_2\t50"])
    message = "line 1: score 'abc' is not a finite number; expected a decimal number"
    with pytest.raises(InputError, match=failure(path, message)):
        read_mirdb(path, species="rno")


@pytest.mark.parametrize(
    ("text", "value"), [(" 60", 60), ("60 ", 60), (" 60 ", 60), ("  1e2  ", 100)]
)
def test_read_mirdb_reads_a_score_with_spaces_around_it(tmp_path, text, value):
    """R 4.5.3 reads each of these as a number. D-049 gives the probe."""
    path = write_lines(tmp_path, ["hsa-a\tNM_1\t50", f"hsa-b\tNM_2\t{text}"])
    assert read_mirdb(path)["target_score"].tolist() == [50, value]


def test_read_mirdb_reads_each_form_of_a_decimal_number(tmp_path):
    forms = ["50", "96.25", "+5", "-5", ".5", "5.", "1e2", "1.5E-1"]
    path = write_lines(tmp_path, [f"hsa-a\tNM_1\t{form}" for form in forms])
    assert read_mirdb(path)["target_score"].tolist() == [50, 96.25, 5, -5, 0.5, 5, 100, 0.15]


def test_read_mirdb_keeps_a_hash_a_quote_and_the_text_na_D049(tmp_path):
    lines = [
        "hsa-a#1\tNM_1\t50",
        "hsa-b'2\tNM_2\t60",
        'hsa-c"3\tNM_3\t70',
        '"hsa-d\tNM_4\t80',
        "'hsa-e\tNM_5\t90",
        "NA\tNM#_6\t95",
    ]
    got = read_mirdb(write_lines(tmp_path, lines))
    assert got["mirdb_id"].tolist() == ["hsa-a#1", "hsa-b'2", 'hsa-c"3', '"hsa-d', "'hsa-e", "NA"]
    assert got["target_refseq_id"].tolist() == ["NM_1", "NM_2", "NM_3", "NM_4", "NM_5", "NM#_6"]
    assert got["species"].tolist() == ["hsa", "hsa", "hsa", '"hsa', "'hsa", "NA"]
    assert got["mirna_name"].tolist() == ["a#1", "b'2", 'c"3', "d", "e", ""]
    assert not got.isna().to_numpy().any()


@pytest.mark.parametrize(
    ("lines", "ids"),
    [
        (
            # The R input of the last note of mirdb.tsv: R drops 2 of these 4 lines.
            ["hsa-a\tNM_1\t50", "hsa-b\tNM'_2\t10", "hsa-c\tNM_3\t11", "hsa-d\tNM_4\t12"],
            ["hsa-a", "hsa-b", "hsa-c", "hsa-d"],
        ),
        (
            # A field that starts with a quote character, and a later line that closes it. A reader
            # that takes a quote character joins the lines 2 and 3 into one row.
            ["hsa-a\tNM_1\t50", '"hsa-b\tNM_2\t10', 'hsa-c"\tNM_3\t11', "hsa-d\tNM_4\t12"],
            ["hsa-a", '"hsa-b', 'hsa-c"', "hsa-d"],
        ),
    ],
    ids=["note", "closing-quote"],
)
def test_read_mirdb_keeps_each_line_that_a_quote_character_would_join_D049(tmp_path, lines, ids):
    assert golden_note("analyses/mirdb.tsv", r"line 2 of (\d+) lines") == str(len(lines))
    assert golden_note("analyses/mirdb.tsv", r"gives (\d+) rows with the IDs") == "2"
    got = read_mirdb(write_lines(tmp_path, lines))
    assert got["mirdb_id"].tolist() == ids
    assert got["target_score"].tolist() == [50, 10, 11, 12]


def test_read_mirdb_gives_the_full_id_as_the_species_of_an_id_without_a_hyphen(tmp_path):
    got = read_mirdb(write_lines(tmp_path, ["lin4\tNM_1\t55"]))
    assert got["species"].tolist() == ["lin4"]
    assert got["mirna_name"].tolist() == [""]


@pytest.mark.parametrize("blanks", [0, 1, 3])
def test_read_mirdb_rejects_a_line_with_2_fields(tmp_path, blanks):
    path = write_lines(tmp_path, ["hsa-a\tNM_1\t50", *[""] * blanks, "hsa-b\tNM_2"])
    message = f"line {2 + blanks} has 2 fields; expected 3: {LINE_FORMAT}"
    with pytest.raises(InputError, match=failure(path, message)):
        read_mirdb(path)


@pytest.mark.parametrize("blanks", [0, 1, 3])
def test_read_mirdb_rejects_a_line_with_4_fields(tmp_path, blanks):
    path = write_lines(tmp_path, ["hsa-a\tNM_1\t50", *[""] * blanks, "hsa-b\tNM_2\t60\textra"])
    message = f"line {2 + blanks} has 4 fields; expected 3: {LINE_FORMAT}"
    with pytest.raises(InputError, match=failure(path, message)):
        read_mirdb(path)


@pytest.mark.parametrize("line", ["hsa-b", "   "])
def test_read_mirdb_rejects_a_line_with_1_field(tmp_path, line):
    path = write_lines(tmp_path, ["hsa-a\tNM_1\t50", line])
    message = f"line 2 has 1 field; expected 3: {LINE_FORMAT}"
    with pytest.raises(InputError, match=failure(path, message)):
        read_mirdb(path)


def test_read_mirdb_rejects_a_url():
    message = "'https://x': the path is a URL; expected a local path"
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        read_mirdb("https://x")


@pytest.mark.parametrize("text", ["", "\n", "\r\n\n\n"])
def test_read_mirdb_rejects_an_empty_file(tmp_path, text):
    path = tmp_path / "x.txt"
    path.write_text(text)
    with pytest.raises(
        InputError, match=failure(path, "the file has no line; expected miRDB lines")
    ):
        read_mirdb(path)


@pytest.mark.parametrize(
    ("data", "number", "reason"),
    [
        (b"hsa-a\tNM_1\t50\nhsa-\xff\tNM_2\t60\n", 2, "invalid start byte"),
        (b"hsa-a\tNM_1\t50\nhsa-b\tNM_2\t6\xc3", 2, "unexpected end of data"),
        (b"\xef\xbbhsa-a\tNM_1\t50\n", 1, "invalid continuation byte"),
        (gzip.compress(b"hsa-a\tNM_1\t50\n\nhsa-\xff\tNM_2\t60\n"), 3, "invalid start byte"),
    ],
    ids=["a bad byte", "the end of the file", "an incomplete byte order mark", "gzip"],
)
def test_read_mirdb_names_the_line_of_text_that_is_not_utf8(tmp_path, data, number, reason):
    path = tmp_path / "x.txt"
    path.write_bytes(data)
    message = f"line {number}: {NOT_UTF8.format(reason)}"
    with pytest.raises(InputError, match=failure(path, message)):
        read_mirdb(path)


@pytest.mark.parametrize("chunk", [None, 1], ids=["default", "1"])
def test_read_mirdb_names_a_bad_line_before_text_that_is_not_utf8_after_many_lines(
    tmp_path, monkeypatch, chunk
):
    """The text reader of Python decodes about 8 KiB ahead of the line that it gives."""
    if chunk is not None:
        monkeypatch.setattr("seedmatch.mirdb._CHUNK_LINES", chunk)
    lines = ["hsa-a\tNM_1\t50", "hsa-b\tNM_2", *["hsa-c\tNM_3\t60"] * 1000, "hsa-\udcff\tNM_4\t70"]
    path = write_lines(tmp_path, lines)
    with pytest.raises(
        InputError, match=failure(path, f"line 2 has 2 fields; expected 3: {LINE_FORMAT}")
    ):
        read_mirdb(path)


def test_read_mirdb_rejects_a_line_above_the_limit_D049(tmp_path, monkeypatch):
    monkeypatch.setattr("seedmatch._lines.LINE_LIMIT", 20)
    path = write_lines(tmp_path, ["hsa-a\tNM_1\t50", "hsa-b\tNM_2\t" + "9" * 20])
    message = "line 2: the line has more than 20 characters; expected at most 20, with the line end"
    with pytest.raises(InputError, match=failure(path, message)):
        read_mirdb(path)


def test_read_mirdb_reads_a_line_of_the_limit_length(tmp_path, monkeypatch):
    line = "hsa-a\tNM_1\t50"
    monkeypatch.setattr("seedmatch._lines.LINE_LIMIT", len(line) + 1)
    assert len(read_mirdb(write_lines(tmp_path, [line, line]))) == 2


@pytest.mark.parametrize("head", [b"", b"\n", b"\r\n\n"], ids=["data", "blank line", "blank lines"])
def test_read_mirdb_skips_a_byte_order_mark_at_the_start_of_the_file(tmp_path, head):
    """R reads these files: R skips the mark, and the line after it is blank or a data line."""
    path = tmp_path / "x.txt"
    path.write_bytes(b"\xef\xbb\xbf" + head + SAMPLE.read_bytes())
    pd.testing.assert_frame_equal(read_mirdb(path), read_mirdb(SAMPLE))


@pytest.mark.parametrize(
    "line", ["hsa-\x00b\tNM_2\t60", "hsa-b\tNM\x00_2\t60", "hsa-b\tNM_2\t6\x000"]
)
def test_read_mirdb_rejects_a_nul_character_D049(tmp_path, line):
    path = write_lines(tmp_path, ["hsa-a\tNM_1\t50", line])
    message = "line 2: the line has a NUL character (U+0000); expected text without NUL characters"
    with pytest.raises(InputError, match=failure(path, message)):
        read_mirdb(path)


@pytest.mark.parametrize("ending", ["\r\n", "\r"])
def test_read_mirdb_reads_each_line_end(tmp_path, ending):
    path = tmp_path / "x.txt"
    path.write_bytes(SAMPLE.read_bytes().replace(b"\n", ending.encode()))
    pd.testing.assert_frame_equal(read_mirdb(path), read_mirdb(SAMPLE))


ENDINGS = pytest.mark.parametrize("ending", ["\n", "\r\n", "\r"], ids=["LF", "CRLF", "CR"])
SPACE_LINES = ["hsa-a\tNM_1\t50", " hsa-b\tNM_2\t60", " hsa-c\tNM_3\t70"]


@ENDINGS
def test_read_mirdb_reads_each_line_end_before_lines_that_start_with_a_space(tmp_path, ending):
    got = read_mirdb(write_lines(tmp_path, SPACE_LINES, ending=ending))
    assert got["mirdb_id"].tolist() == ["hsa-a", " hsa-b", " hsa-c"]


@ENDINGS
def test_read_mirdb_names_a_bad_line_after_lines_that_start_with_a_space(tmp_path, ending):
    path = write_lines(tmp_path, [*SPACE_LINES, "hsa-d\tNM_4\tabc"], ending=ending)
    message = "line 4: score 'abc' is not a finite number; expected a decimal number"
    with pytest.raises(InputError, match=failure(path, message)):
        read_mirdb(path)


@pytest.mark.parametrize(
    ("text", "chunk", "ids"),
    [
        ("hsa-a\tNM_1\t50\n\ufeffhsa-b\tNM_2\t60\n", 1, ["hsa-a", "\ufeffhsa-b"]),
        ("\ufeff\ufeffhsa-a\tNM_1\t50\nhsa-b\tNM_2\t60\n", None, ["\ufeffhsa-a", "hsa-b"]),
        ("\n\ufeffhsa-a\tNM_1\t50\nhsa-b\tNM_2\t60\n", None, ["\ufeffhsa-a", "hsa-b"]),
    ],
    ids=["first line of a chunk", "after the byte order mark", "after a blank line"],
)
def test_read_mirdb_keeps_a_u_feff_after_the_start_of_the_file(
    tmp_path, monkeypatch, text, chunk, ids
):
    if chunk is not None:
        monkeypatch.setattr("seedmatch.mirdb._CHUNK_LINES", chunk)
    path = tmp_path / "x.txt"
    path.write_text(text, encoding="utf-8")
    assert read_mirdb(path)["mirdb_id"].tolist() == ids


def test_read_mirdb_reads_a_last_line_without_a_line_end(tmp_path):
    path = tmp_path / "x.txt"
    path.write_bytes(SAMPLE.read_bytes().rstrip(b"\n"))
    pd.testing.assert_frame_equal(read_mirdb(path), read_mirdb(SAMPLE))


@pytest.mark.parametrize("chunk", [1, 2, 5, 11, 24, 25])
def test_read_mirdb_gives_the_same_table_for_each_chunk_size(tmp_path, monkeypatch, chunk):
    path = tmp_path / "spaced.txt"
    path.write_text("".join(f"{line}\n\n" for line in SAMPLE.read_text().splitlines()))
    expected = read_mirdb(SAMPLE)
    monkeypatch.setattr("seedmatch.mirdb._CHUNK_LINES", chunk)
    pd.testing.assert_frame_equal(read_mirdb(path), expected)
    rat = expected[expected["species"] == "rno"].reset_index(drop=True)
    pd.testing.assert_frame_equal(read_mirdb(path, species="rno"), rat)


@pytest.mark.parametrize("chunk", [1, 2, 3, 4, 5, 6, 7])
def test_read_mirdb_numbers_the_lines_across_chunks(tmp_path, monkeypatch, chunk):
    lines = ["hsa-a\tNM_1\t50", "", "hsa-b\tNM_2\t60", "", "", "hsa-c\tNM_3\tabc"]
    path = write_lines(tmp_path, lines)
    monkeypatch.setattr("seedmatch.mirdb._CHUNK_LINES", chunk)
    message = "line 6: score 'abc' is not a finite number; expected a decimal number"
    with pytest.raises(InputError, match=failure(path, message)):
        read_mirdb(path)


@pytest.mark.parametrize("chunk", [1, 2, 3, 4, 5, 6, 7])
@pytest.mark.parametrize(
    ("line", "message"),
    [
        ("hsa-c\tNM_3", f"line 6 has 2 fields; expected 3: {LINE_FORMAT}"),
        ("hsa-c\tNM_3\t70\textra", f"line 6 has 4 fields; expected 3: {LINE_FORMAT}"),
        (
            "hsa-c\tNM_3\t" + "9" * 20,
            "line 6: the line has more than 25 characters; expected at most 25, with the line end",
        ),
        (
            "hsa-\x00c\tNM_3\t70",
            "line 6: the line has a NUL character (U+0000); expected text without NUL characters",
        ),
    ],
    ids=["2 fields", "4 fields", "above the limit", "NUL"],
)
def test_read_mirdb_numbers_a_bad_line_across_chunks(tmp_path, monkeypatch, chunk, line, message):
    lines = ["hsa-a\tNM_1\t50", "", "hsa-b\tNM_2\t60", "", "", line]
    path = write_lines(tmp_path, lines)
    monkeypatch.setattr("seedmatch.mirdb._CHUNK_LINES", chunk)
    monkeypatch.setattr("seedmatch._lines.LINE_LIMIT", 25)
    with pytest.raises(InputError, match=failure(path, message)):
        read_mirdb(path)


SCORE_ABC = "line 3: score 'abc' is not a finite number; expected a decimal number"
NUL_3 = "line 3: the line has a NUL character (U+0000); expected text without NUL characters"
FIELDS_3 = "line 3 has 2 fields; expected 3: the miRDB ID, the RefSeq ID, and the score"
LONG_3 = "line 3: the line has more than 25 characters; expected at most 25, with the line end"
UTF8_3 = f"line 3: {NOT_UTF8.format('invalid start byte')}"


@pytest.mark.parametrize("chunk", [None, 1, 2], ids=["default", "1", "2"])
@pytest.mark.parametrize(
    ("bad_lines", "message"),
    [
        (["hsa-c\tNM_3\tabc", "hsa-d\tNM_4"], SCORE_ABC),
        (["hsa-c\tNM_3\tabc", "hsa-d\tNM_4\t4\x000"], SCORE_ABC),
        (["hsa-c\tNM_3\tabc", "hsa-d\tNM_4\t" + "4" * 30], SCORE_ABC),
        (["hsa-c\tNM_3\t3\x000", "hsa-d\tNM_4"], NUL_3),
        (["hsa-c\tNM_3\t3\x000", "hsa-d\tNM_4\t" + "4" * 30], NUL_3),
        (["hsa-c\tNM_3", "hsa-d\tNM_4\tabc"], FIELDS_3),
        (["hsa-c\tNM_3\t3\x000", "hsa-d\tNM_4\tabc"], NUL_3),
        (["hsa-c\tNM_3\t" + "3" * 30, "hsa-d\tNM_4\tabc"], LONG_3),
        (["hsa-c\tNM_3\tabc", "hsa-\udcffd\tNM_4\t40"], SCORE_ABC),
        (["hsa-c\tNM_3", "hsa-\udcffd\tNM_4\t40"], FIELDS_3),
        (["hsa-\udcffc\tNM_3\t30", "hsa-d\tNM_4\tabc"], UTF8_3),
        (["hsa-\udcffc\tNM_3", "hsa-d\tNM_4\tabc"], UTF8_3),
    ],
    ids=[
        "value, fields",
        "value, NUL",
        "value, long",
        "NUL, fields",
        "NUL, long",
        "fields, value",
        "NUL, value",
        "long, value",
        "value, UTF-8",
        "fields, UTF-8",
        "UTF-8, value",
        "UTF-8 and fields, value",
    ],
)
def test_read_mirdb_names_the_first_bad_line_for_each_chunk_size(
    tmp_path, monkeypatch, chunk, bad_lines, message
):
    if chunk is not None:
        monkeypatch.setattr("seedmatch.mirdb._CHUNK_LINES", chunk)
    monkeypatch.setattr("seedmatch._lines.LINE_LIMIT", 25)
    path = write_lines(tmp_path, ["hsa-a\tNM_1\t50", "hsa-b\tNM_2\t60", *bad_lines])
    with pytest.raises(InputError, match=failure(path, message)):
        read_mirdb(path)


@OBJECT_TEXT_ONLY
@pytest.mark.parametrize("chunk", [None, 1])
def test_equal_text_values_share_one_object(tmp_path, monkeypatch, chunk):
    lines = [
        "hsa-let-7\tNM_1\t50",
        "hsa-let-7\tNM_2\t60",
        "hsa-miR-1\tNM_1\t70",
        "mmu-let-7\tNM_1\t80",
    ]
    if chunk is not None:
        monkeypatch.setattr("seedmatch.mirdb._CHUNK_LINES", chunk)
    got = read_mirdb(write_lines(tmp_path, lines))
    assert got["mirdb_id"][0] is got["mirdb_id"][1]
    assert got["mirna_name"][0] is got["mirna_name"][1]
    assert got["species"][0] is got["species"][2]
    assert got["mirna_name"][0] is got["mirna_name"][3]
    assert got["target_refseq_id"][0] is got["target_refseq_id"][2]
    assert got["target_refseq_id"][0] is got["target_refseq_id"][3]


def test_top_level_package_exports_read_mirdb():
    assert seedmatch.read_mirdb is read_mirdb
