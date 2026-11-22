import dataclasses
import gzip
import os
import warnings

import pandas as pd
import pytest
from helpers import (
    DATA,
    NOT_UTF8,
    OBJECT_TEXT_ONLY,
    TEXT_DTYPE,
    exactly,
    golden,
    golden_note,
    golden_rows,
)

import seedmatch
from seedmatch.errors import InputError
from seedmatch.network import StringLinks, read_string_links, stringdb

SAMPLE = DATA / "string_links_sample.txt"
CUTOFFS = "network/string_cutoffs.tsv"
LINKS = "network/string_links.tsv"
HEADER = "protein1 protein2 combined_score"
FIELDS = "protein1, protein2, and combined_score"
NOT_A_SCORE = "is not an integer from 0 to 1000; expected a STRING combined score"
PROTEIN_1 = "32630.ENSSYNP00000000001"
PROTEIN_6 = "32630.ENSSYNP00000000006"
PREFIX = [HEADER, "9.A 9.B 500"]


def write_lines(tmp_path, lines, name="x.txt", ending="\n"):
    path = tmp_path / name
    # surrogateescape writes a "\udcff" in a line as the byte 0xFF, which is not UTF-8.
    path.write_bytes(
        "".join(f"{line}{ending}" for line in lines).encode("utf-8", "surrogateescape")
    )
    return path


def failure(path, text):
    return exactly(f"{os.fspath(path)!r}: {text}")


def sample_lines():
    return SAMPLE.read_text().splitlines()


def assert_same_links(got, expected):
    pd.testing.assert_frame_equal(got.links, expected.links)
    assert got.cutoff == expected.cutoff


@pytest.mark.parametrize(
    "case", golden_rows(CUTOFFS), ids=[row.case for row in golden_rows(CUTOFFS)]
)
def test_read_string_links_match_r(case):
    got = read_string_links(SAMPLE, percentile=float(case.percentile))
    table = golden(LINKS)
    expected = table[table["case"] == case.case]
    assert got.cutoff == float(case.cutoff)
    assert len(got.links) == int(case.n_above) == len(expected)
    assert got.links["protein1"].tolist() == expected["protein1"].tolist()
    assert got.links["protein2"].tolist() == expected["protein2"].tolist()
    assert got.links["combined_score"].tolist() == expected["combined_score"].astype(int).tolist()


def test_read_string_links_give_the_table_of_the_documentation():
    got = read_string_links(SAMPLE)
    assert isinstance(got, StringLinks)
    assert got.cutoff == 897.4999999999999
    assert type(got.cutoff) is float
    assert list(got.links.columns) == ["protein1", "protein2", "combined_score"]
    assert [str(dtype) for dtype in got.links.dtypes] == [TEXT_DTYPE, TEXT_DTYPE, "int64"]
    assert isinstance(got.links.index, pd.RangeIndex)
    assert got.links.to_numpy().tolist() == [
        [PROTEIN_1, PROTEIN_6, 900],
        [PROTEIN_6, PROTEIN_1, 900],
    ]
    with pytest.raises(dataclasses.FrozenInstanceError):
        got.cutoff = 1.0


def test_read_string_links_apply_the_threshold_before_the_percentile():
    got = read_string_links(SAMPLE, score_threshold=0, percentile=0)
    assert got.cutoff == 150.0
    assert len(got.links) == 24
    default = read_string_links(SAMPLE, percentile=0)
    assert default.cutoff == 400.0
    assert len(default.links) == 20


def test_read_string_links_accept_the_scores_0_and_1000(tmp_path):
    path = write_lines(tmp_path, [HEADER, "9.A 9.B 0", "9.B 9.A 1000", "9.C 9.D 500"])
    got = read_string_links(path, score_threshold=0, percentile=0)
    assert got.cutoff == 0.0
    assert got.links["combined_score"].tolist() == [1000, 500]


@pytest.mark.parametrize(
    ("line", "message"),
    [
        ("9.C 9.D abc", f"line 3: score 'abc' {NOT_A_SCORE}"),
        ("9.C 9.D", f"line 3 has 2 fields; expected 3: {FIELDS}"),
        ("9.C 9.D 700 9.E", f"line 3 has 4 fields; expected 3: {FIELDS}"),
        ("9.C 9.D 450.5", f"line 3: score '450.5' {NOT_A_SCORE}"),
        ("9.C 9.D -5", f"line 3: score '-5' {NOT_A_SCORE}"),
        ("9.C 9.D 5000", f"line 3: score '5000' {NOT_A_SCORE}"),
        ("9.C 9.D ", f"line 3: score '' {NOT_A_SCORE}"),
        ("9.C 9.D inf", f"line 3: score 'inf' {NOT_A_SCORE}"),
        ("9.C 9.D nan", f"line 3: score 'nan' {NOT_A_SCORE}"),
        ("9.C 9.D 1e400", f"line 3: score '1e400' {NOT_A_SCORE}"),
        ("9.C 9.D 99999999999999999999", f"line 3: score '99999999999999999999' {NOT_A_SCORE}"),
        (" 9.D 700", "line 3: protein1 is empty; expected a STRING protein ID"),
        ("9.C  700", "line 3: protein2 is empty; expected a STRING protein ID"),
        ("  ", "line 3: protein1 is empty; expected a STRING protein ID"),
        (" 9.D abc", "line 3: protein1 is empty; expected a STRING protein ID"),
    ],
    ids=[
        "text score",
        "2 fields",
        "4 fields",
        "decimal score",
        "negative score",
        "score above 1000",
        "empty score",
        "inf",
        "nan",
        "exponent too large",
        "score above 2^64",
        "empty protein1",
        "empty protein2",
        "spaces only",
        "empty protein1 before bad score",
    ],
)
def test_read_string_links_reject_a_bad_line_D059(tmp_path, line, message):
    path = write_lines(tmp_path, [*PREFIX, line])
    with warnings.catch_warnings():
        # pandas warns about a cast of inf, and the reader must not show that warning.
        warnings.simplefilter("error")
        with pytest.raises(InputError, match=failure(path, message)):
            read_string_links(path)


def test_read_string_links_name_the_first_bad_line_D059(tmp_path):
    path = write_lines(tmp_path, [*PREFIX, "9.C 9.D 4e2", "9.E 9.F abc", " 9.D 700"])
    with pytest.raises(InputError, match=failure(path, f"line 4: score 'abc' {NOT_A_SCORE}")):
        read_string_links(path)


def test_read_string_links_keep_hash_and_quote_text_D059(tmp_path):
    path = write_lines(tmp_path, [*PREFIX, "9.C#x 9.D 700", "9.C'x 9.D 700", '9."C 9.D 600'])
    got = read_string_links(path, percentile=0)
    assert got.cutoff == 500.0
    assert got.links["protein1"].tolist() == ["9.C#x", "9.C'x", '9."C']
    assert got.links["protein2"].tolist() == ["9.D", "9.D", "9.D"]
    assert got.links["combined_score"].tolist() == [700, 700, 600]


@pytest.mark.parametrize("pandas_rejects_the_integer_read", [False, True])
def test_read_string_links_accept_a_score_with_an_integer_value(
    tmp_path, monkeypatch, pandas_rejects_the_integer_read
):
    if pandas_rejects_the_integer_read:
        # pandas 2.3 reads these scores as int64. Another version can raise ValueError for them.
        read_text = stringdb.parse_data_lines

        def parse_data_lines(*args, dtype=str, **kwargs):
            if dtype is not str:
                message = "invalid literal for int() with base 10: '400.0'"
                raise ValueError(message)
            return read_text(*args, dtype=dtype, **kwargs)

        monkeypatch.setattr(stringdb, "parse_data_lines", parse_data_lines)
    lines = [HEADER, "9.A 9.B 400.0", "9.C 9.D 4e2", "9.E 9.F +400", "9.G 9.H 300"]
    got = read_string_links(write_lines(tmp_path, lines), score_threshold=0, percentile=0)
    assert got.cutoff == 300.0
    assert got.links["combined_score"].tolist() == [400, 400, 400]
    assert str(got.links["combined_score"].dtype) == "int64"


@pytest.mark.parametrize(
    ("first", "shown"),
    [
        ("a b c", "'a b c'"),
        ("9.A 9.B 500", "'9.A 9.B 500'"),
        ("protein1 protein2", "'protein1 protein2'"),
        (f"{HEADER} ", f"'{HEADER} '"),
        ("", "''"),
    ],
)
def test_read_string_links_reject_a_bad_header(tmp_path, first, shown):
    path = write_lines(tmp_path, [first, "9.A 9.B 500", "9.B 9.A 500"])
    message = f"line 1 is {shown}; expected the header '{HEADER}'"
    with pytest.raises(InputError, match=failure(path, message)):
        read_string_links(path)


def test_read_string_links_reject_an_empty_file(tmp_path):
    path = tmp_path / "x.txt"
    path.write_bytes(b"")
    message = f"the file has no line; expected the header '{HEADER}'"
    with pytest.raises(InputError, match=failure(path, message)):
        read_string_links(path)


def test_read_string_links_reject_a_file_without_a_link(tmp_path):
    message = "the file has no line after the header; expected STRING links"
    for name, lines in (("header.txt", [HEADER]), ("blank.txt", [HEADER, "", ""])):
        path = write_lines(tmp_path, lines, name)
        with pytest.raises(InputError, match=failure(path, message)):
            read_string_links(path)
    low = write_lines(tmp_path, [HEADER, "9.A 9.B 399", "9.B 9.A 399"], "low.txt")
    message = "no link has a score of 400 or more; expected at least one link"
    with pytest.raises(InputError, match=failure(low, message)):
        read_string_links(low)
    message = "no link has a score of 1000 or more; expected at least one link"
    with pytest.raises(InputError, match=failure(SAMPLE, message)):
        read_string_links(SAMPLE, score_threshold=1000)


def test_read_string_links_reject_a_long_line_D059(tmp_path):
    path = write_lines(tmp_path, [*PREFIX, "9" * 1030 + " 9.D 700"])
    message = (
        "line 3: the line has more than 1024 characters; expected at most 1024, with the line end"
    )
    with pytest.raises(InputError, match=failure(path, message)):
        read_string_links(path)


def test_read_string_links_read_a_line_of_the_limit_length(tmp_path):
    line = "9" * 1015 + " 9.D 700"
    assert len(line) + 1 == 1024
    got = read_string_links(write_lines(tmp_path, [*PREFIX, line]), percentile=0)
    assert got.links["protein1"].tolist() == ["9" * 1015]


def test_read_string_links_reject_a_nul_character_D059(tmp_path):
    path = write_lines(tmp_path, [*PREFIX, "9.C 9.\x00D 700"])
    message = "line 3: the line has a NUL character (U+0000); expected text without NUL characters"
    with pytest.raises(InputError, match=failure(path, message)):
        read_string_links(path)


def test_read_string_links_parse_no_line_from_a_nul_line_on(tmp_path, monkeypatch):
    reads = []
    read_text = stringdb.parse_data_lines

    def parse_data_lines(chunk, is_data, **kwargs):
        reads.append(is_data.tolist())
        return read_text(chunk, is_data, **kwargs)

    monkeypatch.setattr(stringdb, "parse_data_lines", parse_data_lines)
    path = write_lines(tmp_path, [*PREFIX, "9.C 9.\x00D 700", "9.G 9.H 600"])
    with pytest.raises(InputError, match="line 3: the line has a NUL character"):
        read_string_links(path)
    assert reads == [[False, True, False, False]]


def test_read_string_links_reject_a_url():
    message = "'https://x': the path is a URL; expected a local path"
    with pytest.raises(InputError, match=exactly(message)):
        read_string_links("https://x")


def test_read_string_links_name_the_line_of_text_that_is_not_utf8(tmp_path):
    path = tmp_path / "x.txt"
    path.write_bytes(HEADER.encode() + b"\n9.A 9.B 500\n9.C \xff 600\n")
    message = f"line 3: {NOT_UTF8.format('invalid start byte')}"
    with pytest.raises(InputError, match=failure(path, message)):
        read_string_links(path)


def test_read_string_links_raise_without_a_link_above_the_cutoff():
    message = (
        "no link has a score above the cutoff 900.0, the percentile 1.0 of the 22 scores of "
        "400 or more; expected at least one link"
    )
    with pytest.raises(InputError, match=failure(SAMPLE, message)):
        read_string_links(SAMPLE, percentile=1)
    note = r"^# percentile_1: .* stops with: (.*?), because"
    assert golden_note(CUTOFFS, note) == "non-character argument"


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        (
            {"score_threshold": -1},
            "score_threshold -1 is out of range; expected an integer from 0 to 1000",
        ),
        (
            {"score_threshold": 1001},
            "score_threshold 1001 is out of range; expected an integer from 0 to 1000",
        ),
        ({"score_threshold": True}, "score_threshold True has type bool; expected an integer"),
        ({"score_threshold": 400.0}, "score_threshold 400.0 has type float; expected an integer"),
        (
            {"percentile": 1.2},
            "percentile 1.2 is out of range; expected a number from 0 to 1",
        ),
        (
            {"percentile": -0.1},
            "percentile -0.1 is out of range; expected a number from 0 to 1",
        ),
        (
            {"percentile": float("nan")},
            "percentile nan is not a number; expected an int or a float, not NaN",
        ),
        ({"percentile": True}, "percentile True has type bool; expected an int or a float"),
        ({"percentile": "0.95"}, "percentile '0.95' has type str; expected an int or a float"),
    ],
)
def test_read_string_links_reject_a_bad_argument(kwargs, message):
    with pytest.raises(InputError, match=exactly(message)):
        read_string_links(SAMPLE, **kwargs)


def test_read_string_links_read_gzip(tmp_path):
    path = tmp_path / "10116.protein.links.v12.5.txt.gz"
    path.write_bytes(gzip.compress(SAMPLE.read_bytes()))
    assert_same_links(read_string_links(path), read_string_links(SAMPLE))


@pytest.mark.parametrize("ending", ["\r\n", "\r"])
def test_read_string_links_read_each_line_end(tmp_path, ending):
    path = tmp_path / "x.txt"
    path.write_bytes(SAMPLE.read_bytes().replace(b"\n", ending.encode()))
    assert_same_links(read_string_links(path), read_string_links(SAMPLE))


def test_read_string_links_skip_a_byte_order_mark(tmp_path):
    path = tmp_path / "x.txt"
    path.write_bytes(b"\xef\xbb\xbf" + SAMPLE.read_bytes())
    assert_same_links(read_string_links(path), read_string_links(SAMPLE))


def test_read_string_links_skip_a_blank_line(tmp_path):
    lines = sample_lines()
    path = write_lines(tmp_path, [lines[0], "", "", *lines[1:], ""])
    assert_same_links(read_string_links(path), read_string_links(SAMPLE))
    bad = write_lines(tmp_path, [HEADER, "", "", "9.C 9.D abc"], "bad.txt")
    with pytest.raises(InputError, match=failure(bad, f"line 4: score 'abc' {NOT_A_SCORE}")):
        read_string_links(bad)


ENDINGS = pytest.mark.parametrize("ending", ["\n", "\r\n", "\r"], ids=["LF", "CRLF", "CR"])
TAB_LINES = [HEADER, "9.A 9.B 500", "\t9.C 9.D 600", "\t9.E 9.F 700"]


@ENDINGS
def test_read_string_links_read_each_line_end_before_lines_that_start_with_a_tab(tmp_path, ending):
    path = write_lines(tmp_path, TAB_LINES, ending=ending)
    got = read_string_links(path, score_threshold=0, percentile=0)
    assert got.links["protein1"].tolist() == ["\t9.C", "\t9.E"]


@ENDINGS
def test_read_string_links_name_a_bad_line_after_lines_that_start_with_a_tab(tmp_path, ending):
    path = write_lines(tmp_path, [*TAB_LINES, "9.G 9.H abc"], ending=ending)
    with pytest.raises(InputError, match=failure(path, f"line 5: score 'abc' {NOT_A_SCORE}")):
        read_string_links(path)


@pytest.mark.parametrize(
    ("text", "chunk", "proteins"),
    [
        (f"{HEADER}\n9.A 9.B 500\n\ufeff9.C 9.D 600\n", 2, ["\ufeff9.C"]),
        (f"{HEADER}\n\ufeff9.A 9.B 600\n9.C 9.D 500\n", None, ["\ufeff9.A"]),
    ],
    ids=["first line of a chunk", "first line of the first chunk"],
)
def test_read_string_links_keep_a_u_feff_after_the_start_of_the_file(
    tmp_path, monkeypatch, text, chunk, proteins
):
    if chunk is not None:
        monkeypatch.setattr("seedmatch.network.stringdb._CHUNK_LINES", chunk)
    path = tmp_path / "x.txt"
    path.write_text(text, encoding="utf-8")
    got = read_string_links(path, score_threshold=0, percentile=0)
    assert got.links["protein1"].tolist() == proteins


@pytest.mark.parametrize("chunk", [1, 5])
def test_read_string_links_read_more_than_one_chunk(tmp_path, monkeypatch, chunk):
    expected = {p: read_string_links(SAMPLE, percentile=p) for p in (0, 0.5, 0.95)}
    monkeypatch.setattr("seedmatch.network.stringdb._CHUNK_LINES", chunk)
    for percentile, result in expected.items():
        assert_same_links(read_string_links(SAMPLE, percentile=percentile), result)
    lines = sample_lines()
    lines[12] = "32630.ENSSYNP00000000008 32630.ENSSYNP00000000009 abc"
    path = write_lines(tmp_path, lines)
    with pytest.raises(InputError, match=failure(path, f"line 13: score 'abc' {NOT_A_SCORE}")):
        read_string_links(path)


SCORE_3 = f"line 3: score 'abc' {NOT_A_SCORE}"
NUL_3 = "line 3: the line has a NUL character (U+0000); expected text without NUL characters"
FIELDS_3 = "line 3 has 2 fields; expected 3: protein1, protein2, and combined_score"
LONG_3 = "line 3: the line has more than 40 characters; expected at most 40, with the line end"
UTF8_3 = f"line 3: {NOT_UTF8.format('invalid start byte')}"


@pytest.mark.parametrize("chunk", [None, 1, 2], ids=["default", "1", "2"])
@pytest.mark.parametrize(
    ("bad_lines", "message"),
    [
        (["9.C 9.D abc", "9.E 9.F"], SCORE_3),
        (["9.C 9.D abc", "9.E 9.F 7\x0000"], SCORE_3),
        (["9.C 9.D abc", "9.E 9.F " + "7" * 40], SCORE_3),
        (["9.C 9.D 6\x0000", "9.E 9.F"], NUL_3),
        (["9.C 9.D 6\x0000", "9.E 9.F " + "7" * 40], NUL_3),
        (["9.C 9.D", "9.E 9.F abc"], FIELDS_3),
        (["9.C 9.D 6\x0000", "9.E 9.F abc"], NUL_3),
        (["9.C 9.D " + "6" * 40, "9.E 9.F abc"], LONG_3),
        (["9.C 9.D abc", "9.\udcffE 9.F 700"], SCORE_3),
        (["9.\udcffC 9.D 600", "9.E 9.F abc"], UTF8_3),
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
        "UTF-8, value",
    ],
)
def test_read_string_links_name_the_first_bad_line_for_each_chunk_size(
    tmp_path, monkeypatch, chunk, bad_lines, message
):
    if chunk is not None:
        monkeypatch.setattr("seedmatch.network.stringdb._CHUNK_LINES", chunk)
    # The header has 33 characters with its line end.
    monkeypatch.setattr("seedmatch._lines.LINE_LIMIT", 40)
    path = write_lines(tmp_path, [*PREFIX, *bad_lines])
    with pytest.raises(InputError, match=failure(path, message)):
        read_string_links(path)


@pytest.mark.parametrize("chunk", [None, 1, 2], ids=["default", "1", "2"])
def test_read_string_links_name_a_bad_header_before_a_later_bad_line(tmp_path, monkeypatch, chunk):
    if chunk is not None:
        monkeypatch.setattr("seedmatch.network.stringdb._CHUNK_LINES", chunk)
    path = write_lines(tmp_path, ["protein1 protein2 score", "9.A 9.B 500", "9.C 9.D"])
    message = f"line 1 is 'protein1 protein2 score'; expected the header '{HEADER}'"
    with pytest.raises(InputError, match=failure(path, message)):
        read_string_links(path)


@OBJECT_TEXT_ONLY
def test_read_string_links_share_equal_texts():
    links = read_string_links(SAMPLE, percentile=0.5).links
    assert links["protein1"][0] is links["protein2"][2]
    assert links["protein1"][0] is links["protein1"][1]


def test_read_string_links_is_in_the_top_level_package():
    assert seedmatch.read_string_links is read_string_links
    assert seedmatch.StringLinks is StringLinks
    assert "read_string_links" in seedmatch.__all__
    assert "StringLinks" in seedmatch.__all__
