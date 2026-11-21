import io
import re

import numpy as np
import pytest

from seedmatch._lines import (
    LINE_LIMIT,
    data_bytes,
    field_lines,
    frame_of,
    not_utf8,
    parse_data_lines,
    read_chunks,
)
from seedmatch.errors import InputError

DESCRIBED = "x, y, and z"
LONG_LINE = (
    f"'f': line 1: the line has more than {LINE_LIMIT} characters; "
    f"expected at most {LINE_LIMIT}, with the line end"
)


def failure(text):
    return f"^{re.escape(text)}$"


def data_flags(chunk, first=1):
    return field_lines(chunk, first, "f", separator=" ", fields=3, described=DESCRIBED)


def data_line(length):
    return "a b " + "1" * (length - 4)


class Pipe(io.StringIO):
    def seekable(self):
        return False

    def seek(self, *_args):
        message = "seek"
        raise io.UnsupportedOperation(message)


def test_field_lines_count_the_fields_of_the_separator():
    lines = data_flags(["a b 1\n", "\n", "c d 2\r\n"])
    assert lines.is_data.tolist() == [True, False, True]
    assert lines.error is None


def test_field_lines_accept_2_fields():
    got = field_lines(["a b\n", "\n"], 1, "f", separator=" ", fields=2, described=DESCRIBED)
    assert got.is_data.tolist() == [True, False]
    assert got.error is None


def test_field_lines_name_the_line_of_a_bad_count():
    lines = data_flags(["a b 1\n", "a b\n"], first=7)
    assert lines.is_data.tolist() == [True, False]
    assert type(lines.error) is InputError
    assert str(lines.error) == f"'f': line 8 has 2 fields; expected 3: {DESCRIBED}"


def test_field_lines_write_one_field_in_the_singular():
    lines = data_flags(["a b 1\n", "\n", "a\n"])
    assert lines.is_data.tolist() == [True, False, False]
    assert str(lines.error) == f"'f': line 3 has 1 field; expected 3: {DESCRIBED}"


def test_field_lines_mark_no_data_line_from_the_first_bad_line_on():
    lines = data_flags(["a b 1\n", "a b\n", "c d 2\n", "\n"])
    assert lines.is_data.tolist() == [True, False, False, False]
    assert str(lines.error) == f"'f': line 2 has 2 fields; expected 3: {DESCRIBED}"


def test_field_lines_reject_a_long_line():
    line = data_line(LINE_LIMIT + 1)
    assert len(line) == LINE_LIMIT + 1
    lines = data_flags(["a b 1\n", line])
    assert lines.is_data.tolist() == [True, False]
    assert str(lines.error) == LONG_LINE.replace("line 1", "line 2")


def test_field_lines_accept_a_line_of_the_limit_length():
    line = data_line(LINE_LIMIT)
    assert len(line) == LINE_LIMIT
    lines = data_flags([line])
    assert lines.is_data.tolist() == [True]
    assert lines.error is None


def test_field_lines_name_the_bad_count_of_a_line_of_the_limit_length():
    line = "a" * (LINE_LIMIT - 1) + "\n"
    assert len(line) == LINE_LIMIT
    lines = data_flags([line])
    assert lines.is_data.tolist() == [False]
    assert str(lines.error) == f"'f': line 1 has 1 field; expected 3: {DESCRIBED}"


@pytest.mark.parametrize(
    ("limit", "line", "message"),
    [
        (
            10,
            "a\0" + "b" * 10 + "\n",
            (
                "'f': line 1: the line has more than 10 characters; expected at most 10, "
                "with the line end"
            ),
        ),
        (LINE_LIMIT, "a\0 b\n", f"'f': line 1 has 2 fields; expected 3: {DESCRIBED}"),
    ],
    ids=["the length first", "the fields before the NUL"],
)
def test_field_lines_check_the_length_then_the_fields_then_the_nul_of_a_line(
    monkeypatch, limit, line, message
):
    monkeypatch.setattr("seedmatch._lines.LINE_LIMIT", limit)
    lines = data_flags([line])
    assert lines.is_data.tolist() == [False]
    assert str(lines.error) == message


@pytest.mark.parametrize("fields", [1, 0, -1])
def test_field_lines_reject_fewer_than_two_fields(fields):
    with pytest.raises(ValueError, match=failure(f"fields is {fields}; expected 2 or more")) as got:
        field_lines(["a\n", "\n"], 1, "f", separator=" ", fields=fields, described=DESCRIBED)
    # InputError is a ValueError too, but a bad call of the helper is not bad input.
    assert got.type is ValueError


def test_data_bytes_join_the_data_lines_in_utf8():
    chunk = ["a b 1\n", "\n", "é b 2\n"]
    is_data = data_flags(chunk).is_data
    assert data_bytes(chunk, is_data) == "\na b 1\né b 2\n".encode()


@pytest.mark.parametrize("ending", ["\n", "\r\n", "\r"])
def test_data_bytes_turn_each_line_end_into_lf_after_a_blank_line(ending):
    chunk = [f"a b 1{ending}", ending, f"\ufeffc d 2{ending}", "e f 3"]
    is_data = data_flags(chunk).is_data
    assert data_bytes(chunk, is_data) == "\na b 1\n\ufeffc d 2\ne f 3".encode()


def test_parse_data_lines_read_a_mix_of_line_ends():
    chunk = ["a b 1\r\n", "c d 2\r", "e f 3\n", "\r", "g h 4\r\n", "\ufeffi j 5\r", "k l 6"]
    is_data = data_flags(chunk).is_data
    expected = "\na b 1\nc d 2\ne f 3\ng h 4\n\ufeffi j 5\nk l 6"
    assert data_bytes(chunk, is_data) == expected.encode()
    frame = parse_data_lines(chunk, is_data, separator=" ", names=["x", "y", "z"])
    assert frame["x"].tolist() == ["a", "c", "e", "g", "\ufeffi", "k"]
    assert frame["z"].tolist() == ["1", "2", "3", "4", "5", "6"]


@pytest.mark.parametrize(
    ("chunk", "flags", "number"),
    [
        (["a b 1\n", "a b 2\n", "a\0 b 3\n"], [True, True, False], 12),
        (["a b 1\n", "\n", "a b 2\n", "a\0 b 3\n"], [True, False, True, False], 13),
    ],
    ids=["third line", "after a blank line"],
)
def test_field_lines_name_the_line_of_a_nul_character(chunk, flags, number):
    lines = data_flags(chunk, first=10)
    message = (
        f"'f': line {number}: the line has a NUL character (U+0000); "
        "expected text without NUL characters"
    )
    assert lines.is_data.tolist() == flags
    assert str(lines.error) == message


def test_not_utf8_name_the_source_the_reason_and_the_kind():
    error = UnicodeDecodeError("utf-8", b"\xff", 0, 1, "invalid start byte")
    got = not_utf8("f", error, "TRRUST")
    assert type(got) is InputError
    assert str(got) == (
        "'f': the file is not valid UTF-8 text (invalid start byte); "
        "expected a text or gzip TRRUST file"
    )


def test_frame_of_keeps_the_memory_and_the_dtype_of_each_array():
    first_texts = np.array(["a", "b"], dtype=object)
    second_texts = np.array(["c", "d"], dtype=object)
    scores = np.array([1.5, 2.5])
    frame = frame_of({"first": first_texts, "second": second_texts, "score": scores})
    assert list(frame.columns) == ["first", "second", "score"]
    assert frame["first"].dtype == np.dtype(object)
    assert frame["score"].dtype == np.dtype(np.float64)
    assert np.shares_memory(frame["first"].to_numpy(), first_texts)
    assert np.shares_memory(frame["second"].to_numpy(), second_texts)
    assert np.shares_memory(frame["score"].to_numpy(), scores)


def test_read_chunks_give_the_first_line_number_of_each_chunk():
    handle = io.StringIO("".join(f"l{i}\n" for i in range(5)))
    got = list(read_chunks(handle, 2))
    assert [first for first, _ in got] == [1, 3, 5]
    assert [chunk for _, chunk in got] == [["l0\n", "l1\n"], ["l2\n", "l3\n"], ["l4\n"]]


@pytest.mark.parametrize("text", ["", "\ufeff"], ids=["empty", "mark only"])
def test_read_chunks_give_no_chunk_for_a_file_without_a_line(text):
    assert list(read_chunks(io.StringIO(text), 10)) == []


@pytest.mark.parametrize("chunk_lines", [0, -1])
def test_read_chunks_reject_fewer_than_one_line_in_a_chunk(chunk_lines):
    message = f"chunk_lines is {chunk_lines}; expected 1 or more"
    with pytest.raises(ValueError, match=failure(message)) as got:
        list(read_chunks(io.StringIO("a\n"), chunk_lines))
    assert got.type is ValueError


@pytest.mark.parametrize("chunk_lines", [1, 2, 3])
def test_read_chunks_remove_the_byte_order_mark_of_line_1(chunk_lines):
    handle = io.StringIO("\ufeffa\n\ufeffb\n\ufeffc\n")
    lines = [line for _, chunk in read_chunks(handle, chunk_lines) for line in chunk]
    assert lines == ["a\n", "\ufeffb\n", "\ufeffc\n"]


@pytest.mark.parametrize("length", [LINE_LIMIT, LINE_LIMIT + 5])
def test_read_chunks_after_a_byte_order_mark_leave_a_long_first_line_long(length):
    handle = io.StringIO("\ufeff" + data_line(length) + "\nc d 2\n")
    ((_, chunk),) = read_chunks(handle, 10)
    assert str(data_flags(chunk).error) == LONG_LINE


@pytest.mark.parametrize("end", ["\n", "\r"])
def test_read_chunks_after_a_byte_order_mark_keep_a_first_line_of_the_limit_length(end):
    line = data_line(LINE_LIMIT - 1) + end
    ((_, chunk),) = read_chunks(io.StringIO("\ufeff" + line + "c d 2\n", newline=""), 10)
    assert chunk == [line, "c d 2\n"]
    assert data_flags(chunk).is_data.tolist() == [True, True]


def test_read_chunks_after_a_byte_order_mark_keep_a_last_line_of_the_limit_length():
    line = data_line(LINE_LIMIT)
    ((_, chunk),) = read_chunks(io.StringIO("\ufeff" + line), 10)
    assert chunk == [line]
    assert data_flags(chunk).is_data.tolist() == [True]


def test_read_chunks_read_a_handle_that_cannot_seek_with_a_byte_order_mark_and_a_long_line():
    handle = Pipe("\ufeff" + data_line(LINE_LIMIT) + "\nc d 2\n")
    ((_, chunk),) = read_chunks(handle, 10)
    assert str(data_flags(chunk).error) == LONG_LINE


def test_read_chunks_read_a_handle_that_cannot_seek_with_a_byte_order_mark_and_a_valid_line():
    line = data_line(LINE_LIMIT - 1) + "\n"
    ((_, chunk),) = read_chunks(Pipe("\ufeff" + line + "c d 2\n"), 10)
    assert chunk == [line, "c d 2\n"]
    assert data_flags(chunk).is_data.tolist() == [True, True]


def outcome(text):
    ((_, chunk),) = read_chunks(io.StringIO(text, newline=""), 10)
    lines = data_flags(chunk)
    if lines.error is not None:
        return "error", str(lines.error)
    return "lines", lines.is_data.tolist(), chunk


@pytest.mark.parametrize("end", ["\n", "\r\n", "\r"])
@pytest.mark.parametrize("length", [LINE_LIMIT - 2, LINE_LIMIT - 1, LINE_LIMIT, LINE_LIMIT + 5])
def test_read_chunks_give_the_same_result_with_and_without_a_byte_order_mark(length, end):
    text = data_line(length) + end + "c d 2\n"
    assert outcome("\ufeff" + text) == outcome(text)


def test_read_chunks_cut_a_long_line_into_parts_of_the_limit_plus_one():
    handle = io.StringIO("a\n" + "x" * (2 * (LINE_LIMIT + 1) + 1) + "\nend\n")
    (_, chunk), *rest = read_chunks(handle, 10)
    assert not rest
    assert [len(line) for line in chunk] == [2, LINE_LIMIT + 1, LINE_LIMIT + 1, 2, 4]
    assert [line.endswith("\n") for line in chunk] == [True, False, False, True, True]


@pytest.mark.parametrize("mark", ["", "\ufeff"])
def test_read_chunks_cut_a_long_first_line_at_the_limit_plus_two_with_the_mark_included(mark):
    handle = io.StringIO(mark + "x" * (2 * LINE_LIMIT) + "\n")
    (_, chunk), *rest = read_chunks(handle, 10)
    assert not rest
    assert len(chunk[0]) == LINE_LIMIT + 2 - len(mark)
    assert not chunk[0].endswith("\n")


def tab_lines(chunk):
    return field_lines(chunk, 1, "f", separator="\t", fields=3, described=DESCRIBED)


def test_parse_data_lines_keep_hash_quote_and_na_as_text():
    chunk = ["a#1\t'b\t\"c\n", "\n", 'NA\t"d\tNA\n', "e\t\t#\r\n", "'f\t\"g\th'\n"]
    frame = parse_data_lines(chunk, tab_lines(chunk).is_data, separator="\t", names=["x", "y", "z"])
    assert frame.to_dict("list") == {
        "x": ["a#1", "NA", "e", "'f"],
        "y": ["'b", '"d', "", '"g'],
        "z": ['"c', "NA", "#", "h'"],
    }
    assert [str(dtype) for dtype in frame.dtypes] == ["object"] * 3


def test_parse_data_lines_give_the_dtype_to_the_parser():
    chunk = ["a\t5\t1.5\n", "b\t6\t2\n"]
    frame = parse_data_lines(
        chunk,
        tab_lines(chunk).is_data,
        separator="\t",
        names=["x", "y", "z"],
        dtype={"x": str, "y": np.int64, "z": np.float64},
    )
    assert [str(dtype) for dtype in frame.dtypes] == ["object", "int64", "float64"]
    assert frame["y"].tolist() == [5, 6]


@pytest.mark.parametrize(("separator", "space"), [("\t", "    "), (" ", "\t\t\t\t")])
@pytest.mark.parametrize("before", range(8))
def test_parse_data_lines_keep_the_leading_whitespace_of_a_line_at_the_end_of_a_parser_buffer(
    separator, space, before
):
    # The C parser of pandas 2.3.3 reads its input in buffers of 262,144 bytes. With the option
    # skip_blank_lines, it lost the whitespace of a line that starts 1 to 4 bytes before the end.
    start = 262_144 - before
    plain = separator.join(["a", "b", "1"]) + "\n"
    count = (start - 1) // len(plain) - 2
    last = start - 1 - count * len(plain)
    chunk = [
        *[plain] * count,
        separator.join(["c", "d", "9" * (last - 5)]) + "\n",
        separator.join([f"{space}e", "f", "2"]) + "\n",
        plain,
    ]
    is_data = field_lines(chunk, 1, "f", separator=separator, fields=3, described=DESCRIBED).is_data
    assert data_bytes(chunk, is_data).index(f"{space}e".encode()) == start
    frame = parse_data_lines(chunk, is_data, separator=separator, names=["x", "y", "z"])
    assert frame["x"].tolist() == ["a"] * count + ["c", f"{space}e", "a"]
