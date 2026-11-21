"""Helpers that read a table file in chunks of lines, shared by the readers of the text tables."""

import csv
import io
from collections.abc import Hashable, Iterable, Iterator, Mapping, Sequence
from functools import partial
from itertools import chain, compress, islice, repeat
from operator import contains
from typing import Any, Final, NamedTuple, TextIO

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from seedmatch.errors import InputError

LINE_LIMIT: Final = 2**10
_MIN_FIELDS: Final = 2


def read_chunks(handle: TextIO, chunk_lines: int) -> Iterator[tuple[int, list[str]]]:
    """Yield the number of the first line of each chunk, and its lines with their line ends.

    A line of more than `LINE_LIMIT` characters comes as parts. Each part but the last has
    `LINE_LIMIT + 1` characters and no line end, thus a chunk holds a bounded number of
    characters. Line 1 is the exception: its first part has `LINE_LIMIT + 2` characters, or
    `LINE_LIMIT + 1` after the byte order mark, which line 1 loses. The handle does not need to
    seek. A `chunk_lines` below 1 raises `ValueError` at the first step of the iteration. The
    iteration decodes the file, thus a `UnicodeDecodeError` comes from the loop, and the loop must
    be inside the `try` that calls `not_utf8`.
    """
    if chunk_lines < 1:
        message = f"chunk_lines is {chunk_lines}; expected 1 or more"
        raise ValueError(message)
    next_line = partial(handle.readline, LINE_LIMIT + 1)
    # R skips a byte order mark, also before a blank first line. The mark takes 1 character of the
    # read, thus line 1 reads 1 more. Without it, a line at the limit loses its line end to a
    # blank line 2, and line 1 joins line 2 in `data_bytes`.
    head = handle.readline(LINE_LIMIT + 2).removeprefix("\ufeff")
    parts = chain([head] if head else [], iter(next_line, ""))
    first = 1
    while chunk := list(islice(parts, chunk_lines)):
        yield first, chunk
        first += len(chunk)


class DataLines(NamedTuple):
    """The data lines of a chunk before its first bad line, and the error of that line.

    The error is `None` if the chunk has no bad line.
    """

    is_data: NDArray[np.bool_]
    error: InputError | None


def field_lines(  # noqa: PLR0913
    chunk: list[str], first: int, source: str, *, separator: str, fields: int, described: str
) -> DataLines:
    """Return the data lines of the chunk and the error of its first bad line.

    A bad line has more than `LINE_LIMIT` characters, another field count, or a NUL character.
    The caller parses the data lines first, and then raises the error. Thus a bad value on an
    earlier line raises first.

    A data line has `fields - 1` separators. A blank line is no data line, and no bad line, thus
    `fields` is 2 or more, and a smaller value raises `ValueError`. `first` is the number of the
    first line of the chunk, and `described` names the fields in the message.
    """
    if fields < _MIN_FIELDS:
        message = f"fields is {fields}; expected {_MIN_FIELDS} or more"
        raise ValueError(message)
    count = len(chunk)
    separators = np.fromiter(map(str.count, chunk, repeat(separator)), dtype=np.int64, count=count)
    lengths = np.fromiter(map(len, chunk), dtype=np.int64, count=count)
    has_nul = np.fromiter(map(contains, chunk, repeat("\0")), dtype=np.bool_, count=count)
    is_data: NDArray[np.bool_] = separators == fields - 1
    for offset in np.flatnonzero(~is_data | has_nul | (lengths > LINE_LIMIT)):
        line, number = chunk[offset], first + int(offset)
        if len(line) > LINE_LIMIT:
            message = (
                f"{source!r}: line {number}: the line has more than {LINE_LIMIT} characters; "
                f"expected at most {LINE_LIMIT}, with the line end"
            )
        elif line.rstrip("\r\n") and separators[offset] != fields - 1:
            found = int(separators[offset]) + 1
            noun = "field" if found == 1 else "fields"
            message = (
                f"{source!r}: line {number} has {found} {noun}; expected {fields}: {described}"
            )
        elif "\0" in line:
            message = (
                f"{source!r}: line {number}: the line has a NUL character (U+0000); "
                "expected text without NUL characters"
            )
        else:
            continue
        is_data[offset:] = False
        return DataLines(is_data, InputError(message))
    return DataLines(is_data, None)


def data_bytes(chunk: list[str], is_data: NDArray[np.bool_]) -> bytes:
    """Return the data lines of the chunk as UTF-8 bytes, after one blank line.

    Each CR becomes LF, thus a CRLF gives a blank line, which pandas skips. pandas 2.3.3 removes a
    U+FEFF at the start of its input, and it misreads a lone CR before a line that starts with a
    space.
    """
    # Bytes, not a StringIO: pandas holds a text buffer in 4-byte characters, and the peak
    # memory of a chunk is then twice as large.
    data = "".join(chain(("\n",), compress(chunk, is_data.tolist()))).encode()
    # `open_text` reads with newline="", thus each CR in a line is a line end.
    return data.replace(b"\r", b"\n")


def parse_data_lines(
    chunk: list[str],
    is_data: NDArray[np.bool_],
    *,
    separator: str,
    names: Sequence[str],
    dtype: str
    | type[str]
    | np.dtype[Any]
    | Mapping[Hashable, str | type[str] | np.dtype[Any]] = str,
) -> pd.DataFrame:
    """Return the data lines of the chunk as a table with the columns `names`.

    A quote character, `#`, and the text `NA` are text, and a field keeps its spaces. `dtype`
    goes to `pandas.read_csv` as it is. `field_lines` must have made `is_data` for the same
    `separator` and the number of `names`.
    """
    return pd.read_csv(
        io.BytesIO(data_bytes(chunk, is_data)),
        sep=separator,
        header=None,
        names=list(names),
        dtype=dtype,
        quoting=csv.QUOTE_NONE,
        na_filter=False,
    )


def not_utf8(source: str, error: UnicodeDecodeError, kind: str) -> InputError:
    """Return the `InputError` for a file that is not valid UTF-8. `kind` names the file."""
    message = (
        f"{source!r}: the file is not valid UTF-8 text ({error.reason}); "
        f"expected a text or gzip {kind} file"
    )
    return InputError(message)


def vocabulary_codes(column: pd.Series, vocabulary: dict[str, int]) -> NDArray[np.int32]:
    """Return the code of each text of the column. A new text takes the next code."""
    codes, texts = pd.factorize(column)
    lookup = np.fromiter(
        (vocabulary.setdefault(text, len(vocabulary)) for text in texts),
        dtype=np.int32,
        count=len(texts),
    )
    return lookup[codes]


def shared_texts(texts: Iterable[str]) -> NDArray[np.object_]:
    """Return the texts as an object array in which equal texts are one object."""
    objects: dict[str, str] = {}
    return np.array([objects.setdefault(text, text) for text in texts], dtype=object)


def frame_of(columns: Mapping[str, NDArray[Any]]) -> pd.DataFrame:
    """Return the table of the columns, without a copy of the arrays."""
    # A dict of object arrays makes pandas scan each array for dates, which allocates a copy of
    # the table. A Series with its dtype set skips the scan.
    return pd.DataFrame(
        {
            name: pd.Series(values, dtype=values.dtype, copy=False)
            for name, values in columns.items()
        },
        copy=False,
    )
