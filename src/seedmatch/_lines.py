"""Helpers that read a table file in chunks of lines, shared by the readers of the text tables."""

from collections.abc import Iterable, Iterator, Mapping
from functools import partial
from itertools import chain, compress, islice, repeat
from typing import Any, Final, TextIO

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


def field_lines(  # noqa: PLR0913
    chunk: list[str], first: int, source: str, *, separator: str, fields: int, described: str
) -> NDArray[np.bool_]:
    """Return which lines of the chunk are data lines. A bad line raises `InputError`.

    A data line has `fields - 1` separators. A blank line is no data line, thus `fields` is 2 or
    more, and a smaller value raises `ValueError`. `first` is the number of the first line of the
    chunk, and `described` names the fields in the message.
    """
    if fields < _MIN_FIELDS:
        message = f"fields is {fields}; expected {_MIN_FIELDS} or more"
        raise ValueError(message)
    count = len(chunk)
    separators = np.fromiter(map(str.count, chunk, repeat(separator)), dtype=np.int64, count=count)
    lengths = np.fromiter(map(len, chunk), dtype=np.int64, count=count)
    is_data: NDArray[np.bool_] = separators == fields - 1
    for offset in np.flatnonzero(~is_data | (lengths > LINE_LIMIT)):
        line, number = chunk[offset], first + int(offset)
        if len(line) > LINE_LIMIT:
            message = (
                f"{source!r}: line {number}: the line has more than {LINE_LIMIT} characters; "
                f"expected at most {LINE_LIMIT}, with the line end"
            )
            raise InputError(message)
        text = line.rstrip("\r\n")
        if text:
            found = text.count(separator) + 1
            noun = "field" if found == 1 else "fields"
            message = (
                f"{source!r}: line {number} has {found} {noun}; expected {fields}: {described}"
            )
            raise InputError(message)
    return is_data


def data_bytes(chunk: list[str], is_data: NDArray[np.bool_], first: int, source: str) -> bytes:
    """Return the data lines of the chunk as UTF-8 bytes. A NUL character raises `InputError`.

    `first` is the number of the first line of the chunk. The message names the line of the NUL
    character.
    """
    flags = is_data.tolist()
    text = "".join(compress(chunk, flags))
    if "\0" in text:
        bad = next(i for i, line in enumerate(compress(chunk, flags)) if "\0" in line)
        number = first + int(np.flatnonzero(is_data)[bad])
        message = (
            f"{source!r}: line {number}: the line has a NUL character (U+0000); "
            "expected text without NUL characters"
        )
        raise InputError(message)
    # Bytes, not a StringIO: pandas holds a text buffer in 4-byte characters, and the peak
    # memory of a chunk is then twice as large.
    return text.encode()


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
