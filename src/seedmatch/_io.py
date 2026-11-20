"""Paths that the readers open, and the text of the files that they read."""

import bz2
import gzip
import io
import lzma
import os
import re
import zlib
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Final, TextIO

from seedmatch.errors import InputError

_URL: Final = re.compile(r"^[A-Za-z][A-Za-z0-9+.-]*://")
_GZIP_MAGIC: Final = b"\x1f\x8b"
_BZIP2_MAGIC: Final = b"BZh"
_XZ_MAGIC: Final = b"\xfd7zXZ\x00"
_DAMAGED_DATA: Final = (EOFError, OSError, lzma.LZMAError, zlib.error)
INT64_MAX: Final = 2**63 - 1
# A genome coordinate is at most 2^62 - 1, so a coordinate + 1, and the width or the sum of two
# coordinates, stay inside a signed 64-bit integer.
COORDINATE_MAX: Final = 2**62 - 1
_INT64_DIGITS: Final = len(str(INT64_MAX))
_SHOWN_DIGITS: Final = 40


def local_path(path: str | os.PathLike[str]) -> Path:
    """Return the path as a `Path`.

    A URL raises `InputError`, because the core modules do not use the network.
    """
    shown = os.fspath(path)
    if _URL.match(shown):
        message = f"{shown!r}: the path is a URL; expected a local path"
        raise InputError(message)
    return Path(shown)


def parse_int64(shown: str, number: int, label: str, digits: str, maximum: int = INT64_MAX) -> int:
    """Return the value of a field of decimal digits on line `number` of the file `shown`.

    A value above `maximum` raises `InputError`. The default is 2^63 - 1, the limit of a signed
    64-bit integer and of an SQLite INTEGER. `maximum` has the form 2^k - 1, at most 2^63 - 1.
    """
    significant = digits.lstrip("0") or "0"
    # Python does not convert a string of more than 4300 digits, so the length comes first.
    if len(significant) <= _INT64_DIGITS and (value := int(significant)) <= maximum:
        return value
    shown_digits = digits if len(digits) <= _SHOWN_DIGITS else f"{digits[:_SHOWN_DIGITS]}..."
    raise too_large_error(shown, number, label, shown_digits, maximum)


def too_large_error(
    shown: str, number: int, what: str, value: str, maximum: int = INT64_MAX
) -> InputError:
    """Return the error for a number on line `number` of the file `shown` above `maximum`."""
    message = (
        f"{shown!r}: line {number}: {what} {value} is too large; "
        f"expected at most {limit_text(maximum)}"
    )
    return InputError(message)


def limit_text(maximum: int) -> str:
    """Return the text of a limit of the form 2^k - 1, for example `2^63 - 1`."""
    return f"2^{maximum.bit_length()} - 1"


def open_binary(path: Path) -> io.BufferedReader:
    """Open the file for binary reads. A directory raises `InputError`."""
    try:
        return path.open("rb")
    except IsADirectoryError as exc:
        message = f"{os.fspath(path)!r}: the path is a directory; expected a file"
        raise InputError(message) from exc


def _decompressed(
    raw: io.BufferedReader,
) -> gzip.GzipFile | bz2.BZ2File | lzma.LZMAFile | io.BufferedReader:
    """Return a reader of the bytes after the compression that the first bytes show, if any."""
    head = raw.peek(len(_XZ_MAGIC))
    if head.startswith(_GZIP_MAGIC):
        return gzip.GzipFile(fileobj=raw)
    if head.startswith(_BZIP2_MAGIC):
        return bz2.BZ2File(raw)
    if head.startswith(_XZ_MAGIC):
        return lzma.LZMAFile(raw)
    return raw


@contextmanager
def open_text(path: Path) -> Iterator[TextIO]:
    """Open the file as UTF-8 text. The first bytes, not the suffix, show a compressed file.

    A directory, and compressed data that is damaged or truncated, raise `InputError`. The error
    for damaged data comes when the reader of the text reaches the damage.
    """
    shown = os.fspath(path)
    raw = open_binary(path)
    with raw:
        binary = _decompressed(raw)
        # pandas also opens a path with newline="", so its parser gets each line end unchanged.
        with io.TextIOWrapper(binary, encoding="utf-8", newline="") as text:
            try:
                yield text
            except _DAMAGED_DATA as exc:
                # A decompressor raises a bare OSError for damaged data, but a failed read of the
                # disk has an errno, and it stays as it is.
                if binary is raw or (isinstance(exc, OSError) and exc.errno is not None):
                    raise
                message = (
                    f"{shown!r}: the compressed data is damaged or truncated ({exc}); "
                    "expected a complete gzip, bzip2, or xz file"
                )
                raise InputError(message) from exc
