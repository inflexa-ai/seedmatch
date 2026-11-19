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


def local_path(path: str | os.PathLike[str]) -> Path:
    """Return the path as a `Path`.

    A URL raises `InputError`, because the core modules do not use the network.
    """
    shown = os.fspath(path)
    if _URL.match(shown):
        message = f"{shown!r}: the path is a URL; expected a local path"
        raise InputError(message)
    return Path(shown)


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
    try:
        raw = path.open("rb")
    except IsADirectoryError as exc:
        message = f"{shown!r}: the path is a directory; expected a file"
        raise InputError(message) from exc
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
