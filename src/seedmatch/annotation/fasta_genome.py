"""The FASTA genome, plain or BGZF: a reader that decodes only the bytes of the requested range."""

import io
import os
import struct
import zlib
from bisect import bisect_right
from collections.abc import Iterator
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path
from types import MappingProxyType
from typing import Final

import numpy as np
from numpy.typing import NDArray

from seedmatch._io import INT64_MAX, open_binary, parse_int64, too_large_error
from seedmatch.annotation.genome import Genome, duplicate_contig_error
from seedmatch.errors import InputError
from seedmatch.sequence.alphabet import MASKS, encode_dna

GZIP_MAGIC: Final = b"\x1f\x8b"
_GZIP_HEADER: Final = struct.Struct("<2sBBIBBH")
_SUBFIELD: Final = struct.Struct("<2sH")
_FEXTRA: Final = 4
_BLOCK_SIZE_TAG: Final = b"BC"
_BLOCK_SIZE_FIELD: Final = 2
_TRAILER_SIZE: Final = 8
_ISIZE: Final = struct.Struct("<I")
# The BGZF specification limits the uncompressed data of a block to 65536 bytes.
_BGZF_BLOCK_LIMIT: Final = 65536
_READ_SIZE: Final = 1 << 20
_FAI_COLUMNS: Final = 5
_FAI_NUMBERS: Final = ("length", "offset", "line bases", "line width")
_TERMINATOR_SIZES: Final = (1, 2)
_LINE_ENDS: Final = b"\r\n"
_LETTERS: Final = "".join(MASKS)


@dataclass(frozen=True, slots=True)
class _Entry:
    """The length of a contig, the byte of its first base, and the layout of its lines."""

    length: int
    offset: int
    line_bases: int
    line_width: int

    def position(self, base: int) -> int:
        """Return the byte of the base. The first base of the contig is base 0."""
        line, column = divmod(base, self.line_bases)
        return self.offset + line * self.line_width + column


class FastaGenome(Genome):
    """A genome in a FASTA file, plain or BGZF, with the index files of samtools.

    The file `<path>.fai` is necessary, and a BGZF file also needs `<path>.gzi`. The file is read
    through its indexes and no file is written. The reader gives the masks of the bases, so each
    soft-masked base is uppercase. A request reads the bytes of its range, and for BGZF it inflates
    the blocks of the range and keeps the last block.
    """

    def __init__(self, path: str | os.PathLike[str]) -> None:
        """Open the file, and read its indexes."""
        super().__init__(path)
        self._shown = os.fspath(self.path)
        self._fai = Path(f"{self._shown}.fai")
        self._gzi = Path(f"{self._shown}.gzi")
        self._compressed: list[int] = []
        self._uncompressed: list[int] = []
        self._cached: tuple[int, bytes] | None = None
        self._handle = open_binary(self.path)
        try:
            self._bgzf = self._handle.read(len(GZIP_MAGIC)) == GZIP_MAGIC
            if self._bgzf:
                self._read_block(0)
            self._entries = self._read_fai()
            if self._bgzf:
                self._read_gzi()
            self._check_extent()
        except Exception:
            self._handle.close()
            raise
        self._contigs = MappingProxyType(
            {name: entry.length for name, entry in self._entries.items()}
        )

    @property
    def contigs(self) -> MappingProxyType[str, int]:
        """Map the name of each contig to its length, in the order of the `.fai`."""
        return self._contigs

    def _codes(self, contig: str, start: int, end: int) -> NDArray[np.uint8]:
        """Return the masks of the range, from the bytes of the range only."""
        entry = self._entries[contig]
        first, last = entry.position(start - 1), entry.position(end - 1)
        chunks = self._bgzf_chunks(first, last) if self._bgzf else self._plain_chunks(first, last)
        total = end - start + 1
        codes = np.empty(total, dtype=np.uint8)
        filled = 0
        for chunk in chunks:
            bases = chunk.translate(None, _LINE_ENDS)
            if filled + len(bases) > total:
                raise self._wrong_length(contig, start, end, f"more than {_plural(total, 'base')}")
            codes[filled : filled + len(bases)] = self._encode(bases, contig, start + filled)
            filled += len(bases)
        if filled < total:
            raise self._wrong_length(contig, start, end, _plural(filled, "base"))
        return codes

    def close(self) -> None:
        """Close the file, and drop the kept block."""
        self._cached = None
        self._handle.close()

    def _encode(self, bases: bytes, contig: str, first: int) -> NDArray[np.uint8]:
        """Encode a piece of a range. The first base of the piece has the number `first`."""
        # Latin-1 maps each byte to one character, so a byte that is not ASCII reaches the encoder
        # as a letter that it rejects, and a position in the text is a position in the bytes.
        text = bases.decode("latin-1")
        try:
            return encode_dna(text)
        except InputError as exc:
            index = _first_invalid(text)
            message = (
                f"{self._shown!r}: contig {contig!r} base {first + index} is {text[index]!r}, "
                f"not an IUPAC DNA letter; expected one of {_LETTERS} in either case"
            )
            raise InputError(message) from exc

    def _wrong_length(self, contig: str, start: int, end: int, found: str) -> InputError:
        message = (
            f"{self._shown!r}: contig {contig!r}: the range {start}..{end} has {found} in the "
            f"file; expected {_plural(end - start + 1, 'base')}, as the .fai says"
        )
        return InputError(message)

    def _plain_chunks(self, first: int, last: int) -> Iterator[bytes]:
        """Yield the bytes `first` to `last` of the file, in pieces of at most `_READ_SIZE`."""
        self._handle.seek(first)
        remaining = last - first + 1
        while remaining > 0:
            chunk = self._handle.read(min(remaining, _READ_SIZE))
            if not chunk:
                return
            remaining -= len(chunk)
            yield chunk

    def _bgzf_chunks(self, first: int, last: int) -> Iterator[bytes]:
        """Yield the uncompressed bytes `first` to `last`, one piece for each block of them."""
        index = bisect_right(self._uncompressed, first) - 1
        while index < len(self._uncompressed) and self._uncompressed[index] <= last:
            start = self._uncompressed[index]
            yield self._inflate(index)[max(first - start, 0) : last - start + 1]
            index += 1

    def _read_block(self, offset: int) -> bytes:
        """Return the deflate data and the trailer of the BGZF block at the compressed offset.

        A block that the file cuts short, and each error of `_read_header`, raise `InputError`.
        """
        size, header_size = self._read_header(offset)
        body_size = size - header_size
        body = self._handle.read(body_size)
        if len(body) < body_size:
            raise self._damaged(offset, "the file ends inside it")
        return body

    def _read_header(self, offset: int) -> tuple[int, int]:
        """Return the size of the BGZF block at the compressed offset, and the size of its header.

        The file stays at the end of the header. A header that the file cuts short, a block
        without the extra field `BC`, and a size too small for the header and the trailer raise
        `InputError`.
        """
        self._handle.seek(offset)
        fixed = self._handle.read(_GZIP_HEADER.size)
        if len(fixed) < _GZIP_HEADER.size:
            raise self._cut_header(offset)
        magic, _method, flags, _time, _flags, _system, extra_size = _GZIP_HEADER.unpack(fixed)
        size = None
        if magic == GZIP_MAGIC and flags & _FEXTRA:
            extra = self._handle.read(extra_size)
            if len(extra) < extra_size:
                raise self._cut_header(offset)
            size = _block_size(extra)
        if size is None:
            message = (
                f"{self._shown!r}: the gzip block at compressed offset {offset} is not BGZF; "
                "expected the extra field BC that `bgzip` writes"
            )
            raise InputError(message)
        header_size = _GZIP_HEADER.size + extra_size
        if size - header_size < _TRAILER_SIZE:
            raise self._damaged(offset, f"its size {size} is too small for its header and trailer")
        return size, header_size

    def _inflate(self, index: int) -> bytes:
        """Return the uncompressed bytes of a block. The last block that it inflated is kept."""
        if self._cached is not None and self._cached[0] == index:
            return self._cached[1]
        offset = self._compressed[index]
        body = self._read_block(offset)
        try:
            block = zlib.decompress(body[:-_TRAILER_SIZE], wbits=-15)
        except zlib.error as exc:
            raise self._damaged(offset, f"its data does not inflate ({exc})") from exc
        (checksum,) = struct.unpack_from("<I", body, len(body) - _TRAILER_SIZE)
        if zlib.crc32(block) != checksum:
            raise self._damaged(offset, "its checksum does not match its data")
        if index + 1 < len(self._uncompressed):
            expected = self._uncompressed[index + 1] - self._uncompressed[index]
            if len(block) != expected:
                problem = (
                    f"the block at compressed offset {offset} has {_plural(len(block), 'byte')}, "
                    f"but the .gzi gives {expected}"
                )
                raise self._gzi_error(problem)
        self._cached = (index, block)
        return block

    def _gzi_error(self, problem: str) -> InputError:
        message = (
            f"{os.fspath(self._gzi)!r}: {problem}; expected the .gzi of this file. "
            f"Make it with `bgzip -r {self._shown}`"
        )
        return InputError(message)

    def _cut_header(self, offset: int) -> InputError:
        message = (
            f"{self._shown!r}: the file ends before the BGZF block at compressed offset {offset} "
            "has a complete header; expected an intact file and the .gzi of this file"
        )
        return InputError(message)

    def _damaged(self, offset: int, reason: str) -> InputError:
        message = (
            f"{self._shown!r}: the BGZF block at compressed offset {offset} is damaged: {reason}; "
            "expected an intact file"
        )
        return InputError(message)

    def _read_fai(self) -> dict[str, _Entry]:
        """Read the `.fai`. Raise `InputError` for a line that is not valid."""
        shown = os.fspath(self._fai)
        entries: dict[str, _Entry] = {}
        try:
            with io.TextIOWrapper(open_binary(self._fai), encoding="utf-8", newline="\n") as handle:
                for number, line in enumerate(handle, start=1):
                    name, entry = _parse_fai_line(shown, number, line.removesuffix("\n"))
                    if name in entries:
                        raise duplicate_contig_error(shown, name)
                    entries[name] = entry
        except FileNotFoundError:
            message = (
                f"{shown!r}: the index file does not exist; "
                f"make it with `samtools faidx {self._shown}`"
            )
            raise FileNotFoundError(message) from None
        except UnicodeDecodeError as exc:
            message = f"{shown!r}: the file is not valid UTF-8; expected UTF-8 contig names"
            raise InputError(message) from exc
        return entries

    def _check_extent(self) -> None:
        """Raise `InputError` for a contig that the `.fai` puts after the end of the data.

        The test is one comparison for each contig and reads no data. For BGZF the end of the data
        is not known before a read, so the limit is the start of the last block plus the largest
        size of a block.
        """
        if self._bgzf:
            limit = self._uncompressed[-1] + _BGZF_BLOCK_LIMIT
        else:
            limit = os.fstat(self._handle.fileno()).st_size
        for name, entry in self._entries.items():
            if not entry.length:
                continue
            last = entry.position(entry.length - 1)
            if last >= limit:
                message = (
                    f"{os.fspath(self._fai)!r}: contig {name!r} has its last base at byte {last}, "
                    f"but the file holds at most {_plural(limit, 'byte')} of data; the .fai does "
                    f"not match the file. Make it with `samtools faidx {self._shown}`"
                )
                raise InputError(message)

    def _read_gzi(self) -> None:
        """Read the `.gzi`: a count, then the compressed and uncompressed offset of each block.

        Offsets that do not increase, and the errors of `_check_uncompressed_offsets`, raise
        `InputError`.
        """
        shown = os.fspath(self._gzi)
        try:
            with open_binary(self._gzi) as handle:
                data = handle.read()
        except FileNotFoundError:
            message = (
                f"{shown!r}: the index file does not exist; make it with `bgzip -r {self._shown}`"
            )
            raise FileNotFoundError(message) from None
        count = int.from_bytes(data[:8], "little")
        expected = 8 + 16 * count
        if len(data) != expected:
            message = (
                f"{shown!r}: the file has {_plural(len(data), 'byte')}; "
                f"expected {_plural(expected, 'byte')} for {_plural(count, 'block')}"
            )
            raise InputError(message)
        values = struct.unpack_from(f"<{2 * count}Q", data, 8)
        offsets = values[0::2]
        size = os.fstat(self._handle.fileno()).st_size
        if offsets and max(offsets) >= size:
            offset = next(offset for offset in offsets if offset >= size)
            problem = (
                f"compressed offset {offset} is not below the size of the file, "
                f"{_plural(size, 'byte')}"
            )
            raise self._gzi_error(problem)
        # The first block has no entry in the file.
        self._compressed = [0, *offsets]
        self._uncompressed = [0, *values[1::2]]
        for column, column_offsets in (
            ("compressed", self._compressed),
            ("uncompressed", self._uncompressed),
        ):
            for previous, offset in pairwise(column_offsets):
                if offset <= previous:
                    problem = (
                        f"{column} offset {offset} is not above the offset before it, {previous}"
                    )
                    raise self._gzi_error(problem)
        self._check_uncompressed_offsets()

    def _check_uncompressed_offsets(self) -> None:
        """Make sure that each uncompressed offset is the sum of the ISIZE of the blocks before it.

        The test reads the header and the ISIZE field of each block before the last `.gzi` entry,
        and inflates no block. A `.gzi` entry inside a block raises `InputError`.
        """
        offset = total = 0
        for listed, expected in zip(self._compressed[1:], self._uncompressed[1:], strict=True):
            while offset < listed:
                size, _header_size = self._read_header(offset)
                end = offset + size
                if end > listed:
                    problem = (
                        f"compressed offset {listed} is inside the block at compressed offset "
                        f"{offset}; the next block starts at {end}"
                    )
                    raise self._gzi_error(problem)
                # The block ends at or before an entry that is below the file size, so the read
                # gets the 4 bytes.
                self._handle.seek(end - _ISIZE.size)
                (isize,) = _ISIZE.unpack(self._handle.read(_ISIZE.size))
                total += isize
                offset = end
            if total != expected:
                problem = (
                    f"uncompressed offset {expected} of the block at compressed offset {listed} is "
                    f"not {total}, the sum of the ISIZE fields of the blocks before it"
                )
                raise self._gzi_error(problem)


def _block_size(extra: bytes) -> int | None:
    """Return the size of a BGZF block from the extra field of its header, or `None` if no `BC`."""
    position = 0
    while position + _SUBFIELD.size <= len(extra):
        tag, size = _SUBFIELD.unpack_from(extra, position)
        data = position + _SUBFIELD.size
        if tag == _BLOCK_SIZE_TAG and size == _BLOCK_SIZE_FIELD and data + size <= len(extra):
            (block_size,) = struct.unpack_from("<H", extra, data)
            return int(block_size) + 1
        position = data + size
    return None


def _plural(count: int, noun: str) -> str:
    return f"{count} {noun}" if count == 1 else f"{count} {noun}s"


def _first_invalid(text: str) -> int:
    """Return the index of the first letter that `encode_dna` rejects. The text must have one."""
    valid, invalid = 0, len(text)
    while invalid - valid > 1:
        middle = (valid + invalid) // 2
        try:
            encode_dna(text[:middle])
        except InputError:
            invalid = middle
        else:
            valid = middle
    return valid


def _count(shown: str, number: int, label: str, text: str) -> int:
    """Return the number in a column of a `.fai` line. It must fit in a signed 64-bit integer."""
    if not (text.isascii() and text.isdecimal()):
        message = (
            f"{shown!r}: line {number}: {label} {text!r} is not a non-negative integer; "
            "expected digits only"
        )
        raise InputError(message)
    return parse_int64(shown, number, label, text)


def _parse_fai_line(shown: str, number: int, line: str) -> tuple[str, _Entry]:
    fields = line.split("\t")
    if len(fields) != _FAI_COLUMNS:
        message = (
            f"{shown!r}: line {number} has {_plural(len(fields), 'column')}; "
            f"expected {_FAI_COLUMNS} columns separated by tabs: "
            "name, length, offset, line bases, line width"
        )
        raise InputError(message)
    name, *texts = fields
    length, offset, line_bases, line_width = (
        _count(shown, number, label, text) for label, text in zip(_FAI_NUMBERS, texts, strict=True)
    )
    if line_bases < 1:
        message = (
            f"{shown!r}: line {number}: line bases {line_bases} is below 1; expected 1 or more"
        )
        raise InputError(message)
    if line_width - line_bases not in _TERMINATOR_SIZES:
        message = (
            f"{shown!r}: line {number}: line width {line_width} minus line bases {line_bases} is "
            f"{line_width - line_bases}; expected 1 (LF) or 2 (CRLF)"
        )
        raise InputError(message)
    entry = _Entry(length, offset, line_bases, line_width)
    last = entry.position(length - 1) if length else 0
    if last > INT64_MAX:
        raise too_large_error(shown, number, "the byte of the last base", str(last))
    return name, entry
