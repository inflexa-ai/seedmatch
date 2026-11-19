"""The 2bit genome of UCSC: a reader that decodes only the bytes of the requested range."""

import mmap
import os
import struct
from dataclasses import dataclass
from types import MappingProxyType
from typing import Final

import numpy as np
from numpy.typing import NDArray

from seedmatch.annotation.genome import Genome
from seedmatch.errors import InputError
from seedmatch.sequence.alphabet import MASKS

_SIGNATURES: Final = {order: struct.pack(f"{order}I", 0x1A412743) for order in "<>"}
_HEADER_SIZE: Final = 16
_BASES_PER_BYTE: Final = 4
# Version 1 has 64-bit offsets in the index. Each other number of the file is 32-bit.
_OFFSET_FORMATS: Final = {0: "I", 1: "Q"}
_N: Final = MASKS["N"]

# The 2-bit codes T, C, A, G are 0 to 3, and the first base of a byte is in the high bits.
_CODE_MASKS: Final = np.array([MASKS[base] for base in "TCAG"], dtype=np.uint8)
_UNPACK: Final = _CODE_MASKS[(np.arange(256)[:, None] >> np.array([6, 4, 2, 0])) & 3]


@dataclass(frozen=True, slots=True)
class _Record:
    """Where the packed bases of a contig start, and its N blocks as 0-based half-open spans."""

    dna: int
    n_starts: NDArray[np.int64]
    n_ends: NDArray[np.int64]


class TwoBitGenome(Genome):
    """A genome in a 2bit file of version 0 or 1, in either byte order.

    The file is mapped read-only, and no file is written. The reader gives the masks of the bases,
    and it ignores the mask blocks, so each soft-masked base is uppercase. The N blocks of a contig
    are read at the first request of the contig.
    """

    def __init__(self, path: str | os.PathLike[str]) -> None:
        """Open the file, and read its header and index."""
        super().__init__(path)
        self._shown = os.fspath(self.path)
        self._records: dict[str, _Record] = {}
        with self.path.open("rb") as handle:
            self._order = self._byte_order(handle.read(4))
            self._data = mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ)
        try:
            self._offsets, lengths = self._read_index()
        except Exception:
            self._data.close()
            raise
        self._contigs = MappingProxyType(lengths)

    @property
    def contigs(self) -> MappingProxyType[str, int]:
        """Map the name of each contig to its length, in the order of the index."""
        return self._contigs

    def _codes(self, contig: str, start: int, end: int) -> NDArray[np.uint8]:
        """Return the masks of the range. The N blocks that overlap it give the mask 15."""
        record = self._record(contig)
        begin, stop = start - 1, end
        first_byte = begin // _BASES_PER_BYTE
        last_byte = (stop - 1) // _BASES_PER_BYTE
        packed = self._data[record.dna + first_byte : record.dna + last_byte + 1]
        skipped = begin - first_byte * _BASES_PER_BYTE
        unpacked: NDArray[np.uint8] = _UNPACK[np.frombuffer(packed, dtype=np.uint8)].reshape(-1)
        codes = unpacked[skipped : skipped + stop - begin]
        low = int(np.searchsorted(record.n_ends, begin, side="right"))
        high = int(np.searchsorted(record.n_starts, stop, side="left"))
        for n_start, n_end in zip(
            record.n_starts[low:high].tolist(), record.n_ends[low:high].tolist(), strict=True
        ):
            codes[max(n_start, begin) - begin : min(n_end, stop) - begin] = _N
        return codes

    def close(self) -> None:
        """Unmap the file."""
        self._data.close()

    def _byte_order(self, head: bytes) -> str:
        """Return the struct prefix of the byte order that the signature shows."""
        for order, signature in _SIGNATURES.items():
            if head == signature:
                return order
        expected = " or ".join(signature.hex() for signature in _SIGNATURES.values())
        message = (
            f"{self._shown!r}: not a 2bit file; the first bytes are {head.hex()}, "
            f"expected {expected}"
        )
        raise InputError(message)

    def _require(self, end: int, what: str) -> None:
        """Raise `InputError` if the file ends before the byte `end`."""
        if end > len(self._data):
            message = (
                f"{self._shown!r}: the file ends before {what}; "
                f"expected at least {end} bytes, the file has {len(self._data)}"
            )
            raise InputError(message)

    def _read_fields(self, fmt: str, offset: int, what: str) -> tuple[int, ...]:
        fmt = self._order + fmt
        self._require(offset + struct.calcsize(fmt), what)
        return struct.unpack_from(fmt, self._data, offset)

    def _words(self, offset: int, count: int, what: str) -> NDArray[np.int64]:
        """Read `count` uint32 values from the offset."""
        self._require(offset + 4 * count, what)
        raw = self._data[offset : offset + 4 * count]
        return np.frombuffer(raw, dtype=np.dtype(self._order + "u4")).astype(np.int64)

    def _read_index(self) -> tuple[dict[str, int], dict[str, int]]:
        """Read the header and the index. Return the offset and the length of each contig."""
        self._require(_HEADER_SIZE, "the header")
        version, count = self._read_fields("II", 4, "the header")
        if version not in _OFFSET_FORMATS:
            message = f"{self._shown!r}: the 2bit file has version {version}; expected 0 or 1"
            raise InputError(message)
        offset_format = _OFFSET_FORMATS[version]
        offsets: dict[str, int] = {}
        lengths: dict[str, int] = {}
        position = _HEADER_SIZE
        for _ in range(count):
            (name_size,) = self._read_fields("B", position, "the index")
            self._require(position + 1 + name_size, "the index")
            raw_name = self._data[position + 1 : position + 1 + name_size]
            try:
                name = raw_name.decode("utf-8")
            except UnicodeDecodeError as exc:
                message = f"{self._shown!r}: contig name {raw_name!r} is not valid UTF-8"
                raise InputError(message) from exc
            if name in offsets:
                message = (
                    f"{self._shown!r}: contig {name!r} occurs twice in the index; "
                    "expected each name once"
                )
                raise InputError(message)
            position += 1 + name_size
            (offset,) = self._read_fields(offset_format, position, "the index")
            position += struct.calcsize(offset_format)
            (length,) = self._read_fields("I", offset, f"the record of contig {name!r}")
            offsets[name] = offset
            lengths[name] = length
        return offsets, lengths

    def _record(self, contig: str) -> _Record:
        """Read the N blocks of the contig at its first request, and check that its bases exist."""
        record = self._records.get(contig)
        if record is None:
            record = self._read_record(contig)
            self._records[contig] = record
        return record

    def _read_record(self, contig: str) -> _Record:
        what = f"the record of contig {contig!r}"
        offset = self._offsets[contig]
        (n_count,) = self._read_fields("I", offset + 4, what)
        blocks = self._words(offset + 8, 2 * n_count, what)
        mask_at = offset + 8 + 8 * n_count
        (mask_count,) = self._read_fields("I", mask_at, what)
        # The mask blocks are a count, the starts, and the sizes. A reserved uint32 follows them.
        dna = mask_at + 4 + 8 * mask_count + 4
        packed_size = -(-self.contigs[contig] // _BASES_PER_BYTE)
        self._require(dna + packed_size, f"the DNA of contig {contig!r}")
        return _Record(dna, blocks[:n_count], blocks[:n_count] + blocks[n_count:])
