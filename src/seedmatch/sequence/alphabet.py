"""The IUPAC DNA alphabet as 4-bit masks: bit 0 is A, bit 1 is C, bit 2 is G, bit 3 is T."""

from collections.abc import Mapping
from types import MappingProxyType
from typing import Final

import numpy as np
from numpy.typing import NDArray

from seedmatch.errors import InputError

_A, _C, _G, _T = 1, 2, 4, 8

MASKS: Final[Mapping[str, int]] = MappingProxyType(
    {
        "A": _A,
        "C": _C,
        "G": _G,
        "T": _T,
        "M": _A | _C,
        "R": _A | _G,
        "W": _A | _T,
        "S": _C | _G,
        "Y": _C | _T,
        "K": _G | _T,
        "V": _A | _C | _G,
        "H": _A | _C | _T,
        "D": _A | _G | _T,
        "B": _C | _G | _T,
        "N": _A | _C | _G | _T,
    }
)

_COMPLEMENT: Final[Mapping[str, str]] = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "M": "K",
    "K": "M",
    "R": "Y",
    "Y": "R",
    "V": "B",
    "B": "V",
    "H": "D",
    "D": "H",
    "W": "W",
    "S": "S",
    "N": "N",
}

_TABLE_SIZE: Final = 256
_LETTERS: Final = "".join(MASKS)


def _byte_table(values: Mapping[str, int]) -> bytes:
    """Map the byte of each letter, in either case, to its value. Each other byte maps to 0."""
    table = bytearray(_TABLE_SIZE)
    for letter, value in values.items():
        table[ord(letter)] = value
        table[ord(letter.lower())] = value
    return bytes(table)


_ENCODE_TABLE: Final = _byte_table(MASKS)
_COMPLEMENT_TABLE: Final = _byte_table({letter: ord(_COMPLEMENT[letter]) for letter in MASKS})


def _decode_table() -> bytes:
    """Map each mask to the byte of its letter. Each other value maps to 0."""
    table = bytearray(_TABLE_SIZE)
    for letter, mask in MASKS.items():
        table[mask] = ord(letter)
    return bytes(table)


_DECODE_TABLE: Final = _decode_table()


def _translate(text: str, table: bytes) -> bytes:
    """Map each letter of the text through the table. Raise for the first letter that is not DNA."""
    # A character outside ASCII becomes "?", which the tables map to 0, so the position stays the
    # position in the text.
    translated = text.encode("ascii", errors="replace").translate(table)
    bad = translated.find(0)
    if bad >= 0:
        message = (
            f"invalid DNA letter {text[bad]!r} at position {bad + 1}; "
            f"expected one of {_LETTERS} in either case"
        )
        raise InputError(message)
    return translated


def encode_dna(text: str) -> NDArray[np.uint8]:
    """Encode IUPAC DNA text as one mask per letter. Lowercase letters give the uppercase mask."""
    return np.frombuffer(_translate(text, _ENCODE_TABLE), dtype=np.uint8).copy()


def decode_dna(codes: NDArray[np.uint8]) -> str:
    """Decode masks to uppercase IUPAC DNA text."""
    translated = codes.tobytes().translate(_DECODE_TABLE)
    bad = translated.find(0)
    if bad >= 0:
        message = f"invalid mask {codes[bad]} at position {bad + 1}; expected 1 to 15"
        raise InputError(message)
    return translated.decode("ascii")


def reverse_complement(text: str) -> str:
    """Return the reverse complement of IUPAC DNA text, in uppercase."""
    return _translate(text, _COMPLEMENT_TABLE)[::-1].decode("ascii")
