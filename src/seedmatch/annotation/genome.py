"""The genome interface: the contigs, and the bases of a range of a contig."""

import os
from abc import ABC, abstractmethod
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Final, Self

import numpy as np
from numpy.typing import NDArray

from seedmatch._checks import as_limit
from seedmatch._io import local_path
from seedmatch.errors import InputError
from seedmatch.sequence.alphabet import decode_dna

_LISTED_CONTIGS: Final = 10
_MITOCHONDRIAL: Final = {"chrM": "MT", "MT": "chrM"}


def duplicate_contig_error(shown: str, contig: str) -> InputError:
    """Return the error for a contig name that occurs twice in the index of a genome file."""
    message = f"{shown!r}: contig {contig!r} occurs twice in the index; expected each name once"
    return InputError(message)


def listed_names(names: Iterable[str]) -> str:
    """Return the first names in sorted order as text for an error message, with `...` for more."""
    ordered = sorted(names)
    listed = ", ".join(ordered[:_LISTED_CONTIGS])
    return f"{listed}, ..." if len(ordered) > _LISTED_CONTIGS else listed


class Genome(ABC):
    """The sequences of the contigs of one genome. A range has 1-based inclusive coordinates."""

    path: Path

    def __init__(self, path: str | os.PathLike[str]) -> None:
        """Set the path of the genome file. A URL raises `InputError`."""
        self.path = local_path(path)

    @property
    @abstractmethod
    def contigs(self) -> Mapping[str, int]:
        """Map the name of each contig to its length."""

    def codes(self, contig: str, start: int, end: int) -> NDArray[np.uint8]:
        """Return the masks of the bases `start` to `end` of the contig.

        The masks do not mark the soft-masked bases of the file. `end == start - 1` gives an empty
        array. An unknown contig, and a range outside the contig, raise `InputError`.
        """
        first, last = self._check_range(contig, start, end)
        if last < first:
            return np.empty(0, dtype=np.uint8)
        return self._codes(contig, first, last)

    @abstractmethod
    def _codes(self, contig: str, start: int, end: int) -> NDArray[np.uint8]:
        """Return the masks of a range that `codes` has checked.

        The contig exists, and `1 <= start <= end <= length` holds, so the range is not empty.
        """

    @abstractmethod
    def close(self) -> None:
        """Release the file of the genome."""

    def sequence(self, contig: str, start: int, end: int) -> str:
        """Return the bases of the range as uppercase text."""
        return decode_dna(self.codes(contig, start, end))

    def resolve(self, name: str) -> str | None:
        """Return the contig of the genome for the name of an annotation contig, or `None`.

        The result is the first of these names that the genome has: the same name, the name without
        `chr`, the name with `chr`, and the other name of the mitochondrial contig (`chrM` and `MT`
        are the same contig).
        """
        candidates = (name, name.removeprefix("chr"), f"chr{name}", _MITOCHONDRIAL.get(name, name))
        contigs = self.contigs
        return next((candidate for candidate in candidates if candidate in contigs), None)

    def __enter__(self) -> Self:
        """Return the genome, so that `with` closes it."""
        return self

    def __exit__(self, *_exc_info: object) -> None:
        """Close the genome."""
        self.close()

    def _check_range(self, contig: str, start: int, end: int) -> tuple[int, int]:
        """Return the bounds as plain ints. Raise `InputError` for an unknown contig or a bad range.

        An empty range, `end == start - 1`, is valid for an `end` from 0 to the length.
        """
        first = as_limit("start", start)
        last = as_limit("end", end)
        length = self.contigs.get(contig)
        if length is None:
            listed = listed_names(self.contigs)
            message = f"contig {contig!r} is not in the genome; expected one of {listed}"
            raise InputError(message)
        if last < first - 1:
            message = (
                f"contig {contig!r}: range {first}..{last} ends before it starts; "
                "expected end >= start - 1 (end == start - 1 is an empty range)"
            )
            raise InputError(message)
        if first < 1 or last > length:
            message = (
                f"contig {contig!r}: range {first}..{last} is outside the contig of length "
                f"{length}; expected 1 <= start and end <= {length}"
            )
            raise InputError(message)
        return first, last
