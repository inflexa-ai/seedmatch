"""A set of sequences in one buffer."""

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Final, Self

import numpy as np
from numpy.typing import NDArray

from seedmatch.errors import InputError
from seedmatch.sequence.alphabet import MASKS, decode_dna, encode_dna

_MAX_MASK: Final = max(MASKS.values())


def _check_masks(codes: NDArray[np.uint8]) -> None:
    """Raise for the first code that is not a mask, because the kernels index tables with codes."""
    # min and max need no temporary array, so a large buffer costs no extra memory.
    if codes.size == 0 or (codes.min() >= 1 and codes.max() <= _MAX_MASK):
        return
    bad = int(np.flatnonzero((codes < 1) | (codes > _MAX_MASK))[0])
    message = f"invalid mask {codes[bad]} at index {bad}; expected 1 to {_MAX_MASK}"
    raise InputError(message)


@dataclass(frozen=True, eq=False)
class SequenceSet:
    """Sequences as masks in one buffer. Sequence i is codes[offsets[i]:offsets[i + 1]].

    The sequences of one gene share a gene ID. The arrays are read-only.
    """

    ids: tuple[str, ...]
    gene_ids: tuple[str, ...]
    codes: NDArray[np.uint8]
    offsets: NDArray[np.int64]

    def __post_init__(self) -> None:
        """Make sure that the parts fit together, because the search kernels trust the offsets."""
        n = len(self.ids)
        if len(self.gene_ids) != n:
            message = (
                f"gene_ids has {len(self.gene_ids)} entries; "
                f"expected one for each of the {n} sequences"
            )
            raise InputError(message)
        for name, array, dtype in (
            ("codes", self.codes, np.uint8),
            ("offsets", self.offsets, np.int64),
        ):
            if array.dtype != dtype:
                message = f"{name} has dtype {array.dtype}; expected {np.dtype(dtype)}"
                raise InputError(message)
            if array.ndim != 1:
                message = f"{name} has {array.ndim} dimensions; expected a 1-D array"
                raise InputError(message)
            if not array.flags.c_contiguous:
                message = f"{name} is not C-contiguous; expected a C-contiguous array"
                raise InputError(message)
        if len(self.offsets) != n + 1:
            message = f"offsets has {len(self.offsets)} entries; expected {n + 1} for {n} sequences"
            raise InputError(message)
        if self.offsets[0] != 0:
            message = f"offsets must start at 0; got {self.offsets[0]}"
            raise InputError(message)
        if self.offsets[-1] != len(self.codes):
            message = (
                f"offsets must end at {len(self.codes)}, the length of codes; "
                f"got {self.offsets[-1]}"
            )
            raise InputError(message)
        if (np.diff(self.offsets) < 0).any():
            message = "offsets must not decrease"
            raise InputError(message)
        _check_masks(self.codes)
        self.codes.flags.writeable = False
        self.offsets.flags.writeable = False

    @classmethod
    def from_records(
        cls,
        records: Iterable[tuple[str, str]],
        gene_ids: Sequence[str] | None = None,
    ) -> Self:
        """Encode (ID, DNA text) records. Without gene_ids, the gene ID of a sequence is its ID."""
        ids: list[str] = []
        ends: list[int] = []
        # One growing buffer keeps the peak memory near the size of the result.
        buffer = bytearray()
        for name, text in records:
            try:
                buffer += encode_dna(text).data
            except InputError as exc:
                message = f"record {name!r}: {exc}"
                raise InputError(message) from exc
            ids.append(name)
            ends.append(len(buffer))
        offsets = np.zeros(len(ids) + 1, dtype=np.int64)
        offsets[1:] = ends
        return cls(
            ids=tuple(ids),
            gene_ids=tuple(ids) if gene_ids is None else tuple(gene_ids),
            codes=np.frombuffer(buffer, dtype=np.uint8),
            offsets=offsets,
        )

    def __len__(self) -> int:
        """Return the number of sequences."""
        return len(self.ids)

    def sequence(self, i: int) -> str:
        """Return sequence i as uppercase text."""
        if not 0 <= i < len(self):
            expected = f"0 to {len(self) - 1}" if len(self) else "none, because the set is empty"
            message = f"sequence index {i} is outside the set; expected {expected}"
            raise InputError(message)
        return decode_dna(self.codes[self.offsets[i] : self.offsets[i + 1]])

    def lengths(self) -> NDArray[np.int64]:
        """Return the length of each sequence."""
        return np.diff(self.offsets)
