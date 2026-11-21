"""A set of sequences in one buffer."""

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Final, Self

import numpy as np
from numpy.typing import NDArray

from seedmatch._checks import check_offsets, check_vector, shown
from seedmatch.errors import InputError
from seedmatch.sequence.alphabet import MASKS, decode_dna, encode_dna_bytes

_MAX_MASK: Final = max(MASKS.values())


def _check_masks(codes: NDArray[np.uint8]) -> None:
    """Raise for the first code that is not a mask, because the kernels index tables with codes."""
    # min and max need no temporary array, so a large buffer costs no extra memory.
    if codes.size == 0 or (codes.min() >= 1 and codes.max() <= _MAX_MASK):
        return
    bad = int(np.flatnonzero((codes < 1) | (codes > _MAX_MASK))[0])
    message = f"invalid mask {codes[bad]} at index {bad}; expected 1 to {_MAX_MASK}"
    raise InputError(message)


def _check_text(field: str, values: tuple[str, ...]) -> None:
    """Raise for the first value that is not a str, because a DE `gene_id` is text."""
    for value in values:
        if not isinstance(value, str):
            message = f"{field} holds {shown(value)} of type {type(value).__name__}; expected a str"
            raise InputError(message)


@dataclass(frozen=True, eq=False)
class SequenceSet:
    """Sequences as masks in one buffer. Sequence i is codes[offsets[i]:offsets[i + 1]].

    The sequences of one gene share a gene ID. Each ID and each gene ID is a `str`. The
    constructor does not copy `codes` and `offsets`, but it makes them read-only. Thus the arrays
    of the caller also become read-only. The search kernels trust the values that the constructor
    checks. Thus a write through a base array after the call can give a wrong result or a crash.
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
        _check_text("ids", self.ids)
        _check_text("gene_ids", self.gene_ids)
        for name, array, dtype in (
            ("codes", self.codes, np.uint8),
            ("offsets", self.offsets, np.int64),
        ):
            check_vector(name, array, dtype)
            if not array.flags.c_contiguous:
                message = f"{name} is not C-contiguous; expected a C-contiguous array"
                raise InputError(message)
        check_offsets(self.offsets, n, "sequences", len(self.codes), "the length of codes")
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
        if isinstance(gene_ids, str):
            message = f"gene_ids {shown(gene_ids)} is a str; expected a sequence of str"
            raise InputError(message)
        ids: list[str] = []
        ends: list[int] = []
        # One growing buffer keeps the peak memory near the size of the result.
        buffer = bytearray()
        for name, text in records:
            try:
                buffer += encode_dna_bytes(text)
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
