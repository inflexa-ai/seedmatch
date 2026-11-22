"""FASTA files."""

import os
import re
from itertools import chain
from pathlib import Path
from typing import Final

import numpy as np

from seedmatch._checks import as_limit, check_at_least
from seedmatch._io import open_text
from seedmatch.errors import InputError
from seedmatch.sequence.alphabet import encode_dna_bytes
from seedmatch.sequence.sequence_set import SequenceSet

_WHITE_SPACE = re.compile(r"\s")
# One encode call for many lines keeps the throughput, and the bound keeps their memory small.
_BATCH_CHARACTERS: Final = 1 << 18


def _add(codes: bytearray, lines: list[str], first: int, path: Path, name: str) -> None:
    """Add the masks of the lines to `codes`, and clear `lines`. Line `first` is `lines[0]`."""
    # The text reader ends each line with LF only, thus LF is the only line end to remove.
    try:
        codes += encode_dna_bytes("".join(lines).replace("\n", ""))
    except InputError:
        # One line at a time, to name the line of the first bad letter.
        for number, line in enumerate(lines, start=first):
            try:
                codes += encode_dna_bytes(line.rstrip("\n"))
            except InputError as exc:
                message = f"{path}, line {number}: record {name!r}: {exc}"
                raise InputError(message) from exc
    lines.clear()


def _read(path: Path) -> SequenceSet:
    """Read the file one line at a time, and add the masks of the sequence lines to one buffer."""
    ids: list[str] = []
    offsets = [0]
    codes = bytearray()
    lines: list[str] = []
    first = size = 0
    with open_text(path, newline=None) as handle:
        # The skip of the mark (D-071). The codec utf-8-sig is not used, because it decodes an
        # incomplete mark at the end of the file to no text, with no error.
        start = handle.readline().removeprefix("\ufeff")
        for number, line in enumerate(chain((start,) if start else (), handle), start=1):
            if line.startswith(">"):
                if ids:
                    _add(codes, lines, first, path, ids[-1])
                    offsets.append(len(codes))
                ids.append(_WHITE_SPACE.split(line[1:], maxsplit=1)[0])
                size = 0
            elif ids:
                if not lines:
                    first = number
                lines.append(line)
                size += len(line)
                if size >= _BATCH_CHARACTERS:
                    _add(codes, lines, first, path, ids[-1])
                    size = 0
            elif line != "\n":
                message = (
                    f"{path}, line {number}: sequence line before the first header; "
                    "expected a header line that starts with '>'"
                )
                raise InputError(message)
    if ids:
        _add(codes, lines, first, path, ids[-1])
        offsets.append(len(codes))
    return SequenceSet(
        ids=tuple(ids),
        gene_ids=tuple(ids),
        codes=np.frombuffer(codes, dtype=np.uint8),
        offsets=np.array(offsets, dtype=np.int64),
    )


def read_fasta(path: str | os.PathLike[str]) -> SequenceSet:
    """Read a DNA FASTA file. The ID of a record is its header text up to the first white space.

    The file is UTF-8 text, or gzip, bzip2, or xz of UTF-8 text. The first bytes of the file, not
    the suffix, show the compression. R reads only gzip, and it reads a truncated gzip file up to
    the damage, with no error. Here, damaged or truncated compressed data raises `InputError`
    (`D-074`).

    The reader skips a byte order mark at the start of the file (`D-071`). A letter that is not
    DNA raises `InputError` with the path, the line, the position in the line, and the record.

    The reader encodes the lines in batches of about 2^18 characters into one buffer of masks.
    Thus the peak memory is the size of the bases and about 1 MiB. A line longer than a batch adds
    about three times its length. For one record of 16 Mi bases in lines of 60, the peak of
    `tracemalloc` is 1.2 times the bases, and the call takes 0.065 seconds. Machine: Apple M2 Max,
    Python 3.12.
    """
    path = Path(path)
    try:
        return _read(path)
    except UnicodeDecodeError as exc:
        message = (
            f"{path} is not valid UTF-8 text ({exc.reason}); "
            "expected UTF-8 text, or gzip, bzip2, or xz of UTF-8 text"
        )
        raise InputError(message) from exc


def write_fasta(
    sequences: SequenceSet,
    path: str | os.PathLike[str],
    *,
    width: int = 80,
) -> None:
    """Write the IDs and bases of the set. Each sequence line has at most width letters.

    A `width` that is not an integer of 1 or more raises `InputError`.
    """
    width = as_limit("width", width)
    check_at_least("width", width, 1)
    for name in sequences.ids:
        if _WHITE_SPACE.search(name):
            message = f"id {name!r} has white space; expected an ID without white space"
            raise InputError(message)
    with Path(path).open("w", encoding="utf-8", newline="\n") as handle:
        for i, name in enumerate(sequences.ids):
            text = sequences.sequence(i)
            handle.write(f">{name}\n")
            handle.writelines(
                f"{text[start : start + width]}\n" for start in range(0, len(text), width)
            )
