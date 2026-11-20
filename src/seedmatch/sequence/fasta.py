"""FASTA files."""

import io
import os
import re
from collections.abc import Iterator
from pathlib import Path

from seedmatch._io import open_binary
from seedmatch.errors import InputError
from seedmatch.sequence.sequence_set import SequenceSet

_WHITE_SPACE = re.compile(r"\s")


def _records(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (ID, text) for each record, and read the file one line at a time."""
    name: str | None = None
    lines: list[str] = []
    with io.TextIOWrapper(open_binary(path), encoding="utf-8") as handle:
        for number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(lines)
                name, lines = _WHITE_SPACE.split(line[1:], maxsplit=1)[0], []
            elif name is None:
                message = (
                    f"{path}, line {number}: sequence line before the first header; "
                    "expected a header line that starts with '>'"
                )
                raise InputError(message)
            else:
                lines.append(line)
    if name is not None:
        yield name, "".join(lines)


def read_fasta(path: str | os.PathLike[str]) -> SequenceSet:
    """Read a DNA FASTA file. The ID of a record is its header text up to the first white space."""
    path = Path(path)
    try:
        return SequenceSet.from_records(_records(path))
    except UnicodeDecodeError as exc:
        message = f"{path} is not valid UTF-8 text ({exc.reason}); expected a text FASTA file"
        raise InputError(message) from exc


def write_fasta(
    sequences: SequenceSet,
    path: str | os.PathLike[str],
    *,
    width: int = 80,
) -> None:
    """Write the IDs and bases of the set. Each sequence line has at most width letters."""
    if width < 1:
        message = f"width {width} is invalid; expected at least 1"
        raise InputError(message)
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
