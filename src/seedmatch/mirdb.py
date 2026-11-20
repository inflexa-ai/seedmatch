"""The reader of the miRDB table of predicted miRNA targets."""

import csv
import io
import os
import re
from dataclasses import dataclass, field
from functools import partial
from itertools import compress, islice, repeat
from typing import Final, TextIO

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from seedmatch._checks import shown
from seedmatch._io import local_path, open_text
from seedmatch.errors import InputError

_FIELDS: Final = ("mirdb_id", "target_refseq_id", "target_score")
_CHUNK_LINES: Final = 2**18
_LINE_LIMIT: Final = 2**10
_DECIMAL: Final = re.compile(r" *[+-]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][+-]?[0-9]+)? *")


def read_mirdb(path: str | os.PathLike[str], *, species: str | None = None) -> pd.DataFrame:
    """Read the miRDB table of predicted miRNA targets, and split each miRNA ID.

    The path is a local path: a URL raises `InputError`, because the core modules do not use the
    network. The file is UTF-8 text, or gzip, bzip2, or xz of UTF-8 text, and the first bytes, not
    the suffix, show the compression. It is the file that `seedmatch-fetch mirdb` downloads. The
    file has no header. Each line has 3 fields with a tab between them: the miRDB ID, the RefSeq
    ID of the target, and the score. A blank line is skipped.

    The columns are `mirdb_id`, `target_refseq_id`, `target_score` (`float64`), `species`, and
    `mirna_name`, as `R/get_data.R` of SeedMatchR names them, and the rows keep the order of the
    file. The text columns have the dtype `object`. `species` is the text of the ID before the
    first `-`, and `mirna_name` is the text after it. An ID without `-` gives the full ID as
    `species` and an empty `mirna_name`. With `species`, the result keeps only the rows of this
    species code. A code without a row raises `InputError` that lists the codes of the file.

    The score is a decimal number with an optional sign and exponent, for example `96.4` or
    `1e2`, and it is finite. Spaces before and after the number are valid, as in R. Another score
    raises `InputError` that names the line, also for a row of another species (`D-049`). The ID
    fields stay as they are: a quote character, a comment character, and the text `NA` have no
    special meaning (`D-049`).

    These inputs also raise `InputError`:

    - a line with another number of fields
    - a line of more than 1024 characters with its line end (`D-049`)
    - a line with a NUL character (`D-049`)
    - text that is not valid UTF-8
    - a file without a data line

    The limit of 1024 characters bounds the memory of a chunk in bytes. R reads a longer line.
    The longest line of the miRDB 6.0 file of `seedmatch-fetch mirdb` has 51 characters with its
    line end. R skips a byte order mark at the start of the file, and so does the reader. The line
    number in a message counts the blank lines.

    The reader reads 2^18 lines at a time, and the time is linear in the number of lines. Equal
    texts in the result are one object, thus the result holds about 40 bytes for each row and the
    distinct strings. The miRDB 6.0 file of `seedmatch-fetch mirdb` has 6,831,595 lines. A read of
    all species takes 6.2 seconds, and the peak memory above the start is about 670 MB. With one
    species code, the read takes 5.8 seconds, and the peak is about 155 MB.
    """
    source = os.fspath(path)
    file = local_path(source)
    wanted = _as_species(species)
    rows = _Rows()
    try:
        with open_text(file) as handle:
            _read_chunks(handle, source, wanted, rows)
    except UnicodeDecodeError as exc:
        message = (
            f"{source!r}: the file is not valid UTF-8 text ({exc.reason}); "
            "expected a text or gzip miRDB file"
        )
        raise InputError(message) from exc
    return rows.table(source, wanted)


def _as_species(species: object) -> str | None:
    if species is None or isinstance(species, str):
        return species
    message = f"species {shown(species)} has type {type(species).__name__}; expected a str or None"
    raise InputError(message)


@dataclass
class _Rows:
    """The kept rows of the chunks read so far.

    A text column holds, for each row, the code of its text in a vocabulary. Thus the chunks keep
    one string for each distinct text, not one for each row. `species` holds the species codes of
    the data lines, and only a read with a species code fills it.
    """

    ids: dict[str, int] = field(default_factory=dict)
    refseqs: dict[str, int] = field(default_factory=dict)
    id_codes: list[NDArray[np.int32]] = field(default_factory=list)
    refseq_codes: list[NDArray[np.int32]] = field(default_factory=list)
    scores: list[NDArray[np.float64]] = field(default_factory=list)
    species: set[str] = field(default_factory=set)
    data_lines: int = 0

    def add(
        self,
        chunk: list[str],
        is_data: NDArray[np.bool_],
        first: int,
        source: str,
        wanted: str | None,
    ) -> None:
        """Parse the data lines of the chunk, and keep the rows of the species `wanted`.

        `first` is the number of the first line of the chunk. `wanted` is `None` for all species.
        """
        numbers = first + np.flatnonzero(is_data)
        if not numbers.size:
            return
        flags = is_data.tolist()
        text = "".join(compress(chunk, flags))
        if "\0" in text:
            bad = next(i for i, line in enumerate(compress(chunk, flags)) if "\0" in line)
            message = (
                f"{source!r}: line {numbers[bad]}: the line has a NUL character (U+0000); "
                "expected text without NUL characters"
            )
            raise InputError(message)
        # Bytes, not a StringIO: pandas holds a text buffer in 4-byte characters, and the peak
        # memory of a chunk is then twice as large.
        data = text.encode()
        del text
        frame = pd.read_csv(
            io.BytesIO(data),
            sep="\t",
            header=None,
            names=list(_FIELDS),
            dtype=str,
            quoting=csv.QUOTE_NONE,
            na_filter=False,
        )
        scores = _parse_scores(frame["target_score"].tolist(), numbers, source)
        self.data_lines += len(frame)
        if wanted is not None:
            codes, texts = pd.factorize(frame["mirdb_id"])
            species_of_id = [text.partition("-")[0] for text in texts]
            self.species.update(species_of_id)
            keep = np.array([code == wanted for code in species_of_id], dtype=np.bool_)[codes]
            frame, scores = frame[keep], scores[keep]
        self.id_codes.append(_vocabulary_codes(frame["mirdb_id"], self.ids))
        self.refseq_codes.append(_vocabulary_codes(frame["target_refseq_id"], self.refseqs))
        self.scores.append(scores)

    def table(self, source: str, wanted: str | None) -> pd.DataFrame:
        """Return the table of the kept rows. A species without a row raises `InputError`."""
        if not self.data_lines:
            message = f"{source!r}: the file has no line; expected miRDB lines"
            raise InputError(message)
        if wanted is not None and wanted not in self.species:
            codes = ", ".join(sorted(self.species))
            message = f"{source!r}: no row has the species {wanted!r}; expected one of {codes}"
            raise InputError(message)
        ids = list(self.ids)
        parts = [text.partition("-") for text in ids]
        id_codes = _joined(self.id_codes)
        mirdb_id = np.array(ids, dtype=object)[id_codes]
        species = _shared([part[0] for part in parts])[id_codes]
        mirna_name = _shared([part[2] for part in parts])[id_codes]
        del id_codes
        refseq_id = np.array(list(self.refseqs), dtype=object)[_joined(self.refseq_codes)]
        columns = {
            "mirdb_id": mirdb_id,
            "target_refseq_id": refseq_id,
            "target_score": _joined(self.scores),
            "species": species,
            "mirna_name": mirna_name,
        }
        # A dict of object arrays makes pandas scan each array for dates, which allocates a copy of
        # the table. A Series with its dtype set skips the scan.
        return pd.DataFrame(
            {
                name: pd.Series(values, dtype=values.dtype, copy=False)
                for name, values in columns.items()
            },
            copy=False,
        )


def _read_chunks(handle: TextIO, source: str, wanted: str | None, rows: _Rows) -> None:
    next_line = partial(handle.readline, _LINE_LIMIT + 1)
    first = 1
    while chunk := list(islice(iter(next_line, ""), _CHUNK_LINES)):
        if first == 1:
            # R skips a byte order mark, also before a blank first line.
            chunk[0] = chunk[0].removeprefix("\ufeff")
        rows.add(chunk, _data_lines(chunk, first, source), first, source, wanted)
        first += len(chunk)


def _data_lines(chunk: list[str], first: int, source: str) -> NDArray[np.bool_]:
    """Return which lines of the chunk are data lines. A bad line raises `InputError`."""
    tabs = np.fromiter(map(str.count, chunk, repeat("\t")), dtype=np.int64, count=len(chunk))
    lengths = np.fromiter(map(len, chunk), dtype=np.int64, count=len(chunk))
    is_data: NDArray[np.bool_] = tabs == len(_FIELDS) - 1
    for offset in np.flatnonzero(~is_data | (lengths > _LINE_LIMIT)):
        line, number = chunk[offset], first + int(offset)
        if len(line) > _LINE_LIMIT:
            message = (
                f"{source!r}: line {number}: the line has more than {_LINE_LIMIT} characters; "
                f"expected at most {_LINE_LIMIT}, with the line end"
            )
            raise InputError(message)
        text = line.rstrip("\r\n")
        if text:
            count = text.count("\t") + 1
            noun = "field" if count == 1 else "fields"
            message = (
                f"{source!r}: line {number} has {count} {noun}; expected {len(_FIELDS)}: "
                "the miRDB ID, the RefSeq ID, and the score"
            )
            raise InputError(message)
    return is_data


def _parse_scores(texts: list[str], numbers: NDArray[np.int64], source: str) -> NDArray[np.float64]:
    """Return the scores. The first text that is not a finite decimal number raises `InputError`."""
    if all(map(_DECIMAL.fullmatch, texts)):
        scores = np.array(texts, dtype=object).astype(np.float64)
        # The pattern accepts an exponent that is too large, and float() gives inf for it.
        infinite = np.flatnonzero(~np.isfinite(scores))
        if not infinite.size:
            return scores
        bad = int(infinite[0])
    else:
        bad = next(i for i, text in enumerate(texts) if _DECIMAL.fullmatch(text) is None)
    message = (
        f"{source!r}: line {numbers[bad]}: score {texts[bad]!r} is not a finite number; "
        "expected a decimal number"
    )
    raise InputError(message)


def _vocabulary_codes(column: pd.Series, vocabulary: dict[str, int]) -> NDArray[np.int32]:
    """Return the code of each text of the column. A new text takes the next code."""
    codes, texts = pd.factorize(column)
    lookup = np.fromiter(
        (vocabulary.setdefault(text, len(vocabulary)) for text in texts),
        dtype=np.int32,
        count=len(texts),
    )
    return lookup[codes]


def _joined[T: np.generic](parts: list[NDArray[T]]) -> NDArray[T]:
    """Return the parts as one array. The list is empty after the call, so the parts are freed."""
    joined = np.concatenate(parts)
    parts.clear()
    return joined


def _shared(texts: list[str]) -> NDArray[np.object_]:
    """Return the texts as an object array in which equal texts are one object."""
    objects: dict[str, str] = {}
    return np.array([objects.setdefault(text, text) for text in texts], dtype=object)
