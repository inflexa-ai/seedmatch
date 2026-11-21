"""The reader of the STRING protein links, and the cutoff of their scores."""

import os
from dataclasses import dataclass, field
from typing import Final

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from seedmatch._checks import as_cutoff, as_limit, concatenate_and_clear
from seedmatch._io import local_path, open_text
from seedmatch._lines import (
    field_lines,
    frame_of,
    not_utf8,
    parse_data_lines,
    read_chunks,
    vocabulary_codes,
)
from seedmatch._rmath import r_quantile
from seedmatch.errors import InputError

STRING_HEADER: Final = "protein1 protein2 combined_score"
STRING_SCORE_MAX: Final = 1000
_COLUMNS: Final = ("protein1", "protein2", "combined_score")
_CHUNK_LINES: Final = 2**16


@dataclass(frozen=True, eq=False)
class StringLinks:
    """The links of the STRING network with a score above the cutoff, and the cutoff.

    `links` has the columns `protein1` and `protein2` (`object`, equal texts are one object) and
    `combined_score` (`int64`), and a `RangeIndex`. The rows keep the order of the file. `cutoff`
    is the percentile of the scores of the links that the file holds with the score threshold or
    more.
    """

    links: pd.DataFrame
    cutoff: float


def read_string_links(
    path: str | os.PathLike[str], *, score_threshold: int = 400, percentile: float = 0.95
) -> StringLinks:
    """Read the STRING protein links, and keep the links with a score above a percentile.

    This is the selection of `get_stringPPI_to_stringGGI` of SeedMatchR, which reads the links
    through STRINGdb 2.22.0 (`R/rstring.R:1064-1067`). The path is a local path: a URL raises
    `InputError`, because the core modules do not use the network. The file is UTF-8 text, or
    gzip, bzip2, or xz of UTF-8 text, and the first bytes, not the suffix, show the compression.
    It is the file that `seedmatch-fetch string` downloads. Line 1 is the header
    `protein1 protein2 combined_score`. Each other line has 3 fields with one space between them:
    two protein IDs and the combined score. A blank line is skipped.

    The reader keeps the links with `combined_score >= score_threshold`, and it takes the
    `percentile` of their scores with `r_quantile`, the quantile of R. Each kept row counts, thus
    a symmetric file counts each pair two times, as R does. The result keeps the links with a
    score above this cutoff, in file order. `score_threshold` is an int from 0 to 1000. The
    default 400 is the default of STRINGdb 2.22.0. `percentile` is a number from 0 to 1. The
    default 0.95 is the default `p` of `get_stringPPI_to_stringGGI`.

    The score is an integer from 0 to 1000. A score with an integer value, for example `400.0`
    or `4e2`, is valid (`D-059`). The protein IDs stay as they are: a quote character and `#` have
    no special meaning (`D-059`). These inputs raise `InputError`:

    - a first line that is not the header
    - a line with another number of fields (`D-059`)
    - a score that is not an integer from 0 to 1000 (`D-059`)
    - an empty protein ID
    - a line of more than 1024 characters with its line end (`D-059`)
    - a line with a NUL character (`D-059`)
    - text that is not valid UTF-8
    - a file without a link, with no link of `score_threshold` or more, or with no link above the
      cutoff. R stops with `non-character argument` for the last input.

    R skips a byte order mark at the start of the file, and so does the reader. The line number
    in a message counts the blank lines.

    The reader reads 2^16 lines at a time, and the time is linear in the number of lines. The work
    memory is one chunk, 10 bytes for each kept row, and the distinct protein IDs. The 10 bytes
    are two 4-byte codes and a 2-byte score. The cutoff needs 8 more bytes for each kept row for
    a moment. The result holds 24 bytes for each link and the distinct strings.

    The STRING rat file of `seedmatch-fetch string 10116` has 14,100,659 lines and 761,435,565
    bytes of text. The longest line has 54 characters with its line end. With the defaults, the
    read takes 8.6 seconds, and the peak memory above the start is about 95 MB. The result has
    98,920 links, the cutoff is 940.0, and the links have 8,809 proteins. With `score_threshold=0`
    and `percentile=0`, the result has 13,894,314 links, the read takes 10.1 seconds, and the peak
    is about 660 MB. Machine: Apple M2 Max, Python 3.12.
    """
    source = os.fspath(path)
    file = local_path(source)
    threshold = _as_threshold(score_threshold)
    probability = _as_percentile(percentile)
    rows = _Rows(threshold)
    try:
        with open_text(file) as handle:
            for first, chunk in read_chunks(handle, _CHUNK_LINES):
                rows.add(chunk, first, source)
    except UnicodeDecodeError as exc:
        raise not_utf8(source, exc, "STRING links") from exc
    return rows.result(source, probability)


def _as_threshold(value: object) -> int:
    threshold = as_limit("score_threshold", value)
    if not 0 <= threshold <= STRING_SCORE_MAX:
        message = (
            f"score_threshold {threshold} is out of range; "
            f"expected an integer from 0 to {STRING_SCORE_MAX}"
        )
        raise InputError(message)
    return threshold


def _as_percentile(value: object) -> float:
    percentile = as_cutoff("percentile", value)
    if not 0.0 <= percentile <= 1.0:
        message = f"percentile {percentile} is out of range; expected a number from 0 to 1"
        raise InputError(message)
    return percentile


@dataclass
class _Rows:
    """The kept links of the chunks read so far.

    A link holds the codes of its two proteins in one vocabulary, and its score. Thus the chunks
    keep one string for each distinct protein ID, not one for each link.
    """

    threshold: int
    proteins: dict[str, int] = field(default_factory=dict)
    first_codes: list[NDArray[np.int32]] = field(default_factory=list)
    second_codes: list[NDArray[np.int32]] = field(default_factory=list)
    scores: list[NDArray[np.int16]] = field(default_factory=list)
    has_header: bool = False
    data_lines: int = 0

    def add(self, chunk: list[str], first: int, source: str) -> None:
        """Parse the data lines of the chunk, and keep the links with a score of the threshold.

        `first` is the number of the first line. A chunk that starts at line 1 starts with the
        header.
        """
        if first == 1:
            _check_header(chunk[0], source)
            self.has_header = True
        is_data = field_lines(
            chunk,
            first,
            source,
            separator=" ",
            fields=len(_COLUMNS),
            described="protein1, protein2, and combined_score",
        )
        if first == 1:
            is_data[0] = False
        numbers = first + np.flatnonzero(is_data)
        if not numbers.size:
            return
        frame, scores = _read_chunk(chunk, is_data, first, source, numbers)
        self.data_lines += len(frame)
        keep = scores >= self.threshold
        if not keep.any():
            return
        frame = frame[keep]
        self.first_codes.append(vocabulary_codes(frame["protein1"], self.proteins))
        self.second_codes.append(vocabulary_codes(frame["protein2"], self.proteins))
        self.scores.append(scores[keep].astype(np.int16))

    def result(self, source: str, percentile: float) -> StringLinks:
        """Return the links above the cutoff. A file without such a link raises `InputError`."""
        if not self.has_header:
            message = f"{source!r}: the file has no line; expected the header {STRING_HEADER!r}"
            raise InputError(message)
        if not self.data_lines:
            message = f"{source!r}: the file has no line after the header; expected STRING links"
            raise InputError(message)
        if not self.scores:
            message = (
                f"{source!r}: no link has a score of {self.threshold} or more; "
                "expected at least one link"
            )
            raise InputError(message)
        scores = concatenate_and_clear(self.scores)
        cutoff = r_quantile(scores, percentile)
        is_above = scores > cutoff
        if not is_above.any():
            message = (
                f"{source!r}: no link has a score above the cutoff {cutoff}, "
                f"the percentile {percentile} of the {len(scores)} scores of "
                f"{self.threshold} or more; expected at least one link"
            )
            raise InputError(message)
        texts = np.array(list(self.proteins), dtype=object)
        links = frame_of(
            {
                "protein1": texts[concatenate_and_clear(self.first_codes)[is_above]],
                "protein2": texts[concatenate_and_clear(self.second_codes)[is_above]],
                "combined_score": scores[is_above].astype(np.int64),
            }
        )
        return StringLinks(links=links, cutoff=cutoff)


def _check_header(line: str, source: str) -> None:
    """Raise `InputError` if line 1 is not the header."""
    text = line.rstrip("\r\n")
    if text != STRING_HEADER:
        message = f"{source!r}: line 1 is {text!r}; expected the header {STRING_HEADER!r}"
        raise InputError(message)


def _read_chunk(
    chunk: list[str],
    is_data: NDArray[np.bool_],
    first: int,
    source: str,
    numbers: NDArray[np.int64],
) -> tuple[pd.DataFrame, NDArray[np.int64]]:
    """Return the table of the data lines and their scores. The first bad line raises `InputError`.

    The read with an integer score is the fast path. If pandas rejects a score, or a score is
    outside the range, the function reads the chunk again with a text score to find the first bad
    line. A score such as `4e2` has an integer value, and it is valid.
    """
    try:
        # pandas casts an infinite score to int64, which warns, before it raises ValueError.
        with np.errstate(invalid="ignore"):
            frame = parse_data_lines(
                chunk,
                is_data,
                first,
                source,
                separator=" ",
                names=_COLUMNS,
                dtype={"protein1": str, "protein2": str, "combined_score": np.dtype(np.int64)},
            )
    except InputError:
        raise
    except (ValueError, OverflowError):
        pass
    else:
        scores = frame["combined_score"].to_numpy()
        if scores.min() >= 0 and scores.max() <= STRING_SCORE_MAX:
            _check_rows(frame, np.zeros(len(frame), dtype=np.bool_), numbers, source)
            return frame, scores
    frame = parse_data_lines(chunk, is_data, first, source, separator=" ", names=_COLUMNS)
    values = pd.to_numeric(frame["combined_score"], errors="coerce").to_numpy(dtype=np.float64)
    is_score = (values == np.floor(values)) & (values >= 0) & (values <= STRING_SCORE_MAX)
    _check_rows(frame, ~is_score, numbers, source)
    return frame, values.astype(np.int64)


def _check_rows(
    frame: pd.DataFrame, is_bad_score: NDArray[np.bool_], numbers: NDArray[np.int64], source: str
) -> None:
    """Raise `InputError` for the first row with an empty protein ID or a bad score."""
    first_empty = frame["protein1"].eq("").to_numpy()
    second_empty = frame["protein2"].eq("").to_numpy()
    bad = first_empty | second_empty | is_bad_score
    if not bad.any():
        return
    row = int(np.argmax(bad))
    if first_empty[row]:
        problem = "protein1 is empty; expected a STRING protein ID"
    elif second_empty[row]:
        problem = "protein2 is empty; expected a STRING protein ID"
    else:
        problem = (
            f"score {frame['combined_score'].iloc[row]!r} is not an integer from 0 to "
            f"{STRING_SCORE_MAX}; expected a STRING combined score"
        )
    message = f"{source!r}: line {numbers[row]}: {problem}"
    raise InputError(message)
