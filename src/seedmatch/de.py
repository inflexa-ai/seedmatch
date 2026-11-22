"""The table of DESeq2 results: the reader for the GEO files and the filter."""

import csv
import os
from collections.abc import Iterator, Sequence
from itertools import batched, repeat
from typing import Final, TextIO

import pandas as pd

from seedmatch._checks import as_cutoff, check_columns, check_numeric_columns
from seedmatch._frames import with_text_dtype
from seedmatch._io import local_path, open_text
from seedmatch.errors import InputError

DE_REQUIRED_COLUMNS: Final = ("gene_id", "baseMean", "log2FoldChange", "pvalue", "padj")

_DESEQ2_COLUMNS: Final = (
    "gene_id",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
    "symbol",
)
_NUMERIC_COLUMNS: Final = ("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
_TEXT_COLUMNS: Final = ("gene_id", "symbol")
_FILTER_COLUMNS: Final = ("baseMean", "log2FoldChange", "pvalue", "padj")
_BATCH_ROWS: Final = 2**15


class _Lines:
    """The lines of a text handle. `ended` is true after the last line."""

    def __init__(self, handle: TextIO) -> None:
        self._handle = handle
        self.ended = False

    def __iter__(self) -> Iterator[str]:
        yield from self._handle
        self.ended = True


def _data_rows(handle: TextIO, shown: str) -> Iterator[list[str]]:
    """Yield the fields of each data row, with an empty field for each field after the last.

    The header is the first line that is not blank. After the header, a line `""` is also no row.
    """
    width = len(_DESEQ2_COLUMNS)
    lines = _Lines(handle)
    # The default dialect reads a quote as read.csv does. In strict mode, the reader also rejects
    # a field such as `"ab"c`, which read.csv and pandas read as `abc`.
    reader = csv.reader(lines, delimiter="\t")
    header = False
    start = 1
    try:
        for row in reader:
            line, start = start, reader.line_num + 1
            # Each row ends at the end of a line. Only a quoted field without its closing quote
            # makes the reader read past the last line, and it then gives the rest of the file as
            # the field, with no error.
            if lines.ended:
                message = (
                    f"{shown!r}: line {line}: the quoted field has no closing quote; "
                    "expected a closing quote before the end of the file"
                )
                raise InputError(message)
            if not row or (header and row == [""]):
                continue
            if not header:
                if len(row) != width:
                    message = (
                        f"{shown!r}: table has {len(row)} columns; "
                        f"expected {width}: {', '.join(_DESEQ2_COLUMNS)}"
                    )
                    raise InputError(message)
                header = True
            elif len(row) > width:
                message = (
                    f"{shown!r}: line {line} has {len(row)} fields; "
                    f"expected at most {width}, as in the header"
                )
                raise InputError(message)
            else:
                row.extend(repeat("", width - len(row)))
                yield row
    except csv.Error as exc:
        message = f"{shown!r}: line {reader.line_num}: {exc}"
        raise InputError(message) from exc
    if not header:
        message = f"{shown!r}: the file has no header line; expected a header line and rows"
        raise InputError(message)


def _table(rows: Sequence[list[str]], shown: str) -> pd.DataFrame:
    """Return the rows as a table with the DESeq2 columns and their types."""
    table = pd.DataFrame(rows, columns=pd.Index(_DESEQ2_COLUMNS), dtype=object)
    for name in _TEXT_COLUMNS:
        table[name] = table[name].where(table[name] != "NA")
    for name in _NUMERIC_COLUMNS:
        field = table[name]
        try:
            # read.csv of R also reads an empty field as missing, but only in a number column.
            table[name] = field.where(~field.isin(("NA", ""))).astype("float64")
        except ValueError as exc:
            message = f"{shown!r}: column {name!r}: {exc}; expected a number, NA, or an empty field"
            raise InputError(message) from exc
    return with_text_dtype(table, _TEXT_COLUMNS)


def read_deseq2_table(path: str | os.PathLike[str]) -> pd.DataFrame:
    """Read a tab-separated DESeq2 table with a header line, and name its columns.

    The path is a local path: a URL raises `InputError`, because the core modules do not use the
    network. The file is UTF-8 text, or gzip, bzip2, or xz of UTF-8 text. As in `read.csv` of R, the
    first bytes of the file, not the suffix, show the compression. The names replace the header, as
    `download_parse_file` of R does. `gene_id` and `symbol` are text, and the other columns are
    `float64`. As in `read.csv` of R, only `NA` is missing in a text column, and `NA` or an empty
    field is missing in a number column.

    The reader obeys the rules of `read.csv` for the lines. A quoted field can hold a tab or a line
    end. A blank line is no row, and after the header, a line `""` is also no row. A line of spaces
    is a row. A row with fewer than 8 fields gets an empty field for each field after its last
    field. A quoted field without its closing quote, a row with more than 8 fields, and a field of
    more than `csv.field_size_limit()` characters raise `InputError`.
    """
    shown = os.fspath(path)
    file = local_path(shown)
    try:
        with open_text(file) as handle:
            # The text of a parsed row takes more memory than its row in the table, thus a batch
            # bounds the memory of the text.
            parts = [
                _table(rows, shown) for rows in batched(_data_rows(handle, shown), _BATCH_ROWS)
            ]
    except UnicodeDecodeError as exc:
        message = (
            f"{shown!r}: the file is not valid UTF-8 text ({exc.reason}); "
            "expected a text or gzip DESeq2 table"
        )
        raise InputError(message) from exc
    if not parts:
        return _table((), shown)
    return pd.concat(parts, ignore_index=True)


def filter_de(
    de: pd.DataFrame,
    *,
    fdr_cutoff: float = 1.0,
    fc_cutoff: float = 0.0,
    basemean_cutoff: float = 30.0,
) -> pd.DataFrame:
    """Return the rows of the DE table that `filter_res` of R keeps, in the same order.

    A row needs a `padj`, a `pvalue`, and a `log2FoldChange`. It also needs `padj` at most
    `fdr_cutoff`, `abs(log2FoldChange)` at least `fc_cutoff`, and `baseMean` at least
    `basemean_cutoff`. The result has a new index. A cutoff is an int or a float: a bool, text,
    `None`, and NaN raise `InputError` (`D-025`).
    """
    fdr = as_cutoff("fdr_cutoff", fdr_cutoff)
    fc = as_cutoff("fc_cutoff", fc_cutoff)
    basemean = as_cutoff("basemean_cutoff", basemean_cutoff)
    check_columns(de, DE_REQUIRED_COLUMNS)
    check_numeric_columns(de, _FILTER_COLUMNS)
    # A missing baseMean fails its comparison, as in R, so the row drops without a separate test.
    keep = (
        de[["padj", "pvalue", "log2FoldChange"]].notna().all(axis=1)
        & (de["padj"] <= fdr)
        & (de["log2FoldChange"].abs() >= fc)
        & (de["baseMean"] >= basemean)
    )
    return de.loc[keep].reset_index(drop=True)
