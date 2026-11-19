"""The table of DESeq2 results: the reader for the GEO files and the filter."""

import bz2
import gzip
import io
import lzma
import os
import warnings
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Final, TextIO

import pandas as pd

from seedmatch._checks import as_cutoff, check_columns, check_numeric_columns
from seedmatch._io import local_path
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
_FILTER_COLUMNS: Final = ("baseMean", "log2FoldChange", "pvalue", "padj")
_GZIP_MAGIC: Final = b"\x1f\x8b"
_BZIP2_MAGIC: Final = b"BZh"
_XZ_MAGIC: Final = b"\xfd7zXZ\x00"


def _decompressed(
    raw: io.BufferedReader,
) -> gzip.GzipFile | bz2.BZ2File | lzma.LZMAFile | io.BufferedReader:
    """Return a reader of the bytes after the compression that the first bytes show, if any."""
    head = raw.peek(len(_XZ_MAGIC))
    if head.startswith(_GZIP_MAGIC):
        return gzip.GzipFile(fileobj=raw)
    if head.startswith(_BZIP2_MAGIC):
        return bz2.BZ2File(raw)
    if head.startswith(_XZ_MAGIC):
        return lzma.LZMAFile(raw)
    return raw


@contextmanager
def _open_text(path: Path) -> Iterator[TextIO]:
    """Open the file as UTF-8 text. The first bytes, not the suffix, show a compressed file."""
    with path.open("rb") as raw:
        binary = _decompressed(raw)
        # pandas also opens a path with newline="", so its parser gets each line end unchanged.
        with io.TextIOWrapper(binary, encoding="utf-8", newline="") as text:
            yield text


def read_deseq2_table(path: str | os.PathLike[str]) -> pd.DataFrame:
    """Read a tab-separated DESeq2 table with a header line, and name its columns.

    The path is a local path: a URL raises `InputError`, because the core modules do not use the
    network. The file is UTF-8 text, or gzip, bzip2, or xz of UTF-8 text. As in `read.csv` of R, the
    first bytes of the file, not the suffix, show the compression. The names replace the header, as
    `download_parse_file` of R does. `gene_id` and `symbol` are text, and the other columns are
    `float64`. As in `read.csv` of R, only `NA` is missing in a text column, and `NA` or an empty
    field is missing in a number column.
    """
    shown = os.fspath(path)
    file = local_path(shown)
    try:
        # pandas gets an open handle and no path, so it cannot open a URL or guess a compression.
        with _open_text(file) as handle, warnings.catch_warnings():
            # With index_col=False, pandas drops the extra fields of a row and only warns.
            warnings.simplefilter("error", pd.errors.ParserWarning)
            table = pd.read_csv(
                handle,
                sep="\t",
                index_col=False,
                dtype=str,
                na_values=["NA"],
                keep_default_na=False,
                compression=None,
            )
    except UnicodeDecodeError as exc:
        message = (
            f"{shown!r}: the file is not valid UTF-8 text ({exc.reason}); "
            "expected a text or gzip DESeq2 table"
        )
        raise InputError(message) from exc
    except pd.errors.EmptyDataError as exc:
        message = f"{shown!r}: the file has no header line; expected a header line and rows"
        raise InputError(message) from exc
    except pd.errors.ParserWarning as exc:
        message = f"{shown!r}: a row has more fields than the header; expected at most as many"
        raise InputError(message) from exc
    except pd.errors.ParserError as exc:
        message = f"{shown!r}: {str(exc).strip()}"
        raise InputError(message) from exc
    if table.shape[1] != len(_DESEQ2_COLUMNS):
        message = (
            f"{shown!r}: table has {table.shape[1]} columns; "
            f"expected {len(_DESEQ2_COLUMNS)}: {', '.join(_DESEQ2_COLUMNS)}"
        )
        raise InputError(message)
    table.columns = pd.Index(_DESEQ2_COLUMNS)
    for name in _NUMERIC_COLUMNS:
        field = table[name]
        try:
            # read.csv of R also reads an empty field as missing, but only in a number column.
            table[name] = field.where(field != "").astype("float64")
        except ValueError as exc:
            message = f"{shown!r}: column {name!r}: {exc}; expected a number, NA, or an empty field"
            raise InputError(message) from exc
    return table


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
