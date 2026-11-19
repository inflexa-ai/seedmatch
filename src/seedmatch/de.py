"""The table of DESeq2 results: the reader for the GEO files and the filter."""

import os
import warnings
from collections.abc import Sequence
from typing import Final

import pandas as pd

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


def read_deseq2_table(path: str | os.PathLike[str]) -> pd.DataFrame:
    """Read a tab-separated DESeq2 table with a header line, and name its columns.

    The file can be gzip. The names replace the header, as `download_parse_file` of R does.
    `gene_id` and `symbol` are text, and the other columns are `float64`. As in `read.csv` of R,
    only `NA` is missing in a text column, and `NA` or an empty field is missing in a number
    column.
    """
    shown = os.fspath(path)
    try:
        with warnings.catch_warnings():
            # With index_col=False, pandas drops the extra fields of a row and only warns.
            warnings.simplefilter("error", pd.errors.ParserWarning)
            table = pd.read_csv(
                path,
                sep="\t",
                index_col=False,
                dtype=str,
                na_values=["NA"],
                keep_default_na=False,
            )
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


def check_columns(de: pd.DataFrame, required: Sequence[str]) -> None:
    """Raise `InputError` that lists the missing columns in sorted order."""
    missing = sorted(set(required) - set(de.columns))
    if missing:
        message = f"DE table has missing columns: {missing}; expected {list(required)}"
        raise InputError(message)


def check_numeric_columns(de: pd.DataFrame, names: Sequence[str]) -> None:
    """Raise `InputError` that names the first column whose dtype is not numeric."""
    for name in names:
        if not pd.api.types.is_numeric_dtype(de[name]):
            message = f"column {name!r} has dtype {de[name].dtype}; expected a numeric dtype"
            raise InputError(message)


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
    `basemean_cutoff`. The result has a new index.
    """
    check_columns(de, DE_REQUIRED_COLUMNS)
    check_numeric_columns(de, _FILTER_COLUMNS)
    # A missing baseMean fails its comparison, as in R, so the row drops without a separate test.
    keep = (
        de[["padj", "pvalue", "log2FoldChange"]].notna().all(axis=1)
        & (de["padj"] <= fdr_cutoff)
        & (de["log2FoldChange"].abs() >= fc_cutoff)
        & (de["baseMean"] >= basemean_cutoff)
    )
    return de.loc[keep].reset_index(drop=True)
