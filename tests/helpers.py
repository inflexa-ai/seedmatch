import bz2
import gzip
import lzma
import math
import zlib
from functools import cache
from pathlib import Path
from typing import Any

import pandas as pd

from seedmatch import (
    FeatureType,
    Guide,
    Reference,
    Seed,
    alignment_text,
    annotate_seed_counts,
    count_matches,
    e_value_null,
    filter_de,
    find_matches,
    full_search_ranges,
    kmer_counts,
    match_distances,
    seed_match_report,
)
from seedmatch.sequence import SequenceSet, read_fasta

DATA = Path(__file__).parent / "data"
GOLDEN = Path(__file__).parent / "golden"
GUIDE = "UUAUAGAGCAAGAACACUGUUUU"


COMMENT_LINES = "".join(f"# comment {i}\n" for i in range(5000)).encode()


def damaged_compressed_files() -> dict[str, tuple[bytes, type[Exception]]]:
    """Map the name of a damaged file to its bytes and to the exception of its decompressor.

    The file holds 5000 comment lines, which a GTF reader skips, so a reader reaches the damage.
    """
    gzipped, bzipped, xzipped = (
        pack(COMMENT_LINES) for pack in (gzip.compress, bz2.compress, lzma.compress)
    )
    return {
        "truncated gzip": (gzipped[: len(gzipped) // 2], EOFError),
        "gzip with a bad header": (b"\x1f\x8b" + bytes(30), gzip.BadGzipFile),
        "gzip with a bad block": (gzipped[:10] + b"\xff" * 40, zlib.error),
        "truncated bzip2": (bzipped[: len(bzipped) // 2], EOFError),
        "bzip2 with a bad block": (b"BZh9" + b"\xff" * 40, OSError),
        "truncated xz": (xzipped[: len(xzipped) // 2], EOFError),
        "xz with a bad block": (b"\xfd7zXZ\x00" + b"\xff" * 40, lzma.LZMAError),
    }


def golden(path: str) -> pd.DataFrame:
    return pd.read_csv(GOLDEN / path, sep="\t", comment="#", dtype=str, keep_default_na=False)


# The field names of each tuple are the columns of the golden file, so no static type can name them.
def golden_rows(path: str) -> list[Any]:
    return list(golden(path).itertuples(index=False))


def golden_row(path: str, case: str, alternative: str) -> Any:
    (row,) = (
        row for row in golden_rows(path) if row.case == case and row.alternative == alternative
    )
    return row


def ecdf_golden_row(input_name: str, test: str, alternative: str) -> Any:
    (row,) = (
        row
        for row in golden_rows("stats/ecdf.tsv")
        if (row.input, row.test, row.alternative) == (input_name, test, alternative)
    )
    return row


# R finds D of the KS test with a cumulative sum of 1/m and -1/n, and it takes the weights of the
# DTS statistic from floating-point divisions. Thus each statistic differs from its exact form by
# rounding errors. A p-value has no such noise.
STATISTIC_TOLERANCE = {"rel_tol": 1e-9, "abs_tol": 1e-12}


def dts_p_value_error(p_r: float, n_boot: int) -> float:
    """Return the largest difference that is plausible between two DTS p-values of n_boot draws.

    Each p-value is a fraction of n draws, so the difference has the variance 2 p (1 - p) / n.
    """
    return 4 * math.sqrt(2 * p_r * (1 - p_r) / n_boot) + 2 / n_boot


# The cache shares one dict of lists, so a test must not change them.
@cache
def stats_samples() -> dict[str, tuple[list[float], list[float]]]:
    cases: dict[str, tuple[list[float], list[float]]] = {}
    for case, sample, value in golden("stats/samples.tsv").itertuples(index=False):
        x, y = cases.setdefault(case, ([], []))
        (x if sample == "x" else y).append(float(value))
    return cases


# A SequenceSet is immutable, so each test can share the one set.
@cache
def sirna_seqs() -> SequenceSet:
    return read_fasta(DATA / "Test_seqs.fa")


# The file has a mer7m8 column, so this reader keeps it. read_deseq2_table would reject the file.
# A DataFrame is mutable, so each call reads the file again.
def sirna_de() -> pd.DataFrame:
    return pd.read_csv(DATA / "DESEQ2_results.tsv", sep="\t")


def run_core_pipeline() -> None:
    seqs = read_fasta(DATA / "Test_seqs.fa")
    de = filter_de(sirna_de(), basemean_cutoff=0)
    annotate_seed_counts(de, seqs, Guide(GUIDE).target(Seed.MER7M8))
    full_search_ranges(Guide(GUIDE), seqs)
    kmer_counts(seqs, k=3)
    alignment_text(find_matches(Guide(GUIDE).target(Seed.MER7M8), seqs).iloc[0], Guide(GUIDE), seqs)
    match_distances(find_matches(Guide(GUIDE).target(Seed.MER6), seqs))
    seed_match_report(de, seqs, Guide(GUIDE))
    e_value_null(Guide(GUIDE), seqs, rng=0, n_random=2)
    annotation = DATA / "annotation"
    with Reference.from_paths(
        ensdb=annotation / "synthetic.sqlite", genome=annotation / "synthetic.2bit"
    ) as ref:
        utr3 = ref.sequences(FeatureType.UTR3)
    count_matches(Guide(GUIDE).target(Seed.MER7M8), utr3)


# The columns that an R golden file names differently from the Python frames. The column `names` of
# a vmatchPattern2 result repeats `seqnames`, so the comparison drops it.
R_COLUMNS = {
    "seqs": "seq_id",
    "seqnames": "seq_id",
    "Category": "category",
    "Group": "group",
    "best.match": "best_match",
    "TargetSeq": "target_seq",
    "EditDistance": "edit_distance",
    "AllowIndels": "indels",
    "IN.DESEQ2": "detected",
}
_R_BOOLEANS = {"TRUE": True, "FALSE": False, "Detected": True, "Undetected": False}
_R_MISSING = "NA"
_FLOAT_RELATIVE_TOLERANCE = 1e-9


def _kind(column: pd.Series) -> str:
    if pd.api.types.is_bool_dtype(column):
        return "bool"
    if pd.api.types.is_integer_dtype(column):
        return "int"
    if pd.api.types.is_float_dtype(column):
        return "float"
    return "str"


_READERS = {"bool": _R_BOOLEANS.__getitem__, "int": int, "float": float, "str": str}


# The kind comes from the dtype of the frame, so a text column that holds digits stays text.
def _equal(value: object, text: str, kind: str) -> bool:
    missing = bool(pd.isna(value))
    if text == _R_MISSING or missing:
        return text == _R_MISSING and missing
    try:
        expected = _READERS[kind](text)
    except (KeyError, ValueError):
        return False
    if kind == "float":
        return math.isclose(float(value), expected, rel_tol=_FLOAT_RELATIVE_TOLERANCE)
    return value == expected


def assert_frame_matches_golden(got: pd.DataFrame, path: str) -> None:
    """Fail at the first row and column where `got` differs from the golden file.

    The comparison reads only the columns of the golden file, under their Python names, so `got` can
    have more columns. An integer, a boolean, and a string must be equal. A float must agree within
    a relative tolerance of 1e-9.
    """
    expected = golden(path).drop(columns=["names"], errors="ignore").rename(columns=R_COLUMNS)
    if len(got) != len(expected):
        message = f"{path}: the frame has {len(got)} rows, expected {len(expected)}"
        raise AssertionError(message)
    for name in expected.columns:
        if name not in got.columns:
            columns = list(got.columns)
            message = f"{path}: the frame has no column {name!r}; its columns are {columns}"
            raise AssertionError(message)
    first: tuple[int, int, str, object, str, str] | None = None
    for position, name in enumerate(expected.columns):
        kind = _kind(got[name])
        for row, (value, text) in enumerate(zip(got[name], expected[name], strict=True)):
            if not _equal(value, text, kind):
                if first is None or (row, position) < first[:2]:
                    first = (row, position, name, value, text, kind)
                break
    if first is not None:
        row, _, name, value, text, kind = first
        shown = repr(text) if kind == "str" else text
        message = f"{path}: row {row} (0-based) column {name!r}: got {value!r}, expected {shown}"
        raise AssertionError(message)


# Values that are not a DataFrame, and the type name that a message gives for each one.
NOT_A_DATA_FRAME = [
    ({"gene_id": ["g1"]}, "dict"),
    ([1, 2], "list"),
    (None, "NoneType"),
    (pd.Series(["g1"]), "Series"),
]
