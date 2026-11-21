import bz2
import gzip
import io
import lzma
import re
import warnings

import pandas as pd
import pytest
from helpers import DATA, NOT_A_DATA_FRAME, damaged_compressed_files, golden_rows, sirna_de

import seedmatch
from seedmatch import de
from seedmatch.de import DE_REQUIRED_COLUMNS, filter_de, read_deseq2_table
from seedmatch.errors import InputError

DESEQ2_COLUMNS = [
    "gene_id",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
    "symbol",
]
ROWS = [
    ["ENSRNOG1", 2138.5, -8.25, 0.5, -23.5, 2.5e-123, 3.5e-119, "Ttr"],
    ["ENSRNOG2", 437.5, -1.5, 0.25, -12.5, 2.0e-36, 1.5e-32, "Kpna6"],
]


def read_cases():
    return pd.read_csv(DATA / "de_filter_cases.tsv", sep="\t")


def golden_filter_ids(cutoffs):
    return [
        row.gene_id
        for row in golden_rows("de/filter_res.tsv")
        if (float(row.fdr_cutoff), float(row.fc_cutoff), float(row.basemean_cutoff)) == cutoffs
    ]


def write_table(path, header):
    text = "\n".join("\t".join(map(str, row)) for row in [header, *ROWS]) + "\n"
    if path.suffix == ".gz":
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)


@pytest.mark.parametrize("cutoffs", [(1, 0, 30), (0.05, 0.5, 10)])
def test_filter_de_matches_r(cutoffs):
    expected = golden_filter_ids(cutoffs)
    assert expected
    got = filter_de(
        read_cases(), fdr_cutoff=cutoffs[0], fc_cutoff=cutoffs[1], basemean_cutoff=cutoffs[2]
    )
    assert got["gene_id"].tolist() == expected


def test_filter_de_defaults_are_the_defaults_of_filter_res():
    assert filter_de(read_cases())["gene_id"].tolist() == golden_filter_ids((1, 0, 30))


def test_filter_de_keeps_the_input_and_resets_the_index():
    cases = read_cases()
    before = cases.copy()
    got = filter_de(cases, fdr_cutoff=0.05, fc_cutoff=0.5, basemean_cutoff=10)
    pd.testing.assert_frame_equal(cases, before)
    assert got.index.tolist() == list(range(len(got)))
    assert list(got.columns) == list(cases.columns)


def test_filter_de_removes_rows_with_a_missing_value():
    got = filter_de(read_cases(), basemean_cutoff=0.0)
    assert not got[["padj", "pvalue", "log2FoldChange", "baseMean"]].isna().to_numpy().any()
    assert not {"g02_na_padj", "g03_na_pvalue", "g04_na_log2fc", "g05_na_basemean"} & set(
        got["gene_id"]
    )


def test_filter_de_names_the_missing_columns():
    with pytest.raises(InputError, match=r"missing columns: \['padj', 'pvalue'\]"):
        filter_de(pd.DataFrame({"gene_id": [], "baseMean": [], "log2FoldChange": []}))


def test_filter_de_rejects_a_repeated_padj_column_D051():
    """R stops too, with `Can't transform a data frame with duplicate names.` (dplyr 1.2.1)."""
    cases = read_cases()
    de = pd.concat([cases, cases[["padj"]]], axis=1)
    message = "DE table has 2 columns 'padj'; expected one column for each name"
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        filter_de(de)


def test_required_columns_are_the_five_columns_of_filter_res():
    assert DE_REQUIRED_COLUMNS == ("gene_id", "baseMean", "log2FoldChange", "pvalue", "padj")


@pytest.mark.parametrize("suffix", [".tsv", ".tsv.gz"])
def test_read_deseq2_table_names_the_eight_columns(tmp_path, suffix):
    path = tmp_path / f"table{suffix}"
    write_table(path, ["Geneid", "mean", "lfc", "se", "wald", "p", "fdr", "name"])
    got = read_deseq2_table(path)
    assert got.columns.tolist() == DESEQ2_COLUMNS
    assert got["gene_id"].tolist() == ["ENSRNOG1", "ENSRNOG2"]
    assert got["padj"].tolist() == [3.5e-119, 1.5e-32]
    assert got["symbol"].tolist() == ["Ttr", "Kpna6"]


def test_read_deseq2_table_accepts_a_str_path(tmp_path):
    path = tmp_path / "table.tsv"
    write_table(path, DESEQ2_COLUMNS)
    assert len(read_deseq2_table(str(path))) == len(ROWS)


def test_read_deseq2_table_rejects_a_table_with_more_than_8_columns_D021():
    with pytest.raises(InputError, match=r"table has 9 columns; expected 8"):
        read_deseq2_table(DATA / "DESEQ2_results.tsv")


def test_read_deseq2_table_names_the_path_of_a_table_with_too_few_columns(tmp_path):
    path = tmp_path / "short.tsv"
    path.write_text("a\tb\tc\n1\t2\t3\n")
    with pytest.raises(InputError, match=r"short\.tsv.*table has 3 columns; expected 8"):
        read_deseq2_table(path)


def test_read_deseq2_table_result_passes_the_filter(tmp_path):
    path = tmp_path / "table.tsv.gz"
    write_table(path, DESEQ2_COLUMNS)
    assert filter_de(read_deseq2_table(path), basemean_cutoff=1000.0)["gene_id"].tolist() == [
        "ENSRNOG1"
    ]


def test_filter_de_accepts_the_test_data():
    assert len(filter_de(sirna_de())) > 0


def test_top_level_package_exports_the_de_functions():
    assert (seedmatch.read_deseq2_table, seedmatch.filter_de) == (read_deseq2_table, filter_de)


def write_text(tmp_path, text):
    path = tmp_path / "table.tsv"
    path.write_text(text)
    return path


def test_read_deseq2_table_reads_missing_values_as_r_does(tmp_path):
    rows = [
        ["g1", "10", "1.5", "0.1", "2.5", "0.01", "", ""],
        ["g2", "NA", "NaN", "Inf", "-Inf", "0.02", "0.5", "null"],
        ["NA", "5", "1.5", "0.1", "2.5", "0.03", "0.6", "None"],
        ["g4", "5", "1.5", "0.1", "2.5", "0.03", "0.6", "NA"],
    ]
    text = "\n".join("\t".join(row) for row in [DESEQ2_COLUMNS, *rows]) + "\n"
    got = read_deseq2_table(write_text(tmp_path, text))
    assert got["symbol"].tolist()[:3] == ["", "null", "None"]
    assert got["symbol"].isna().tolist() == [False, False, False, True]
    assert got["gene_id"].isna().tolist() == [False, False, True, False]
    assert got["padj"].isna().tolist() == [True, False, False, False]
    assert got["baseMean"].isna().tolist() == [False, True, False, False]
    assert got["log2FoldChange"].isna().tolist() == [False, True, False, False]
    assert (got["lfcSE"][1], got["stat"][1]) == (float("inf"), float("-inf"))
    assert {str(got[name].dtype) for name in DESEQ2_COLUMNS[1:7]} == {"float64"}


def test_read_deseq2_table_reads_gene_id_and_symbol_as_text(tmp_path):
    rows = [["0012", "1", "1", "1", "1", "1", "1", "5"], ["34", "1", "1", "1", "1", "1", "1", "6"]]
    text = "\n".join("\t".join(row) for row in [DESEQ2_COLUMNS, *rows]) + "\n"
    got = read_deseq2_table(write_text(tmp_path, text))
    assert got["gene_id"].tolist() == ["0012", "34"]
    assert got["symbol"].tolist() == ["5", "6"]


def test_read_deseq2_table_names_the_column_with_a_value_that_is_not_a_number_D021(tmp_path):
    row = ["g1", "abc", "1", "1", "1", "1", "1", "x"]
    text = "\n".join("\t".join(r) for r in [DESEQ2_COLUMNS, row]) + "\n"
    with pytest.raises(InputError, match=r"column 'baseMean'.*'abc'"):
        read_deseq2_table(write_text(tmp_path, text))


def test_read_deseq2_table_reads_a_table_without_rows(tmp_path):
    got = read_deseq2_table(write_text(tmp_path, "\t".join(DESEQ2_COLUMNS) + "\n"))
    assert got.empty
    assert got.columns.tolist() == DESEQ2_COLUMNS


def test_read_deseq2_table_rejects_an_empty_file(tmp_path):
    with pytest.raises(InputError, match=r"table\.tsv.*no header line"):
        read_deseq2_table(write_text(tmp_path, ""))


@pytest.mark.parametrize("extra_rows", [1, 2])
def test_read_deseq2_table_rejects_a_row_with_more_fields_than_the_header_D021(
    tmp_path, extra_rows
):
    row = ["g1", "1", "1", "1", "1", "1", "1", "x"]
    rows = [[*row, "y"], *[row] * (extra_rows - 1)]
    text = "\n".join("\t".join(r) for r in [DESEQ2_COLUMNS, *rows]) + "\n"
    with pytest.raises(InputError, match=r"table\.tsv.*Expected 8 fields in line 2, saw 9"):
        read_deseq2_table(write_text(tmp_path, text))


class FiltersResetOnRead(io.StringIO):
    """A handle whose read resets the warning filters, as another thread can do at that time.

    `warnings.catch_warnings` of Python 3.12 restores the filters of the process at its exit.
    """

    def read(self, size=-1):
        warnings.resetwarnings()
        return super().read(size)


def test_read_deseq2_table_rejects_a_row_with_more_fields_when_another_thread_resets_the_filters(
    tmp_path, monkeypatch
):
    row = ["g1", "1", "1", "1", "1", "1", "1", "x"]
    text = "\n".join("\t".join(r) for r in [DESEQ2_COLUMNS, [*row, "y"], row]) + "\n"
    path = write_text(tmp_path, text)
    monkeypatch.setattr(de, "open_text", lambda _file: FiltersResetOnRead(text))
    with (
        warnings.catch_warnings(),
        pytest.raises(InputError, match=r"table\.tsv.*Expected 8 fields in line 2, saw 9"),
    ):
        read_deseq2_table(path)


def test_read_deseq2_table_rejects_a_later_row_with_more_fields_than_the_header(tmp_path):
    row = ["g1", "1", "1", "1", "1", "1", "1", "x"]
    text = "\n".join("\t".join(r) for r in [DESEQ2_COLUMNS, row, [*row, "y"]]) + "\n"
    with pytest.raises(InputError, match=r"table\.tsv.*Expected 8 fields in line 3, saw 9"):
        read_deseq2_table(write_text(tmp_path, text))


@pytest.mark.parametrize("column", ["baseMean", "log2FoldChange", "pvalue", "padj"])
def test_filter_de_names_a_filter_column_that_is_not_numeric(column):
    cases = read_cases()
    cases[column] = cases[column].astype(str)
    with pytest.raises(InputError, match=rf"column '{column}' has dtype (object|str)"):
        filter_de(cases)


@pytest.mark.parametrize("column", ["baseMean", "log2FoldChange", "pvalue", "padj"])
@pytest.mark.parametrize("dtype", ["bool", "complex128"])
def test_filter_de_rejects_a_bool_and_a_complex_column_D024(column, dtype):
    cases = read_cases()
    cases[column] = cases[column].astype(dtype)
    message = f"column '{column}' has dtype {dtype}; expected a numeric dtype"
    with pytest.raises(InputError, match=f"^{message}$"):
        filter_de(cases)


@pytest.mark.parametrize(("de", "type_name"), NOT_A_DATA_FRAME)
def test_filter_de_rejects_a_table_that_is_not_a_data_frame(de, type_name):
    message = f"de has type {type_name}; expected a pandas.DataFrame"
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        filter_de(de)


def test_read_deseq2_table_refuses_a_url_and_opens_no_connection(monkeypatch):
    url = "https://example.org/de.tsv"

    def refuse(*_args, **_kwargs):
        message = "urlopen was called"
        raise AssertionError(message)

    monkeypatch.setattr("urllib.request.urlopen", refuse)
    message = f"{url!r}: the path is a URL; expected a local path"
    with pytest.raises(InputError, match=re.escape(message)):
        read_deseq2_table(url)


def test_read_deseq2_table_finds_gzip_from_the_content_not_the_suffix(tmp_path):
    packed = tmp_path / "table.gz"
    write_table(packed, DESEQ2_COLUMNS)
    renamed = packed.rename(tmp_path / "table.tsv")
    assert read_deseq2_table(renamed)["gene_id"].tolist() == ["ENSRNOG1", "ENSRNOG2"]


def test_read_deseq2_table_rejects_text_that_is_not_utf8(tmp_path):
    path = tmp_path / "table.tsv"
    row = [*ROWS[0][:-1], "Caf\xe9"]
    text = "\t".join(DESEQ2_COLUMNS) + "\n" + "\t".join(map(str, row)) + "\n"
    path.write_bytes(text.encode("latin-1"))
    with pytest.raises(InputError, match=r"table\.tsv': the file is not valid UTF-8 text"):
        read_deseq2_table(path)


@pytest.mark.parametrize(
    ("cutoff", "message"),
    [
        (float("nan"), "fdr_cutoff nan is not a number; expected an int or a float, not NaN"),
        (None, "fdr_cutoff None has type NoneType; expected an int or a float"),
        ("0.05", "fdr_cutoff '0.05' has type str; expected an int or a float"),
    ],
)
def test_filter_de_rejects_a_cutoff_that_is_not_a_number_D025(cutoff, message):
    with pytest.raises(InputError, match=re.escape(message)):
        filter_de(sirna_de(), fdr_cutoff=cutoff)


@pytest.mark.parametrize("name", ["fdr_cutoff", "fc_cutoff", "basemean_cutoff"])
def test_filter_de_checks_each_cutoff_D025(name):
    with pytest.raises(InputError, match=f"^{name} nan is not a number"):
        filter_de(sirna_de(), **{name: float("nan")})


@pytest.mark.parametrize("compress", [bz2.compress, lzma.compress], ids=["bzip2", "xz"])
def test_read_deseq2_table_reads_bzip2_and_xz_from_the_content_as_r_does(tmp_path, compress):
    plain = tmp_path / "plain.tsv"
    write_table(plain, DESEQ2_COLUMNS)
    packed = tmp_path / "table.tsv"
    packed.write_bytes(compress(plain.read_bytes()))
    assert read_deseq2_table(packed)["gene_id"].tolist() == ["ENSRNOG1", "ENSRNOG2"]


@pytest.mark.parametrize("name", damaged_compressed_files())
def test_read_deseq2_table_names_the_file_of_a_damaged_compressed_table(tmp_path, name):
    path = tmp_path / "table.tsv"
    path.write_bytes(damaged_compressed_files()[name][0])
    prefix = f"{str(path)!r}: the compressed data is damaged or truncated ("
    with pytest.raises(InputError, match=re.escape(prefix)) as raised:
        read_deseq2_table(path)
    assert str(raised.value).endswith("); expected a complete gzip, bzip2, or xz file")


def test_read_deseq2_table_names_the_path_of_a_directory(tmp_path):
    message = f"{str(tmp_path)!r}: the path is a directory; expected a file"
    with pytest.raises(InputError, match=re.escape(message)):
        read_deseq2_table(tmp_path)
