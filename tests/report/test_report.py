import dataclasses
import math
import re

import numpy as np
import pandas as pd
import pytest
from helpers import (
    GUIDE,
    NOT_A_DATA_FRAME,
    TEXT_DTYPE,
    assert_frame_matches_golden,
    sirna_de,
    sirna_seqs,
)

import seedmatch
from seedmatch import Guide, SequenceSet
from seedmatch.errors import InputError
from seedmatch.report import grid
from seedmatch.report.grid import FullSearchOptions
from seedmatch.report.report import (
    REPORT_CATEGORIES,
    ReportOptions,
    SeedMatchReport,
    seed_match_report,
)

# The order of expected.names in R/SeedMatchReport.R:1-5.
R_CATEGORIES = (
    *(f"Full-D{d}" for d in range(5)),
    *(f"18mer-D{d}" for d in range(5)),
    *(f"15mer-D{d}" for d in (4, 3, 2, 1, 0)),
    "8mer",
    "7mer-m8",
    "7mer-A1",
    "6mer",
)
ROW_NAMES = [
    "In silico predictions",
    "Expressed predictions",
    "Off-target predictions",
    "% off-target",
]
# The 8mer site of the guide, and a subject that holds only the 6mer site of the guide.
SITE_8MER = "CTCTATAA"
SITE_6MER = "GTCTATAG"
FOLD_CHANGE = -2.0
PADJ = 0.01


def de_of(*rows):
    """Return a DE table from (gene_id, log2FoldChange, padj, baseMean) rows."""
    return pd.DataFrame(rows, columns=["gene_id", "log2FoldChange", "padj", "baseMean"])


def empty_de():
    """Return a DE table without rows and with the dtypes of a table that has rows."""
    return de_of().astype(
        {"gene_id": object, "log2FoldChange": "float64", "padj": "float64", "baseMean": "float64"}
    )


def seqs_of(*records, gene_ids=None):
    return SequenceSet.from_records(records, gene_ids=gene_ids)


def report_row(got, name):
    return got.report.set_index("category").loc[name]


@pytest.mark.parametrize(
    ("variant", "options"),
    [("default", ReportOptions()), ("wobble", ReportOptions(wobble=True, fixed_pattern=False))],
)
def test_report_matches_r(variant, options):
    got = seed_match_report(sirna_de(), sirna_seqs(), Guide(GUIDE), options)
    assert_frame_matches_golden(got.report, f"seedmatchr/report_{variant}.tsv")
    assert_frame_matches_golden(got.matches, f"seedmatchr/report_matches_{variant}.tsv")


def test_report_known_snapshot_values():
    row = report_row(
        seed_match_report(sirna_de(), sirna_seqs(), Guide(GUIDE)), "In silico predictions"
    )
    assert (row["Total"], row["15mer-D4"], row["8mer"], row["6mer"]) == (125, 73, 1, 4)


def test_report_counts_genes_not_sequences_D003():
    seqs = seqs_of(("t1", SITE_8MER), ("t2", SITE_8MER), gene_ids=["g", "g"])
    de = de_of(("g", FOLD_CHANGE, PADJ, 100.0))
    got = seed_match_report(de, seqs, Guide(GUIDE))
    for name in ROW_NAMES[:3]:
        assert report_row(got, name)["8mer"] == 1
        assert report_row(got, name)["Total"] == 1
    assert len(got.matches) == 2


@pytest.mark.parametrize("order", [(0, 1), (1, 0)])
def test_the_best_match_of_a_gene_is_the_best_over_its_sequences_D003(order):
    records = [("t1", SITE_6MER), ("t2", SITE_8MER)]
    seqs = seqs_of(*(records[i] for i in order), gene_ids=["g", "g"])
    got = seed_match_report(de_of(("g", FOLD_CHANGE, PADJ, 1.0)), seqs, Guide(GUIDE))
    row = report_row(got, "In silico predictions")
    assert (row["8mer"], row["6mer"], row["Total"]) == (1, 0, 1)
    assert sorted(got.matches["best_match"]) == ["6mer", "8mer"]


def test_percent_is_zero_when_no_gene_is_expressed():
    de = de_of(("other", 0.0, 1.0, 100.0))
    pct = report_row(seed_match_report(de, sirna_seqs(), Guide(GUIDE)), "% off-target")
    assert (pct[list(REPORT_CATEGORIES)] == 0).all()
    assert pct["Total"] == 0


def test_report_categories_are_the_r_names_in_the_r_order():
    assert R_CATEGORIES == REPORT_CATEGORIES


def test_report_categories_are_the_labels_of_the_best_match():
    # The ranking is the one source of the labels, so the test reads it and spells out no copy.
    labels = {label for _, label in grid._BEST_MATCH_RANKING}  # noqa: SLF001
    assert set(REPORT_CATEGORIES) == labels


def test_report_columns_rows_and_dtypes():
    got = seed_match_report(sirna_de(), sirna_seqs(), Guide(GUIDE), group="si-Ttr")
    assert got.report.columns.tolist() == ["category", *REPORT_CATEGORIES, "Total", "group"]
    assert got.report["category"].tolist() == ROW_NAMES
    assert got.report["group"].unique().tolist() == ["si-Ttr"]
    assert got.report.index.tolist() == [0, 1, 2, 3]
    numbers = got.report[[*REPORT_CATEGORIES, "Total"]]
    assert (numbers.dtypes == np.float64).all()


def test_the_total_column_is_the_sum_of_the_categories_except_in_the_percent_row():
    got = seed_match_report(sirna_de(), sirna_seqs(), Guide(GUIDE)).report.set_index("category")
    counts = got.loc[ROW_NAMES[:3], list(REPORT_CATEGORIES)]
    assert counts.sum(axis=1).tolist() == got.loc[ROW_NAMES[:3], "Total"].tolist()
    total = got["Total"]
    assert (
        total["% off-target"]
        == total["Off-target predictions"] / total["Expressed predictions"] * 100
    )


def test_percent_is_the_share_of_off_targets_among_the_expressed_genes():
    seqs = seqs_of(("a", SITE_8MER), ("b", SITE_8MER), ("c", SITE_8MER))
    de = de_of(("a", FOLD_CHANGE, PADJ, 1.0), ("b", 0.0, PADJ, 1.0), ("d", 0.0, PADJ, 1.0))
    got = seed_match_report(de, seqs, Guide(GUIDE))
    assert report_row(got, "In silico predictions")["8mer"] == 3
    assert report_row(got, "Expressed predictions")["8mer"] == 2
    assert report_row(got, "Off-target predictions")["8mer"] == 1
    assert report_row(got, "% off-target")["8mer"] == 50
    assert report_row(got, "% off-target")["Total"] == 50


@pytest.mark.parametrize(
    ("options", "search"),
    [
        (
            ReportOptions(),
            FullSearchOptions(
                seed_edit_distances=(0,), indels=True, fixed_pattern=True, fixed_subject=True
            ),
        ),
        (
            ReportOptions(indels=False, fixed_pattern=False, fixed_subject=False),
            FullSearchOptions(seed_edit_distances=(0,)),
        ),
        (
            ReportOptions(wobble=True, fixed_pattern=False, seed_edit_distances=(0, 1)),
            FullSearchOptions(wobble=True, indels=True, fixed_pattern=False, fixed_subject=True),
        ),
    ],
)
def test_matches_is_the_full_search_table_then_the_de_columns(options, search):
    sequences = sirna_seqs()
    got = seed_match_report(sirna_de(), sequences, Guide(GUIDE), options).matches
    searched = seedmatch.full_search(Guide(GUIDE), sequences, search, "SeedMatchReport")
    de_columns = ["detected", "log2FoldChange", "padj", "baseMean"]
    assert got.columns.tolist() == [*searched.columns, *de_columns]
    assert got[searched.columns.tolist()].equals(searched)
    assert got["detected"].dtype == bool
    assert (got[de_columns[1:]].dtypes == np.float64).all()


def test_a_gene_without_a_de_row_is_undetected_with_missing_values():
    seqs = seqs_of(("a", SITE_8MER), ("b", SITE_8MER))
    got = seed_match_report(de_of(("a", FOLD_CHANGE, PADJ, 5.0)), seqs, Guide(GUIDE))
    assert got.matches["detected"].tolist() == [True, False]
    assert got.matches.loc[0, ["log2FoldChange", "padj", "baseMean"]].tolist() == [
        FOLD_CHANGE,
        PADJ,
        5.0,
    ]
    assert got.matches.loc[1, ["log2FoldChange", "padj", "baseMean"]].isna().all()
    assert report_row(got, "Expressed predictions")["8mer"] == 1


def test_a_gene_in_the_de_table_without_a_sequence_changes_nothing():
    seqs = seqs_of(("a", SITE_8MER))
    de = de_of(("a", FOLD_CHANGE, PADJ, 5.0), ("z", FOLD_CHANGE, PADJ, 5.0))
    got = seed_match_report(de, seqs, Guide(GUIDE))
    assert len(got.matches) == 1
    assert report_row(got, "Expressed predictions")["Total"] == 1


@pytest.mark.parametrize(
    ("fold_change", "padj"),
    [(math.nan, PADJ), (FOLD_CHANGE, math.nan), (math.nan, math.nan)],
)
def test_a_missing_fold_change_or_padj_is_not_an_off_target(fold_change, padj):
    de = de_of(("a", fold_change, padj, 5.0))
    got = seed_match_report(de, seqs_of(("a", SITE_8MER)), Guide(GUIDE))
    assert report_row(got, "Expressed predictions")["8mer"] == 1
    assert report_row(got, "Off-target predictions")["Total"] == 0


@pytest.mark.parametrize(
    ("fold_change", "padj", "off"),
    [
        (-1.0, 0.05, 1),
        (-0.999, 0.05, 0),
        (-1.0, 0.0500001, 0),
        (-1.5, 0.01, 1),
        (1.5, 0.01, 0),
    ],
)
def test_the_cutoffs_are_inclusive_and_the_fold_change_must_be_a_decrease(fold_change, padj, off):
    de = de_of(("a", fold_change, padj, 5.0))
    options = ReportOptions(log2fc_cut=-1.0, padj_cut=0.05)
    got = seed_match_report(de, seqs_of(("a", SITE_8MER)), Guide(GUIDE), options)
    assert report_row(got, "Off-target predictions")["8mer"] == off


def test_an_integer_cutoff_is_stored_as_a_float():
    options = ReportOptions(log2fc_cut=-1, padj_cut=1)
    assert (options.log2fc_cut, options.padj_cut) == (-1.0, 1.0)
    assert type(options.log2fc_cut) is float
    assert type(options.padj_cut) is float


def test_a_numpy_number_is_a_valid_cutoff_and_is_stored_as_a_float():
    options = ReportOptions(log2fc_cut=np.float32(-1), padj_cut=np.int64(1))
    assert (options.log2fc_cut, options.padj_cut) == (-1.0, 1.0)
    assert type(options.log2fc_cut) is float
    assert type(options.padj_cut) is float


def test_a_gene_with_no_match_is_in_no_row():
    got = seed_match_report(
        de_of(("a", FOLD_CHANGE, PADJ, 5.0)), seqs_of(("a", "AAAAAAAA")), Guide(GUIDE)
    )
    assert got.matches["best_match"].tolist() == ["No.matches"]
    assert (got.report[[*REPORT_CATEGORIES, "Total"]].iloc[:3] == 0).all().all()


def test_an_empty_sequence_set_and_a_typed_empty_de_table_give_zero_counts():
    got = seed_match_report(empty_de(), seqs_of(), Guide(GUIDE))
    assert got.matches.empty
    text_columns = got.matches[["seq_id", "gene_id", "group", "best_match"]]
    assert [str(dtype) for dtype in text_columns.dtypes] == [TEXT_DTYPE] * 4
    assert (got.report[[*REPORT_CATEGORIES, "Total"]] == 0).all().all()


def test_a_typed_empty_de_table_with_sequences_leaves_every_gene_undetected():
    got = seed_match_report(empty_de(), seqs_of(("a", SITE_8MER)), Guide(GUIDE))
    assert got.matches["detected"].tolist() == [False]
    assert report_row(got, "In silico predictions")["8mer"] == 1
    assert report_row(got, "Expressed predictions")["Total"] == 0


def test_an_empty_sequence_set_gives_zero_counts():
    got = seed_match_report(de_of(("a", FOLD_CHANGE, PADJ, 5.0)), seqs_of(), Guide(GUIDE))
    assert got.matches.empty
    assert (got.report[[*REPORT_CATEGORIES, "Total"]] == 0).all().all()
    assert got.report["category"].tolist() == ROW_NAMES


def test_a_duplicate_gene_in_the_de_table_uses_its_first_row():
    seqs = seqs_of(("a", SITE_8MER), ("b", SITE_8MER))
    de = de_of(("a", 0.0, 0.5, 1.0), ("b", FOLD_CHANGE, PADJ, 2.0), ("a", FOLD_CHANGE, PADJ, 3.0))
    got = seed_match_report(de, seqs, Guide(GUIDE))
    assert got.matches["seq_id"].tolist() == ["a", "b"]
    assert got.matches["baseMean"].tolist() == [1.0, 2.0]
    assert report_row(got, "Off-target predictions")["8mer"] == 1


def test_a_duplicate_gene_in_the_de_table_does_not_copy_a_sequence_row():
    seqs = seqs_of(("t1", SITE_8MER), ("t2", SITE_6MER), gene_ids=["a", "a"])
    de = de_of(("a", FOLD_CHANGE, PADJ, 1.0), ("a", FOLD_CHANGE, PADJ, 2.0), ("a", 0.0, 1.0, 3.0))
    got = seed_match_report(de, seqs, Guide(GUIDE))
    assert got.matches["seq_id"].tolist() == ["t1", "t2"]
    assert report_row(got, "In silico predictions")["Total"] == 1


@pytest.mark.parametrize("dtype", ["object", "string", "category"])
def test_gene_ids_of_a_text_dtype_match(dtype):
    de = de_of(("a", FOLD_CHANGE, PADJ, 5.0)).astype({"gene_id": dtype})
    got = seed_match_report(de, seqs_of(("a", SITE_8MER)), Guide(GUIDE))
    assert got.matches["detected"].tolist() == [True]
    assert str(got.matches["gene_id"].dtype) == TEXT_DTYPE


def test_the_de_columns_of_integers_become_float64():
    de = pd.DataFrame({"gene_id": ["a"], "log2FoldChange": [-2], "padj": [0], "baseMean": [7]})
    got = seed_match_report(de, seqs_of(("a", SITE_8MER)), Guide(GUIDE))
    assert (got.matches[["log2FoldChange", "padj", "baseMean"]].dtypes == np.float64).all()
    assert report_row(got, "Off-target predictions")["8mer"] == 1


def test_the_report_does_not_change_the_de_table():
    de = sirna_de()
    before = de.copy()
    seed_match_report(de, sirna_seqs(), Guide(GUIDE))
    pd.testing.assert_frame_equal(de, before)


def test_workers_does_not_change_the_result():
    alone = seed_match_report(sirna_de(), sirna_seqs(), Guide(GUIDE), workers=1)
    pooled = seed_match_report(sirna_de(), sirna_seqs(), Guide(GUIDE), workers=4)
    pd.testing.assert_frame_equal(alone.report, pooled.report)
    pd.testing.assert_frame_equal(alone.matches, pooled.matches)


def test_the_wobble_adds_in_silico_predictions_D026():
    plain = seed_match_report(sirna_de(), sirna_seqs(), Guide(GUIDE))
    wobble = seed_match_report(
        sirna_de(), sirna_seqs(), Guide(GUIDE), ReportOptions(wobble=True, fixed_pattern=False)
    )
    total = "In silico predictions"
    assert report_row(wobble, total)["Total"] > report_row(plain, total)["Total"]


@pytest.mark.parametrize("make", [ReportOptions, FullSearchOptions])
def test_the_wobble_with_a_fixed_pattern_causes_an_error_D026(make):
    message = "wobble=True needs fixed_pattern=False"
    with pytest.raises(InputError, match=re.escape(message)):
        make(wobble=True, fixed_pattern=True)


def test_the_report_is_a_frozen_dataclass_of_two_frames():
    got = seed_match_report(sirna_de(), sirna_seqs(), Guide(GUIDE))
    assert isinstance(got, SeedMatchReport)
    assert dataclasses.is_dataclass(got)
    with pytest.raises(dataclasses.FrozenInstanceError):
        got.report = got.matches


def test_the_results_of_two_calls_are_not_equal_and_each_result_has_a_hash():
    first = seed_match_report(sirna_de(), sirna_seqs(), Guide(GUIDE))
    second = seed_match_report(sirna_de(), sirna_seqs(), Guide(GUIDE))
    assert first != second
    assert first == first  # noqa: PLR0124  a result equals itself, because it compares by identity
    assert hash(first) != hash(second)
    assert len({first, second}) == 2


def test_the_names_are_exported():
    assert seedmatch.seed_match_report is seed_match_report
    assert seedmatch.SeedMatchReport is SeedMatchReport
    assert seedmatch.ReportOptions is ReportOptions
    assert seedmatch.REPORT_CATEGORIES is REPORT_CATEGORIES
    assert seedmatch.report.seed_match_report is seed_match_report


def test_report_options_defaults_copy_r():
    options = ReportOptions()
    assert dataclasses.astuple(options) == (True, False, -0.7655347, 0.05, (0,), True, True)


def test_report_options_are_frozen_and_normalized():
    options = ReportOptions(seed_edit_distances=[0, 1])
    assert options.seed_edit_distances == (0, 1)
    with pytest.raises(dataclasses.FrozenInstanceError):
        options.wobble = True


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"log2fc_cut": "-1"}, "log2fc_cut '-1' has type str"),
        ({"padj_cut": None}, "padj_cut None has type NoneType"),
        ({"seed_edit_distances": ()}, "seed_edit_distances is empty"),
        ({"seed_edit_distances": (-1,)}, "seed_edit_distances -1 is negative"),
        ({"seed_edit_distances": (0, 0)}, "seed_edit_distances repeats 0"),
        ({"seed_edit_distances": {0, 1}}, "seed_edit_distances .* has type set"),
        ({"indels": 1}, "indels 1 has type int; expected a bool"),
        ({"wobble": "yes"}, "wobble 'yes' has type str; expected a bool"),
        ({"fixed_pattern": 0}, "fixed_pattern 0 has type int; expected a bool"),
        ({"fixed_subject": None}, "fixed_subject None has type NoneType; expected a bool"),
    ],
)
def test_report_options_reject_a_bad_value(kwargs, message):
    with pytest.raises(InputError, match=message):
        ReportOptions(**kwargs)


@pytest.mark.parametrize("name", ["log2fc_cut", "padj_cut"])
def test_a_nan_cutoff_is_an_error_D024(name):
    with pytest.raises(InputError, match=f"{name} nan is not a number"):
        ReportOptions(**{name: math.nan})


@pytest.mark.parametrize(
    ("name", "value", "kind"),
    [("log2fc_cut", True, "bool"), ("padj_cut", False, "bool"), ("log2fc_cut", "-1", "str")],
)
def test_a_cutoff_that_is_not_an_int_or_a_float_is_an_error_D024(name, value, kind):
    message = f"{name} {value!r} has type {kind}; expected an int or a float"
    with pytest.raises(InputError, match=message):
        ReportOptions(**{name: value})


def test_a_cutoff_that_a_float_cannot_hold_is_an_error_D024():
    with pytest.raises(InputError, match=r"padj_cut 1000.* is too large; expected a number"):
        ReportOptions(padj_cut=10**400)


def test_a_missing_column_is_an_error_that_lists_the_missing_columns_in_sorted_order():
    de = pd.DataFrame({"gene_id": ["a"], "padj": [0.01]})
    with pytest.raises(InputError, match=r"missing columns: \['baseMean', 'log2FoldChange'\]"):
        seed_match_report(de, seqs_of(("a", SITE_8MER)), Guide(GUIDE))


def test_a_gene_id_column_that_is_not_text_is_an_error_that_names_the_dtype():
    de = de_of((1, FOLD_CHANGE, PADJ, 5.0))
    with pytest.raises(InputError, match=r"'gene_id' has dtype int64; expected text"):
        seed_match_report(de, seqs_of(("1", SITE_8MER)), Guide(GUIDE))


@pytest.mark.parametrize("column", ["log2FoldChange", "padj", "baseMean"])
def test_a_column_that_is_not_numeric_is_an_error_that_names_the_column_and_the_dtype_D024(column):
    de = de_of(("a", FOLD_CHANGE, PADJ, 5.0)).astype({column: "object"})
    de[column] = "x"
    message = f"column '{column}' has dtype {TEXT_DTYPE}; expected a numeric"
    with pytest.raises(InputError, match=message):
        seed_match_report(de, seqs_of(("a", SITE_8MER)), Guide(GUIDE))


@pytest.mark.parametrize("column", ["log2FoldChange", "padj", "baseMean"])
@pytest.mark.parametrize("dtype", ["bool", "complex128"])
def test_a_bool_and_a_complex_column_are_not_numeric_D024(column, dtype):
    de = de_of(("a", FOLD_CHANGE, PADJ, 5.0)).astype({column: dtype})
    message = f"column '{column}' has dtype {dtype}; expected a numeric dtype"
    with pytest.raises(InputError, match=f"^{message}$"):
        seed_match_report(de, seqs_of(("a", SITE_8MER)), Guide(GUIDE))


@pytest.mark.parametrize(("de", "type_name"), NOT_A_DATA_FRAME)
def test_a_table_that_is_not_a_data_frame_is_an_error(de, type_name):
    message = f"de has type {type_name}; expected a pandas.DataFrame"
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        seed_match_report(de, seqs_of(("a", SITE_8MER)), Guide(GUIDE))


def test_a_guide_that_is_too_short_for_a_seed_is_an_error_before_any_search():
    with pytest.raises(InputError, match="18mer"):
        seed_match_report(sirna_de(), sirna_seqs(), Guide("UUAUAGAGCAAGAACACU"))
