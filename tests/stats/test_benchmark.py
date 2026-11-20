import dataclasses
import math

import numpy as np
import pandas as pd
import pytest
from helpers import (
    GUIDE,
    NOT_A_DATA_FRAME,
    STATISTIC_TOLERANCE,
    assert_curve_matches_golden,
    exactly,
    golden,
    golden_note,
    sirna_de,
)

import seedmatch
import seedmatch.stats
from seedmatch import Guide, Seed, SequenceSet, annotate_seed_counts
from seedmatch.errors import InputError
from seedmatch.stats import PrCurve, RocCurve, pr_curve, roc_curve
from seedmatch.stats.benchmark import Benchmark, benchmark_predictions

BENCHMARK = golden("analyses/benchmark.tsv")
BENCHMARK_FILE = "analyses/benchmark.tsv"
UNPATCHED = r"the unpatched benchmark_preds\(res, padj = 1e-10\) "


def small_de(**columns):
    """Return 6 genes. The gold standard of the default cutoffs is g1 and g2."""
    table = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3", "g4", "g5", "g6"],
            "padj": [0.01, 0.1, 0.01, 0.01, 0.5, 0.9],
            "log2FoldChange": [-1.0, -2.0, 0.5, 1.0, -1.0, -1.0],
            "mer7m8": [3, 2, 1, 0, 1, 0],
        }
    )
    for name, values in columns.items():
        table[name] = values
    return table


def case_de(case):
    """Return the table of a case of the golden file, as the R script makes it."""
    de = sirna_de()
    if case == "benchmark_na":
        passing = de.loc[(de["padj"] <= 1e-10) & (de["log2FoldChange"] < 0), "gene_id"]
        de.loc[de["gene_id"].isin(passing.head(2)), "padj"] = np.nan
    return de


def close(value, reference):
    return math.isclose(value, float(reference), **STATISTIC_TOLERANCE)


@pytest.mark.parametrize("case", BENCHMARK["case"])
def test_benchmark_matches_r(case):
    row = BENCHMARK[BENCHMARK["case"] == case].iloc[0]
    result = benchmark_predictions(
        case_de(case),
        column=row["column"],
        padj_cutoff=float(row["padj"]),
        l2fc_cutoff=float(row["l2fc"]),
    )
    expected = (int(row["n_positive"]), int(row["n_negative"]))
    assert (result.pr.n_positive, result.pr.n_negative) == expected
    assert (result.roc.n_positive, result.roc.n_negative) == expected
    assert close(result.pr.auc_integral, row["auc_integral"])
    assert close(result.pr.auc_davis_goadrich, row["auc_davis_goadrich"])
    assert close(result.roc.auc, row["roc_auc"])
    assert_curve_matches_golden(
        result.pr.curve, "analyses/pr_curve.tsv", case, ["recall", "precision", "threshold"]
    )
    assert_curve_matches_golden(
        result.roc.curve, "analyses/roc_curve.tsv", case, ["fpr", "tpr", "threshold"]
    )


def test_the_padj_cutoff_applies_D006():
    """D-006: the filter of R keeps each row with a padj. The note of benchmark.tsv gives 63."""
    unpatched = int(golden_note(BENCHMARK_FILE, UNPATCHED + r"keeps (\d+) genes"))
    result = benchmark_predictions(sirna_de(), padj_cutoff=1e-10)
    assert unpatched == 63
    assert result.pr.n_positive == 24


def test_the_roc_curve_uses_the_labels_D006():
    """D-006: R gives the labels as scores.class1. The notes of benchmark.tsv give its ROC area.

    The unpatched call has the 63 labels of benchmark_1, because the PR areas are equal.
    """
    unpatched_pr = float(
        golden_note(BENCHMARK_FILE, UNPATCHED + r"gives the PR auc\.integral (\S+) ")
    )
    unpatched_roc = float(
        golden_note(BENCHMARK_FILE, UNPATCHED + r"gives .* and the ROC auc (\S+)$")
    )
    result = benchmark_predictions(sirna_de())
    row = BENCHMARK[BENCHMARK["case"] == "benchmark_1"].iloc[0]
    assert close(result.pr.auc_integral, unpatched_pr)
    assert close(result.roc.auc, row["roc_auc"])
    assert abs(result.roc.auc - unpatched_roc) > 0.1


def test_the_curves_are_the_curves_of_the_scores_and_the_labels_of_the_gold_standard():
    de = small_de()
    result = benchmark_predictions(de)
    labels = [True, True, False, False, False, False]
    expected_pr = pr_curve(de["mer7m8"], labels)
    expected_roc = roc_curve(de["mer7m8"], labels)
    assert result.pr.auc_integral == expected_pr.auc_integral
    assert result.pr.curve.equals(expected_pr.curve)
    assert result.roc.auc == expected_roc.auc
    assert result.roc.curve.equals(expected_roc.curve)


def test_each_row_of_a_gene_in_the_gold_standard_is_a_positive():
    de = pd.DataFrame(
        {
            "gene_id": ["g1", "g1", "g2", "g3"],
            "padj": [0.01, 0.9, 0.9, 0.01],
            "log2FoldChange": [-1.0, 1.0, -1.0, 1.0],
            "mer7m8": [2, 1, 1, 0],
        }
    )
    result = benchmark_predictions(de)
    assert (result.pr.n_positive, result.pr.n_negative) == (2, 2)
    assert (result.roc.n_positive, result.roc.n_negative) == (2, 2)


def test_a_row_without_padj_is_a_negative():
    de = small_de(padj=[np.nan, 0.1, 0.01, 0.01, 0.5, 0.9])
    result = benchmark_predictions(de)
    assert (result.pr.n_positive, result.pr.n_negative) == (1, 5)


def test_a_row_without_a_fold_change_is_a_negative():
    de = small_de(log2FoldChange=[-1.0, np.nan, 0.5, 1.0, -1.0, -1.0])
    result = benchmark_predictions(de)
    assert (result.pr.n_positive, result.pr.n_negative) == (1, 5)


def test_a_gene_with_a_row_without_padj_and_a_row_that_passes_is_a_positive():
    de = pd.DataFrame(
        {
            "gene_id": ["g1", "g1", "g2"],
            "padj": [np.nan, 0.01, 0.9],
            "log2FoldChange": [-1.0, -1.0, 1.0],
            "mer7m8": [1, 2, 0],
        }
    )
    assert benchmark_predictions(de).pr.n_positive == 2


def test_the_cutoffs_include_padj_and_exclude_the_fold_change_limit():
    de = small_de(
        padj=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1], log2FoldChange=[-1.0, 0.0, 1.0, 0.0, 1.0, 0.0]
    )
    result = benchmark_predictions(de, padj_cutoff=0.1, l2fc_cutoff=0.0)
    assert result.pr.n_positive == 1
    result = benchmark_predictions(de, padj_cutoff=0.1, l2fc_cutoff=0.5)
    assert result.pr.n_positive == 4


def test_the_cutoffs_and_the_column_are_keyword_only():
    with pytest.raises(TypeError):
        benchmark_predictions(small_de(), "mer7m8")


def test_benchmark_takes_the_scores_of_the_column_that_the_caller_names():
    de = small_de(count=[0, 0, 5, 4, 3, 2])
    named = benchmark_predictions(de, column="count")
    assert named.roc.auc == 0.0
    assert benchmark_predictions(de).roc.auc == 1.0


def test_benchmark_accepts_an_infinite_score():
    de = small_de(mer7m8=[np.inf, 2.0, 1.0, 0.0, 1.0, -np.inf])
    assert benchmark_predictions(de).roc.auc == 1.0


def test_benchmark_accepts_the_integer_and_float_cutoffs_and_numpy_scalars():
    de = small_de()
    expected = benchmark_predictions(de, padj_cutoff=1, l2fc_cutoff=0)
    got = benchmark_predictions(de, padj_cutoff=np.float64(1.0), l2fc_cutoff=np.int64(0))
    assert got.pr.n_positive == expected.pr.n_positive == 4


def test_benchmark_gives_a_gene_without_a_sequence_the_score_0():
    """A count column of annotate_seed_counts with shared_genes=False has 0 for such a gene."""
    target = Guide(GUIDE).target(Seed.MER7M8)
    sequences = SequenceSet.from_records(
        [("t1", "AACTCTATAAA"), ("t2", "AACTCTATACTCTATA"), ("t3", "GGGGGGGG")],
        gene_ids=["g1", "g2", "g3"],
    )
    counts = annotate_seed_counts(
        small_de().drop(columns="mer7m8"), sequences, target, shared_genes=False
    )
    assert counts["mer7m8"].tolist() == [1, 2, 0, 0, 0, 0]
    result = benchmark_predictions(counts)
    assert (result.pr.n_positive, result.pr.n_negative) == (2, 4)
    assert result.roc.curve["threshold"].tolist()[-1] == 0.0


def test_benchmark_gives_the_same_result_for_a_repeated_index_and_a_categorical_gene_id():
    expected = benchmark_predictions(small_de())
    de = small_de(gene_id=pd.Categorical(["g1", "g2", "g3", "g4", "g5", "g6"]))
    de.index = [5, 5, 1, 1, 0, 0]
    got = benchmark_predictions(de)
    assert got.pr.curve.equals(expected.pr.curve)
    assert got.roc.curve.equals(expected.roc.curve)
    assert got.pr.n_positive == 2


MISSING_IDS = [
    pytest.param(None, None, id="None-None"),
    pytest.param(np.nan, np.nan, id="NaN-NaN"),
    pytest.param(pd.NA, pd.NA, id="NA-NA"),
    pytest.param(None, np.nan, id="None-NaN"),
    pytest.param(None, pd.NA, id="None-NA"),
    pytest.param(np.nan, pd.NA, id="NaN-NA"),
    pytest.param(pd.NA, None, id="NA-None"),
]


def missing_ids_de(first, second, passing):
    """Return 6 rows. Rows 0 and 1 have a missing gene ID, and only row `passing` passes."""
    ids = pd.Series([first, second, "g3", "g4", "g5", "g6"], dtype=object)
    padj = [0.5] * 6
    l2fc = [1.0] * 6
    padj[passing], l2fc[passing] = 0.01, -1.0
    return small_de(gene_id=ids, padj=padj, log2FoldChange=l2fc)


@pytest.mark.parametrize("passing", [0, 1])
@pytest.mark.parametrize(("first", "second"), MISSING_IDS)
def test_all_missing_gene_ids_are_one_gene_as_na_is_in_r(first, second, passing):
    """In R, NA %in% c(NA) is TRUE, and R has one NA. Thus None, NaN, and pd.NA match each other."""
    result = benchmark_predictions(missing_ids_de(first, second, passing))
    assert (result.pr.n_positive, result.pr.n_negative) == (2, 4)


@pytest.mark.parametrize(("first", "second"), MISSING_IDS)
def test_a_missing_gene_id_is_a_negative_if_no_row_with_a_missing_id_passes(first, second):
    de = missing_ids_de(first, second, passing=2)
    result = benchmark_predictions(de)
    assert (result.pr.n_positive, result.pr.n_negative) == (1, 5)


def test_benchmark_does_not_change_its_input():
    de = small_de(padj=[np.nan, 0.1, 0.01, 0.01, 0.5, 0.9])
    before = de.copy()
    benchmark_predictions(de)
    pd.testing.assert_frame_equal(de, before)


def test_the_result_is_a_frozen_record_of_the_two_curves():
    result = benchmark_predictions(small_de())
    assert isinstance(result, Benchmark)
    assert [field.name for field in dataclasses.fields(result)] == ["pr", "roc"]
    assert isinstance(result.pr, PrCurve)
    assert isinstance(result.roc, RocCurve)
    with pytest.raises(dataclasses.FrozenInstanceError):
        result.pr = result.pr


@pytest.mark.parametrize("name", ["Benchmark", "benchmark_predictions"])
def test_a_benchmark_name_is_exported_from_the_package_and_from_stats(name):
    assert name in seedmatch.__all__
    assert name in seedmatch.stats.__all__
    assert getattr(seedmatch, name) is getattr(seedmatch.stats, name)


def test_benchmark_rejects_a_nan_score_D048():
    """D-048: the patched benchmark_preds stops: NAs are not allowed in subscripted assignments."""
    de = small_de(mer7m8=[3.0, 2.0, 1.0, np.nan, 1.0, 0.0])
    message = "column 'mer7m8' of de holds NaN for gene 'g4'; expected a score for each row"
    with pytest.raises(InputError, match=exactly(message)):
        benchmark_predictions(de)


def test_benchmark_rejects_a_missing_value_of_a_nullable_integer_score_D048():
    de = small_de(mer7m8=pd.array([3, 2, 1, 0, None, 0], dtype="Int64"))
    message = "column 'mer7m8' of de holds NaN for gene 'g5'; expected a score for each row"
    with pytest.raises(InputError, match=exactly(message)):
        benchmark_predictions(de)


@pytest.mark.parametrize("name", ["padj_cutoff", "l2fc_cutoff"])
def test_benchmark_rejects_a_nan_cutoff_D048(name):
    """D-048: R keeps no gene for a NaN cutoff and gives NaN areas (the probe in the entry)."""
    message = f"{name} nan is not a number; expected an int or a float, not NaN"
    with pytest.raises(InputError, match=exactly(message)):
        benchmark_predictions(small_de(), **{name: math.nan})


@pytest.mark.parametrize("name", ["padj_cutoff", "l2fc_cutoff"])
@pytest.mark.parametrize(
    ("value", "type_name"), [(True, "bool"), ("0.15", "str"), (None, "NoneType")]
)
def test_benchmark_rejects_a_cutoff_that_is_not_a_number_D048(name, value, type_name):
    message = f"{name} {value!r} has type {type_name}; expected an int or a float"
    with pytest.raises(InputError, match=exactly(message)):
        benchmark_predictions(small_de(), **{name: value})


@pytest.mark.parametrize(
    ("de", "cutoff", "shown"),
    [
        pytest.param(sirna_de(), 0, "0.0", id="r-note"),
        pytest.param(small_de(), 1e-10, "1e-10", id="small"),
    ],
)
def test_benchmark_rejects_a_gold_standard_without_a_positive_D048(de, cutoff, shown):
    """D-048: for padj = 0 (the case r-note), the notes of benchmark.tsv give the areas NaN."""
    message = (
        f"no gene has padj <= {shown} and log2FoldChange < 0.0; expected at least one positive gene"
    )
    with pytest.raises(InputError, match=exactly(message)):
        benchmark_predictions(de, padj_cutoff=cutoff)


@pytest.mark.parametrize(("l2fc_cutoff", "shown"), [(10, "10.0"), (math.inf, "inf")])
def test_benchmark_rejects_a_gold_standard_without_a_negative_D048(l2fc_cutoff, shown):
    """D-048: for padj = 1 and l2fc = Inf, the notes of benchmark.tsv give the PR areas 1 and 1.

    They give the ROC area NaN. The cutoff 10 keeps the same labels on this data.
    """
    message = (
        f"each row is a positive with padj <= 1.0 and log2FoldChange < {shown}; "
        "expected at least one negative row"
    )
    with pytest.raises(InputError, match=exactly(message)):
        benchmark_predictions(sirna_de(), padj_cutoff=1, l2fc_cutoff=l2fc_cutoff)


def test_benchmark_rejects_a_table_without_rows():
    message = (
        "no gene has padj <= 0.15 and log2FoldChange < 0.0; expected at least one positive gene"
    )
    with pytest.raises(InputError, match=exactly(message)):
        benchmark_predictions(small_de().iloc[:0])


@pytest.mark.parametrize("missing", ["gene_id", "padj", "log2FoldChange", "mer7m8"])
def test_benchmark_rejects_a_missing_column(missing):
    message = (
        f"DE table has missing columns: ['{missing}']; "
        "expected ['gene_id', 'padj', 'log2FoldChange', 'mer7m8']"
    )
    with pytest.raises(InputError, match=exactly(message)):
        benchmark_predictions(small_de().drop(columns=missing))


def test_benchmark_rejects_a_column_that_the_table_lacks():
    message = (
        "DE table has missing columns: ['count']; "
        "expected ['gene_id', 'padj', 'log2FoldChange', 'count']"
    )
    with pytest.raises(InputError, match=exactly(message)):
        benchmark_predictions(small_de(), column="count")


def test_benchmark_rejects_a_repeated_column():
    de = pd.concat([small_de(), small_de()[["padj"]]], axis=1)
    message = "DE table has 2 columns 'padj'; expected one column for each name"
    with pytest.raises(InputError, match=exactly(message)):
        benchmark_predictions(de)


@pytest.mark.parametrize("name", ["padj", "log2FoldChange", "mer7m8"])
def test_benchmark_rejects_a_text_column(name):
    de = small_de()
    de[name] = pd.Series([str(value) for value in de[name]], dtype=object)
    message = f"column '{name}' has dtype object; expected a numeric dtype"
    with pytest.raises(InputError, match=exactly(message)):
        benchmark_predictions(de)


def test_benchmark_rejects_gene_ids_that_are_not_text():
    message = "'gene_id' has dtype int64; expected text values. Convert the column with astype(str)"
    with pytest.raises(InputError, match=exactly(message)):
        benchmark_predictions(small_de(gene_id=[1, 2, 3, 4, 5, 6]))


def test_benchmark_rejects_a_column_name_that_is_not_text():
    message = "column 5 has type int; expected a str"
    with pytest.raises(InputError, match=exactly(message)):
        benchmark_predictions(small_de(), column=5)


@pytest.mark.parametrize(("de", "type_name"), NOT_A_DATA_FRAME)
def test_benchmark_rejects_a_de_that_is_not_a_data_frame(de, type_name):
    message = f"de has type {type_name}; expected a pandas.DataFrame"
    with pytest.raises(InputError, match=exactly(message)):
        benchmark_predictions(de)
