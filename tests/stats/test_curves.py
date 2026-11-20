import dataclasses
import math
import re
from decimal import Decimal, localcontext
from functools import cache

import numpy as np
import pandas as pd
import pytest
from helpers import STATISTIC_TOLERANCE, golden
from hypothesis import given, settings
from hypothesis import strategies as st
from scipy.stats import rankdata

import seedmatch
import seedmatch.stats
from seedmatch.errors import InputError
from seedmatch.stats import PrCurve, RocCurve, pr_curve, roc_curve

CURVES = [pytest.param(pr_curve, id="pr"), pytest.param(roc_curve, id="roc")]
LABEL_VALUES = "0, 1, False, or True"
LABEL_MESSAGE = f"expected {LABEL_VALUES}"


@cache
def curve_samples() -> dict[str, tuple[list[float], list[int]]]:
    cases: dict[str, tuple[list[float], list[int]]] = {}
    for case, score, label in golden("analyses/curve_samples.tsv").itertuples(index=False):
        scores, labels = cases.setdefault(case, ([], []))
        scores.append(float(score))
        labels.append(int(label))
    return cases


CASES = [case for case in curve_samples() if not case.startswith("benchmark_")]


def golden_curve(path: str, case: str) -> pd.DataFrame:
    table = golden(path)
    return table[table["case"] == case].reset_index(drop=True)


def assert_columns_match(got: pd.DataFrame, expected: pd.DataFrame, columns: list[str]) -> None:
    assert len(got) == len(expected)
    for name in columns:
        for row, (value, reference) in enumerate(zip(got[name], expected[name], strict=True)):
            assert math.isclose(value, float(reference), **STATISTIC_TOLERANCE), (name, row)


def golden_row(path: str, case: str) -> pd.Series:
    table = golden(path)
    return table[table["case"] == case].iloc[0]


@pytest.mark.parametrize("case", CASES)
def test_pr_curve_matches_r(case):
    scores, labels = curve_samples()[case]
    result = pr_curve(scores, labels)
    areas = golden_row("analyses/pr.tsv", case)
    assert math.isclose(result.auc_integral, float(areas.auc_integral), **STATISTIC_TOLERANCE)
    assert math.isclose(
        result.auc_davis_goadrich, float(areas.auc_davis_goadrich), **STATISTIC_TOLERANCE
    )
    expected = golden_curve("analyses/pr_curve.tsv", case)
    assert_columns_match(result.curve, expected, ["recall", "precision", "threshold"])
    assert result.curve["threshold"].tolist() == expected["threshold"].astype(float).tolist()
    assert (result.n_positive, result.n_negative) == (sum(labels), len(labels) - sum(labels))


@pytest.mark.parametrize("case", CASES)
def test_roc_curve_matches_r(case):
    scores, labels = curve_samples()[case]
    result = roc_curve(scores, labels)
    assert math.isclose(
        result.auc, float(golden_row("analyses/roc.tsv", case).auc), **STATISTIC_TOLERANCE
    )
    expected = golden_curve("analyses/roc_curve.tsv", case)
    assert_columns_match(result.curve, expected, ["fpr", "tpr", "threshold"])
    assert result.curve["threshold"].tolist() == expected["threshold"].astype(float).tolist()
    assert (result.n_positive, result.n_negative) == (sum(labels), len(labels) - sum(labels))


SCORES = st.one_of(
    st.integers(-3, 3).map(float),
    st.floats(-5, 5),
    st.sampled_from([math.inf, -math.inf]),
)


@st.composite
def samples(draw, *, with_negative):
    """Draw scores with many ties and labels. The first label is a positive."""
    size = draw(st.integers(2, 30))
    scores = draw(st.lists(SCORES, min_size=size, max_size=size))
    labels = [True, *draw(st.lists(st.booleans(), min_size=size - 1, max_size=size - 1))]
    if with_negative:
        labels[1] = False
    return scores, labels


def pairs(scores, labels):
    positives = [s for s, label in zip(scores, labels, strict=True) if label]
    negatives = [s for s, label in zip(scores, labels, strict=True) if not label]
    return positives, negatives


@settings(deadline=None)
@given(samples(with_negative=True))
def test_roc_auc_is_the_probability_that_a_positive_scores_higher(sample):
    scores, labels = sample
    positives, negatives = pairs(scores, labels)
    higher = sum(p > n for p in positives for n in negatives)
    equal = sum(p == n for p in positives for n in negatives)
    expected = (higher + equal / 2) / (len(positives) * len(negatives))
    assert math.isclose(roc_curve(scores, labels).auc, expected, **STATISTIC_TOLERANCE)


def precision_lines(scores, labels):
    """Return (r_prev, r, a, b) of each step that adds a positive, from the definitions."""
    positives = sum(labels)
    tp_prev = fp_prev = 0
    lines = []
    for threshold in sorted(set(scores), reverse=True):
        tp = sum(1 for s, label in zip(scores, labels, strict=True) if s >= threshold and label)
        fp = sum(1 for s, label in zip(scores, labels, strict=True) if s >= threshold and not label)
        if tp > tp_prev:
            h = (fp - fp_prev) / (tp - tp_prev)
            b = (fp_prev - h * tp_prev) / positives
            lines.append((tp_prev / positives, tp / positives, 1 + h, b))
        tp_prev, fp_prev = tp, fp
    return lines


def midpoint_integral(low, high, a, b, points=10**4):
    width = (high - low) / points
    middle = low + (np.arange(points) + 0.5) * width
    return float((middle / (a * middle + b)).sum() * width)


@settings(deadline=None)
@given(samples(with_negative=False))
def test_pr_auc_integral_equals_a_numeric_integral(sample):
    scores, labels = sample
    expected = sum(midpoint_integral(*line) for line in precision_lines(scores, labels))
    assert math.isclose(pr_curve(scores, labels).auc_integral, expected, abs_tol=1e-6)


def exact_auc_integral(scores, labels):
    """Return the area of the formula of `auc_integral` in decimal arithmetic with 50 digits."""
    order = np.argsort(-scores, kind="stable")
    positive = labels[order]
    tp_all = np.cumsum(positive)
    fp_all = np.arange(1, len(scores) + 1) - tp_all
    ends = np.flatnonzero(np.append(scores[order][1:] != scores[order][:-1], True))
    with localcontext() as context:
        context.prec = 50
        positives = Decimal(int(positive.sum()))
        total = Decimal(0)
        tp_prev = fp_prev = 0
        for tp, fp in zip(tp_all[ends].tolist(), fp_all[ends].tolist(), strict=True):
            if tp > tp_prev:
                h = Decimal(fp - fp_prev) / (tp - tp_prev)
                a = 1 + h
                b = (fp_prev - h * tp_prev) / positives
                recall, recall_prev = Decimal(tp) / positives, Decimal(tp_prev) / positives
                if b == 0:
                    total += (recall - recall_prev) / a
                else:
                    logs = (a * recall + b).ln() - (a * recall_prev + b).ln()
                    total += ((recall - recall_prev) - b / a * logs) / a
            tp_prev, fp_prev = tp, fp
    return float(total)


@pytest.mark.parametrize("n_positive", [11, 2044])
def test_auc_integral_keeps_its_precision_for_few_positives_D050(n_positive):
    generator = np.random.default_rng(20261125)
    labels = np.zeros(200_000, dtype=bool)
    labels[generator.choice(len(labels), n_positive, replace=False)] = True
    scores = generator.normal(size=len(labels)) + 0.8 * labels
    got = pr_curve(scores, labels).auc_integral
    assert math.isclose(got, exact_auc_integral(scores, labels), rel_tol=1e-12, abs_tol=0.0)


def test_auc_integral_has_a_small_error_for_one_positive_in_a_million_scores_D050():
    scores = np.random.default_rng(4711).normal(size=1_000_000)
    labels = np.zeros(len(scores), dtype=bool)
    labels[np.argmin(scores)] = True
    got = pr_curve(scores, labels).auc_integral
    # The relative error is 3.2e-10 here. PRROC has 3.7e-4.
    assert math.isclose(got, exact_auc_integral(scores, labels), rel_tol=1e-9, abs_tol=0.0)


@settings(deadline=None)
@given(samples(with_negative=False))
def test_pr_curve_goes_from_recall_1_to_0(sample):
    scores, labels = sample
    recall = pr_curve(scores, labels).curve["recall"].to_numpy()
    assert recall[0] == 1.0
    assert recall[-1] == 0.0
    assert (np.diff(recall) <= 0).all()


@settings(deadline=None)
@given(samples(with_negative=False))
def test_pr_curve_has_at_most_2k_plus_2_plus_max_100_p_rows(sample):
    scores, labels = sample
    bound = 2 * len(set(scores)) + 2 + max(100, sum(labels))
    assert len(pr_curve(scores, labels).curve) <= bound


def test_large_uses_the_recall_step_1_over_p():
    scores, labels = curve_samples()["large"]
    positives = sum(labels)
    assert positives >= 100
    curve = pr_curve(scores, labels).curve.iloc[1:-1]
    steps = [
        np.diff(np.sort(group["recall"].to_numpy()))
        for _, group in curve.groupby("threshold")
        if len(group) > 1
    ]
    differences = np.concatenate(steps)
    near_zero = differences < 1e-12
    near_step = np.isclose(differences, 1 / positives, rtol=1e-9, atol=0)
    assert (near_zero | near_step).all()
    assert near_step.sum() > 100


def assert_rows_stay_inside_their_step(scores, labels):
    curve = pr_curve(scores, labels).curve.iloc[1:-1]
    scores_array, positive = np.asarray(scores, dtype=float), np.asarray(labels, dtype=bool)
    for threshold, rows in curve.groupby("threshold"):
        reached = positive[scores_array >= threshold].sum() / positive.sum()
        before = positive[scores_array > threshold].sum() / positive.sum()
        assert rows["recall"].max() == reached
        assert rows["recall"].min() >= before


@settings(deadline=None)
@given(samples(with_negative=False))
def test_pr_rows_of_a_step_stay_between_its_recalls(sample):
    assert_rows_stay_inside_their_step(*sample)


def test_pr_rows_of_a_step_stay_between_its_recalls_for_alternating_labels():
    # Without the limit, the last inner row of 5 steps is 1 ulp above the recall of its step.
    assert_rows_stay_inside_their_step(np.arange(40.0, 0, -1), np.tile([1, 0], 20))


def test_pr_curve_has_no_inner_rows_if_each_step_adds_one_positive():
    generator = np.random.default_rng(0)
    labels = np.zeros(200, dtype=bool)
    labels[generator.choice(200, 150, replace=False)] = True
    result = pr_curve(np.arange(200.0, 0, -1), labels)
    assert result.n_positive == 150
    assert len(result.curve) == 200 + 2


def test_a_small_example_has_the_areas_of_its_definitions():
    result = pr_curve([0.9, 0.8, 0.7, 0.6], [1, 0, 1, 0])
    assert math.isclose(result.auc_integral, 1 - 0.5 * math.log(1.5), rel_tol=1e-12)
    assert math.isclose(result.auc_davis_goadrich, 19 / 24, rel_tol=1e-12)
    assert roc_curve([0.9, 0.8, 0.7, 0.6], [1, 0, 1, 0]).auc == 0.75


def test_the_roc_threshold_of_a_row_is_the_score_of_the_points_that_the_next_row_adds():
    result = roc_curve([0.9, 0.8, 0.7, 0.6], [1, 0, 1, 0])
    assert result.curve["fpr"].tolist() == [0.0, 0.0, 0.5, 0.5, 1.0]
    assert result.curve["tpr"].tolist() == [0.0, 0.5, 0.5, 1.0, 1.0]
    assert result.curve["threshold"].tolist() == [0.9, 0.8, 0.7, 0.6, 0.6]


def test_roc_auc_equals_the_rank_sum_formula_on_a_large_input():
    generator = np.random.default_rng(3)
    scores = generator.integers(0, 300, size=20_000).astype(float)
    labels = generator.random(20_000) < 1 / (1 + np.exp(-(scores - 150) / 100))
    positives = labels.sum()
    ranks = rankdata(scores)
    expected = (ranks[labels].sum() - positives * (positives + 1) / 2) / (
        positives * (len(labels) - positives)
    )
    assert math.isclose(roc_curve(scores, labels).auc, expected, **STATISTIC_TOLERANCE)


ZERO_LABELS = [1, 0, 1, 0, 1, 0, 0, 1]


# PRROC 1.4 (R 4.5.3) shows the sign of the last zero score in the order of the input.
@pytest.mark.parametrize("curve", CURVES)
@pytest.mark.parametrize(
    ("scores", "labels", "negative"),
    [
        pytest.param([-0.0, 0.0, 1.0, -1.0, 0.0, -0.0, 0.0, 2.0], ZERO_LABELS, False, id="last+0"),
        pytest.param([0.0, -0.0, 1.0, -1.0, -0.0, 0.0, -0.0, 2.0], ZERO_LABELS, True, id="last-0"),
        pytest.param([0.0, -0.0, 0.0, -0.0], [1, 0, 0, 1], True, id="only-zeros"),
    ],
)
def test_a_zero_threshold_has_the_sign_of_the_last_zero_score_as_in_r(
    curve, scores, labels, negative
):
    thresholds = curve(scores, labels).curve["threshold"].to_numpy()
    zeros = thresholds[thresholds == 0]
    assert len(zeros) > 0
    assert (np.signbit(zeros) == negative).all()


def test_an_infinite_score_is_valid():
    result = roc_curve([math.inf, 1.0, -math.inf], [1, 0, 0])
    assert result.auc == 1.0
    assert result.curve["threshold"].tolist() == [math.inf, 1.0, -math.inf, -math.inf]
    assert pr_curve([math.inf, 1.0, -math.inf], [1, 0, 0]).auc_integral == 1.0


def test_pr_curve_accepts_labels_without_a_negative():
    result = pr_curve([3.0, 2.0, 2.0, 1.0], [1, 1, 1, 1])
    assert (result.n_positive, result.n_negative) == (4, 0)
    assert result.auc_integral == 1.0
    assert result.auc_davis_goadrich == 1.0
    assert (result.curve["precision"] == 1.0).all()


@pytest.mark.parametrize("curve", CURVES)
def test_a_curve_takes_bool_labels_a_series_and_a_tuple_as_it_takes_ints(curve):
    scores, labels = curve_samples()["ties_small"]
    expected = curve(scores, labels)
    for variant in (
        [bool(label) for label in labels],
        np.array(labels, dtype=bool),
        np.array(labels, dtype=float),
        pd.Series(labels),
        tuple(labels),
        [np.int64(label) for label in labels],
        np.array(labels, dtype=object),
        pd.array([bool(label) for label in labels], dtype="boolean"),
    ):
        got = curve(pd.Series(scores), variant)
        assert got.curve.equals(expected.curve)
        assert (got.n_positive, got.n_negative) == (expected.n_positive, expected.n_negative)


@pytest.mark.parametrize("curve", CURVES)
def test_a_curve_does_not_change_its_inputs(curve):
    scores = np.array([2.0, 1.0, 2.0, 0.0])
    labels = np.array([1, 0, 0, 1])
    curve(scores, labels)
    assert scores.tolist() == [2.0, 1.0, 2.0, 0.0]
    assert labels.tolist() == [1, 0, 0, 1]


def test_the_results_are_frozen_records_with_float64_frames():
    pr = pr_curve([2.0, 1.0, 0.0], [1, 0, 1])
    roc = roc_curve([2.0, 1.0, 0.0], [1, 0, 1])
    assert isinstance(pr, PrCurve)
    assert isinstance(roc, RocCurve)
    assert [field.name for field in dataclasses.fields(pr)] == [
        "auc_integral",
        "auc_davis_goadrich",
        "curve",
        "n_positive",
        "n_negative",
    ]
    assert [field.name for field in dataclasses.fields(roc)] == [
        "auc",
        "curve",
        "n_positive",
        "n_negative",
    ]
    assert list(pr.curve.columns) == ["recall", "precision", "threshold"]
    assert list(roc.curve.columns) == ["fpr", "tpr", "threshold"]
    assert set(pr.curve.dtypes) == set(roc.curve.dtypes) == {np.dtype("float64")}
    assert {type(pr.auc_integral), type(pr.auc_davis_goadrich), type(roc.auc)} == {float}
    assert {type(pr.n_positive), type(pr.n_negative), type(roc.n_positive)} == {int}
    with pytest.raises(dataclasses.FrozenInstanceError):
        pr.n_positive = 0


@pytest.mark.parametrize("name", ["PrCurve", "RocCurve", "pr_curve", "roc_curve"])
def test_a_curve_name_is_exported_from_the_package_and_from_stats(name):
    assert name in seedmatch.__all__
    assert name in seedmatch.stats.__all__
    assert getattr(seedmatch, name) is getattr(seedmatch.stats, name)


@pytest.mark.parametrize("curve", CURVES)
def test_curves_reject_a_nan_score(curve):
    with pytest.raises(InputError, match=r"^scores\[1\] is NaN; expected a number$"):
        curve([1.0, math.nan, 2.0], [1, 0, 1])


@pytest.mark.parametrize("curve", CURVES)
def test_curves_reject_empty_scores(curve):
    with pytest.raises(InputError, match=r"^scores is empty; expected at least 1 value$"):
        curve([], [])


@pytest.mark.parametrize("curve", CURVES)
def test_curves_reject_a_label_of_2(curve):
    message = rf"^labels\[3\] is 2; {LABEL_MESSAGE}$"
    with pytest.raises(InputError, match=message):
        curve([4.0, 3.0, 2.0, 1.0], [1, 0, 1, 2])


@pytest.mark.parametrize("curve", CURVES)
@pytest.mark.parametrize(
    ("labels", "shown"),
    [
        ([1, 0, -1, 1], "-1"),
        ([1, 0, 0.5, 1], "0.5"),
        ([1, 0, math.nan, 1], "nan"),
        ([1, 0, "1", 1], "'1'"),
        ([1, 0, None, 1], "None"),
        (pd.array([True, False, None, True], dtype="boolean"), "<NA>"),
        (np.array([1.0, 0.0, 2.0, 1.0]), "2.0"),
        (pd.Series([1, 0, 3, 1]), "3"),
    ],
)
def test_curves_reject_a_label_that_is_not_0_or_1(curve, labels, shown):
    message = rf"^labels\[2\] is {re.escape(shown)}; {LABEL_MESSAGE}$"
    with pytest.raises(InputError, match=message):
        curve([4.0, 3.0, 2.0, 1.0], labels)


@pytest.mark.parametrize("curve", CURVES)
def test_curves_reject_labels_of_another_length(curve):
    message = r"^scores has 5 values and labels has 4; expected the same number$"
    with pytest.raises(InputError, match=message):
        curve([5.0, 4.0, 3.0, 2.0, 1.0], [1, 0, 1, 0])


@pytest.mark.parametrize("curve", CURVES)
def test_curves_reject_labels_of_rows_with_different_lengths(curve):
    message = rf"^labels\[0\] is \[1, 0\]; {LABEL_MESSAGE}$"
    with pytest.raises(InputError, match=message):
        curve([4.0, 3.0], [[1, 0], [1]])


@pytest.mark.parametrize("curve", CURVES)
def test_curves_reject_labels_that_are_not_a_vector(curve):
    with pytest.raises(InputError, match=r"^labels has 2 dimensions; expected a 1-D array$"):
        curve([4.0, 3.0, 2.0, 1.0], [[1, 0], [1, 0]])
    with pytest.raises(
        InputError, match=rf"^labels has type int; expected a 1-D array of {LABEL_VALUES}$"
    ):
        curve([4.0, 3.0, 2.0, 1.0], 1)


def test_pr_curve_rejects_labels_without_a_positive():
    message = r"^labels have no positive; expected at least one label 1 or True$"
    with pytest.raises(InputError, match=message):
        pr_curve([3.0, 2.0, 1.0], [0, 0, 0])


def test_roc_curve_rejects_labels_without_a_positive():
    message = r"^labels have no positive; expected at least one label 1 or True$"
    with pytest.raises(InputError, match=message):
        roc_curve([3.0, 2.0, 1.0], [False, False, False])


def test_roc_curve_rejects_labels_without_a_negative():
    message = r"^labels have no negative; expected at least one label 0 or False$"
    with pytest.raises(InputError, match=message):
        roc_curve([3.0, 2.0, 1.0], [1, 1, 1])
