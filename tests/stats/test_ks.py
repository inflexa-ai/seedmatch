import itertools
import math
import re
import warnings
from collections.abc import Sequence

import numpy as np
import pytest
from helpers import STATISTIC_TOLERANCE, golden_row, golden_rows, stats_samples
from hypothesis import given, settings
from hypothesis import strategies as st
from scipy.special import kolmogorov

from seedmatch.errors import InputError
from seedmatch.stats import Alternative, PValueMethod, ks_test

P_VALUE_TOLERANCE = {"rel_tol": 1e-9, "abs_tol": 0.0}
D032_ROWS = {
    ("ks_asym_cut_series", "two.sided"),
    ("ks_asym_small_p", "greater"),
    ("ks_100x100_shift", "two.sided"),
    ("ks_100x100_shift", "greater"),
    ("ks_120x150_shift", "greater"),
    ("ks_400x300_shift", "greater"),
}


def run_case(case, alternative):
    x, y = stats_samples()[case]
    return ks_test(x, y, alternative=Alternative(alternative))


def limit_inputs(case, alternative):
    """Return (en, D) of the asymptotic formulas, with D from the lattice value that R rounds."""
    x, y = stats_samples()[case]
    m, n = len(x), len(y)
    d = round(float(golden_row("stats/ks.tsv", case, alternative).statistic) * m * n)
    return m * n / (m + n), d / (m * n)


@pytest.mark.parametrize(
    "row",
    [
        pytest.param(row, id=f"{row.case}-{row.alternative}")
        for row in golden_rows("stats/ks.tsv")
        if (row.case, row.alternative) not in D032_ROWS
    ],
)
def test_ks_matches_r(row):
    result = run_case(row.case, row.alternative)
    assert math.isclose(result.statistic, float(row.statistic), **STATISTIC_TOLERANCE)
    assert math.isclose(result.p_value, float(row.p_value), **P_VALUE_TOLERANCE)
    assert result.method.value == row.method


def test_ks_result_carries_the_sizes_and_the_alternative():
    result = ks_test([1.0, 2.0, 3.0], [1.5, 2.5], alternative=Alternative.LESS)
    assert (result.n_x, result.n_y, result.alternative) == (3, 2, Alternative.LESS)


def test_ks_default_alternative_is_the_default_of_r():
    x, y = stats_samples()["ks_5x7_plain"]
    assert ks_test(x, y) == ks_test(x, y, alternative=Alternative.TWO_SIDED)
    assert ks_test(x, y).alternative is Alternative.TWO_SIDED


def test_ks_rejects_the_alternative_as_text():
    message = (
        "alternative 'less' is not valid; "
        "expected one of Alternative.TWO_SIDED, Alternative.LESS, Alternative.GREATER"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        ks_test([1.0], [2.0], alternative="less")


def test_ks_rejects_an_unknown_alternative():
    with pytest.raises(
        InputError,
        match=r"alternative 'sideways' is not valid; expected one of Alternative\.TWO_SIDED",
    ):
        ks_test([1.0], [2.0], alternative="sideways")


@pytest.mark.parametrize(
    ("m", "n", "method"),
    [
        (1, 9999, PValueMethod.EXACT),
        (9999, 1, PValueMethod.EXACT),
        (1, 10000, PValueMethod.ASYMPTOTIC),
        (2, 5000, PValueMethod.ASYMPTOTIC),
    ],
)
def test_ks_uses_the_exact_p_value_only_below_10000_pairs(m, n, method):
    generator = np.random.default_rng(7)
    result = ks_test(generator.normal(size=m), generator.normal(size=n))
    assert result.method is method


def test_ks_gives_no_warning_for_ties_D037():
    # R 4.5.3 warns "p-value will be approximate in the presence of ties" for these samples.
    x, y = np.repeat([1.0, 2.0, 3.0], 40), np.repeat([2.0, 3.0, 4.0], 50)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        result = ks_test(x, y)
    assert result.method is PValueMethod.ASYMPTOTIC


def test_ks_with_equal_samples_has_the_statistic_0_and_the_p_value_1():
    for alternative in Alternative:
        result = ks_test([1.0, 2.0, 2.0], [1.0, 2.0, 2.0], alternative=alternative)
        assert (result.statistic, result.p_value) == (0.0, 1.0)


def test_ks_less_of_two_samples_is_greater_of_the_samples_in_the_other_order():
    generator = np.random.default_rng(11)
    x, y = generator.integers(0, 6, size=9), generator.integers(0, 6, size=8)
    less = ks_test(x, y, alternative=Alternative.LESS)
    greater = ks_test(y, x, alternative=Alternative.GREATER)
    assert (less.statistic, less.p_value) == (greater.statistic, greater.p_value)


def test_ks_does_not_change_its_inputs():
    x, y = np.array([3.0, 1.0, 2.0]), np.array([2.5, 0.5])
    ks_test(x, y)
    assert x.tolist() == [3.0, 1.0, 2.0]
    assert y.tolist() == [2.5, 0.5]


@pytest.mark.parametrize("alternative", list(Alternative))
def test_ks_accepts_an_infinite_value(alternative):
    # An infinite value is the greatest or the least value of the pooled sample.
    with_inf = ks_test([1.0, np.inf, -np.inf], [0.0, 2.0], alternative=alternative)
    with_numbers = ks_test([1.0, 100.0, -100.0], [0.0, 2.0], alternative=alternative)
    assert with_inf == with_numbers
    assert not math.isnan(with_inf.statistic)
    assert ks_test([np.inf], [np.inf], alternative=alternative).p_value == 1.0


def test_ks_rejects_an_empty_sample():
    with pytest.raises(InputError, match=r"^x is empty"):
        ks_test([], [1.0])
    with pytest.raises(InputError, match=r"^y is empty"):
        ks_test([1.0], [])


def test_ks_rejects_nan():
    with pytest.raises(InputError, match=r"^x\[1\] is NaN"):
        ks_test([1.0, np.nan], [1.0])
    with pytest.raises(InputError, match=r"^y\[0\] is NaN"):
        ks_test([1.0], [np.nan, 2.0])


def numerators(x: Sequence[float], y: Sequence[float]) -> dict[Alternative, int]:
    """Return d of each alternative from the ECDFs at each pooled value, with integer arithmetic."""
    m, n = len(x), len(y)
    gaps = [sum(v <= t for v in x) * n - sum(v <= t for v in y) * m for t in sorted({*x, *y})]
    return {
        Alternative.GREATER: max(gaps),
        Alternative.LESS: max(-gap for gap in gaps),
        Alternative.TWO_SIDED: max(abs(gap) for gap in gaps),
    }


@st.composite
def small_samples_with_ties(draw):
    m = draw(st.integers(min_value=1, max_value=11))
    n = draw(st.integers(min_value=1, max_value=12 - m))
    values = st.sampled_from([-2.5, 0.0, 1.0, 4.0])
    return draw(st.lists(values, min_size=m, max_size=m)), draw(
        st.lists(values, min_size=n, max_size=n)
    )


@settings(deadline=None, max_examples=40)
@given(small_samples_with_ties())
def test_exact_p_value_equals_the_count_of_the_permutations(pair):
    x, y = pair
    m = len(x)
    observed = numerators(x, y)
    pooled = sorted(x + y)
    at_least = dict.fromkeys(Alternative, 0)
    total = 0
    for chosen in itertools.combinations(range(len(pooled)), m):
        in_x = set(chosen)
        permuted = numerators(
            [pooled[k] for k in chosen], [pooled[k] for k in range(len(pooled)) if k not in in_x]
        )
        total += 1
        for alternative in Alternative:
            at_least[alternative] += permuted[alternative] >= observed[alternative]
    for alternative in Alternative:
        result = ks_test(x, y, alternative=alternative)
        assert result.method is PValueMethod.EXACT
        assert result.statistic == observed[alternative] / (m * len(y))
        assert math.isclose(result.p_value, at_least[alternative] / total, rel_tol=0, abs_tol=1e-12)


def test_ks_asymptotic_two_sided_uses_the_limit_distribution_D032():
    row = golden_row("stats/ks.tsv", "ks_asym_cut_series", "two.sided")
    result = run_case(row.case, row.alternative)
    en, d = limit_inputs(row.case, row.alternative)
    expected = float(kolmogorov(math.sqrt(en) * d))
    assert 0.9 < math.sqrt(en) * d < 1
    assert result.method is PValueMethod.ASYMPTOTIC
    assert math.isclose(result.p_value, expected, rel_tol=1e-12)
    # The cut series of R differs from the limit distribution in the fifth digit.
    assert abs(float(row.p_value) - result.p_value) / result.p_value > 1e-7


def test_ks_asymptotic_one_sided_keeps_a_small_p_value_D032():
    row = golden_row("stats/ks.tsv", "ks_asym_small_p", "greater")
    result = run_case(row.case, row.alternative)
    en, d = limit_inputs(row.case, row.alternative)
    assert float(row.p_value) == 0
    assert 0 < result.p_value < 1e-17
    assert math.isclose(result.p_value, math.exp(-2 * en * d**2), rel_tol=1e-12)


@pytest.mark.parametrize(
    ("case", "alternative"),
    [
        ("ks_100x100_shift", "two.sided"),
        ("ks_100x100_shift", "greater"),
        ("ks_120x150_shift", "greater"),
        ("ks_400x300_shift", "greater"),
    ],
)
def test_ks_keeps_a_p_value_that_r_rounds_to_zero_D032(case, alternative):
    row = golden_row("stats/ks.tsv", case, alternative)
    result = run_case(case, alternative)
    en, d = limit_inputs(case, alternative)
    if alternative == "two.sided":
        expected = float(kolmogorov(math.sqrt(en) * d))
    else:
        expected = math.exp(-2 * en * d**2)
    assert float(row.p_value) == 0
    assert 0 < result.p_value < 1e-40
    assert math.isclose(result.p_value, expected, rel_tol=1e-12)
    assert math.isclose(result.statistic, float(row.statistic), **STATISTIC_TOLERANCE)


def test_ks_statistic_cannot_exceed_1_where_the_sum_of_r_does_D032():
    x, y = stats_samples()["ks_100x100_shift"]
    assert float(golden_row("stats/ks.tsv", "ks_100x100_shift", "two.sided").statistic) > 1
    assert ks_test(x, y).statistic == 1.0
