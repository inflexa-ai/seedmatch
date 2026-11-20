import itertools
import math
import re
import sys
import warnings

import numpy as np
import pytest
from helpers import golden_row, golden_rows, stats_samples
from hypothesis import given, settings
from hypothesis import strategies as st

from seedmatch.errors import InputError
from seedmatch.stats import Alternative, PValueMethod, wilcoxon_test

# W is a multiple of 0.5, and each p-value comes from a few float operations or from sums of
# non-negative counts, so a relative tolerance of 1e-9 leaves a wide margin.
TOLERANCE = {"rel_tol": 1e-9, "abs_tol": 0.0}
ALL_EQUAL = "wilcoxon_all_equal"


@pytest.mark.parametrize(
    "row",
    [
        pytest.param(row, id=f"{row.case}-{row.alternative}")
        for row in golden_rows("stats/wilcoxon.tsv")
        if row.case != ALL_EQUAL
    ],
)
def test_wilcoxon_matches_r(row):
    x, y = stats_samples()[row.case]
    result = wilcoxon_test(x, y, alternative=Alternative(row.alternative))
    assert math.isclose(result.statistic, float(row.statistic), **TOLERANCE)
    assert math.isclose(result.p_value, float(row.p_value), **TOLERANCE)
    assert result.method.value == row.method


def test_wilcoxon_result_carries_the_sizes_and_the_alternative():
    result = wilcoxon_test([1.0, 2.0, 3.0], [1.5, 2.5], alternative=Alternative.LESS)
    assert (result.n_x, result.n_y, result.alternative) == (3, 2, Alternative.LESS)


def test_wilcoxon_default_alternative_is_the_default_of_r():
    x, y = stats_samples()["wilcoxon_10x12_plain"]
    assert wilcoxon_test(x, y) == wilcoxon_test(x, y, alternative=Alternative.TWO_SIDED)
    assert wilcoxon_test(x, y).alternative is Alternative.TWO_SIDED


def test_wilcoxon_rejects_the_alternative_as_text():
    message = (
        "alternative 'less' is not valid; "
        "expected one of Alternative.TWO_SIDED, Alternative.LESS, Alternative.GREATER"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        wilcoxon_test([1.0], [2.0], alternative="less")


def test_wilcoxon_rejects_an_unknown_alternative():
    with pytest.raises(
        InputError,
        match=r"alternative 'sideways' is not valid; expected one of Alternative\.TWO_SIDED",
    ):
        wilcoxon_test([1.0], [2.0], alternative="sideways")


def test_wilcoxon_rejects_an_empty_sample():
    with pytest.raises(InputError, match=r"^x is empty"):
        wilcoxon_test([], [1.0])
    with pytest.raises(InputError, match=r"^y is empty"):
        wilcoxon_test([1.0], [])


def test_wilcoxon_rejects_nan():
    with pytest.raises(InputError, match=r"^x\[1\] is NaN"):
        wilcoxon_test([1.0, np.nan], [1.0])
    with pytest.raises(InputError, match=r"^y\[0\] is NaN"):
        wilcoxon_test([1.0], [np.nan, 2.0])


@pytest.mark.parametrize("alternative", list(Alternative))
def test_wilcoxon_accepts_an_infinite_value(alternative):
    # Only the order of the values matters, thus an infinite value is the greatest or the least.
    with_inf = wilcoxon_test([1.0, np.inf, -np.inf], [0.0, 2.0], alternative=alternative)
    with_numbers = wilcoxon_test([1.0, 100.0, -100.0], [0.0, 2.0], alternative=alternative)
    assert with_inf == with_numbers
    assert wilcoxon_test([np.inf], [np.inf], alternative=alternative).p_value == 1.0


def test_wilcoxon_gives_no_warning_for_ties_D037():
    # R 4.5.3 warns "cannot compute exact p-value with ties" for these samples.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        result = wilcoxon_test([1.0, 2.0, 3.0, 3.0], [3.0, 4.0, 5.0])
    assert result.method is PValueMethod.ASYMPTOTIC


def test_wilcoxon_does_not_change_its_inputs():
    x, y = np.array([3.0, 1.0, 2.0]), np.array([2.5, 0.5])
    wilcoxon_test(x, y)
    assert x.tolist() == [3.0, 1.0, 2.0]
    assert y.tolist() == [2.5, 0.5]


@pytest.mark.parametrize(
    ("x", "y", "method"),
    [
        ([1.0, 2.0, 3.0], [3.5, 4.0], PValueMethod.EXACT),
        ([1.0, 2.0, 3.0], [3.0, 4.0], PValueMethod.ASYMPTOTIC),
        ([1.0, 1.0, 3.0], [3.5, 4.0], PValueMethod.ASYMPTOTIC),
        (np.arange(49.0), np.arange(49.0) + 0.5, PValueMethod.EXACT),
        (np.arange(50.0), np.arange(49.0) + 0.5, PValueMethod.ASYMPTOTIC),
        (np.arange(49.0), np.arange(50.0) + 0.5, PValueMethod.ASYMPTOTIC),
    ],
)
def test_wilcoxon_uses_the_exact_p_value_only_for_small_samples_without_ties(x, y, method):
    assert wilcoxon_test(x, y).method is method


@pytest.mark.parametrize(
    ("decimals", "method"), [(None, PValueMethod.EXACT), (0, PValueMethod.ASYMPTOTIC)]
)
def test_wilcoxon_less_of_two_samples_is_greater_of_the_samples_in_the_other_order(
    decimals, method
):
    generator = np.random.default_rng(11)
    x, y = generator.normal(size=9), generator.normal(size=8)
    if decimals is not None:
        x, y = np.round(x, decimals), np.round(y, decimals)
    less = wilcoxon_test(x, y, alternative=Alternative.LESS)
    greater = wilcoxon_test(y, x, alternative=Alternative.GREATER)
    assert less.method is greater.method is method
    assert math.isclose(less.p_value, greater.p_value, **TOLERANCE)
    assert less.statistic + greater.statistic == 9 * 8


def test_wilcoxon_exact_tail_keeps_its_relative_precision():
    # The 49 + 49 ranks have C(98, 49) arrangements, more than 2^63. The lowest rank sum has one.
    below, above = np.arange(49.0), np.arange(49.0) + 100
    one_in = math.comb(98, 49)
    less = wilcoxon_test(below, above, alternative=Alternative.LESS)
    assert less.method is PValueMethod.EXACT
    assert less.statistic == 0
    assert math.isclose(less.p_value, 1 / one_in, rel_tol=1e-12)
    two_sided = wilcoxon_test(below, above)
    assert math.isclose(two_sided.p_value, 2 / one_in, rel_tol=1e-12)
    assert wilcoxon_test(below, above, alternative=Alternative.GREATER).p_value == 1.0
    most = wilcoxon_test(above, below, alternative=Alternative.GREATER)
    assert most.statistic == 49 * 49
    assert math.isclose(most.p_value, 1 / one_in, rel_tol=1e-12)


def test_wilcoxon_asymptotic_upper_tail_keeps_a_small_p_value():
    # Each sample is one tie group. 1 - ndtr(z) would give 0 here, but the p-value is about 1e-45.
    x, y = [1.0] * 100, [0.0] * 100
    result = wilcoxon_test(x, y, alternative=Alternative.GREATER)
    sigma = math.sqrt(10000 / 12 * (201 - 2 * (100**3 - 100) / (200 * 199)))
    z = (5000 - 0.5) / sigma
    expected = 0.5 * math.erfc(z / math.sqrt(2))
    assert result.method is PValueMethod.ASYMPTOTIC
    assert 0 < expected < 1e-40
    assert math.isclose(result.p_value, expected, rel_tol=1e-9)
    two_sided = wilcoxon_test(x, y)
    assert math.isclose(two_sided.p_value, 2 * expected, rel_tol=1e-9)


# R 4.5.3: wilcox.test(seq_len(948) + 1000, seq_len(948), alternative = "greater")$p.value
R_SEPARATED_GREATER = 2.5372015660847956e-311
# R 4.5.3: wilcox.test(seq_len(948) + 1000, seq_len(948))$p.value
R_SEPARATED_TWO_SIDED = 5.0744031321695912e-311


def test_wilcoxon_asymptotic_p_value_keeps_the_subnormal_floats_of_r():
    # ndtr(-z) is 0.0 for z above about 37.5, but a float64 holds p-values down to 5e-324.
    x, y = np.arange(948.0) + 1000, np.arange(948.0)
    greater = wilcoxon_test(x, y, alternative=Alternative.GREATER)
    two_sided = wilcoxon_test(x, y)
    assert greater.method is PValueMethod.ASYMPTOTIC
    assert 0 < greater.p_value < sys.float_info.min
    assert math.isclose(greater.p_value, R_SEPARATED_GREATER, **TOLERANCE)
    assert math.isclose(two_sided.p_value, R_SEPARATED_TWO_SIDED, **TOLERANCE)
    assert wilcoxon_test(x, y, alternative=Alternative.LESS).p_value == 1.0


def test_wilcoxon_p_value_does_not_pass_1():
    # The float sum of all counts over C(64, 20) is 1.0000000000000002 without the clamp.
    x, y = np.arange(20.0), np.arange(44.0) + 100
    result = wilcoxon_test(x, y, alternative=Alternative.GREATER)
    assert result.method is PValueMethod.EXACT
    assert result.p_value == 1.0


@st.composite
def small_samples_without_ties(draw):
    m = draw(st.integers(min_value=1, max_value=11))
    n = draw(st.integers(min_value=1, max_value=12 - m))
    values = draw(
        st.lists(st.integers(-1000, 1000), min_size=m + n, max_size=m + n, unique=True).map(
            lambda drawn: [float(value) for value in drawn]
        )
    )
    return values[:m], values[m:]


@settings(deadline=None, max_examples=60)
@given(small_samples_without_ties())
def test_exact_p_value_equals_the_count_of_the_subsets(pair):
    x, y = pair
    m, n = len(x), len(y)
    pooled = sorted(x + y)
    observed = sum(pooled.index(value) + 1 for value in x) - m * (m + 1) // 2
    sums = [
        sum(chosen) - m * (m + 1) // 2 for chosen in itertools.combinations(range(1, m + n + 1), m)
    ]
    at_most = sum(w <= observed for w in sums) / len(sums)
    at_least = sum(w >= observed for w in sums) / len(sums)
    expected = {
        Alternative.LESS: at_most,
        Alternative.GREATER: at_least,
        Alternative.TWO_SIDED: min(1.0, 2 * (at_least if observed > m * n / 2 else at_most)),
    }
    assert len(sums) == math.comb(m + n, m)
    for alternative in Alternative:
        result = wilcoxon_test(x, y, alternative=alternative)
        assert result.method is PValueMethod.EXACT
        assert result.statistic == observed
        assert math.isclose(result.p_value, expected[alternative], rel_tol=1e-12, abs_tol=0.0)


@pytest.mark.parametrize("alternative", list(Alternative))
def test_wilcoxon_with_all_values_equal_gives_1_D033(alternative):
    x, y = stats_samples()[ALL_EQUAL]
    golden = golden_row("stats/wilcoxon.tsv", ALL_EQUAL, alternative.value)
    result = wilcoxon_test(x, y, alternative=alternative)
    # R gives NaN for the two-sided test and 1 for the one-sided tests.
    assert golden.p_value == ("NA" if alternative is Alternative.TWO_SIDED else "1")
    assert math.isclose(result.statistic, float(golden.statistic), **TOLERANCE)
    assert result.statistic == 17.5
    assert result.p_value == 1.0
    assert result.method is PValueMethod.ASYMPTOTIC
    assert result.method.value == golden.method
