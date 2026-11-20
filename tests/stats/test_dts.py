import itertools
import math
import time
import tracemalloc

import numpy as np
import pytest
from helpers import golden_rows, stats_samples
from hypothesis import given, settings
from hypothesis import strategies as st

from seedmatch.errors import InputError
from seedmatch.stats import PValueMethod, dts_test
from seedmatch.stats.dts import _gaps, _statistics

# R takes E - F and the weights from floating-point divisions and adds the terms in its own order,
# so its statistic differs from the integer form by rounding errors.
STATISTIC_TOLERANCE = {"rel_tol": 1e-9, "abs_tol": 1e-12}
R_ROWS = [pytest.param(row, id=row.case) for row in golden_rows("stats/dts.tsv")]


def definition(x, y):
    """Return the DTS statistic from the three ECDFs at the left value of each gap."""
    pooled = np.sort(np.concatenate((x, y)))
    distinct = np.unique(pooled)
    left, width = distinct[:-1], np.diff(distinct)
    e = np.searchsorted(np.sort(x), left, side="right") / len(x)
    f = np.searchsorted(np.sort(y), left, side="right") / len(y)
    g = np.searchsorted(pooled, left, side="right") / pooled.size
    return float((width * np.abs(e - f) / np.sqrt(2 * g * (1 - g) / pooled.size)).sum())


@pytest.mark.parametrize("row", R_ROWS)
def test_dts_statistic_matches_r(row):
    x, y = stats_samples()[row.case]
    result = dts_test(x, y, rng=1, n_boot=1)
    assert math.isclose(result.statistic, float(row.statistic), **STATISTIC_TOLERANCE)


@pytest.mark.parametrize("row", R_ROWS)
def test_the_definition_of_the_statistic_matches_r(row):
    x, y = stats_samples()[row.case]
    assert math.isclose(definition(x, y), float(row.statistic), **STATISTIC_TOLERANCE)


@st.composite
def samples(draw):
    values = st.one_of(
        st.sampled_from([-2.5, 0.0, 1.0, 4.0, 4.5]),
        st.floats(min_value=-1e3, max_value=1e3, allow_nan=False, allow_infinity=False),
    )
    return (
        draw(st.lists(values, min_size=1, max_size=15)),
        draw(st.lists(values, min_size=1, max_size=15)),
    )


@settings(deadline=None, max_examples=100)
@given(samples())
def test_dts_statistic_equals_the_definition(pair):
    x, y = pair
    result = dts_test(x, y, rng=0, n_boot=1)
    assert math.isclose(result.statistic, definition(x, y), **STATISTIC_TOLERANCE)


@settings(deadline=None, max_examples=100)
@given(samples())
def test_dts_statistic_is_symmetric(pair):
    x, y = pair
    forward = dts_test(x, y, rng=0, n_boot=1).statistic
    backward = dts_test(y, x, rng=0, n_boot=1).statistic
    assert math.isclose(forward, backward, **STATISTIC_TOLERANCE)


def test_dts_with_equal_values_has_the_statistic_0_and_the_p_value_1():
    result = dts_test([1.0, 1.0], [1.0, 1.0, 1.0], rng=0)
    assert (result.statistic, result.p_value) == (0.0, 1.0)


def test_dts_result_carries_the_sizes_and_has_no_alternative():
    result = dts_test([1.0, 2.0, 3.0], [1.5, 2.5], rng=0, n_boot=10)
    assert result.method is PValueMethod.PERMUTATION
    assert result.alternative is None
    assert (result.n_x, result.n_y) == (3, 2)


@pytest.mark.parametrize("row", R_ROWS)
def test_dts_p_value_agrees_with_r_within_the_permutation_error(row):
    x, y = stats_samples()[row.case]
    n_boot, p_r = int(row.n_boot), float(row.p_value)
    result = dts_test(x, y, rng=1, n_boot=n_boot)
    # Each p-value is a fraction of n draws, so the difference has the variance 2 p (1 - p) / n.
    tolerance = 4 * math.sqrt(2 * p_r * (1 - p_r) / n_boot) + 2 / n_boot
    assert abs(result.p_value - p_r) <= tolerance


def test_dts_p_value_is_the_same_for_the_same_seed():
    x, y = stats_samples()["dts_30x40_ties"]
    first = dts_test(x, y, rng=7)
    assert dts_test(x, y, rng=7) == first
    assert dts_test(x, y, rng=np.random.default_rng(7)) == first
    assert dts_test(x, y, rng=np.int64(7)) == first


def test_dts_draws_the_permutations_from_the_generator_that_the_caller_gives():
    x, y = stats_samples()["dts_30x40_ties"]
    generator = np.random.default_rng(7)
    before = generator.bit_generator.state
    dts_test(x, y, rng=generator, n_boot=10)
    assert generator.bit_generator.state != before


def test_dts_p_value_is_never_0():
    # Only 2 of the C(40, 20) arrangements, a fraction of 1e-11, give the observed statistic. Thus
    # no seed reaches it.
    x, y = np.arange(20.0), np.arange(20.0) + 100
    assert dts_test(x, y, rng=1).p_value == 1 / 4000
    assert dts_test(x, y, rng=1, n_boot=10).p_value == 1 / 20


MINIMUM = ([6.0, 0.0, 10.0, 4.0, 2.0, 7.0], [3.0, 8.0, 6.0, 5.0, 1.0])


def test_the_statistic_of_the_minimum_case_is_the_least_of_the_arrangements():
    x, y = MINIMUM
    pooled = sorted(x + y)
    statistics = [
        definition(
            [pooled[k] for k in chosen],
            [pooled[k] for k in range(len(pooled)) if k not in set(chosen)],
        )
        for chosen in itertools.combinations(range(len(pooled)), len(x))
    ]
    assert math.isclose(definition(x, y), min(statistics), rel_tol=1e-9)


@pytest.mark.parametrize("n_boot", [1, 63, 64, 65, 128, 130])
def test_dts_p_value_counts_each_permutation_once_for_each_block_size(n_boot):
    # The observed statistic is the least of all arrangements, so each permutation counts. The count
    # must be n_boot, not less and not more.
    x, y = MINIMUM
    for seed in range(5):
        assert dts_test(x, y, rng=seed, n_boot=n_boot).p_value == 1.0


@st.composite
def labelled_samples(draw):
    x, y = draw(samples())
    return x, y, draw(st.permutations([True] * len(x) + [False] * len(y)))


@settings(deadline=None, max_examples=200)
@given(labelled_samples())
def test_the_statistic_of_the_same_labels_is_bit_equal_alone_and_in_a_block(case):
    x, y, labels = case
    ends, expected, scale = _gaps(np.sort(np.concatenate((x, y))), len(x), len(y))
    row = np.array(labels)[np.newaxis, :]
    alone = _statistics(row, ends, expected, scale)[0]
    for size in (2, 5, 64):
        block = np.tile(row, (size, 1))
        assert (_statistics(block, ends, expected, scale) == alone).all()


def test_the_statistic_of_the_same_labels_is_bit_equal_alone_and_in_a_block_for_a_long_sample():
    generator = np.random.default_rng(4)
    values = np.sort(generator.normal(size=700))
    ends, expected, scale = _gaps(values, 300, 400)
    row = (generator.permutation(700) < 300)[np.newaxis, :]
    alone = _statistics(row, ends, expected, scale)[0]
    for size in (2, 5, 64):
        block = np.tile(row, (size, 1))
        assert (_statistics(block, ends, expected, scale) == alone).all()


EXACT_CASES = {
    "no ties": ([0.1, 1.6, 2.9, 4.1], [0.7, 1.0, 3.5, 3.9, 5.2]),
    "ties": ([1.0, 2.0, 2.0, 4.0], [2.0, 3.0, 3.0, 4.0, 5.0]),
}


@pytest.mark.parametrize("case", EXACT_CASES)
def test_dts_p_value_estimates_the_fraction_of_the_arrangements_that_are_not_less(case):
    x, y = EXACT_CASES[case]
    pooled = sorted(x + y)
    observed = definition(x, y)
    statistics = np.array(
        [
            definition(
                [pooled[k] for k in chosen],
                [pooled[k] for k in range(len(pooled)) if k not in set(chosen)],
            )
            for chosen in itertools.combinations(range(len(pooled)), len(x))
        ]
    )
    # An arrangement with the same statistic in exact arithmetic can fall on either side in floats.
    most = float((statistics >= observed * (1 - 1e-9)).mean())
    least = float((statistics >= observed * (1 + 1e-9)).mean())
    n_boot = 20000
    margin = 5 * math.sqrt(0.25 / n_boot)
    estimate = dts_test(x, y, rng=5, n_boot=n_boot).p_value
    assert 0.05 < least <= most < 0.95
    assert least - margin <= estimate <= most + margin


def test_dts_does_not_change_its_inputs():
    x, y = np.array([3.0, 1.0, 2.0]), np.array([2.5, 0.5])
    dts_test(x, y, rng=0, n_boot=10)
    assert x.tolist() == [3.0, 1.0, 2.0]
    assert y.tolist() == [2.5, 0.5]


def test_dts_rejects_an_infinite_value_D034():
    with pytest.raises(InputError, match=r"^x\[2\] inf is infinite; expected a finite number"):
        dts_test([1.0, 2.0, np.inf], [1.0], rng=0)
    with pytest.raises(InputError, match=r"^y\[0\] -inf is infinite; expected a finite number"):
        dts_test([1.0], [-np.inf, 2.0], rng=0)


def test_dts_rejects_an_empty_sample_D034():
    with pytest.raises(InputError, match=r"^x is empty; expected at least 1 value"):
        dts_test([], [1.0], rng=0)
    with pytest.raises(InputError, match=r"^y is empty; expected at least 1 value"):
        dts_test([1.0], [], rng=0)


def test_dts_rejects_nan():
    with pytest.raises(InputError, match=r"^x\[1\] is NaN"):
        dts_test([1.0, np.nan], [1.0], rng=0)
    with pytest.raises(InputError, match=r"^y\[0\] is NaN"):
        dts_test([1.0], [np.nan], rng=0)


@pytest.mark.parametrize(
    ("n_boot", "message"),
    [
        (0, r"n_boot 0 is less than 1; expected 1 or more"),
        (-5, r"n_boot -5 is less than 1; expected 1 or more"),
        (True, r"n_boot True has type bool; expected an integer"),
        (2000.0, r"n_boot 2000\.0 has type float; expected an integer"),
        ("2000", r"n_boot '2000' has type str; expected an integer"),
        (2**63, r"n_boot 9223372036854775808 is too large; expected at most"),
    ],
)
def test_dts_rejects_a_bad_number_of_permutations(n_boot, message):
    with pytest.raises(InputError, match=f"^{message}"):
        dts_test([1.0, 2.0], [3.0], rng=0, n_boot=n_boot)


def test_dts_accepts_a_numpy_integer_for_the_number_of_permutations():
    assert dts_test([1.0, 2.0], [3.0], rng=0, n_boot=np.int64(5)).method is PValueMethod.PERMUTATION


def test_dts_rejects_a_negative_seed():
    with pytest.raises(InputError, match=r"^rng -1 is negative; expected 0 or more"):
        dts_test([1.0, 2.0], [3.0], rng=-1)


@pytest.mark.parametrize(
    "rng", [None, 1.5, "1", True, np.random.RandomState(1), [1], np.array([1, 2])]
)
def test_dts_rejects_an_rng_that_is_neither_a_generator_nor_an_int(rng):
    message = r"^rng .* has type \w+; expected a numpy\.random\.Generator or an int"
    with pytest.raises(InputError, match=message):
        dts_test([1.0, 2.0], [3.0], rng=rng)


def test_dts_accepts_an_int_seed_of_any_size():
    x, y = stats_samples()["dts_30x40_ties"]
    seed = 2**200
    assert dts_test(x, y, rng=seed) == dts_test(x, y, rng=np.random.default_rng(seed))


def peak_bytes(x, y, n_boot):
    tracemalloc.start()
    try:
        dts_test(x, y, rng=1, n_boot=n_boot)
        return tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


def test_dts_memory_does_not_grow_with_the_number_of_permutations():
    generator = np.random.default_rng(2)
    x, y = generator.normal(size=500), generator.normal(size=4500)
    # One block of 64 permutations of 5000 values holds a few arrays of 2.5 MB. All permutations
    # together would hold 51 MB for each such array.
    few, many = peak_bytes(x, y, 128), peak_bytes(x, y, 1280)
    assert many <= 1.2 * few
    assert many < 40e6


def test_dts_runs_on_a_large_table():
    generator = np.random.default_rng(3)
    x, y = generator.normal(0.3, 1, size=2000), generator.normal(0, 1, size=18000)
    start = time.perf_counter()
    result = dts_test(x, y, rng=1, n_boot=2000)
    elapsed = time.perf_counter() - start
    assert elapsed < 30
    assert 0 < result.p_value < 0.01
    assert (result.n_x, result.n_y) == (2000, 18000)
