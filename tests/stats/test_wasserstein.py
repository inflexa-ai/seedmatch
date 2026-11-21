import dataclasses
import math
import time
import tracemalloc
import warnings
from bisect import bisect_left
from fractions import Fraction

import numpy as np
import pandas as pd
import pytest
from helpers import (
    WASSERSTEIN_CASES,
    exactly,
    golden,
    golden_rows,
    stats_samples,
    wasserstein_samples,
)
from hypothesis import given, settings
from hypothesis import strategies as st

from seedmatch.errors import InputError
from seedmatch.stats import WassersteinAreas, wasserstein_areas, wasserstein_difference
from seedmatch.stats.wasserstein import _unrounded_areas

R_ROWS = [
    pytest.param(row, id=f"{row.case}-{row.n_interp}")
    for row in golden_rows("stats/wasserstein.tsv")
]
CASES = sorted({row.case for row in golden_rows("stats/wasserstein.tsv")})


@pytest.mark.parametrize("row", R_ROWS)
def test_wasserstein_matches_r(row):
    seed, non_seed = stats_samples()[row.case]
    areas = wasserstein_areas(seed, non_seed, n_interp=int(row.n_interp))
    assert (areas.absolute, areas.negative_censored, areas.positive_censored) == (
        float(row.abs_auc),
        float(row.neg_cens_auc),
        float(row.pos_cens_auc),
    )


def test_wasserstein_uses_10000_grid_points_by_default():
    seed, non_seed = stats_samples()["wass_40x300_ties"]
    assert wasserstein_areas(seed, non_seed) == wasserstein_areas(seed, non_seed, n_interp=10000)


@pytest.mark.parametrize("case", CASES)
def test_the_censored_areas_swap_with_the_samples(case):
    seed, non_seed = stats_samples()[case]
    forward = wasserstein_areas(seed, non_seed, n_interp=1000)
    backward = wasserstein_areas(non_seed, seed, n_interp=1000)
    assert backward.negative_censored == forward.positive_censored
    assert backward.positive_censored == forward.negative_censored
    assert backward.absolute == forward.absolute


def test_the_areas_are_plain_floats_in_a_frozen_record():
    seed, non_seed = stats_samples()["wass_10x10_plain"]
    areas = wasserstein_areas(seed, non_seed, n_interp=100)
    assert isinstance(areas, WassersteinAreas)
    assert {type(value) for value in dataclasses.astuple(areas)} == {float}
    with pytest.raises(dataclasses.FrozenInstanceError):
        areas.absolute = 0.0


def test_the_censored_areas_add_up_to_the_absolute_area():
    seed, non_seed = stats_samples()["wass_40x300_plain"]
    unrounded = _unrounded_areas(np.array(seed), np.array(non_seed), 10000)
    assert math.isclose(unrounded[0], unrounded[1] + unrounded[2], rel_tol=1e-12)


def exact_areas(seed, non_seed, n_interp):
    """Return the three areas in exact rational arithmetic, from the definition."""
    pooled = [Fraction(value) for value in (*seed, *non_seed)]
    low, high = min(pooled), max(pooled)
    grid = [low + k * (high - low) / (n_interp - 1) for k in range(n_interp)]

    def ecdf(values):
        knots = sorted({Fraction(value) for value in values})
        heights = [
            Fraction(sum(Fraction(value) <= knot for value in values), len(values))
            for knot in knots
        ]

        def at(v):
            if v < knots[0]:
                return Fraction(0)
            if v > knots[-1]:
                return Fraction(1)
            j = bisect_left(knots, v)
            if knots[j] == v:
                return heights[j]
            i = j - 1
            return heights[i] + (heights[j] - heights[i]) * (v - knots[i]) / (knots[j] - knots[i])

        return at

    seed_at, non_seed_at = ecdf(seed), ecdf(non_seed)
    diff = [seed_at(v) - non_seed_at(v) for v in grid]

    def trapezoid(heights):
        return float(
            sum((grid[k] - grid[k - 1]) * (heights[k] + heights[k - 1]) for k in range(1, n_interp))
            / 2
        )

    return (
        trapezoid([abs(d) for d in diff]),
        trapezoid([max(d, 0) for d in diff]),
        trapezoid([max(-d, 0) for d in diff]),
    )


@st.composite
def dyadic_samples(draw):
    # A multiple of 1/8 and a grid with 2^k + 1 points keep the grid exact, thus no grid point falls
    # next to a knot by rounding. The small range of values makes many equal values.
    values = st.integers(min_value=-12, max_value=12).map(lambda k: k / 8)
    seed = draw(st.lists(values, min_size=10, max_size=25).filter(lambda v: len(set(v)) > 1))
    non_seed = draw(st.lists(values, min_size=10, max_size=25).filter(lambda v: len(set(v)) > 1))
    return seed, non_seed, draw(st.sampled_from([2, 3, 5, 9, 17, 33, 65]))


@settings(deadline=None, max_examples=300)
@given(dyadic_samples())
def test_wasserstein_areas_equal_the_definition_in_exact_arithmetic(case):
    seed, non_seed, n_interp = case
    got = _unrounded_areas(np.array(seed), np.array(non_seed), n_interp)
    expected = exact_areas(seed, non_seed, n_interp)
    for value, reference in zip(got, expected, strict=True):
        assert math.isclose(value, reference, rel_tol=1e-9, abs_tol=1e-12)


def test_wasserstein_gives_no_warning_for_ties_D037():
    # R 4.5.3 warns "collapsing to unique 'x' values" for each sample with equal values.
    seed = [1.0, 2.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
    non_seed = [1.0, 1.0, 2.0, 3.0, 3.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        areas = wasserstein_areas(seed, non_seed, n_interp=100)
    assert areas.absolute > 0


def test_wasserstein_does_not_change_its_inputs():
    seed, non_seed = np.arange(10.0)[::-1].copy(), np.arange(12.0)[::-1] + 0.5
    wasserstein_areas(seed, non_seed, n_interp=50)
    assert seed.tolist() == list(range(9, -1, -1))
    assert non_seed.tolist() == [11.5 - k for k in range(12)]


def test_wasserstein_accepts_a_series_and_a_tuple():
    seed, non_seed = stats_samples()["wass_10x10_plain"]
    expected = wasserstein_areas(seed, non_seed, n_interp=100)
    assert wasserstein_areas(pd.Series(seed), tuple(non_seed), n_interp=100) == expected


def test_wasserstein_rejects_fewer_than_10_values():
    nine = [float(k) for k in range(9)]
    ten = [float(k) for k in range(10)]
    with pytest.raises(InputError, match=r"^seed has 9 values; expected at least 10$"):
        wasserstein_areas(nine, ten)
    with pytest.raises(InputError, match=r"^non_seed has 9 values; expected at least 10$"):
        wasserstein_areas(ten, nine)


@pytest.mark.parametrize(
    ("case", "name", "count"),
    [("wass_seed_one_value", "seed", 10), ("wass_non_seed_one_value", "non_seed", 20)],
)
def test_wasserstein_rejects_one_distinct_value(case, name, count):
    seed, non_seed = stats_samples()[case]
    message = rf"^{name} has {count} values that are all 0\.5; expected at least 2 distinct values$"
    with pytest.raises(InputError, match=message):
        wasserstein_areas(seed, non_seed)


def test_wasserstein_treats_negative_zero_and_zero_as_one_value():
    with pytest.raises(
        InputError, match=r"^seed has 10 values that are all 0\.0; expected at least 2"
    ):
        wasserstein_areas([0.0, -0.0] * 5, [float(k) for k in range(10)])


@pytest.mark.parametrize("n_interp", [1, 0, -5])
def test_wasserstein_rejects_n_interp_below_2_D038(n_interp):
    seed, non_seed = stats_samples()["wass_10x10_plain"]
    message = rf"^n_interp {n_interp} is less than 2; expected 2 or more$"
    with pytest.raises(InputError, match=message):
        wasserstein_areas(seed, non_seed, n_interp=n_interp)


@pytest.mark.parametrize(
    ("n_interp", "message"),
    [
        (True, r"n_interp True has type bool; expected an integer"),
        (100.0, r"n_interp 100\.0 has type float; expected an integer"),
        ("100", r"n_interp '100' has type str; expected an integer"),
        (2**63, r"n_interp 9223372036854775808 is too large; expected at most"),
    ],
)
def test_wasserstein_rejects_an_n_interp_that_is_not_a_usable_integer(n_interp, message):
    seed, non_seed = stats_samples()["wass_10x10_plain"]
    with pytest.raises(InputError, match=f"^{message}"):
        wasserstein_areas(seed, non_seed, n_interp=n_interp)


def test_wasserstein_accepts_a_numpy_integer_for_n_interp():
    seed, non_seed = stats_samples()["wass_10x10_plain"]
    assert wasserstein_areas(seed, non_seed, n_interp=np.int64(100)) == wasserstein_areas(
        seed, non_seed, n_interp=100
    )


def test_wasserstein_rejects_nan_and_an_infinite_value():
    ten = [float(k) for k in range(10)]
    with pytest.raises(InputError, match=r"^seed\[2\] is NaN"):
        wasserstein_areas([0.0, 1.0, np.nan, *ten], ten)
    with pytest.raises(InputError, match=r"^non_seed\[0\] -inf is infinite; expected a finite"):
        wasserstein_areas(ten, [-np.inf, *ten])


EXPECTED_RANGE = r"expected a range of at most 8\.988465674311579e\+307$"


def test_wasserstein_rejects_a_range_that_a_float_cannot_hold_D038():
    ten = [float(k) for k in range(10)]
    message = rf"^seed and non_seed span from -1e\+308 to 1e\+308; {EXPECTED_RANGE}"
    with pytest.raises(InputError, match=message):
        wasserstein_areas([-1e308, *ten], [1e308, *ten])


def test_wasserstein_rejects_a_finite_range_whose_sum_overflows_D038():
    # The range is 1.32e308, below the largest float 1.797e308, but the sum of the products of the
    # trapezoid rule is twice the area, and it is 2.4e308.
    seed = np.linspace(-6e307, -6.6e307, 10)
    non_seed = np.linspace(6e307, 6.6e307, 10)
    message = rf"^seed and non_seed span from -6.6e\+307 to 6.6e\+307; {EXPECTED_RANGE}"
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(InputError, match=message):
            wasserstein_areas(seed, non_seed, n_interp=100)


def test_wasserstein_accepts_a_range_just_below_the_limit_D038():
    # The range is 8.8e307. The sum of the products is at most 1.76e308, below the largest float.
    seed = np.linspace(-4.4e307, -4.0e307, 10)
    non_seed = np.linspace(4.0e307, 4.4e307, 10)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        areas = wasserstein_areas(seed, non_seed, n_interp=100)
    assert math.isfinite(areas.absolute)
    assert areas.absolute == areas.negative_censored > 8e307
    assert areas.positive_censored == 0.0


def test_wasserstein_rejects_an_empty_sample():
    ten = [float(k) for k in range(10)]
    with pytest.raises(InputError, match=r"^seed is empty; expected at least 1 value$"):
        wasserstein_areas([], ten)
    with pytest.raises(InputError, match=r"^non_seed is empty; expected at least 1 value$"):
        wasserstein_areas(ten, [])


CHUNKS = [1, 2, 3, 7]


def bits(areas):
    return tuple(value.hex() for value in areas)


@st.composite
def float_samples(draw):
    values = st.floats(min_value=-5, max_value=5, allow_nan=False, allow_infinity=False)
    seed = draw(st.lists(values, min_size=10, max_size=30).filter(lambda v: len(set(v)) > 1))
    non_seed = draw(st.lists(values, min_size=10, max_size=30).filter(lambda v: len(set(v)) > 1))
    return seed, non_seed, draw(st.integers(min_value=2, max_value=70))


@settings(deadline=None, max_examples=150)
@given(st.one_of(dyadic_samples(), float_samples()), st.sampled_from(CHUNKS))
def test_the_areas_are_bit_equal_for_each_chunk_size(case, chunk):
    seed, non_seed, n_interp = case
    whole = _unrounded_areas(np.array(seed), np.array(non_seed), n_interp)
    chunked = _unrounded_areas(np.array(seed), np.array(non_seed), n_interp, chunk=chunk)
    assert bits(chunked) == bits(whole)


@pytest.mark.parametrize("chunk", CHUNKS)
def test_the_areas_of_a_golden_case_are_bit_equal_for_each_chunk_size(chunk):
    seed, non_seed = stats_samples()["wass_10x50_ties"]
    args = (np.array(seed), np.array(non_seed), 10000)
    assert bits(_unrounded_areas(*args, chunk=chunk)) == bits(_unrounded_areas(*args))


def test_the_areas_are_bit_equal_for_a_grid_of_more_than_one_default_chunk():
    generator = np.random.default_rng(8)
    seed, non_seed = generator.normal(0.2, 1, 3000), np.round(generator.normal(0, 1, 4000), 2)
    for n_interp in (65536, 65537, 65538, 150_001):
        args = (seed, non_seed, n_interp)
        assert bits(_unrounded_areas(*args)) == bits(_unrounded_areas(*args, chunk=10**9))


def peak_bytes(seed, non_seed, n_interp):
    tracemalloc.start()
    try:
        wasserstein_areas(seed, non_seed, n_interp=n_interp)
        return tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


def test_wasserstein_memory_does_not_grow_with_n_interp():
    generator = np.random.default_rng(6)
    seed, non_seed = generator.normal(0.2, 1, 20000), generator.normal(0, 1, 20000)
    # One chunk of the grid holds a few arrays of 0.5 MB. The whole grid of 2 million points would
    # hold 16 MB for each such array.
    few, many = peak_bytes(seed, non_seed, 200_000), peak_bytes(seed, non_seed, 2_000_000)
    assert many <= 1.2 * few
    assert many < 10e6


def test_wasserstein_runs_on_a_large_grid_and_large_samples():
    generator = np.random.default_rng(5)
    seed, non_seed = generator.normal(0.2, 1, size=2000), generator.normal(0, 1, size=18000)
    start = time.perf_counter()
    areas = wasserstein_areas(seed, non_seed, n_interp=1_000_000)
    elapsed = time.perf_counter() - start
    assert elapsed < 10
    coarse = wasserstein_areas(seed, non_seed, n_interp=10_000)
    assert math.isclose(areas.absolute, coarse.absolute, abs_tol=0.002)


DIFFERENCE = "plots/wasserstein_difference.tsv"


@pytest.mark.parametrize("case", WASSERSTEIN_CASES)
def test_wasserstein_difference_matches_r(case):
    seed, non_seed = wasserstein_samples(case)
    got = wasserstein_difference(seed, non_seed, n_interp=100)
    table = golden(DIFFERENCE)
    expected = table[(table["case"] == case) & (table["n_interp"] == "100")]
    assert len(got) == len(expected) == 100
    # The x values of R differ from the port in the last bit, thus no comparison with ==.
    for column in ("x", "difference"):
        np.testing.assert_allclose(
            got[column], expected[column].astype(float), rtol=1e-9, atol=1e-12
        )


def test_wasserstein_difference_gives_the_areas():
    seed, non_seed = wasserstein_samples("test_data")
    got = wasserstein_difference(seed, non_seed)
    area = np.trapezoid(np.abs(got["difference"].to_numpy()), got["x"].to_numpy())
    assert abs(area - wasserstein_areas(seed, non_seed).absolute) <= 5e-4


def test_wasserstein_difference_has_the_grid_of_the_areas():
    seed, non_seed = wasserstein_samples("test_data")
    got = wasserstein_difference(seed, non_seed, n_interp=50)
    assert list(got.columns) == ["x", "seed", "non_seed", "difference"]
    assert set(got.dtypes) == {np.dtype("float64")}
    assert isinstance(got.index, pd.RangeIndex)
    assert len(got) == 50
    assert got["x"].iloc[0] == min(*seed, *non_seed)
    assert got["x"].iloc[-1] == max(*seed, *non_seed)
    assert got["seed"].iloc[0] == 1 / 18
    assert got["non_seed"].iloc[0] == 0.0
    assert got["seed"].iloc[-1] == got["non_seed"].iloc[-1] == 1.0
    np.testing.assert_array_equal(got["difference"], got["seed"] - got["non_seed"])
    assert len(wasserstein_difference(seed, non_seed)) == 10000


def test_wasserstein_difference_gives_the_two_ecdf_curves_of_the_areas():
    seed, non_seed = stats_samples()["wass_10x10_plain"]
    got = wasserstein_difference(seed, non_seed, n_interp=100)
    assert (np.diff(got["seed"]) >= 0).all()
    assert (np.diff(got["non_seed"]) >= 0).all()
    area = np.trapezoid(np.abs(got["difference"].to_numpy()), got["x"].to_numpy())
    assert abs(area - wasserstein_areas(seed, non_seed, n_interp=100).absolute) <= 5e-4


TEN = [float(k) for k in range(10)]
REJECTED = [
    pytest.param(
        {"seed": TEN[:9], "non_seed": TEN},
        "seed has 9 values; expected at least 10",
        id="nine-values",
    ),
    pytest.param(
        {"seed": TEN, "non_seed": [1.0] * 10},
        "non_seed has 10 values that are all 1.0; expected at least 2 distinct values",
        id="one-distinct-value",
    ),
    pytest.param(
        {"seed": [0.0, 1.0, np.nan, *TEN], "non_seed": TEN},
        "seed[2] is NaN; expected a number",
        id="nan",
    ),
    pytest.param(
        {"seed": TEN, "non_seed": TEN, "n_interp": 1},
        "n_interp 1 is less than 2; expected 2 or more",
        id="one-point",
    ),
    pytest.param(
        {"seed": TEN, "non_seed": TEN, "n_interp": True},
        "n_interp True has type bool; expected an integer",
        id="bool-points",
    ),
    pytest.param(
        {"seed": [-1e308, *TEN], "non_seed": [1e308, *TEN]},
        "seed and non_seed span from -1e+308 to 1e+308; "
        "expected a range of at most 8.988465674311579e+307",
        id="range",
    ),
]


@pytest.mark.parametrize(("kwargs", "message"), REJECTED)
def test_wasserstein_difference_rejects_what_the_areas_reject(kwargs, message):
    with pytest.raises(InputError, match=exactly(message)):
        wasserstein_areas(**kwargs)
    with pytest.raises(InputError, match=exactly(message)):
        wasserstein_difference(**kwargs)


def test_wasserstein_difference_accepts_a_series_a_tuple_and_a_numpy_integer():
    seed, non_seed = stats_samples()["wass_10x10_plain"]
    expected = wasserstein_difference(seed, non_seed, n_interp=100)
    got = wasserstein_difference(pd.Series(seed), tuple(non_seed), n_interp=np.int64(100))
    pd.testing.assert_frame_equal(got, expected)


def test_wasserstein_difference_does_not_change_its_inputs():
    seed, non_seed = np.arange(10.0)[::-1].copy(), np.arange(12.0)[::-1] + 0.5
    wasserstein_difference(seed, non_seed, n_interp=50)
    assert seed.tolist() == list(range(9, -1, -1))
    assert non_seed.tolist() == [11.5 - k for k in range(12)]


def test_wasserstein_difference_holds_a_few_times_its_result_in_memory():
    generator = np.random.default_rng(7)
    seed, non_seed = generator.normal(0.2, 1, 2000), generator.normal(0, 1, 20000)
    n_interp = 200_000
    tracemalloc.start()
    try:
        wasserstein_difference(seed, non_seed, n_interp=n_interp)
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    assert peak < 4 * 32 * n_interp
