import math
from fractions import Fraction

import numpy as np
import pytest
from helpers import golden_rows
from hypothesis import given, settings
from hypothesis import strategies as st

from seedmatch._rmath import r_quantile, r_round
from seedmatch.errors import InputError


def test_r_round_matches_r():
    mismatches = [
        (row.x, row.digits, row.r, r_round(float(row.x), int(row.digits)))
        for row in golden_rows("stats/round.tsv")
        if r_round(float(row.x), int(row.digits)) != float(row.r)
    ]
    assert mismatches == []


def test_r_round_differs_from_python_round():
    assert r_round(0.0025, 3) == 0.002
    assert round(0.0025, 3) == 0.003
    assert r_round(0.1235, 3) == 0.124


@pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf, 0.0, -0.0])
def test_r_round_keeps_nan_infinity_and_zero(value):
    for digits in (0, 3):
        result = r_round(value, digits)
        assert math.copysign(1.0, result) == math.copysign(1.0, value)
        assert result == value or (math.isnan(result) and math.isnan(value))


@pytest.mark.parametrize(
    ("value", "expected"),
    [(0.5, 0.0), (1.5, 2.0), (2.5, 2.0), (3.5, 4.0), (-0.5, -0.0), (-1.5, -2.0), (-2.5, -2.0)],
)
def test_r_round_with_0_digits_rounds_a_tie_to_the_even_integer(value, expected):
    result = r_round(value, 0)
    assert result == expected
    assert math.copysign(1.0, result) == math.copysign(1.0, expected)


def test_r_round_with_0_digits_keeps_the_sign_of_a_small_negative_value():
    assert math.copysign(1.0, r_round(-0.3, 0)) == -1.0


def test_r_round_gives_a_plain_float():
    result = r_round(np.float64(0.1235), 3)
    assert type(result) is float
    assert type(r_round(2.5, 0)) is float


@pytest.mark.parametrize(
    ("value", "digits"),
    [(1234567890.12345, 6), (123456789012.3456, 6), (2.675, 15), (0.1 + 0.2, 22), (1e300, 3)],
)
def test_r_round_keeps_a_value_that_needs_no_rounding(value, digits):
    # R 4.5.3 gives the value back: it has no more than 15 significant digits to round away.
    assert r_round(value, digits) == value


def test_r_round_rounds_a_value_that_has_digits_to_remove_even_for_22_digits():
    assert r_round(1.2e-22, 22) == 1e-22
    assert r_round(-1.2e-22, 22) == -1e-22


@pytest.mark.parametrize("digits", [-1, 23, 1000, -(10**30)])
def test_r_round_rejects_digits_outside_0_to_22(digits):
    with pytest.raises(
        InputError, match=r"^digits .* is out of range; expected an integer from 0 to 22"
    ):
        r_round(1.5, digits)


@pytest.mark.parametrize("digits", [True, 3.0, "3", None])
def test_r_round_rejects_digits_that_are_not_an_integer(digits):
    with pytest.raises(InputError, match=r"^digits .* has type \w+; expected an integer"):
        r_round(1.5, digits)


def oracle(x: float, digits: int) -> float:
    """Return the nearer of xd and xu, with exact distances and the float xd and xu of R.

    R computes x * 10^digits in float64, then takes the floor and the ceiling of it, then divides
    each by 10^digits in float64. A tie goes to the even floor. This code reads no value from the
    code under test.
    """
    if x < 0:
        return -oracle(-x, digits)
    p = 10**digits
    scaled = float(Fraction(x) * p)
    below, above = math.floor(scaled), math.ceil(scaled)
    lower = float(Fraction(below, p))
    upper = float(Fraction(above, p))
    d_lower = Fraction(x) - Fraction(lower)
    d_upper = Fraction(upper) - Fraction(x)
    if d_upper < d_lower or (d_upper == d_lower and below % 2 == 1):
        return upper
    return lower


@st.composite
def rounding_cases(draw):
    digits = draw(st.integers(min_value=1, max_value=6))
    if draw(st.booleans()):
        # A midpoint between two decimals, moved by up to 3 steps of the float spacing.
        k = draw(st.integers(min_value=0, max_value=10**8))
        x = float(Fraction(2 * k + 1, 2 * 10**digits))
        steps = draw(st.integers(min_value=-3, max_value=3))
        x = math.nextafter(x, math.inf if steps > 0 else -math.inf, steps=abs(steps))
    else:
        x = draw(st.floats(min_value=0, max_value=1e8, allow_nan=False, allow_infinity=False))
    return (-x if draw(st.booleans()) else x), digits


@settings(deadline=None, max_examples=2000)
@given(rounding_cases())
def test_r_round_gives_the_nearer_of_the_floor_and_the_ceiling(case):
    x, digits = case
    assert r_round(x, digits) == oracle(x, digits)


@pytest.mark.parametrize("case", golden_rows("network/r_quantiles.tsv"), ids=lambda row: row.case)
def test_r_quantile_matches_r(case):
    values = np.array([float(text) for text in case.values.split(",")])
    assert r_quantile(values, float(case.probability)) == float(case.quantile)


def test_r_quantile_differs_from_numpy_in_the_last_bit():
    (row,) = (row for row in golden_rows("network/r_quantiles.tsv") if row.case == "links_p95")
    values = np.array([float(text) for text in row.values.split(",")])
    assert r_quantile(values, 0.95) == 897.4999999999999
    assert np.quantile(values, 0.95) == 897.5


def test_r_quantile_keeps_a_repeated_value_exactly():
    value = 7.188239240658031
    assert r_quantile(np.full(27, value), 0.7141294836112025) == value


def test_r_quantile_accepts_integer_values_and_keeps_its_input():
    values = np.array([800, 400, 900, 400, 650], dtype=np.int16)
    before = values.copy()
    assert r_quantile(values, 0.5) == 650.0
    assert type(r_quantile(values, 0.5)) is float
    assert values.tolist() == before.tolist()


@pytest.mark.parametrize("probability", [0.0, 0.1, 0.25, 0.5, 0.73, 0.95, 1.0])
@pytest.mark.parametrize("size", [1, 2, 3, 10, 101])
def test_r_quantile_agrees_with_numpy_to_rounding(size, probability):
    values = np.random.default_rng(size).normal(size=size)
    assert r_quantile(values, probability) == pytest.approx(
        float(np.quantile(values, probability)), rel=1e-12, abs=1e-12
    )
