"""Arithmetic of R that Python does not share."""

import math
from typing import Any, Final

import numpy as np
from numpy.typing import NDArray

from seedmatch._checks import as_limit
from seedmatch.errors import InputError

# 10.0 ** digits is exact in float64 up to 10 ** 22, and no caller needs more digits.
_MAX_DIGITS: Final = 22
_SIGNIFICANT_DIGITS: Final = 15
_LOG10_2: Final = 0.301029995663981195213738894724


def r_round(x: float, digits: int) -> float:
    """Round `x` to `digits` decimals, as `round(x, digits)` of R 4.0.0 and later does.

    R takes xd = floor(x * 10^digits) / 10^digits and xu = ceil(x * 10^digits) / 10^digits in
    float64, and gives the nearer one by the float differences. If the two differences are equal,
    R gives the candidate whose integer numerator is even: floor(x * 10^digits) for xd, and
    ceil(x * 10^digits) for xu. Thus xu wins if the floor is odd. `r_round(0.0025, 3)` is 0.002,
    but Python `round(0.0025, 3)` is 0.003, because Python rounds the exact binary value. NaN, an
    infinite value, and 0 stay the same. A negative `x` is rounded as `-x`. Python and R agree for
    0 digits: a tie goes to the even integer. `digits` is an int from 0 to 22.
    """
    places = as_limit("digits", digits)
    if not 0 <= places <= _MAX_DIGITS:
        message = f"digits {places} is out of range; expected an integer from 0 to {_MAX_DIGITS}"
        raise InputError(message)
    value = float(x)
    if not math.isfinite(value) or value == 0.0:
        return value
    if places == 0:
        return math.copysign(round(value), value)
    if value < 0.0:
        return -_round_magnitude(-value, places)
    return _round_magnitude(value, places)


def _round_magnitude(x: float, places: int) -> float:
    # If the requested decimals reach past 15 significant digits, x stays the same. R takes log10(x)
    # from the binary exponent, and this estimate decides differently from log10(x) near 10^k.
    estimate = _LOG10_2 * (0.5 + (math.frexp(x)[1] - 1))
    if estimate + places > _SIGNIFICANT_DIGITS:
        return x
    scale = float(10**places)
    scaled = x * scale
    below = math.floor(scaled)
    lower = below / scale
    upper = math.ceil(scaled) / scale
    if upper - x < x - lower or (upper - x == x - lower and below % 2 == 1):
        return upper
    return lower


def r_quantile(values: NDArray[np.number[Any]], probability: float) -> float:
    """Return the quantile of `values` as `stats::quantile(values, probability)` of R gives it.

    This is type 7, the default of R. With n values, R takes index = 1 + (n - 1) * probability,
    lo = floor(index), and hi = ceiling(index). The result is the value of rank lo. If index is
    above lo and the value of rank hi is different, the result is (1 - h) * x[lo] + h * x[hi] for
    h = index - lo. The function does the same float64 operations in the same order, thus it
    differs from `numpy.quantile` in the last bit for some inputs. `numpy.quantile` gives 897.5
    and R gives 897.4999999999999 for the 22 scores of the STRING sample at 0.95.

    `values` is a one-dimensional array of finite numbers, and it has at least one value.
    `probability` is from 0 to 1. The caller makes sure of these conditions, and the function
    does not check them. The time is O(n) and the memory is 8 bytes for each value. A float64 copy
    of `values` goes to `numpy.ndarray.partition`, which puts the values of the ranks lo and hi in
    their sorted places and does not order the other values. The input stays as it is.
    """
    partitioned = np.array(values, dtype=np.float64)
    index = 1 + (len(partitioned) - 1) * float(probability)
    low_rank, high_rank = math.floor(index), math.ceil(index)
    partitioned.partition(sorted({low_rank - 1, high_rank - 1}))
    low, high = float(partitioned[low_rank - 1]), float(partitioned[high_rank - 1])
    if index > low_rank and high != low:
        weight = index - low_rank
        return (1 - weight) * low + weight * high
    return low
