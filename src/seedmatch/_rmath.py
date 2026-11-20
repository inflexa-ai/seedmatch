"""Arithmetic of R that Python does not share."""

import math
from typing import Final

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
