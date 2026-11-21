"""The Gaussian density of a sample, as the density plot of `wass_dist` draws it."""

import math
from typing import Final

import numpy as np
import numpy.typing as npt

from seedmatch._checks import shown
from seedmatch._rmath import r_quantile
from seedmatch.errors import InputError

DENSITY_POINTS: Final = 512
_BLOCK: Final = 1024
# `bw.nrd0` rounds the interquartile range of the standard normal, 1.349, down to 1.34.
_IQR_DIVISOR: Final = 1.34


def nrd0_bandwidth(name: str, values: npt.NDArray[np.float64]) -> float:
    """Return the bandwidth of `bw.nrd0` of R: `0.9 * min(sd, iqr / 1.34) * n ** -0.2`.

    `sd` is the standard deviation with n - 1 in the divisor, and `iqr` is the interquartile range
    of the type 7 quantiles of R. If `iqr` is 0, the rule uses `sd`. `values` is a
    one-dimensional array of finite numbers with at least 2 distinct values. `name` is the name of
    the sample in the message of an error. The caller makes sure of these conditions.

    A result that is 0 or not a finite number raises `InputError`. The result is 0 if the square
    of the spread underflows. It is `inf` or NaN if the sums of squares of `np.std` overflow and
    `sd` sets the result, because the interquartile range is 0 or `sd` is NaN. If `sd` is `inf` and
    the interquartile range term is finite, that term wins `min` and the call succeeds. The time
    is O(n).
    """
    # The sums of squares can overflow. The check of the result decides, as the docstring says.
    with np.errstate(over="ignore", invalid="ignore"):
        sd = float(np.std(values, ddof=1))
    iqr = r_quantile(values, 0.75) - r_quantile(values, 0.25)
    spread = min(sd, iqr / _IQR_DIVISOR) or sd
    bandwidth = 0.9 * spread * math.pow(values.size, -0.2)
    if math.isfinite(bandwidth) and bandwidth > 0:
        return bandwidth
    reason = "close together" if bandwidth == 0 else "far apart"
    message = (
        f"the bandwidth of the {name} sample is {shown(bandwidth)}; "
        f"expected a finite number above 0, because the values are too {reason}"
    )
    raise InputError(message)


def gaussian_density(
    values: npt.NDArray[np.float64], grid: npt.NDArray[np.float64], bandwidth: float
) -> npt.NDArray[np.float64]:
    """Return the Gaussian kernel density of `values` at each point of `grid`.

    For a grid point g, the density is the sum over the values v of
    `exp(-0.5 * ((g - v) / h) ** 2)`, divided by `n * h * sqrt(2 * pi)`, with the bandwidth h.
    This is the exact sum. `density()` of R bins the values and uses an FFT, thus it differs by a
    fraction of the maximum.

    The function reads the values in blocks of `_BLOCK` = 1024, and it reads `_BLOCK` at each call.
    Thus the work memory is `grid.size * 1024 * 8` bytes, 4 MiB for 512 points. The time is
    O(`grid.size` * n) for n values.
    """
    total = np.zeros(grid.size)
    buffer = np.empty((grid.size, min(_BLOCK, values.size)))
    # A distance that overflows to infinity gives the weight exp(-inf) = 0, the correct weight.
    with np.errstate(over="ignore"):
        for start in range(0, values.size, _BLOCK):
            block = values[start : start + _BLOCK]
            work = buffer[:, : block.size]
            np.subtract(grid[:, np.newaxis], block, out=work)
            work /= bandwidth
            np.square(work, out=work)
            work *= -0.5
            np.exp(work, out=work)
            total += work.sum(axis=1)
    total /= values.size * bandwidth * math.sqrt(2 * math.pi)
    return total
