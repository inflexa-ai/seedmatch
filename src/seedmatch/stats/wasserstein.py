"""The Wasserstein areas between the ECDF of a seed sample and the ECDF of a non-seed sample."""

import sys
from dataclasses import dataclass
from typing import Final, NamedTuple

import numpy as np
import numpy.typing as npt
from numpy.typing import ArrayLike

from seedmatch._checks import as_limit, check_at_least, shown
from seedmatch._rmath import r_round
from seedmatch.errors import InputError
from seedmatch.stats.types import as_sample, check_range

_MIN_VALUES: Final = 10
_DECIMALS: Final = 3
_GRID_CHUNK: Final = 2**16
# Each product of the trapezoid rule is a width times a sum of two heights in [0, 1], so a sum of
# products is at most twice the range of the pooled values. The sum overflows above this range.
_MAX_RANGE: Final = sys.float_info.max / 2


@dataclass(frozen=True)
class WassersteinAreas:
    """The areas of `wasserstein_areas`, each rounded to 3 decimals with the rule of R.

    Let diff be the ECDF of the seed sample minus the ECDF of the non-seed sample. `absolute` is
    the area of |diff|. `negative_censored` is the area of diff where it is positive, and
    `positive_censored` is the area of -diff where diff is negative. The names follow R, which sets
    the negative or the positive values of diff to 0 before it integrates.
    """

    absolute: float
    negative_censored: float
    positive_censored: float


def wasserstein_areas(
    seed: ArrayLike, non_seed: ArrayLike, *, n_interp: int = 10000
) -> WassersteinAreas:
    """Return the areas between two ECDFs, as `wass_dist` of R does (`R/wass_dist.R:61-81`).

    Each ECDF is evaluated at its own values, and equal values merge to one point, as in
    `stats::approx`. The function gives no warning for equal values (`D-037`). The linear
    interpolation of each ECDF is read at `n_interp` points from the minimum to the maximum of the
    pooled values, with 0 below the first value of the sample and 1 above its last value. The
    trapezoid rule integrates the three differences.

    Each sample needs at least 10 values and at least 2 distinct values. NaN and an infinite value
    raise `InputError`, as the failure of R does. `n_interp` is at least 2, and the pooled range is
    at most half of the largest float (`D-038`). The work is O(n log n + n_interp log n), because
    `np.unique` sorts each sample and `np.searchsorted` searches for each grid point. The memory
    does not grow with `n_interp`, because the grid goes through in chunks.
    """
    seed_values = _as_wasserstein_sample("seed", seed)
    non_seed_values = _as_wasserstein_sample("non_seed", non_seed)
    check_wasserstein_range("seed", "non_seed", seed_values, non_seed_values)
    points = as_limit("n_interp", n_interp)
    check_at_least("n_interp", points, 2)
    absolute, negative, positive = _unrounded_areas(seed_values, non_seed_values, points)
    return WassersteinAreas(
        absolute=r_round(absolute, _DECIMALS),
        negative_censored=r_round(negative, _DECIMALS),
        positive_censored=r_round(positive, _DECIMALS),
    )


def _as_wasserstein_sample(name: str, values: ArrayLike) -> npt.NDArray[np.float64]:
    sample = as_sample(name, values, finite=True)
    check_wasserstein_sample(name, sample)
    return sample


def check_wasserstein_sample(name: str, sample: npt.NDArray[np.float64]) -> None:
    """Raise `InputError` unless the sample has 10 values or more, and 2 distinct values or more."""
    if sample.size < _MIN_VALUES:
        message = f"{name} has {sample.size} values; expected at least {_MIN_VALUES}"
        raise InputError(message)
    if sample.min() == sample.max():
        message = (
            f"{name} has {sample.size} values that are all {shown(float(sample[0]))}; "
            "expected at least 2 distinct values"
        )
        raise InputError(message)


def check_wasserstein_range(
    first: str,
    second: str,
    x: npt.NDArray[np.float64],
    y: npt.NDArray[np.float64],
) -> None:
    """Raise `InputError` if the pooled range is above half of the largest float (`D-038`)."""
    check_range(first, second, x, y, _MAX_RANGE)


def _unrounded_areas(
    seed: npt.NDArray[np.float64],
    non_seed: npt.NDArray[np.float64],
    n_interp: int,
    chunk: int = _GRID_CHUNK,
) -> tuple[float, float, float]:
    """Return the absolute, negative-censored, and positive-censored areas before rounding.

    The grid goes through in chunks of `chunk` trapezoids, so the memory does not grow with
    `n_interp`. A chunk reads its first point again from the end of the previous chunk. Each point
    comes from its own index, thus a chunk does not change a point or a product.
    """
    low = float(min(seed.min(), non_seed.min()))
    high = float(max(seed.max(), non_seed.max()))
    seed_ecdf, non_seed_ecdf = _ecdf(seed), _ecdf(non_seed)
    totals = [0.0, 0.0, 0.0]
    for start in range(0, n_interp - 1, chunk):
        stop = min(start + chunk, n_interp - 1)
        grid = _grid(low, high, n_interp, start, stop + 1)
        widths = np.diff(grid)
        diff = _interpolate(seed_ecdf, grid)
        diff -= _interpolate(non_seed_ecdf, grid)
        heights = (
            np.abs(diff),
            np.where(diff < 0, 0.0, diff),
            np.where(diff > 0, 0.0, np.abs(diff)),
        )
        totals = [_add_products(total, widths, y) for total, y in zip(totals, heights, strict=True)]
    absolute, negative, positive = (total / 2 for total in totals)
    return absolute, negative, positive


def _grid(low: float, high: float, n_interp: int, first: int, stop: int) -> npt.NDArray[np.float64]:
    """Return the points `first` to `stop - 1` of `seq(low, high, length.out = n_interp)` of R.

    Point k is `low + k * step`, and the last point is `high`.
    """
    grid = np.arange(first, stop, dtype=np.float64)
    grid *= (high - low) / (n_interp - 1)
    grid += low
    if stop == n_interp:
        grid[-1] = high
    return grid


class _Ecdf(NamedTuple):
    """The merged points of an ECDF, and the lengths of its segments."""

    knots: npt.NDArray[np.float64]
    heights: npt.NDArray[np.float64]
    widths: npt.NDArray[np.float64]
    steps: npt.NDArray[np.float64]


def _ecdf(sample: npt.NDArray[np.float64]) -> _Ecdf:
    knots, counts = np.unique(sample, return_counts=True)
    heights = np.cumsum(counts) / sample.size
    return _Ecdf(knots, heights, np.diff(knots), np.diff(heights))


def _interpolate(ecdf: _Ecdf, grid: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Return the interpolated ECDF at each grid point.

    Segment i joins `knots[i]` and `knots[i + 1]`. A point takes the last segment that starts at or
    below it, and the expression of `approx1` of R reads the line.
    """
    knots = ecdf.knots
    segment = np.searchsorted(knots, grid, side="right")
    segment -= 1
    np.clip(segment, 0, knots.size - 2, out=segment)
    values = grid - knots[segment]
    values /= ecdf.widths[segment]
    values *= ecdf.steps[segment]
    values += ecdf.heights[segment]
    values[grid < knots[0]] = 0.0
    # R returns the height of the last knot at that knot without the expression, which can differ
    # from the height in the last bit. The height is 1.
    values[grid >= knots[-1]] = 1.0
    return values


def _add_products(
    total: float, widths: npt.NDArray[np.float64], y: npt.NDArray[np.float64]
) -> float:
    """Return `total` plus the sum of widths[i] * (y[i + 1] + y[i]), as `caTools::trapz` adds.

    R adds the products in order. `np.sum` adds in pairs, in another order, and a different last
    bit can move an area across a rounding boundary of `r_round`. The total goes into the first
    product, so that the chunks of a grid add in the order of one sum over all products.
    """
    products = y[1:] + y[:-1]
    products *= widths
    products[0] += total
    return float(np.cumsum(products)[-1])
