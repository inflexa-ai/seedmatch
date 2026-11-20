"""The two-sample Wilcoxon rank-sum test."""

import math
from typing import Final, assert_never

import numpy as np
import numpy.typing as npt
from numpy.typing import ArrayLike
from scipy.special import log_ndtr

from seedmatch.stats.types import (
    Alternative,
    PValueMethod,
    TwoSampleResult,
    as_alternative,
    as_sample,
)

_EXACT_SIZE_LIMIT: Final = 50


def wilcoxon_test(
    x: ArrayLike, y: ArrayLike, *, alternative: Alternative = Alternative.TWO_SIDED
) -> TwoSampleResult:
    """Do the two-sample Wilcoxon rank-sum test, as `stats::wilcox.test` of R does.

    The statistic W is the rank sum of `x` minus n_x (n_x + 1) / 2, with the average rank for tied
    values. `GREATER` tests `x` above `y`. An infinite value is the greatest or the least value.
    The p-value is exact if both sizes are below 50 and no two values are equal. Otherwise, it comes
    from the normal approximation with the tie correction and the continuity correction. If each
    value is equal, the p-value is 1 (`D-033`). The test gives no warning for ties (`D-037`).
    """
    chosen = as_alternative(alternative)
    sample_x = as_sample("x", x, finite=False)
    sample_y = as_sample("y", y, finite=False)
    m, n = sample_x.size, sample_y.size
    statistic, tie_sizes = _rank_statistic(sample_x, sample_y)
    no_ties = tie_sizes.size == m + n
    if m < _EXACT_SIZE_LIMIT and n < _EXACT_SIZE_LIMIT and no_ties:
        method = PValueMethod.EXACT
        p_value = _exact_p_value(m, n, int(statistic), chosen)
    else:
        method = PValueMethod.ASYMPTOTIC
        p_value = _asymptotic_p_value(m, n, statistic, tie_sizes, chosen)
    return TwoSampleResult(
        statistic=statistic,
        p_value=min(max(p_value, 0.0), 1.0),
        method=method,
        alternative=chosen,
        n_x=m,
        n_y=n,
    )


def _rank_statistic(
    x: npt.NDArray[np.float64], y: npt.NDArray[np.float64]
) -> tuple[float, npt.NDArray[np.int64]]:
    """Return W and the sizes of the groups of equal pooled values."""
    m = x.size
    _, group, sizes = np.unique(np.concatenate((x, y)), return_inverse=True, return_counts=True)
    average_rank = np.cumsum(sizes) - (sizes - 1) / 2
    return float(average_rank[group[:m]].sum()) - m * (m + 1) / 2, sizes


def _rank_sum_counts(m: int, n: int) -> npt.NDArray[np.float64]:
    """Return the number of choices of the positions of `x` that give W = 0, ..., m * n.

    Let ways[k, s] count the subsets of k ranks with the sum s, for the ranks 0, ..., i - 1. The
    rank i adds ways[k - 1, s - i] to ways[k, s]. A subset of m ranks has W = s - m (m - 1) / 2. A
    sum above the greatest W cannot lead to a smaller one, thus the table stops there. The counts
    pass 2^63, thus they are float64. Each step adds non-negative numbers, thus the relative error
    stays near 1e-15.
    """
    offset = m * (m - 1) // 2
    top = offset + m * n
    ways = np.zeros((m + 1, top + 1))
    ways[0, 0] = 1.0
    for rank in range(m + n):
        ways[1:, rank:] += ways[:-1, : top + 1 - rank]
    return ways[m, offset:]


def _exact_p_value(m: int, n: int, w: int, alternative: Alternative) -> float:
    """Return a tail of W when all choices of the positions of `x` are equally likely.

    Each tail is a sum of counts, not 1 minus the other tail, thus a small p-value keeps its
    relative precision.
    """
    counts = _rank_sum_counts(m, n)
    total = math.comb(m + n, m)
    at_most = float(counts[: w + 1].sum()) / total
    at_least = float(counts[w:].sum()) / total
    match alternative:
        case Alternative.LESS:
            return at_most
        case Alternative.GREATER:
            return at_least
        case Alternative.TWO_SIDED:
            return min(1.0, 2 * (at_least if 2 * w > m * n else at_most))
        case _:
            assert_never(alternative)


def _asymptotic_p_value(
    m: int, n: int, w: float, tie_sizes: npt.NDArray[np.int64], alternative: Alternative
) -> float:
    pooled = m + n
    if tie_sizes.size == 1:
        # One group of ties makes the variance 0 (D-033).
        return 1.0
    ties = tie_sizes.astype(np.float64)
    variance = m * n / 12 * (pooled + 1 - float((ties**3 - ties).sum()) / (pooled * (pooled - 1)))
    z = w - m * n / 2
    match alternative:
        case Alternative.LESS:
            correction = -0.5
        case Alternative.GREATER:
            correction = 0.5
        case Alternative.TWO_SIDED:
            correction = 0.5 * ((z > 0) - (z < 0))
        case _:
            assert_never(alternative)
    z = (z - correction) / math.sqrt(variance)
    match alternative:
        case Alternative.LESS:
            return _normal_lower_tail(z)
        case Alternative.GREATER:
            return _normal_lower_tail(-z)
        case Alternative.TWO_SIDED:
            return 2 * min(_normal_lower_tail(z), _normal_lower_tail(-z))
        case _:
            assert_never(alternative)


def _normal_lower_tail(z: float) -> float:
    # ndtr(z) is 0.0 below about -37.5, but exp(log_ndtr(z)) keeps the subnormal floats down to
    # 5e-324, as pnorm of R does.
    return math.exp(float(log_ndtr(z)))
