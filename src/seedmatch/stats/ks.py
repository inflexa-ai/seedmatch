"""The two-sample Kolmogorov-Smirnov test."""

import math
from typing import Final

import numpy as np
import numpy.typing as npt
from numpy.typing import ArrayLike
from scipy.special import kolmogorov

from seedmatch.stats.types import (
    Alternative,
    PValueMethod,
    TwoSampleResult,
    as_alternative,
    as_sample,
)

_EXACT_PAIRS_LIMIT: Final = 10_000


def ks_test(
    x: ArrayLike, y: ArrayLike, *, alternative: Alternative = Alternative.TWO_SIDED
) -> TwoSampleResult:
    """Do the two-sample Kolmogorov-Smirnov test, as `stats::ks.test` of R does.

    The statistic is D, `D+` for `GREATER`, or `D-` for `LESS`, from the ECDFs F of `x` and G of
    `y`. `GREATER` tests F above G. An infinite value is the greatest or the least value. The
    p-value is exact if the product of the two sizes is below 10000, also with ties, and
    asymptotic otherwise (`D-032`). The test gives no warning for ties (`D-037`).
    """
    chosen = as_alternative(alternative)
    sample_x = as_sample("x", x, finite=False)
    sample_y = as_sample("y", y, finite=False)
    m, n = sample_x.size, sample_y.size
    gap, checkpoint = _lattice(sample_x, sample_y)
    d = int(_bound(gap[checkpoint], chosen).max())
    method = PValueMethod.EXACT if m * n < _EXACT_PAIRS_LIMIT else PValueMethod.ASYMPTOTIC
    if d < 1:
        p_value = 1.0
    elif method is PValueMethod.EXACT:
        p_value = _exact_p_value(m, n, d, checkpoint, chosen)
    else:
        p_value = _asymptotic_p_value(m, n, d, chosen)
    return TwoSampleResult(
        statistic=d / (m * n),
        p_value=min(max(p_value, 0.0), 1.0),
        method=method,
        alternative=chosen,
        n_x=m,
        n_y=n,
    )


def _lattice(
    x: npt.NDArray[np.float64], y: npt.NDArray[np.float64]
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.bool_]]:
    """Return `i * n - j * m` and the checkpoint flag at each position of the pooled order.

    Entry k - 1 describes the position k, where i and j count the values of `x` and `y` up to k. A
    position is a checkpoint if it ends a run of equal values, because an ECDF jumps only there.
    The integer gap is exact, so the statistic d / (m * n) has no rounding of a running sum.
    """
    m, n = x.size, y.size
    pooled = np.concatenate((x, y))
    order = np.argsort(pooled)
    i = np.cumsum(order < m, dtype=np.int64)
    j = np.arange(1, m + n + 1, dtype=np.int64) - i
    values = pooled[order]
    return i * n - j * m, np.append(values[1:] > values[:-1], True)


def _bound(gap: npt.NDArray[np.int64], alternative: Alternative) -> npt.NDArray[np.int64]:
    match alternative:
        case Alternative.GREATER:
            return gap
        case Alternative.LESS:
            return -gap
        case Alternative.TWO_SIDED:
            return np.abs(gap)


def _exact_p_value(
    m: int, n: int, d: int, checkpoint: npt.NDArray[np.bool_], alternative: Alternative
) -> float:
    """Return P(D >= d / (m * n)) when all choices of the positions of `x` are equally likely.

    A path from (0, 0) to (m, n) fails at a checkpoint (i, j) if the bound of the alternative is at
    least d there. C(i, j) is the probability that a path to (i, j) has failed. It is 1 at a point
    that fails, and else the average of its two predecessors, weighted by i and j (Viehmann, 2021).
    The tail itself has no cancellation, thus a small p-value keeps its relative precision.
    """
    i = np.arange(m + 1)[:, np.newaxis]
    j = np.arange(n + 1)[np.newaxis, :]
    on_checkpoint = np.concatenate(([False], checkpoint))[i + j]
    fails = (on_checkpoint & (_bound(i * n - j * m, alternative) >= d)).tolist()
    previous = [0.0] * (n + 1)
    for row in range(m + 1):
        current = [0.0] * (n + 1)
        for column in range(n + 1):
            if fails[row][column]:
                current[column] = 1.0
            elif row + column:
                left = current[column - 1] if column else 0.0
                current[column] = (row * previous[column] + column * left) / (row + column)
        previous = current
    return previous[n]


def _asymptotic_p_value(m: int, n: int, d: int, alternative: Alternative) -> float:
    en = m * n / (m + n)
    statistic = d / (m * n)
    if alternative is Alternative.TWO_SIDED:
        return float(kolmogorov(math.sqrt(en) * statistic))
    return math.exp(-2 * en * statistic**2)
