"""The DTS two-sample test."""

import operator
from typing import Final

import numpy as np
import numpy.typing as npt
from numpy.typing import ArrayLike

from seedmatch._checks import as_limit, check_at_least, check_not_negative, shown
from seedmatch.errors import InputError
from seedmatch.stats.types import PValueMethod, TwoSampleResult, as_sample

_BLOCK_SIZE: Final = 64


def dts_test(
    x: ArrayLike, y: ArrayLike, *, rng: np.random.Generator | int, n_boot: int = 2000
) -> TwoSampleResult:
    """Do the DTS two-sample test, as `twosamples::dts_test` of R does with `p = 1`.

    Let z be the sorted pooled values, with n values. The statistic is the sum over the gaps
    between consecutive values of the width of the gap times `|E - F| / sqrt(2 G (1 - G) / n)`.
    E and F are the ECDFs of `x` and `y`, and G is the pooled ECDF, each at the left value of the
    gap. The p-value is the fraction of `n_boot` random choices of the positions of `x` in z whose
    statistic is not less than the observed statistic. If no choice reaches it, the p-value is
    `1 / (2 * n_boot)`. The permutations come from `rng`, a `Generator` or an int seed, so the
    p-value agrees with R only within the error of the permutations. An infinite value and an empty
    sample raise `InputError` (`D-034`).
    """
    sample_x = as_sample("x", x, finite=True)
    sample_y = as_sample("y", y, finite=True)
    generator = _as_generator(rng)
    boots = as_limit("n_boot", n_boot)
    check_at_least("n_boot", boots, 1)
    m, n = sample_x.size, sample_y.size
    pooled = np.concatenate((sample_x, sample_y))
    order = np.argsort(pooled, kind="stable")
    ends, expected, scale = _gaps(pooled[order], m, n)
    observed = _statistics((order < m)[np.newaxis, :], ends, expected, scale)[0]
    labels = np.arange(m + n) < m
    at_least = 0
    for start in range(0, boots, _BLOCK_SIZE):
        block = np.tile(labels, (min(_BLOCK_SIZE, boots - start), 1))
        generator.permuted(block, axis=1, out=block)
        at_least += int((_statistics(block, ends, expected, scale) >= observed).sum())
    return TwoSampleResult(
        statistic=float(observed),
        p_value=at_least / boots if at_least else 1 / (2 * boots),
        method=PValueMethod.PERMUTATION,
        alternative=None,
        n_x=m,
        n_y=n,
    )


def _as_generator(rng: object) -> np.random.Generator:
    if isinstance(rng, np.random.Generator):
        return rng
    if isinstance(rng, bool | np.bool_) or not isinstance(rng, int | np.integer):
        message = (
            f"rng {shown(rng)} has type {type(rng).__name__}; "
            "expected a numpy.random.Generator or an int"
        )
        raise InputError(message)
    seed = operator.index(rng)
    check_not_negative("rng", seed)
    return np.random.default_rng(seed)


def _gaps(
    values: npt.NDArray[np.float64], m: int, n: int
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.int64], npt.NDArray[np.float64]]:
    """Return `ends`, `expected`, and `scale` for the gaps of the sorted pooled values.

    A gap with width 0 adds nothing to the statistic, so only a position that ends a run of equal
    values counts. `ends[g]` is the index of the left value of gap g, and k = `ends[g] + 1` values
    lie at or below it. With c of them from `x`, E - F = (c N - k m) / (m n), where N = m + n, and
    G = k / N. Thus the term of the gap is `scale[g] * |c N - expected[g]|`, with `expected` = k m.
    """
    pooled_size = m + n
    width = np.diff(values)
    ends = np.flatnonzero(width > 0)
    below = ends + 1
    scale = width[ends] * pooled_size**1.5 / (m * n * np.sqrt(2.0 * below * (pooled_size - below)))
    return ends, below * m, scale


def _statistics(
    labels: npt.NDArray[np.bool_],
    ends: npt.NDArray[np.int64],
    expected: npt.NDArray[np.int64],
    scale: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Return the statistic of each row of `labels`, which marks the `x` positions in sorted order.

    The same labels give the same bits alone and in a block of any size, thus the comparison `>=`
    of the p-value counts a permutation that equals the observed arrangement. This holds because
    `np.take` gives a C-contiguous array, and then each row adds its terms in the same order. The
    selection `[:, ends]` gives an F-contiguous array, and the sum of such an array with 2 or more
    rows adds column by column, in another order than the sum of 1 row.
    """
    from_x = np.take(np.cumsum(labels, axis=1, dtype=np.int64), ends, axis=1)
    deviation = np.abs(from_x * labels.shape[1] - expected)
    statistics: npt.NDArray[np.float64] = (deviation * scale).sum(axis=1)
    return statistics
