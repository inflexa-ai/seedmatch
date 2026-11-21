import math
import tracemalloc
import warnings

import numpy as np
import pytest
from helpers import WASSERSTEIN_CASES, exactly, golden_note, wasserstein_samples

from seedmatch.errors import InputError
from seedmatch.stats import density
from seedmatch.stats.density import DENSITY_POINTS, gaussian_density, nrd0_bandwidth

DENSITY = "plots/wasserstein_density.tsv"


@pytest.mark.parametrize("group", ["non-seed", "seed"])
@pytest.mark.parametrize("case", WASSERSTEIN_CASES)
def test_nrd0_bandwidth_matches_r(case, group):
    seed, non_seed = wasserstein_samples(case)
    sample = np.array(seed if group == "seed" else non_seed)
    expected = float(golden_note(DENSITY, rf"^# bandwidth {case} {group}: (.*)$"))
    assert math.isclose(nrd0_bandwidth(group, sample), expected, rel_tol=1e-15, abs_tol=0)


@pytest.mark.parametrize(
    ("values", "spread"),
    [
        pytest.param([0.0] * 8 + [-1.0, 1.0], math.sqrt(2 / 9), id="zero-iqr-uses-the-sd"),
        pytest.param(
            [-100.0, -1.0, -0.5, 0.0, 0.5, 1.0, 100.0],
            (0.75 - (-0.75)) / 1.34,
            id="narrow-iqr-uses-the-iqr",
        ),
        pytest.param([0.0] * 5 + [1.0] * 5, math.sqrt(2.5 / 9), id="smaller-sd-is-used"),
    ],
)
def test_nrd0_bandwidth_takes_the_smaller_spread(values, spread):
    sample = np.array(values)
    assert math.isclose(
        nrd0_bandwidth("seed", sample), 0.9 * spread * sample.size**-0.2, rel_tol=1e-12
    )


@pytest.mark.parametrize(
    "sample",
    [np.arange(10) * 5e-324, np.arange(10) * 1e-300],
    ids=["denormal-spread", "spread-below-the-smallest-normal-square"],
)
def test_nrd0_bandwidth_rejects_a_bandwidth_of_zero(sample):
    message = (
        "the bandwidth of the seed sample is 0.0; expected a finite number above 0, "
        "because the values are too close together"
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(InputError, match=exactly(message)):
            nrd0_bandwidth("seed", sample)


def test_nrd0_bandwidth_rejects_a_bandwidth_that_is_not_finite():
    # The sum of squares of the deviations overflows, and the interquartile range is 0.
    sample = np.array([0.0] * 8 + [4e307, -4e307])
    message = (
        "the bandwidth of the non_seed sample is inf; expected a finite number above 0, "
        "because the values are too far apart"
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(InputError, match=exactly(message)):
            nrd0_bandwidth("non_seed", sample)


def test_nrd0_bandwidth_rejects_a_bandwidth_that_is_nan():
    # Two partial sums of the mean overflow to inf and -inf, and their sum is NaN.
    sample = np.zeros(40)
    sample[0:40:8] = 4.4e307
    sample[1:40:8] = -4.4e307
    message = (
        "the bandwidth of the seed sample is nan; expected a finite number above 0, "
        "because the values are too far apart"
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(InputError, match=exactly(message)):
            nrd0_bandwidth("seed", sample)


def exact_density(values, grid, bandwidth):
    return [
        math.fsum(math.exp(-0.5 * ((point - value) / bandwidth) ** 2) for value in values)
        / (len(values) * bandwidth * math.sqrt(2 * math.pi))
        for point in grid
    ]


def test_gaussian_density_is_the_exact_sum(monkeypatch):
    values, grid = np.array([0.0, 1.0, 3.0]), np.array([-1.0, 0.0, 2.0])
    expected = exact_density(values, grid, 0.5)
    np.testing.assert_allclose(gaussian_density(values, grid, 0.5), expected, rtol=0, atol=1e-15)
    monkeypatch.setattr(density, "_BLOCK", 2)
    np.testing.assert_allclose(gaussian_density(values, grid, 0.5), expected, rtol=0, atol=1e-15)


def test_gaussian_density_is_the_exact_sum_for_more_values_than_a_block():
    values = np.random.default_rng(3).normal(size=2500)
    grid = np.linspace(-4.0, 4.0, 37)
    got = gaussian_density(values, grid, 0.3)
    np.testing.assert_allclose(got, exact_density(values, grid, 0.3), rtol=1e-13, atol=0)
    assert got.dtype == np.float64
    assert got.shape == grid.shape


def test_gaussian_density_gives_a_far_value_the_weight_zero():
    values, grid = np.array([0.0, 1.0, 1e300]), np.array([0.0])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        got = gaussian_density(values, grid, 1.0)
    expected = (1.0 + math.exp(-0.5)) / (3 * math.sqrt(2 * math.pi))
    np.testing.assert_allclose(got, [expected], rtol=1e-15, atol=0)


def test_gaussian_density_integrates_to_one():
    values = np.random.default_rng(4).normal(size=500)
    grid = np.linspace(-8.0, 8.0, 4001)
    area = np.trapezoid(gaussian_density(values, grid, 0.4), grid)
    assert math.isclose(area, 1.0, abs_tol=1e-9)


def test_gaussian_density_does_not_change_its_inputs():
    values, grid = np.array([3.0, 0.0, 1.0]), np.array([2.0, -1.0, 0.0])
    gaussian_density(values, grid, 0.5)
    assert values.tolist() == [3.0, 0.0, 1.0]
    assert grid.tolist() == [2.0, -1.0, 0.0]


def peak_bytes(values, grid):
    tracemalloc.start()
    try:
        gaussian_density(values, grid, 0.3)
        return tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


def test_gaussian_density_work_memory_is_one_block_of_the_grid(monkeypatch):
    values = np.random.default_rng(5).normal(size=5000)
    grid = np.linspace(-4.0, 4.0, DENSITY_POINTS)
    block_bytes = grid.size * 1024 * 8
    default = peak_bytes(values, grid)
    assert block_bytes <= default < 1.25 * block_bytes
    monkeypatch.setattr(density, "_BLOCK", 16)
    assert peak_bytes(values, grid) < default / 10
