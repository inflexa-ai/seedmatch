import io
import math
import time
import tracemalloc
import warnings

import numpy as np
import pandas as pd
import pytest
from helpers import WASSERSTEIN_CASES, exactly, golden_note, golden_rows, wasserstein_samples
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import to_hex

from seedmatch import wasserstein_difference
from seedmatch.errors import InputError
from seedmatch.plot import plot_wasserstein
from seedmatch.plot._style import linewidth_points

DENSITY = "plots/wasserstein_density.tsv"
ECDF = "plots/wasserstein_ecdf.tsv"
GROUPS = ["non-seed", "seed"]
RELATIVE = 1e-12
TEST_DATA_LIMITS = (-8.743981512594443, 4.002092292512264)
TEST_DATA_TITLE = "Wass(- censored) = 2.196   Wass(+ censored) = 0   Wass (total) = 2.196"


def draw(case="test_data", **options):
    seed, non_seed = wasserstein_samples(case)
    return plot_wasserstein(seed, non_seed, **options)


def r_density(case, group):
    rows = [row for row in golden_rows(DENSITY) if (row.case, row.group) == (case, group)]
    return np.array([float(row.x) for row in rows]), np.array([float(row.density) for row in rows])


def r_ecdf(case, group):
    return [row for row in golden_rows(ECDF) if (row.case, row.group) == (case, group)]


def legend_texts(axes):
    legend = axes.get_legend()
    return [text.get_text() for text in legend.get_texts()], legend.get_title().get_text()


@pytest.mark.parametrize("case", WASSERSTEIN_CASES)
def test_plot_wasserstein_density_agrees_with_r(case):
    lines = draw(case).axes[0].get_lines()
    assert len(lines) == 2
    for line, group in zip(lines, GROUPS, strict=True):
        r_x, r_y = r_density(case, group)
        x, y = line.get_xydata().T
        assert x.size == r_x.size == 512
        np.testing.assert_allclose(x, r_x, rtol=RELATIVE, atol=RELATIVE)
        # R computes the density with an FFT, and the port sums the exact kernels.
        assert np.abs(y - r_y).max() <= 0.005 * r_y.max(), group


@pytest.mark.parametrize("case", WASSERSTEIN_CASES)
def test_plot_wasserstein_draws_the_ecdf_of_the_areas_D065(case):
    ecdf_axes = draw(case).axes[1]
    x_low, x_high = ecdf_axes.get_xlim()
    lines = ecdf_axes.get_lines()
    assert len(lines) == 2
    for line, group in zip(lines, GROUPS, strict=True):
        rows = r_ecdf(case, group)
        assert (rows[0].x, rows[-1].x) == ("-Inf", "Inf")
        finite = [(float(row.x), float(row.y)) for row in rows[1:-1]]
        points = line.get_xydata()
        assert points.shape == (len(finite) + 3, 2)
        np.testing.assert_allclose(points[2:-1], finite, rtol=RELATIVE)
        assert points[0, 1] == points[1, 1] == 0.0
        assert points[-1, 1] == 1.0
        assert points[1, 0] == points[2, 0]
        assert points[0, 0] == x_low
        assert points[-1, 0] == x_high


def test_plot_wasserstein_draws_the_ecdf_of_the_seed_values_from_the_lowest_value():
    seed_line = draw().axes[1].get_lines()[1]
    points = seed_line.get_xydata()
    np.testing.assert_allclose(
        points[:3],
        [[TEST_DATA_LIMITS[0], 0.0], [-8.16461452145323, 0.0], [-8.16461452145323, 1 / 18]],
        rtol=RELATIVE,
    )
    np.testing.assert_allclose(points[-1], [TEST_DATA_LIMITS[1], 1.0], rtol=RELATIVE)


def test_plot_wasserstein_draws_the_difference():
    seed, non_seed = wasserstein_samples("test_data")
    frame = wasserstein_difference(seed, non_seed)
    (line,) = draw().axes[2].get_lines()
    np.testing.assert_array_equal(line.get_xydata(), frame[["x", "difference"]].to_numpy())
    assert to_hex(line.get_color()) == "#000000"
    assert line.get_linewidth() == pytest.approx(linewidth_points(0.5))
    (short,) = draw(n_interp=100).axes[2].get_lines()
    assert short.get_xydata().shape == (100, 2)


def test_plot_wasserstein_gives_the_r_titles_and_labels():
    figure = draw()
    density_axes, ecdf_axes, difference_axes = figure.axes
    assert figure.get_suptitle() == TEST_DATA_TITLE
    assert figure.get_suptitle() == golden_note(
        DENSITY, r"^# title test_data n_interp 10000: (.*)$"
    )
    r_low, r_high = (
        float(limit)
        for limit in golden_note(DENSITY, r"^# panel x range test_data: (.*)$").split(",")
    )
    assert (r_low, r_high) == pytest.approx(TEST_DATA_LIMITS, rel=RELATIVE)
    for axes in figure.axes:
        assert axes.get_xlim() == pytest.approx((r_low, r_high), rel=RELATIVE, abs=RELATIVE)
        assert axes.get_xlim() == pytest.approx(TEST_DATA_LIMITS, rel=RELATIVE, abs=RELATIVE)
        assert axes.get_title() == ""
    x_label = golden_note(DENSITY, r"^# density plot: .* / x (.*) / y .* / title none$")
    assert [axes.get_xlabel() for axes in figure.axes] == [x_label] * 3
    assert density_axes.get_ylabel() == golden_note(
        DENSITY, r"^# density plot: .* / y (.*) / title none$"
    )
    assert ecdf_axes.get_ylabel() == golden_note(
        DENSITY, r"^# ecdf plot: .* / y (.*) / title none$"
    )
    assert difference_axes.get_ylabel() == golden_note(DENSITY, r"^# difference plot: .* / y (.*)$")
    assert (x_label, density_axes.get_ylabel(), ecdf_axes.get_ylabel()) == (
        "log2(FC)",
        "Density",
        "eCDF",
    )
    assert difference_axes.get_ylabel() == "Seed - NonSeed eCDF"
    for axes in (density_axes, ecdf_axes):
        assert legend_texts(axes) == (["non-seed", "seed"], "group")
    assert difference_axes.get_legend() is None


def test_plot_wasserstein_fills_each_density_down_to_zero():
    density_axes = draw().axes[0]
    for fill, line in zip(density_axes.collections, density_axes.get_lines(), strict=True):
        vertices = fill.get_paths()[0].vertices
        np.testing.assert_array_equal(vertices[1:513], line.get_xydata())
        assert (vertices[513:, 1] == 0.0).all()


def test_plot_wasserstein_gives_the_r_colors():
    density_axes, ecdf_axes, _ = draw().axes
    fills = golden_note(DENSITY, r"^# density plot: .* fills non-seed (#\w+), seed #\w+ /")
    seed_fill = golden_note(DENSITY, r"^# density plot: .* fills non-seed #\w+, seed (#\w+) /")
    expected_fills = [f"{fill.lower()}40" for fill in (fills, seed_fill)]
    assert expected_fills == ["#f8766d40", "#00bfc440"]
    got_fills = [
        to_hex(collection.get_facecolor()[0], keep_alpha=True)
        for collection in density_axes.collections
    ]
    assert got_fills == expected_fills
    outline = golden_note(DENSITY, r"^# density plot: .* outline (\w+) /")
    assert [to_hex(line.get_color()) for line in density_axes.get_lines()] == [to_hex(outline)] * 2
    assert [line.get_linewidth() for line in density_axes.get_lines()] == pytest.approx(
        [linewidth_points(0.5)] * 2
    )
    non_seed_color = golden_note(DENSITY, r"^# ecdf plot: .* colours non-seed (\w+), seed \w+ /")
    seed_color = golden_note(DENSITY, r"^# ecdf plot: .* colours non-seed \w+, seed (\w+) /")
    assert [to_hex(line.get_color()) for line in ecdf_axes.get_lines()] == [
        to_hex(non_seed_color),
        to_hex(seed_color),
    ]
    assert [line.get_linewidth() for line in ecdf_axes.get_lines()] == pytest.approx(
        [linewidth_points(1.25)] * 2
    )


@pytest.mark.parametrize("n_interp", [100, 10000])
@pytest.mark.parametrize("case", WASSERSTEIN_CASES)
def test_plot_wasserstein_gives_the_r_title_and_range_of_each_case(case, n_interp):
    figure = draw(case, n_interp=n_interp)
    assert figure.get_suptitle() == golden_note(
        DENSITY, rf"^# title {case} n_interp {n_interp}: (.*)$"
    )
    r_low, r_high = (
        float(limit)
        for limit in golden_note(DENSITY, rf"^# panel x range {case}: (.*)$").split(",")
    )
    for axes in figure.axes:
        # The note of R has 15 significant digits.
        assert axes.get_xlim() == pytest.approx((r_low, r_high), rel=1e-12)


def test_plot_wasserstein_gives_the_areas_of_n_interp_100_in_its_title():
    assert draw(n_interp=100).get_suptitle() == (
        "Wass(- censored) = 2.197   Wass(+ censored) = 0   Wass (total) = 2.197"
    )


def test_plot_wasserstein_styles_the_axes_with_the_panel_of_theme_grey():
    figure = draw()
    assert [to_hex(axes.get_facecolor()).lower() for axes in figure.axes] == ["#ebebeb"] * 3
    assert not any(spine.get_visible() for axes in figure.axes for spine in axes.spines.values())


def test_plot_wasserstein_puts_each_legend_right_of_its_axes_inside_the_figure():
    figure = draw()
    FigureCanvasAgg(figure).draw()
    inside = figure.bbox
    for axes in figure.axes[:2]:
        legend = axes.get_legend().get_window_extent()
        panel = axes.get_window_extent()
        assert legend.x0 >= panel.x1
        assert inside.x0 <= legend.x0
        assert legend.x1 <= inside.x1
        assert inside.y0 <= legend.y0
        assert legend.y1 <= inside.y1


def test_plot_wasserstein_takes_the_size():
    assert tuple(draw().get_size_inches()) == (15.0, 5.0)
    assert tuple(draw(size=(9, 3)).get_size_inches()) == (9.0, 3.0)


def test_plot_wasserstein_draws_a_figure_that_saves_without_a_warning():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        draw().savefig(io.BytesIO(), format="png")


def test_plot_wasserstein_accepts_a_series_and_an_array():
    seed, non_seed = wasserstein_samples("test_data")
    expected = draw()
    got = plot_wasserstein(pd.Series(seed), np.array(non_seed))
    assert got.get_suptitle() == expected.get_suptitle()
    for axes, expected_axes in zip(got.axes, expected.axes, strict=True):
        for line, expected_line in zip(axes.get_lines(), expected_axes.get_lines(), strict=True):
            np.testing.assert_array_equal(line.get_xydata(), expected_line.get_xydata())


def test_plot_wasserstein_does_not_change_its_inputs():
    seed, non_seed = np.arange(10.0)[::-1].copy(), np.arange(12.0)[::-1] + 0.5
    plot_wasserstein(seed, non_seed, n_interp=50)
    assert seed.tolist() == list(range(9, -1, -1))
    assert non_seed.tolist() == [11.5 - k for k in range(12)]


TEN = [float(k) for k in range(10)]
BAD_ARGUMENTS = [
    pytest.param({"seed": TEN[:9]}, "seed has 9 values; expected at least 10", id="nine-values"),
    pytest.param(
        {"non_seed": [1.0] * 10},
        "non_seed has 10 values that are all 1.0; expected at least 2 distinct values",
        id="one-distinct-value",
    ),
    pytest.param(
        {"seed": [0.0, 1.0, math.nan, *TEN]}, "seed[2] is NaN; expected a number", id="nan"
    ),
    pytest.param({"n_interp": 1}, "n_interp 1 is less than 2; expected 2 or more", id="one-point"),
    pytest.param(
        {"seed": [-1e308, *TEN], "non_seed": [1e308, *TEN]},
        "the values of seed and non_seed go from -1e+308 to 1e+308; "
        "expected values from -9e306 to 9e306 with a range of at most 9e306",
        id="range",
    ),
    pytest.param(
        {"size": (0, 5)},
        "size (0, 5) is not valid; expected (width, height) in inches, two finite numbers above 0",
        id="zero-size",
    ),
]


@pytest.mark.parametrize(("kwargs", "message"), BAD_ARGUMENTS)
def test_plot_wasserstein_rejects_bad_input(kwargs, message):
    arguments = {"seed": TEN, "non_seed": TEN} | kwargs
    with pytest.raises(InputError, match=exactly(message)):
        plot_wasserstein(**arguments)


def test_plot_wasserstein_rejects_a_range_wider_than_the_bound():
    seed, non_seed = np.linspace(-4.6e306, 0, 10), np.linspace(0, 4.6e306, 12)
    message = (
        "the values of seed and non_seed go from -4.6e+306 to 4.6e+306; "
        "expected values from -9e306 to 9e306 with a range of at most 9e306"
    )
    with pytest.raises(InputError, match=exactly(message)):
        plot_wasserstein(seed, non_seed, n_interp=100)


def test_plot_wasserstein_rejects_values_outside_the_bound():
    seed, non_seed = np.linspace(1e308, 1.0001e308, 10), np.linspace(1e308, 1.0001e308, 12)
    message = (
        "the values of seed and non_seed go from 1e+308 to 1.0001e+308; "
        "expected values from -9e306 to 9e306 with a range of at most 9e306"
    )
    with pytest.raises(InputError, match=exactly(message)):
        plot_wasserstein(seed, non_seed, n_interp=100)


@pytest.mark.parametrize("smallest", [1e-308, 1e-310])
def test_plot_wasserstein_rejects_a_bandwidth_below_the_bound(smallest):
    seed = [0.0] * 7 + [smallest, 1.0, 1.0]
    message = (
        r"the bandwidth of the seed sample is \S+; expected at least 4\.4326920044603633e-308, "
    )
    with pytest.raises(InputError, match=message):
        plot_wasserstein(seed, np.linspace(0, 1, 12), n_interp=100)


def test_plot_wasserstein_draws_a_bandwidth_at_the_bound():
    seed = [0.0] * 7 + [2e-307, 1.0, 1.0]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        figure = plot_wasserstein(seed, np.linspace(0, 1, 12), n_interp=100)
        figure.savefig(io.BytesIO(), format="png")


@pytest.mark.parametrize("size", [(15.0, 5.0), (6.0, 2.0)], ids=["default-size", "small-size"])
def test_plot_wasserstein_draws_a_range_at_the_bound(size):
    seed, non_seed = np.linspace(-4.5e306, 0, 10), np.linspace(0, 4.5e306, 12)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        figure = plot_wasserstein(seed, non_seed, n_interp=100, size=size)
        figure.savefig(io.BytesIO(), format="png")


# The ECDF of a sample with a denormal spread divides by a denormal width, and numpy warns about the
# overflow. `wasserstein_areas` has the same limit, and the error of the bandwidth follows it.
@pytest.mark.filterwarnings("ignore:overflow encountered in divide:RuntimeWarning")
@pytest.mark.parametrize("name", ["seed", "non_seed"])
def test_plot_wasserstein_rejects_a_sample_whose_bandwidth_is_zero(name):
    message = (
        f"the bandwidth of the {name} sample is 0.0; expected a finite number above 0, "
        "because the values are too close together"
    )
    arguments = {"seed": TEN, "non_seed": TEN} | {name: np.arange(10) * 5e-324}
    with pytest.raises(InputError, match=exactly(message)):
        plot_wasserstein(**arguments)


def test_plot_wasserstein_runs_on_large_samples_in_bounded_memory():
    generator = np.random.default_rng(1)
    non_seed, seed = generator.normal(size=20000), generator.normal(size=2000) - 0.1
    tracemalloc.start()
    try:
        start = time.perf_counter()
        figure = plot_wasserstein(seed, non_seed, n_interp=10000)
        elapsed = time.perf_counter() - start
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    assert elapsed < 10
    assert peak < 16 * 2**20
    assert figure.axes[0].get_lines()[0].get_xydata().shape == (512, 2)
