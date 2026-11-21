import io
import math
import time
import tracemalloc
import warnings

import numpy as np
import pytest
from helpers import exactly, golden, golden_note, sirna_de
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import LineCollection, QuadMesh
from matplotlib.colors import hsv_to_rgb, to_hex

from seedmatch import benchmark_predictions
from seedmatch.errors import InputError
from seedmatch.plot import plot_benchmark
from seedmatch.stats import Benchmark, pr_curve, roc_curve

BENCHMARK = "plots/benchmark_plot.tsv"
SIZE_MESSAGE = (
    "size (0, 5) is not valid; expected (width, height) in inches, two finite numbers above 0"
)
LABELS = [1, 0] * 10
THRESHOLD_MESSAGE = (
    "the {} curve has the finite threshold {}; "
    "expected a finite threshold from -9e305 to 9e305 for the colors"
)
R_ERROR = "PR error increasing 'x' and 'y' values expected"
# Points of a line in base graphics: `lwd = 3` is 3 * 0.75 points.
LWD_POINTS = 2.25


def benchmark():
    return benchmark_predictions(sirna_de())


def of_scores(scores):
    """Return the benchmark of the PR curve and the ROC curve of the scores with 10 positives."""
    return Benchmark(pr=pr_curve(scores, LABELS), roc=roc_curve(scores, LABELS))


def lines_of(axes):
    (collection,) = axes.collections
    assert isinstance(collection, LineCollection)
    return collection


def colors_of(collection):
    return [to_hex(color) for color in collection.get_colors()]


def colors_by_threshold(axes, curve):
    """Map the threshold of each segment, the threshold of its first row, to its set of colors."""
    colors: dict[float, set[str]] = {}
    for threshold, color in zip(curve["threshold"][:-1], colors_of(lines_of(axes)), strict=True):
        colors.setdefault(threshold, set()).add(color)
    return colors


def r_segments(case, curve):
    """Return the segments of the case and the curve in the R golden, and their colors."""
    table = golden(BENCHMARK)
    rows = table[(table["case"] == case) & (table["curve"] == curve)]
    coordinates = rows[["x0", "y0", "x1", "y1"]].astype(float).to_numpy()
    segments = coordinates.reshape(-1, 2, 2)
    return segments, [color.lower() for color in rows["color"]]


def drawn_color_bar(figure, index):
    """Return the face colors and the range of the color bar of the curve, after a draw."""
    FigureCanvasAgg(figure).draw()
    axes = figure.axes[2 + index]
    (mesh,) = (collection for collection in axes.collections if isinstance(collection, QuadMesh))
    return [to_hex(color) for color in mesh.get_facecolor()], axes.get_ylim()


def test_plot_benchmark_draws_the_r_segments():
    result = benchmark()
    figure = plot_benchmark(result)
    pr_axes, roc_axes, *bars = figure.axes
    assert len(bars) == 2
    pr_lines, roc_lines = lines_of(pr_axes), lines_of(roc_axes)
    assert (len(pr_lines.get_segments()), len(roc_lines.get_segments())) == (105, 3)

    for lines, curve, columns, name in (
        (pr_lines, result.pr.curve, ["recall", "precision"], "PR"),
        (roc_lines, result.roc.curve, ["fpr", "tpr"], "ROC"),
    ):
        points = curve[columns].to_numpy()
        port = np.array(lines.get_segments())
        np.testing.assert_array_equal(port[:, 0], points[:-1])
        np.testing.assert_array_equal(port[:, 1], points[1:])
        r, r_colors = r_segments("benchmark_1", name)
        assert port.shape == r.shape
        np.testing.assert_allclose(port, r, rtol=1e-12, atol=0)
        assert colors_of(lines) == r_colors

    pr_segments = np.array(pr_lines.get_segments())
    np.testing.assert_allclose(pr_segments[0], [(1.0, 63 / 179), (1.0, 63 / 179)], rtol=1e-12)
    pr_colors = colors_of(pr_lines)
    assert (pr_colors[0], pr_colors[-1]) == ("#ff0000", "#cc00ff")
    assert set(pr_colors) == {"#ff0000", "#00ff60", "#cc00ff"}
    np.testing.assert_allclose(
        np.array(roc_lines.get_segments()),
        [
            [(0, 0), (2 / 116, 0)],
            [(2 / 116, 0), (4 / 116, 14 / 63)],
            [(4 / 116, 14 / 63), (1, 1)],
        ],
        rtol=1e-12,
        atol=0,
    )
    assert colors_of(roc_lines) == ["#cc00ff", "#00ff60", "#ff0000"]


def test_plot_benchmark_gives_the_r_titles_and_labels():
    figure = plot_benchmark(benchmark())
    pr_axes, roc_axes = figure.axes[:2]
    for axes, name in ((pr_axes, "PR"), (roc_axes, "ROC")):
        note = golden_note(BENCHMARK, rf"^# title {name}: (.*)$")
        assert axes.get_title() == note.replace(" | ", "\n")
        assert f"xlab {axes.get_xlabel()} / ylab {axes.get_ylabel()}" == golden_note(
            BENCHMARK, rf"^# labels {name}: (.*)$"
        )
    assert pr_axes.get_title() == "PR curve\nAUC = 0.4723177"
    assert roc_axes.get_title() == "ROC curve\nAUC = 0.591954"

    r_width = float(golden_note(BENCHMARK, r"^# lwd: (.*)$"))
    limits = [float(value) for value in golden_note(BENCHMARK, r"^# usr: (.*)$").split(",")]
    for axes in (pr_axes, roc_axes):
        assert list(lines_of(axes).get_linewidth()) == [r_width * 0.75] == [LWD_POINTS]
        assert axes.get_xlim() == pytest.approx(limits[:2], rel=1e-12)
        assert axes.get_ylim() == pytest.approx(limits[2:], rel=1e-12)
    assert pr_axes.get_xlim() == pytest.approx((-0.04, 1.04), rel=1e-12)
    assert tuple(figure.get_size_inches()) == (12.0, 5.0)

    # R gives the labels to `scores.class1`, thus its AUC is not the AUC of the port (`D-006`).
    unpatched = golden_note(BENCHMARK, r"^# title ROC unpatched: (.*)$")
    assert unpatched == "ROC curve | AUC = 0.3762679"
    assert roc_axes.get_title().replace("\n", " | ") != unpatched


def test_plot_benchmark_has_a_color_bar_for_each_curve():
    figure = plot_benchmark(benchmark())
    assert len(figure.axes) == 4
    scale = golden_note(BENCHMARK, r"^# scale: (.*)$")
    for index, name in enumerate(["PR", "ROC"]):
        colors, limits = drawn_color_bar(figure, index)
        assert len(colors) == 100
        assert colors[0] == "#ff0000"
        assert colors[-1] == "#cc00ff"
        for entry in scale.split(" / "):
            position, color = entry.split(" ")
            assert colors[int(position) - 1] == color.lower()
        note = golden_note(BENCHMARK, rf"^# color bar {name}: (.*)$")
        assert note == "first y 0 / last y 2 / colors 100"
        assert limits == (0, 2)


@pytest.mark.parametrize(
    ("score", "r_note", "r_text"),
    [
        (0.0, "one score 0", f"{R_ERROR} / segment colors NA"),
        (-3.0, "one score -3", f"{R_ERROR} / segment colors #00FF85"),
        (2.0, "one_score_2 color bar", "first y 1.8 / last y 2.18 / colors 100"),
    ],
)
def test_plot_benchmark_colors_one_distinct_score_D069(score, r_note, r_text):
    assert golden_note(BENCHMARK, rf"^# {r_note}: (.*)$") == r_text
    figure = plot_benchmark(of_scores(np.full(20, score)))
    for index, axes in enumerate(figure.axes[:2]):
        lines = lines_of(axes)
        assert len(colors_of(lines)) == len(lines.get_segments())
        assert set(colors_of(lines)) == {"#00ff60"}
        colors, limits = drawn_color_bar(figure, index)
        assert limits == pytest.approx((score - 0.5, score + 0.5))
        assert len(colors) == 100

    r, r_colors = r_segments("one_score_2", "PR")
    np.testing.assert_allclose(np.array(lines_of(figure.axes[0]).get_segments()), r, rtol=1e-12)
    assert set(r_colors) == {"#00ff85"}


def test_plot_benchmark_colors_an_infinite_score_D069():
    scores = np.array([math.inf] + [1.0] * 9 + [0.0] * 10)
    result = of_scores(scores)
    figure = plot_benchmark(result)
    for index, (axes, curve) in enumerate(
        zip(figure.axes[:2], [result.pr.curve, result.roc.curve], strict=True)
    ):
        assert colors_by_threshold(axes, curve) == {
            0.0: {"#ff0000"},
            1.0: {"#cc00ff"},
            math.inf: {"#cc00ff"},
        }
        assert drawn_color_bar(figure, index)[1] == (0, 1)

    r_note = golden_note(BENCHMARK, r"^# an Inf score: (.*)$")
    assert r_note == "PR error 'to' must be a finite number / segment colors #FF0000,NA"
    assert r_note.split("segment colors ")[1].split(",")[0].lower() == "#ff0000"


@pytest.mark.parametrize("curve", ["PR", "ROC"])
@pytest.mark.parametrize(
    ("low", "high", "first_in_pr", "first_in_roc"),
    [
        (-1.7e308, 1.7e308, "-1.7e+308", "1.7e+308"),
        (9e307, 9.1e307, "9e+307", "9.1e+307"),
        (0.0, 1e306, "1e+306", "1e+306"),
        (-9.1e305, 0.0, "-9.1e+305", "-9.1e+305"),
        (
            1.7976931348623157e308,
            1.7976931348623157e308,
            "1.7976931348623157e+308",
            "1.7976931348623157e+308",
        ),
    ],
)
def test_plot_benchmark_rejects_a_threshold_that_is_too_large_D069(
    curve, low, high, first_in_pr, first_in_roc
):
    wide, narrow = np.array([low, high] * 10), np.array([0.0, 1.0] * 10)
    pr_scores, roc_scores = (wide, narrow) if curve == "PR" else (narrow, wide)
    result = Benchmark(pr=pr_curve(pr_scores, LABELS), roc=roc_curve(roc_scores, LABELS))
    message = THRESHOLD_MESSAGE.format(curve, first_in_pr if curve == "PR" else first_in_roc)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(InputError, match=exactly(message)):
            plot_benchmark(result)


@pytest.mark.parametrize(
    ("low", "high", "expected"),
    [
        (-9e305, 9e305, {"#ff0000", "#cc00ff"}),
        (0.0, 9e305, {"#ff0000", "#cc00ff"}),
        (9e305, 9e305, {"#00ff60"}),
        (-9e305, -9e305, {"#00ff60"}),
    ],
)
def test_plot_benchmark_takes_a_threshold_just_inside_the_bound_D069(low, high, expected):
    result = of_scores(np.array([low, high] * 10))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        figure = plot_benchmark(result)
        out = io.BytesIO()
        figure.savefig(out, format="png")
    assert out.getvalue().startswith(b"\x89PNG")
    for axes in figure.axes[:2]:
        assert set(colors_of(lines_of(axes))) == expected


def test_plot_benchmark_colors_a_negative_infinite_score_D069():
    result = of_scores(np.array([-math.inf] + [1.0] * 9 + [0.0] * 10))
    figure = plot_benchmark(result)
    for axes, curve in zip(figure.axes[:2], [result.pr.curve, result.roc.curve], strict=True):
        assert colors_by_threshold(axes, curve) == {
            -math.inf: {"#ff0000"},
            0.0: {"#ff0000"},
            1.0: {"#cc00ff"},
        }


def test_plot_benchmark_colors_only_infinite_scores():
    figure = plot_benchmark(of_scores(np.full(20, math.inf)))
    for index, axes in enumerate(figure.axes[:2]):
        assert set(colors_of(lines_of(axes))) == {"#00ff60"}
        assert drawn_color_bar(figure, index)[1] == (-0.5, 0.5)


def test_plot_benchmark_rounds_the_color_index_half_to_even():
    # The thresholds 0, 0.5, 1 give the positions 0, 49.5, and 99. R rounds 50.5 to 50, thus the
    # index from 0 is 49, and the color is the 50th color of the scale.
    scores = np.array([0.0, 0.5, 1.0] * 6 + [0.0, 0.5])
    labels = [1, 0] * 10
    result = Benchmark(pr=pr_curve(scores, labels), roc=roc_curve(scores, labels))
    figure = plot_benchmark(result)
    colors = colors_of(lines_of(figure.axes[1]))
    thresholds = result.roc.curve["threshold"][:-1]
    assert dict(zip(thresholds, colors, strict=True))[0.5] == "#00ff60"


def test_plot_benchmark_gives_each_color_bar_the_range_of_its_curve():
    result = Benchmark(
        pr=pr_curve(np.array([0.0, 1.0] * 10), LABELS),
        roc=roc_curve(np.array([0.0, 5.0] * 10), LABELS),
    )
    figure = plot_benchmark(result)
    assert drawn_color_bar(figure, 0)[1] == (0, 1)
    assert drawn_color_bar(figure, 1)[1] == (0, 5)


def test_plot_benchmark_puts_each_color_bar_to_the_right_of_its_axes():
    figure = plot_benchmark(benchmark())
    FigureCanvasAgg(figure).draw()
    pr, roc, pr_bar, roc_bar = (axes.get_position() for axes in figure.axes)
    assert pr.x1 <= pr_bar.x0 <= pr_bar.x1 <= roc.x0
    assert roc.x1 <= roc_bar.x0


def test_plot_benchmark_computes_the_position_in_the_order_of_r():
    # R gives round(1 + 99 * 101 / 198) = round(51.5) = 52, thus the index from 0 is 51. The
    # product after the division gives a position just below 50.5, and then the index 50.
    scores = np.array([0.0, 101.0, 198.0] * 6 + [0.0, 101.0])
    result = of_scores(scores)
    colors = colors_by_threshold(plot_benchmark(result).axes[1], result.roc.curve)
    assert colors[101.0] == {to_hex(hsv_to_rgb((0.8 * 51 / 99, 1.0, 1.0)))}


def test_plot_benchmark_takes_a_size():
    assert tuple(plot_benchmark(benchmark(), size=(8, 3)).get_size_inches()) == (8.0, 3.0)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"benchmark": "x"}, "benchmark has type str; expected a Benchmark"),
        ({"benchmark": benchmark().pr}, "benchmark has type PrCurve; expected a Benchmark"),
        ({"size": (0, 5)}, SIZE_MESSAGE),
    ],
)
def test_plot_benchmark_rejects_bad_input(kwargs, message):
    arguments = {"benchmark": benchmark()} | kwargs
    with pytest.raises(InputError, match=exactly(message)):
        plot_benchmark(**arguments)


def test_plot_benchmark_draws_to_a_file():
    figure = plot_benchmark(benchmark())
    out = io.BytesIO()
    figure.savefig(out, format="png")
    assert out.getvalue().startswith(b"\x89PNG")


def test_plot_benchmark_runs_on_many_rows_in_bounded_time_and_memory():
    rng = np.random.default_rng(0)
    scores = rng.permutation(50_000).astype(float)
    labels = rng.random(50_000) < 0.3
    result = Benchmark(pr=pr_curve(scores, labels), roc=roc_curve(scores, labels))
    rows = len(result.pr.curve) + len(result.roc.curve)
    tracemalloc.start()
    try:
        begin = time.perf_counter()
        figure = plot_benchmark(result)
        elapsed = time.perf_counter() - begin
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    segments = sum(len(lines_of(axes).get_segments()) for axes in figure.axes[:2])
    assert segments == rows - 2
    assert elapsed < 10
    assert peak < 128 * 2**20
