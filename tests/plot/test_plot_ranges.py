import io
import math
import time
import tracemalloc
import warnings

import numpy as np
import pandas as pd
import pytest
from helpers import GUIDE, exactly, golden, golden_note, sirna_seqs
from matplotlib.colors import to_hex

from seedmatch import Guide, Seed, find_matches
from seedmatch.errors import InputError
from seedmatch.plot import plot_ranges

RANGES = "plots/ranges_plot.tsv"
SIZE_MESSAGE = (
    "size (0, 5) is not valid; expected (width, height) in inches, "
    "two finite numbers above 0 and below 50"
)
BAR_HEIGHT_MESSAGE = "bar_height {} is out of range; expected a number above 0 and at most 1"
POSITION_MESSAGE = "ranges table has {} {} in row {}; expected a position from -2^52 to 2^52"
LIMIT = 2**52


def eight():
    """Return the table `eight` of the R golden: the starts 10, 30, ..., 150, and start + 6 ends."""
    starts = np.arange(10, 151, 20)
    return pd.DataFrame({"start": starts, "end": starts + 6})


def boxes(figure):
    """Return the four vertices of each bar of the only Axes, as an array of shape (n, 4, 2)."""
    (axes,) = figure.axes
    (collection,) = axes.collections
    return np.array([path.vertices[:4] for path in collection.get_paths()])


def r_rows(case):
    """Return the layer rows of the case in the R golden, with NaN for NA."""
    table = golden(RANGES)
    rows = table.loc[table["case"] == case, ["xmin", "xmax", "ymin", "ymax"]]
    return rows.replace("NA", "nan").astype(float).reset_index(drop=True)


def test_plot_ranges_draws_each_row_D068():
    figure = plot_ranges(eight())
    (axes,) = figure.axes
    (collection,) = axes.collections
    paths = collection.get_paths()
    assert len(paths) == 8
    np.testing.assert_allclose(
        paths[0].vertices[:4],
        [(9.5, 0.875), (16.5, 0.875), (16.5, 1.125), (9.5, 1.125)],
        rtol=0,
        atol=1e-12,
    )
    assert axes.get_xlim() == pytest.approx((2.15, 163.85), rel=1e-12)
    assert axes.get_ylim() == (0, 9)
    assert to_hex(collection.get_facecolor()[0]) == "#0000ff"
    assert to_hex(collection.get_edgecolor()[0]) == "#0000ff"
    assert list(collection.get_linewidth()) == [0.0]
    assert (axes.get_title(), axes.get_xlabel(), axes.get_ylabel()) == (
        "Ranges",
        "Position",
        "Ranges",
    )
    assert tuple(figure.get_size_inches()) == (7.0, 5.0)

    r = r_rows("eight_rows")
    assert len(r) == 6
    port = boxes(figure)[:6]
    np.testing.assert_allclose(port[:, 0, 0] + 0.5, r["xmin"])
    np.testing.assert_allclose(port[:, 1, 0] - 0.5, r["xmax"])
    np.testing.assert_allclose(port[:, 0, 1], r["ymin"])
    np.testing.assert_allclose(port[:, 2, 1], r["ymax"])
    assert len(boxes(figure)) == 8 > len(r)


def test_plot_ranges_has_one_collection_and_no_other_artist():
    (axes,) = plot_ranges(eight()).axes
    assert len(axes.collections) == 1
    assert (len(axes.patches), len(axes.lines), len(axes.texts)) == (0, 0, 0)


def test_plot_ranges_uses_the_options_D068():
    assert golden_note(RANGES, r"^# options: (.*)$") == (
        "title none / coord limits none / layer rows 6"
    )
    assert golden_note(RANGES, r"^# options call: (.*)$") == (
        'plotRangesDF(eight, xlim = c(0, 1000), main = "My title", col = "red")'
    )
    assert golden_note(RANGES, r"^# bar_height 0.8: ymin ([\d.]+) ") == "0.6"
    assert golden_note(RANGES, r"^# bar_height 0.8: ymin [\d.]+ / ymax ([\d.]+)$") == "1.4"

    figure = plot_ranges(
        eight(), x_limits=(0.0, 1000.0), title="My title", color="red", bar_height=0.8, size=(4, 3)
    )
    (axes,) = figure.axes
    (collection,) = axes.collections
    assert axes.get_xlim() == (0, 1000)
    assert axes.get_title() == "My title"
    assert to_hex(collection.get_facecolor()[0]) == "#ff0000"
    assert to_hex(collection.get_edgecolor()[0]) == "#ff0000"
    assert tuple(figure.get_size_inches()) == (4.0, 3.0)
    first = boxes(figure)[0]
    np.testing.assert_allclose([first[0, 1], first[2, 1]], [0.6, 1.4])
    assert axes.get_ylim() == (0, 9)


def test_plot_ranges_puts_each_row_at_its_position_D068():
    subset = eight().iloc[[2, 6, 7]]
    centres = boxes(plot_ranges(subset))[:, :, 1].mean(axis=1)
    np.testing.assert_allclose(centres, [1, 2, 3])
    named = subset.set_axis(["a", "b", "c"])
    np.testing.assert_allclose(boxes(plot_ranges(named))[:, :, 1].mean(axis=1), [1, 2, 3])

    r = r_rows("subset")
    np.testing.assert_allclose(r["ymin"], [2.875, 6.875, 7.875])
    np.testing.assert_allclose(boxes(plot_ranges(subset))[:, 0, 0] + 0.5, r["xmin"])
    assert r_rows("text_row_names")["ymin"].isna().all()
    assert golden_note(RANGES, r"^# text_row_names warning: (.*)$") == "NAs introduced by coercion"


def test_plot_ranges_shows_a_range_of_one_base_D068():
    table = pd.DataFrame({"start": [5, 20], "end": [10, 20]})
    r = r_rows("one_base")
    assert r["xmax"][1] - r["xmin"][1] == 0
    assert golden_note(RANGES, r"^# plotIRanges one_base warning: (.*)$") == (
        "zero-length arrow is of indeterminate angle and so skipped"
    )
    port = boxes(plot_ranges(table))
    assert (port[1, 0, 0], port[1, 1, 0]) == (19.5, 20.5)
    assert (port[0, 0, 0], port[0, 1, 0]) == (4.5, 10.5)


def test_plot_ranges_draws_a_find_matches_result():
    matches = find_matches(Guide(GUIDE).target(Seed.MER7M8), sirna_seqs())
    assert len(matches) == 20
    first = matches.iloc[0]
    assert (first["gene_id"], first["start"], first["end"]) == ("ENSRNOG00000016275", 73, 79)
    port = boxes(plot_ranges(matches))
    assert port.shape == (20, 4, 2)
    assert (port[0, 0, 0], port[0, 1, 0]) == (72.5, 79.5)
    np.testing.assert_allclose(port[:, 0, 0], matches["start"] - 0.5)
    np.testing.assert_allclose(port[:, 1, 0], matches["end"] + 0.5)


def test_plot_ranges_takes_a_nullable_integer_table():
    table = pd.DataFrame(
        {"start": pd.array([1, 5], dtype="Int64"), "end": pd.array([3, 9], dtype="Int64")}
    )
    port = boxes(plot_ranges(table))
    np.testing.assert_allclose(port[:, 0, 0], [0.5, 4.5])
    np.testing.assert_allclose(port[:, 1, 0], [3.5, 9.5])


def test_plot_ranges_takes_the_limits_from_the_ranges_with_an_expansion():
    table = pd.DataFrame({"start": [-4, 10], "end": [0, 50]})
    (axes,) = plot_ranges(table).axes
    low, high = -4.5, 50.5
    assert axes.get_xlim() == pytest.approx((low - 0.05 * (high - low), high + 0.05 * (high - low)))
    assert axes.get_ylim() == (0, 3)


def test_plot_ranges_has_the_gray_panel_of_the_default_theme():
    (axes,) = plot_ranges(eight()).axes
    assert to_hex(axes.get_facecolor()).lower() == "#ebebeb"
    assert not any(spine.get_visible() for spine in axes.spines.values())


@pytest.mark.parametrize("title", ["$x$", r"cost $5 and $6", r"$\badcommand$"])
def test_plot_ranges_draws_the_title_as_written(title):
    figure = plot_ranges(eight(), title=title)
    assert figure.axes[0].title.get_parse_math() is False
    figure.savefig(io.BytesIO(), format="png")


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"ranges": [1, 2]}, "ranges has type list; expected a pandas.DataFrame"),
        (
            {"ranges": eight().drop(columns="end")},
            "ranges table has missing columns: ['end']; expected ['start', 'end']",
        ),
        (
            {"ranges": eight().astype({"start": "float64"})},
            "ranges table column 'start' has dtype float64; expected an integer dtype",
        ),
        (
            {"ranges": eight().iloc[:0]},
            "ranges table has no row; expected at least one range",
        ),
        (
            {"ranges": pd.DataFrame({"start": [1, 3, 9, 20], "end": [2, 4, 7, 10]})},
            "ranges table has start 9 and end 7 in row 2; expected end >= start",
        ),
        (
            {"x_limits": (1.0, -1.0)},
            (
                "x_limits (1.0, -1.0) is not valid; expected (low, high), "
                "two finite numbers from -9e306 to 9e306 with low < high "
                "and high - low at most 9e306"
            ),
        ),
        (
            {"ranges": pd.DataFrame({"start": [LIMIT + 1], "end": [LIMIT + 1]})},
            POSITION_MESSAGE.format("start", LIMIT + 1, 0),
        ),
        (
            {"ranges": pd.DataFrame({"start": [5, -LIMIT - 1], "end": [6, 7]})},
            POSITION_MESSAGE.format("start", -LIMIT - 1, 1),
        ),
        (
            {"ranges": pd.DataFrame({"start": [1, 2], "end": [3, LIMIT + 1]})},
            POSITION_MESSAGE.format("end", LIMIT + 1, 1),
        ),
        (
            {"ranges": pd.DataFrame({"start": [np.iinfo(np.int64).min], "end": [0]})},
            POSITION_MESSAGE.format("start", np.iinfo(np.int64).min, 0),
        ),
        (
            {"ranges": pd.DataFrame({"start": [1], "end": np.array([2**64 - 1], dtype=np.uint64)})},
            POSITION_MESSAGE.format("end", 2**64 - 1, 0),
        ),
        (
            {"ranges": pd.DataFrame({"start": [0], "end": [-LIMIT - 1]})},
            POSITION_MESSAGE.format("end", -LIMIT - 1, 0),
        ),
        (
            {
                "ranges": pd.DataFrame(
                    {"start": [0, LIMIT + 1, LIMIT + 5], "end": [LIMIT + 1, LIMIT + 1, LIMIT + 5]}
                )
            },
            POSITION_MESSAGE.format("start", LIMIT + 1, 1),
        ),
        (
            {"ranges": pd.DataFrame({"start": [0, 1, 2], "end": [5, LIMIT + 1, LIMIT + 9]})},
            POSITION_MESSAGE.format("end", LIMIT + 1, 1),
        ),
        ({"title": 5}, "title 5 has type int; expected a str"),
        (
            {"color": "blu"},
            (
                "color 'blu' is not a matplotlib color; expected a color name or a hex code, "
                "for example '#3498db'"
            ),
        ),
        ({"color": 5}, "color 5 has type int; expected a str"),
        ({"bar_height": 0.0}, BAR_HEIGHT_MESSAGE.format("0.0")),
        ({"bar_height": 1.5}, BAR_HEIGHT_MESSAGE.format("1.5")),
        ({"bar_height": math.nan}, BAR_HEIGHT_MESSAGE.format("nan")),
        ({"bar_height": True}, "bar_height True has type bool; expected an int or a float"),
        ({"size": (0, 5)}, SIZE_MESSAGE),
    ],
)
def test_plot_ranges_rejects_bad_input(kwargs, message):
    arguments = {"ranges": eight()} | kwargs
    with pytest.raises(InputError, match=exactly(message)):
        plot_ranges(**arguments)


@pytest.mark.parametrize(
    "limits", [(-1e308, 1e308), (-4.6e306, 4.6e306), (1e308, 1.0001e308), (8.9e306, 9.1e306)]
)
def test_plot_ranges_rejects_x_limits_wider_than_the_bound(limits):
    message = (
        f"x_limits {limits!r} is not valid; expected (low, high), "
        "two finite numbers from -9e306 to 9e306 with low < high and high - low at most 9e306"
    )
    with pytest.raises(InputError, match=exactly(message)):
        plot_ranges(eight(), x_limits=limits)


@pytest.mark.parametrize("size", [(7.0, 5.0), (1.0, 1.0)], ids=["default-size", "small-size"])
def test_plot_ranges_draws_x_limits_at_the_bound(size):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        figure = plot_ranges(eight(), x_limits=(-4.5e306, 4.5e306), size=size)
        figure.savefig(io.BytesIO(), format="png")


@pytest.mark.parametrize("position", [LIMIT, -LIMIT])
def test_plot_ranges_takes_a_range_of_one_base_at_the_float_limit(position):
    table = pd.DataFrame({"start": [position], "end": [position]})
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        figure = plot_ranges(table)
        figure.savefig(io.BytesIO(), format="png")
    (axes,) = figure.axes
    low, high = axes.get_xlim()
    assert low < high


def test_plot_ranges_keeps_the_vertices_of_a_large_position_exact():
    table = pd.DataFrame({"start": [2**40, 10**9, LIMIT - 1], "end": [2**40 + 5, 10**9, LIMIT - 1]})
    port = boxes(plot_ranges(table))
    np.testing.assert_array_equal(port[:, 0, 0], [2**40 - 0.5, 10**9 - 0.5, LIMIT - 1.5])
    np.testing.assert_array_equal(port[:, 1, 0], [2**40 + 5.5, 10**9 + 0.5, LIMIT - 0.5])


def test_plot_ranges_takes_the_whole_bar_height_range():
    first_row = boxes(plot_ranges(eight(), bar_height=1))[0]
    np.testing.assert_allclose([first_row[0, 1], first_row[2, 1]], [0.5, 1.5])


def test_plot_ranges_does_not_change_ranges():
    table = eight()
    before = table.copy()
    plot_ranges(table)
    pd.testing.assert_frame_equal(table, before)


def test_plot_ranges_runs_on_many_rows_in_bounded_time_and_memory():
    starts = np.arange(100_000) * 10
    table = pd.DataFrame({"start": starts, "end": starts + 6})
    tracemalloc.start()
    try:
        begin = time.perf_counter()
        figure = plot_ranges(table)
        elapsed = time.perf_counter() - begin
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    (axes,) = figure.axes
    assert len(axes.collections) == 1
    paths = axes.collections[0].get_paths()
    assert len(paths) == 100_000
    # A 3-D array of vertices gives paths that share their codes. A list of arrays gives each path
    # its own codes, and then the call takes more time and memory.
    assert paths[0].codes is paths[1].codes
    assert elapsed < 10
    assert peak < 64 * 2**20
