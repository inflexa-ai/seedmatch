import io
import math
import re
import subprocess
import sys
import warnings
from functools import partial

import matplotlib as mpl
import numpy as np
import pandas as pd
import pytest
from helpers import exactly, sirna_de, sirna_seqs, wasserstein_samples
from matplotlib.colors import to_hex

from seedmatch import benchmark_predictions, ecdf_test, kmer_counts, kmer_gene_counts
from seedmatch.errors import InputError
from seedmatch.plot import (
    PALETTE,
    plot_alignment,
    plot_benchmark,
    plot_ecdf,
    plot_kmer_histogram,
    plot_ranges,
    plot_wasserstein,
)
from seedmatch.plot._style import (
    TEXT_COLOR,
    apply_layout,
    as_color,
    as_limits,
    as_size,
    check_text_pixels,
    expanded,
    linewidth_points,
    new_figure,
    pasted,
    style_gray_axes,
    style_theme_axes,
)

PLOT_SETUP = """
import io
import sys

import matplotlib
import numpy as np
import pandas as pd

import seedmatch as sm
from seedmatch.plot import (
    plot_alignment,
    plot_benchmark,
    plot_ecdf,
    plot_kmer_histogram,
    plot_ranges,
    plot_wasserstein,
)

rng = np.random.default_rng(0)
de = pd.DataFrame(
    {"gene_id": [f"G{i}" for i in range(40)], "log2FoldChange": rng.normal(size=40)}
)
target, background = list(de["gene_id"][:10]), list(de["gene_id"][10:])
benchmark_de = de.assign(padj=[0.01] * 10 + [1.0] * 30, mer7m8=np.arange(40) % 4)
benchmark_de.loc[:9, "log2FoldChange"] = -1.0
result = sm.ecdf_test(de, target=target, background=background)
PLOT_CALLS = {
    "plot_ecdf": lambda: plot_ecdf(de, {"Background": background, "target": target}, result),
    "plot_wasserstein": lambda: plot_wasserstein(
        de["log2FoldChange"][:10], de["log2FoldChange"][10:], n_interp=100
    ),
    "plot_kmer_histogram": lambda: plot_kmer_histogram(
        pd.Series(rng.integers(0, 30, size=64), index=[f"K{i}" for i in range(64)]), "K0"
    ),
    "plot_alignment": lambda: plot_alignment("3'-UUU-5'\\n    |  \\n5'-...AAA...-3'"),
    "plot_ranges": lambda: plot_ranges(pd.DataFrame({"start": [1, 5, 20], "end": [3, 9, 20]})),
    "plot_benchmark": lambda: plot_benchmark(sm.benchmark_predictions(benchmark_de)),
}
"""

# A fresh process gives the snapshot of the unchanged rc parameters, whatever the order of the
# tests. The forms matplotlib.get_backend(), dict(rcParams), and rcParams == other resolve the
# backend, and that loads pyplot. Thus the snapshot uses dict.items and leaves out the backend key.
PYPLOT_SCRIPT = (
    PLOT_SETUP
    + """
for call in PLOT_CALLS.values():
    call().savefig(io.BytesIO(), format="png")
print("matplotlib.pyplot" in sys.modules)
"""
)

RC_PARAMS_SCRIPT = (
    PLOT_SETUP
    + """
def snapshot():
    return {k: v for k, v in dict.items(matplotlib.rcParams) if k != "backend"}

changed = {}
for name, call in PLOT_CALLS.items():
    before = snapshot()
    call().savefig(io.BytesIO(), format="png")
    after = snapshot()
    changed[name] = sorted(key for key in before if before[key] != after[key])
print(changed)
"""
)


def run_script(code: str) -> subprocess.CompletedProcess[str]:
    # The arguments are the interpreter path and a constant script, so no input reaches the command.
    return subprocess.run(  # noqa: S603
        [sys.executable, "-c", code], check=False, capture_output=True, text=True
    )


def test_plot_package_names_the_extra_without_matplotlib():
    code = (
        "import sys\n"
        "sys.modules['matplotlib'] = None\n"
        "try:\n"
        "    import seedmatch.plot\n"
        "except ImportError as error:\n"
        "    print(type(error.__cause__).__name__)\n"
        "    raise\n"
    )
    done = run_script(code)
    assert done.returncode != 0
    assert done.stdout.strip() == "ModuleNotFoundError"
    assert "ImportError" in done.stderr
    assert "install the extra 'plot'" in done.stderr


def test_plot_functions_load_no_pyplot():
    done = run_script(PYPLOT_SCRIPT)
    assert done.returncode == 0, done.stderr
    assert done.stdout.strip() == "False"


def test_plot_functions_keep_the_rc_params():
    done = run_script(RC_PARAMS_SCRIPT)
    assert done.returncode == 0, done.stderr
    assert done.stdout.strip() == (
        "{'plot_ecdf': [], 'plot_wasserstein': [], 'plot_kmer_histogram': [], "
        "'plot_alignment': [], 'plot_ranges': [], 'plot_benchmark': []}"
    )


def test_the_palette_is_the_palette_of_r():
    assert PALETTE == (
        "black",
        "#3498db",
        "#2ecc71",
        "#f1c40f",
        "#e74c3c",
        "#9b59b6",
        "#1abc9c",
        "#f39c12",
        "#d35400",
    )


@pytest.mark.parametrize(
    ("value", "text"),
    [
        (0.51, "0.51"),
        (748.0, "748"),
        (1e-05, "1e-05"),
        (1000.01, "1000.01"),
        (-0.0, "0"),
        (0.1 + 0.2, "0.3"),
        (1 / 3, "0.333333333333333"),
        (1234567.12345678, "1234567.12345678"),
        (math.nan, "NaN"),
        (math.inf, "Inf"),
        (-math.inf, "-Inf"),
    ],
)
def test_pasted_gives_the_text_of_r(value, text):
    assert pasted(value) == text


@pytest.mark.parametrize(("linewidth", "points"), [(1.5, 3.2009), (1.25, 2.6674), (0.5, 1.0670)])
def test_linewidth_points_converts_the_ggplot2_width(linewidth, points):
    assert linewidth_points(linewidth) == pytest.approx(points, abs=1e-4)


def test_expanded_adds_five_percent_of_the_range_on_each_side():
    assert expanded(-1.0, 1.0) == pytest.approx((-1.1, 1.1), rel=1e-12)
    assert expanded(0.0, 10.0) == pytest.approx((-0.5, 10.5), rel=1e-12)


def test_as_size_returns_floats():
    assert as_size((5, 6.5)) == (5.0, 6.5)
    assert all(isinstance(side, float) for side in as_size([5, 6]))


def test_as_size_accepts_numpy_numbers():
    size = as_size((np.int64(5), np.float32(6)))
    assert size == (5.0, 6.0)
    assert all(type(side) is float for side in size)


@pytest.mark.parametrize(
    "size",
    [
        (0, 5),
        (5, -1),
        (5, math.inf),
        (True, 5),
        (5,),
        (1, 2, 3),
        (10**400, 5),
        "5x5",
        (50, 5),
        (5, 50.0),
        (10**6, 5),
    ],
)
def test_as_size_rejects_a_bad_size(size):
    message = (
        f"size {size!r} is not valid; expected (width, height) in inches, "
        "two finite numbers above 0 and below 50"
    )
    with pytest.raises(InputError, match=exactly(message)):
        as_size(size)


def test_as_size_accepts_a_side_just_below_50():
    assert as_size((49.99, 0.01)) == (49.99, 0.01)


def test_as_limits_returns_floats():
    assert as_limits("x_limits", (-1, 2)) == (-1.0, 2.0)


def test_as_limits_accepts_numpy_numbers():
    limits = as_limits("x_limits", (np.int64(-1), np.float32(2.5)))
    assert limits == (-1.0, 2.5)
    assert all(type(limit) is float for limit in limits)


@pytest.mark.parametrize(
    "limits",
    [
        (1.0, -1.0),
        (1, 1),
        (0.0, math.inf),
        (math.nan, 1.0),
        (False, 1),
        (0,),
        (1, 2, 3),
        (-(10**400), 1),
        (-1e308, 1e308),
        (1e308, 1.0001e308),
        (8.9e306, 9.1e306),
        "01",
    ],
)
def test_as_limits_rejects_a_bad_pair(limits):
    message = (
        f"x_limits {limits!r} is not valid; expected (low, high), "
        "two finite numbers from -9e306 to 9e306 with low < high and high - low at most 9e306"
    )
    with pytest.raises(InputError, match=exactly(message)):
        as_limits("x_limits", limits)


def test_as_color_returns_a_color_text():
    assert as_color("palette[0]", "#3498db") == "#3498db"
    assert as_color("palette[0]", "black") == "black"


def test_as_color_rejects_a_text_that_is_not_a_color():
    message = (
        "palette[1] 'blu' is not a matplotlib color; expected a color name or a hex code, "
        "for example '#3498db'"
    )
    with pytest.raises(InputError, match=exactly(message)):
        as_color("palette[1]", "blu")


def test_as_color_rejects_a_value_that_is_not_text():
    with pytest.raises(InputError, match=exactly("palette[0] 5 has type int; expected a str")):
        as_color("palette[0]", 5)


def test_new_figure_has_the_size_and_the_constrained_layout():
    figure = new_figure((4.0, 3.0))
    assert tuple(figure.get_size_inches()) == (4.0, 3.0)
    assert type(figure.get_layout_engine()).__name__ == "ConstrainedLayoutEngine"
    assert to_hex(figure.get_facecolor()) == "#ffffff"


@pytest.mark.parametrize("size", [(1e-320, 5.0), (5.0, 1e-320), (0.005, 5.0)])
def test_new_figure_rejects_a_side_below_one_pixel_before_an_artist_exists(size):
    # On Linux, np.dot warns for the transforms of the first artist of such a figure.
    with pytest.raises(InputError, match=exactly(too_small(size))):
        new_figure(size)


def test_style_theme_axes_applies_the_theme_of_r():
    axes = new_figure((5.0, 5.0)).add_subplot()
    style_theme_axes(axes)
    assert to_hex(axes.get_facecolor()) == "#ffffff"
    assert {to_hex(spine.get_edgecolor()) for spine in axes.spines.values()} == {"#2c3e50"}
    for axis in (axes.xaxis, axes.yaxis):
        assert to_hex(axis.label.get_color()) == TEXT_COLOR
        assert axis.label.get_fontsize() == 12
        assert {tick.label1.get_fontsize() for tick in axis.get_major_ticks()} == {9.6}
        grid = axis.get_gridlines()[0]
        assert to_hex(grid.get_color()) == "#6b8cae"
        assert grid.get_linestyle() == "--"
        assert grid.get_linewidth() == pytest.approx(linewidth_points(0.2))


def test_style_gray_axes_applies_the_panel_of_theme_grey():
    axes = new_figure((5.0, 5.0)).add_subplot()
    style_gray_axes(axes)
    assert to_hex(axes.get_facecolor()).lower() == "#ebebeb"
    assert not any(spine.get_visible() for spine in axes.spines.values())
    for axis in (axes.xaxis, axes.yaxis):
        grid = axis.get_gridlines()[0]
        assert to_hex(grid.get_color()) == "#ffffff"
        assert grid.get_linewidth() == pytest.approx(linewidth_points(0.5))


def sized_calls():
    """Return a call of each plot function that takes `size`, on the test data."""
    de = sirna_de()
    groups = {
        "Background": de.loc[de["mer7m8"] == 0, "gene_id"],
        "mer7m8": de.loc[de["mer7m8"] >= 1, "gene_id"],
    }
    result = ecdf_test(de, target=groups["mer7m8"], background=groups["Background"])
    seed, non_seed = wasserstein_samples("test_data")
    counts = kmer_gene_counts(kmer_counts(sirna_seqs()))
    ranges = pd.DataFrame({"start": [1, 5, 9], "end": [3, 9, 20]})
    return {
        "plot_ranges": lambda **options: plot_ranges(ranges, **options),
        "plot_ecdf": lambda **options: plot_ecdf(de, groups, result, **options),
        "plot_kmer_histogram": lambda **options: plot_kmer_histogram(counts, "CTCTATA", **options),
        "plot_wasserstein": lambda **options: plot_wasserstein(
            seed, non_seed, n_interp=100, **options
        ),
        "plot_benchmark": lambda **options: plot_benchmark(benchmark_predictions(de), **options),
    }


SIZED = list(sized_calls())


def each_call(size=None):
    """Return a call of each plot function on the test data, with `size` if it takes a size."""
    options = {} if size is None else {"size": size}
    calls = {name: partial(call, **options) for name, call in sized_calls().items()}
    calls["plot_alignment"] = partial(plot_alignment, "3'-UUU-5'\n    |  \n5'-...AAA...-3'")
    return calls


TOO_LOW = re.compile(
    r"figure dpi 3 is too low for the font size ([\d.]+) of a text of the figure; "
    r"expected a figure dpi of at least (\d+), from matplotlib\.rcParams\['figure\.dpi'\]"
)


@pytest.mark.parametrize("name", list(each_call()))
def test_each_plot_function_rejects_a_dpi_too_low_for_its_texts_and_accepts_the_bound(name):
    with mpl.rc_context({"figure.dpi": 3}), pytest.raises(InputError) as error:
        each_call()[name]()
    found = TOO_LOW.fullmatch(str(error.value))
    assert found, str(error.value)
    font_size, bound = float(found[1]), int(found[2])
    assert bound == math.ceil(72 / font_size)
    # A large figure keeps the layout away from the collapse at a low dpi.
    with mpl.rc_context({"figure.dpi": bound}):
        figure = each_call(size=(40.0, 40.0))[name]()
        figure.savefig(io.BytesIO(), format="png")


def test_a_plot_function_rejects_a_font_size_of_the_rc_params_too_small_for_the_dpi():
    message = (
        "figure dpi 30 is too low for the font size 1 of a text of the figure; "
        "expected a figure dpi of at least 72, from matplotlib.rcParams['figure.dpi']"
    )
    with (
        mpl.rc_context({"figure.dpi": 30, "font.size": 1}),
        pytest.raises(InputError, match=exactly(message)),
    ):
        sized_calls()["plot_ranges"]()


def test_apply_layout_ignores_a_hidden_text():
    with mpl.rc_context({"figure.dpi": 30}):
        figure = new_figure((5.0, 5.0))
        figure.add_subplot()
        hidden = figure.text(0.5, 0.5, "x", fontsize=1, visible=False)
        assert apply_layout(figure) is figure
        hidden.set_visible(True)
        with pytest.raises(InputError, match=r"font size 1 of a text"):
            check_text_pixels(figure)


def test_a_plot_function_ignores_an_empty_title_with_a_font_size_too_small_for_the_dpi():
    ranges = pd.DataFrame({"start": [1, 5], "end": [3, 9]})
    with mpl.rc_context({"figure.dpi": 30, "axes.titlesize": 1}):
        plot_ranges(ranges, title="").savefig(io.BytesIO(), format="png")


@pytest.mark.parametrize(
    ("x", "labelsize"),
    [((0.0, 1.0), 1.0), ((1e6, 1e6 + 1), 20.0)],
    ids=["tick labels", "offset text"],
)
def test_check_text_pixels_counts_an_axis_text_that_a_draw_fills(x, labelsize):
    with mpl.rc_context({"figure.dpi": 30, "xtick.labelsize": 1}):
        figure = new_figure((5.0, 5.0))
        axes = figure.add_subplot()
        axes.plot(x, (0.0, 1.0))
        axes.xaxis.set_tick_params(labelsize=labelsize)
        with pytest.raises(InputError, match=r"font size 1 of a text"):
            check_text_pixels(figure)
        with pytest.raises(RuntimeError, match="invalid ppem value"):
            figure.savefig(io.BytesIO(), format="png")


def test_check_text_pixels_ignores_an_empty_offset_text():
    with mpl.rc_context({"figure.dpi": 30, "xtick.labelsize": 1}):
        figure = new_figure((5.0, 5.0))
        axes = figure.add_subplot()
        axes.plot((0.0, 1.0), (0.0, 1.0))
        axes.xaxis.set_tick_params(labelsize=20)
        check_text_pixels(figure)
        figure.savefig(io.BytesIO(), format="png")


def test_check_text_pixels_rejects_a_dpi_with_a_fraction_below_the_bound_of_its_message():
    message = (
        "figure dpi 3.9 is too low for the font size 20 of a text of the figure; "
        "expected a figure dpi of at least 4, from matplotlib.rcParams['figure.dpi']"
    )
    with mpl.rc_context({"figure.dpi": 3.9}):
        figure = new_figure((5.0, 5.0))
        figure.text(0.5, 0.5, "x", fontsize=20)
        with pytest.raises(InputError, match=exactly(message)):
            check_text_pixels(figure)
        figure.set_dpi(4)
        check_text_pixels(figure)


@pytest.mark.parametrize("dpi", [math.nan, math.inf, 0.0, -100.0])
def test_new_figure_rejects_a_dpi_that_is_not_finite_and_above_0(dpi):
    message = (
        f"figure dpi {dpi:g} is not valid; "
        "expected a finite figure dpi above 0, from matplotlib.rcParams['figure.dpi']"
    )
    with (
        mpl.rc_context({"figure.dpi": dpi}),
        pytest.raises(InputError, match=exactly(message)),
    ):
        new_figure((5.0, 5.0))


def too_small(size):
    return (
        f"size {size!r} is too small for the texts of the figure; expected a larger width or height"
    )


@pytest.mark.parametrize("name", SIZED)
def test_each_plot_function_rejects_a_size_too_small_for_its_texts(name):
    with pytest.raises(InputError, match=exactly(too_small((0.3, 0.3)))):
        sized_calls()[name](size=(0.3, 0.3))


def test_plot_ecdf_rejects_a_legend_too_wide_for_the_default_size():
    de = sirna_de()
    groups = {"g" * 60: de.loc[de["mer7m8"] >= 1, "gene_id"]}
    groups["Background"] = de.loc[de["mer7m8"] == 0, "gene_id"]
    result = ecdf_test(de, target=groups["g" * 60], background=groups["Background"])
    with pytest.raises(InputError, match=exactly(too_small((5.0, 5.0)))):
        plot_ecdf(de, groups, result)


def test_the_too_small_message_shows_the_size_as_floats():
    with pytest.raises(InputError, match=exactly(too_small((1.0, 1.0)))):
        sized_calls()["plot_ecdf"](size=(1, 1))


@pytest.mark.parametrize("name", SIZED)
@pytest.mark.parametrize("size", [(1e-320, 5.0), (5.0, 1e-320)])
def test_each_plot_function_rejects_a_side_below_one_pixel(name, size):
    call = sized_calls()[name]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(InputError, match=exactly(too_small(size))):
            call(size=size)


def test_a_layout_warning_other_than_the_collapse_stays_a_warning():
    ranges = pd.DataFrame({"start": [1, 5], "end": [3, 9]})
    with pytest.warns(UserWarning, match=r"Glyph 12354 .* missing from font"):
        plot_ranges(ranges, title="\u3042")


def test_a_layout_warning_that_an_outer_filter_raises_is_not_a_size_error():
    ranges = pd.DataFrame({"start": [1, 5], "end": [3, 9]})
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(UserWarning, match=r"Glyph 12354 .* missing from font"):
            plot_ranges(ranges, title="\u3042")


def test_a_plot_function_restores_the_warning_filters_after_a_collapse():
    call = sized_calls()["plot_ranges"]
    before = list(warnings.filters)
    with pytest.raises(InputError, match=exactly(too_small((0.3, 0.3)))):
        call(size=(0.3, 0.3))
    assert warnings.filters == before


def test_a_plot_function_restores_the_warning_filters_after_a_layout():
    call = sized_calls()["plot_ranges"]
    before = list(warnings.filters)
    call()
    assert warnings.filters == before


# The color bars of plot_benchmark take their width at the first draw, thus the test draws one
# time before it reads the positions.
@pytest.mark.parametrize("name", SIZED)
def test_each_plot_function_keeps_its_layout_for_each_save(name):
    figure = sized_calls()[name]()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        figure.savefig(io.BytesIO(), format="png", dpi=100)
        positions = [axes.get_position().bounds for axes in figure.axes]
        for file_format, dpi in (("png", 50), ("png", 300), ("pdf", 72), ("svg", 72)):
            figure.savefig(io.BytesIO(), format=file_format, dpi=dpi)
            assert [axes.get_position().bounds for axes in figure.axes] == positions


# A color bar keeps its aspect ratio, thus its drawn width follows the size of the figure. The
# original position is the position that the layout gave.
@pytest.mark.parametrize("name", SIZED)
def test_each_plot_function_keeps_its_axes_when_the_size_changes(name):
    figure = sized_calls()[name]()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        figure.savefig(io.BytesIO(), format="png", dpi=100)
        positions = [axes.get_position(original=True).bounds for axes in figure.axes]
        figure.set_size_inches(10, 8)
        figure.savefig(io.BytesIO(), format="png", dpi=100)
        assert [axes.get_position(original=True).bounds for axes in figure.axes] == positions


@pytest.mark.parametrize("name", SIZED)
def test_each_plot_function_lays_out_again_after_set_layout_engine(name):
    figure = sized_calls()[name]()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        figure.savefig(io.BytesIO(), format="png", dpi=100)
        positions = [axes.get_position().bounds for axes in figure.axes]
        figure.set_size_inches(10, 8)
        figure.set_layout_engine("constrained")
        figure.savefig(io.BytesIO(), format="png", dpi=100)
        assert [axes.get_position().bounds for axes in figure.axes] != positions
