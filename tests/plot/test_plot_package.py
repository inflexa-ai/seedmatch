import math
import subprocess
import sys

import numpy as np
import pytest
from helpers import exactly
from matplotlib.colors import to_hex

from seedmatch.errors import InputError
from seedmatch.plot import PALETTE
from seedmatch.plot._style import (
    TEXT_COLOR,
    as_color,
    as_limits,
    as_size,
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
from seedmatch.plot import plot_ecdf

rng = np.random.default_rng(0)
de = pd.DataFrame(
    {"gene_id": [f"G{i}" for i in range(40)], "log2FoldChange": rng.normal(size=40)}
)
target, background = list(de["gene_id"][:10]), list(de["gene_id"][10:])
result = sm.ecdf_test(de, target=target, background=background)
PLOT_CALLS = {
    "plot_ecdf": lambda: plot_ecdf(de, {"Background": background, "target": target}, result),
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
    assert done.stdout.strip() == "{'plot_ecdf': []}"


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
    "size", [(0, 5), (5, -1), (5, math.inf), (True, 5), (5,), (1, 2, 3), (10**400, 5), "5x5"]
)
def test_as_size_rejects_a_bad_size(size):
    message = (
        f"size {size!r} is not valid; expected (width, height) in inches, "
        "two finite numbers above 0"
    )
    with pytest.raises(InputError, match=exactly(message)):
        as_size(size)


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
        "01",
    ],
)
def test_as_limits_rejects_a_bad_pair(limits):
    message = (
        f"x_limits {limits!r} is not valid; expected (low, high), "
        "two finite numbers with low < high"
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
