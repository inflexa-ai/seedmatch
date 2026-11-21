"""The style that each figure of `seedmatch.plot` shares: the palette, the theme, and the checks."""

import math
from typing import Final

import numpy as np
from matplotlib.axes import Axes
from matplotlib.colors import is_color_like
from matplotlib.figure import Figure
from matplotlib.text import Text

from seedmatch._checks import as_text, shown
from seedmatch.errors import InputError

PALETTE: Final = (
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
TEXT_COLOR: Final = "#34495e"

PANEL_EDGE: Final = "#2c3e50"
_GRID_COLOR: Final = "#6b8cae"
_AXIS_LABEL_SIZE: Final = 12
# The size of `axis.text` in `theme_grey` is rel(0.8), and `SeedMatchR.theme` sets the base size 12.
_TICK_LABEL_SIZE: Final = 9.6
_GRAY_PANEL: Final = "#EBEBEB"
_PANEL_EXPANSION: Final = 0.05
_PAIR_LENGTH: Final = 2


def linewidth_points(linewidth: float) -> float:
    """Return the width in points of a `linewidth` of ggplot2, which counts in units of 0.75 mm."""
    return linewidth * 72.27 / 25.4 * 72 / 96


def pasted(value: float) -> str:
    """Return the text of `paste(value)` of R: 15 significant digits, and `0` for `-0`.

    NaN gives `NaN`, infinity gives `Inf`, and minus infinity gives `-Inf`. R gives `1e+05` for
    100000, but this function gives `100000`.
    """
    if math.isnan(value):
        return "NaN"
    if math.isinf(value):
        return "Inf" if value > 0 else "-Inf"
    text = format(value, ".15g")
    return "0" if text == "-0" else text


def expanded(low: float, high: float) -> tuple[float, float]:
    """Return the limits of an axis with the default expansion of ggplot2, 5% on each side."""
    margin = _PANEL_EXPANSION * (high - low)
    return low - margin, high + margin


def _finite_pair(value: object) -> tuple[float, float] | None:
    """Return the two finite numbers of a `tuple` or a `list`, or `None` for another value."""
    if not isinstance(value, tuple | list) or len(value) != _PAIR_LENGTH:
        return None
    pair: list[float] = []
    for item in value:
        if isinstance(item, bool | np.bool_) or not isinstance(
            item, int | float | np.integer | np.floating
        ):
            return None
        try:
            number = float(item)
        except OverflowError:
            return None
        if not math.isfinite(number):
            return None
        pair.append(number)
    return pair[0], pair[1]


def as_size(size: object) -> tuple[float, float]:
    """Return `size` as two floats, the width and the height in inches. A bool is not a number."""
    pair = _finite_pair(size)
    if pair is None or min(pair) <= 0:
        message = (
            f"size {shown(size)} is not valid; expected (width, height) in inches, "
            "two finite numbers above 0"
        )
        raise InputError(message)
    return pair


def as_limits(name: str, value: object) -> tuple[float, float]:
    """Return `value` as two floats, low and high, with low < high. A bool is not a number."""
    pair = _finite_pair(value)
    if pair is None or pair[0] >= pair[1]:
        message = (
            f"{name} {shown(value)} is not valid; expected (low, high), "
            "two finite numbers with low < high"
        )
        raise InputError(message)
    return pair


def as_color(name: str, value: object) -> str:
    """Return `value` if it is a `str` that matplotlib reads as a color."""
    color = as_text(name, value)
    if not is_color_like(color):
        message = (
            f"{name} {shown(color)} is not a matplotlib color; "
            "expected a color name or a hex code, for example '#3498db'"
        )
        raise InputError(message)
    return color


def as_written[T: Text](text: T) -> T:
    """Make the text draw as written, and return it. A `$` is not the start of a formula.

    Use it for a text from an argument. matplotlib reads each pair of `$` as math, and a bad
    formula raises when the figure draws.
    """
    text.set_parse_math(False)
    return text


def new_figure(size: tuple[float, float]) -> Figure:
    """Return an empty `Figure` of the size in inches. It does not use `matplotlib.pyplot`."""
    return Figure(figsize=size, layout="constrained", facecolor="white")


def style_theme_axes(axes: Axes) -> None:
    """Apply the values of `SeedMatchR.theme` to the Axes: its colors, its grid, and its sizes.

    The sizes are near the sizes of R, not equal. The function changes only the Axes.
    """
    for axis in (axes.xaxis, axes.yaxis):
        axis.label.set_color(TEXT_COLOR)
        axis.label.set_fontsize(_AXIS_LABEL_SIZE)
    axes.tick_params(labelsize=_TICK_LABEL_SIZE)
    axes.set_facecolor("white")
    for spine in axes.spines.values():
        spine.set_edgecolor(PANEL_EDGE)
    axes.grid(visible=True, color=_GRID_COLOR, linestyle="--", linewidth=linewidth_points(0.2))
    axes.set_axisbelow(True)


def style_gray_axes(axes: Axes) -> None:
    """Apply the panel of `theme_grey` of ggplot2 to the Axes: gray, white grid, and no edge.

    The function changes only the Axes.
    """
    axes.set_facecolor(_GRAY_PANEL)
    for spine in axes.spines.values():
        spine.set_visible(False)
    axes.grid(visible=True, color="white", linestyle="-", linewidth=linewidth_points(0.5))
    axes.set_axisbelow(True)
