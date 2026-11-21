"""The style that each figure of `seedmatch.plot` shares: the palette, the theme, and the checks."""

import math
import warnings
from typing import TYPE_CHECKING, Final, cast

import matplotlib as mpl
import numpy as np
from matplotlib.axes import Axes
from matplotlib.axis import Axis
from matplotlib.colors import is_color_like
from matplotlib.figure import Figure
from matplotlib.text import Text

from seedmatch._checks import as_text, shown
from seedmatch.errors import InputError

if TYPE_CHECKING:
    from matplotlib.layout_engine import LayoutEngine

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
# The `MaxNLocator` of matplotlib 3.11.1 overflows when the width of a view divided by its bin
# count reaches 1e307, and the bin count can be 1. Near the largest float, the ticks after the end
# of a view also overflow. The expansion of 5% on each side of a range of 9e306 gives 9.9e306.
AXIS_BOUND: Final = 9e306
# ggsave of ggplot2 4.0.2 rejects a side of 50 inches or more. The bound also keeps the canvas of
# apply_layout below about 95 MiB at 100 dpi.
_MAX_SIDE: Final = 50.0
# The start of the UserWarning of matplotlib 3.11.1 for a constrained layout that collapses.
_COLLAPSE: Final = "constrained_layout not applied"
_POINTS_PER_INCH: Final = 72
_DPI_FROM: Final = "from matplotlib.rcParams['figure.dpi']"


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
    if pair is None or min(pair) <= 0 or max(pair) >= _MAX_SIDE:
        message = (
            f"size {shown(size)} is not valid; expected (width, height) in inches, "
            "two finite numbers above 0 and below 50"
        )
        raise InputError(message)
    return pair


def _inside_axis_bound(low: float, high: float) -> bool:
    """Return whether `low` and `high` are inside +-`AXIS_BOUND`, at most `AXIS_BOUND` apart."""
    return low >= -AXIS_BOUND and high <= AXIS_BOUND and high - low <= AXIS_BOUND


def check_axis_range(name: str, low: float, high: float) -> None:
    """Raise `InputError` if the values from `low` to `high` cannot be on one x axis."""
    if not _inside_axis_bound(low, high):
        message = (
            f"the values of {name} go from {low} to {high}; "
            "expected values from -9e306 to 9e306 with a range of at most 9e306"
        )
        raise InputError(message)


def as_limits(name: str, value: object) -> tuple[float, float]:
    """Return `value` as two floats, low and high, with low < high.

    Each limit is inside -`AXIS_BOUND` to `AXIS_BOUND`, and high - low is at most `AXIS_BOUND`.
    A bool is not a number.
    """
    pair = _finite_pair(value)
    if pair is None or pair[0] >= pair[1] or not _inside_axis_bound(*pair):
        message = (
            f"{name} {shown(value)} is not valid; expected (low, high), "
            "two finite numbers from -9e306 to 9e306 with low < high and high - low at most 9e306"
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


def _too_small(figure: Figure) -> InputError:
    size = tuple(figure.get_size_inches().tolist())
    message = (
        f"size {shown(size)} is too small for the texts of the figure; "
        "expected a larger width or height"
    )
    return InputError(message)


def new_figure(size: tuple[float, float]) -> Figure:
    """Return an empty `Figure` of the size in inches. It does not use `matplotlib.pyplot`.

    The dpi is `matplotlib.rcParams["figure.dpi"]`. A dpi that is not finite and above 0 raises
    `InputError`, and so does a side below one pixel at the dpi.
    """
    dpi = float(mpl.rcParams["figure.dpi"])
    if not (math.isfinite(dpi) and dpi > 0):
        message = (
            f"figure dpi {dpi:g} is not valid; expected a finite figure dpi above 0, {_DPI_FROM}"
        )
        raise InputError(message)
    figure = Figure(figsize=size, dpi=dpi, layout="constrained", facecolor="white")
    # Such a side makes the transforms of matplotlib divide to infinity. Then the first artist
    # warns in np.dot on some platforms, and the layout raises a raw error.
    if figure.get_size_inches().min() * figure.dpi < 1:
        raise _too_small(figure)
    return figure


def check_text_pixels(figure: Figure) -> None:
    """Raise `InputError` if the font size of a drawn text is below 1 pixel at the figure dpi.

    A drawn text is visible and not empty. The function first fills the tick labels and the offset
    text of each axis, as a draw does.
    """
    # matplotlib does not draw an empty text, for example the default titles of an Axes. A draw
    # fills the tick labels and the offset text of an axis, thus an empty one can draw later.
    for axis in figure.findobj(Axis):
        axis.get_ticklabels(which="both")
        axis.get_offset_text().set_text(axis.get_major_formatter().get_offset())
    sizes = [
        text.get_fontsize()
        for text in figure.findobj(Text)
        if text.get_visible() and text.get_text()
    ]
    if not sizes:
        return
    # matplotlib makes each font size at least 1 point, thus the division is safe.
    smallest = float(min(sizes))
    # FreeType gets the dpi as an integer, and it rejects a size below half a pixel. The bound of
    # 1 pixel keeps a margin.
    if smallest * math.floor(figure.dpi) < _POINTS_PER_INCH:
        needed = math.ceil(_POINTS_PER_INCH / smallest)
        message = (
            f"figure dpi {figure.dpi:g} is too low for the font size {smallest:g} of a text of the "
            f"figure; expected a figure dpi of at least {needed}, {_DPI_FROM}"
        )
        raise InputError(message)


def apply_layout(figure: Figure) -> Figure:
    """Apply the constrained layout of the figure one time, and keep the positions of the Axes.

    The positions stay for each later draw. Thus a save at each dpi gives no layout warning. A size
    too small for the texts raises `InputError`, and so does a figure dpi too low for a text
    (`check_text_pixels`).
    """
    check_text_pixels(figure)
    # new_figure gives each figure of a plot function a constrained layout engine.
    engine = cast("LayoutEngine", figure.get_layout_engine())
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message=_COLLAPSE, category=UserWarning)
        try:
            engine.execute(figure)
        except UserWarning as warning:
            # An outer filter can also raise a warning that is not the collapse.
            if not str(warning).startswith(_COLLAPSE):
                raise
            raise _too_small(figure) from None
    figure.set_layout_engine("none")
    return figure


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
