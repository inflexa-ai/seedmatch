"""The PR curve and the ROC curve of a benchmark (`plot.PRROC` of PRROC 1.4)."""

from typing import Final

import numpy as np
import numpy.typing as npt
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, Normalize, hsv_to_rgb
from matplotlib.figure import Figure

from seedmatch._checks import check_type, shown
from seedmatch.errors import InputError
from seedmatch.plot._style import apply_layout, as_size, new_figure
from seedmatch.stats import Benchmark

_COLORS: Final = 100
_LAST: Final = _COLORS - 1
# The scale of PRROC: `hsv(h = seq(0, 1, length = 100) * 0.8, s = 1, v = 1)`.
_SCALE: Final = hsv_to_rgb(
    np.column_stack((0.8 * np.arange(_COLORS) / _LAST, np.ones(_COLORS), np.ones(_COLORS)))
)
# The color of a curve without a range of thresholds (`D-069`).
_MIDDLE: Final = 49
# Base graphics of R: `lwd = 3` of `segments` is 3 * 0.75 points.
_LINE_WIDTH: Final = 2.25
# The default `xaxs = "r"` of R for the limits 0 and 1 adds 4% of the range on each side.
_VIEW: Final = (-0.04, 1.04)
# About the largest float divided by 198. Inside this bound, `99 * (high - low)` is finite, and so
# are the sums of two neighbor boundaries in the color bar of matplotlib.
_MAX_THRESHOLD: Final = 9e305


def _color_indexes(
    name: str, thresholds: npt.NDArray[np.float64]
) -> tuple[npt.NDArray[np.intp], tuple[float, float]]:
    """Return the index in `_SCALE` of each threshold, and the range of the color bar.

    R uses the extremes of the thresholds and stops for an infinite one. Here the range comes from
    the finite thresholds (`D-069`). Without two distinct finite thresholds, the range is c - 0.5 to
    c + 0.5. If c - 0.5 and c + 0.5 round to c, matplotlib widens the color bar by 10% of |c| on
    each side (0.9 c to 1.1 c for c above 0, and 1.1 c to 0.9 c for c below 0). A finite threshold
    outside -9e305 to 9e305 raises `InputError`.
    """
    finite = thresholds[np.isfinite(thresholds)]
    too_large = np.flatnonzero(np.abs(finite) > _MAX_THRESHOLD)
    if too_large.size:
        message = (
            f"the {name} curve has the finite threshold {shown(float(finite[too_large[0]]))}; "
            "expected a finite threshold from -9e305 to 9e305 for the colors"
        )
        raise InputError(message)
    if finite.size and finite.min() < finite.max():
        low, high = float(finite.min()), float(finite.max())
        position = np.clip((_LAST * (thresholds - low)) / (high - low), 0, _LAST)
        # R rounds `1 + position` half to even, and its index counts from 1. Rounding `position`
        # alone gives another color at a tie.
        return np.rint(1 + position).astype(np.intp) - 1, (low, high)
    centre = float(finite[0]) if finite.size else 0.0
    return np.full(thresholds.size, _MIDDLE, dtype=np.intp), (centre - 0.5, centre + 0.5)


def _draw_curve(
    axes: Axes, curve: pd.DataFrame, columns: tuple[str, str], indexes: npt.NDArray[np.intp]
) -> None:
    """Draw a segment from each row of the curve to the next row, with its color from `indexes`."""
    points = curve[list(columns)].to_numpy(dtype=np.float64)
    segments = np.stack((points[:-1], points[1:]), axis=1)
    lines = LineCollection(list(segments), colors=_SCALE[indexes[:-1]], linewidths=_LINE_WIDTH)
    axes.add_collection(lines, autolim=False)
    axes.set_xlim(*_VIEW)
    axes.set_ylim(*_VIEW)


def plot_benchmark(benchmark: Benchmark, *, size: tuple[float, float] = (12.0, 5.0)) -> Figure:
    """Return the PR curve and the ROC curve of a `Benchmark`, as `plot.PRROC` of R draws them.

    `benchmark` is the result of `benchmark_predictions`. The figure has the PR Axes and the ROC
    Axes in one row, with a color bar to the right of each. `Figure.axes` lists the PR Axes, the ROC
    Axes, and then the two color bars. Each Axes has one `LineCollection` with a segment from each
    row of its curve to the next row, in the order of the rows. The points are `recall` and
    `precision` for the PR curve, and `fpr` and `tpr` for the ROC curve. The line width is 2.25
    points, the width of `lwd = 3`. The view is -0.04 to 1.04 on each axis, as in R.

    A segment has the color of its first row, as in R. The color is one of the 100 colors of the
    scale of PRROC, from red to violet. The position on the scale is the threshold of the row, from
    the lowest to the highest finite threshold of the curve. The color bar shows this scale.
    Where R stops, the port draws a figure (`D-069`). An infinite threshold gets the color of the
    nearest end of the scale. A curve without two distinct finite thresholds gets the middle color,
    and its color bar goes from c - 0.5 to c + 0.5 for the finite threshold c, or for 0 if there
    is none. If c - 0.5 and c + 0.5 round to c, matplotlib widens the color bar by 10% of |c| on
    each side (0.9 c to 1.1 c for c above 0, and 1.1 c to 0.9 c for c below 0). A finite threshold
    outside -9e305 to 9e305 raises `InputError` before the function draws. The message names the
    first such threshold in the order of the rows. Inside the bound, each position on the scale and
    each boundary of the color bar is finite.

    The titles give the area, with 7 significant digits as in R. The PR area is `auc_integral`,
    the default of `plot.PRROC`. For many negatives, it differs from PRROC (`D-050`). The Axes
    keep the default style of matplotlib, because the plots of R use base graphics.

    `size` is the width and the height in inches, each above 0 and below 50, as in `ggsave`. A size
    too small for the texts, the legends, and the labels raises `InputError`. The layout does not
    count the width of a title, thus the figure cuts a title that is wider than the figure. The
    figure keeps its layout, thus a later change of the figure does not move the Axes of the curves.
    A color bar keeps its aspect ratio, thus `set_size_inches` changes its drawn width. To apply the
    layout again, call `figure.set_layout_engine("constrained")`. A `benchmark` that is not a
    `Benchmark`, and a `size` that is not valid, also raise `InputError`.

    The function returns a `Figure`. It does not use `matplotlib.pyplot`. It writes no file, and
    it does not change `benchmark` or `rcParams`.

    The time and the memory are O(R) for R rows of the two curves. The PR curve has at most
    2K + 2 + max(100, P) rows for K distinct scores and P positives. The function has no row limit,
    because a limit drops rows without a message. For 10^5 distinct scores (200,003 rows), the
    call takes 0.36 seconds, 1.8 µs for each row. With `Figure.savefig` to PNG at 100 dpi, it takes
    0.71 seconds, 3.5 µs for each row. The peak of `tracemalloc` is 69.1 MiB, 362 bytes for each
    row, and 71.4 MiB, 374 bytes for each row, with `savefig`. For 10^6 scores of 8 values with
    30% positives (262,411 rows), the call takes 0.47 seconds, 1.8 µs for each row, and 0.89
    seconds with `savefig`, 3.4 µs for each row. The peak is 94.7 MiB, 379 bytes for each row, and
    98.3 MiB, 393 bytes for each row, with `savefig`. Machine: Apple M2 Max, Python 3.12.
    """
    check_type("benchmark", benchmark, Benchmark)
    figure_size = as_size(size)
    pr_indexes, pr_limits = _color_indexes(
        "PR", benchmark.pr.curve["threshold"].to_numpy(dtype=np.float64)
    )
    roc_indexes, roc_limits = _color_indexes(
        "ROC", benchmark.roc.curve["threshold"].to_numpy(dtype=np.float64)
    )
    figure = new_figure(figure_size)
    pr_axes = figure.add_subplot(1, 2, 1)
    roc_axes = figure.add_subplot(1, 2, 2)
    _draw_curve(pr_axes, benchmark.pr.curve, ("recall", "precision"), pr_indexes)
    _draw_curve(roc_axes, benchmark.roc.curve, ("fpr", "tpr"), roc_indexes)
    pr_axes.set_title(f"PR curve\nAUC = {benchmark.pr.auc_integral:.7g}")
    pr_axes.set_xlabel("Recall")
    pr_axes.set_ylabel("Precision")
    roc_axes.set_title(f"ROC curve\nAUC = {benchmark.roc.auc:.7g}")
    roc_axes.set_xlabel("FPR")
    roc_axes.set_ylabel("Sensitivity")
    for axes, limits in ((pr_axes, pr_limits), (roc_axes, roc_limits)):
        mappable = ScalarMappable(norm=Normalize(*limits), cmap=ListedColormap(_SCALE))
        figure.colorbar(mappable, ax=axes)
    return apply_layout(figure)
