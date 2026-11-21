"""The ECDF figure of the fold changes of gene groups (`.plot_ecdf` of R)."""

import math
from collections.abc import Iterable, Mapping, Sequence
from typing import TYPE_CHECKING, Final

import numpy as np
import numpy.typing as npt
import pandas as pd
from matplotlib.figure import Figure

from seedmatch._checks import as_items, as_text, check_columns, check_numeric_columns, check_type
from seedmatch._rmath import r_round
from seedmatch.errors import InputError, warn
from seedmatch.plot._style import (
    PALETTE,
    PANEL_EDGE,
    TEXT_COLOR,
    apply_layout,
    as_color,
    as_limits,
    as_size,
    as_written,
    check_axis_range,
    expanded,
    linewidth_points,
    new_figure,
    pasted,
    style_theme_axes,
)
from seedmatch.stats.ecdf import EcdfTestResult, assign_groups
from seedmatch.stats.wasserstein import ecdf

if TYPE_CHECKING:
    from matplotlib.lines import Line2D

_Y_VIEW: Final = (-0.05, 1.05)
_TWO_DECIMALS_FROM: Final = 0.01
_THREE_DECIMALS_FROM: Final = 0.001


# Each keyword-only argument is an independent option of the figure of `.plot_ecdf`.
def plot_ecdf(  # noqa: PLR0913
    de: pd.DataFrame,
    groups: Mapping[str, Iterable[str]],
    result: EcdfTestResult,
    *,
    title: str = "ECDF",
    palette: Sequence[str] = PALETTE,
    x_limits: tuple[float, float] = (-1.0, 1.0),
    size: tuple[float, float] = (5.0, 5.0),
) -> Figure:
    """Return the figure of `.plot_ecdf` of R: one ECDF line of the fold changes for each group.

    `de` needs the columns `gene_id` and `log2FoldChange`. `groups` maps a group name to gene IDs,
    as for `assign_groups`, and it needs at least one group. The figure has one line for each key
    of `groups`, in the order of the keys, and line i has the color `palette[i]`. A row in no group
    has no line (`D-064`). The values of a line are the finite `log2FoldChange` values of the rows
    in its group, and each group needs at least one. A group with a row that has no finite value
    warns with `SeedMatchWarning`, and the figure omits that row, as R does. The legend label of a
    line is its name and the count of its values (`D-064`). `result` is the result of
    `ecdf_test`, and the left title of the Axes gives its statistic and its p-value (`D-063`).
    `title` is the title of the figure. `palette` needs a color for each group. `x_limits` are the
    limits of the x axis before the expansion of 5% on each side.

    `size` is the width and the height in inches, each above 0 and below 50, as in `ggsave`. A size
    too small for the texts, the legends, and the labels raises `InputError`. The figure keeps its
    layout, thus a later change of the figure does not move the Axes. To apply the layout again,
    call `figure.set_layout_engine("constrained")`.

    A value that is not valid raises `InputError`. A value of `x_limits` or of the fold changes
    outside -9e306 to 9e306, or a range wider than 9e306, raises `InputError`.

    Each line is a step function, `drawstyle="steps-post"`, and it extends from the lowest to the
    highest of the view and of all values. The function returns a `Figure` and does not use
    `matplotlib.pyplot`. It writes no file, and it does not change `de` or `rcParams`.

    The time is O(n log n) for n rows, because `numpy.unique` sorts the values of each group. The
    figure holds two floats for each distinct value. For 22,000 rows in two groups, the call takes
    0.033 seconds, and the call with `Figure.savefig` to PNG at 100 dpi takes 0.07 seconds. The
    peak memory above the start of the two is 3.0 MiB. Machine: Apple M2 Max, Python 3.12.
    """
    check_columns(de, ("gene_id", "log2FoldChange"))
    check_numeric_columns(de, ["log2FoldChange"])
    if isinstance(groups, Mapping) and not groups:
        message = "groups has no group; expected at least one group"
        raise InputError(message)
    group_of_row = assign_groups(de, groups)
    check_type("result", result, EcdfTestResult)
    heading = as_text("title", title)
    colors = as_items("palette", palette, as_color, "str")
    names = [str(name) for name in group_of_row.cat.categories]
    if len(colors) < len(names):
        message = (
            f"palette has {len(colors)} color(s) for {len(names)} groups; "
            "expected a color for each group"
        )
        raise InputError(message)
    x_view = expanded(*as_limits("x_limits", x_limits))
    figure_size = as_size(size)
    samples, row_counts = _group_values(de, group_of_row, len(names))
    _check_samples(names, samples, row_counts)
    lowest, highest = _check_value_range(samples)
    left = min(x_view[0], lowest)
    right = max(x_view[1], highest)

    figure = new_figure(figure_size)
    axes = figure.add_subplot()
    lines: list[Line2D] = []
    labels: list[str] = []
    # The palette can have more colors than the figure has groups.
    for name, color, sample in zip(names, colors, samples, strict=False):
        curve = ecdf(sample)
        labels.append(f"{name}: {sample.size:,}")
        (line,) = axes.plot(
            np.concatenate(([left], curve.knots, [right])),
            np.concatenate(([0.0], curve.heights, [1.0])),
            drawstyle="steps-post",
            color=color,
            linewidth=linewidth_points(1.5),
            label=labels[-1],
        )
        lines.append(line)
    axes.set_xlim(*x_view)
    axes.set_ylim(*_Y_VIEW)
    axes.set_xlabel("Log2(Fold Change)")
    axes.set_ylabel("Cumulative Fraction")
    style_theme_axes(axes)
    as_written(figure.suptitle(heading, fontsize=14, color=TEXT_COLOR))
    subtitle = (
        f"Statistic: {pasted(r_round(result.statistic, 3))}\n"
        f"p-value: {_p_value_label(result.p_value)}"
    )
    axes.set_title(subtitle, loc="left", fontsize=12, color=TEXT_COLOR)
    # Without the labels, matplotlib leaves out the entry of a label that starts with an underscore.
    legend = axes.legend(
        lines,
        labels,
        title="Group",
        loc="upper left",
        fontsize=12,
        labelcolor=TEXT_COLOR,
        facecolor="white",
        edgecolor=PANEL_EDGE,
        framealpha=1.0,
    )
    legend.get_title().set(fontsize=12, fontweight="bold", color=TEXT_COLOR)
    for text in legend.get_texts():
        as_written(text)
    return apply_layout(figure)


def _p_value_label(p_value: float | None) -> str:
    """Return the p-value text of the subtitle. R stops for NaN and for 0.01 (`D-063`)."""
    if p_value is None:
        return "Undefined"
    if math.isnan(p_value):
        return "NaN"
    if p_value >= _TWO_DECIMALS_FROM:
        return pasted(r_round(p_value, 2))
    if p_value >= _THREE_DECIMALS_FROM:
        return pasted(r_round(p_value, 3))
    return "p<0.001"


def _group_values(
    de: pd.DataFrame, group_of_row: pd.Series, count: int
) -> tuple[list[npt.NDArray[np.float64]], npt.NDArray[np.intp]]:
    """Return the finite `log2FoldChange` values and the number of rows of each group.

    Both come in the order of the groups.
    """
    values = de["log2FoldChange"].to_numpy(dtype=np.float64, na_value=np.nan)
    codes = group_of_row.cat.codes.to_numpy()
    in_group = codes >= 0
    kept = in_group & np.isfinite(values)
    kept_codes = codes[kept]
    order = np.argsort(kept_codes, kind="stable")
    ordered = values[kept][order]
    bounds = np.searchsorted(kept_codes[order], np.arange(count + 1))
    samples = [ordered[bounds[i] : bounds[i + 1]] for i in range(count)]
    return samples, np.bincount(codes[in_group], minlength=count)


def _check_value_range(samples: list[npt.NDArray[np.float64]]) -> tuple[float, float]:
    """Return the lowest and the highest drawn value. Raise `InputError` if no x axis holds them."""
    lowest = min(float(sample.min()) for sample in samples)
    highest = max(float(sample.max()) for sample in samples)
    check_axis_range("log2FoldChange", lowest, highest)
    return lowest, highest


def _check_samples(
    names: list[str], samples: list[npt.NDArray[np.float64]], row_counts: npt.NDArray[np.intp]
) -> None:
    """Raise `InputError` for a group without a value, and warn for each group that loses rows."""
    for name, sample in zip(names, samples, strict=True):
        if sample.size == 0:
            message = (
                f"group {name!r} has no row with a finite log2FoldChange in the DE table; "
                "expected at least 1 row"
            )
            raise InputError(message)
    for name, sample, rows in zip(names, samples, row_counts, strict=True):
        if rows > sample.size:
            warn(
                f"group {name!r} has {rows - sample.size} row(s) without a finite "
                "log2FoldChange; the figure omits them"
            )
