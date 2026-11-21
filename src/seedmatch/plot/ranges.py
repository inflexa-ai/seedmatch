"""The bars of the ranges of a table (`plotIRanges` and `plotRangesDF` of R)."""

from typing import Final

import numpy as np
import pandas as pd
from matplotlib.collections import PolyCollection
from matplotlib.figure import Figure

from seedmatch._checks import (
    as_number,
    as_text,
    check_columns,
    check_integer_columns,
    shown,
)
from seedmatch.errors import InputError
from seedmatch.plot._style import (
    as_color,
    as_limits,
    as_size,
    as_written,
    expanded,
    new_figure,
    style_gray_axes,
)

# A bar edge is a position plus or minus 0.5, and a float above 2^52 cannot hold the half. At
# exactly +-2^52, one edge rounds, and the bar is 0.5 wide.
_MAX_POSITION: Final = 2**52


# Each keyword-only argument is an option of `plotIRanges` or of `plotRangesDF`.
def plot_ranges(  # noqa: PLR0913
    ranges: pd.DataFrame,
    *,
    x_limits: tuple[float, float] | None = None,
    title: str = "Ranges",
    color: str = "blue",
    bar_height: float = 0.25,
    size: tuple[float, float] = (7.0, 5.0),
) -> Figure:
    """Return the bars of the ranges of `plotIRanges` and `plotRangesDF` of R, one for each row.

    `ranges` is a table with the integer columns `start` and `end`, for example a result of
    `find_matches`. Another column and the labels of the rows have no effect. Row i, from 1 in the
    order of the table, gets a bar from `start - 0.5` to `end + 0.5`, at the height i. Thus a range
    of one base has the width 1. R draws only the first 6 rows of a data.frame, at the positions of
    its row names, and it draws a range of one base with the width 0 (`D-068`). A value that is not
    a table, a missing column, a column that is not an integer column, a table without a row, a
    `start` or an `end` outside -2^52 to 2^52, and a row with `end < start` raise `InputError`.

    `x_limits` is the range of the x axis. Its default is the range of the bars with 5% on each
    side. `title` is the title of the Axes, `color` is the color of the bars, and `size` is the
    width and the height in inches. `bar_height` is the height of a bar, above 0 and at most 1.
    The value 1 is the distance between two rows, and a higher bar covers the bar of the next row.
    A value that is not valid raises `InputError`. A limit outside -9e306 to 9e306, or an
    `x_limits` range wider than 9e306, raises `InputError`. The y axis goes from 0 to n + 1 for n
    rows.

    The function returns a `Figure` with one Axes and one `PolyCollection`. It does not use
    `matplotlib.pyplot`. It writes no file, and it does not change `ranges` or `rcParams`.

    The time and the memory are O(n) for n rows. One `PolyCollection` holds the bars with 4
    vertices for each row, and no row has an artist of its own. The function has no row limit,
    because a limit drops rows without a message, as `head()` of R does. For 10^5 ranges, the call
    takes 0.065 seconds, 0.65 µs for each row. With `Figure.savefig` to PNG at 100 dpi, it takes
    0.24 seconds, 2.4 µs for each row. The peak of `tracemalloc` is 40.7 MiB, 427 bytes for each
    row. Machine: Apple M2 Max, Python 3.12.
    """
    check_columns(ranges, ("start", "end"), name="ranges", label="ranges table")
    check_integer_columns(ranges, ["start", "end"], label="ranges table")
    count = len(ranges)
    if count == 0:
        message = "ranges table has no row; expected at least one range"
        raise InputError(message)
    for column in ("start", "end"):
        values = ranges[column]
        if int(values.min()) < -_MAX_POSITION or int(values.max()) > _MAX_POSITION:
            outside = (values < -_MAX_POSITION) | (values > _MAX_POSITION)
            row = int(np.argmax(outside.to_numpy(dtype=bool)))
            message = (
                f"ranges table has {column} {int(values.iloc[row])} in row {row}; "
                "expected a position from -2^52 to 2^52"
            )
            raise InputError(message)
    reversed_rows = np.flatnonzero((ranges["end"] < ranges["start"]).to_numpy(dtype=bool))
    if reversed_rows.size:
        row = int(reversed_rows[0])
        message = (
            f"ranges table has start {int(ranges['start'].iloc[row])} "
            f"and end {int(ranges['end'].iloc[row])} in row {row}; expected end >= start"
        )
        raise InputError(message)
    limits = None if x_limits is None else as_limits("x_limits", x_limits)
    heading = as_text("title", title)
    face = as_color("color", color)
    height = as_number("bar_height", bar_height)
    if not 0 < height <= 1:
        message = (
            f"bar_height {shown(bar_height)} is out of range; "
            "expected a number above 0 and at most 1"
        )
        raise InputError(message)
    figure_size = as_size(size)

    left = ranges["start"].to_numpy(dtype=np.float64) - 0.5
    right = ranges["end"].to_numpy(dtype=np.float64) + 0.5
    centres = np.arange(1, count + 1, dtype=np.float64)
    vertices = np.empty((count, 4, 2))
    vertices[:, [0, 3], 0] = left[:, None]
    vertices[:, [1, 2], 0] = right[:, None]
    vertices[:, [0, 1], 1] = (centres - height / 2)[:, None]
    vertices[:, [2, 3], 1] = (centres + height / 2)[:, None]

    figure = new_figure(figure_size)
    axes = figure.add_subplot()
    # PolyCollection reads a 3-D array through a path that shares the codes of all bars. For 10^5
    # rows, the call takes 0.07 seconds and 41 MiB, against 0.29 seconds and 65 MiB with a list of
    # arrays. The stub names a Sequence.
    bars = PolyCollection(
        vertices,  # type: ignore[arg-type]
        facecolors=face,
        edgecolors=face,
        linewidths=0,
    )
    axes.add_collection(bars, autolim=False)
    axes.set_xlim(*(limits or expanded(float(left.min()), float(right.max()))))
    axes.set_ylim(0, count + 1)
    as_written(axes.set_title(heading))
    axes.set_xlabel("Position")
    axes.set_ylabel("Ranges")
    style_gray_axes(axes)
    return figure
