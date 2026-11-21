"""The alignment of a guide with a seed match as a figure (`plotAlignment` of R, ggplot output)."""

from typing import Final

from matplotlib.figure import Figure

from seedmatch._checks import as_text
from seedmatch.errors import InputError
from seedmatch.plot._style import as_written, check_text_pixels, new_figure

_NAMES: Final = ("siRNA ", " " * 6, "target")
_FONT_SIZE: Final = 11
# DejaVu Sans Mono is the first font of `font.monospace` in the default rcParams, and it ships with
# matplotlib. Its bold face at 11 points advances about 0.09 inch for each character, so 0.1 inch
# leaves room.
_INCHES_PER_CHARACTER: Final = 0.1
_MARGIN: Final = 0.4
_HEIGHT: Final = 0.9


def plot_alignment(text: str) -> Figure:
    """Return the figure of `plotAlignment` of R (`output.type = "ggplot"`), the lines of a text.

    `text` is the result of `alignment_text`, three lines with a line end between them. A text
    with another number of lines raises `InputError`. The figure has no Axes and one `Text`. The
    text is the three lines in order, each after its name and one space. The names are `siRNA `,
    six spaces, and `target` (`D-067`). R names the first line `miRNA` in the figure, but the
    data.frame names it `siRNA `. The font is bold and monospace, with the size 11 of `geom_text`,
    and the block of text is at the center of the figure. A `$` in `text` draws as written. A
    figure dpi too low for the text raises `InputError`.

    The label is the text with the names. The width of the figure is 0.1 inch for each character
    of the longest line of the label, plus 0.4 inch. The height is 0.9 inch. The face is white.
    The function returns a `Figure` and does not use `matplotlib.pyplot`. It writes no file, and
    it does not change `rcParams`.

    The time and the memory are O(W) for a label whose longest line has W characters.
    `Figure.savefig` at 100 dpi draws a canvas of 3,600 bytes for each character of W. For the
    label of the first `mer7m8` match of the test data with `left_pad=20` and `right_pad=20`, W is
    84. The call takes 0.1 ms, and the call with `Figure.savefig` to PNG at 100 dpi takes 13 ms.
    For W = 2,000, the call with `Figure.savefig` takes 0.27 seconds. Machine: Apple M2 Max,
    Python 3.12.
    """
    lines = as_text("text", text).split("\n")
    if len(lines) != len(_NAMES):
        message = (
            f"text has {len(lines)} line(s); expected the {len(_NAMES)} lines of alignment_text"
        )
        raise InputError(message)
    rows = [f"{name} {line}" for name, line in zip(_NAMES, lines, strict=True)]
    label = "\n".join(rows)
    width = _INCHES_PER_CHARACTER * max(map(len, rows)) + _MARGIN
    figure = new_figure((width, _HEIGHT))
    as_written(
        figure.text(
            0.5,
            0.5,
            label,
            ha="center",
            va="center",
            multialignment="left",
            fontfamily="monospace",
            fontweight="bold",
            fontsize=_FONT_SIZE,
        )
    )
    check_text_pixels(figure)
    return figure
