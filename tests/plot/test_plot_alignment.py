import io
import re

import pytest
from helpers import GUIDE, exactly, golden, golden_note, sirna_seqs
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import to_hex

from seedmatch import Guide, Seed, alignment_text, find_matches
from seedmatch.errors import InputError
from seedmatch.plot import plot_alignment

ALIGNMENT_PLOT = "plots/alignment_plot.tsv"
NAMES = ["siRNA ", " " * 6, "target"]


def first_match_text(**pads):
    guide = Guide(GUIDE)
    match = find_matches(guide.target(Seed.MER7M8), sirna_seqs()).iloc[0]
    return alignment_text(match, guide, sirna_seqs(), **pads)


def test_plot_alignment_shows_the_r_lines_D067():
    r_names = re.findall(
        r"\[([^\]]*)\]", golden_note(ALIGNMENT_PLOT, r"^# data.frame row names: (.*)$")
    )
    assert r_names == NAMES
    r_lines = golden(ALIGNMENT_PLOT)["text"].tolist()
    assert r_lines[0].startswith("miRNA ")

    figure = plot_alignment(first_match_text())
    assert figure.axes == []
    (label,) = figure.texts
    lines = label.get_text().split("\n")
    assert lines == ["siRNA " + r_lines[0][len("miRNA ") :], r_lines[1], r_lines[2]]
    assert [line[:6] for line in lines] == r_names
    assert to_hex(figure.get_facecolor()) == "#ffffff"


def test_plot_alignment_uses_a_bold_monospace_font():
    (label,) = plot_alignment(first_match_text()).texts
    assert label.get_fontfamily() == ["monospace"]
    assert label.get_fontweight() == "bold"
    assert label.get_fontsize() == 11
    assert label.get_position() == (0.5, 0.5)
    assert label.get_horizontalalignment() == "center"
    assert label.get_verticalalignment() == "center"
    # Text has no public getter for the alignment of its lines.
    assert label._get_multialignment() == "left"  # noqa: SLF001


def test_plot_alignment_sizes_the_figure_by_its_longest_line():
    text = first_match_text()
    longest = max(
        len(name) + 1 + len(line) for name, line in zip(NAMES, text.split("\n"), strict=True)
    )
    assert tuple(plot_alignment(text).get_size_inches()) == pytest.approx(
        (0.1 * longest + 0.4, 0.9)
    )
    wide = first_match_text(left_pad=20, right_pad=20)
    assert tuple(plot_alignment(wide).get_size_inches()) == pytest.approx(
        (0.1 * (longest + 40) + 0.4, 0.9)
    )


@pytest.mark.parametrize("longest", [0, 1, 2])
def test_plot_alignment_takes_the_width_from_the_longest_line(longest):
    lengths = [4, 6, 5]
    lengths[longest] = 12
    text = "\n".join("x" * length for length in lengths)
    width = max(len(name) + 1 + length for name, length in zip(NAMES, lengths, strict=True))
    assert tuple(plot_alignment(text).get_size_inches()) == pytest.approx((0.1 * width + 0.4, 0.9))


def test_plot_alignment_fits_the_text():
    figure = plot_alignment(first_match_text(left_pad=20, right_pad=20))
    (label,) = figure.texts
    extent = label.get_window_extent(FigureCanvasAgg(figure).get_renderer())
    box = figure.bbox
    assert box.x0 <= extent.x0
    assert extent.x1 <= box.x1
    assert box.y0 <= extent.y0
    assert extent.y1 <= box.y1


@pytest.mark.parametrize("text", ["$x$\n$y\nz$", "cost $5 and $6\n$a$\nb", "$\\badcommand$\n.\n."])
def test_plot_alignment_draws_the_text_as_written(text):
    figure = plot_alignment(text)
    assert figure.texts[0].get_parse_math() is False
    figure.savefig(io.BytesIO(), format="png")


@pytest.mark.parametrize(
    ("text", "message"),
    [
        (5, "text 5 has type int; expected a str"),
        ("a\nb", "text has 2 line(s); expected the 3 lines of alignment_text"),
        ("a\nb\nc\nd", "text has 4 line(s); expected the 3 lines of alignment_text"),
        ("a\nb\nc\n", "text has 4 line(s); expected the 3 lines of alignment_text"),
        ("", "text has 1 line(s); expected the 3 lines of alignment_text"),
        # A bare CR is not a line end: `str.splitlines` would read this text as 3 lines.
        ("a\nb\rc", "text has 2 line(s); expected the 3 lines of alignment_text"),
    ],
)
def test_plot_alignment_rejects_bad_text(text, message):
    with pytest.raises(InputError, match=exactly(message)):
        plot_alignment(text)
