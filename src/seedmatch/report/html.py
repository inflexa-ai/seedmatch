"""The HTML of the table that `SeedMatchReport` of R builds with `gt` (`.make_table`)."""

import html
import math
from typing import Final, NamedTuple

import numpy as np

from seedmatch._checks import (
    check_columns,
    check_numeric_columns,
    check_text_columns,
    check_type,
)
from seedmatch.report.report import REPORT_CATEGORIES, SeedMatchReport
from seedmatch.seed import Guide, Seed, seed_window

_TITLE: Final = "In-silico siRNA Binding Prediction"
_SUBTITLE: Final = "Identifying siRNA hits in the transcriptome"
_PERCENT_ROW: Final = "% off-target"
_DISTANCES: Final = 5
# The minus sign of gt is U+2212. The result is ASCII, thus it is a character reference.
_MINUS: Final = "&#8722;"


class _Window(NamedTuple):
    seed: Seed
    label: str
    prefix: str


_WINDOWS: Final = (
    _Window(Seed.FULL, "Full Guide Strand", "Full"),
    _Window(Seed.MER18, "18-mer", "18mer"),
    _Window(Seed.MER15, "15-mer", "15mer"),
)
_SITES: Final = ("8mer", "7mer-m8", "7mer-A1", "6mer")
_COLUMNS: Final = (
    *(f"{window.prefix}-D{distance}" for window in _WINDOWS for distance in range(_DISTANCES)),
    *_SITES,
    "Total",
)
_WIDTH: Final = len(_COLUMNS) + 1

_CSS: Final = """
table.seedmatch-report { border-collapse: collapse; color: #000000; font-size: 14px; }
table.seedmatch-report th, table.seedmatch-report td { padding: 4px 8px; }
table.seedmatch-report td { text-align: right; font-variant-numeric: tabular-nums; }
table.seedmatch-report thead th { text-align: center; font-weight: bold; }
table.seedmatch-report th.title { text-align: left; font-size: 125%; }
table.seedmatch-report th.subtitle { text-align: left; font-size: 85%; font-weight: normal; }
table.seedmatch-report th.group { text-align: left; border-top: 2px solid #d3d3d3; }
table.seedmatch-report th.stub { text-align: left; font-weight: normal; }
table.seedmatch-report tr.striped { background-color: rgba(128, 128, 128, 0.05); }
"""


def _text(value: str) -> str:
    return html.escape(value).encode("ascii", "xmlcharrefreplace").decode("ascii")


def _background(percent: float) -> str:
    """Return the blend from white (0) to red (100) of a percent value, clipped to 0 to 100."""
    shade = round(2.55 * (100 - min(max(percent, 0.0), 100.0)))
    return f"#FF{shade:02X}{shade:02X}"


def _cell_text(value: float, *, percent: bool) -> str:
    """Return the text of `fmt_number` (no decimals) or `fmt_percent` (two decimals) of `gt`."""
    spec = ",.2f" if percent else ",.0f"
    suffix = "%" if percent else ""
    if 0 <= value < math.inf:
        # abs gives -0.0 no sign, as gt does: gt adds the sign only for a value below 0.
        return format(abs(value), spec) + suffix
    if math.isnan(value):
        return "NaN"
    if math.isinf(value):
        digits, negative = "Inf", value < 0
    else:
        digits = format(-value, spec)
        # gt drops the sign of a count that rounds to 0, but it keeps the sign of a percent.
        negative = percent or digits != "0"
    return f"{_MINUS if negative else ''}{digits}{suffix}"


def _cell(value: float, *, percent: bool, colored: bool) -> str:
    text = _cell_text(value, percent=percent)
    if percent and colored and math.isfinite(value):
        return f'<td style="background-color: {_background(value)}">{text}</td>'
    return f"<td>{text}</td>"


def _row(name: str, values: list[float], *, striped: bool) -> str:
    percent = name == _PERCENT_ROW
    *categories, total = values
    cells = [_cell(value, percent=percent, colored=True) for value in categories]
    cells.append(_cell(total, percent=percent, colored=False))
    start = '<tr class="striped">' if striped else "<tr>"
    return f'{start}<th scope="row" class="stub">{_text(name)}</th>{"".join(cells)}</tr>'


def _header(guide: Guide) -> list[str]:
    length = len(guide.sequence)
    spanners = []
    for window in _WINDOWS:
        first, last = seed_window(window.seed, length)
        spanners.append(
            f'<th colspan="{_DISTANCES}" scope="colgroup">{window.label} (g{first}:g{last})</th>'
        )
    labels = "".join(f'<th scope="col">D{distance}</th>' for distance in range(_DISTANCES))
    sites = "".join(f'<th rowspan="2" scope="col">{name}</th>' for name in (*_SITES, "Total"))
    return [
        "<thead>",
        f'<tr><th class="title" colspan="{_WIDTH}">{_TITLE}</th></tr>',
        f'<tr><th class="subtitle" colspan="{_WIDTH}">{_SUBTITLE}</th></tr>',
        f'<tr><th rowspan="2"></th>{"".join(spanners)}{sites}</tr>',
        f"<tr>{labels * len(_WINDOWS)}</tr>",
        "</thead>",
    ]


def report_html(report: SeedMatchReport, guide: Guide) -> str:
    """Return the HTML of the `gt` table of `SeedMatchReport` of R: a `<style>` and a `<table>`.

    `guide` must be the guide of the call of `seed_match_report` that made `report`. The labels of
    the windows come from `guide`. R gives fixed texts, which are wrong for the 15-mer window
    and for a guide of another length than 23 bases (`D-070`). The columns, the title, and the
    subtitle are those of `gt`, in its order. The table has the class `seedmatch-report` and no
    `id`, and each rule of the style selects this class.

    Each group of `report.report`, in the order of its first row, has a `<tbody>`, and the even
    rows of a group have the class `striped`. A cell of the row `% off-target` has two decimals and
    a percent sign. Another cell has no decimals. Both have a comma between thousands.

    The texts of `NaN`, `Inf`, and a negative value are those of `gt`. The minus sign is a
    character reference. The port omits the space that `gt` writes before an `Inf` in a count
    column that also has `-Inf`. A pandas missing value has the text of `NaN`, but `gt` writes
    `NA` for an NA of R. A cell that is not finite has no color.

    The other category cells of the row `% off-target` have a linear blend from white (0) to red
    (100), with black text. `gt` blends in another color space and gives white text on a dark
    cell. 0 and 100 give the colors of `gt`.

    `report.report` must have the columns of `seed_match_report`. Each column of counts must be
    numeric, and each `category` and `group` value must be a `str`. A value that is not valid
    raises `InputError`. Each text is escaped, and each character outside ASCII is a character
    reference, thus the result is ASCII. The time and the memory are O(R) for R rows of
    `report.report`. For 10^5 rows, the call takes 1.0 seconds, and the peak of `tracemalloc` is
    119 MiB, for a result of 46 MiB.
    """
    check_type("report", report, SeedMatchReport)
    check_type("guide", guide, Guide)
    frame = report.report
    check_columns(
        frame,
        ("category", *REPORT_CATEGORIES, "Total", "group"),
        name="report.report",
        label="report table",
    )
    check_numeric_columns(frame, (*REPORT_CATEGORIES, "Total"))
    check_text_columns(frame, ("category", "group"), label="report table")
    matrix = frame[list(_COLUMNS)].to_numpy(dtype=np.float64, na_value=np.nan)
    names = frame["category"].tolist()
    rows_of_group: dict[str, list[int]] = {}
    for position, group in enumerate(frame["group"].tolist()):
        rows_of_group.setdefault(group, []).append(position)
    parts = [f"<style>{_CSS}</style>", '<table class="seedmatch-report">', *_header(guide)]
    for group, positions in rows_of_group.items():
        parts.append("<tbody>")
        parts.append(
            f'<tr><th class="group" colspan="{_WIDTH}" scope="colgroup">{_text(group)}</th></tr>'
        )
        parts.extend(
            _row(names[position], matrix[position].tolist(), striped=number % 2 == 0)
            for number, position in enumerate(positions, start=1)
        )
        parts.append("</tbody>")
    parts.extend(["</table>", ""])
    return "\n".join(parts)
