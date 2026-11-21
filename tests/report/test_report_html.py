import dataclasses
import re
import time
import tracemalloc
from dataclasses import dataclass, field
from functools import cache
from html.parser import HTMLParser

import numpy as np
import pandas as pd
import pytest
from helpers import GUIDE, exactly, golden, golden_note, sirna_de, sirna_seqs

from seedmatch import Guide, SeedMatchReport, report_html, seed_match_report
from seedmatch.errors import InputError
from seedmatch.report import REPORT_CATEGORIES

TABLE = "plots/report_table.tsv"
ROW_NAMES = [
    "In silico predictions",
    "Expressed predictions",
    "Off-target predictions",
    "% off-target",
]
PERCENT_ROW = "% off-target"
MINUS = "\u2212"
SITE_COLUMNS = ["8mer", "7mer-m8", "7mer-A1", "6mer", "Total"]
# The backgrounds of the percent row, by the text of the cell. The port blends from white to red.
# gt gives #FFBFAA, #FFA68B, #FFC8B5, #FFD5C6, and #FF6846 for the same cells.
PORT_BACKGROUNDS = {
    "0.00%": "#FFFFFF",
    "100.00%": "#FF0000",
    "33.33%": "#FFAAAA",
    "45.83%": "#FF8A8A",
    "28.77%": "#FFB6B6",
    "22.22%": "#FFC6C6",
    "75.00%": "#FF4040",
}


@dataclass
class Cell:
    tag: str
    attrs: dict[str, str | None]
    text: str = ""


@dataclass
class Row:
    attrs: dict[str, str | None]
    cells: list[Cell] = field(default_factory=list)


class Page(HTMLParser):
    """The rows of `thead` and of each `tbody`, with the tag, attributes, and text of each cell."""

    def __init__(self, text):
        super().__init__()
        self.tags: list[tuple[str, dict[str, str | None]]] = []
        self.head: list[Row] = []
        self.bodies: list[list[Row]] = []
        self.css = ""
        self._rows: list[Row] | None = None
        self._row: Row | None = None
        self._cell: Cell | None = None
        self._in_style = False
        self.feed(text)
        self.close()

    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))
        if tag == "thead":
            self._rows = self.head
        elif tag == "tbody":
            self._rows = []
            self.bodies.append(self._rows)
        elif tag == "tr":
            assert self._rows is not None
            self._row = Row(dict(attrs))
            self._rows.append(self._row)
        elif tag in ("th", "td"):
            assert self._row is not None
            self._cell = Cell(tag, dict(attrs))
            self._row.cells.append(self._cell)
        elif tag == "style":
            self._in_style = True

    def handle_endtag(self, tag):
        if tag in ("th", "td"):
            self._cell = None
        elif tag == "style":
            self._in_style = False

    def handle_data(self, data):
        if self._in_style:
            self.css += data
        elif self._cell is not None:
            self._cell.text += data

    def data_rows(self, body=0):
        """Return the rows of a `tbody` after its group row."""
        group, *rows = self.bodies[body]
        assert [cell.attrs["class"] for cell in group.cells] == ["group"]
        return rows


def rules(css):
    """Return the (selectors, declarations) of each rule of the CSS."""
    found = []
    for block in css.split("}"):
        if "{" in block:
            selectors, declarations = block.split("{")
            found.append(([part.strip() for part in selectors.split(",")], declarations))
    return found


@cache
def default_report() -> SeedMatchReport:
    return seed_match_report(sirna_de(), sirna_seqs(), Guide(GUIDE))


@cache
def default_text() -> str:
    return report_html(default_report(), Guide(GUIDE))


@cache
def default_page() -> Page:
    return Page(default_text())


def html_of(frame, guide=GUIDE):
    return report_html(dataclasses.replace(default_report(), report=frame), Guide(guide))


def frame_copy():
    return default_report().report.copy()


def golden_columns():
    return golden_note(TABLE, r"columns: (.+)").split(",")


def table_cells(page):
    """Return (row, column, cell) for each data cell of the first `tbody`, row by row."""
    columns = golden_columns()
    return [
        (row.cells[0].text, column, cell)
        for row in page.data_rows()
        for column, cell in zip(columns, row.cells[1:], strict=True)
    ]


def cell_of(page, row, column):
    (cell,) = (cell for r, c, cell in table_cells(page) if (r, c) == (row, column))
    return cell


def gt_text(key):
    """Return the text of a note. The note writes the minus sign U+2212 as the text `U+2212`."""
    return golden_note(TABLE, rf"^# {re.escape(key)}: (.*)$").replace("U+2212", MINUS)


def mixed_frame():
    """Return the report frame with the rows of the groups B and A alternating, B first."""
    both = pd.concat([frame_copy().assign(group="B"), frame_copy().assign(group="A")])
    return both.sort_index(kind="stable").reset_index(drop=True)


def test_report_html_gives_the_r_cell_texts():
    got = [(row, column, cell.text) for row, column, cell in table_cells(default_page())]
    expected = list(golden(TABLE)[["row", "column", "text"]].itertuples(index=False, name=None))
    assert got == expected
    assert (
        golden_columns()
        == [
            f"{window}-D{distance}"
            for window in ("Full", "18mer", "15mer")
            for distance in range(5)
        ]
        + SITE_COLUMNS
    )
    assert cell_of(default_page(), "In silico predictions", "Total").text == "125"
    assert cell_of(default_page(), PERCENT_ROW, "18mer-D4").text == "45.83%"
    assert cell_of(default_page(), PERCENT_ROW, "Total").text == "33.60%"


def test_report_html_gives_the_r_headers():
    page = default_page()
    colspan = golden_note(TABLE, r"colspan: (\d+)")
    title, subtitle, spanners, labels = page.head
    assert title.cells == [
        Cell("th", {"class": "title", "colspan": colspan}, golden_note(TABLE, r"title: (.+)"))
    ]
    assert subtitle.cells == [
        Cell("th", {"class": "subtitle", "colspan": colspan}, golden_note(TABLE, r"subtitle: (.+)"))
    ]
    stub, *windows = spanners.cells[:4]
    assert stub == Cell("th", {"rowspan": "2"}, "")
    assert [cell.text for cell in windows[:2]] == [
        golden_note(TABLE, r"spanner Full: label (.+?) / columns"),
        golden_note(TABLE, r"spanner 18mer: label (.+?) / columns"),
    ]
    assert all(
        cell.tag == "th" and cell.attrs == {"colspan": "5", "scope": "colgroup"} for cell in windows
    )
    assert spanners.cells[4:] == [
        Cell("th", {"rowspan": "2", "scope": "col"}, name) for name in SITE_COLUMNS
    ]
    assert labels.cells == [Cell("th", {"scope": "col"}, f"D{d}") for d in range(5)] * 3
    group = page.bodies[0][0]
    assert group.cells == [
        Cell(
            "th",
            {"class": "group", "colspan": colspan, "scope": "colgroup"},
            golden_note(TABLE, r"group: (.+)"),
        )
    ]
    assert len(page.bodies) == 1
    assert [row.cells[0].text for row in page.data_rows()] == ROW_NAMES
    assert all(
        row.cells[0].tag == "th" and row.cells[0].attrs == {"scope": "row", "class": "stub"}
        for row in page.data_rows()
    )
    assert all(
        cell.tag == "td" and "colspan" not in cell.attrs for cell in page.data_rows()[0].cells[1:]
    )
    assert len(page.data_rows()[0].cells) == int(colspan)


def test_report_html_is_one_style_and_one_table():
    page = default_page()
    tags = [tag for tag, _ in page.tags]
    assert tags.count("style") == 1
    assert tags.count("table") == 1
    assert [attrs for tag, attrs in page.tags if tag == "table"] == [{"class": "seedmatch-report"}]
    assert default_text().isascii()
    (base,) = (
        declarations
        for selectors, declarations in rules(page.css)
        if selectors == ["table.seedmatch-report"]
    )
    assert "color: #000000;" in base


def test_report_html_labels_the_windows_from_the_guide_D070():
    spanners = default_page().head[2].cells[1:4]
    assert [cell.text for cell in spanners] == [
        golden_note(TABLE, r"spanner Full: label (.+?) / columns"),
        golden_note(TABLE, r"spanner 18mer: label (.+?) / columns"),
        "15-mer (g2:g16)",
    ]
    # R labels the window of the 15mer columns g2:g19, but the window is g2 to g16.
    assert golden_note(TABLE, r"spanner 15mer: label (.+?) / columns") == "15-mer (g2:g19)"
    short = GUIDE[:21]
    report = seed_match_report(sirna_de(), sirna_seqs(), Guide(short))
    page = Page(report_html(report, Guide(short)))
    assert [cell.text for cell in page.head[2].cells[1:4]] == [
        "Full Guide Strand (g2:g21)",
        "18-mer (g2:g19)",
        "15-mer (g2:g16)",
    ]
    assert (
        golden_note(TABLE, r"spanner Full of a 21-base guide: (.+)") == "Full Guide Strand (g2:g23)"
    )


def test_report_html_has_no_fixed_id_D070():
    page = default_page()
    assert "id=" not in default_text()
    assert all("id" not in attrs for _, attrs in page.tags)
    assert golden_note(TABLE, r"id: (\w+)") == "mygt"
    assert "mygt" not in default_text()
    selectors = [selector for group, _ in rules(page.css) for selector in group]
    assert selectors
    assert all(selector.startswith("table.seedmatch-report") for selector in selectors)
    assert "#" not in "".join(selectors)


def test_report_html_colors_the_percent_row():
    page = default_page()
    colored = 0
    for row, column, background in golden(TABLE)[["row", "column", "background"]].itertuples(
        index=False, name=None
    ):
        cell = cell_of(page, row, column)
        if background == "NA":
            assert "style" not in cell.attrs, (row, column)
            continue
        colored += 1
        assert cell.attrs == {"style": f"background-color: {PORT_BACKGROUNDS[cell.text]}"}, column
        if cell.text in ("0.00%", "100.00%"):
            assert PORT_BACKGROUNDS[cell.text] == background
    assert colored == len(REPORT_CATEGORIES)
    assert "style" not in cell_of(page, PERCENT_ROW, "Total").attrs


def test_report_html_clips_the_percent_color():
    frame = frame_copy()
    row = frame.index[frame["category"] == PERCENT_ROW][0]
    frame.loc[row, ["Full-D0", "Full-D1"]] = [150.0, -5.0]
    page = Page(html_of(frame))
    high, low = (cell_of(page, PERCENT_ROW, name) for name in ("Full-D0", "Full-D1"))
    assert (high.text, high.attrs["style"]) == ("150.00%", "background-color: #FF0000")
    assert (low.text, low.attrs["style"]) == (
        gt_text("fmt_percent -5"),
        "background-color: #FFFFFF",
    )


def test_report_html_rounds_the_percent_color_of_a_half_to_even():
    # 2.55 * 70 is 178.5 and 2.55 * 30 is 76.5. Round half up gives #FFB3B3 and #FF4D4D.
    frame = frame_copy()
    row = frame.index[frame["category"] == PERCENT_ROW][0]
    frame.loc[row, ["Full-D0", "Full-D1"]] = [30.0, 70.0]
    page = Page(html_of(frame))
    got = [cell_of(page, PERCENT_ROW, name) for name in ("Full-D0", "Full-D1")]
    assert [(cell.text, cell.attrs["style"]) for cell in got] == [
        ("30.00%", "background-color: #FFB2B2"),
        ("70.00%", "background-color: #FF4C4C"),
    ]


# The texts of gt 1.4.0 are the notes `fmt_number <value>` and `fmt_percent <value>` of the golden
# file. gt writes the minus sign U+2212, and it writes a space before an Inf in a count column that
# also has -Inf, which the port omits. A count that rounds to 0 has no sign, but a percent keeps
# its sign. gt gives no sign to -0.
SPECIAL_VALUES = [
    pytest.param(ROW_NAMES[0], np.inf, "fmt_number Inf", None, id="count-inf"),
    pytest.param(ROW_NAMES[0], -np.inf, "fmt_number -Inf", None, id="count-minus-inf"),
    pytest.param(ROW_NAMES[0], -0.4, "fmt_number -0.4", None, id="count-rounds-to-zero"),
    pytest.param(ROW_NAMES[0], -5.0, "fmt_number -5", None, id="count-negative"),
    pytest.param(ROW_NAMES[0], -1234.4, "fmt_number -1234.4", None, id="count-negative-with-comma"),
    pytest.param(ROW_NAMES[0], -0.0, "fmt_number -0", None, id="count-negative-zero-value"),
    pytest.param(PERCENT_ROW, np.inf, "fmt_percent Inf", None, id="percent-inf"),
    pytest.param(PERCENT_ROW, -np.inf, "fmt_percent -Inf", None, id="percent-minus-inf"),
    pytest.param(PERCENT_ROW, -0.4, "fmt_percent -0.4", "#FFFFFF", id="percent-negative"),
    pytest.param(PERCENT_ROW, -5.0, "fmt_percent -5", "#FFFFFF", id="percent-negative-whole"),
    pytest.param(PERCENT_ROW, 1e-9, "fmt_percent 1e-9", "#FFFFFF", id="percent-rounds-to-zero"),
    pytest.param(PERCENT_ROW, -0.004, "fmt_percent -0.004", "#FFFFFF", id="percent-negative-zero"),
    pytest.param(PERCENT_ROW, 12345.678, "fmt_percent 12345.678", "#FF0000", id="percent-comma"),
    pytest.param(PERCENT_ROW, -0.0, "fmt_percent -0", "#FFFFFF", id="percent-negative-zero-value"),
]


@pytest.mark.parametrize(("row", "value", "key", "background"), SPECIAL_VALUES)
def test_report_html_gives_the_gt_text_of_an_infinite_or_negative_cell(row, value, key, background):
    frame = frame_copy()
    frame.loc[frame["category"] == row, "Full-D0"] = value
    page_text = html_of(frame)
    cell = cell_of(Page(page_text), row, "Full-D0")
    text = gt_text(key)
    assert cell.text == text
    # A cell without a color has no style. A finite percent has the clipped color.
    expected = {} if background is None else {"style": f"background-color: {background}"}
    assert cell.attrs == expected
    assert page_text.isascii()
    assert ("&#8722;" in page_text) == (MINUS in text)


@pytest.mark.parametrize(("dtype", "missing"), [("float64", np.nan), ("Float64", pd.NA)])
def test_report_html_gives_a_missing_cell_the_text_nan_and_no_style(dtype, missing):
    frame = frame_copy().astype(dict.fromkeys(["Full-D0", "6mer"], dtype))
    percent = frame.index[frame["category"] == PERCENT_ROW][0]
    frame.loc[percent, "Full-D0"] = missing
    frame.loc[0, "6mer"] = missing
    page = Page(html_of(frame))
    # gt gives NaN for a NaN, but NA for an NA (the notes `fmt_number NA` and `fmt_percent NA`).
    # A pandas frame has one missing value, and the port gives it the text of the NaN.
    for row, column, key in (
        (PERCENT_ROW, "Full-D0", "fmt_percent NaN"),
        (ROW_NAMES[0], "6mer", "fmt_number NaN"),
    ):
        cell = cell_of(page, row, column)
        assert (cell.text, cell.attrs) == (gt_text(key), {})
    assert cell_of(page, PERCENT_ROW, "Full-D1").text == "100.00%"


def test_report_html_stripes_rows_2_and_4():
    page = Page(html_of(mixed_frame()))
    positions = []
    for body in range(2):
        classes = [row.attrs.get("class") for row in page.data_rows(body)]
        assert classes == [None, "striped", None, "striped"]
        positions.append(",".join(str(i) for i, name in enumerate(classes, start=1) if name))
    assert positions == [golden_note(TABLE, r"striped rows: (.+)")] * 2


def test_report_html_stripes_each_group_from_its_own_first_row():
    # A has 3 rows and B has 6 rows. A stripe count over all groups would start B on a stripe.
    frame = frame_copy()
    order = [
        (0, "A"),
        (1, "A"),
        (0, "B"),
        (1, "B"),
        (2, "B"),
        (3, "B"),
        (0, "B"),
        (1, "B"),
        (2, "A"),
    ]
    rows = [frame.iloc[[position]].assign(group=group) for position, group in order]
    page = Page(html_of(pd.concat(rows, ignore_index=True)))
    assert [body[0].cells[0].text for body in page.bodies] == ["A", "B"]
    classes = [[row.attrs.get("class") for row in page.data_rows(body)] for body in range(2)]
    assert classes == [[None, "striped", None], [None, "striped", None, "striped", None, "striped"]]


def test_report_html_puts_each_group_in_its_own_body():
    page = Page(html_of(mixed_frame()))
    assert [body[0].cells[0].text for body in page.bodies] == ["B", "A"]
    for body in range(2):
        assert [row.cells[0].text for row in page.data_rows(body)] == ROW_NAMES


def test_report_html_puts_a_comma_in_a_large_count():
    frame = frame_copy()
    frame.loc[0, "Total"] = 12345.0
    assert cell_of(Page(html_of(frame)), ROW_NAMES[0], "Total").text == golden_note(
        TABLE, r"fmt_number 12345: (.+)"
    )
    assert golden_note(TABLE, r"fmt_number 12345: (.+)") == "12,345"


def test_report_html_gives_a_percent_with_two_decimals():
    frame = frame_copy()
    percent = frame.index[frame["category"] == PERCENT_ROW][0]
    frame.loc[percent, "Total"] = 0.125
    cell = cell_of(Page(html_of(frame)), PERCENT_ROW, "Total")
    assert cell.text == golden_note(TABLE, r"fmt_percent 0.125: (.+)")


def test_report_html_escapes_text():
    frame = frame_copy().assign(group="<b>&é")
    frame.loc[0, "category"] = "<i>"
    text = html_of(frame)
    assert "&lt;b&gt;&amp;&#233;" in text
    assert "&lt;i&gt;" in text
    assert text.isascii()
    page = Page(text)
    assert page.bodies[0][0].cells[0].text == "<b>&é"
    assert page.data_rows()[0].cells[0].text == "<i>"


def test_report_html_gives_a_table_without_rows():
    page = Page(html_of(frame_copy().iloc[:0]))
    assert page.bodies == []
    assert len(page.head) == 4


def frame_without(column):
    return frame_copy().drop(columns=column)


def with_frame(frame):
    return dataclasses.replace(default_report(), report=frame)


# Each case builds its input when the test runs, because a frame needs a call of seed_match_report.
BAD_INPUT = [
    pytest.param(
        lambda: ("x", Guide(GUIDE)),
        "report has type str; expected a SeedMatchReport",
        id="report-type",
    ),
    pytest.param(
        lambda: (default_report(), "UUA"), "guide has type str; expected a Guide", id="guide-type"
    ),
    pytest.param(
        lambda: (with_frame("x"), Guide(GUIDE)),
        "report.report has type str; expected a pandas.DataFrame",
        id="frame-type",
    ),
    pytest.param(
        lambda: (with_frame(frame_without("Total")), Guide(GUIDE)),
        "report table has missing columns: ['Total']; expected "
        f"{['category', *REPORT_CATEGORIES, 'Total', 'group']}",
        id="missing-column",
    ),
    pytest.param(
        lambda: (with_frame(frame_copy().astype({"Total": object})), Guide(GUIDE)),
        "column 'Total' has dtype object; expected a numeric dtype",
        id="text-count",
    ),
    pytest.param(
        lambda: (with_frame(frame_copy().astype({"Full-D3": bool})), Guide(GUIDE)),
        "column 'Full-D3' has dtype bool; expected a numeric dtype",
        id="bool-count",
    ),
    pytest.param(
        lambda: (with_frame(frame_copy().assign(category=[1, 2, 3, 4])), Guide(GUIDE)),
        "report table column 'category' holds 1 in row 0; expected a str in each row",
        id="number-category",
    ),
    pytest.param(
        lambda: (with_frame(frame_copy().assign(group=None)), Guide(GUIDE)),
        "report table column 'group' holds None in row 0; expected a str in each row",
        id="missing-group",
    ),
]


@pytest.mark.parametrize(("build", "message"), BAD_INPUT)
def test_report_html_rejects_bad_input(build, message):
    report, guide = build()
    with pytest.raises(InputError, match=exactly(message)):
        report_html(report, guide)


def test_report_html_runs_on_many_rows_in_bounded_time_and_memory():
    groups = 2_500
    frame = pd.concat([frame_copy()] * groups, ignore_index=True)
    frame["group"] = np.repeat(np.arange(groups), len(ROW_NAMES)).astype(str)
    report = dataclasses.replace(default_report(), report=frame)
    tracemalloc.start()
    try:
        begin = time.perf_counter()
        text = report_html(report, Guide(GUIDE))
        elapsed = time.perf_counter() - begin
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    assert text.count("<tbody>") == groups
    assert text.count("<tr") == 4 + groups * (1 + len(ROW_NAMES))
    assert elapsed < 10
    assert peak < 32 * 2**20
