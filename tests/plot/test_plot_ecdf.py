import dataclasses
import io
import math
import re
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
import pytest
from helpers import exactly, golden_note, golden_rows, sirna_de
from matplotlib.colors import to_hex
from matplotlib.text import Text

from seedmatch import EcdfTest, EcdfTestResult, SeedMatchWarning, ecdf_test
from seedmatch.errors import InputError
from seedmatch.plot import plot_ecdf

LINES = "plots/ecdf_lines.tsv"
LABELS = "plots/ecdf_labels.tsv"
RELATIVE = 1e-12


def groups_of(de: pd.DataFrame) -> dict[str, pd.Series]:
    return {
        "Background": de.loc[de["mer7m8"] == 0, "gene_id"],
        "mer7m8": de.loc[de["mer7m8"] >= 1, "gene_id"],
    }


def ks_result(de: pd.DataFrame) -> EcdfTestResult:
    groups = groups_of(de)
    return ecdf_test(de, target=groups["mer7m8"], background=groups["Background"])


def draw(de=None, groups=None, result=None, **options):
    table = sirna_de() if de is None else de
    return plot_ecdf(
        table,
        groups_of(table) if groups is None else groups,
        ks_result(table) if result is None else result,
        **options,
    )


def legend_of(figure) -> tuple[list[str], str]:
    legend = figure.axes[0].get_legend()
    return [text.get_text() for text in legend.get_texts()], legend.get_title().get_text()


def line_of(figure, label: str):
    (line,) = (line for line in figure.axes[0].get_lines() if line.get_label() == label)
    return line


def points_of(figure) -> list[np.ndarray]:
    return [line.get_xydata() for line in figure.axes[0].get_lines()]


def note_legend(key: str) -> list[tuple[str, str]]:
    """Return the label and the color of each legend entry of a note of the golden file."""
    note = golden_note(LINES, rf"^# legend {re.escape(key)}: (.*)$")
    return [
        (match.group(1), match.group(2))
        for entry in note.split(" / ")
        if (match := re.fullmatch(r"(.*) \[(.*)\]", entry))
    ]


def r_layer(case: str) -> dict[str, list[tuple[float, float]]]:
    """Return the layer rows of the R plot for each group, in the order of the file."""
    rows: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for row in golden_rows(LINES):
        if row.case == case:
            rows[row.group].append((float(row.x), float(row.y)))
    return rows


def test_plot_ecdf_draws_the_r_steps():
    de = sirna_de()
    figure = draw(de)
    (axes,) = figure.axes
    background, seed = axes.get_lines()
    assert axes.get_xlim() == pytest.approx((-1.1, 1.1), rel=RELATIVE)
    assert axes.get_ylim() == pytest.approx((-0.05, 1.05), rel=RELATIVE)
    low, high = -8.16461452145323, 3.42272530137105
    assert seed.get_xydata().shape == (20, 2)
    assert seed.get_xydata()[0] == pytest.approx((low, 0.0), rel=RELATIVE)
    assert seed.get_xydata()[1] == pytest.approx(
        (-8.16461452145323, 0.05555555555555555), rel=RELATIVE
    )
    assert seed.get_xydata()[-1] == pytest.approx((high, 1.0), rel=RELATIVE)
    assert background.get_xydata().shape == (163, 2)
    assert background.get_xydata()[1] == pytest.approx(
        (-3.44283360609722, 0.006211180124223602), rel=RELATIVE
    )
    assert [to_hex(line.get_color()) for line in (background, seed)] == ["#000000", "#3498db"]
    assert [line.get_linewidth() for line in (background, seed)] == pytest.approx(
        [3.2009, 3.2009], abs=1e-4
    )
    assert {line.get_drawstyle() for line in (background, seed)} == {"steps-post"}
    assert legend_of(figure) == (["Background: 161", "mer7m8: 18"], "Group")
    assert figure.get_suptitle() == "ECDF"
    assert axes.get_title(loc="left") == "Statistic: 0.51\np-value: p<0.001"
    r_rows = r_layer("test_data")
    for group, line in (("Background", background), ("mer7m8", seed)):
        finite = [(x, y) for x, y in r_rows[group] if math.isfinite(x)]
        got = line.get_xydata()[1:-1]
        assert len(got) == len(finite)
        for (x, y), (expected_x, expected_y) in zip(got, finite, strict=True):
            assert math.isclose(x, expected_x, rel_tol=RELATIVE), group
            assert math.isclose(y, expected_y, rel_tol=RELATIVE), group


def test_plot_ecdf_gives_the_r_titles_and_labels():
    figure = draw()
    (axes,) = figure.axes
    assert figure.get_suptitle() == golden_note(LINES, r"^# labels: title=(.*) / x=")
    assert axes.get_xlabel() == golden_note(LINES, r"^# labels: .* / x=(.*) / y=")
    assert axes.get_ylabel() == golden_note(LINES, r"^# labels: .* / y=(.*)$")
    subtitle = golden_note(LINES, r"^# subtitle KS: (.*)$")
    assert axes.get_title(loc="left").replace("\n", " | ") == subtitle
    low, high = (float(limit) for limit in golden_note(LINES, r"^# coord limits: (.*)$").split(","))
    r_x, r_y = golden_note(LINES, r"^# panel range: (.*)$").split(" / ")
    r_x_range = [float(limit) for limit in r_x.removeprefix("x ").split(",")]
    r_y_range = [float(limit) for limit in r_y.removeprefix("y ").split(",")]
    assert axes.get_xlim() == pytest.approx(r_x_range, rel=RELATIVE)
    assert axes.get_ylim() == pytest.approx(r_y_range, rel=RELATIVE)
    assert (low, high) == (-1.0, 1.0)
    assert r_x_range == pytest.approx([low - 0.1, high + 0.1], rel=RELATIVE)


def test_plot_ecdf_gives_the_r_legend():
    figure = draw()
    r_legend = note_legend("test_data")
    texts, title = legend_of(figure)
    assert texts == [label for label, _ in r_legend]
    assert title == "Group"
    colors = [to_hex(line.get_color()) for line in figure.axes[0].get_lines()]
    assert colors == [to_hex(color) for _, color in r_legend]


def label_cases():
    for row in golden_rows(LABELS):
        if row.p_value == "NA":
            continue
        if row.label.startswith("ERROR"):
            # R stops for these two p-values, and the port gives them a label (D-063).
            label = {"0.01": "0.01", "NaN": "NaN"}[row.p_value]
        else:
            label = row.label
        yield pytest.param(row.p_value, label, id=f"p={row.p_value}")


@pytest.mark.parametrize(("p_value", "label"), list(label_cases()))
def test_plot_ecdf_labels_each_p_value_D063(p_value, label):
    de = sirna_de()
    p = None if p_value == "NULL" else float(p_value)
    result = dataclasses.replace(ks_result(de), p_value=p)
    subtitle = draw(de, result=result).axes[0].get_title(loc="left")
    assert subtitle.split("\n")[1] == f"p-value: {label}"


@pytest.mark.parametrize(
    "statistic", ["0.00050000000000000001", "0.0025000000000000001", "0.1235", "0.0115"]
)
def test_plot_ecdf_rounds_the_statistic_as_r_does(statistic):
    row = next(
        row for row in golden_rows("stats/round.tsv") if (row.digits, row.x) == ("3", statistic)
    )
    # Python rounds the binary value of these numbers, thus its text differs from R.
    assert round(float(statistic), 3) != float(row.r)
    result = dataclasses.replace(ks_result(sirna_de()), statistic=float(statistic))
    assert draw(result=result).axes[0].get_title(loc="left").split("\n")[0] == (
        f"Statistic: {row.r}"
    )


def test_plot_ecdf_labels_the_p_value_of_exactly_one_hundredth_D063():
    row = next(row for row in golden_rows(LABELS) if row.p_value == "0.01")
    assert row.label == "ERROR: object 'pvalue_label' not found"
    result = dataclasses.replace(ks_result(sirna_de()), p_value=20 / 2000)
    assert draw(result=result).axes[0].get_title(loc="left").endswith("p-value: 0.01")


@pytest.mark.parametrize(
    ("test", "note", "statistic", "p_value"),
    [
        (EcdfTest.WILCOXON, "Wilcoxen", "748", "1"),
        (EcdfTest.WASSERSTEIN, "Wass", "2.196", "Undefined"),
        (EcdfTest.DTS, "DTS", "74.912", "p<0.001"),
    ],
)
def test_plot_ecdf_labels_the_statistic_of_each_test(test, note, statistic, p_value):
    de = sirna_de()
    groups = groups_of(de)
    options = {"rng": 0} if test is EcdfTest.DTS else {}
    result = ecdf_test(
        de, target=groups["mer7m8"], background=groups["Background"], test=test, **options
    )
    subtitle = draw(de, result=result).axes[0].get_title(loc="left")
    assert subtitle == f"Statistic: {statistic}\np-value: {p_value}"
    assert subtitle.replace("\n", " | ") == golden_note(LINES, rf"^# subtitle {note}: (.*)$")


def test_plot_ecdf_draws_only_the_groups_D064():
    de = sirna_de()
    extra = de.iloc[[0]].assign(gene_id="ENSNOGROUP", log2FoldChange=5.0)
    with_extra = pd.concat([extra, de], ignore_index=True)
    figure = draw(with_extra, groups_of(de), ks_result(de))
    assert len(figure.axes[0].get_lines()) == 2
    assert legend_of(figure)[0] == ["Background: 161", "mer7m8: 18"]
    assert max(line.get_xydata()[:, 0].max() for line in figure.axes[0].get_lines()) < 5
    assert set(r_layer("no_group_r")) == {"Background", "mer7m8", "NA"}
    assert note_legend("no_group_r")[-1] == ("NA", "grey50")


def test_plot_ecdf_counts_the_values_of_each_line_D064():
    de = sirna_de()
    groups = groups_of(de)
    padded = [*groups["Background"], "ABSENT1", "ABSENT2", groups["Background"].iloc[0]]
    figure = draw(de, {"Background": padded, "mer7m8": groups["mer7m8"]}, ks_result(de))
    assert legend_of(figure)[0] == ["Background: 161", "mer7m8: 18"]
    assert note_legend("absent_and_repeated")[0] == ("Background: 164", "black")
    first_seed = groups["mer7m8"].iloc[0]
    missing = de.copy()
    missing.loc[missing["gene_id"] == first_seed, "log2FoldChange"] = np.nan
    message = "group 'mer7m8' has 1 row(s) without a finite log2FoldChange; the figure omits them"
    with pytest.warns(SeedMatchWarning, match=exactly(message)):
        figure = draw(missing, groups_of(missing), ks_result(de))
    assert legend_of(figure)[0] == ["Background: 161", "mer7m8: 17"]
    assert line_of(figure, "mer7m8: 17").get_xydata().shape == (19, 2)
    assert note_legend("na_fold_change")[1] == ("mer7m8: 18", "#3498db")
    assert golden_note(LINES, r"^# na_fold_change rows of mer7m8: (\d+)$") == "19"


def test_plot_ecdf_follows_the_order_of_the_keys_D064():
    de = sirna_de()
    groups = groups_of(de)
    figure = draw(de, {"mer7m8": groups["mer7m8"], "Background": groups["Background"]})
    assert legend_of(figure)[0] == ["mer7m8: 18", "Background: 161"]
    black, blue = figure.axes[0].get_lines()
    assert to_hex(black.get_color()) == "#000000"
    assert black.get_xydata().shape == (20, 2)
    assert blue.get_xydata().shape == (163, 2)
    assert note_legend("factor_order_reversed")[0] == ("mer7m8: 161", "black")


def test_plot_ecdf_puts_a_comma_in_a_large_count():
    rng = np.random.default_rng(0)
    de = pd.DataFrame(
        {
            "gene_id": [f"G{i}" for i in range(12345)],
            "log2FoldChange": rng.normal(size=12345),
        }
    )
    groups = {"Background": de["gene_id"][:12000], "target": de["gene_id"][12000:]}
    result = ecdf_test(de, target=groups["target"], background=groups["Background"])
    figure = draw(de, groups, result)
    expected = [label for label, _ in note_legend("big_counts")]
    assert expected == ["Background: 12,000", "target: 345"]
    assert legend_of(figure)[0] == expected


def test_plot_ecdf_keeps_the_legend_entry_of_a_name_that_starts_with_an_underscore():
    de = sirna_de()
    groups = groups_of(de)
    figure = draw(de, {"_hidden": groups["Background"], "mer7m8": groups["mer7m8"]}, ks_result(de))
    assert legend_of(figure) == (["_hidden: 161", "mer7m8: 18"], "Group")


@pytest.mark.parametrize(
    ("name", "title"),
    [("$x$", "cost $5 and $6"), (r"$\badcommand$", r"a $\badcommand$ b")],
    ids=["dollar-signs", "bad-formula"],
)
def test_plot_ecdf_draws_a_user_text_as_written(name, title):
    de = sirna_de()
    groups = groups_of(de)
    figure = draw(de, {name: groups["Background"], "mer7m8": groups["mer7m8"]}, title=title)
    assert figure.get_suptitle() == title
    assert legend_of(figure) == ([f"{name}: 161", "mer7m8: 18"], "Group")
    (suptitle,) = (text for text in figure.findobj(Text) if text.get_text() == title)
    user_texts = [suptitle, *figure.axes[0].get_legend().get_texts()]
    assert [text.get_parse_math() for text in user_texts] == [False] * 3
    figure.savefig(io.BytesIO(), format="png")


def test_plot_ecdf_extends_the_steps_to_the_view():
    figure = draw(x_limits=(-10.0, 10.0))
    assert figure.axes[0].get_xlim() == pytest.approx((-11.0, 11.0), rel=RELATIVE)
    for points in points_of(figure):
        assert points[0, 0] == -11.0
        assert points[-1, 0] == 11.0


def test_plot_ecdf_extends_the_steps_to_the_lowest_and_the_highest_value():
    figure = draw(x_limits=(-1.0, 1.0))
    for points in points_of(figure):
        assert points[0, 0] == pytest.approx(-8.16461452145323, rel=RELATIVE)
        assert points[-1, 0] == pytest.approx(3.42272530137105, rel=RELATIVE)


def test_plot_ecdf_steps_up_by_the_count_of_a_repeated_value():
    de = pd.DataFrame({"gene_id": list("abcdef"), "log2FoldChange": [1.0, 1.0, 2.0, 3.0, 3.0, 3.0]})
    groups = {"x": ["a", "b", "c"], "y": ["d", "e", "f"]}
    figure = draw(de, groups, ks_result(sirna_de()), x_limits=(0.0, 4.0))
    first, second = points_of(figure)
    np.testing.assert_allclose(
        first, [[-0.2, 0.0], [1.0, 2 / 3], [2.0, 1.0], [4.2, 1.0]], rtol=RELATIVE
    )
    np.testing.assert_allclose(second, [[-0.2, 0.0], [3.0, 1.0], [4.2, 1.0]], rtol=RELATIVE)


def test_plot_ecdf_drops_a_non_finite_fold_change():
    de = sirna_de()
    groups = groups_of(de)
    changed = de.copy()
    rows = changed["gene_id"].isin(groups["mer7m8"].iloc[:2])
    changed.loc[rows, "log2FoldChange"] = [np.inf, -np.inf]
    message = "group 'mer7m8' has 2 row(s) without a finite log2FoldChange; the figure omits them"
    with pytest.warns(SeedMatchWarning, match=exactly(message)):
        figure = draw(changed, groups_of(changed), ks_result(de))
    assert legend_of(figure)[0] == ["Background: 161", "mer7m8: 16"]
    assert np.isfinite(figure.axes[0].get_lines()[1].get_xydata()).all()


def test_plot_ecdf_gives_no_warning_for_finite_values():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        draw()


def test_plot_ecdf_warns_for_a_gene_in_two_groups():
    de = sirna_de()
    groups = groups_of(de)
    shared = {
        "Background": pd.concat([groups["Background"], groups["mer7m8"].iloc[:1]]),
        "mer7m8": groups["mer7m8"],
    }
    with pytest.warns(SeedMatchWarning, match="share the gene"):
        draw(de, shared, ks_result(de))


BAD_ARGUMENTS = [
    pytest.param(
        {"palette": ("black",)},
        "palette has 1 color(s) for 2 groups; expected a color for each group",
        id="short-palette",
    ),
    pytest.param(
        {"palette": ("black", "blu")},
        "palette[1] 'blu' is not a matplotlib color; expected a color name or a hex code, "
        "for example '#3498db'",
        id="not-a-color",
    ),
    pytest.param(
        {"palette": "black"},
        "palette 'black' is a str; expected a sequence of str, for example ('black',)",
        id="palette-text",
    ),
    pytest.param(
        {"palette": ("black", 5)},
        "palette[1] 5 has type int; expected a str",
        id="palette-number",
    ),
    pytest.param(
        {"x_limits": (1.0, -1.0)},
        "x_limits (1.0, -1.0) is not valid; expected (low, high), "
        "two finite numbers from -9e306 to 9e306 with low < high and high - low at most 9e306",
        id="reversed-limits",
    ),
    pytest.param(
        {"x_limits": (0.0, math.inf)},
        "x_limits (0.0, inf) is not valid; expected (low, high), "
        "two finite numbers from -9e306 to 9e306 with low < high and high - low at most 9e306",
        id="infinite-limit",
    ),
    pytest.param(
        {"size": (0, 5)},
        "size (0, 5) is not valid; expected (width, height) in inches, "
        "two finite numbers above 0 and below 50",
        id="zero-size",
    ),
    pytest.param({"title": 5}, "title 5 has type int; expected a str", id="title"),
    pytest.param({"result": "x"}, "result has type str; expected a EcdfTestResult", id="result"),
    pytest.param({"groups": {}}, "groups has no group; expected at least one group", id="no-group"),
    pytest.param(
        {"groups": ["a"]},
        "groups ['a'] has type list; expected a mapping from a group name to gene IDs",
        id="groups-list",
    ),
]


@pytest.mark.parametrize(("kwargs", "message"), BAD_ARGUMENTS)
def test_plot_ecdf_rejects_a_bad_argument(kwargs, message):
    with pytest.raises(InputError, match=exactly(message)):
        draw(**kwargs)


WIDE_GROUPS = {"x": ["a", "b", "c"], "y": ["d", "e", "f"]}
RANGE = "expected values from -9e306 to 9e306 with a range of at most 9e306"


def wide_de(low: float, high: float) -> pd.DataFrame:
    """Return a table whose groups have the lowest value `low` and the highest value `high`."""
    values = [low, 0.0, 1.0, 2.0, 3.0, high]
    return pd.DataFrame({"gene_id": list("abcdef"), "log2FoldChange": values})


@pytest.mark.parametrize(
    "limits", [(-1e308, 1e308), (-4.6e306, 4.6e306), (1e308, 1.0001e308), (8.9e306, 9.1e306)]
)
def test_plot_ecdf_rejects_x_limits_wider_than_the_bound(limits):
    message = (
        f"x_limits {limits!r} is not valid; expected (low, high), "
        "two finite numbers from -9e306 to 9e306 with low < high and high - low at most 9e306"
    )
    with pytest.raises(InputError, match=exactly(message)):
        draw(x_limits=limits)


def test_plot_ecdf_rejects_values_wider_than_the_bound():
    message = f"the values of log2FoldChange go from -1e+308 to 1e+308; {RANGE}"
    with pytest.raises(InputError, match=exactly(message)):
        draw(wide_de(-1e308, 1e308), WIDE_GROUPS, ks_result(sirna_de()))


def test_plot_ecdf_rejects_values_outside_the_bound():
    de = pd.DataFrame({"gene_id": list("abcdef"), "log2FoldChange": [1.75e308] * 6})
    message = f"the values of log2FoldChange go from 1.75e+308 to 1.75e+308; {RANGE}"
    with pytest.raises(InputError, match=exactly(message)):
        draw(de, WIDE_GROUPS, ks_result(sirna_de()))


@pytest.mark.parametrize("size", [(5.0, 5.0), (2.0, 2.0)], ids=["default-size", "small-size"])
def test_plot_ecdf_draws_at_the_bound(size):
    de = wide_de(-9e306, 0.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        figure = draw(de, WIDE_GROUPS, ks_result(sirna_de()), x_limits=(0.0, 9e306), size=size)
        figure.savefig(io.BytesIO(), format="png")


@pytest.mark.parametrize("missing", ["absent", "nan"])
def test_plot_ecdf_rejects_a_group_without_a_value(missing):
    de = sirna_de()
    groups = groups_of(de)
    genes = ["ENSABSENT"] if missing == "absent" else list(groups["mer7m8"])
    if missing == "nan":
        de = de.copy()
        de.loc[de["gene_id"].isin(genes), "log2FoldChange"] = np.nan
    message = (
        "group 'mer7m8' has no row with a finite log2FoldChange in the DE table; "
        "expected at least 1 row"
    )
    with pytest.raises(InputError, match=exactly(message)):
        draw(de, {"Background": groups["Background"], "mer7m8": genes}, ks_result(sirna_de()))


def test_plot_ecdf_rejects_a_de_without_log2_fold_change():
    de = sirna_de().drop(columns="log2FoldChange")
    message = (
        "DE table has missing columns: ['log2FoldChange']; expected ['gene_id', 'log2FoldChange']"
    )
    with pytest.raises(InputError, match=exactly(message)):
        plot_ecdf(de, {"a": ["x"]}, ks_result(sirna_de()))


def test_plot_ecdf_rejects_a_text_log2_fold_change():
    de = sirna_de().astype({"log2FoldChange": str})
    message = "column 'log2FoldChange' has dtype object; expected a numeric dtype"
    with pytest.raises(InputError, match=exactly(message)):
        plot_ecdf(de, groups_of(sirna_de()), ks_result(sirna_de()))


def test_plot_ecdf_does_not_change_de():
    de = sirna_de()
    before = de.copy()
    draw(de)
    pd.testing.assert_frame_equal(de, before)


def test_plot_ecdf_takes_the_size_and_the_title_and_the_palette():
    figure = draw(title="My title", size=(4.0, 3.0), palette=("red", "blue"))
    assert figure.get_suptitle() == "My title"
    assert tuple(figure.get_size_inches()) == (4.0, 3.0)
    colors = [to_hex(line.get_color()) for line in figure.axes[0].get_lines()]
    assert colors == ["#ff0000", "#0000ff"]


def test_plot_ecdf_styles_the_axes_with_the_theme_of_r():
    (axes,) = draw().axes
    grid = golden_note(LINES, r"panel\.grid\.major colour (#[0-9a-f]{6})")
    panel = golden_note(LINES, r"panel\.background fill white colour (#[0-9a-f]{6})")
    assert to_hex(axes.xaxis.get_gridlines()[0].get_color()) == grid
    assert {to_hex(spine.get_edgecolor()) for spine in axes.spines.values()} == {panel}
    legend = axes.get_legend()
    # Legend has no public getter of its location code.
    assert legend._loc == 2  # noqa: SLF001
    assert to_hex(legend.get_frame().get_edgecolor()) == panel
    assert to_hex(legend.get_frame().get_facecolor()) == "#ffffff"
    assert legend.get_title().get_fontweight() == "bold"
    assert legend.get_title().get_fontsize() == 12
    (subtitle,) = (
        child
        for child in axes.get_children()
        if isinstance(child, Text) and child.get_text().startswith("Statistic")
    )
    assert subtitle.get_fontsize() == 12
