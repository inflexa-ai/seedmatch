import io
import re
import time
import tracemalloc

import numpy as np
import pandas as pd
import pytest
from helpers import exactly, golden_note, sirna_seqs
from matplotlib.colors import to_hex

from seedmatch import SequenceSet, kmer_counts, kmer_gene_counts, kmer_histogram
from seedmatch.errors import InputError
from seedmatch.plot import plot_kmer_histogram

KMER_HIST = "plots/kmer_hist.tsv"
TITLE = "K-mer gene frequencies"
SIZE_MESSAGE = (
    "size (0, 5) is not valid; expected (width, height) in inches, two finite numbers above 0"
)


def data_gene_counts():
    return kmer_gene_counts(kmer_counts(sirna_seqs()))


def r_pages():
    """Return the (main, xlab, ylab) of each page that R draws, from the note `title calls`."""
    note = golden_note(KMER_HIST, r"^# title calls: (.*)$")
    return [
        re.fullmatch(r"main=(.*), xlab=(.*), ylab=(.*)", call).groups()
        for call in note.split(" / ")
    ]


def draw(target="CTCTATA", **options):
    return plot_kmer_histogram(data_gene_counts(), target, **options)


def test_plot_kmer_histogram_draws_the_bins_and_the_target_D066():
    counts = data_gene_counts()
    histogram = kmer_histogram(counts)
    figure = plot_kmer_histogram(counts, "CTCTATA")
    (axes,) = figure.axes
    bars = axes.patches
    assert len(bars) == len(histogram) == 44
    first = bars[0]
    assert (first.get_x(), first.get_width(), first.get_height()) == (-0.5, 2, 2142)
    np.testing.assert_array_equal([bar.get_x() for bar in bars], histogram["low"] - 0.5)
    np.testing.assert_array_equal(
        [bar.get_width() for bar in bars], histogram["high"] - histogram["low"] + 1
    )
    np.testing.assert_array_equal([bar.get_height() for bar in bars], histogram["kmers"])
    assert {to_hex(bar.get_facecolor()) for bar in bars} == {"#e5e5e5"}
    assert {to_hex(bar.get_edgecolor()) for bar in bars} == {"#000000"}
    assert {bar.get_linewidth() for bar in bars} == {0.75}
    (line,) = axes.get_lines()
    assert list(line.get_xdata()) == [18, 18]
    assert counts["CTCTATA"] == int(golden_note(KMER_HIST, r"^# target CTCTATA: (\d+)$")) == 18
    assert to_hex(line.get_color()) == to_hex("firebrick")
    assert line.get_linewidth() == 2.25
    assert line.get_linewidth() == float(golden_note(KMER_HIST, r"^# abline width: ([\d.]+) w$"))


def test_plot_kmer_histogram_titles_the_first_page_of_r_D066():
    (axes,) = draw().axes
    pages = r_pages()
    assert golden_note(KMER_HIST, r"^# pages: (\d+)$") == "2"
    assert pages[1][0] == "Histogram of gene.counts"
    texts = (axes.get_title(), axes.get_xlabel(), axes.get_ylabel())
    assert texts == pages[0]
    assert texts == (TITLE, "Number of genes with k-mer", "Frequency")
    assert axes.title.get_fontsize() == 14.4
    assert axes.title.get_fontweight() == "bold"


def test_plot_kmer_histogram_has_no_top_and_no_right_edge():
    (axes,) = draw().axes
    visible = {name: spine.get_visible() for name, spine in axes.spines.items()}
    assert visible == {"left": True, "bottom": True, "top": False, "right": False}


def test_plot_kmer_histogram_places_the_target_line_by_its_count():
    counts = pd.Series([7, 12, 12, 30], index=["A", "C", "G", "T"])
    for target, count in zip(counts.index, counts, strict=True):
        (line,) = plot_kmer_histogram(counts, target).axes[0].get_lines()
        assert list(line.get_xdata()) == [count, count]


def test_plot_kmer_histogram_has_the_default_bins_and_size():
    figure = draw()
    assert len(figure.axes[0].patches) == 44
    assert tuple(figure.get_size_inches()) == (7.0, 5.0)
    for rows, bars in ((50, 50), (51, 26)):
        counts = pd.Series(range(rows), index=[f"K{i}" for i in range(rows)])
        assert len(plot_kmer_histogram(counts, "K0").axes[0].patches) == bars


def test_plot_kmer_histogram_takes_the_title_the_bins_and_the_size():
    figure = draw(title="Counts", max_bins=10, size=(4, 3))
    (axes,) = figure.axes
    assert axes.get_title() == "Counts"
    assert len(axes.patches) == 9
    assert axes.patches[0].get_width() == 10
    assert tuple(figure.get_size_inches()) == (4.0, 3.0)


@pytest.mark.parametrize("title", ["$x$", r"cost $5 and $6", r"$\badcommand$"])
def test_plot_kmer_histogram_draws_the_title_as_written(title):
    figure = draw(title=title)
    assert figure.axes[0].title.get_parse_math() is False
    figure.savefig(io.BytesIO(), format="png")


@pytest.mark.parametrize("target", ["CUCUAUA", "ctctata", "CTCTATAA", "NOTAKMER"])
def test_plot_kmer_histogram_rejects_a_target_that_is_not_a_kmer_D066(target):
    assert golden_note(KMER_HIST, rf"^# target {target}: (.*)$") == "NA"
    message = (
        f"target {target!r} is not a k-mer of gene_counts; "
        "expected a k-mer of 7 letters, for example 'AAAAAAA'"
    )
    with pytest.raises(InputError, match=exactly(message)):
        draw(target)


def test_plot_kmer_histogram_names_the_first_label_of_another_index():
    counts = pd.Series([7, 12], index=["ACG", "TTT"])
    message = (
        "target 'AC' is not a k-mer of gene_counts; expected a k-mer of 3 letters, "
        "for example 'ACG'"
    )
    with pytest.raises(InputError, match=exactly(message)):
        plot_kmer_histogram(counts, "AC")


def test_plot_kmer_histogram_rejects_a_target_that_occurs_more_than_one_time():
    counts = pd.Series([7, 12, 9], index=["ACG", "TTT", "ACG"])
    message = "target 'ACG' occurs 2 times in gene_counts; expected one row for each k-mer"
    with pytest.raises(InputError, match=exactly(message)):
        plot_kmer_histogram(counts, "ACG")
    (line,) = plot_kmer_histogram(counts, "TTT").axes[0].get_lines()
    assert list(line.get_xdata()) == [12, 12]


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"title": 5}, "title 5 has type int; expected a str"),
        ({"target": 5}, "target 5 has type int; expected a str"),
        ({"size": (0, 5)}, SIZE_MESSAGE),
        ({"max_bins": 0}, "max_bins 0 is less than 1; expected 1 or more"),
        ({"gene_counts": [1, 2]}, "gene_counts has type list; expected a Series"),
    ],
)
def test_plot_kmer_histogram_rejects_a_bad_argument(kwargs, message):
    arguments = {"gene_counts": data_gene_counts(), "target": "CTCTATA"} | kwargs
    with pytest.raises(InputError, match=exactly(message)):
        plot_kmer_histogram(**arguments)


def test_plot_kmer_histogram_does_not_change_the_counts():
    counts = data_gene_counts()
    before = counts.copy()
    draw()
    pd.testing.assert_series_equal(counts, before)


def test_plot_kmer_histogram_runs_on_the_kmers_of_k_8_in_bounded_time_and_memory():
    labels = kmer_counts(SequenceSet.from_records([]), k=8).kmers()
    generator = np.random.default_rng(1)
    counts = pd.Series(generator.integers(0, 5000, size=len(labels)), index=labels)
    tracemalloc.start()
    try:
        start = time.perf_counter()
        figure = plot_kmer_histogram(counts, labels[-1])
        elapsed = time.perf_counter() - start
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    assert elapsed < 10
    assert peak < 16 * 2**20
    assert len(figure.axes[0].patches) <= 50
