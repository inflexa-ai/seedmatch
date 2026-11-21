import re
from collections import Counter

import numpy as np
import pandas as pd
import pytest
from helpers import DATA, GUIDE, exactly, golden, golden_note, sirna_seqs
from hypothesis import given, settings
from hypothesis import strategies as st

from seedmatch import Guide, Seed, SequenceSet, kmer_counts, kmer_gene_counts, kmer_histogram
from seedmatch.errors import InputError
from seedmatch.kmers import MAX_K

KMER_COUNTS = golden("analyses/kmer_counts.tsv")
KMER_GENE_COUNTS = golden("analyses/kmer_gene_counts.tsv")
KMER_HIST = "plots/kmer_hist.tsv"
CASES = list(dict.fromkeys(KMER_COUNTS["case"]))
CASE_PARAMETERS = re.compile(r"_k(\d+)_s(\d+)$")


def case_set(case):
    if case.startswith("test_data"):
        everything = sirna_seqs()
        size = 10 if "first10" in case else len(everything)
        return SequenceSet.from_records(
            (everything.ids[i], everything.sequence(i)) for i in range(size)
        )
    fixture = pd.read_csv(DATA / "kmer_cases.tsv", sep="\t", dtype=str, keep_default_na=False)
    return SequenceSet.from_records(
        zip(fixture["seq_id"], fixture["sequence"], strict=True), gene_ids=fixture["gene_id"]
    )


def case_counts(case):
    k, step = (int(value) for value in CASE_PARAMETERS.search(case).groups())
    return kmer_counts(case_set(case), k=k, step=step)


def cells(counts):
    """Return the entries of the matrix as (row, k-mer): count."""
    labels = counts.kmers()
    coo = counts.matrix.tocoo()
    return {
        (int(row), labels[column]): int(value)
        for row, column, value in zip(coo.row, coo.col, coo.data, strict=True)
    }


def listed_cells(counts):
    """Return the entries as (seq_id, gene_id, k-mer, count), by row and then by column."""
    matrix = counts.matrix
    rows = np.repeat(np.arange(matrix.shape[0]), np.diff(matrix.indptr))
    labels = counts.kmers()
    return [
        (counts.ids[row], counts.gene_ids[row], labels[column], int(value))
        for row, column, value in zip(rows, matrix.indices, matrix.data, strict=True)
    ]


@pytest.mark.parametrize("case", CASES)
def test_kmer_counts_match_r(case):
    rows = KMER_COUNTS[KMER_COUNTS["case"] == case]
    expected = list(
        zip(
            rows["seq_id"],
            rows["gene_id"],
            rows["kmer"],
            rows["count"].astype(int),
            strict=True,
        )
    )
    assert listed_cells(case_counts(case)) == expected


def test_to_frame_has_the_r_layout():
    frame = case_counts("cases_k2_s1").to_frame()
    expected = golden("analyses/kmer_frame.tsv")
    assert list(frame.columns) == ["seq_id", *expected.columns[1:]]
    assert frame["seq_id"].tolist() == expected["tx_id"].tolist()
    assert frame["gene_id"].tolist() == expected["gene_id"].tolist()
    assert (frame.dtypes.iloc[2:] == "int64").all()
    counts = expected.iloc[:, 2:].astype(int)
    assert frame.iloc[:, 2:].to_numpy().tolist() == counts.to_numpy().tolist()
    assert frame["GG"].eq(0).all()
    assert isinstance(frame.index, pd.RangeIndex)


def test_the_columns_are_in_the_order_of_biostrings():
    labels = kmer_counts(SequenceSet.from_records([]), k=2).kmers()
    assert labels[:6] == ("AA", "AC", "AG", "AT", "CA", "CC")
    assert (len(labels), labels[-1]) == (16, "TT")


CHUNK = st.one_of(
    st.text(alphabet="ACGT", max_size=10),
    st.text(alphabet="acgt", max_size=6),
    st.text(alphabet="MRWSYKVHDBN", min_size=1, max_size=3),
    st.integers(1, 6).map(lambda size: "N" * size),
)
SEQUENCES = st.lists(st.lists(CHUNK, max_size=5).map("".join), max_size=6)


def oracle(texts, k, step):
    expected = Counter()
    for row, text in enumerate(texts):
        upper = text.upper()
        for start in range(0, len(upper) - k + 1, step):
            window = upper[start : start + k]
            if set(window) <= set("ACGT"):
                expected[row, window] += 1
    return dict(expected)


@settings(deadline=None)
@given(texts=SEQUENCES, k=st.integers(1, 4), step=st.integers(1, 4))
def test_kmer_counts_equal_an_oracle(texts, k, step):
    sequences = SequenceSet.from_records(
        [(f"s{i}", text) for i, text in enumerate(texts)],
        gene_ids=[f"g{i % 3}" for i in range(len(texts))],
    )
    # A block of 7 starts splits most sequences, and a window can cross the end of a block.
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr("seedmatch.kmers._BLOCK_STARTS", 7)
        counts = kmer_counts(sequences, k=k, step=step)
    matrix = counts.matrix
    assert cells(counts) == oracle(texts, k, step)
    assert matrix.shape == (len(texts), 4**k)
    assert matrix.dtype == np.int64
    assert matrix.has_canonical_format
    assert (matrix.data > 0).all()
    assert (counts.k, counts.step) == (k, step)
    assert (counts.ids, counts.gene_ids) == (sequences.ids, sequences.gene_ids)


def test_kmer_counts_skip_a_window_with_an_ambiguous_letter():
    sequences = SequenceSet.from_records([("t", "ACGTNACGTA")])
    assert cells(kmer_counts(sequences, k=2)) == {
        (0, "AC"): 2,
        (0, "CG"): 2,
        (0, "GT"): 2,
        (0, "TA"): 1,
    }
    assert cells(kmer_counts(sequences, k=2, step=2)) == {
        (0, "AC"): 1,
        (0, "GT"): 1,
        (0, "CG"): 1,
        (0, "TA"): 1,
    }


def test_kmer_counts_accept_k_12():
    counts = kmer_counts(SequenceSet.from_records([("t", "ACGTACGTACGT")]), k=MAX_K)
    assert counts.matrix.shape == (1, 4**12)
    assert counts.matrix.indices.tolist() == [int("012301230123", 4)]


def test_kmer_gene_counts_match_r_for_one_sequence_per_gene():
    result = kmer_gene_counts(kmer_counts(sirna_seqs(), k=7))
    rows = KMER_GENE_COUNTS[KMER_GENE_COUNTS["case"] == "test_data_k7_s1"]
    assert len(rows) == 4**7
    assert result.index.tolist() == rows["kmer"].tolist()
    assert result.tolist() == rows["n"].astype(int).tolist()
    assert (result.name, result.index.name, result.dtype) == ("genes", "kmer", "int64")


def test_kmer_gene_counts_give_the_counts_of_kmer_hist_D010():
    note = golden_note(
        "analyses/kmer_gene_counts.tsv",
        r"n of the mer7m8 target CTCTATA for test_data_k7_s1 is (\d+)",
    )
    target = Guide(GUIDE).target(Seed.MER7M8)
    result = kmer_gene_counts(kmer_counts(sirna_seqs(), k=7))
    assert target.sequence == "CTCTATA"
    assert result[target.sequence] == int(note)


def test_kmer_gene_counts_count_genes_not_sequences_D042():
    result = kmer_gene_counts(case_counts("cases_k2_s1"))
    rows = KMER_COUNTS[KMER_COUNTS["case"] == "cases_k2_s1"]
    genes = rows.groupby("kmer")["gene_id"].nunique().reindex(result.index, fill_value=0)
    in_r = KMER_GENE_COUNTS[KMER_GENE_COUNTS["case"] == "cases_k2_s1"]
    r_value = in_r.set_index("kmer")["n"].astype(int).reindex(result.index)
    assert result.tolist() == genes.tolist()
    assert (r_value >= result).all()
    assert (r_value > result).any()


def test_kmer_gene_counts_of_an_empty_set_are_0_for_each_kmer():
    result = kmer_gene_counts(kmer_counts(SequenceSet.from_records([]), k=3))
    assert (len(result), result.sum(), result.dtype) == (64, 0, "int64")


def test_kmer_counts_reject_k_above_12_D043():
    sequences = SequenceSet.from_records([("t", "ACGT")])
    message = "k 13 is out of range; expected an integer from 1 to 12"
    with pytest.raises(InputError, match=exactly(message)):
        kmer_counts(sequences, k=13)


def test_kmer_counts_reject_k_below_1():
    sequences = SequenceSet.from_records([("t", "ACGT")])
    message = "k 0 is out of range; expected an integer from 1 to 12"
    with pytest.raises(InputError, match=exactly(message)):
        kmer_counts(sequences, k=0)


def test_kmer_counts_reject_a_bool_k():
    sequences = SequenceSet.from_records([("t", "ACGT")])
    message = "k True has type bool; expected an integer"
    with pytest.raises(InputError, match=exactly(message)):
        kmer_counts(sequences, k=True)


def test_kmer_counts_reject_step_below_1():
    sequences = SequenceSet.from_records([("t", "ACGT")])
    message = "step 0 is less than 1; expected 1 or more"
    with pytest.raises(InputError, match=exactly(message)):
        kmer_counts(sequences, step=0)


def test_kmer_counts_reject_a_list_of_sequences():
    message = "sequences has type list; expected a SequenceSet"
    with pytest.raises(InputError, match=exactly(message)):
        kmer_counts(["ACGT"])


def test_kmer_gene_counts_reject_a_data_frame():
    message = "counts has type DataFrame; expected a KmerCounts"
    with pytest.raises(InputError, match=exactly(message)):
        kmer_gene_counts(pd.DataFrame())


def test_kmer_counts_keep_each_sequence_with_a_repeated_id_D044():
    counts = kmer_counts(SequenceSet.from_records([("t1", "ACGT"), ("t1", "AAAA")]), k=2)
    assert counts.ids == ("t1", "t1")
    assert cells(counts) == {(0, "AC"): 1, (0, "CG"): 1, (0, "GT"): 1, (1, "AA"): 3}
    assert len(counts.to_frame()) == 2


def test_kmer_counts_of_an_empty_set():
    counts = kmer_counts(SequenceSet.from_records([]), k=3)
    frame = counts.to_frame()
    assert counts.matrix.shape == (0, 64)
    assert list(frame.columns) == ["seq_id", "gene_id", *counts.kmers()]
    assert len(frame) == 0


def data_gene_counts():
    return kmer_gene_counts(kmer_counts(sirna_seqs()))


def r_gene_count_table():
    """Return the number of k-mers for each gene count, as `table(gene.counts)` of R gives it."""
    table = golden(KMER_HIST)
    return pd.Series(table["kmers"].astype(int).to_numpy(), index=table["gene_count"].astype(int))


def test_kmer_histogram_bins_the_integer_counts_D066():
    histogram = kmer_histogram(data_gene_counts())
    assert list(histogram.columns) == ["low", "high", "kmers"]
    assert (histogram.dtypes == "int64").all()
    assert isinstance(histogram.index, pd.RangeIndex)
    assert histogram.iloc[:2].to_numpy().tolist() == [[0, 1, 2142], [2, 3, 2242]]
    assert histogram.iloc[-1].tolist() == [86, 87, 1]
    assert len(histogram) == 44
    assert histogram["kmers"].sum() == 16384
    # The first bin of R holds the counts 0, 1, and 2, and the first bin of the port holds 0 and 1.
    r_table = r_gene_count_table()
    first_bar_of_r = int(golden_note(KMER_HIST, r"^# hist counts: (\d+),"))
    assert first_bar_of_r == r_table.loc[:2].sum() == 3394
    assert histogram["kmers"].iloc[0] == r_table.loc[:1].sum() == 2142


def test_kmer_histogram_defaults_to_at_most_50_bins():
    assert len(kmer_histogram(data_gene_counts())) == 44
    fifty = kmer_histogram(pd.Series(range(50)))
    assert len(fifty) == 50
    assert (fifty["high"] == fifty["low"]).all()
    wider = kmer_histogram(pd.Series(range(51)))
    assert len(wider) == 26
    assert (wider["high"] - wider["low"] == 1).all()


def test_kmer_histogram_sums_the_counts_of_r():
    r_table = r_gene_count_table()
    histogram = kmer_histogram(data_gene_counts())
    assert r_table.sum() == 16384
    for row in histogram.itertuples():
        assert r_table.loc[row.low : row.high].sum() == row.kmers


def test_kmer_histogram_widens_the_bins_to_max_bins():
    counts = data_gene_counts()
    ten = kmer_histogram(counts, max_bins=10)
    assert len(ten) == 9
    assert (ten["high"] - ten["low"] + 1 == 10).all()
    assert ten["kmers"].tolist()[:5] == [8484, 6021, 1570, 225, 59]
    assert ten["kmers"].tolist()[-3:] == [4, 1, 1]
    assert kmer_histogram(counts, max_bins=3).to_numpy().tolist() == [
        [0, 49, 16359],
        [50, 99, 25],
    ]


def test_kmer_histogram_starts_at_the_bin_of_the_lowest_count():
    counts = pd.Series([7, 12, 12, 30], index=["A", "C", "G", "T"])
    fifty = kmer_histogram(counts)
    assert len(fifty) == 24
    assert fifty.iloc[0].tolist() == [7, 7, 1]
    assert fifty.iloc[-1].tolist() == [30, 30, 1]
    assert fifty.loc[fifty["low"] == 12, "kmers"].tolist() == [2]
    assert kmer_histogram(counts, max_bins=5).to_numpy().tolist() == [
        [0, 9, 1],
        [10, 19, 2],
        [20, 29, 0],
        [30, 39, 1],
    ]


def test_kmer_histogram_gives_one_bin_for_equal_counts():
    assert kmer_histogram(pd.Series([5, 5, 5])).to_numpy().tolist() == [[5, 5, 3]]
    assert kmer_histogram(pd.Series([0, 3]), max_bins=1).to_numpy().tolist() == [[0, 4, 2]]


@pytest.mark.parametrize("dtype", ["uint8", "int32", "Int64", "UInt16"])
def test_kmer_histogram_accepts_each_integer_dtype(dtype):
    counts = pd.Series([1, 2, 2, 9], index=list("ACGT"), dtype=dtype)
    histogram = kmer_histogram(counts, max_bins=3)
    assert histogram.to_numpy().tolist() == [[0, 4, 3], [5, 9, 1]]
    assert (histogram.dtypes == "int64").all()


def test_kmer_histogram_computes_the_bins_in_int64_for_a_narrow_dtype():
    counts = pd.Series([0, 200], dtype="uint8")
    assert kmer_histogram(counts, max_bins=1).to_numpy().tolist() == [[0, 499, 2]]


def test_kmer_histogram_counts_the_bins_from_the_lowest_count():
    counts = pd.Series([2**60, 2**60 + 3])
    assert kmer_histogram(counts).to_numpy().tolist() == [
        [2**60, 2**60, 1],
        [2**60 + 1, 2**60 + 1, 0],
        [2**60 + 2, 2**60 + 2, 0],
        [2**60 + 3, 2**60 + 3, 1],
    ]


def test_kmer_histogram_accepts_the_largest_count():
    counts = pd.Series([0, 2**61], index=["A", "C"])
    histogram = kmer_histogram(counts, max_bins=2)
    assert histogram["kmers"].sum() == 2
    assert histogram["high"].iloc[-1] >= 2**61
    assert (histogram.dtypes == "int64").all()
    assert kmer_histogram(counts, max_bins=1).to_numpy().tolist() == [[0, 5 * 10**18 - 1, 2]]


@settings(deadline=None)
@given(
    values=st.lists(st.integers(0, 3000), min_size=1, max_size=60),
    max_bins=st.integers(1, 80),
)
def test_kmer_histogram_equals_an_oracle(values, max_bins):
    counts = pd.Series(values, index=[f"k{i}" for i in range(len(values))])
    histogram = kmer_histogram(counts, max_bins=max_bins)
    lo, hi = min(values), max(values)
    widths = {factor * 10**power for power in range(5) for factor in (1, 2, 5)}
    fits = [w for w in sorted(widths) if hi // w - lo // w + 1 <= max_bins]
    width = fits[0]
    assert (histogram["high"] - histogram["low"] + 1 == width).all()
    assert histogram["low"].tolist() == [b * width for b in range(lo // width, hi // width + 1)]
    expected = Counter(value // width for value in values)
    assert histogram["kmers"].tolist() == [expected[b] for b in range(lo // width, hi // width + 1)]
    assert len(histogram) <= max_bins


@pytest.mark.parametrize(
    ("counts", "kwargs", "message"),
    [
        ([1, 2], {}, "gene_counts has type list; expected a Series"),
        (
            pd.Series([1.0, 2.0], index=["A", "C"]),
            {},
            "gene_counts has dtype float64; expected an integer dtype, as kmer_gene_counts gives",
        ),
        (
            pd.Series([True, False], index=["A", "C"]),
            {},
            "gene_counts has dtype bool; expected an integer dtype, as kmer_gene_counts gives",
        ),
        (
            pd.Series([1, None, 3, None], index=["A", "C", "G", "T"], dtype="Int64"),
            {},
            "gene_counts has a missing value for k-mer 'C'; expected a count",
        ),
        (
            pd.Series([], dtype="int64"),
            {},
            "gene_counts has no value; expected at least one k-mer",
        ),
        (
            pd.Series([-1, 5], index=["AAA", "AAC"]),
            {},
            "gene_counts holds -1 for k-mer 'AAA'; expected 0 or more",
        ),
        (
            pd.Series([4, 0, -3, -1], index=["AAA", "AAC", "AAG", "AAT"]),
            {},
            "gene_counts holds -3 for k-mer 'AAG'; expected 0 or more",
        ),
        (
            pd.Series([0, 2**61 + 1], index=["AAA", "AAC"]),
            {},
            "gene_counts holds 2305843009213693953 for k-mer 'AAC'; expected at most 2^61",
        ),
        (pd.Series([1, 2]), {"max_bins": 0}, "max_bins 0 is less than 1; expected 1 or more"),
        (pd.Series([1, 2]), {"max_bins": True}, "max_bins True has type bool; expected an integer"),
    ],
    ids=[
        "list",
        "float64",
        "bool",
        "missing",
        "empty",
        "negative",
        "later negative",
        "too large",
        "max_bins 0",
        "max_bins bool",
    ],
)
def test_kmer_histogram_rejects_bad_input(counts, kwargs, message):
    with pytest.raises(InputError, match=exactly(message)):
        kmer_histogram(counts, **kwargs)


def test_kmer_histogram_does_not_change_the_counts():
    counts = data_gene_counts()
    before = counts.copy()
    kmer_histogram(counts)
    pd.testing.assert_series_equal(counts, before)
