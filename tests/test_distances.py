import itertools
from typing import assert_never

import numpy as np
import pandas as pd
import pytest
from helpers import GUIDE, NOT_A_DATA_FRAME, TEXT_DTYPE, exactly, golden, sirna_seqs
from hypothesis import given, settings
from hypothesis import strategies as st

from seedmatch import Guide, MatchAnchor, Seed, find_matches, match_distances
from seedmatch.errors import InputError

RANGES = golden("analyses/distances_ranges.tsv")
DISTANCES = golden("analyses/distances.tsv")
CASES = list(dict.fromkeys(RANGES["case"]))
COLUMNS = ["seq_id", "gene_id", "distance"]


def case_matches(case):
    rows = RANGES[RANGES["case"] == case]
    return pd.DataFrame(
        {
            "seq_id": rows["seq_id"].to_numpy(),
            "gene_id": ("gene of " + rows["seq_id"]).to_numpy(),
            "start": rows["start"].astype(np.int64).to_numpy(),
            "end": rows["end"].astype(np.int64).to_numpy(),
            "width": rows["width"].astype(np.int64).to_numpy(),
        }
    )


def golden_distances(case, anchor):
    return DISTANCES[(DISTANCES["case"] == case) & (DISTANCES["anchor"] == anchor.value)]


def matches_of(rows):
    """Return a matches table from (seq_id, start) pairs. The width of a match is 1 to 7."""
    return pd.DataFrame(
        {
            "seq_id": [seq_id for seq_id, _ in rows],
            "gene_id": [f"gene of {seq_id}" for seq_id, _ in rows],
            "start": np.array([start for _, start in rows], dtype=np.int64),
            "end": np.array([start + start % 7 for _, start in rows], dtype=np.int64),
            "width": np.array([1 + start % 7 for _, start in rows], dtype=np.int64),
        }
    )


def coordinate(start, anchor):
    width = 1 + start % 7
    match anchor:
        case MatchAnchor.FIVE_PRIME_END:
            return start + width - 1
        case MatchAnchor.THREE_PRIME_END:
            return start
        case MatchAnchor.CENTER:
            return start + width / 2
        case _:
            assert_never(anchor)


@pytest.mark.parametrize("anchor", list(MatchAnchor))
@pytest.mark.parametrize("case", CASES)
def test_match_distances_match_r(case, anchor):
    matches = case_matches(case)
    sizes = matches["seq_id"].value_counts()
    expected = golden_distances(case, anchor)
    expected = expected[expected["seq_id"].map(sizes) >= 2]
    got = match_distances(matches, anchor=anchor)
    assert "NA" not in expected["distance"].tolist()
    assert got["seq_id"].tolist() == expected["seq_id"].tolist()
    assert got["distance"].tolist() == pytest.approx(expected["distance"].astype(float), rel=1e-9)


def test_a_sequence_with_one_match_gives_no_distance_D011():
    matches = case_matches("mer6")
    sizes = matches["seq_id"].value_counts()
    single = list(dict.fromkeys(matches["seq_id"][matches["seq_id"].map(sizes) == 1]))
    assert single
    rows = golden_distances("mer6", MatchAnchor.FIVE_PRIME_END)
    for seq_id in single:
        assert rows.loc[rows["seq_id"] == seq_id, "distance"].tolist() == ["NA", "0"]
    got = match_distances(matches)
    assert not set(single) & set(got["seq_id"])
    assert set(got["seq_id"]) == set(sizes[sizes >= 2].index)


def test_a_table_of_one_match_gives_no_row_D011():
    got = match_distances(case_matches("synthetic").query("seq_id == 'b'"))
    assert list(got.columns) == COLUMNS
    assert len(got) == 0


def test_center_is_the_start_plus_half_the_width():
    matches = case_matches("synthetic")
    got = match_distances(matches, anchor=MatchAnchor.CENTER)
    assert got.loc[got["seq_id"] == "a", "distance"].tolist() == [20, 25.5, 5.5]


STARTS = st.lists(st.tuples(st.sampled_from("abcd"), st.integers(-20, 200)), max_size=24)


@settings(max_examples=200, deadline=None)
@given(starts=STARTS, anchor=st.sampled_from(list(MatchAnchor)))
def test_each_pair_occurs_once(starts, anchor):
    groups = {}
    for seq_id, start in starts:
        groups.setdefault(seq_id, []).append(coordinate(start, anchor))
    expected = [
        (seq_id, f"gene of {seq_id}", abs(first - second))
        for seq_id, coordinates in groups.items()
        for first, second in itertools.combinations(coordinates, 2)
    ]
    got = match_distances(matches_of(starts), anchor=anchor)
    assert list(got.itertuples(index=False, name=None)) == expected


def test_the_rows_of_one_sequence_need_not_be_adjacent():
    matches = matches_of([("a", 10), ("b", 50), ("a", 30), ("c", 5), ("a", 31)])
    got = match_distances(matches, anchor=MatchAnchor.THREE_PRIME_END)
    assert got["seq_id"].tolist() == ["a", "a", "a"]
    assert got["distance"].tolist() == [20, 21, 1]


def test_all_rows_with_a_missing_seq_id_are_one_sequence():
    matches = matches_of([("a", 10), (None, 20), ("a", 30), (None, 25)])
    got = match_distances(matches, anchor=MatchAnchor.THREE_PRIME_END)
    assert got["seq_id"].iloc[0] == "a"
    assert pd.isna(got["seq_id"].iloc[1])
    assert got["distance"].tolist() == [20, 5]


def test_gene_id_is_the_gene_of_the_first_row_of_the_sequence():
    matches = matches_of([("a", 10), ("b", 20), ("a", 30)])
    matches["gene_id"] = ["first", "other", "second"]
    assert match_distances(matches)["gene_id"].tolist() == ["first"]


def test_rows_keep_the_order_of_the_sequences():
    sequences = sirna_seqs()
    matches = find_matches(Guide(GUIDE).target(Seed.MER6), sequences)
    got = match_distances(matches)
    position = {seq_id: index for index, seq_id in enumerate(sequences.ids)}
    seen = [position[seq_id] for seq_id in got["seq_id"].drop_duplicates()]
    assert len(seen) > 1
    assert seen == sorted(seen)
    assert list(got.columns) == COLUMNS
    assert got["distance"].dtype == np.float64
    assert isinstance(got.index, pd.RangeIndex)
    assert (got["gene_id"] == got["seq_id"]).all()


@pytest.mark.parametrize(("matches", "type_name"), NOT_A_DATA_FRAME)
def test_match_distances_rejects_a_table_that_is_not_a_data_frame(matches, type_name):
    message = f"matches has type {type_name}; expected a pandas.DataFrame"
    with pytest.raises(InputError, match=exactly(message)):
        match_distances(matches)


def test_match_distances_rejects_a_missing_column():
    matches = case_matches("synthetic").drop(columns="width")
    message = (
        "matches table has missing columns: ['width']; "
        "expected ['seq_id', 'gene_id', 'start', 'end', 'width']"
    )
    with pytest.raises(InputError, match=exactly(message)):
        match_distances(matches)


def test_match_distances_rejects_a_float_start():
    matches = case_matches("synthetic")
    matches["start"] = matches["start"].astype(float)
    message = "column 'start' has dtype float64; expected an integer dtype"
    with pytest.raises(InputError, match=exactly(message)):
        match_distances(matches)


def test_match_distances_rejects_a_bool_end():
    matches = case_matches("synthetic")
    matches["end"] = True
    message = "column 'end' has dtype bool; expected an integer dtype"
    with pytest.raises(InputError, match=exactly(message)):
        match_distances(matches)


def test_match_distances_rejects_a_missing_width():
    matches = case_matches("synthetic")
    matches["width"] = pd.array([7, None, 8, 7, 7, 7], dtype="Int64")
    message = "column 'width' has missing values; expected an integer in each row"
    with pytest.raises(InputError, match=exactly(message)):
        match_distances(matches)


def test_match_distances_rejects_an_anchor_as_text():
    message = (
        "anchor 'center' is not valid; expected one of "
        "MatchAnchor.FIVE_PRIME_END, MatchAnchor.THREE_PRIME_END, MatchAnchor.CENTER"
    )
    with pytest.raises(InputError, match=exactly(message)):
        match_distances(case_matches("synthetic"), anchor="center")


def test_ids_that_are_not_text_keep_their_values():
    matches = pd.DataFrame(
        {"seq_id": [7, 7], "gene_id": [70, 70], "start": [1, 9], "end": [5, 13], "width": [5, 5]}
    )
    got = match_distances(matches)
    assert got["seq_id"].tolist() == [7]
    assert got["gene_id"].tolist() == [70]
    assert type(got["seq_id"][0]) is int


def test_an_empty_table_gives_an_empty_result():
    got = match_distances(case_matches("synthetic").iloc[:0])
    assert list(got.columns) == COLUMNS
    assert len(got) == 0
    assert [str(got[name].dtype) for name in ("seq_id", "gene_id")] == [TEXT_DTYPE] * 2
    assert got["distance"].dtype == np.float64
    assert isinstance(got.index, pd.RangeIndex)
