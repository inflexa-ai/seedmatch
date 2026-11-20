import re
from typing import Any

import pandas as pd
import pytest
from helpers import GUIDE, golden_rows, sirna_seqs

from seedmatch import (
    CustomSeed,
    Guide,
    SearchOptions,
    Seed,
    SequenceSet,
    alignment_text,
    find_matches,
)
from seedmatch.errors import InputError
from seedmatch.seed import seed_window

GUIDE_LENGTH = len(GUIDE)
LINE_FILL = 14
FIRST_MATCH = {"seq_id": "ENSRNOG00000016275", "start": 73, "end": 79}
GUIDE_LINE = "3'-UUUUGUCACAAGAACGAGAUAUU-5'"


def alignment_cases() -> list[Any]:
    cases: dict[str, list[Any]] = {}
    for row in golden_rows("analyses/alignment.tsv"):
        cases.setdefault(row.case, []).append(row)
    return [sorted(rows, key=lambda row: int(row.line)) for rows in cases.values()]


@pytest.mark.parametrize("rows", alignment_cases(), ids=lambda rows: rows[0].case)
def test_alignment_text_matches_r(rows):
    first = rows[0]
    assert [row.line for row in rows] == ["1", "2", "3"]
    if first.seed == "custom":
        seed = CustomSeed("custom", int(first.start_pos), int(first.stop_pos))
    else:
        seed = Seed(first.seed)
    assert seed_window(seed, GUIDE_LENGTH) == (int(first.start_pos), int(first.stop_pos))
    match = {"seq_id": first.seq_id, "start": int(first.start), "end": int(first.end)}
    text = alignment_text(
        match,
        Guide(GUIDE),
        sirna_seqs(),
        seed=seed,
        left_pad=int(first.left_pad),
        right_pad=int(first.right_pad),
    )
    assert text.split("\n") == [row.text for row in rows]


def test_alignment_text_shows_the_full_guide_D045():
    """D-045: the unpatched R stops for this match with `nchar(s1) == nchar(s2) is not TRUE`.

    The cause is a guide target that starts at g2, thus the window has one base more. The port
    shows g1 to g23, in reverse order.
    """
    guide_line = alignment_text(FIRST_MATCH, Guide(GUIDE), sirna_seqs()).split("\n")[0]
    assert guide_line.strip() == GUIDE_LINE
    assert guide_line.strip()[3:-3] == GUIDE[::-1]
    assert len(guide_line.strip()[3:-3]) == 23


def test_alignment_text_puts_spaces_outside_the_sequence_D046():
    """D-046: R stops with a `subseq` error for a window that is not inside the sequence."""
    target_at_3 = "GG" + "CTCTATA" + "AAA"
    target_at_end = "C" * 20 + "CTCTATA"
    seqs = SequenceSet.from_records([("s1", target_at_3), ("s2", target_at_end)])
    guide = Guide(GUIDE)
    before = alignment_text({"seq_id": "s1", "start": 3, "end": 9}, guide, seqs).split("\n")
    after = alignment_text({"seq_id": "s2", "start": 21, "end": 27}, guide, seqs).split("\n")
    assert before == [
        f"{' ' * 5}{GUIDE_LINE}{' ' * 3}",
        f"{' ' * 8}{' ' * 14}{'|' * 9}{' ' * 6}",
        f"  5'-...{' ' * 13}GGCUCUAUAA...-3'",
    ]
    assert after == [
        f"{' ' * 5}{GUIDE_LINE}{' ' * 3}",
        f"{' ' * 8}{' ' * 4}|{' ' * 6}|{' ' * 3}{'|' * 7} {' ' * 6}",
        f"  5'-...{'C' * 15}CUCUAUA ...-3'",
    ]
    assert {len(line) for line in before + after} == {GUIDE_LENGTH + LINE_FILL}


def test_alignment_text_pairs_nothing_in_the_padding_of_a_window_past_the_start_D046():
    seqs = SequenceSet.from_records([("s1", "GG" + "CTCTATA" + "AAA")])
    lines = alignment_text(
        {"seq_id": "s1", "start": 3, "end": 9}, Guide(GUIDE), seqs, left_pad=3, right_pad=2
    ).split("\n")
    assert lines == [
        f"{' ' * 8}{GUIDE_LINE}{' ' * 5}",
        f"{' ' * 8}{' ' * 3}{' ' * 14}{'|' * 9}{' ' * 2}{' ' * 6}",
        f"  5'-...{' ' * 16}GGCUCUAUAAAA...-3'",
    ]
    assert {len(line) for line in lines} == {GUIDE_LENGTH + 3 + 2 + LINE_FILL}


@pytest.mark.parametrize(("start", "end"), [(-30, -24), (40, 46)])
def test_alignment_text_gives_only_spaces_to_a_window_with_no_base_D046(start, end):
    seqs = SequenceSet.from_records([("s1", "A" * 10)])
    lines = alignment_text({"seq_id": "s1", "start": start, "end": end}, Guide(GUIDE), seqs).split(
        "\n"
    )
    assert lines[1:] == [" " * (GUIDE_LENGTH + LINE_FILL), f"  5'-...{' ' * GUIDE_LENGTH}...-3'"]


def test_alignment_text_accepts_a_row_of_find_matches():
    seqs = sirna_seqs()
    guide = Guide(GUIDE)
    row = find_matches(guide.target(Seed.MER7M8), seqs).iloc[0]
    expected = alignment_text(FIRST_MATCH, guide, seqs)
    assert alignment_text(row, guide, seqs) == expected
    assert alignment_text(dict(row), guide, seqs) == expected
    assert expected.split("\n")[2] == "  5'-...GAAACAGUGUUCUUGCUCUAUAA...-3'"


def test_alignment_text_accepts_a_mismatch_match_before_the_start_D046():
    guide = Guide(GUIDE)
    seqs = SequenceSet.from_records([("s1", "TCTATAGGGG")])
    row = find_matches(guide.target(Seed.MER7M8), seqs, SearchOptions(max_mismatch=1)).iloc[0]
    assert row["start"] == 0
    assert alignment_text(row, guide, seqs).split("\n") == [
        f"{' ' * 5}{GUIDE_LINE}{' ' * 3}",
        f"{' ' * 8}{' ' * 16}{'|' * 6} {' ' * 6}",
        f"  5'-...{' ' * 16}UCUAUAG...-3'",
    ]


def test_alignment_text_accepts_a_mismatch_match_after_the_end_D046():
    guide = Guide(GUIDE)
    seqs = SequenceSet.from_records([("s1", "CTCTAT")])
    row = find_matches(guide.target(Seed.MER7M8), seqs, SearchOptions(max_mismatch=1)).iloc[0]
    assert (row["start"], row["end"]) == (1, 7)
    assert alignment_text(row, guide, seqs).split("\n") == [
        f"{' ' * 5}{GUIDE_LINE}{' ' * 3}",
        f"{' ' * 8}{' ' * 15}{'|' * 6}{' ' * 2}{' ' * 6}",
        f"  5'-...{' ' * 15}CUCUAU{' ' * 2}...-3'",
    ]


def ids_of(*names):
    return SequenceSet.from_records([(name, "A" * 100) for name in names])


@pytest.mark.parametrize("end", [78, 80])
def test_alignment_text_rejects_an_indel_match(end):
    message = (
        f"match 'g1' 73-{end} has width {end - 72}, but the window of mer7m8 has width 7; "
        "expected a match without an insertion or a deletion"
    )
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        alignment_text({"seq_id": "g1", "start": 73, "end": end}, Guide(GUIDE), ids_of("g1"))


def test_alignment_text_names_the_seed_of_a_custom_window_in_the_width_message():
    message = (
        "match 'g1' 73-79 has width 7, but the window of g3-g12 has width 10; "
        "expected a match without an insertion or a deletion"
    )
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        alignment_text(
            {"seq_id": "g1", "start": 73, "end": 79},
            Guide(GUIDE),
            ids_of("g1"),
            seed=CustomSeed("g3-g12", 3, 12),
        )


def test_alignment_text_rejects_an_unknown_seq_id():
    message = "seq_id 'g9' is not an ID of the sequences; expected the ID of one sequence"
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        alignment_text({"seq_id": "g9", "start": 73, "end": 79}, Guide(GUIDE), ids_of("g1"))


def test_alignment_text_rejects_a_repeated_seq_id():
    message = "seq_id 'g1' names 2 sequences; expected the ID of one sequence"
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        alignment_text(
            {"seq_id": "g1", "start": 73, "end": 79}, Guide(GUIDE), ids_of("g1", "g2", "g1")
        )


@pytest.mark.parametrize("name", ["left_pad", "right_pad"])
def test_alignment_text_rejects_a_negative_pad(name):
    message = f"{name} -1 is negative; expected 0 or more"
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        alignment_text(FIRST_MATCH, Guide(GUIDE), sirna_seqs(), **{name: -1})


SHORT_SEQUENCE = "GG" + "CTCTATA" + "AAA"


@pytest.mark.parametrize("name", ["left_pad", "right_pad"])
def test_alignment_text_accepts_a_pad_equal_to_the_length_of_the_sequence(name):
    seqs = SequenceSet.from_records([("s1", SHORT_SEQUENCE)])
    size = len(SHORT_SEQUENCE)
    text = alignment_text(
        {"seq_id": "s1", "start": 3, "end": 9}, Guide(GUIDE), seqs, **{name: size}
    )
    assert {len(line) for line in text.split("\n")} == {GUIDE_LENGTH + size + LINE_FILL}


@pytest.mark.parametrize(
    ("name", "pad", "sequence"),
    [
        ("left_pad", 13, SHORT_SEQUENCE),
        ("right_pad", 13, SHORT_SEQUENCE),
        ("left_pad", 500, SHORT_SEQUENCE),
        ("right_pad", 2**62, SHORT_SEQUENCE),
        ("left_pad", 1, ""),
    ],
)
def test_alignment_text_rejects_a_pad_longer_than_the_sequence_D046(name, pad, sequence):
    seqs = SequenceSet.from_records([("s1", sequence)])
    size = len(sequence)
    message = (
        f"{name} {pad} is greater than the length {size} of sequence 's1'; expected at most {size}"
    )
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        alignment_text({"seq_id": "s1", "start": 3, "end": 9}, Guide(GUIDE), seqs, **{name: pad})


@pytest.mark.parametrize("key", ["seq_id", "start", "end"])
def test_alignment_text_rejects_a_match_without_start(key):
    match = {name: value for name, value in FIRST_MATCH.items() if name != key}
    message = f"match has no {key!r}; expected a row of find_matches with seq_id, start, and end"
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        alignment_text(match, Guide(GUIDE), sirna_seqs())


@pytest.mark.parametrize(
    ("match", "type_name"),
    [
        (["ENSRNOG00000016275", 73, 79], "list"),
        (pd.DataFrame([FIRST_MATCH]), "DataFrame"),
        (None, "NoneType"),
    ],
)
def test_alignment_text_rejects_a_match_that_is_not_a_row(match, type_name):
    message = (
        f"match has type {type_name}; expected a row of find_matches, a pandas.Series or a mapping"
    )
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        alignment_text(match, Guide(GUIDE), sirna_seqs())


def test_alignment_text_rejects_a_start_that_is_not_an_integer():
    match = FIRST_MATCH | {"start": 73.0}
    message = "start 73.0 has type float; expected an integer"
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        alignment_text(match, Guide(GUIDE), sirna_seqs())


def test_alignment_text_rejects_a_seed_that_does_not_fit_the_guide():
    message = "18mer needs a guide of at least 19 nt, got 9"
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        alignment_text(FIRST_MATCH, Guide("UUAUAGAGC"), sirna_seqs(), seed=Seed.MER18)
