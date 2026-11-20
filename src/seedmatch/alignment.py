"""The text of the alignment of a guide with a seed match, as `plotAlignment` gives it."""

from collections.abc import Mapping
from typing import Final

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from seedmatch._checks import as_limit, as_text, check_not_negative, shown
from seedmatch.errors import InputError
from seedmatch.seed import Guide, Seed, SeedSpec, seed_window
from seedmatch.sequence import SequenceSet
from seedmatch.sequence.alphabet import decode_dna

_PAIRING_LETTER: Final = str.maketrans("ACGU", "TGCA")


def _field(match: pd.Series | Mapping[str, object], key: str) -> object:
    if key not in match:
        message = (
            f"match has no {key!r}; expected a row of find_matches with seq_id, start, and end"
        )
        raise InputError(message)
    return match[key]


def _position(sequences: SequenceSet, seq_id: str) -> int:
    count = sequences.ids.count(seq_id)
    if count != 1:
        problem = "is not an ID of the sequences" if count == 0 else f"names {count} sequences"
        message = f"seq_id {shown(seq_id)} {problem}; expected the ID of one sequence"
        raise InputError(message)
    return sequences.ids.index(seq_id)


def _check_pads(seq_id: str, size: int, pads: Mapping[str, int]) -> None:
    for name, pad in pads.items():
        if pad > size:
            message = (
                f"{name} {pad} is greater than the length {size} of sequence {shown(seq_id)}; "
                f"expected at most {size}"
            )
            raise InputError(message)


def _window(codes: NDArray[np.uint8], first: int, last: int) -> str:
    """Return the letters of a sequence from `first` to `last`, 1-based and inclusive.

    A position outside the sequence gives a space.
    """
    inside_first, inside_last = max(first, 1), min(last, len(codes))
    if inside_first > inside_last:
        return " " * (last - first + 1)
    letters = decode_dna(codes[inside_first - 1 : inside_last])
    return " " * (inside_first - first) + letters + " " * (last - inside_last)


# Each argument is an independent input of the text, so no parameter object fits.
def alignment_text(  # noqa: PLR0913
    match: pd.Series | Mapping[str, object],
    guide: Guide,
    sequences: SequenceSet,
    *,
    seed: SeedSpec = Seed.MER7M8,
    left_pad: int = 0,
    right_pad: int = 0,
) -> str:
    """Return the alignment of the guide with the sequence at a match, as three lines of text.

    `match` is a row of `find_matches`, a `Series` or a mapping with the keys `seq_id`, `start`, and
    `end`. `seed` is the seed of the search that found the match. Its window of the guide, from
    `seed_window`, must have the width of the match. A match with another width raises
    `InputError`. `find_matches` with `SearchOptions(indels=True)` can give such a match.

    The three lines are the guide from 3' to 5', the pairs, and the sequence at the match, with `U`
    for `T`. The sequence part has the whole guide, g1 to the last base (`D-045`), and it starts
    `left_pad` bases before and ends `right_pad` bases after the part that pairs with the guide. A
    position outside the sequence is a space (`D-046`). A `|` shows a Watson-Crick pair. A G:U pair,
    an IUPAC letter, and a space show no `|`. For a guide of N bases, each line has
    N + `left_pad` + `right_pad` + 14 characters. The lines are joined with a line end, and the text
    has no final line end.

    Each pad must be at most the length of the sequence of the match, because a larger pad shows no
    base and its text would grow without a bound. The time and the memory are O(W) for a window of W
    bases, and W is at most N + 2L for a sequence of length L. The function also scans the IDs for
    `seq_id` in O(n) for a set of n sequences, thus a loop over m matches costs O(m · n).
    """
    if not isinstance(match, pd.Series | Mapping):
        message = (
            f"match has type {type(match).__name__}; "
            "expected a row of find_matches, a pandas.Series or a mapping"
        )
        raise InputError(message)
    seq_id = as_text("seq_id", _field(match, "seq_id"))
    start = as_limit("start", _field(match, "start"))
    end = as_limit("end", _field(match, "end"))
    left = as_limit("left_pad", left_pad)
    right = as_limit("right_pad", right_pad)
    check_not_negative("left_pad", left)
    check_not_negative("right_pad", right)
    name = guide.target(seed).name
    length = len(guide.sequence)
    first, last = seed_window(seed, length)
    if end - start + 1 != last - first + 1:
        message = (
            f"match {shown(seq_id)} {start}-{end} has width {end - start + 1}, "
            f"but the window of {name} has width {last - first + 1}; "
            "expected a match without an insertion or a deletion"
        )
        raise InputError(message)
    position = _position(sequences, seq_id)
    codes = sequences.codes[sequences.offsets[position] : sequences.offsets[position + 1]]
    _check_pads(seq_id, len(codes), {"left_pad": left, "right_pad": right})
    window = _window(codes, start - (length - last) - left, end + (first - 1) + right)
    reverse = guide.sequence[::-1]
    pairs = "".join(
        "|" if letter == pairing else " "
        for letter, pairing in zip(
            window[left : left + length], reverse.translate(_PAIRING_LETTER), strict=True
        )
    )
    return "\n".join(
        (
            f"{' ' * (5 + left)}3'-{reverse}-5'{' ' * (3 + right)}",
            f"{' ' * 8}{' ' * left}{pairs}{' ' * right}{' ' * 6}",
            f"  5'-...{window.replace('T', 'U')}...-3'",
        )
    )
