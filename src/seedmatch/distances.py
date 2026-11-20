"""The distances between the matches of a sequence, as `match_distances` of R gives them."""

from enum import StrEnum
from typing import Final, assert_never

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from seedmatch._checks import as_data_frame, as_member, check_columns, check_integer_columns

_COLUMNS: Final = ("seq_id", "gene_id", "start", "end", "width")
_INTEGER_COLUMNS: Final = ("start", "end", "width")


class MatchAnchor(StrEnum):
    """The position of a match that the distance uses, with the values of `feature.pos` of R."""

    FIVE_PRIME_END = "fivePrimeEnd"
    THREE_PRIME_END = "threePrimeEnd"
    CENTER = "center"


def _coordinates(matches: pd.DataFrame, anchor: MatchAnchor) -> NDArray[np.float64]:
    match anchor:
        case MatchAnchor.FIVE_PRIME_END:
            return matches["end"].to_numpy(dtype=np.float64)
        case MatchAnchor.THREE_PRIME_END:
            return matches["start"].to_numpy(dtype=np.float64)
        case MatchAnchor.CENTER:
            starts = matches["start"].to_numpy(dtype=np.float64)
            return starts + matches["width"].to_numpy(dtype=np.float64) / 2
        case _:
            assert_never(anchor)


def match_distances(
    matches: pd.DataFrame, *, anchor: MatchAnchor = MatchAnchor.FIVE_PRIME_END
) -> pd.DataFrame:
    """Return the distance between each pair of matches in a sequence, as `match_distances` of R.

    `matches` has the columns `seq_id`, `gene_id`, `start`, `end`, and `width`, for example the
    result of `find_matches`. `start`, `end`, and `width` have an integer dtype. The rows with one
    `seq_id` are one sequence, also if they are not adjacent. All rows with a missing `seq_id` are
    one sequence. The sequences come in the order of their first row, and the matches of a sequence
    keep the order of the rows.

    `anchor` sets the coordinate of a match: `end` for `FIVE_PRIME_END`, `start` for
    `THREE_PRIME_END`, and `start + width / 2` for `CENTER`. For the coordinates c_1 to c_n of a
    sequence, the pairs come in the R order: (1, 2), (1, 3), ..., (1, n), (2, 3), and so on. The
    distance of a pair is |c_i - c_j|.

    The result has the columns `seq_id`, `gene_id` (of the first row of the sequence), and
    `distance`, and it has a `RangeIndex`. `distance` is `float64`, because `CENTER` gives half
    values. A sequence with one match has no pair, thus it gives no row. R gives the values `NA`
    and 0 for it (`D-011`). An empty `matches` gives an empty result.

    A sequence with n matches gives n(n - 1)/2 rows. No bound below the size of the result is
    possible: 10,000 matches in one sequence give 49,995,000 rows. A row takes 24 bytes, because
    `distance` and each of the two ID columns take 8 bytes. The function counts the pairs of each
    sequence first and allocates the result one time. It fills the pairs of one match with one
    vectorized difference, and it writes the difference into the result. Thus the work memory,
    without the result, is O(m) for m matches. The time is O(m log m + P) for P pairs, because of
    the stable sort of the rows, and there is one numpy call for each match that has a later match
    in its sequence.

    Measured on an Apple M2 Max: 10,000 matches in one sequence take 0.3 s. The result holds
    1.2 GB, and the peak memory of the process is 1.4 GB. 1,000,000 matches in 500,000 sequences of
    2 matches give 500,000 rows in 0.6 s, with a peak memory of 0.35 GB.
    """
    table = as_data_frame(matches, name="matches")
    check_columns(table, _COLUMNS, name="matches", label="matches table")
    check_integer_columns(table, _INTEGER_COLUMNS)
    kind = as_member("anchor", anchor, MatchAnchor)

    # All rows with a missing seq_id are one sequence, so each row has a code of 0 or more.
    # np.bincount rejects the code -1 that the default of factorize gives to a missing value.
    codes, ids = pd.factorize(table["seq_id"], use_na_sentinel=False)
    sizes = np.bincount(codes, minlength=len(ids))
    order = np.argsort(codes, kind="stable")
    coordinates = _coordinates(table, kind)[order]
    # In the sorted rows, a match pairs with the later matches of its sequence.
    later = np.repeat(np.cumsum(sizes), sizes) - np.arange(len(order)) - 1
    pairs = sizes * (sizes - 1) // 2

    distance = np.empty(int(pairs.sum()), dtype=np.float64)
    position = 0
    rows = np.flatnonzero(later)
    for row, count in zip(rows.tolist(), later[rows].tolist(), strict=True):
        slot = distance[position : position + count]
        np.subtract(coordinates[row + 1 : row + 1 + count], coordinates[row], out=slot)
        np.abs(slot, out=slot)
        position += count

    # The sort is stable, so the first sorted row of a sequence is its first row in the table.
    first_rows = order[np.cumsum(sizes) - sizes]
    # copy=False keeps the peak memory at the size of the result. A copy doubles it.
    return pd.DataFrame(
        {
            "seq_id": np.repeat(ids.to_numpy(dtype=object), pairs),
            "gene_id": np.repeat(table["gene_id"].iloc[first_rows].to_numpy(dtype=object), pairs),
            "distance": distance,
        },
        copy=False,
    )
