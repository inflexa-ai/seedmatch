"""The letter compatibility of a Biostrings search, for each pair of `fixed` flags."""

import numpy as np
from numpy.typing import NDArray

_MASK_COUNT = 16


def compat_table(*, fixed_pattern: bool, fixed_subject: bool) -> NDArray[np.bool_]:
    """Return the (16, 16) table: `[p, s]` is true if pattern mask `p` matches subject mask `s`.

    A fixed side takes its letter as a literal. A side that is not fixed takes its letter as the
    set of its bases:

    - both sides fixed: the letters are equal
    - only the pattern fixed: the bases of `p` are a subset of the bases of `s`
    - only the subject fixed: the bases of `s` are a subset of the bases of `p`
    - no side fixed: the letters have a base in common

    Mask 0 is not a letter, so its row and its column are false.
    """
    p = np.arange(_MASK_COUNT)[:, np.newaxis]
    s = np.arange(_MASK_COUNT)[np.newaxis, :]
    table: NDArray[np.bool_]
    if fixed_pattern and fixed_subject:
        table = p == s
    elif fixed_pattern:
        table = (p & ~s) == 0
    elif fixed_subject:
        table = (s & ~p) == 0
    else:
        table = (p & s) != 0
    table[0, :] = False
    table[:, 0] = False
    return table
