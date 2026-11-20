"""The public search: count the matches of a pattern, or find their ranges.

The functions obey the search rules of `docs/design.md`. The counts and the ranges of one call come
from the same kernel, so they agree.
"""

from dataclasses import dataclass
from typing import Final

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from seedmatch._checks import as_flag, as_limit, check_not_negative
from seedmatch.errors import InputError
from seedmatch.search import kernels
from seedmatch.search.compat import compat_table
from seedmatch.seed import Target
from seedmatch.sequence import SequenceSet
from seedmatch.sequence.alphabet import encode_dna

# Biostrings 2.78.0 stops for a longer pattern in each search mode. The limit also bounds the cost
# of a search with a large max_mismatch.
_MAX_PATTERN_LENGTH: Final = 20000


@dataclass(frozen=True)
class SearchOptions:
    """The options of a search, with the names of `matchPattern` in snake case.

    `fixed_pattern` and `fixed_subject` are the two entries of the `fixed` argument of R. With
    `indels` on and `max_mismatch` 0, the search is exact. `min_mismatch` applies only if `indels`
    is off.
    """

    max_mismatch: int = 0
    min_mismatch: int = 0
    indels: bool = False
    fixed_pattern: bool = True
    fixed_subject: bool = True

    def __post_init__(self) -> None:
        """Make sure that the options are valid, because the kernels trust them."""
        # Plain int limits give each kernel one type signature, thus one numba cache entry.
        object.__setattr__(self, "max_mismatch", as_limit("max_mismatch", self.max_mismatch))
        object.__setattr__(self, "min_mismatch", as_limit("min_mismatch", self.min_mismatch))
        object.__setattr__(self, "indels", as_flag("indels", self.indels))
        object.__setattr__(self, "fixed_pattern", as_flag("fixed_pattern", self.fixed_pattern))
        object.__setattr__(self, "fixed_subject", as_flag("fixed_subject", self.fixed_subject))
        for name, value in (
            ("max_mismatch", self.max_mismatch),
            ("min_mismatch", self.min_mismatch),
        ):
            check_not_negative(name, value)
        if self.min_mismatch > self.max_mismatch:
            message = (
                f"min_mismatch {self.min_mismatch} is greater than "
                f"max_mismatch {self.max_mismatch}; expected at most {self.max_mismatch}"
            )
            raise InputError(message)
        if self.indels and self.min_mismatch > 0:
            message = f"min_mismatch must be 0 when indels is on; got {self.min_mismatch}"
            raise InputError(message)


def _encode_pattern(pattern: str | Target) -> NDArray[np.uint8]:
    text = pattern.sequence if isinstance(pattern, Target) else pattern
    if len(text) > _MAX_PATTERN_LENGTH:
        message = f"pattern has {len(text)} letters; expected at most {_MAX_PATTERN_LENGTH}"
        raise InputError(message)
    masks = encode_dna(text)
    if masks.size == 0:
        message = "empty pattern; expected a pattern of one letter or more"
        raise InputError(message)
    return masks


def _compat(options: SearchOptions) -> NDArray[np.bool_]:
    return compat_table(fixed_pattern=options.fixed_pattern, fixed_subject=options.fixed_subject)


def _uses_indel_kernel(options: SearchOptions) -> bool:
    return options.indels and options.max_mismatch > 0


def count_matches(
    pattern: str | Target,
    sequences: SequenceSet,
    options: SearchOptions = SearchOptions(),
) -> NDArray[np.int64]:
    """Return the number of matches of the pattern in each sequence, in the order of the set."""
    masks = _encode_pattern(pattern)
    compat = _compat(options)
    if _uses_indel_kernel(options):
        counts = kernels.indel_counts(
            sequences.codes, sequences.offsets, masks, compat, options.max_mismatch
        )
    else:
        counts = kernels.hamming_counts(
            sequences.codes,
            sequences.offsets,
            masks,
            compat,
            options.max_mismatch,
            options.min_mismatch,
        )
    return np.asarray(counts, dtype=np.int64)


def match_ranges(
    pattern: str | Target,
    sequences: SequenceSet,
    options: SearchOptions = SearchOptions(),
) -> tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.int64]]:
    """Return the sequence index, the start, and the end of each match, 1-based and inclusive.

    The arrays have the rows of `find_matches`, in its order. The sequence index is the position of
    the sequence in the set, thus a sequence ID that occurs more than once does not hide which
    sequence has the match. The ranges follow the rules of `find_matches`.
    """
    masks = _encode_pattern(pattern)
    compat = _compat(options)
    if _uses_indel_kernel(options):
        index, starts, ends = kernels.indel_ranges(
            sequences.codes, sequences.offsets, masks, compat, options.max_mismatch
        )
    else:
        index, starts = kernels.hamming_ranges(
            sequences.codes,
            sequences.offsets,
            masks,
            compat,
            options.max_mismatch,
            options.min_mismatch,
        )
        ends = starts + (masks.size - 1)
    return (
        np.asarray(index, dtype=np.int64),
        np.asarray(starts, dtype=np.int64),
        np.asarray(ends, dtype=np.int64),
    )


def find_matches(
    pattern: str | Target,
    sequences: SequenceSet,
    options: SearchOptions = SearchOptions(),
) -> pd.DataFrame:
    """Return the ranges of the matches, 1-based and inclusive.

    The columns are `seq_id`, `gene_id`, `start`, `end`, and `width`. The rows come in the order of
    the sequences, and then in the order of Biostrings. A mismatch match can reach past either end
    of a sequence, so `start` can be less than 1 and `end` can be greater than the length.
    """
    index, starts, ends = match_ranges(pattern, sequences, options)
    seq_ids = np.array(sequences.ids, dtype=object)
    gene_ids = np.array(sequences.gene_ids, dtype=object)
    return pd.DataFrame(
        {
            "seq_id": seq_ids[index],
            "gene_id": gene_ids[index],
            "start": starts,
            "end": ends,
            "width": ends - starts + 1,
        }
    )
