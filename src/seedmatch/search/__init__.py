"""Pattern search with the matching rules of Biostrings."""

from seedmatch.search.api import (
    MAX_PATTERN_LENGTH,
    SearchOptions,
    count_matches,
    find_matches,
    match_ranges,
)

__all__ = ["MAX_PATTERN_LENGTH", "SearchOptions", "count_matches", "find_matches", "match_ranges"]
