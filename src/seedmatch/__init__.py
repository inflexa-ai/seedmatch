"""seedmatch: a Python port of SeedMatchR."""

from seedmatch.de import filter_de, read_deseq2_table
from seedmatch.report import (
    REPORT_CATEGORIES,
    FullSearchOptions,
    FullSearchRangesOptions,
    ReportOptions,
    SeedMatchReport,
    annotate_seed_counts,
    best_match,
    full_search,
    full_search_ranges,
    seed_match_report,
)
from seedmatch.search import SearchOptions, count_matches, find_matches
from seedmatch.seed import CustomSeed, Guide, Seed, SeedSpec, Target
from seedmatch.sequence import SequenceSet, read_fasta, reverse_complement, write_fasta

__all__ = [
    "REPORT_CATEGORIES",
    "CustomSeed",
    "FullSearchOptions",
    "FullSearchRangesOptions",
    "Guide",
    "ReportOptions",
    "SearchOptions",
    "Seed",
    "SeedMatchReport",
    "SeedSpec",
    "SequenceSet",
    "Target",
    "annotate_seed_counts",
    "best_match",
    "count_matches",
    "filter_de",
    "find_matches",
    "full_search",
    "full_search_ranges",
    "read_deseq2_table",
    "read_fasta",
    "reverse_complement",
    "seed_match_report",
    "write_fasta",
]
