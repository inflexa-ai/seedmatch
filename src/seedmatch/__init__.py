"""seedmatch: a Python port of SeedMatchR."""

from seedmatch.annotation import Features, FeatureType, PerGene, Reference, TranscriptFilter
from seedmatch.de import filter_de, read_deseq2_table
from seedmatch.errors import InputError, SeedMatchError
from seedmatch.report import (
    DEFAULT_SEEDS,
    NONSEED_SEEDS,
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
    "DEFAULT_SEEDS",
    "NONSEED_SEEDS",
    "REPORT_CATEGORIES",
    "CustomSeed",
    "FeatureType",
    "Features",
    "FullSearchOptions",
    "FullSearchRangesOptions",
    "Guide",
    "InputError",
    "PerGene",
    "Reference",
    "ReportOptions",
    "SearchOptions",
    "Seed",
    "SeedMatchError",
    "SeedMatchReport",
    "SeedSpec",
    "SequenceSet",
    "Target",
    "TranscriptFilter",
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
