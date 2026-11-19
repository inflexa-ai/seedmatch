"""The SeedMatchR modes on a DE table, and the searches over many seeds."""

from seedmatch.report.annotate import annotate_seed_counts
from seedmatch.report.grid import (
    DEFAULT_SEEDS,
    NONSEED_SEEDS,
    FullSearchOptions,
    FullSearchRangesOptions,
    best_match,
    full_search,
    full_search_ranges,
)
from seedmatch.report.report import (
    REPORT_CATEGORIES,
    ReportOptions,
    SeedMatchReport,
    seed_match_report,
)

__all__ = [
    "DEFAULT_SEEDS",
    "NONSEED_SEEDS",
    "REPORT_CATEGORIES",
    "FullSearchOptions",
    "FullSearchRangesOptions",
    "ReportOptions",
    "SeedMatchReport",
    "annotate_seed_counts",
    "best_match",
    "full_search",
    "full_search_ranges",
    "seed_match_report",
]
