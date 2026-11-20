"""seedmatch: a Python port of SeedMatchR."""

from seedmatch.alignment import alignment_text
from seedmatch.annotation import Features, FeatureType, PerGene, Reference, TranscriptFilter
from seedmatch.de import filter_de, read_deseq2_table
from seedmatch.distances import MatchAnchor, match_distances
from seedmatch.errors import InputError, SeedMatchError, SeedMatchWarning
from seedmatch.kmers import KmerCounts, kmer_counts, kmer_gene_counts
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
from seedmatch.stats import (
    Alternative,
    EcdfTest,
    EcdfTestResult,
    PValueMethod,
    WassersteinAreas,
    assign_groups,
    ecdf_test,
    wasserstein_areas,
)

__all__ = [
    "DEFAULT_SEEDS",
    "NONSEED_SEEDS",
    "REPORT_CATEGORIES",
    "Alternative",
    "CustomSeed",
    "EcdfTest",
    "EcdfTestResult",
    "FeatureType",
    "Features",
    "FullSearchOptions",
    "FullSearchRangesOptions",
    "Guide",
    "InputError",
    "KmerCounts",
    "MatchAnchor",
    "PValueMethod",
    "PerGene",
    "Reference",
    "ReportOptions",
    "SearchOptions",
    "Seed",
    "SeedMatchError",
    "SeedMatchReport",
    "SeedMatchWarning",
    "SeedSpec",
    "SequenceSet",
    "Target",
    "TranscriptFilter",
    "WassersteinAreas",
    "alignment_text",
    "annotate_seed_counts",
    "assign_groups",
    "best_match",
    "count_matches",
    "ecdf_test",
    "filter_de",
    "find_matches",
    "full_search",
    "full_search_ranges",
    "kmer_counts",
    "kmer_gene_counts",
    "match_distances",
    "read_deseq2_table",
    "read_fasta",
    "reverse_complement",
    "seed_match_report",
    "wasserstein_areas",
    "write_fasta",
]
