"""seedmatch: a Python port of SeedMatchR.

A `D-NNN` ID in the code names an entry of `divergences.md` in this package, a change of behavior
from SeedMatchR. `r-map.md` in this package gives the Python name of each SeedMatchR function. Read
them with `importlib.resources.files("seedmatch").joinpath("divergences.md").read_text()`.
"""

from seedmatch.alignment import alignment_text
from seedmatch.annotation import Features, FeatureType, PerGene, Reference, TranscriptFilter
from seedmatch.de import filter_de, read_deseq2_table
from seedmatch.distances import MatchAnchor, match_distances
from seedmatch.errors import InputError, SeedMatchError, SeedMatchWarning
from seedmatch.evalue import e_value_null, e_values
from seedmatch.kmers import KmerCounts, kmer_counts, kmer_gene_counts, kmer_histogram
from seedmatch.mirdb import read_mirdb
from seedmatch.network import (
    Aggregator,
    GeneNetwork,
    StringLinks,
    read_string_links,
    read_trrust,
    string_gene_network,
    string_scores,
    trrust_scores,
)
from seedmatch.noncanonical import (
    CompensatoryOptions,
    WobbleBulgeOptions,
    annotate_wobble_bulge_counts,
    compensatory_sites,
    wobble_bulge_matches,
    wobble_bulge_null,
    wobble_bulge_patterns,
)
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
    report_html,
    seed_match_report,
)
from seedmatch.search import SearchOptions, count_matches, find_matches
from seedmatch.seed import CustomSeed, Guide, Seed, SeedSpec, Target
from seedmatch.sequence import SequenceSet, read_fasta, reverse_complement, write_fasta
from seedmatch.stats import (
    Alternative,
    Benchmark,
    EcdfTest,
    EcdfTestResult,
    PrCurve,
    PValueMethod,
    RocCurve,
    WassersteinAreas,
    assign_groups,
    benchmark_predictions,
    ecdf_test,
    pr_curve,
    roc_curve,
    wasserstein_areas,
    wasserstein_difference,
)

__all__ = [
    "DEFAULT_SEEDS",
    "NONSEED_SEEDS",
    "REPORT_CATEGORIES",
    "Aggregator",
    "Alternative",
    "Benchmark",
    "CompensatoryOptions",
    "CustomSeed",
    "EcdfTest",
    "EcdfTestResult",
    "FeatureType",
    "Features",
    "FullSearchOptions",
    "FullSearchRangesOptions",
    "GeneNetwork",
    "Guide",
    "InputError",
    "KmerCounts",
    "MatchAnchor",
    "PValueMethod",
    "PerGene",
    "PrCurve",
    "Reference",
    "ReportOptions",
    "RocCurve",
    "SearchOptions",
    "Seed",
    "SeedMatchError",
    "SeedMatchReport",
    "SeedMatchWarning",
    "SeedSpec",
    "SequenceSet",
    "StringLinks",
    "Target",
    "TranscriptFilter",
    "WassersteinAreas",
    "WobbleBulgeOptions",
    "alignment_text",
    "annotate_seed_counts",
    "annotate_wobble_bulge_counts",
    "assign_groups",
    "benchmark_predictions",
    "best_match",
    "compensatory_sites",
    "count_matches",
    "e_value_null",
    "e_values",
    "ecdf_test",
    "filter_de",
    "find_matches",
    "full_search",
    "full_search_ranges",
    "kmer_counts",
    "kmer_gene_counts",
    "kmer_histogram",
    "match_distances",
    "pr_curve",
    "read_deseq2_table",
    "read_fasta",
    "read_mirdb",
    "read_string_links",
    "read_trrust",
    "report_html",
    "reverse_complement",
    "roc_curve",
    "seed_match_report",
    "string_gene_network",
    "string_scores",
    "trrust_scores",
    "wasserstein_areas",
    "wasserstein_difference",
    "wobble_bulge_matches",
    "wobble_bulge_null",
    "wobble_bulge_patterns",
    "write_fasta",
]
