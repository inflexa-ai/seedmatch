"""The two-sample tests of an ECDF comparison."""

from seedmatch.stats.ks import ks_test
from seedmatch.stats.types import Alternative, PValueMethod, TwoSampleResult, as_sample
from seedmatch.stats.wilcoxon import wilcoxon_test

__all__ = [
    "Alternative",
    "PValueMethod",
    "TwoSampleResult",
    "as_sample",
    "ks_test",
    "wilcoxon_test",
]
