"""The two-sample tests of an ECDF comparison."""

from seedmatch.stats.dts import dts_test
from seedmatch.stats.ks import ks_test
from seedmatch.stats.types import Alternative, PValueMethod, TwoSampleResult, as_sample
from seedmatch.stats.wasserstein import WassersteinAreas, wasserstein_areas
from seedmatch.stats.wilcoxon import wilcoxon_test

__all__ = [
    "Alternative",
    "PValueMethod",
    "TwoSampleResult",
    "WassersteinAreas",
    "as_sample",
    "dts_test",
    "ks_test",
    "wasserstein_areas",
    "wilcoxon_test",
]
