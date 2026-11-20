"""The two-sample tests of an ECDF comparison."""

from seedmatch.stats.dts import dts_test
from seedmatch.stats.ecdf import EcdfTest, EcdfTestResult, assign_groups, ecdf_test
from seedmatch.stats.ks import ks_test
from seedmatch.stats.types import Alternative, PValueMethod, TwoSampleResult, as_sample
from seedmatch.stats.wasserstein import WassersteinAreas, wasserstein_areas
from seedmatch.stats.wilcoxon import wilcoxon_test

__all__ = [
    "Alternative",
    "EcdfTest",
    "EcdfTestResult",
    "PValueMethod",
    "TwoSampleResult",
    "WassersteinAreas",
    "as_sample",
    "assign_groups",
    "dts_test",
    "ecdf_test",
    "ks_test",
    "wasserstein_areas",
    "wilcoxon_test",
]
