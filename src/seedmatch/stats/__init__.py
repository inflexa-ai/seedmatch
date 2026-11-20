"""The two-sample tests of an ECDF comparison."""

from seedmatch.stats.ks import ks_test
from seedmatch.stats.types import Alternative, PValueMethod, TwoSampleResult, as_sample

__all__ = ["Alternative", "PValueMethod", "TwoSampleResult", "as_sample", "ks_test"]
