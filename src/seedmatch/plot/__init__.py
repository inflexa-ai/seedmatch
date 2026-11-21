"""Figures of seedmatch results. The extra `plot` installs matplotlib."""

try:
    # The import only makes sure that matplotlib is present, so the module does not use the name.
    import matplotlib as mpl  # noqa: F401
except ImportError as error:
    message = (
        "seedmatch.plot requires matplotlib; "
        "install the extra 'plot', for example: uv pip install 'seedmatch[plot]'"
    )
    raise ImportError(message) from error

from seedmatch.plot._style import PALETTE
from seedmatch.plot.alignment import plot_alignment
from seedmatch.plot.benchmark import plot_benchmark
from seedmatch.plot.ecdf import plot_ecdf
from seedmatch.plot.kmers import plot_kmer_histogram
from seedmatch.plot.ranges import plot_ranges
from seedmatch.plot.wasserstein import plot_wasserstein

__all__ = [
    "PALETTE",
    "plot_alignment",
    "plot_benchmark",
    "plot_ecdf",
    "plot_kmer_histogram",
    "plot_ranges",
    "plot_wasserstein",
]
