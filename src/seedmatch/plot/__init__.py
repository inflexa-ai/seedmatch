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
from seedmatch.plot.ecdf import plot_ecdf

__all__ = ["PALETTE", "plot_ecdf"]
