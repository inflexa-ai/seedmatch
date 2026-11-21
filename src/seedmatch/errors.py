"""Errors that seedmatch raises, and the warning that it gives."""

import os
import warnings
from pathlib import Path


class SeedMatchError(Exception):
    """Base class of each error that seedmatch raises."""


class InputError(SeedMatchError, ValueError):
    """Invalid input. The message names the bad value and the expected value."""


class FetchError(SeedMatchError):
    """A download failed. The message names the URL or the path, and the value that was expected."""


class SeedMatchWarning(UserWarning):
    """Base class of each warning that seedmatch gives."""


_PACKAGE_DIRECTORY = str(Path(__file__).parent) + os.sep


def warn(message: str) -> None:
    """Warn with `SeedMatchWarning`. The warning names the first line outside the package.

    The line is the same for each depth of the calls inside seedmatch, thus a caller of a public
    function sees its own line, also if a plot function calls the function that warns.
    """
    warnings.warn(message, SeedMatchWarning, skip_file_prefixes=(_PACKAGE_DIRECTORY,))
