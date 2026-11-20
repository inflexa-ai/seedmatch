"""Errors that seedmatch raises."""


class SeedMatchError(Exception):
    """Base class of each error that seedmatch raises."""


class InputError(SeedMatchError, ValueError):
    """Invalid input. The message names the bad value and the expected value."""


class FetchError(SeedMatchError):
    """A download failed. The message names the URL or the path, and the value that was expected."""


class SeedMatchWarning(UserWarning):
    """Base class of each warning that seedmatch gives."""
