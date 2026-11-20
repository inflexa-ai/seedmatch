"""Errors that seedmatch raises."""


class SeedMatchError(Exception):
    """Base class of each error that seedmatch raises."""


class InputError(SeedMatchError, ValueError):
    """Invalid input. The message names the bad value and the expected value."""


class SeedMatchWarning(UserWarning):
    """Base class of each warning that seedmatch gives."""
