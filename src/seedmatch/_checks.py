"""Checks of one option value or one column, shared by the option classes and the DE functions."""

import math
import operator
from collections.abc import Sequence
from typing import SupportsIndex

import numpy as np
import pandas as pd

from seedmatch.errors import InputError

_KERNEL_INT_MAX = np.iinfo(np.int64).max


def shown(value: object) -> str:
    """Return the repr of the value for an error message.

    Python refuses the text form of an int above its digit limit, and a message must not raise
    another error. Such an int is named by its size.
    """
    try:
        return repr(value)
    except ValueError:
        if isinstance(value, int):
            return f"<an int of {value.bit_length()} bits>"
        raise


def as_limit(name: str, value: object) -> int:
    """Return the value as a plain int that the kernels can hold. A bool is not a limit."""
    if isinstance(value, bool | np.bool_) or not isinstance(value, SupportsIndex):
        raise InputError(_not_an_integer(name, value))
    try:
        limit = operator.index(value)
    except TypeError:
        # Each ndarray has __index__, but only an integer array with one element converts.
        raise InputError(_not_an_integer(name, value)) from None
    if limit > _KERNEL_INT_MAX:
        message = f"{name} {shown(limit)} is too large; expected at most 2^63 - 1"
        raise InputError(message)
    return limit


def _not_an_integer(name: str, value: object) -> str:
    return f"{name} {shown(value)} has type {type(value).__name__}; expected an integer"


def check_not_negative(name: str, value: int) -> None:
    """Raise `InputError` if the value is below 0."""
    if value < 0:
        message = f"{name} {shown(value)} is negative; expected 0 or more"
        raise InputError(message)


def as_flag(name: str, value: object) -> bool:
    """Return the value as a plain bool. An int is not a flag."""
    if isinstance(value, bool | np.bool_):
        return bool(value)
    message = f"{name} {shown(value)} has type {type(value).__name__}; expected a bool"
    raise InputError(message)


def as_cutoff(name: str, value: object) -> float:
    """Return the value as a plain float. A bool is not a cutoff, and neither is NaN."""
    if isinstance(value, bool | np.bool_) or not isinstance(
        value, int | float | np.integer | np.floating
    ):
        message = (
            f"{name} {shown(value)} has type {type(value).__name__}; expected an int or a float"
        )
        raise InputError(message)
    try:
        cutoff = float(value)
    except OverflowError:
        message = f"{name} {shown(value)} is too large; expected a number that a float can hold"
        raise InputError(message) from None
    if math.isnan(cutoff):
        message = f"{name} {shown(value)} is not a number; expected an int or a float, not NaN"
        raise InputError(message)
    return cutoff


def check_columns(de: pd.DataFrame, required: Sequence[str]) -> None:
    """Raise `InputError` that lists the missing columns in sorted order."""
    missing = sorted(set(required) - set(de.columns))
    if missing:
        message = f"DE table has missing columns: {missing}; expected {list(required)}"
        raise InputError(message)


def check_numeric_columns(de: pd.DataFrame, names: Sequence[str]) -> None:
    """Raise `InputError` that names the first column whose dtype is not numeric."""
    for name in names:
        if not pd.api.types.is_numeric_dtype(de[name]):
            message = f"column {name!r} has dtype {de[name].dtype}; expected a numeric dtype"
            raise InputError(message)


def _is_text(values: pd.Series) -> bool:
    """Return True if each value that is present is a str. A missing value is permitted."""
    present = values.dropna()
    if present.empty:
        return True
    if isinstance(present.dtype, pd.CategoricalDtype):
        present = present.cat.categories.to_series()
    return pd.api.types.infer_dtype(present) == "string"


def check_gene_ids_are_text(de: pd.DataFrame) -> None:
    """Raise `InputError` if the `gene_id` column of the DE table has values that are not text.

    The IDs of a `SequenceSet` are text, so a number column would match no gene in silence.
    """
    if not _is_text(de["gene_id"]):
        message = (
            f"'gene_id' has dtype {de['gene_id'].dtype}; expected text values. "
            "Convert the column with astype(str)"
        )
        raise InputError(message)
