"""Checks of the values, columns, and arrays that cross the public API into seedmatch."""

import math
import operator
from collections import Counter
from collections.abc import Callable, Iterable, Sequence
from enum import Enum
from typing import SupportsIndex

import numpy as np
import pandas as pd
from numpy.typing import NDArray

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


def check_at_least(name: str, value: int, minimum: int) -> None:
    """Raise `InputError` if the value is below the minimum."""
    if value < minimum:
        message = f"{name} {shown(value)} is less than {minimum}; expected {minimum} or more"
        raise InputError(message)


def as_flag(name: str, value: object) -> bool:
    """Return the value as a plain bool. An int is not a flag."""
    if isinstance(value, bool | np.bool_):
        return bool(value)
    message = f"{name} {shown(value)} has type {type(value).__name__}; expected a bool"
    raise InputError(message)


def as_member[M: Enum](name: str, value: object, members: type[M]) -> M:
    """Return the value if it is a member of the enum. Its text is not a member."""
    if isinstance(value, members):
        return value
    expected = ", ".join(f"{members.__name__}.{member.name}" for member in members)
    message = f"{name} {shown(value)} is not valid; expected one of {expected}"
    raise InputError(message)


def _article_names(types: tuple[type, ...]) -> str:
    """Return the names of the types for a message: `a A`, `a A or a B`, or `a A, a B, or a C`."""
    *head, last = (f"a {kind.__name__}" for kind in types)
    if not head:
        return last
    if len(head) == 1:
        return f"{head[0]} or {last}"
    return f"{', '.join(head)}, or {last}"


def check_type(name: str, value: object, expected: type | tuple[type, ...]) -> None:
    """Raise `InputError` unless the value is an instance of the expected type.

    `expected` can be a tuple of types, and then the message names each type.
    """
    if not isinstance(value, expected):
        types = expected if isinstance(expected, tuple) else (expected,)
        message = f"{name} has type {type(value).__name__}; expected {_article_names(types)}"
        raise InputError(message)


def as_text(name: str, value: object) -> str:
    """Return the value if it is a `str`."""
    if not isinstance(value, str):
        message = f"{name} {shown(value)} has type {type(value).__name__}; expected a str"
        raise InputError(message)
    return value


# Each argument is an independent part of the message or of the check, so no parameter object fits.
def as_items[T](  # noqa: PLR0913
    name: str,
    value: Iterable[object],
    as_item: Callable[[str, object], T],
    item_type: str,
    *,
    container: str = "a sequence",
    example: str | None = None,
) -> tuple[T, ...]:
    """Return the items of the value as a tuple, each item checked by `as_item`.

    A `str` is one value and not a set of values, thus it raises `InputError`, and so does a value
    that is not iterable. The message for a `str` gives the `str` in a tuple as an example, unless
    the caller gives another `example`. `as_item` gets `name[index]` and the item. The caller keeps
    its own rule for an empty set.
    """
    if isinstance(value, str):
        message = (
            f"{name} {shown(value)} is a str; expected {container} of {item_type}, "
            f"for example ({example or repr(value)},)"
        )
        raise InputError(message)
    try:
        items = tuple(value)
    except TypeError:
        message = f"{name} {shown(value)} has type {type(value).__name__}; expected {container}"
        raise InputError(message) from None
    return tuple(as_item(f"{name}[{position}]", item) for position, item in enumerate(items))


def as_generator(rng: object) -> np.random.Generator:
    """Return `rng` if it is a `Generator`, or a `Generator` that an int seed makes."""
    if isinstance(rng, np.random.Generator):
        return rng
    if isinstance(rng, bool | np.bool_) or not isinstance(rng, int | np.integer):
        message = (
            f"rng {shown(rng)} has type {type(rng).__name__}; "
            "expected a numpy.random.Generator or an int"
        )
        raise InputError(message)
    seed = operator.index(rng)
    check_not_negative("rng", seed)
    return np.random.default_rng(seed)


def as_number(name: str, value: object) -> float:
    """Return the value as a plain float. A bool is not a number. NaN and infinity are floats."""
    if isinstance(value, bool | np.bool_) or not isinstance(
        value, int | float | np.integer | np.floating
    ):
        message = (
            f"{name} {shown(value)} has type {type(value).__name__}; expected an int or a float"
        )
        raise InputError(message)
    try:
        return float(value)
    except OverflowError:
        message = f"{name} {shown(value)} is too large; expected a number that a float can hold"
        raise InputError(message) from None


def as_cutoff(name: str, value: object) -> float:
    """Return the value as a plain float. A bool is not a cutoff, and neither is NaN."""
    cutoff = as_number(name, value)
    if math.isnan(cutoff):
        message = f"{name} {shown(value)} is not a number; expected an int or a float, not NaN"
        raise InputError(message)
    return cutoff


def as_data_frame(
    value: object, *, name: str = "de", expected: str = "a pandas.DataFrame"
) -> pd.DataFrame:
    """Return `value` if it is a DataFrame, and raise `InputError` otherwise.

    `name` is the argument that the message names, and `expected` is the end of the message.
    """
    if not isinstance(value, pd.DataFrame):
        message = f"{name} has type {type(value).__name__}; expected {expected}"
        raise InputError(message)
    return value


def check_columns(
    table: object, required: Sequence[str], *, name: str = "de", label: str = "DE table"
) -> None:
    """Raise `InputError` if `table` is not a DataFrame, or if it lacks a required column.

    A required column whose label occurs more than once raises `InputError` first, because the
    label would select a table and not a column. The message for a missing column lists the missing
    columns in sorted order. `name` is the argument that the message for a wrong type names, and
    `label` names the table in the other messages.
    """
    frame = as_data_frame(table, name=name)
    check_unique_columns(frame, required, label=label)
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        message = f"{label} has missing columns: {missing}; expected {list(required)}"
        raise InputError(message)


def check_unique_columns(
    table: pd.DataFrame, names: Sequence[str], *, label: str = "DE table"
) -> None:
    """Raise `InputError` for the first of `names` that labels more than one column of `table`.

    A name that no column has is valid. `label` names the table in the message.
    """
    occurrences = Counter(table.columns)
    for column in dict.fromkeys(names):
        if occurrences[column] > 1:
            message = (
                f"{label} has {occurrences[column]} columns {column!r}; "
                "expected one column for each name"
            )
            raise InputError(message)


def check_text_columns(table: pd.DataFrame, names: Sequence[str], *, label: str) -> None:
    """Raise `InputError` for the first value of the columns `names` that is not a `str`.

    The columns come in the order of `names`, and the row is the position. `None`, NaN, and
    `pd.NA` are not a `str`. The time is O(n) for each column of n rows. The table must have the
    columns, and `label` names the table in the message.
    """
    for name in names:
        column = table[name]
        if pd.api.types.infer_dtype(column, skipna=False) == "string" and not column.hasnans:
            continue
        values = column.to_numpy(dtype=object)
        bad = next((row for row, value in enumerate(values) if not isinstance(value, str)), None)
        if bad is not None:
            message = (
                f"{label} column {name!r} holds {shown(values[bad])} in row {bad}; "
                "expected a str in each row"
            )
            raise InputError(message)


def number_text(value: float) -> str:
    """Return the value for a message: NaN as `NaN`, and a whole number without a decimal point."""
    if math.isnan(value):
        return "NaN"
    return str(int(value)) if value.is_integer() else repr(value)


def concatenate_and_clear[T: np.generic](parts: list[NDArray[T]]) -> NDArray[T]:
    """Return the parts as one array. The list is empty after the call, so the parts are freed."""
    joined = np.concatenate(parts)
    parts.clear()
    return joined


def column_values(
    de: pd.DataFrame,
    column: str,
    is_valid: Callable[[NDArray[np.float64]], NDArray[np.bool_]],
    expected: str,
) -> NDArray[np.float64]:
    """Return the numeric column as a `float64` array, with NaN for a missing value.

    `is_valid` marks the valid values. The first row that it rejects raises `InputError`. The
    message names the column, the value, and the `gene_id` of the row, and `expected` ends it. The
    table must have the column `gene_id`.
    """
    values = de[column].to_numpy(dtype=np.float64, na_value=np.nan)
    invalid = np.flatnonzero(~is_valid(values))
    if invalid.size:
        row = int(invalid[0])
        message = (
            f"column {column!r} of de holds {number_text(float(values[row]))} "
            f"for gene {shown(de['gene_id'].iloc[row])}; expected {expected}"
        )
        raise InputError(message)
    return values


def check_numeric_columns(de: pd.DataFrame, names: Sequence[str]) -> None:
    """Raise `InputError` that names the first column whose dtype is not numeric.

    As `is.numeric` of R does, this check rejects a bool column and a complex column.
    """
    for name in names:
        dtype = de[name].dtype
        is_number = (
            pd.api.types.is_numeric_dtype(dtype)
            and not pd.api.types.is_bool_dtype(dtype)
            and not pd.api.types.is_complex_dtype(dtype)
        )
        if not is_number:
            message = f"column {name!r} has dtype {de[name].dtype}; expected a numeric dtype"
            raise InputError(message)


def check_integer_columns(table: pd.DataFrame, names: Sequence[str]) -> None:
    """Raise `InputError` that names the first column that is not an integer column.

    A bool column is not an integer column. A column of a nullable integer dtype must have a value
    in each row.
    """
    for name in names:
        column = table[name]
        if not pd.api.types.is_integer_dtype(column.dtype):
            message = f"column {name!r} has dtype {column.dtype}; expected an integer dtype"
            raise InputError(message)
        if column.hasnans:
            message = f"column {name!r} has missing values; expected an integer in each row"
            raise InputError(message)


def is_text(values: pd.Series) -> bool:
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
    if not is_text(de["gene_id"]):
        message = (
            f"'gene_id' has dtype {de['gene_id'].dtype}; expected text values. "
            "Convert the column with astype(str)"
        )
        raise InputError(message)


def check_one_dimensional(name: str, array: NDArray[np.generic]) -> None:
    """Raise `InputError` if the array is not 1-D."""
    if array.ndim != 1:
        message = f"{name} has {array.ndim} dimensions; expected a 1-D array"
        raise InputError(message)


def check_vector(name: str, array: NDArray[np.generic], dtype: type[np.generic]) -> None:
    """Raise `InputError` if the array is not a 1-D array of the dtype."""
    if array.dtype != dtype:
        message = f"{name} has dtype {array.dtype}; expected {np.dtype(dtype)}"
        raise InputError(message)
    check_one_dimensional(name, array)


def check_offsets(
    offsets: NDArray[np.int64], count: int, unit: str, end: int, end_name: str
) -> None:
    """Raise `InputError` unless item i of the `count` items is `offsets[i]:offsets[i + 1]`.

    The offsets must go from 0 to `end` without a decrease. `unit` names the items, and `end_name`
    names the value of `end` in the message.
    """
    if len(offsets) != count + 1:
        message = f"offsets has {len(offsets)} entries; expected {count + 1} for {count} {unit}"
        raise InputError(message)
    if offsets[0] != 0:
        message = f"offsets must start at 0; got {offsets[0]}"
        raise InputError(message)
    if offsets[-1] != end:
        message = f"offsets must end at {end}, {end_name}; got {offsets[-1]}"
        raise InputError(message)
    decreases = np.flatnonzero(np.diff(offsets) < 0)
    if decreases.size:
        i = int(decreases[0]) + 1
        message = (
            f"offsets[{i}] {offsets[i]} is less than offsets[{i - 1}] {offsets[i - 1]}; "
            "expected offsets that do not decrease"
        )
        raise InputError(message)
