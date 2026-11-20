"""The types and the input checks that the tests of `seedmatch.stats` share."""

from dataclasses import dataclass
from enum import StrEnum
from typing import Any

import numpy as np
import numpy.typing as npt
from numpy.typing import ArrayLike

from seedmatch._checks import as_member, as_number, check_one_dimensional
from seedmatch.errors import InputError


class Alternative(StrEnum):
    """The alternative hypothesis, with the values of the `alternative` argument of R."""

    TWO_SIDED = "two.sided"
    LESS = "less"
    GREATER = "greater"


class PValueMethod(StrEnum):
    """How a p-value comes from the statistic."""

    EXACT = "exact"
    ASYMPTOTIC = "asymptotic"
    PERMUTATION = "permutation"


@dataclass(frozen=True)
class TwoSampleResult:
    """The result of a two-sample test. `alternative` is `None` for a test that has no direction."""

    statistic: float
    p_value: float
    method: PValueMethod
    alternative: Alternative | None
    n_x: int
    n_y: int


def as_alternative(value: object) -> Alternative:
    """Return the value if it is a member of `Alternative`. Its text is not a member."""
    return as_member("alternative", value, Alternative)


def check_range(
    first: str,
    second: str,
    x: npt.NDArray[np.float64],
    y: npt.NDArray[np.float64],
    maximum: float,
) -> None:
    """Raise `InputError` if the two samples together span more than `maximum`.

    A range that overflows to infinity is above each maximum. The message names the two samples.
    """
    low = float(min(x.min(), y.min()))
    high = float(max(x.max(), y.max()))
    if high - low > maximum:
        message = (
            f"{first} and {second} span from {low} to {high}; expected a range of at most {maximum}"
        )
        raise InputError(message)


def as_sample(name: str, values: ArrayLike, *, finite: bool) -> npt.NDArray[np.float64]:
    """Return the values as a new 1-D `float64` array.

    An empty sample, a value that is not a number, and NaN raise `InputError`, and so does an
    infinite value if `finite` is true. A message names the sample and the index of the bad value.
    A string is not a number, even if it holds digits.
    """
    sample = _as_vector(name, values)
    if sample.size == 0:
        message = f"{name} is empty; expected at least 1 value"
        raise InputError(message)
    missing = np.isnan(sample)
    if missing.any():
        message = f"{name}[{missing.argmax()}] is NaN; expected a number"
        raise InputError(message)
    infinite = np.isinf(sample)
    if finite and infinite.any():
        index = infinite.argmax()
        message = f"{name}[{index}] {sample[index]} is infinite; expected a finite number"
        raise InputError(message)
    return sample


def as_vector(name: str, values: ArrayLike, *, kinds: str, expected: str) -> npt.NDArray[Any]:
    """Return the values as a 1-D array of a dtype kind in `kinds`, or else as a 1-D object array.

    A value with no 1-D form (a set, a generator, a dict, a str, or a scalar) raises `InputError`.
    The message names the type of the value, and `expected` names the items of a 1-D array.
    """
    try:
        array = np.asarray(values)
    except ValueError:
        # A nested sequence with rows of different lengths has no array form without a dtype.
        array = None
    if array is not None and array.ndim == 0:
        # A set, a generator, a dict, a str, and a scalar each make a 0-D array.
        message = f"{name} has type {type(values).__name__}; expected a 1-D array of {expected}"
        raise InputError(message)
    if array is None or array.dtype.kind not in kinds:
        # numpy turns a list of numbers and strings into strings, so only objects name the index.
        array = np.asarray(values, dtype=object)
    check_one_dimensional(name, array)
    return array


def _as_vector(name: str, values: ArrayLike) -> npt.NDArray[np.float64]:
    array = as_vector(name, values, kinds="iuf", expected="numbers")
    if array.dtype.kind in "iuf" and not _holds_a_bool(values):
        return array.astype(np.float64)
    if array.dtype.kind != "O":
        # numpy made floats of a list with a bool, so only the objects name the bool.
        array = np.asarray(values, dtype=object)
        check_one_dimensional(name, array)
    numbers = [as_number(f"{name}[{index}]", value) for index, value in enumerate(array.tolist())]
    return np.array(numbers, dtype=np.float64)


def _holds_a_bool(values: ArrayLike) -> bool:
    """Return True if a list or a tuple holds a bool. numpy turns `[1.0, True]` into floats."""
    return isinstance(values, list | tuple) and any(
        isinstance(value, bool | np.bool_) for value in values
    )
