"""The types and the input checks that the tests of `seedmatch.stats` share."""

from dataclasses import dataclass
from enum import StrEnum

import numpy as np
import numpy.typing as npt
from numpy.typing import ArrayLike

from seedmatch._checks import as_number, check_one_dimensional, shown
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
    if isinstance(value, Alternative):
        return value
    expected = ", ".join(f"Alternative.{member.name}" for member in Alternative)
    message = f"alternative {shown(value)} is not valid; expected one of {expected}"
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


def _as_vector(name: str, values: ArrayLike) -> npt.NDArray[np.float64]:
    try:
        array = np.asarray(values)
    except ValueError:
        # A nested sequence with rows of different lengths has no array form without a dtype.
        array = None
    if array is not None and array.dtype.kind in "iuf" and not _holds_a_bool(values):
        check_one_dimensional(name, array)
        return array.astype(np.float64)
    # numpy turns a list of numbers and strings into strings, so only the objects name the index.
    objects = np.asarray(values, dtype=object)
    check_one_dimensional(name, objects)
    numbers = [as_number(f"{name}[{index}]", value) for index, value in enumerate(objects.tolist())]
    return np.array(numbers, dtype=np.float64)


def _holds_a_bool(values: ArrayLike) -> bool:
    """Return True if a list or a tuple holds a bool. numpy turns `[1.0, True]` into floats."""
    return isinstance(values, list | tuple) and any(
        isinstance(value, bool | np.bool_) for value in values
    )
