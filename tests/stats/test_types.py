import dataclasses
import re

import numpy as np
import pandas as pd
import pytest

from seedmatch.errors import InputError
from seedmatch.stats import Alternative, PValueMethod, TwoSampleResult, as_sample
from seedmatch.stats.types import as_alternative


def test_alternatives_have_the_names_of_r():
    assert {member.name: member.value for member in Alternative} == {
        "TWO_SIDED": "two.sided",
        "LESS": "less",
        "GREATER": "greater",
    }
    assert Alternative("two.sided") is Alternative.TWO_SIDED


def test_p_value_methods():
    assert {member.name: member.value for member in PValueMethod} == {
        "EXACT": "exact",
        "ASYMPTOTIC": "asymptotic",
        "PERMUTATION": "permutation",
    }


def test_two_sample_result_is_a_frozen_record_with_a_missing_alternative_for_dts():
    result = TwoSampleResult(
        statistic=0.5, p_value=0.25, method=PValueMethod.PERMUTATION, alternative=None, n_x=3, n_y=4
    )
    assert [field.name for field in dataclasses.fields(result)] == [
        "statistic",
        "p_value",
        "method",
        "alternative",
        "n_x",
        "n_y",
    ]
    assert result.alternative is None
    with pytest.raises(dataclasses.FrozenInstanceError):
        result.p_value = 0.5


def test_as_alternative_accepts_a_member():
    assert as_alternative(Alternative.LESS) is Alternative.LESS


@pytest.mark.parametrize("value", ["less", "two.sided", "Less", None, 1])
def test_as_alternative_rejects_text_and_other_values(value):
    message = (
        f"alternative {value!r} is not valid; "
        "expected one of Alternative.TWO_SIDED, Alternative.LESS, Alternative.GREATER"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        as_alternative(value)


@pytest.mark.parametrize(
    "values",
    [
        [1, 2, 3],
        (1.5, 2.5),
        np.array([1, 2], dtype=np.int32),
        np.array([1.5, 2.5], dtype=np.float32),
        pd.Series([1.5, 2.5]),
        pd.Series([1, 2], dtype="Float64"),
        range(3),
    ],
)
def test_as_sample_returns_a_float64_vector(values):
    sample = as_sample("x", values, finite=True)
    assert sample.dtype == np.float64
    assert sample.ndim == 1
    assert sample.tolist() == [float(value) for value in values]


def test_as_sample_converts_objects_that_are_numbers():
    sample = as_sample("x", pd.Series([3, 1.5], dtype=object), finite=True)
    assert sample.dtype == np.float64
    assert sample.tolist() == [3.0, 1.5]


def test_as_sample_converts_an_int_that_numpy_keeps_as_an_object():
    assert as_sample("x", [1, 2**70], finite=True).tolist() == [1.0, float(2**70)]


def test_as_sample_rejects_a_bool_in_a_list_of_numbers():
    with pytest.raises(
        InputError, match=r"^x\[1\] True has type bool; expected an int or a float$"
    ):
        as_sample("x", [1.0, True], finite=False)


def test_as_sample_rejects_a_numpy_bool_in_a_tuple_of_numbers():
    with pytest.raises(
        InputError, match=r"^x\[1\] np\.True_ has type bool; expected an int or a float$"
    ):
        as_sample("x", (2, np.True_), finite=False)


def test_as_sample_returns_a_new_array():
    values = np.array([3.0, 1.0, 2.0])
    sample = as_sample("x", values, finite=True)
    sample[0] = 99.0
    assert values.tolist() == [3.0, 1.0, 2.0]


def test_as_sample_keeps_an_infinite_value_if_it_is_not_finite_only():
    assert as_sample("x", [1.0, np.inf, -np.inf], finite=False).tolist() == [1.0, np.inf, -np.inf]


@pytest.mark.parametrize("values", [[], np.array([]), pd.Series([], dtype=float)])
def test_as_sample_rejects_an_empty_sample(values):
    with pytest.raises(InputError, match=r"^y is empty; expected at least 1 value$"):
        as_sample("y", values, finite=False)


@pytest.mark.parametrize("finite", [True, False])
def test_as_sample_rejects_nan_and_gives_its_index(finite):
    with pytest.raises(InputError, match=r"^x\[2\] is NaN; expected a number$"):
        as_sample("x", [1.0, 2.0, np.nan, 4.0], finite=finite)


def test_as_sample_rejects_a_missing_value_of_a_nullable_column():
    with pytest.raises(InputError, match=r"^x\[1\] is NaN; expected a number$"):
        as_sample("x", pd.Series([1.0, None], dtype="Float64"), finite=False)


@pytest.mark.parametrize(
    ("values", "text"),
    [
        ([1.0, "a"], r"x\[1\] 'a' has type str"),
        (["1.5", "2.5"], r"x\[0\] '1.5' has type str"),
        ([None, 1.0], r"x\[0\] None has type NoneType"),
        ([True, False], r"x\[0\] True has type bool"),
        (np.array([True, False]), r"x\[0\] .* has type bool"),
        ([1.0, 2 + 1j], r"x\[1\] \(2\+1j\) has type complex"),
        ([[1.0, 2.0], [3.0]], r"x\[0\] \[1\.0, 2\.0\] has type list"),
    ],
)
def test_as_sample_rejects_a_value_that_is_not_a_number_and_gives_its_index(values, text):
    with pytest.raises(InputError, match=f"^{text}.*; expected an int or a float"):
        as_sample("x", values, finite=False)


def test_as_sample_rejects_an_int_that_a_float_cannot_hold():
    with pytest.raises(InputError, match=r"^x\[1\] .* is too large"):
        as_sample("x", [1, 10**400], finite=False)


def test_as_sample_rejects_an_infinite_value_if_it_must_be_finite():
    with pytest.raises(InputError, match=r"^x\[1\] inf is infinite; expected a finite number$"):
        as_sample("x", [1.0, np.inf], finite=True)
    with pytest.raises(InputError, match=r"^x\[0\] -inf is infinite; expected a finite number$"):
        as_sample("x", [-np.inf, 1.0], finite=True)


@pytest.mark.parametrize(
    ("values", "dimensions"),
    [(np.zeros((2, 2)), 2), (5.0, 0), ([[1.0, 2.0], [3.0, 4.0]], 2)],
)
def test_as_sample_rejects_an_array_that_is_not_1d(values, dimensions):
    message = f"^x has {dimensions} dimensions; expected a 1-D array$"
    with pytest.raises(InputError, match=message):
        as_sample("x", values, finite=False)
