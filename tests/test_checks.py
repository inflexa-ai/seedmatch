import math
import re
import sys
import warnings
from enum import StrEnum

import numpy as np
import pandas as pd
import pytest

from seedmatch._checks import (
    as_cutoff,
    as_flag,
    as_items,
    as_limit,
    as_member,
    as_number,
    as_text,
    check_at_least,
    check_not_negative,
    check_numeric_columns,
)
from seedmatch.errors import InputError
from seedmatch.report import FullSearchOptions

# Python refuses the text form of an int with more digits than the limit, so repr raises ValueError.
HUGE = 10**5000


@pytest.fixture(autouse=True)
def digit_limit():
    previous = sys.get_int_max_str_digits()
    sys.set_int_max_str_digits(4300)
    yield
    sys.set_int_max_str_digits(previous)


def test_the_digit_limit_makes_the_text_form_of_a_huge_int_fail():
    with pytest.raises(ValueError, match="digits"):
        repr(HUGE)


def test_as_limit_gives_input_error_for_an_int_above_the_range_of_the_kernels():
    with pytest.raises(InputError, match=r"distance .*\d+ bits.* is too large; expected at most"):
        as_limit("distance", HUGE)


def test_a_huge_negative_int_is_a_limit_that_check_not_negative_rejects():
    limit = as_limit("distance", -HUGE)
    with pytest.raises(InputError, match=r"distance .*\d+ bits.* is negative; expected 0 or more"):
        check_not_negative("distance", limit)


def test_check_at_least_rejects_a_value_below_the_minimum_and_accepts_the_minimum():
    check_at_least("n_boot", 1, 1)
    check_at_least("n_interp", 2, 2)
    check_at_least("n_interp", 3, 2)
    for value in (0, -3):
        with pytest.raises(
            InputError, match=rf"^n_boot {value} is less than 1; expected 1 or more$"
        ):
            check_at_least("n_boot", value, 1)
    for value in (1, 0, -3):
        with pytest.raises(
            InputError, match=rf"^n_interp {value} is less than 2; expected 2 or more$"
        ):
            check_at_least("n_interp", value, 2)


def test_a_huge_negative_int_is_a_limit_that_check_at_least_rejects():
    limit = as_limit("n_boot", -HUGE)
    with pytest.raises(InputError, match=r"n_boot .*\d+ bits.* is less than 1; expected 1 or more"):
        check_at_least("n_boot", limit, 1)


def test_as_flag_gives_input_error_for_a_huge_int():
    with pytest.raises(InputError, match=r"indels .*\d+ bits.* has type int; expected a bool"):
        as_flag("indels", HUGE)


@pytest.mark.parametrize("sign", [1, -1])
def test_as_cutoff_gives_input_error_for_an_int_that_a_float_cannot_hold(sign):
    with pytest.raises(InputError, match=r"padj_cut .*\d+ bits.* is too large; expected a number"):
        as_cutoff("padj_cut", sign * HUGE)


def test_as_number_returns_a_float_and_keeps_nan_and_infinity():
    assert as_number("x", 3) == 3.0
    assert isinstance(as_number("x", 3), float)
    assert as_number("x", float("-inf")) == -math.inf
    assert math.isnan(as_number("x", math.nan))
    with pytest.raises(InputError, match=r"x True has type bool; expected an int or a float"):
        as_number("x", value=True)
    with pytest.raises(InputError, match=r"x '1' has type str; expected an int or a float"):
        as_number("x", "1")


def test_a_value_with_a_text_form_is_named_by_its_repr():
    with pytest.raises(InputError, match=r"distance 9223372036854775808 is too large"):
        as_limit("distance", 2**63)
    with pytest.raises(InputError, match=r"indels 'yes' has type str"):
        as_flag("indels", "yes")
    with pytest.raises(InputError, match=r"padj_cut '1' has type str"):
        as_cutoff("padj_cut", "1")


def test_the_options_name_a_huge_int_by_its_size():
    with pytest.raises(InputError, match=r"seeds .*\d+ bits.* has type int; expected a tuple"):
        FullSearchOptions(seeds=HUGE)
    with pytest.raises(InputError, match=r"seeds holds .*\d+ bits.* of type int; expected a Seed"):
        FullSearchOptions(seeds=(HUGE,))
    with pytest.raises(InputError, match=r"seed_edit_distances .*\d+ bits.* is too large"):
        FullSearchOptions(seed_edit_distances=(HUGE,))


class Colour(StrEnum):
    RED = "red"
    DARK_BLUE = "dark.blue"


def test_as_member_returns_a_member_of_the_enum():
    assert as_member("colour", Colour.DARK_BLUE, Colour) is Colour.DARK_BLUE


@pytest.mark.parametrize("value", ["red", "dark.blue", "Red", None, 1])
def test_as_member_rejects_the_text_of_a_member_and_other_values(value):
    message = f"colour {value!r} is not valid; expected one of Colour.RED, Colour.DARK_BLUE"
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        as_member("colour", value, Colour)


def test_as_text_returns_a_str_and_names_another_type():
    assert as_text("symbols[0]", "Ttr") == "Ttr"
    with pytest.raises(InputError, match=r"^symbols\[1\] 5 has type int; expected a str$"):
        as_text("symbols[1]", 5)


def test_as_items_returns_a_tuple_of_the_checked_items_in_order():
    assert as_items("ids", iter(["b", "a", "b"]), as_text, "str") == ("b", "a", "b")
    assert as_items("ids", [], as_text, "str") == ()


def test_as_items_names_a_bare_str_with_an_example():
    message = "ids 'g1' is a str; expected a sequence of str, for example ('g1',)"
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        as_items("ids", "g1", as_text, "str")
    message = "ids 'g1' is a str; expected an iterable of gene IDs, for example ('g7',)"
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        as_items("ids", "g1", as_text, "gene IDs", container="an iterable", example="'g7'")


def test_as_items_names_a_value_that_is_not_iterable_and_the_index_of_a_bad_item():
    with pytest.raises(InputError, match=r"^ids 5 has type int; expected a sequence$"):
        as_items("ids", 5, as_text, "str")
    with pytest.raises(InputError, match=r"^ids 5 has type int; expected an iterable$"):
        as_items("ids", 5, as_text, "str", container="an iterable")
    with pytest.raises(InputError, match=r"^ids\[2\] None has type NoneType; expected a str$"):
        as_items("ids", ["a", "b", None], as_text, "str")


@pytest.mark.parametrize(
    ("column", "dtype"),
    [
        (np.array([True, False]), "bool"),
        (np.array([1 + 2j, 3 + 0j]), "complex128"),
        (pd.array([True, None], dtype="boolean"), "boolean"),
        (["a", "b"], "object"),
    ],
)
def test_check_numeric_columns_rejects_a_bool_column_a_complex_column_and_text(column, dtype):
    de = pd.DataFrame({"value": column})
    message = rf"^column 'value' has dtype {dtype}; expected a numeric dtype$"
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(InputError, match=message):
            check_numeric_columns(de, ["value"])


@pytest.mark.parametrize("dtype", ["int64", "float64", "uint8", "Int64", "Float64"])
def test_check_numeric_columns_accepts_an_integer_and_a_float_column(dtype):
    check_numeric_columns(pd.DataFrame({"value": pd.array([1, 2], dtype=dtype)}), ["value"])
