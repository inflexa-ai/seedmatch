import math
import sys

import pytest

from seedmatch._checks import (
    as_cutoff,
    as_flag,
    as_limit,
    as_number,
    check_not_negative,
    check_positive,
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


def test_check_positive_rejects_0_and_a_negative_int_and_accepts_1():
    check_positive("n_boot", 1)
    for value in (0, -3):
        with pytest.raises(
            InputError, match=rf"^n_boot {value} is less than 1; expected 1 or more$"
        ):
            check_positive("n_boot", value)


def test_a_huge_negative_int_is_a_limit_that_check_positive_rejects():
    limit = as_limit("n_boot", -HUGE)
    with pytest.raises(InputError, match=r"n_boot .*\d+ bits.* is less than 1; expected 1 or more"):
        check_positive("n_boot", limit)


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
