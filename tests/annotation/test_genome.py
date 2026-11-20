import re
from pathlib import Path

import numpy as np
import pytest

from seedmatch.annotation.genome import Genome
from seedmatch.errors import InputError


class Names(Genome):
    """A genome of contigs with length 10. Its bases are all N. It records each call of `_codes`."""

    def __init__(self, names, path="names.2bit"):
        super().__init__(path)
        self._lengths = dict.fromkeys(names, 10)
        self.closed = 0
        self.calls = []

    @property
    def contigs(self):
        return self._lengths

    def _codes(self, contig, start, end):
        self.calls.append((contig, start, end))
        return np.full(end - start + 1, 15, dtype=np.uint8)

    def close(self):
        self.closed += 1


@pytest.mark.parametrize(
    ("names", "name", "expected"),
    [
        (["1", "X"], "1", "1"),
        (["1", "X"], "chr1", "1"),
        (["1", "X"], "chrX", "X"),
        (["1", "X"], "Z", None),
        (["1", "X"], "chrZ", None),
        (["chr1", "chrX"], "1", "chr1"),
        (["chr1", "chrX"], "chr1", "chr1"),
        (["chr1", "chrX"], "X", "chrX"),
        (["1", "chr1"], "1", "1"),
        (["1", "chr1"], "chr1", "chr1"),
        (["MT"], "MT", "MT"),
        (["MT"], "chrM", "MT"),
        (["chrM"], "chrM", "chrM"),
        (["chrM"], "MT", "chrM"),
        (["chrM", "MT"], "MT", "MT"),
        (["chrM", "MT"], "chrM", "chrM"),
        (["M"], "chrM", "M"),
        (["1", "X"], "MT", None),
        (["1", "X"], "chrM", None),
    ],
)
def test_resolve_maps_chr_names_and_mt(names, name, expected):
    assert Names(names).resolve(name) == expected


def test_path_is_a_path():
    assert Names(["1"]).path == Path("names.2bit")


def test_genome_refuses_a_url_as_the_path():
    message = "'https://example.org/genome.2bit': the path is a URL; expected a local path"
    with pytest.raises(InputError, match=re.escape(message)):
        Names(["1"], "https://example.org/genome.2bit")


def test_sequence_decodes_the_codes():
    assert Names(["1"]).sequence("1", 2, 5) == "NNNN"


def test_with_block_closes_the_genome_when_it_exits():
    genome = Names(["1"])
    with genome as entered:
        assert entered is genome
        assert genome.closed == 0
    assert genome.closed == 1


def test_with_block_closes_the_genome_when_the_body_raises():
    genome = Names(["1"])
    message = "body"
    with pytest.raises(RuntimeError, match=message), genome:
        raise RuntimeError(message)
    assert genome.closed == 1


def test_unknown_contig_lists_the_contigs_in_sorted_order():
    message = "contig 'Y' is not in the genome; expected one of 1, 2, X"
    with pytest.raises(InputError, match=re.escape(message)):
        Names(["X", "2", "1"]).codes("Y", 1, 1)


def test_unknown_contig_lists_only_the_first_ten_contigs():
    names = [f"c{number:02d}" for number in range(12, 0, -1)]
    listed = ", ".join(f"c{number:02d}" for number in range(1, 11))
    message = f"contig 'Y' is not in the genome; expected one of {listed}, ..."
    with pytest.raises(InputError, match=re.escape(message)) as raised:
        Names(names).codes("Y", 1, 1)
    assert "c11" not in str(raised.value)


def test_unknown_contig_with_ten_contigs_lists_all_of_them_without_an_ellipsis():
    names = [f"c{number:02d}" for number in range(1, 11)]
    with pytest.raises(InputError) as raised:
        Names(names).codes("Y", 1, 1)
    assert str(raised.value).endswith("c09, c10")


def test_unknown_contig_raises_for_an_empty_range_too():
    with pytest.raises(InputError, match=r"contig 'Y' is not in the genome"):
        Names(["1"]).codes("Y", 1, 0)


@pytest.mark.parametrize(("start", "end"), [(1, 10), (10, 10), (1, 1), (4, 7)])
def test_range_inside_the_contig_is_accepted(start, end):
    assert len(Names(["1"]).codes("1", start, end)) == end - start + 1


@pytest.mark.parametrize(("start", "end"), [(1, 0), (5, 4), (10, 9), (11, 10)])
def test_empty_range_inside_the_boundaries_gives_an_empty_array(start, end):
    codes = Names(["1"]).codes("1", start, end)
    assert codes.dtype == np.uint8
    assert codes.shape == (0,)
    assert Names(["1"]).sequence("1", start, end) == ""


@pytest.mark.parametrize(
    ("start", "end", "shown"),
    [(0, 5, "0..5"), (-3, 5, "-3..5"), (5, 11, "5..11"), (11, 11, "11..11"), (1, 100, "1..100")],
)
def test_range_outside_the_contig_names_the_range_and_the_length(start, end, shown):
    message = (
        f"contig '1': range {shown} is outside the contig of length 10; "
        "expected 1 <= start and end <= 10"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        Names(["1"]).codes("1", start, end)


@pytest.mark.parametrize(("start", "end"), [(0, -1), (12, 11), (-5, -6)])
def test_empty_range_outside_the_boundaries_raises(start, end):
    with pytest.raises(InputError, match=r"is outside the contig of length 10"):
        Names(["1"]).codes("1", start, end)


@pytest.mark.parametrize(("start", "end"), [(6, 4), (10, 1)])
def test_range_that_ends_before_it_starts_raises(start, end):
    message = (
        f"contig '1': range {start}..{end} ends before it starts; "
        "expected end >= start - 1 (end == start - 1 is an empty range)"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        Names(["1"]).codes("1", start, end)


@pytest.mark.parametrize(
    ("start", "end", "name", "value", "kind"),
    [
        ("1", 5, "start", "'1'", "str"),
        (1, 5.0, "end", "5.0", "float"),
        (True, 5, "start", "True", "bool"),
        (1, None, "end", "None", "NoneType"),
    ],
)
def test_range_bound_that_is_not_an_integer_raises(start, end, name, value, kind):
    message = f"{name} {value} has type {kind}; expected an integer"
    with pytest.raises(InputError, match=re.escape(message)):
        Names(["1"]).codes("1", start, end)


def test_range_bounds_can_be_numpy_integers():
    assert Names(["1"]).sequence("1", np.int64(2), np.uint16(4)) == "NNN"


def test_reader_gets_only_checked_ranges_with_plain_ints():
    genome = Names(["1"])
    genome.codes("1", np.int64(2), np.uint16(4))
    assert genome.calls == [("1", 2, 4)]
    assert {type(value) for value in genome.calls[0][1:]} == {int}


@pytest.mark.parametrize(
    ("contig", "start", "end", "message"),
    [
        ("Y", 1, 5, "contig 'Y' is not in the genome"),
        ("Y", 1, 0, "contig 'Y' is not in the genome"),
        ("1", 0, 5, r"range 0\.\.5 is outside the contig"),
        ("1", 5, 11, r"range 5\.\.11 is outside the contig"),
        ("1", 6, 4, r"range 6\.\.4 ends before it starts"),
        ("1", "1", 5, "start '1' has type str; expected an integer"),
    ],
)
def test_reader_is_not_called_for_an_unknown_contig_or_a_bad_range(contig, start, end, message):
    genome = Names(["1"])
    with pytest.raises(InputError, match=message):
        genome.codes(contig, start, end)
    assert genome.calls == []


def test_reader_is_not_called_for_an_empty_range():
    genome = Names(["1"])
    assert genome.codes("1", 4, 3).shape == (0,)
    assert genome.calls == []
