import os
import re
import shutil
import struct
from functools import cache

import numpy as np
import pytest
from annotation_cases import ANNOTATION
from helpers import golden_rows
from hypothesis import given
from hypothesis import strategies as st

from seedmatch.annotation.twobit import TwoBitGenome
from seedmatch.errors import InputError
from seedmatch.sequence import read_fasta
from seedmatch.sequence.alphabet import encode_dna

TWOBIT = ANNOTATION / "synthetic.2bit"
FILES = ["synthetic.2bit", "synthetic.v1.2bit"]
LENGTHS = {"1": 150002, "X": 5003, "SCAF1.1": 2001}
HEADER_SIZE = 16


def index_entries(data):
    """List (name, position of the offset, offset) for each index entry of a little-endian file."""
    (version,) = struct.unpack_from("<I", data, 4)
    (count,) = struct.unpack_from("<I", data, 8)
    offset_format = "<I" if version == 0 else "<Q"
    entries = []
    position = HEADER_SIZE
    for _ in range(count):
        size = data[position]
        name = data[position + 1 : position + 1 + size].decode()
        position += 1 + size
        (offset,) = struct.unpack_from(offset_format, data, position)
        entries.append((name, position, offset))
        position += struct.calcsize(offset_format)
    return entries


def to_big_endian(data):
    """Write each number of the header, the index, and the record headers in big-endian order."""
    out = bytearray(data)

    def swap(position, words, width=4):
        for start in range(position, position + width * words, width):
            out[start : start + width] = out[start : start + width][::-1]

    (version,) = struct.unpack_from("<I", data, 4)
    swap(0, 4)
    for _name, position, offset in index_entries(data):
        swap(position, 1, width=4 if version == 0 else 8)
        (n_blocks,) = struct.unpack_from("<I", data, offset + 4)
        mask_position = offset + 8 + 8 * n_blocks
        (mask_blocks,) = struct.unpack_from("<I", data, mask_position)
        swap(offset, 2 + 2 * n_blocks)
        swap(mask_position, 2 + 2 * mask_blocks)
    return bytes(out)


@cache
def fasta_sequences():
    seqs = read_fasta(ANNOTATION / "synthetic.fa")
    return {name: seqs.sequence(i) for i, name in enumerate(seqs.ids)}


@st.composite
def ranges(draw):
    contig = draw(st.sampled_from(sorted(LENGTHS)))
    start = draw(st.integers(1, LENGTHS[contig] + 1))
    end = draw(st.integers(start - 1, LENGTHS[contig]))
    return contig, start, end


def write(tmp_path, data, name="genome.2bit"):
    path = tmp_path / name
    path.write_bytes(data)
    return path


@pytest.mark.parametrize("name", FILES)
@pytest.mark.parametrize("row", golden_rows("annotation/genome_ranges.tsv"), ids=lambda r: r.label)
def test_twobit_ranges_match_r(name, row):
    with TwoBitGenome(ANNOTATION / name) as genome:
        assert genome.sequence(row.contig, int(row.start), int(row.end)) == row.twobit


@pytest.mark.parametrize("name", FILES)
def test_contigs_give_the_lengths(name):
    with TwoBitGenome(ANNOTATION / name) as genome:
        assert dict(genome.contigs) == LENGTHS


def test_contigs_keep_the_order_of_the_index():
    with TwoBitGenome(TWOBIT) as genome:
        assert list(genome.contigs) == ["1", "X", "SCAF1.1"]


def test_contigs_cannot_be_changed():
    with TwoBitGenome(TWOBIT) as genome, pytest.raises(TypeError):
        genome.contigs["Y"] = 1


def test_path_is_the_path_of_the_file_as_a_path():
    with TwoBitGenome(str(TWOBIT)) as genome:
        assert genome.path == TWOBIT


@pytest.mark.parametrize("name", FILES)
def test_whole_contigs_match_the_fasta_file(name):
    with TwoBitGenome(ANNOTATION / name) as genome:
        for contig, length in LENGTHS.items():
            assert genome.sequence(contig, 1, length) == fasta_sequences()[contig]


@pytest.mark.parametrize("name", FILES)
@given(ranges())
def test_any_range_equals_the_slice_of_the_contig(name, query):
    contig, start, end = query
    with TwoBitGenome(ANNOTATION / name) as genome:
        assert genome.sequence(contig, start, end) == fasta_sequences()[contig][start - 1 : end]


def n_runs(contig):
    return [(m.start() + 1, m.end()) for m in re.finditer("N+", fasta_sequences()[contig])]


@pytest.mark.parametrize(("contig", "run"), [(c, r) for c in LENGTHS for r in n_runs(c)])
def test_ranges_that_start_or_end_at_an_n_block_edge_equal_the_slice(contig, run):
    first, last = run
    with TwoBitGenome(TWOBIT) as genome:
        for start in (first - 3, first - 1, first, first + 1, last, last + 1):
            for end in (first - 1, first, last - 1, last, last + 1, last + 3):
                if 1 <= start <= end + 1 and end <= LENGTHS[contig]:
                    expected = fasta_sequences()[contig][start - 1 : end]
                    assert genome.sequence(contig, start, end) == expected


def test_n_runs_of_the_fixture_include_one_that_starts_at_the_first_base():
    assert n_runs("1")[0][0] == 1
    assert len(n_runs("1")) >= 2


def test_codes_are_the_uint8_masks_of_the_bases():
    with TwoBitGenome(TWOBIT) as genome:
        codes = genome.codes("1", 22398, 22442)
    assert codes.dtype == np.uint8
    assert np.array_equal(codes, encode_dna(fasta_sequences()["1"][22397:22442]))


def test_codes_can_be_changed_by_the_caller():
    with TwoBitGenome(TWOBIT) as genome:
        codes = genome.codes("X", 1, 8)
        codes[:] = 0
        assert genome.sequence("X", 1, 8) == "CGTCCCGG"


def test_n_block_gives_mask_15_in_the_codes():
    with TwoBitGenome(TWOBIT) as genome:
        codes = genome.codes("1", 98, 103)
    assert codes[:3].tolist() == [15, 15, 15]
    assert 15 not in codes[3:].tolist()


@pytest.mark.parametrize("name", FILES)
def test_big_endian_file_gives_the_same_sequences(tmp_path, name):
    data = (ANNOTATION / name).read_bytes()
    big = to_big_endian(data)
    assert big[:4] == bytes.fromhex("1a412743")
    assert big != data
    path = write(tmp_path, big)
    with TwoBitGenome(path) as genome:
        assert dict(genome.contigs) == LENGTHS
        for contig, length in LENGTHS.items():
            assert genome.sequence(contig, 1, length) == fasta_sequences()[contig]
        for row in golden_rows("annotation/genome_ranges.tsv"):
            got = genome.sequence(row.contig, int(row.start), int(row.end))
            assert got == row.twobit, row.label


def test_empty_range_gives_an_empty_array():
    with TwoBitGenome(TWOBIT) as genome:
        codes = genome.codes("X", 10, 9)
        assert codes.dtype == np.uint8
        assert codes.shape == (0,)
        assert genome.sequence("X", 5004, 5003) == ""


def test_range_after_the_contig_end_raises():
    with TwoBitGenome(TWOBIT) as genome, pytest.raises(InputError, match=r"5004.*length 5003"):
        genome.codes("X", 5003, 5004)


def test_unknown_contig_raises_and_names_the_contigs():
    with (
        TwoBitGenome(TWOBIT) as genome,
        pytest.raises(InputError, match=r"'Y'.*1, SCAF1.1, X"),
    ):
        genome.codes("Y", 1, 5)


def test_file_that_is_not_2bit_raises(tmp_path):
    path = write(tmp_path, b"This is a text file, not a genome.\n")
    message = "not a 2bit file; the first bytes are 54686973, expected 4327411a or 1a412743"
    with pytest.raises(InputError, match=re.escape(message)):
        TwoBitGenome(path)


@pytest.mark.parametrize("data", [b"", b"C'A"])
def test_file_that_is_too_short_for_the_signature_is_not_2bit(tmp_path, data):
    with pytest.raises(InputError, match="not a 2bit file"):
        TwoBitGenome(write(tmp_path, data))


@pytest.mark.parametrize("order", ["<", ">"])
def test_unknown_version_raises(tmp_path, order):
    data = bytearray(TWOBIT.read_bytes() if order == "<" else to_big_endian(TWOBIT.read_bytes()))
    data[4:8] = struct.pack(f"{order}I", 2)
    with pytest.raises(InputError, match=r"version 2; expected 0 or 1"):
        TwoBitGenome(write(tmp_path, bytes(data)))


@pytest.mark.parametrize(
    ("keep", "what"),
    [
        (10, "the header"),
        (20, "the index"),
        (0x2A, "the record of contig '1'"),
    ],
)
def test_truncated_file_raises(tmp_path, keep, what):
    path = write(tmp_path, TWOBIT.read_bytes()[:keep])
    with pytest.raises(InputError, match=rf"ends before {what}"):
        TwoBitGenome(path)


def test_file_that_ends_inside_the_blocks_of_a_record_raises_at_the_request_of_the_contig(tmp_path):
    data = TWOBIT.read_bytes()
    offset = {name: offset for name, _, offset in index_entries(data)}["SCAF1.1"]
    path = write(tmp_path, data[: offset + 6])
    with TwoBitGenome(path) as genome:
        assert genome.sequence("X", 1, 4) == "CGTC"
        with pytest.raises(InputError, match=r"ends before the record of contig 'SCAF1.1'"):
            genome.codes("SCAF1.1", 1, 4)


def test_file_that_ends_inside_the_dna_raises_at_the_request_of_the_contig_for_each_range(tmp_path):
    path = write(tmp_path, TWOBIT.read_bytes()[:-1])
    with TwoBitGenome(path) as genome:
        assert genome.sequence("1", 1, 4) == "NNNN"
        for start, end in [(1, 4), (1998, 2001), (1, 2001)]:
            with pytest.raises(InputError, match=r"ends before the DNA of contig 'SCAF1.1'"):
                genome.codes("SCAF1.1", start, end)


def test_contig_name_that_is_not_utf8_raises(tmp_path):
    data = bytearray(TWOBIT.read_bytes())
    (_, position, _), _x, _scaffold = index_entries(bytes(data))
    data[position - 1] = 0xFF
    message = r"contig name b'\\xff' is not valid UTF-8"
    with pytest.raises(InputError, match=message):
        TwoBitGenome(write(tmp_path, bytes(data)))


@pytest.mark.parametrize("name", FILES)
def test_duplicate_contig_name_raises(tmp_path, name):
    data = bytearray((ANNOTATION / name).read_bytes())
    _first, (_, position, _), _scaffold = index_entries(bytes(data))
    data[position - 1] = ord("1")
    path = write(tmp_path, bytes(data))
    message = f"{os.fspath(path)!r}: contig '1' occurs twice in the index; expected each name once"
    with pytest.raises(InputError, match=re.escape(message)):
        TwoBitGenome(path)


def test_missing_file_raises_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        TwoBitGenome(tmp_path / "missing.2bit")


def test_url_raises():
    with pytest.raises(InputError, match="the path is a URL; expected a local path"):
        TwoBitGenome("https://example.org/genome.2bit")


def test_reads_from_a_read_only_directory(tmp_path):
    directory = tmp_path / "refs"
    directory.mkdir()
    file = directory / "synthetic.2bit"
    shutil.copy(TWOBIT, file)
    file.chmod(0o444)
    directory.chmod(0o555)
    try:
        with TwoBitGenome(file) as genome:
            assert genome.sequence("X", 1, 5003) == fasta_sequences()["X"]
        assert [entry.name for entry in directory.iterdir()] == ["synthetic.2bit"]
    finally:
        directory.chmod(0o755)


def test_close_releases_the_file_and_a_later_request_raises():
    genome = TwoBitGenome(TWOBIT)
    assert genome.sequence("X", 1, 4) == "CGTC"
    genome.close()
    genome.close()
    with pytest.raises(ValueError, match="mmap closed"):
        genome.codes("X", 1, 4)


def test_with_block_closes_the_file():
    with TwoBitGenome(TWOBIT) as genome:
        pass
    with pytest.raises(ValueError, match="mmap closed"):
        genome.codes("X", 1, 4)


def test_resolve_maps_the_contigs_of_the_file():
    with TwoBitGenome(TWOBIT) as genome:
        assert genome.resolve("chr1") == "1"
        assert genome.resolve("chrX") == "X"
        assert genome.resolve("SCAF1.1") == "SCAF1.1"
        assert genome.resolve("Z") is None
