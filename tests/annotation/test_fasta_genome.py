import gzip
import os
import re
import shutil
import struct
import tracemalloc
import zlib
from functools import cache

import numpy as np
import pytest
from annotation_cases import ANNOTATION
from helpers import golden_rows
from hypothesis import given, settings
from hypothesis import strategies as st

from seedmatch.annotation.fasta_genome import FastaGenome
from seedmatch.annotation.genome_file import open_genome
from seedmatch.annotation.twobit import TwoBitGenome
from seedmatch.errors import InputError
from seedmatch.sequence import read_fasta
from seedmatch.sequence.alphabet import encode_dna

NAMES = ["synthetic.fa", "synthetic.fa.gz"]
PLAIN = ANNOTATION / "synthetic.fa"
BGZF = ANNOTATION / "synthetic.fa.gz"
LENGTHS = {"1": 150002, "X": 5003, "SCAF1.1": 2001}
LAYOUTS = {
    "lf": ("\n", None),
    "crlf": ("\r\n", None),
    "lf_bgzf": ("\n", 211),
    "crlf_bgzf": ("\r\n", 1000),
}
BASES_PER_LINE = 60


@cache
def fasta_sequences():
    seqs = read_fasta(PLAIN)
    return {name: seqs.sequence(i) for i, name in enumerate(seqs.ids)}


def fasta_bytes(sequences, ending):
    """Return the FASTA bytes of the sequences and the text of their .fai."""
    body = bytearray()
    fai = []
    for name, sequence in sequences.items():
        body += f">{name}{ending}".encode()
        offset = len(body)
        for start in range(0, len(sequence), BASES_PER_LINE):
            body += sequence[start : start + BASES_PER_LINE].encode() + ending.encode()
        fai.append(
            f"{name}\t{len(sequence)}\t{offset}\t{BASES_PER_LINE}\t{BASES_PER_LINE + len(ending)}\n"
        )
    return bytes(body), "".join(fai)


def bgzf_block(data, subfields=b""):
    compressor = zlib.compressobj(wbits=-15)
    body = compressor.compress(data) + compressor.flush()
    total = 12 + len(subfields) + 6 + len(body) + 8
    extra = subfields + struct.pack("<2sHH", b"BC", 2, total - 1)
    header = struct.pack("<2sBBIBBH", b"\x1f\x8b", 8, 4, 0, 0, 255, len(extra))
    return header + extra + body + struct.pack("<II", zlib.crc32(data), len(data))


def bgzf_bytes(data, block_size, subfields=b""):
    """Return the BGZF bytes of the data in blocks of `block_size` bytes, and its .gzi bytes."""
    out = bytearray()
    pairs = []
    for start in range(0, len(data), block_size):
        if start:
            pairs.append((len(out), start))
        out += bgzf_block(data[start : start + block_size], subfields)
    out += bgzf_block(b"", subfields)
    gzi = struct.pack("<Q", len(pairs)) + b"".join(struct.pack("<QQ", *pair) for pair in pairs)
    return bytes(out), gzi


def write_genome(directory, sequences, ending="\n", block_size=None, subfields=b""):
    """Write `genome.fa` with its index files. A block size gives BGZF. Return the path."""
    data, fai = fasta_bytes(sequences, ending)
    name = "genome.fa"
    path = directory / name
    if block_size is None:
        path.write_bytes(data)
    else:
        packed, gzi = bgzf_bytes(data, block_size, subfields)
        path.write_bytes(packed)
        (directory / f"{name}.gzi").write_bytes(gzi)
    (directory / f"{name}.fai").write_text(fai)
    return path


def copy_as(directory, name, new_name):
    """Copy a fixture file and its index files. Return the path of the copy."""
    for suffix in ("", ".fai", ".gzi"):
        source = ANNOTATION / f"{name}{suffix}"
        if source.exists():
            shutil.copy(source, directory / f"{new_name}{suffix}")
    return directory / new_name


def read_gzi(path):
    data = path.read_bytes()
    (count,) = struct.unpack_from("<Q", data)
    values = struct.unpack_from(f"<{2 * count}Q", data, 8)
    return list(zip(values[0::2], values[1::2], strict=True))


@pytest.fixture(scope="module")
def layouts(tmp_path_factory):
    directory = tmp_path_factory.mktemp("layouts")
    paths = {}
    for layout, (ending, block_size) in LAYOUTS.items():
        (directory / layout).mkdir()
        paths[layout] = write_genome(directory / layout, fasta_sequences(), ending, block_size)
    return paths


@st.composite
def ranges(draw):
    contig = draw(st.sampled_from(sorted(LENGTHS)))
    start = draw(st.integers(1, LENGTHS[contig]))
    end = draw(st.integers(start, LENGTHS[contig]))
    return contig, start, end


def edges(length):
    low = [1, 2, 59, 60, 61, 62, 119, 120, 121]
    return sorted(
        {position for position in [*low, *(length - p + 1 for p in low)] if 1 <= position <= length}
    )


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("row", golden_rows("annotation/genome_ranges.tsv"), ids=lambda r: r.label)
def test_fasta_ranges_match_r(name, row):
    with FastaGenome(ANNOTATION / name) as genome:
        assert genome.sequence(row.contig, int(row.start), int(row.end)) == row.bgzf


@pytest.mark.parametrize("name", NAMES)
def test_contigs_give_the_lengths_in_the_order_of_the_fai(name):
    with FastaGenome(ANNOTATION / name) as genome:
        assert dict(genome.contigs) == LENGTHS
        assert list(genome.contigs) == ["1", "X", "SCAF1.1"]


def test_contigs_cannot_be_changed():
    with FastaGenome(PLAIN) as genome, pytest.raises(TypeError):
        genome.contigs["Y"] = 1


def test_path_is_the_path_of_the_file_as_a_path():
    with FastaGenome(str(PLAIN)) as genome:
        assert genome.path == PLAIN


@pytest.mark.parametrize("name", NAMES)
def test_fasta_and_twobit_give_the_same_contigs(name):
    with (
        FastaGenome(ANNOTATION / name) as fasta,
        TwoBitGenome(ANNOTATION / "synthetic.2bit") as twobit,
    ):
        assert dict(fasta.contigs) == dict(twobit.contigs)
        for contig, length in LENGTHS.items():
            assert fasta.sequence(contig, 1, length) == twobit.sequence(contig, 1, length)


def test_the_helper_writes_the_fai_of_the_fixture():
    assert fasta_bytes(fasta_sequences(), "\n")[1] == (ANNOTATION / "synthetic.fa.fai").read_text()


@pytest.mark.parametrize("layout", LAYOUTS)
def test_whole_contigs_match_the_fasta_file(layouts, layout):
    with FastaGenome(layouts[layout]) as genome:
        assert dict(genome.contigs) == LENGTHS
        for contig, length in LENGTHS.items():
            assert genome.sequence(contig, 1, length) == fasta_sequences()[contig]


@pytest.mark.parametrize("layout", LAYOUTS)
@settings(deadline=None)
@given(ranges())
def test_any_range_equals_the_slice_of_the_contig(layouts, layout, query):
    contig, start, end = query
    with FastaGenome(layouts[layout]) as genome:
        assert genome.sequence(contig, start, end) == fasta_sequences()[contig][start - 1 : end]


@pytest.mark.parametrize("layout", LAYOUTS)
@pytest.mark.parametrize("contig", LENGTHS)
def test_ranges_at_the_line_ends_and_the_contig_ends_equal_the_slice(layouts, layout, contig):
    sequence = fasta_sequences()[contig]
    positions = edges(LENGTHS[contig])
    with FastaGenome(layouts[layout]) as genome:
        for start in positions:
            for end in positions:
                if start <= end < start + 300:
                    assert genome.sequence(contig, start, end) == sequence[start - 1 : end]


@pytest.mark.parametrize("name", NAMES)
def test_codes_are_the_uint8_masks_of_the_bases(name):
    with FastaGenome(ANNOTATION / name) as genome:
        codes = genome.codes("1", 22398, 22442)
    assert codes.dtype == np.uint8
    assert np.array_equal(codes, encode_dna(fasta_sequences()["1"][22397:22442]))


@pytest.mark.parametrize("name", NAMES)
def test_codes_can_be_changed_by_the_caller(name):
    with FastaGenome(ANNOTATION / name) as genome:
        codes = genome.codes("X", 1, 8)
        codes[:] = 0
        assert genome.sequence("X", 1, 8) == "CGTCCCGG"


def test_soft_masked_and_iupac_letters_give_uppercase_masks(tmp_path):
    path = write_genome(tmp_path, {"a": "acgtnRYKMswbdhvN"})
    with FastaGenome(path) as genome:
        assert genome.sequence("a", 1, 16) == "ACGTNRYKMSWBDHVN"


def test_last_line_without_a_line_end_is_read(tmp_path):
    data, fai = fasta_bytes(fasta_sequences(), "\n")
    assert data.endswith(b"\n")
    path = tmp_path / "genome.fa"
    path.write_bytes(data.removesuffix(b"\n"))
    (tmp_path / "genome.fa.fai").write_text(fai)
    with FastaGenome(path) as genome:
        last = LENGTHS["SCAF1.1"]
        assert genome.sequence("SCAF1.1", 1, last) == fasta_sequences()["SCAF1.1"]
        assert genome.sequence("SCAF1.1", last, last) == fasta_sequences()["SCAF1.1"][-1]


@pytest.mark.parametrize("block_size", [None, 20])
def test_letter_that_is_not_dna_raises_and_names_the_file_and_the_contig(tmp_path, block_size):
    path = write_genome(tmp_path, {"a": "ACGT", "b": "AC*TACGT"}, block_size=block_size)
    with FastaGenome(path) as genome:
        assert genome.sequence("a", 1, 4) == "ACGT"
        message = (
            f"{os.fspath(path)!r}: contig 'b' base 3 is '*', not an IUPAC DNA letter; "
            "expected one of ACGTMRWSYKVHDBN in either case"
        )
        with pytest.raises(InputError, match=re.escape(message)):
            genome.codes("b", 1, 8)


def test_letter_that_is_not_dna_names_its_base_in_the_contig(tmp_path):
    path = write_genome(tmp_path, {"a": "ACGTACGT", "b": "ACGTAC*T"})
    message = f"{os.fspath(path)!r}: contig 'b' base 7 is '*', not an IUPAC DNA letter"
    with FastaGenome(path) as genome, pytest.raises(InputError, match=re.escape(message)):
        genome.codes("b", 3, 8)


def test_letter_after_the_first_piece_of_a_plain_file_names_its_base_in_the_contig(tmp_path):
    length = 3 * 2**20
    sequence = "A" * 2097168 + "*" + "A" * (length - 2097169)
    path = write_genome(tmp_path, {"1": sequence})
    message = f"{os.fspath(path)!r}: contig '1' base 2097169 is '*', not an IUPAC DNA letter"
    with FastaGenome(path) as genome, pytest.raises(InputError, match=re.escape(message)):
        genome.codes("1", 1, length)


@pytest.mark.parametrize(("start", "end"), [(1, 400), (200, 300), (250, 250), (240, 399)])
def test_letter_in_a_later_block_names_its_base_in_the_contig(tmp_path, start, end):
    sequence = "ACGT" * 62 + "AC*TA" + "ACGT" * 37
    path = write_genome(tmp_path, {"a": sequence}, block_size=50)
    message = f"{os.fspath(path)!r}: contig 'a' base 251 is '*', not an IUPAC DNA letter"
    with FastaGenome(path) as genome:
        if start <= 251 <= end:
            with pytest.raises(InputError, match=re.escape(message)):
                genome.codes("a", start, end)
        else:
            assert genome.sequence("a", start, end) == sequence[start - 1 : end]


def test_letter_at_the_first_base_names_base_1(tmp_path):
    path = write_genome(tmp_path, {"a": "*ACGT"})
    message = f"{os.fspath(path)!r}: contig 'a' base 1 is '*', not an IUPAC DNA letter"
    with FastaGenome(path) as genome, pytest.raises(InputError, match=re.escape(message)):
        genome.codes("a", 1, 5)


@pytest.mark.parametrize(("start", "end"), [(1, 177), (97, 120), (97, 97)])
def test_letter_at_the_first_byte_of_a_later_block_names_its_base(tmp_path, start, end):
    # With blocks of 50 bytes, the byte of base 97 is the first byte of the third block.
    sequence = "ACGT" * 24 + "*" + "ACGT" * 20
    path = write_genome(tmp_path, {"a": sequence}, block_size=50)
    assert read_gzi(tmp_path / "genome.fa.gzi")[1][1] == 100
    message = f"{os.fspath(path)!r}: contig 'a' base 97 is '*', not an IUPAC DNA letter"
    with FastaGenome(path) as genome, pytest.raises(InputError, match=re.escape(message)):
        genome.codes("a", start, end)


def test_letter_at_the_first_byte_of_the_second_piece_of_a_plain_file_names_its_base(tmp_path):
    read_size = 2**20
    line, column = divmod(read_size, 61)
    assert column < 60
    base = line * 60 + column + 1
    sequence = "A" * (base - 1) + "*" + "A" * 1000
    path = write_genome(tmp_path, {"a": sequence})
    message = f"{os.fspath(path)!r}: contig 'a' base {base} is '*', not an IUPAC DNA letter"
    with FastaGenome(path) as genome, pytest.raises(InputError, match=re.escape(message)):
        genome.codes("a", 1, len(sequence))


def test_byte_that_is_not_ascii_is_a_letter_that_is_not_dna(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_bytes(b">a\nAC\xffT\n")
    (tmp_path / "genome.fa.fai").write_text("a\t4\t3\t4\t5\n")
    message = f"{os.fspath(path)!r}: contig 'a' base 3 is 'ÿ', not an IUPAC DNA letter"
    with FastaGenome(path) as genome, pytest.raises(InputError, match=re.escape(message)):
        genome.codes("a", 1, 4)


def test_bytes_outside_the_range_are_not_read(tmp_path):
    path = write_genome(tmp_path, {"a": "ACGTACGT" + "*" * 8})
    with FastaGenome(path) as genome:
        assert genome.sequence("a", 1, 8) == "ACGTACGT"
        with pytest.raises(InputError, match=r"base 9 is '\*'"):
            genome.codes("a", 1, 9)


def test_file_with_shorter_lines_than_its_fai_says_raises(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_bytes(b">a\nACGTA\nCGTAC\n")
    (tmp_path / "genome.fa.fai").write_text("a\t10\t3\t10\t11\n")
    message = (
        f"{os.fspath(path)!r}: contig 'a': the range 1..10 has 9 bases in the file; "
        "expected 10 bases, as the .fai says"
    )
    with FastaGenome(path) as genome:
        assert genome.sequence("a", 1, 5) == "ACGTA"
        with pytest.raises(InputError, match=re.escape(message)):
            genome.codes("a", 1, 10)


def test_file_that_becomes_shorter_after_open_raises(tmp_path):
    path = write_genome(tmp_path, {"a": "ACGT" * 5000})
    with FastaGenome(path) as genome:
        os.truncate(path, 12000)
        message = r"the range 15000\.\.15010 has 0 bases in the file; expected 11 bases"
        with pytest.raises(InputError, match=message):
            genome.codes("a", 15000, 15010)
        with pytest.raises(InputError, match=r"the range 11000\.\.16000 has \d+ bases in the file"):
            genome.codes("a", 11000, 16000)


def test_message_of_a_range_of_one_base_that_the_file_lacks_has_the_singular(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_bytes(b">a\nAC\nGT\n")
    (tmp_path / "genome.fa.fai").write_text("a\t4\t3\t4\t5\n")
    message = (
        f"{os.fspath(path)!r}: contig 'a': the range 3..3 has 0 bases in the file; "
        "expected 1 base, as the .fai says"
    )
    with FastaGenome(path) as genome, pytest.raises(InputError, match=re.escape(message)):
        genome.codes("a", 3, 3)


def test_contig_that_ends_at_the_last_byte_of_a_plain_file_is_valid(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_bytes(b">a\nACGT")
    (tmp_path / "genome.fa.fai").write_text("a\t4\t3\t10\t11\n")
    with FastaGenome(path) as genome:
        assert genome.sequence("a", 1, 4) == "ACGT"


def test_contig_that_ends_after_the_last_byte_of_a_plain_file_raises(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_bytes(b">a\nACGT")
    (tmp_path / "genome.fa.fai").write_text("a\t5\t3\t10\t11\n")
    message = (
        f"{os.fspath(path)}.fai': contig 'a' has its last base at byte 7, but the file holds at "
        "most 7 bytes of data; the .fai does not match the file. "
        f"Make it with `samtools faidx {os.fspath(path)}`"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        FastaGenome(path)


@pytest.mark.parametrize("name", NAMES)
def test_fai_of_a_longer_genome_raises_at_open_and_names_the_contig(tmp_path, name):
    path = copy_as(tmp_path, name, name)
    (tmp_path / f"{name}.fai").write_text("1\t150002\t3\t60\t61\nX\t5003\t900000\t60\t61\n")
    message = f"{os.fspath(path)}.fai': contig 'X' has its last base at byte 905085"
    with pytest.raises(
        InputError, match=re.escape(message) + r".*the \.fai does not match the file"
    ):
        FastaGenome(path)


@pytest.mark.parametrize("name", NAMES)
def test_fai_that_claims_a_contig_of_a_billion_bases_raises_and_not_memory_error(tmp_path, name):
    path = copy_as(tmp_path, name, name)
    (tmp_path / f"{name}.fai").write_text("1\t1000000000\t3\t60\t61\n")
    with pytest.raises(InputError, match=r"contig '1' has its last base at byte"):
        FastaGenome(path)


@pytest.mark.parametrize(("length", "valid"), [(196076, True), (196077, False)])
def test_bgzf_contig_must_end_before_the_last_block_start_plus_65536(tmp_path, length, valid):
    path = copy_as(tmp_path, "synthetic.fa.gz", "synthetic.fa.gz")
    (tmp_path / "synthetic.fa.gz.fai").write_text(f"1\t{length}\t3\t1000000\t1000001\n")
    if valid:
        with FastaGenome(path) as genome:
            assert genome.sequence("1", 1, 8) == fasta_sequences()["1"][:8]
    else:
        message = r"last base at byte 196079, but the file holds at most 196079 bytes of data"
        with pytest.raises(InputError, match=message):
            FastaGenome(path)


def test_contig_of_length_zero_needs_no_data(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_bytes(b"")
    (tmp_path / "genome.fa.fai").write_text("a\t0\t3\t60\t61\n")
    with FastaGenome(path) as genome:
        assert dict(genome.contigs) == {"a": 0}


def test_file_whose_lines_are_longer_than_its_fai_says_raises(tmp_path):
    path = write_genome(tmp_path, {"a": "ACGTACGTAC"})
    path.write_bytes(b">a\nACGTACGTACGTACGTACGT\n")
    (tmp_path / "genome.fa.fai").write_text("a\t10\t3\t5\t6\n")
    with (
        FastaGenome(path) as genome,
        pytest.raises(InputError, match=r"has more than 10 bases.*expected 10 bases"),
    ):
        genome.codes("a", 1, 10)


def test_missing_fai_raises_and_names_the_command(tmp_path):
    shutil.copy(PLAIN, tmp_path / "synthetic.fa")
    with pytest.raises(
        FileNotFoundError, match=r"synthetic\.fa\.fai.*samtools faidx .*synthetic\.fa"
    ):
        FastaGenome(tmp_path / "synthetic.fa")


def test_missing_gzi_raises_and_names_the_command(tmp_path):
    path = copy_as(tmp_path, "synthetic.fa.gz", "synthetic.fa.gz")
    (tmp_path / "synthetic.fa.gz.gzi").unlink()
    with pytest.raises(
        FileNotFoundError, match=r"synthetic\.fa\.gz\.gzi.*bgzip -r .*synthetic\.fa\.gz"
    ):
        FastaGenome(path)


def test_missing_file_raises_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        FastaGenome(tmp_path / "missing.fa")


def test_url_raises():
    with pytest.raises(InputError, match="the path is a URL; expected a local path"):
        FastaGenome("https://example.org/genome.fa")


def test_gzip_that_is_not_bgzf_raises(tmp_path):
    path = tmp_path / "plain.fa.gz"
    path.write_bytes(gzip.compress(b">a\nACGT\n"))
    with pytest.raises(InputError, match=r"compressed offset 0 is not BGZF.*bgzip"):
        FastaGenome(path)


@pytest.mark.parametrize(
    ("fai", "message"),
    [
        ("a\t10\t3\t60\n", "line 1 has 4 columns; expected 5 columns separated by tabs"),
        ("a\t10\t3\t60\t61\t0\n", "line 1 has 6 columns; expected 5 columns separated by tabs"),
        ("a 10 3 60 61\n", "line 1 has 1 column; expected 5 columns separated by tabs"),
        ("a\t10\t3\t60\t61\n\n", "line 2 has 1 column; expected 5 columns separated by tabs"),
        ("a\tten\t3\t60\t61\n", "line 1: length 'ten' is not a non-negative integer"),
        ("a\t10\t-3\t60\t61\n", "line 1: offset '-3' is not a non-negative integer"),
        ("a\t10\t3\t6.5\t61\n", "line 1: line bases '6.5' is not a non-negative integer"),
        ("a\t10\t3\t60\t61\r\n", "line 1: line width '61\\r' is not a non-negative integer"),
        ("a\t10\t3\t0\t1\n", "line 1: line bases 0 is below 1; expected 1 or more"),
        (
            "a\t10\t3\t60\t60\n",
            "line 1: line width 60 minus line bases 60 is 0; expected 1 (LF) or 2 (CRLF)",
        ),
        (
            "a\t10\t3\t60\t63\n",
            "line 1: line width 63 minus line bases 60 is 3; expected 1 (LF) or 2 (CRLF)",
        ),
        (
            "a\t10\t3\t60\t59\n",
            "line 1: line width 59 minus line bases 60 is -1; expected 1 (LF) or 2 (CRLF)",
        ),
        (
            "a\t10\t3\t60\t61\nb\t10\t3\t60\t6\n",
            "line 2: line width 6 minus line bases 60 is -54; expected",
        ),
        (
            "a\t9223372036854775808\t3\t60\t61\n",
            "line 1: length 9223372036854775808 is too large; expected at most 2^63 - 1",
        ),
        (
            "a\t10\t1000000000000000000000000000000\t60\t61\n",
            (
                "line 1: offset 1000000000000000000000000000000 is too large; "
                "expected at most 2^63 - 1"
            ),
        ),
        (
            "a\t10\t3\t9223372036854775808\t9223372036854775809\n",
            "line 1: line bases 9223372036854775808 is too large; expected at most 2^63 - 1",
        ),
        (
            "a\t10\t3\t60\t9223372036854775808\n",
            "line 1: line width 9223372036854775808 is too large; expected at most 2^63 - 1",
        ),
        (
            "a\t10\t9223372036854775807\t60\t61\n",
            (
                "line 1: the byte of the last base 9223372036854775816 is too large; "
                "expected at most 2^63 - 1"
            ),
        ),
        (
            "a\t10\t3\t60\t61\nb\t" + "9" * 5000 + "\t3\t60\t61\n",
            "line 2: length " + "9" * 40 + "... is too large; expected at most 2^63 - 1",
        ),
    ],
)
def test_invalid_fai_raises_and_names_the_line(tmp_path, fai, message):
    path = tmp_path / "genome.fa"
    path.write_bytes(b">a\nACGTACGTAC\n")
    (tmp_path / "genome.fa.fai").write_bytes(fai.encode())
    with pytest.raises(InputError, match=re.escape(f"{os.fspath(path)}.fai': {message}")):
        FastaGenome(path)


def test_fai_with_too_few_columns_raises(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_bytes(b">a\nACGT\n")
    (tmp_path / "genome.fa.fai").write_text("a\t4\t3\n")
    with pytest.raises(InputError, match=r"line 1.*5 columns"):
        FastaGenome(path)


def test_duplicate_contig_name_in_the_fai_raises(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_bytes(b">a\nACGT\n")
    (tmp_path / "genome.fa.fai").write_text("a\t4\t3\t4\t5\nb\t4\t3\t4\t5\na\t4\t3\t4\t5\n")
    fai = f"{os.fspath(path)}.fai"
    message = f"{fai!r}: contig 'a' occurs twice in the index; expected each name once"
    with pytest.raises(InputError, match=re.escape(message)):
        FastaGenome(path)


def test_fai_that_is_not_utf8_raises(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_bytes(b">a\nACGT\n")
    (tmp_path / "genome.fa.fai").write_bytes(b"\xff\t4\t3\t4\t5\n")
    with pytest.raises(InputError, match=r"genome\.fa\.fai.*is not valid UTF-8"):
        FastaGenome(path)


def test_short_last_line_is_valid(tmp_path):
    path = write_genome(tmp_path, {"a": "ACGTACGT" * 8 + "A"})
    with FastaGenome(path) as genome:
        assert genome.sequence("a", 61, 65) == ("ACGTACGT" * 8 + "A")[60:65]


@pytest.mark.parametrize("block_size", [None, 5])
@pytest.mark.parametrize(
    ("ending", "fai"), [("\n", "a\t9\t3\t1\t2\n"), ("\r\n", "a\t9\t4\t1\t3\n")]
)
def test_one_base_for_each_line(tmp_path, ending, fai, block_size):
    sequence = "ACGTRYKMN"
    data = f">a{ending}{''.join(base + ending for base in sequence)}".encode()
    path = tmp_path / "genome.fa"
    if block_size is None:
        path.write_bytes(data)
    else:
        packed, gzi = bgzf_bytes(data, block_size)
        path.write_bytes(packed)
        (tmp_path / "genome.fa.gzi").write_bytes(gzi)
    (tmp_path / "genome.fa.fai").write_text(fai)
    with FastaGenome(path) as genome:
        for start in range(1, 10):
            for end in range(start, 10):
                assert genome.sequence("a", start, end) == sequence[start - 1 : end]


def test_t10_range_crosses_a_block_boundary():
    (row,) = (r for r in golden_rows("annotation/genome_ranges.tsv") if r.label.startswith("t10_"))
    start, end = int(row.start), int(row.end)
    entry = {
        line.split("\t")[0]: [int(field) for field in line.split("\t")[1:]]
        for line in (ANNOTATION / "synthetic.fa.gz.fai").read_text().splitlines()
    }[row.contig]
    _length, offset, bases, width = entry

    def byte(base):
        return offset + (base - 1) // bases * width + (base - 1) % bases

    starts = [
        uncompressed for _compressed, uncompressed in read_gzi(ANNOTATION / "synthetic.fa.gz.gzi")
    ]
    assert [u for u in starts if byte(start) < u <= byte(end)] == [65273]
    with FastaGenome(BGZF) as genome:
        assert genome.sequence(row.contig, start, end) == row.bgzf


@pytest.fixture
def inflations(monkeypatch):
    """Record the length of the data of each block that the reader inflates."""
    calls = []
    real = zlib.decompress

    def decompress(data, *args, **kwargs):
        block = real(data, *args, **kwargs)
        calls.append(len(block))
        return block

    monkeypatch.setattr(zlib, "decompress", decompress)
    return calls


def test_t10_range_inflates_the_two_blocks_that_it_touches(inflations):
    with FastaGenome(BGZF) as genome:
        genome.codes("1", 64000, 64450)
    assert inflations == [65273, 65270]


def test_a_range_inflates_only_the_blocks_from_its_first_byte_to_its_last_byte(
    tmp_path, inflations
):
    path = write_genome(tmp_path, {"a": "ACGT" * 100}, block_size=50)
    with FastaGenome(path) as genome:
        genome.codes("a", 61, 70)
        assert inflations == [50]
        genome.codes("a", 61, 150)
        assert inflations == [50, 50, 50]


def test_the_last_block_is_kept_for_the_next_request(tmp_path, inflations):
    path = write_genome(tmp_path, {"a": "ACGT" * 100}, block_size=50)
    with FastaGenome(path) as genome:
        genome.codes("a", 1, 10)
        genome.codes("a", 1, 10)
        genome.codes("a", 11, 20)
        assert len(inflations) == 1
        genome.codes("a", 1, 100)
        assert len(inflations) == 3
        genome.codes("a", 97, 110)
        assert len(inflations) == 3


def test_blocks_after_the_range_are_not_read(tmp_path):
    path = write_genome(tmp_path, {"a": "ACGT" * 100}, block_size=50)
    compressed, uncompressed = read_gzi(tmp_path / "genome.fa.gzi")[2]
    assert uncompressed == 150
    data = bytearray(path.read_bytes())
    data[compressed + 20] ^= 0xFF
    path.write_bytes(bytes(data))
    with FastaGenome(path) as genome:
        assert genome.sequence("a", 1, 60) == ("ACGT" * 100)[:60]
        message = rf"compressed offset {compressed} is damaged"
        with pytest.raises(InputError, match=message):
            genome.codes("a", 1, 150)


@pytest.mark.parametrize(
    ("edit", "request_", "offset", "reason"),
    [
        ("flip", ("1", 64000, 64450), 18906, r"(its data does not inflate|its checksum does not)"),
        ("cut", ("X", 1, 100), 39090, "the file ends inside it"),
    ],
)
def test_damaged_block_raises_and_names_the_offset(tmp_path, edit, request_, offset, reason):
    path = copy_as(tmp_path, "synthetic.fa.gz", "synthetic.fa.gz")
    data = bytearray(path.read_bytes())
    if edit == "flip":
        data[25000] ^= 0xFF
    else:
        del data[45000:]
    path.write_bytes(bytes(data))
    with FastaGenome(path) as genome:
        assert genome.sequence("1", 1, 100) == fasta_sequences()["1"][:100]
        with pytest.raises(InputError, match=rf"compressed offset {offset} is damaged: {reason}"):
            genome.codes(*request_)


def gzip_block_header(size):
    extra = struct.pack("<2sHH", b"BC", 2, size - 1)
    return struct.pack("<2sBBIBBH", b"\x1f\x8b", 8, 4, 0, 0, 255, len(extra)) + extra


@pytest.mark.parametrize("size", [16, 17, 18, 25])
def test_block_size_below_its_header_and_trailer_raises_without_a_read_of_the_file(tmp_path, size):
    path = tmp_path / "genome.fa"
    path.write_bytes(gzip_block_header(size) + bytes(8 << 20))
    message = (
        f"{os.fspath(path)!r}: the BGZF block at compressed offset 0 is damaged: its size {size} "
        "is too small for its header and trailer; expected an intact file"
    )
    tracemalloc.start()
    try:
        with pytest.raises(InputError, match=re.escape(message)):
            FastaGenome(path)
        _current, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert peak < 1 << 20


@pytest.mark.parametrize(
    "subfields", [b"XX\x03\x00abc", b"XX\x00\x00", b"RG\x02\x00ab" + b"ZZ\x01\x00z"]
)
def test_extra_fields_before_bc_are_skipped(tmp_path, subfields):
    path = write_genome(tmp_path, {"a": "ACGT" * 100}, block_size=50, subfields=subfields)
    with FastaGenome(path) as genome:
        assert genome.sequence("a", 1, 400) == "ACGT" * 100


@pytest.mark.parametrize(
    "extra", [b"XX\x02\x00\x00\x00", b"BC\x02\x00\x00", b"BC\x03\x00\x00\x00\x00", b"BC\x02"]
)
def test_gzip_with_an_extra_field_but_no_bc_is_not_bgzf(tmp_path, extra):
    header = struct.pack("<2sBBIBBH", b"\x1f\x8b", 8, 4, 0, 0, 255, len(extra))
    path = tmp_path / "genome.fa"
    path.write_bytes(header + extra + b"\x03\x00" + bytes(8))
    with pytest.raises(InputError, match=r"compressed offset 0 is not BGZF.*bgzip"):
        FastaGenome(path)


@pytest.mark.parametrize(
    "head",
    [
        b"\x1f\x8b",
        b"\x1f\x8b\x08\x04\x00",
        struct.pack("<2sBBIBBH", b"\x1f\x8b", 8, 4, 0, 0, 255, 6) + b"BC\x02",
    ],
)
def test_gzip_that_ends_inside_the_header_raises(tmp_path, head):
    path = tmp_path / "genome.fa"
    path.write_bytes(head)
    message = (
        f"{os.fspath(path)!r}: the file ends before the BGZF block at compressed offset 0 has a "
        "complete header; expected an intact file and the .gzi of this file"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        FastaGenome(path)


def test_file_that_ends_inside_the_header_of_a_later_block_raises(tmp_path):
    path = copy_as(tmp_path, "synthetic.fa.gz", "synthetic.fa.gz")
    path.write_bytes(path.read_bytes()[: 39090 + 5])
    message = (
        f"{os.fspath(path)!r}: the file ends before the BGZF block at compressed offset 39090 "
        "has a complete header; expected an intact file and the .gzi of this file"
    )
    with FastaGenome(path) as genome:
        assert genome.sequence("1", 1, 4) == fasta_sequences()["1"][:4]
        with pytest.raises(InputError, match=re.escape(message)):
            genome.codes("X", 1, 100)


@pytest.mark.parametrize("offset", [47388, 10**8, 2**63, 2**64 - 1])
def test_gzi_with_an_offset_that_is_not_below_the_file_size_raises_at_open(tmp_path, offset):
    path = copy_as(tmp_path, "synthetic.fa.gz", "synthetic.fa.gz")
    assert path.stat().st_size == 47388
    gzi = bytearray((ANNOTATION / "synthetic.fa.gz.gzi").read_bytes())
    gzi[8:16] = struct.pack("<Q", offset)
    (tmp_path / "synthetic.fa.gz.gzi").write_bytes(bytes(gzi))
    message = (
        f"{os.fspath(path)}.gzi': compressed offset {offset} is not below the size of the file, "
        f"47388 bytes; expected the .gzi of this file. Make it with `bgzip -r {os.fspath(path)}`"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        FastaGenome(path)


@pytest.mark.parametrize(
    ("edit", "reason"),
    [
        ("data", r"its data does not inflate \(.*invalid block type\)"),
        ("checksum", "its checksum does not match its data"),
    ],
)
def test_block_with_bad_data_or_a_bad_checksum_raises(tmp_path, edit, reason):
    path = write_genome(tmp_path, {"a": "ACGT" * 10}, block_size=1000)
    data = bytearray(path.read_bytes())
    (total,) = struct.unpack_from("<H", data, 16)
    if edit == "data":
        data[18:22] = b"\xff" * 4
    else:
        data[total + 1 - 8] ^= 0xFF
    path.write_bytes(bytes(data))
    message = rf"compressed offset 0 is damaged: {reason}; expected an intact file"
    with FastaGenome(path) as genome, pytest.raises(InputError, match=message):
        genome.codes("a", 1, 4)


def test_gzi_that_gives_another_block_length_raises(tmp_path):
    path = copy_as(tmp_path, "synthetic.fa.gz", "synthetic.fa.gz")
    gzi = bytearray((ANNOTATION / "synthetic.fa.gz.gzi").read_bytes())
    gzi[16:24] = struct.pack("<Q", 65272)
    (tmp_path / "synthetic.fa.gz.gzi").write_bytes(bytes(gzi))
    message = (
        r"synthetic\.fa\.gz\.gzi': the block at compressed offset 0 has 65273 bytes, "
        r"but the \.gzi gives 65272"
    )
    with FastaGenome(path) as genome, pytest.raises(InputError, match=message):
        genome.codes("1", 1, 4)


def test_gzi_that_points_into_a_block_raises(tmp_path):
    path = copy_as(tmp_path, "synthetic.fa.gz", "synthetic.fa.gz")
    gzi = bytearray((ANNOTATION / "synthetic.fa.gz.gzi").read_bytes())
    gzi[8:16] = struct.pack("<Q", 18907)
    (tmp_path / "synthetic.fa.gz.gzi").write_bytes(bytes(gzi))
    with FastaGenome(path) as genome:
        assert genome.sequence("1", 1, 4) == fasta_sequences()["1"][:4]
        with pytest.raises(InputError, match=r"compressed offset 18907 is not BGZF.*bgzip"):
            genome.codes("1", 64000, 64450)


@pytest.mark.parametrize(("cut", "size"), [(3, "37 bytes"), (12, "28 bytes"), (39, "1 byte")])
def test_gzi_with_a_wrong_size_raises(tmp_path, cut, size):
    path = copy_as(tmp_path, "synthetic.fa.gz", "synthetic.fa.gz")
    gzi = (ANNOTATION / "synthetic.fa.gz.gzi").read_bytes()
    (tmp_path / "synthetic.fa.gz.gzi").write_bytes(gzi[:-cut])
    message = rf"synthetic\.fa\.gz\.gzi': the file has {size}; expected 40 bytes for 2 blocks"
    with pytest.raises(InputError, match=message):
        FastaGenome(path)


def test_gzi_of_a_file_with_one_block_has_no_entries(tmp_path):
    path = write_genome(tmp_path, {"a": "ACGT" * 10}, block_size=1000)
    assert read_gzi(tmp_path / "genome.fa.gzi") == []
    with FastaGenome(path) as genome:
        assert genome.sequence("a", 1, 40) == "ACGT" * 10


@pytest.mark.parametrize("name", NAMES)
def test_reads_from_a_read_only_directory(tmp_path, name):
    directory = tmp_path / "refs"
    directory.mkdir()
    path = copy_as(directory, name, name)
    expected = sorted(entry.name for entry in directory.iterdir())
    for entry in directory.iterdir():
        entry.chmod(0o444)
    directory.chmod(0o555)
    try:
        with FastaGenome(path) as genome:
            assert genome.sequence("X", 1, 5003) == fasta_sequences()["X"]
        assert sorted(entry.name for entry in directory.iterdir()) == expected
    finally:
        directory.chmod(0o755)


@pytest.mark.parametrize("name", NAMES)
def test_close_releases_the_file_and_a_later_request_raises(name):
    genome = FastaGenome(ANNOTATION / name)
    assert genome.sequence("X", 1, 4) == "CGTC"
    genome.close()
    genome.close()
    with pytest.raises(ValueError, match="closed file"):
        genome.codes("X", 1, 4)


@pytest.mark.parametrize("name", NAMES)
def test_with_block_closes_the_file(name):
    with FastaGenome(ANNOTATION / name) as genome:
        pass
    with pytest.raises(ValueError, match="closed file"):
        genome.codes("X", 1, 4)


@pytest.mark.parametrize("name", NAMES)
def test_resolve_maps_the_contigs_of_the_file(name):
    with FastaGenome(ANNOTATION / name) as genome:
        assert genome.resolve("chr1") == "1"
        assert genome.resolve("chrX") == "X"
        assert genome.resolve("SCAF1.1") == "SCAF1.1"
        assert genome.resolve("Z") is None


@pytest.mark.parametrize(
    ("name", "reader"),
    [
        ("synthetic.2bit", TwoBitGenome),
        ("synthetic.v1.2bit", TwoBitGenome),
        ("synthetic.fa", FastaGenome),
        ("synthetic.fa.gz", FastaGenome),
    ],
)
def test_open_genome_selects_the_reader_from_the_first_bytes(tmp_path, name, reader):
    path = copy_as(tmp_path, name, "genome")
    with open_genome(path) as genome:
        assert type(genome) is reader
        assert genome.sequence("X", 1, 8) == "CGTCCCGG"


def test_open_genome_opens_a_big_endian_2bit_file(tmp_path):
    path = tmp_path / "genome"
    path.write_bytes(struct.pack(">IIII", 0x1A412743, 0, 0, 0))
    with open_genome(path) as genome:
        assert type(genome) is TwoBitGenome
        assert dict(genome.contigs) == {}


def test_open_genome_accepts_a_string_path():
    with open_genome(str(PLAIN)) as genome:
        assert type(genome) is FastaGenome


def test_open_genome_refuses_other_bytes(tmp_path):
    path = tmp_path / "genome"
    path.write_bytes(b"This is a text file, not a genome.\n")
    message = (
        f"{os.fspath(path)!r}: not a 2bit or FASTA file; the first bytes are 54686973, "
        "expected 4327411a or 1a412743 (2bit), 1f8b (gzip), or 3e (FASTA)"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        open_genome(path)


def test_open_genome_refuses_an_empty_file(tmp_path):
    path = tmp_path / "genome"
    path.write_bytes(b"")
    message = f"{os.fspath(path)!r}: the file is empty; expected a 2bit or FASTA file"
    with pytest.raises(InputError, match=re.escape(message)):
        open_genome(path)


def test_open_genome_of_a_gzip_file_that_is_not_bgzf_raises(tmp_path):
    path = tmp_path / "genome"
    path.write_bytes(gzip.compress(b">a\nACGT\n"))
    with pytest.raises(InputError, match=r"is not BGZF.*bgzip"):
        open_genome(path)


def test_open_genome_of_a_fasta_file_without_a_fai_names_the_command(tmp_path):
    path = tmp_path / "genome"
    shutil.copy(PLAIN, path)
    with pytest.raises(FileNotFoundError, match=r"genome\.fai.*samtools faidx"):
        open_genome(path)


def test_open_genome_missing_file_raises_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        open_genome(tmp_path / "missing")


def test_open_genome_refuses_a_url():
    with pytest.raises(InputError, match="the path is a URL; expected a local path"):
        open_genome("https://example.org/genome.2bit")
