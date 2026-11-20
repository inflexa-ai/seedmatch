import re
from pathlib import Path

import numpy as np
import pytest
from helpers import DATA, sirna_seqs

from seedmatch.errors import InputError
from seedmatch.sequence import SequenceSet, read_fasta, write_fasta


def naive_fasta_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            records.append((line[1:].split()[0], ""))
        elif line:
            name, body = records[-1]
            records[-1] = (name, body + line)
    return records


def test_read_fasta_keeps_each_record_of_the_test_data():
    seqs = read_fasta(DATA / "Test_seqs.fa")
    expected = naive_fasta_records(DATA / "Test_seqs.fa")
    assert list(seqs.ids) == [name for name, _ in expected]
    assert [seqs.sequence(i) for i in range(len(seqs))] == [body for _, body in expected]
    assert seqs.gene_ids == seqs.ids


def test_sirna_seqs_is_the_test_data_and_is_shared():
    expected = naive_fasta_records(DATA / "Test_seqs.fa")
    assert sirna_seqs().ids == tuple(name for name, _ in expected)
    assert sirna_seqs() is sirna_seqs()


def test_empty_set_has_no_sequences():
    s = SequenceSet.from_records([])
    assert len(s) == 0
    assert s.offsets.tolist() == [0]


def test_empty_set_has_typed_empty_buffers():
    s = SequenceSet.from_records([])
    assert s.ids == ()
    assert s.gene_ids == ()
    assert s.codes.dtype == np.uint8
    assert s.codes.shape == (0,)
    assert s.lengths().tolist() == []


def test_from_records_builds_one_buffer_with_n_plus_one_offsets():
    s = SequenceSet.from_records([("a", "ACGT"), ("b", ""), ("c", "nnRY")])
    assert len(s) == 3
    assert s.ids == ("a", "b", "c")
    assert s.offsets.tolist() == [0, 4, 4, 8]
    assert s.lengths().tolist() == [4, 0, 4]
    assert s.sequence(0) == "ACGT"
    assert s.sequence(1) == ""
    assert s.sequence(2) == "NNRY"


def test_buffers_have_the_dtypes_and_layout_that_the_kernels_need():
    s = SequenceSet.from_records([("a", "ACGT"), ("b", "TT")])
    assert s.codes.dtype == np.uint8
    assert s.codes.flags.c_contiguous
    assert s.offsets.dtype == np.int64
    assert s.offsets.flags.c_contiguous
    assert s.lengths().dtype == np.int64


def test_soft_masked_bases_are_stored_as_uppercase():
    s = SequenceSet.from_records([("a", "acgtNnRy")])
    assert s.sequence(0) == "ACGTNNRY"


def test_buffers_are_read_only():
    s = SequenceSet.from_records([("a", "ACGT")])
    with pytest.raises(ValueError, match="read-only"):
        s.codes[0] = 1
    with pytest.raises(ValueError, match="read-only"):
        s.offsets[0] = 1


def test_set_cannot_be_reassigned():
    s = SequenceSet.from_records([("a", "ACGT")])
    with pytest.raises(AttributeError):
        s.ids = ("b",)


def test_from_records_accepts_a_generator():
    s = SequenceSet.from_records((name, "AC") for name in "xyz")
    assert s.ids == ("x", "y", "z")


def test_duplicate_ids_are_kept():
    s = SequenceSet.from_records([("a", "AC"), ("a", "GT")])
    assert s.ids == ("a", "a")
    assert [s.sequence(0), s.sequence(1)] == ["AC", "GT"]


def test_gene_ids_default_to_the_ids():
    s = SequenceSet.from_records([("t1", "AC"), ("t2", "GT")])
    assert s.gene_ids == s.ids


def test_gene_ids_group_several_sequences_under_one_gene():
    s = SequenceSet.from_records([("t1", "AC"), ("t2", "GT")], gene_ids=["g1", "g1"])
    assert s.ids == ("t1", "t2")
    assert s.gene_ids == ("g1", "g1")


def test_gene_ids_with_the_wrong_length_name_both_lengths():
    with pytest.raises(InputError, match=r"gene_ids has 1 entries.*the 3 sequences"):
        SequenceSet.from_records([("a", "A"), ("b", "C"), ("c", "G")], gene_ids=["g"])


def test_uracil_in_a_record_names_the_record_and_the_position():
    with pytest.raises(InputError, match=r"record 'tx2'.*'U' at position 3"):
        SequenceSet.from_records([("tx1", "ACGT"), ("tx2", "ACUG")])


@pytest.mark.parametrize("index", [-1, 2, 10])
def test_sequence_rejects_an_index_outside_the_set(index):
    s = SequenceSet.from_records([("a", "AC"), ("b", "GT")])
    with pytest.raises(InputError, match=rf"index {index}.*0 to 1"):
        s.sequence(index)


def test_sequence_of_an_empty_set_rejects_each_index():
    with pytest.raises(InputError, match=r"index 0.*empty"):
        SequenceSet.from_records([]).sequence(0)


def make_parts():
    codes = np.array([1, 2, 4, 8], dtype=np.uint8)
    offsets = np.array([0, 2, 4], dtype=np.int64)
    return {"ids": ("a", "b"), "gene_ids": ("a", "b"), "codes": codes, "offsets": offsets}


def test_direct_construction_with_consistent_parts_gives_a_read_only_set():
    s = SequenceSet(**make_parts())
    assert [s.sequence(0), s.sequence(1)] == ["AC", "GT"]
    assert not s.codes.flags.writeable
    assert not s.offsets.flags.writeable


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"gene_ids": ("a",)}, r"gene_ids has 1 entries.*the 2 sequences"),
        (
            {"codes": np.array([1, 2, 4, 8], dtype=np.int64)},
            r"codes has dtype int64; expected uint8",
        ),
        (
            {"offsets": np.array([0, 2, 4], dtype=np.int32)},
            r"offsets has dtype int32; expected int64",
        ),
        ({"codes": np.array([[1, 2], [4, 8]], dtype=np.uint8)}, r"codes.*1-D"),
        (
            {"codes": np.array([1, 0, 2, 0, 4, 0, 8, 0], dtype=np.uint8)[::2]},
            r"codes.*C-contiguous",
        ),
        ({"offsets": np.array([0, 4], dtype=np.int64)}, r"offsets has 2 entries.*3"),
        ({"offsets": np.array([1, 2, 4], dtype=np.int64)}, r"offsets must start at 0"),
        ({"offsets": np.array([0, 2, 3], dtype=np.int64)}, r"offsets must end at 4"),
        (
            {"offsets": np.array([0, 5, 4], dtype=np.int64)},
            r"offsets\[2\] 4 is less than offsets\[1\] 5; expected offsets that do not decrease",
        ),
    ],
)
def test_direct_construction_rejects_parts_that_do_not_fit_together(change, message):
    with pytest.raises(InputError, match=message):
        SequenceSet(**{**make_parts(), **change})


@pytest.mark.parametrize("bad", [0, 16, 200])
def test_direct_construction_rejects_a_code_that_is_not_a_mask(bad):
    codes = np.array([1, 2, bad, 8], dtype=np.uint8)
    with pytest.raises(InputError, match=rf"invalid mask {bad} at index 2; expected 1 to 15"):
        SequenceSet(**{**make_parts(), "codes": codes})


def test_direct_construction_accepts_each_mask_from_1_to_15():
    codes = np.arange(1, 16, dtype=np.uint8)
    offsets = np.array([0, 15], dtype=np.int64)
    s = SequenceSet(ids=("a",), gene_ids=("a",), codes=codes, offsets=offsets)
    assert len(s.sequence(0)) == 15


def test_d016_read_fasta_takes_the_id_up_to_the_first_white_space(tmp_path):
    """D-016: R keeps the whole header line as the name, but the port ends the ID at white space."""
    ensembl = ">ENSRNOG00000016275 chromosome:mRatBN7.2:1:1:100:1 gene_biotype:protein_coding"
    (tmp_path / "x.fa").write_text(
        f">a desc here\nAC\n>b\tdesc\nGT\n>  c\nTT\n>\nAA\n{ensembl}\nCC\n"
    )
    assert read_fasta(tmp_path / "x.fa").ids == ("a", "b", "", "", "ENSRNOG00000016275")


def test_read_fasta_joins_the_lines_and_ignores_empty_lines_and_carriage_returns(tmp_path):
    (tmp_path / "x.fa").write_bytes(b"\n>a\r\nAC\r\n\r\nGT\r\n\n>b\r\n>c\nacgn\n\n")
    seqs = read_fasta(tmp_path / "x.fa")
    assert seqs.ids == ("a", "b", "c")
    assert [seqs.sequence(i) for i in range(3)] == ["ACGT", "", "ACGN"]


def test_read_fasta_accepts_a_str_path(tmp_path):
    (tmp_path / "x.fa").write_text(">a\nAC\n")
    assert read_fasta(str(tmp_path / "x.fa")).ids == ("a",)


def test_read_fasta_of_an_empty_file_gives_an_empty_set(tmp_path):
    (tmp_path / "x.fa").write_text("")
    assert len(read_fasta(tmp_path / "x.fa")) == 0


def test_read_fasta_of_a_directory_raises_input_error(tmp_path):
    message = f"{str(tmp_path)!r}: the path is a directory; expected a file"
    with pytest.raises(InputError, match=re.escape(message)):
        read_fasta(tmp_path)


def test_read_fasta_rejects_a_sequence_line_before_the_first_header(tmp_path):
    (tmp_path / "x.fa").write_text("\nACGT\n>a\nAC\n")
    with pytest.raises(InputError, match=r"x\.fa.*line 2.*before the first header"):
        read_fasta(tmp_path / "x.fa")


def test_read_fasta_names_the_record_of_a_bad_letter(tmp_path):
    (tmp_path / "x.fa").write_text(">a\nACGT\n>b desc\nAC\nGUT\n")
    with pytest.raises(InputError, match=r"record 'b'.*'U' at position 4"):
        read_fasta(tmp_path / "x.fa")


@pytest.mark.parametrize("letter", ["-", "+", "."])
def test_d017_read_fasta_rejects_gap_mask_and_unknown_letters(tmp_path, letter):
    """D-017: readDNAStringSet keeps these letters, but the port stops with InputError."""
    (tmp_path / "x.fa").write_text(f">x\nAC{letter}GT\n")
    with pytest.raises(InputError, match=rf"record 'x'.*'{re.escape(letter)}' at position 3"):
        read_fasta(tmp_path / "x.fa")


def test_read_fasta_rejects_white_space_inside_a_sequence_line(tmp_path):
    (tmp_path / "x.fa").write_text(">a\nAC GT\n")
    with pytest.raises(InputError, match=r"record 'a'.*' ' at position 3"):
        read_fasta(tmp_path / "x.fa")


def test_read_fasta_rejects_a_file_that_is_not_utf8(tmp_path):
    (tmp_path / "x.fa").write_bytes(b">a\xff\nAC\n")
    with pytest.raises(InputError, match=r"x\.fa.*UTF-8"):
        read_fasta(tmp_path / "x.fa")


def test_read_fasta_of_a_missing_file_raises_the_os_error(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_fasta(tmp_path / "missing.fa")


def test_write_then_read_gives_the_same_set(tmp_path):
    seqs = SequenceSet.from_records([("a", "ACGT" * 50), ("b", ""), ("c", "NNRY")])
    write_fasta(seqs, tmp_path / "x.fa", width=7)
    back = read_fasta(tmp_path / "x.fa")
    assert back.ids == seqs.ids
    assert np.array_equal(back.codes, seqs.codes)
    assert max(len(line) for line in (tmp_path / "x.fa").read_text().splitlines()) == 7


def test_write_fasta_layout(tmp_path):
    seqs = SequenceSet.from_records([("a", "ACGTAC"), ("b", ""), ("c", "NNRYA")])
    write_fasta(seqs, tmp_path / "x.fa", width=3)
    assert (tmp_path / "x.fa").read_text() == ">a\nACG\nTAC\n>b\n>c\nNNR\nYA\n"


def test_write_fasta_wraps_at_80_by_default(tmp_path):
    seqs = SequenceSet.from_records([("a", "A" * 200)])
    write_fasta(seqs, tmp_path / "x.fa")
    lines = (tmp_path / "x.fa").read_text().splitlines()
    assert [len(line) for line in lines] == [2, 80, 80, 40]


def test_write_fasta_of_the_test_data_reads_back_the_same(tmp_path):
    seqs = sirna_seqs()
    write_fasta(seqs, tmp_path / "x.fa")
    back = read_fasta(tmp_path / "x.fa")
    assert back.ids == seqs.ids
    assert np.array_equal(back.codes, seqs.codes)
    assert np.array_equal(back.offsets, seqs.offsets)


def test_write_fasta_of_an_empty_set_makes_an_empty_file(tmp_path):
    write_fasta(SequenceSet.from_records([]), tmp_path / "x.fa")
    assert (tmp_path / "x.fa").read_text() == ""


@pytest.mark.parametrize("width", [0, -5])
def test_write_fasta_rejects_a_width_below_one(tmp_path, width):
    seqs = SequenceSet.from_records([("a", "AC")])
    with pytest.raises(InputError, match=rf"width {width}.*at least 1"):
        write_fasta(seqs, tmp_path / "x.fa", width=width)
    assert not (tmp_path / "x.fa").exists()


@pytest.mark.parametrize("bad_id", ["a b", "a\tb", "a\nb", " a"])
def test_write_fasta_rejects_an_id_with_white_space(tmp_path, bad_id):
    seqs = SequenceSet.from_records([("ok", "AC"), (bad_id, "GT")])
    with pytest.raises(InputError, match=r"id .*white space"):
        write_fasta(seqs, tmp_path / "x.fa")
    assert not (tmp_path / "x.fa").exists()


@pytest.mark.parametrize(
    ("ids", "gene_ids", "message"),
    [
        ([1], ["g"], "ids holds 1 of type int; expected a str"),
        (["a"], [1], "gene_ids holds 1 of type int; expected a str"),
    ],
)
def test_sequence_set_rejects_an_id_that_is_not_text(ids, gene_ids, message):
    records = [(name, "ACGT") for name in ids]
    with pytest.raises(InputError, match=re.escape(message)):
        SequenceSet.from_records(records, gene_ids=gene_ids)


def test_from_records_rejects_one_str_as_the_gene_ids():
    records = [("a", "A"), ("b", "C"), ("c", "G")]
    message = "gene_ids 'abc' is a str; expected a sequence of str"
    with pytest.raises(InputError, match=re.escape(message)):
        SequenceSet.from_records(records, gene_ids="abc")
