import os
import re

import numpy as np
import pytest
from annotation_cases import ALL, CASES, ENSDB, NOT_CODING, copy_ensdb, fixture_reference
from helpers import GOLDEN, golden

import seedmatch.annotation.ensdb
from seedmatch.annotation.ensdb import AnnotationDb, ExonRows
from seedmatch.annotation.features import Features, derive_features
from seedmatch.annotation.filters import FeatureType, PerGene, TranscriptFilter
from seedmatch.errors import InputError

G01 = "ENSSYNG00000000001"
G02 = "ENSSYNG00000000002"
G03 = "ENSSYNG00000000003"
G08 = "ENSSYNG00000000008"
T01 = "ENSSYNT00000000001"
T02 = "ENSSYNT00000000002"
T03 = "ENSSYNT00000000003"
T04 = "ENSSYNT00000000004"
T05 = "ENSSYNT00000000005"
T07 = "ENSSYNT00000000007"
T12 = "ENSSYNT00000000012"
EVERYTHING = CASES["exons_everything"].transcript_filter


def ranges_of(frame, name):
    rows = frame[frame["id"] == name]
    return list(zip(rows["start"].tolist(), rows["end"].tolist(), strict=True))


def golden_ranges(frame, name):
    rows = frame[frame["name"] == name]
    return [(int(start), int(end)) for start, end in zip(rows["start"], rows["end"], strict=True)]


def hand_rows(*transcripts):
    """Return the ExonRows of (tx_id, gene_id, strand, cds, exons) records, in the given order.

    The exons are in transcript order, and `cds` is None for a transcript without a CDS.
    """
    columns = {name: [] for name in ExonRows.__dataclass_fields__}
    for tx_row, (tx_id, gene_id, strand, cds, exons) in enumerate(transcripts, start=1):
        cds_start, cds_end = cds or (0, 0)
        for start, end in exons:
            for name, value in (
                ("tx_id", tx_id),
                ("gene_id", gene_id),
                ("seq_name", "c"),
                ("strand", strand),
                ("cds_start", cds_start),
                ("cds_end", cds_end),
                ("exon_start", start),
                ("exon_end", end),
                ("tx_seq_start", min(s for s, _ in exons)),
                ("tx_row", tx_row),
            ):
                columns[name].append(value)
    dtypes = {"tx_id": object, "gene_id": object, "seq_name": object, "strand": np.int8}
    return ExonRows(
        **{
            name: np.array(values, dtype=dtypes.get(name, np.int64))
            for name, values in columns.items()
        }
    )


# T1 and T2 have the CDS on exactly one exon, so each exon edge is a CDS edge. T3 has one exon,
# and T4 has no CDS.
EDGES = hand_rows(
    ("T1", "G1", 1, (30, 40), ((10, 20), (30, 40), (50, 60))),
    ("T2", "G2", -1, (30, 40), ((50, 60), (30, 40), (10, 20))),
    ("T3", "G3", 1, (120, 180), ((100, 200),)),
    ("T4", "G4", 1, None, ((300, 400),)),
)


@pytest.mark.parametrize("config", CASES)
def test_features_match_r(config):
    case = CASES[config]
    with fixture_reference() as ref:
        got = ref.features(case.feature, case.transcript_filter, per_gene=case.per_gene).to_frame()
    expected = golden(f"annotation/{config}_features.tsv")
    assert got[["id", "rank", "contig", "start", "end", "strand"]].astype(str).values.tolist() == (
        expected[["name", "rank", "seqnames", "start", "end", "strand"]].values.tolist()
    )


def test_union_gives_the_minus_strand_ranges_in_transcript_order_D027():
    r_ranges = golden_ranges(golden("annotation/utr3_union_r_features.tsv"), G02)
    assert len(r_ranges) > 1
    assert r_ranges == sorted(r_ranges)
    with fixture_reference() as ref:
        got = ref.features(FeatureType.UTR3, ALL, per_gene=PerGene.UNION).to_frame()
    g02 = got[got["id"] == G02]
    assert set(g02["strand"]) == {"-"}
    assert ranges_of(got, G02) == r_ranges[::-1]


def test_longest_skips_a_gene_without_a_coding_transcript_D028():
    path = "annotation/utr3_longest_noncoding_features.tsv"
    notes = (GOLDEN / path).read_text()
    assert "# unpatched get_longest_utr stops: subscript contains invalid names" in notes
    tx_gene = golden("annotation/transcripts_utr3_longest_noncoding.tsv").set_index("tx_id")
    r_genes = set(tx_gene.loc[golden(path)["name"], "gene_id"])
    with fixture_reference() as ref:
        listed = ref.transcripts(NOT_CODING).set_index("tx_id")
        got = ref.features(FeatureType.UTR3, NOT_CODING, per_gene=PerGene.LONGEST)
    # The filter keeps G03, but its only transcript has no CDS.
    assert listed.loc[T05, "gene_id"] == G03
    assert listed["cds_start"][listed["gene_id"] == G03].isna().all()
    assert G03 not in got.gene_ids
    assert r_genes <= set(got.gene_ids)
    assert list(got.gene_ids) == sorted(set(listed["gene_id"]) - {G03})


def test_per_gene_applies_to_the_exons_D029():
    r_longest = golden("annotation/exons_longest_r_features.tsv")
    r_all = golden("annotation/exons_all_features.tsv")
    assert r_longest.equals(r_all)
    tx_gene = golden("annotation/transcripts_exons_all.tsv").set_index("tx_id")["gene_id"]
    width = {
        tx_id: sum(end - start + 1 for start, end in golden_ranges(r_all, tx_id))
        for tx_id in tx_gene.index
    }
    with fixture_reference() as ref:
        got = ref.features(FeatureType.EXONS, ALL, per_gene=PerGene.LONGEST)
    frame = got.to_frame()
    assert list(got.gene_ids) == sorted(set(tx_gene))
    for tx_id, gene_id in zip(got.ids, got.gene_ids, strict=True):
        assert tx_gene[tx_id] == gene_id
        assert width[tx_id] == max(width[other] for other in tx_gene.index[tx_gene == gene_id])
        assert ranges_of(frame, tx_id) == golden_ranges(r_all, tx_id)


def test_longest_tie_goes_to_the_lower_tx_row():
    g01 = TranscriptFilter(canonical=False, gene_ids=(G01,))
    with fixture_reference() as ref:
        each = ref.features(FeatureType.UTR3, g01)
        longest = ref.features(FeatureType.UTR3, g01, per_gene=PerGene.LONGEST)
        rows = ref.annotation.exon_rows(g01, coding_only=True)
    width = dict(zip(each.ids, each.widths().tolist(), strict=True))
    tx_row = dict(zip(rows.tx_id, rows.tx_row.tolist(), strict=True))
    tx_seq_start = dict(zip(rows.tx_id, rows.tx_seq_start.tolist(), strict=True))
    assert width[T01] == width[T02]
    assert tx_seq_start[T01] == tx_seq_start[T02]
    assert tx_row[T02] < tx_row[T01]
    assert longest.ids == (T02,)


def test_union_merges_touching_ranges_but_not_a_gap_of_one_base():
    with fixture_reference() as ref:
        each = ref.features(FeatureType.UTR3, ALL).to_frame()
        union = ref.features(FeatureType.UTR3, ALL, per_gene=PerGene.UNION).to_frame()
    assert (22000, 22800) in ranges_of(each, T03)
    assert (22801, 23000) in ranges_of(each, T04)
    # Base 21999 is in no range of G02.
    assert (21500, 21998) in ranges_of(each, T04)
    assert ranges_of(union, G02) == [(25000, 25700), (22000, 23000), (21500, 21998)]


def merged(ranges):
    """Return the union of the ranges in ascending order, with touching ranges as one range."""
    union = []
    for start, end in sorted(ranges):
        if union and start <= union[-1][1] + 1:
            union[-1] = (union[-1][0], max(union[-1][1], end))
        else:
            union.append((start, end))
    return union


@pytest.mark.parametrize("transcript_filter", [ALL, NOT_CODING], ids=["coding", "not coding"])
@pytest.mark.parametrize("feature", FeatureType)
def test_longest_and_union_agree_with_the_transcripts_for_each_feature_type_D029(
    feature, transcript_filter
):
    with fixture_reference() as ref:
        each = ref.features(feature, transcript_filter)
        longest = ref.features(feature, transcript_filter, per_gene=PerGene.LONGEST)
        union = ref.features(feature, transcript_filter, per_gene=PerGene.UNION)
        rows = ref.annotation.exon_rows(transcript_filter, coding_only=False)
    tie_key = {
        tx_id: (start, tx_row)
        for tx_id, start, tx_row in zip(
            rows.tx_id, rows.tx_seq_start.tolist(), rows.tx_row.tolist(), strict=True
        )
    }
    frame = each.to_frame()
    ranges = {tx_id: ranges_of(frame, tx_id) for tx_id in each.ids}
    width = dict(zip(each.ids, each.widths().tolist(), strict=True))
    strand = dict(zip(each.gene_ids, each.strands.tolist(), strict=True))
    transcripts = {}
    for tx_id, gene_id in zip(each.ids, each.gene_ids, strict=True):
        transcripts.setdefault(gene_id, []).append(tx_id)
    genes = sorted(transcripts)
    assert len(genes) > 1
    assert any(len(tx_ids) > 1 for tx_ids in transcripts.values())

    wanted = []
    for gene_id in genes:
        in_order = sorted(transcripts[gene_id], key=tie_key.__getitem__)
        maximum = max(width[tx_id] for tx_id in in_order)
        wanted.append(next(tx_id for tx_id in in_order if width[tx_id] == maximum))
    longest_frame = longest.to_frame()
    assert list(longest.ids) == wanted
    assert list(longest.gene_ids) == genes
    assert {tx_id: ranges_of(longest_frame, tx_id) for tx_id in wanted} == {
        tx_id: ranges[tx_id] for tx_id in wanted
    }

    union_frame = union.to_frame()
    assert list(union.ids) == genes
    assert list(union.gene_ids) == genes
    for gene_id in genes:
        ascending = merged(r for tx_id in transcripts[gene_id] for r in ranges[tx_id])
        expected = ascending if strand[gene_id] == 1 else ascending[::-1]
        assert ranges_of(union_frame, gene_id) == expected, gene_id


def test_union_merges_a_range_inside_an_earlier_range_and_keeps_the_genes_apart():
    rows = hand_rows(
        ("T1", "GA", 1, None, ((1, 1000),)),
        ("T2", "GB", -1, None, ((400, 500), (40, 50), (5, 10))),
        ("T3", "GB", -1, None, ((302, 310), (20, 300))),
        ("T4", "GB", -1, None, ((60, 70),)),
    )
    got = derive_features(rows, FeatureType.EXONS, PerGene.UNION)
    frame = got.to_frame()
    assert got.ids == ("GA", "GB")
    assert ranges_of(frame, "GA") == [(1, 1000)]
    assert ranges_of(frame, "GB") == [(400, 500), (302, 310), (20, 300), (5, 10)]


def test_longest_tie_goes_to_the_lower_start_before_the_lower_tx_row():
    rows = hand_rows(
        ("T1", "G", 1, None, ((50, 59),)),
        ("T2", "G", 1, None, ((10, 19),)),
        ("T3", "G", 1, None, ((5, 13),)),
    )
    assert derive_features(rows, FeatureType.EXONS, PerGene.LONGEST).ids == ("T2",)


def test_transcript_without_a_3utr_has_no_range():
    with fixture_reference() as ref:
        got = ref.features(FeatureType.UTR3, ALL)
    counts = dict(zip(got.ids, np.diff(got.offsets).tolist(), strict=True))
    widths = dict(zip(got.ids, got.widths().tolist(), strict=True))
    assert counts[T07] == counts[T12] == 0
    assert widths[T07] == widths[T12] == 0
    assert all(count > 0 for tx_id, count in counts.items() if tx_id not in {T07, T12})
    assert not {T07, T12} & set(got.to_frame()["id"])


@pytest.mark.parametrize(("per_gene", "name"), [(PerGene.LONGEST, T12), (PerGene.UNION, G08)])
def test_a_gene_whose_transcripts_lack_the_feature_has_a_feature_without_ranges(per_gene, name):
    with fixture_reference() as ref:
        got = ref.features(FeatureType.UTR3, ALL, per_gene=per_gene)
    position = got.gene_ids.index(G08)
    assert got.ids[position] == name
    assert got.offsets[position] == got.offsets[position + 1]


def test_widths_are_the_lengths_of_the_r_sequences():
    expected = golden("annotation/utr3_min0_sequences.tsv")
    with fixture_reference() as ref:
        got = ref.features(FeatureType.UTR3, ALL)
    assert list(got.ids) == expected["name_tx"].tolist()
    assert got.widths().tolist() == expected["width"].astype(int).tolist()


@pytest.mark.parametrize(
    ("feature", "expected"),
    [
        (FeatureType.UTR5, {"T1": [(10, 20)], "T2": [(50, 60)], "T3": [(100, 119)]}),
        (FeatureType.CDS, {"T1": [(30, 40)], "T2": [(30, 40)], "T3": [(120, 180)]}),
        (FeatureType.UTR3, {"T1": [(50, 60)], "T2": [(10, 20)], "T3": [(181, 200)]}),
        (
            FeatureType.EXONS,
            {
                "T1": [(10, 20), (30, 40), (50, 60)],
                "T2": [(50, 60), (30, 40), (10, 20)],
                "T3": [(100, 200)],
                "T4": [(300, 400)],
            },
        ),
    ],
)
def test_cut_keeps_the_part_of_each_exon_in_the_feature(feature, expected):
    got = derive_features(EDGES, feature, PerGene.EACH_TRANSCRIPT)
    frame = got.to_frame()
    assert list(got.ids) == list(expected)
    assert {tx_id: ranges_of(frame, tx_id) for tx_id in got.ids} == expected


def test_cut_features_skip_a_transcript_without_a_cds():
    with fixture_reference() as ref:
        everything = ref.annotation.exon_rows(NOT_CODING, coding_only=False)
        coding = ref.annotation.exon_rows(NOT_CODING, coding_only=True)
    assert T05 in set(everything.tx_id)
    for feature in (FeatureType.UTR3, FeatureType.UTR5, FeatureType.CDS):
        for per_gene in PerGene:
            got = derive_features(everything, feature, per_gene).to_frame()
            wanted = derive_features(coding, feature, per_gene).to_frame()
            assert got.equals(wanted), (feature, per_gene)
            assert T05 not in set(got["id"])


def test_exon_rows_come_in_transcript_order_with_the_tables_of_the_ensdb():
    t03 = TranscriptFilter(canonical=False, tx_ids=(T03,))
    with fixture_reference() as ref:
        rows = ref.annotation.exon_rows(t03, coding_only=True)
    assert rows.tx_id.tolist() == [T03] * 4
    assert rows.gene_id.tolist() == [G02] * 4
    assert rows.seq_name.tolist() == ["1"] * 4
    assert rows.strand.tolist() == [-1] * 4
    assert rows.exon_start.tolist() == [29000, 27000, 25000, 22000]
    assert rows.exon_end.tolist() == [29200, 27400, 25600, 22800]
    assert rows.cds_start.tolist() == [25599] * 4
    assert rows.cds_end.tolist() == [27201] * 4
    assert rows.tx_seq_start.tolist() == [22000] * 4
    assert rows.tx_row.tolist() == [3] * 4
    dtypes = {name: str(getattr(rows, name).dtype) for name in ExonRows.__dataclass_fields__}
    assert dtypes == {
        "tx_id": "object",
        "gene_id": "object",
        "seq_name": "object",
        "strand": "int8",
        **dict.fromkeys(
            (
                "cds_start",
                "cds_end",
                "exon_start",
                "exon_end",
                "tx_seq_start",
                "tx_row",
            ),
            "int64",
        ),
    }


def test_exon_rows_sort_on_tx_id_and_coding_only_keeps_the_transcripts_with_a_cds():
    with fixture_reference() as ref:
        listed = ref.transcripts(EVERYTHING)
        everything = ref.annotation.exon_rows(EVERYTHING, coding_only=False)
        coding = ref.annotation.exon_rows(EVERYTHING, coding_only=True)
    with_cds = set(listed["tx_id"][listed["cds_start"].notna()])
    assert list(dict.fromkeys(everything.tx_id)) == sorted(listed["tx_id"])
    assert list(dict.fromkeys(coding.tx_id)) == sorted(with_cds)
    assert T05 not in with_cds
    assert everything.tx_id.size == len(golden("annotation/exons_everything_features.tsv"))


def test_exon_rows_of_many_chunks_equal_the_rows_of_one_chunk(monkeypatch):
    with fixture_reference() as ref:
        whole = ref.annotation.exon_rows(EVERYTHING, coding_only=False)
        monkeypatch.setattr(seedmatch.annotation.ensdb, "_CHUNK_ROWS", 7)
        chunked = ref.annotation.exon_rows(EVERYTHING, coding_only=False)
    for name in ExonRows.__dataclass_fields__:
        assert getattr(chunked, name).tolist() == getattr(whole, name).tolist(), name
        assert getattr(chunked, name).dtype == getattr(whole, name).dtype, name


def test_exon_rows_with_a_strand_that_is_not_plus_or_minus_raise(tmp_path):
    path = copy_ensdb(
        tmp_path, "UPDATE gene SET seq_strand = 0 WHERE gene_id = 'ENSSYNG00000000008'"
    )
    db = AnnotationDb.open(path)
    message = f"{os.fspath(path)!r}: transcript 'ENSSYNT00000000012' has strand 0; expected 1 or -1"
    with pytest.raises(InputError, match=re.escape(message)):
        db.exon_rows(EVERYTHING, coding_only=False)
    db.close()


@pytest.mark.parametrize("per_gene", PerGene)
@pytest.mark.parametrize("feature", FeatureType)
def test_a_filter_without_a_match_gives_no_feature(feature, per_gene):
    with fixture_reference() as ref:
        got = ref.features(feature, TranscriptFilter(gene_ids=("none",)), per_gene=per_gene)
    assert len(got) == 0
    assert got.offsets.tolist() == [0]
    assert got.widths().tolist() == []
    frame = got.to_frame()
    assert len(frame) == 0
    assert list(frame.columns) == ["id", "gene_id", "contig", "strand", "rank", "start", "end"]


def test_to_frame_gives_one_row_for_each_range():
    with fixture_reference() as ref:
        got = ref.features(FeatureType.UTR3, ALL, per_gene=PerGene.UNION)
    frame = got.to_frame()
    assert {name: str(dtype) for name, dtype in frame.dtypes.items()} == {
        "id": "object",
        "gene_id": "object",
        "contig": "object",
        "strand": "object",
        "rank": "int64",
        "start": "int64",
        "end": "int64",
    }
    assert len(frame) == got.offsets[-1]
    assert (frame["id"] == frame["gene_id"]).all()
    assert frame[frame["id"] == G02]["rank"].tolist() == [1, 2, 3]
    assert list(frame.index) == list(range(len(frame)))


def feature_parts(**changes):
    parts = {
        "ids": ("A", "B"),
        "gene_ids": ("GA", "GB"),
        "contigs": ("1", "1"),
        "strands": np.array([1, -1], dtype=np.int8),
        "starts": np.array([10, 20, 40], dtype=np.int64),
        "ends": np.array([15, 30, 50], dtype=np.int64),
        "offsets": np.array([0, 1, 3], dtype=np.int64),
    }
    return {**parts, **changes}


def test_features_with_parts_that_fit_together_have_their_ranges():
    got = Features(**feature_parts())
    assert got.widths().tolist() == [6, 22]
    assert ranges_of(got.to_frame(), "B") == [(20, 30), (40, 50)]


# B has no range, C is on the minus strand with its ranges in transcript order, and D has a range
# of exactly 8 bases.
AT_LEAST_PARTS = {
    "ids": ("A", "B", "C", "D"),
    "gene_ids": ("GA", "GB", "GC", "GD"),
    "contigs": ("1", "2", "1", "3"),
    "strands": np.array([1, 1, -1, -1], dtype=np.int8),
    "starts": np.array([10, 40, 20, 100], dtype=np.int64),
    "ends": np.array([15, 50, 30, 107], dtype=np.int64),
    "offsets": np.array([0, 1, 1, 3, 4], dtype=np.int64),
}


@pytest.mark.parametrize(
    ("min_width", "ids", "starts", "ends", "offsets"),
    [
        (0, ("A", "B", "C", "D"), [10, 40, 20, 100], [15, 50, 30, 107], [0, 1, 1, 3, 4]),
        (7, ("C", "D"), [40, 20, 100], [50, 30, 107], [0, 2, 3]),
        (8, ("C", "D"), [40, 20, 100], [50, 30, 107], [0, 2, 3]),
        (9, ("C",), [40, 20], [50, 30], [0, 2]),
        (23, (), [], [], [0]),
    ],
)
def test_at_least_keeps_the_features_of_that_width_with_their_ranges(
    min_width, ids, starts, ends, offsets
):
    features = Features(**AT_LEAST_PARTS)
    got = features.at_least(min_width)
    kept = [AT_LEAST_PARTS["ids"].index(name) for name in ids]
    assert got.ids == ids
    assert got.gene_ids == tuple(AT_LEAST_PARTS["gene_ids"][i] for i in kept)
    assert got.contigs == tuple(AT_LEAST_PARTS["contigs"][i] for i in kept)
    assert got.strands.tolist() == [int(AT_LEAST_PARTS["strands"][i]) for i in kept]
    assert got.starts.tolist() == starts
    assert got.ends.tolist() == ends
    assert got.offsets.tolist() == offsets
    assert got.widths().tolist() == [int(features.widths()[i]) for i in kept]


def test_at_least_gives_arrays_that_are_read_only():
    got = Features(**AT_LEAST_PARTS).at_least(8)
    for array in (got.strands, got.starts, got.ends, got.offsets):
        with pytest.raises(ValueError, match="read-only"):
            array[0] = 0


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"gene_ids": ("GA",)}, "gene_ids has 1 entries; expected one for each of the 2 features"),
        ({"contigs": ("1", "1", "X")}, "contigs has 3 entries; expected one for each of the 2"),
        ({"strands": np.array([1, -1], dtype=np.int64)}, "strands has dtype int64; expected int8"),
        ({"starts": np.array([10, 20, 40], dtype=np.int32)}, "starts has dtype int32; expected"),
        (
            {"ends": np.array([[15, 30, 50]], dtype=np.int64)},
            "ends has 2 dimensions; expected a 1-D",
        ),
        (
            {"offsets": np.array([0, 1, 3], dtype=np.int32)},
            "offsets has dtype int32; expected int64",
        ),
        (
            {"strands": np.array([1], dtype=np.int8)},
            "strands has 1 entries; expected one for each of the 2 features",
        ),
        (
            {"ends": np.array([15, 30], dtype=np.int64)},
            "ends has 2 entries; expected one for each of the 3 ranges",
        ),
        ({"offsets": np.array([0, 3], dtype=np.int64)}, "offsets has 2 entries; expected 3 for 2"),
        ({"offsets": np.array([1, 2, 3], dtype=np.int64)}, "offsets must start at 0; got 1"),
        (
            {"offsets": np.array([0, 1, 2], dtype=np.int64)},
            "offsets must end at 3, the number of ranges; got 2",
        ),
        (
            {"offsets": np.array([0, 4, 3], dtype=np.int64)},
            "offsets[2] 3 is less than offsets[1] 4; expected offsets that do not decrease",
        ),
        ({"strands": np.array([1, 0], dtype=np.int8)}, "strands[1] is 0; expected 1 or -1"),
        ({"starts": np.array([0, 20, 40], dtype=np.int64)}, "starts[0] is 0; expected 1 or more"),
        (
            {"ends": np.array([15, 19, 50], dtype=np.int64)},
            "range 1 has start 20 after end 19; expected a start that is not after the end",
        ),
    ],
)
def test_features_reject_parts_that_do_not_fit_together(change, message):
    with pytest.raises(InputError, match=re.escape(message)):
        Features(**feature_parts(**change))


def test_features_arrays_are_read_only():
    with fixture_reference() as ref:
        got = ref.features(FeatureType.UTR3, ALL)
    assert isinstance(got, Features)
    for array in (got.strands, got.starts, got.ends, got.offsets):
        with pytest.raises(ValueError, match="read-only"):
            array[0] = 0


def test_reference_features_defaults_are_the_defaults_of_load_annotations():
    expected = golden("annotation/utr3_features.tsv")
    with fixture_reference() as ref:
        got = ref.features(FeatureType.UTR3)
    assert got.to_frame()["id"].tolist() == expected["name"].tolist()


def test_features_reject_a_feature_of_another_type():
    message = "feature '3UTR' has type str; expected a FeatureType"
    with fixture_reference() as ref, pytest.raises(InputError, match=re.escape(message)):
        ref.features("3UTR")


def test_features_reject_a_per_gene_of_another_type():
    message = "per_gene 'union' has type str; expected a PerGene"
    with fixture_reference() as ref, pytest.raises(InputError, match=re.escape(message)):
        ref.features(FeatureType.UTR3, per_gene="union")


def test_exon_rows_give_0_as_the_cds_of_a_transcript_without_a_cds():
    db = AnnotationDb.open(ENSDB)
    rows = db.exon_rows(
        TranscriptFilter(canonical=False, tx_ids=(T05,), protein_coding=False), coding_only=False
    )
    db.close()
    assert rows.tx_id.tolist() == [T05, T05]
    assert rows.cds_start.tolist() == [0, 0]
    assert rows.cds_end.tolist() == [0, 0]
