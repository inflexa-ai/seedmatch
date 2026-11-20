import re

import numpy as np
import pandas as pd
import pytest
from helpers import GUIDE, NOT_A_DATA_FRAME, golden_rows, sirna_de, sirna_seqs

import seedmatch
from seedmatch import Guide, SearchOptions, Seed, SequenceSet
from seedmatch.errors import InputError
from seedmatch.report import annotate_seed_counts

MER7M8 = Guide(GUIDE).target(Seed.MER7M8)


def golden_counts(seed):
    return [
        [row.gene_id, int(row.count)]
        for row in golden_rows("seedmatchr/deseq2_mode.tsv")
        if row.seed == seed
    ]


def seqs_with(*records, gene_ids=None):
    return SequenceSet.from_records(records, gene_ids=gene_ids)


@pytest.mark.parametrize("seed", ["mer7m8", "mer8", "mer6", "mer7A1", "Full"])
def test_seed_counts_match_r(seed):
    target = Guide(GUIDE).target(Seed(seed))
    got = annotate_seed_counts(sirna_de(), sirna_seqs(), target)
    assert got[["gene_id", seed]].to_numpy().tolist() == golden_counts(seed)


def test_gene_count_is_the_maximum_over_its_sequences_D003():
    seqs = SequenceSet.from_records(
        [("t1", "CTCTATA"), ("t2", "CTCTATACTCTATACTCTATA")], gene_ids=["g", "g"]
    )
    de = pd.DataFrame({"gene_id": ["g"]})
    assert annotate_seed_counts(de, seqs, MER7M8)["mer7m8"].tolist() == [3]


@pytest.mark.parametrize("order", [(0, 1), (1, 0)])
def test_the_maximum_does_not_depend_on_the_order_of_the_sequences_D003(order):
    records = [("t1", "CTCTATA"), ("t2", "CTCTATACTCTATA")]
    seqs = seqs_with(*(records[i] for i in order), gene_ids=["g", "g"])
    de = pd.DataFrame({"gene_id": ["g"]})
    assert annotate_seed_counts(de, seqs, MER7M8)["mer7m8"].tolist() == [2]


def test_existing_column_is_replaced_not_duplicated_D015():
    got = annotate_seed_counts(sirna_de(), sirna_seqs(), MER7M8)
    assert list(got.columns).count("mer7m8") == 1


def test_replaced_column_keeps_its_position_and_gets_the_new_values_D015():
    de = pd.DataFrame({"gene_id": ["g", "h"], "mer7m8": [99, 99], "other": ["x", "y"]})
    seqs = seqs_with(("g", "CTCTATA"), ("h", "AAAA"))
    got = annotate_seed_counts(de, seqs, MER7M8)
    assert got.columns.tolist() == ["gene_id", "mer7m8", "other"]
    assert got["mer7m8"].tolist() == [1, 0]
    assert de["mer7m8"].tolist() == [99, 99]


def test_new_column_is_the_last_column_with_the_target_name_and_the_dtype_int64():
    de = pd.DataFrame({"gene_id": ["g", "h"], "padj": [0.5, 0.25]})
    got = annotate_seed_counts(de, seqs_with(("g", "CTCTATA"), ("h", "AAAA")), MER7M8)
    assert got.columns.tolist() == ["gene_id", "padj", "mer7m8"]
    assert got["mer7m8"].dtype == np.int64


def test_column_argument_names_the_count_column():
    de = pd.DataFrame({"gene_id": ["g"]})
    got = annotate_seed_counts(de, seqs_with(("g", "CTCTATA")), MER7M8, column="hits")
    assert got.columns.tolist() == ["gene_id", "hits"]
    assert got["hits"].tolist() == [1]


def test_custom_target_name_is_the_default_column_name():
    target = Guide(GUIDE).target(seedmatch.CustomSeed("g2-g8", 2, 8))
    got = annotate_seed_counts(pd.DataFrame({"gene_id": ["g"]}), seqs_with(("g", "AAAA")), target)
    assert got.columns.tolist() == ["gene_id", "g2-g8"]


def test_gene_without_sequence_counts_zero_when_not_shared():
    de = pd.DataFrame({"gene_id": ["absent"]})
    got = annotate_seed_counts(de, sirna_seqs(), MER7M8, shared_genes=False)
    assert got["mer7m8"].tolist() == [0]


def test_not_shared_keeps_each_row_and_the_index():
    de = pd.DataFrame({"gene_id": ["absent", "g", "absent"]}, index=[7, 3, 5])
    got = annotate_seed_counts(de, seqs_with(("g", "CTCTATA")), MER7M8, shared_genes=False)
    assert got.index.tolist() == [7, 3, 5]
    assert got["mer7m8"].tolist() == [0, 1, 0]


def test_shared_genes_keeps_the_rows_with_a_sequence_in_order_and_resets_the_index():
    de = pd.DataFrame({"gene_id": ["absent", "h", "absent", "g", "h"]}, index=[7, 3, 5, 9, 1])
    seqs = seqs_with(("g", "CTCTATA"), ("h", "AAAA"))
    got = annotate_seed_counts(de, seqs, MER7M8)
    assert got["gene_id"].tolist() == ["h", "g", "h"]
    assert got["mer7m8"].tolist() == [0, 1, 0]
    assert got.index.tolist() == [0, 1, 2]


def de_with_absent_gene():
    de = sirna_de()
    absent = de.iloc[[0]].assign(gene_id="ENSRNOG_NO_SEQUENCE")
    return pd.concat([de.iloc[:50], absent, de.iloc[50:]], ignore_index=True)


def golden_absent_gene(shared_genes):
    return [
        [row.gene_id, int(row.count)]
        for row in golden_rows("seedmatchr/deseq2_mode_absent_gene.tsv")
        if row.shared_genes == shared_genes
    ]


def test_shared_genes_drops_the_gene_without_a_sequence_as_r_does():
    got = annotate_seed_counts(de_with_absent_gene(), sirna_seqs(), MER7M8)
    expected = golden_absent_gene("TRUE")
    assert len(expected) == len(sirna_de())
    assert got[["gene_id", "mer7m8"]].to_numpy().tolist() == expected


def test_not_shared_fills_the_gene_without_a_sequence_with_zero_as_r_does():
    got = annotate_seed_counts(de_with_absent_gene(), sirna_seqs(), MER7M8, shared_genes=False)
    expected = golden_absent_gene("FALSE")
    assert ["ENSRNOG_NO_SEQUENCE", 0] in expected
    assert got[["gene_id", "mer7m8"]].to_numpy().tolist() == expected


def test_input_is_not_changed():
    de = sirna_de()
    before = de.copy()
    annotate_seed_counts(de, sirna_seqs(), MER7M8)
    annotate_seed_counts(de, sirna_seqs(), MER7M8, shared_genes=False, column="hits")
    pd.testing.assert_frame_equal(de, before)


def test_options_control_the_search():
    seqs = seqs_with(("g", "CTCTAGA"))
    de = pd.DataFrame({"gene_id": ["g"]})
    assert annotate_seed_counts(de, seqs, MER7M8)["mer7m8"].tolist() == [0]
    one_mismatch = SearchOptions(max_mismatch=1)
    assert annotate_seed_counts(de, seqs, MER7M8, options=one_mismatch)["mer7m8"].tolist() == [1]


def test_a_de_table_without_gene_id_names_the_column():
    message = "DE table has no 'gene_id' column; found ['symbol']"
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        annotate_seed_counts(pd.DataFrame({"symbol": ["Ttr"]}), sirna_seqs(), MER7M8)


@pytest.mark.parametrize(("de", "type_name"), NOT_A_DATA_FRAME)
def test_a_table_that_is_not_a_data_frame_is_an_error(de, type_name):
    message = f"de has type {type_name}; expected a pandas.DataFrame"
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        annotate_seed_counts(de, sirna_seqs(), MER7M8)


def test_empty_sequence_set_gives_zero_counts_or_no_rows():
    de = pd.DataFrame({"gene_id": ["g"]})
    empty = SequenceSet.from_records([])
    assert annotate_seed_counts(de, empty, MER7M8, shared_genes=False)["mer7m8"].tolist() == [0]
    got = annotate_seed_counts(de, empty, MER7M8)
    assert got.empty
    assert got.columns.tolist() == ["gene_id", "mer7m8"]
    assert got["mer7m8"].dtype == np.int64


def test_empty_de_table_gives_an_empty_int64_column():
    got = annotate_seed_counts(pd.DataFrame({"gene_id": []}), sirna_seqs(), MER7M8)
    assert got.empty
    assert got["mer7m8"].dtype == np.int64


def test_empty_sequence_and_short_sequence_count_zero():
    seqs = seqs_with(("g", ""), ("h", "CTC"))
    de = pd.DataFrame({"gene_id": ["g", "h"]})
    assert annotate_seed_counts(de, seqs, MER7M8)["mer7m8"].tolist() == [0, 0]


def test_a_gene_in_two_rows_gets_the_same_count_in_both():
    de = pd.DataFrame({"gene_id": ["g", "g"]})
    got = annotate_seed_counts(de, seqs_with(("g", "CTCTATA")), MER7M8)
    assert got["mer7m8"].tolist() == [1, 1]


def test_top_level_package_exports_annotate_seed_counts():
    assert seedmatch.annotate_seed_counts is annotate_seed_counts


@pytest.mark.parametrize("dtype", ["object", "category", "string"])
def test_text_dtypes_of_gene_id_give_the_same_counts(dtype):
    de = pd.DataFrame({"gene_id": pd.Series(["g", "absent", "h", "g"], dtype=dtype)})
    seqs = seqs_with(("g", "CTCTATA"), ("h", "AAAA"))
    assert annotate_seed_counts(de, seqs, MER7M8, shared_genes=False)["mer7m8"].tolist() == [
        1,
        0,
        0,
        1,
    ]
    assert annotate_seed_counts(de, seqs, MER7M8)["gene_id"].tolist() == ["g", "h", "g"]


def test_a_missing_gene_id_is_permitted_and_counts_zero():
    de = pd.DataFrame({"gene_id": ["g", None, float("nan")]})
    got = annotate_seed_counts(de, seqs_with(("g", "CTCTATA")), MER7M8, shared_genes=False)
    assert got["mer7m8"].tolist() == [1, 0, 0]


@pytest.mark.parametrize(
    ("values", "dtype"),
    [
        ([1, 2], "int64"),
        ([1.5, 2.5], "float64"),
        ([1, "g"], "object"),
    ],
)
def test_gene_id_that_is_not_text_names_the_dtype_and_the_conversion(values, dtype):
    de = pd.DataFrame({"gene_id": pd.Series(values, dtype=dtype)})
    with pytest.raises(InputError, match=rf"'gene_id' has dtype {dtype}.*astype\(str\)"):
        annotate_seed_counts(de, seqs_with(("1", "CTCTATA")), MER7M8)


def test_gene_id_with_numeric_categories_is_not_text():
    de = pd.DataFrame({"gene_id": pd.Categorical([1, 2])})
    with pytest.raises(InputError, match=r"'gene_id' has dtype category.*astype\(str\)"):
        annotate_seed_counts(de, seqs_with(("1", "CTCTATA")), MER7M8)


def test_the_column_argument_cannot_be_gene_id():
    de = pd.DataFrame({"gene_id": ["g"]})
    with pytest.raises(InputError, match=r"column 'gene_id' would overwrite the gene IDs"):
        annotate_seed_counts(de, seqs_with(("g", "CTCTATA")), MER7M8, column="gene_id")


def test_a_target_named_gene_id_cannot_name_the_count_column():
    de = pd.DataFrame({"gene_id": ["g"]})
    target = Guide(GUIDE).target(seedmatch.CustomSeed("gene_id", 2, 8))
    with pytest.raises(InputError, match=r"column 'gene_id' would overwrite the gene IDs"):
        annotate_seed_counts(de, seqs_with(("g", "CTCTATA")), target)


def test_annotate_rejects_a_shared_genes_flag_that_is_not_a_bool():
    with pytest.raises(InputError, match="shared_genes 'no' has type str; expected a bool"):
        annotate_seed_counts(sirna_de(), sirna_seqs(), MER7M8, shared_genes="no")


def test_annotate_rejects_a_column_name_that_is_not_text():
    with pytest.raises(InputError, match="column 5 has type int; expected a str"):
        annotate_seed_counts(sirna_de(), sirna_seqs(), MER7M8, column=5)


def test_annotate_rejects_a_wobble_target_with_a_fixed_pattern_D026():
    wobble = Guide(GUIDE).target(Seed.MER7M8, wobble=True)
    message = "target 'mer7m8' has the ambiguity letter 'R'; expected fixed_pattern=False"
    with pytest.raises(InputError, match=message):
        annotate_seed_counts(sirna_de(), sirna_seqs(), wobble)


def test_annotate_counts_a_wobble_target_with_a_pattern_that_is_not_fixed_D026():
    wobble = Guide(GUIDE).target(Seed.MER7M8, wobble=True)
    options = SearchOptions(indels=True, fixed_pattern=False)
    got = annotate_seed_counts(sirna_de(), sirna_seqs(), wobble, options=options)["mer7m8"]
    plain = annotate_seed_counts(sirna_de(), sirna_seqs(), MER7M8)["mer7m8"]
    assert (got.to_numpy() >= plain.to_numpy()).all()
    assert got.sum() > plain.sum()
