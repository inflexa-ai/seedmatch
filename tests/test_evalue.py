import contextlib
import inspect

import numpy as np
import pandas as pd
import pytest
from helpers import GUIDE, NOT_A_DATA_FRAME, TEXT_DTYPE, exactly, golden, sirna_seqs

from seedmatch import (
    Guide,
    SearchOptions,
    Seed,
    SequenceSet,
    annotate_seed_counts,
    e_value_null,
    e_values,
)
from seedmatch import evalue as evalue_module
from seedmatch.errors import InputError
from seedmatch.evalue import _shuffled_guides

EVALUE = golden("analyses/evalue.tsv")
NULL = golden("analyses/evalue_null.tsv")
OBSERVED = golden("analyses/evalue_observed.tsv")
SHUFFLES = NULL["shuffle"].astype(int).max()
# The 4 counts of the gene a and of the gene b in the notes of evalue.tsv.
NOTE_NULL = {"a": [1, 2, 3, 0], "b": [0, 1, 0, 0]}


def make_null(rows):
    """Return a null from a dict of gene ID to counts, or from a list of (gene ID, counts) pairs."""
    pairs = list(rows.items()) if isinstance(rows, dict) else rows
    return pd.DataFrame(
        [counts for _, counts in pairs],
        index=pd.Index([gene for gene, _ in pairs], name="gene_id"),
        columns=range(len(pairs[0][1])),
    )


def make_de(genes, counts, column="mer7m8"):
    return pd.DataFrame({"gene_id": genes, column: counts})


def r_null():
    """Return the null of the golden file, with the rows in the R order (g40 to g01)."""
    counts = NULL["count"].astype(np.int32).to_numpy().reshape(-1, SHUFFLES)
    genes = NULL["gene_id"].to_numpy()[::SHUFFLES]
    return pd.DataFrame(counts, index=pd.Index(genes, name="gene_id"), columns=range(SHUFFLES))


def r_observed():
    return OBSERVED.assign(overall_matches=OBSERVED["overall_matches"].astype(int))


def multi_sequence_set():
    """Return 21 random sequences that share 7 genes, so that a gene has sequences apart."""
    letters = np.random.default_rng(5).choice(list("ACGT"), size=(21, 1500))
    records = [(f"t{i}", "".join(row)) for i, row in enumerate(letters)]
    return SequenceSet.from_records(records, gene_ids=[f"g{(i * 5) % 7}" for i in range(21)])


def oracle_shuffles(guide, count, seed):
    """Draw the shuffles with the rule of D-008, as a reader of the docstring does."""
    generator = np.random.default_rng(seed)
    letters = guide.sequence
    return [
        Guide("".join(letters[i] for i in generator.permutation(len(letters))))
        for _ in range(count)
    ]


def test_e_values_match_r():
    de = r_observed()
    result = e_values(de, r_null(), column="overall_matches")
    assert result["gene_id"].tolist() == EVALUE["gene_id"].tolist()
    assert result["overall_matches"].tolist() == EVALUE["overall_matches"].astype(int).tolist()
    assert result["e_value"].tolist() == EVALUE["e_value"].astype(float).tolist()
    assert result["e_value"].dtype == np.float64


def test_e_values_match_genes_by_gene_id_D008():
    de = r_observed()
    null = r_null()
    expected = e_values(de, null, column="overall_matches")
    order = np.random.default_rng(1).permutation(len(null))
    shuffled = e_values(de, null.iloc[order], column="overall_matches")
    pd.testing.assert_frame_equal(shuffled, expected)
    # The rows of the null are in the reverse order of the genes, so a match by position differs.
    assert null.index.tolist() == sorted(null.index, reverse=True)


def test_e_value_null_counts_each_shuffle_like_annotate_seed_counts_D008():
    sequences = multi_sequence_set()
    guide = Guide(GUIDE)
    options = SearchOptions(max_mismatch=1)
    null = e_value_null(guide, sequences, rng=7, n_random=5, options=options)
    genes = list(dict.fromkeys(sequences.gene_ids))
    table = pd.DataFrame({"gene_id": genes})
    expected = {}
    for j, shuffle in enumerate(oracle_shuffles(guide, 5, 7)):
        annotated = annotate_seed_counts(
            table,
            sequences,
            shuffle.target(Seed.MER7M8),
            column="count",
            options=options,
            shared_genes=False,
        )
        expected[j] = annotated["count"].to_numpy()
    assert null.index.tolist() == genes
    assert null.index.name == "gene_id"
    assert null.columns.tolist() == [0, 1, 2, 3, 4]
    assert (null.dtypes == np.int32).all()
    assert null.to_numpy().max() > 0
    for j, counts in expected.items():
        assert null[j].tolist() == counts.tolist()


def test_e_value_null_uses_the_default_options_of_annotate_seed_counts_D008():
    sequences = sirna_seqs()
    guide = Guide(GUIDE)
    null = e_value_null(guide, sequences, rng=3, n_random=4)
    expected = {
        j: annotate_seed_counts(
            pd.DataFrame({"gene_id": list(dict.fromkeys(sequences.gene_ids))}),
            sequences,
            shuffle.target(Seed.MER7M8),
            column="count",
            shared_genes=False,
        )["count"].tolist()
        for j, shuffle in enumerate(oracle_shuffles(guide, 4, 3))
    }
    assert {j: null[j].tolist() for j in range(4)} == expected
    # An exact search gives the same counts for some other defaults, so the defaults need a test.
    assert (
        inspect.signature(e_value_null).parameters["options"].default
        == inspect.signature(annotate_seed_counts).parameters["options"].default
    )


def test_e_value_null_uses_the_seed_that_the_caller_gives():
    sequences = multi_sequence_set()
    guide = Guide(GUIDE)
    options = SearchOptions(max_mismatch=1)
    null = e_value_null(guide, sequences, rng=2, n_random=3, seed=Seed.MER8, options=options)
    expected = [
        annotate_seed_counts(
            pd.DataFrame({"gene_id": list(dict.fromkeys(sequences.gene_ids))}),
            sequences,
            shuffle.target(Seed.MER8),
            column="count",
            options=options,
            shared_genes=False,
        )["count"].tolist()
        for shuffle in oracle_shuffles(guide, 3, 2)
    ]
    assert [null[j].tolist() for j in range(3)] == expected


def test_each_shuffle_keeps_the_letters_of_the_guide_D008():
    guide = Guide(GUIDE)
    shuffles = list(_shuffled_guides(guide, 50, np.random.default_rng(3)))
    assert len(shuffles) == 50
    assert all(sorted(shuffle.sequence) == sorted(guide.sequence) for shuffle in shuffles)
    assert any(shuffle.sequence != guide.sequence for shuffle in shuffles)


def test_the_shuffles_are_the_permutations_of_the_generator_in_order_D008():
    guide = Guide(GUIDE)
    shuffles = list(_shuffled_guides(guide, 6, np.random.default_rng(11)))
    assert shuffles == oracle_shuffles(guide, 6, 11)


def test_e_value_null_is_the_same_for_the_same_seed():
    sequences = multi_sequence_set()
    guide = Guide(GUIDE)
    options = SearchOptions(max_mismatch=1)
    first = e_value_null(guide, sequences, rng=7, n_random=6, options=options)
    again = e_value_null(guide, sequences, rng=7, n_random=6, options=options)
    from_generator = e_value_null(
        guide, sequences, rng=np.random.default_rng(7), n_random=6, options=options
    )
    other = e_value_null(guide, sequences, rng=8, n_random=6, options=options)
    pd.testing.assert_frame_equal(again, first)
    pd.testing.assert_frame_equal(from_generator, first)
    assert not other.equals(first)


def test_e_value_null_does_not_depend_on_workers():
    sequences = multi_sequence_set()
    guide = Guide(GUIDE)
    options = SearchOptions(max_mismatch=1)
    one = e_value_null(guide, sequences, rng=7, n_random=40, options=options, workers=1)
    four = e_value_null(guide, sequences, rng=7, n_random=40, options=options, workers=4)
    pd.testing.assert_frame_equal(four, one)


def test_e_value_null_draws_from_the_generator_of_the_caller():
    generator = np.random.default_rng(7)
    before = generator.bit_generator.state
    e_value_null(Guide(GUIDE), multi_sequence_set(), rng=generator, n_random=3)
    assert generator.bit_generator.state != before


def test_e_value_null_of_a_set_without_sequences_has_no_row():
    empty = SequenceSet.from_records([])
    null = e_value_null(Guide(GUIDE), empty, rng=1, n_random=3)
    assert null.shape == (0, 3)
    assert null.index.name == "gene_id"
    assert str(null.index.dtype) == TEXT_DTYPE
    assert (null.dtypes == np.int32).all()


def test_e_value_null_gives_one_int32_block():
    null = e_value_null(Guide(GUIDE), sirna_seqs(), rng=1, n_random=4)
    assert np.shares_memory(null.to_numpy(), null[0].to_numpy())
    assert null.to_numpy().dtype == np.int32


def test_e_value_null_gives_the_frame_the_memory_of_the_table_under_copy_on_write(monkeypatch):
    """pandas 3 always has copy-on-write, and then a frame made without `copy=False` copies."""
    sequences = sirna_seqs()
    table = np.empty((4, len(set(sequences.gene_ids))), dtype=np.int32)
    monkeypatch.setattr(evalue_module, "_result_table", lambda *_: table)
    # pandas 2 has copy-on-write only with the option, and pandas 3 deprecates the option.
    copy_on_write = (
        pd.option_context("mode.copy_on_write", True)  # noqa: FBT003 - the API takes pairs
        if pd.__version__.startswith("2.")
        else contextlib.nullcontext()
    )
    with copy_on_write:
        null = e_value_null(Guide(GUIDE), sequences, rng=1, n_random=4)
    assert np.shares_memory(null.to_numpy(), table)


def test_e_value_null_gives_each_gene_one_row_for_its_sequences_apart():
    sequences = SequenceSet.from_records(
        [("s1", "ACGTACGT"), ("s2", "TTTT"), ("s3", "ACGTACGT")], gene_ids=["x", "y", "x"]
    )
    null = e_value_null(Guide(GUIDE), sequences, rng=1, n_random=2)
    assert null.index.tolist() == ["x", "y"]


def test_e_value_null_rejects_a_seed_that_does_not_fit_the_guide_before_any_draw():
    generator = np.random.default_rng(1)
    before = generator.bit_generator.state
    with pytest.raises(InputError, match=r"^18mer needs a guide of at least 19 nt, got 8$"):
        e_value_null(Guide("ACGUACGU"), sirna_seqs(), rng=generator, seed=Seed.MER18)
    assert generator.bit_generator.state == before


def test_e_value_null_allocates_the_result_before_the_first_draw(monkeypatch):
    def refuse(*_):
        raise MemoryError

    monkeypatch.setattr(evalue_module, "_result_table", refuse)
    generator = np.random.default_rng(1)
    before = generator.bit_generator.state
    with pytest.raises(MemoryError):
        e_value_null(Guide(GUIDE), sirna_seqs(), rng=generator, n_random=3)
    assert generator.bit_generator.state == before


def test_e_value_null_draws_the_shuffles_while_the_searches_run(monkeypatch):
    """A search starts after at most 2 * workers + 1 draws, because the pool holds 2 * workers."""
    workers = 2
    drawn = []
    seen = []
    draw = _shuffled_guides
    search = evalue_module.count_matches

    def counting_draw(guide, count, generator):
        for shuffle in draw(guide, count, generator):
            drawn.append(shuffle)
            yield shuffle

    def counting_search(target, sequences, options):
        seen.append(len(drawn))
        return search(target, sequences, options)

    monkeypatch.setattr(evalue_module, "_shuffled_guides", counting_draw)
    monkeypatch.setattr(evalue_module, "count_matches", counting_search)
    null = e_value_null(Guide(GUIDE), sirna_seqs(), rng=1, n_random=60, workers=workers)
    assert len(drawn) == 60
    assert min(seen) <= 2 * workers + 1
    monkeypatch.undo()
    pd.testing.assert_frame_equal(
        null, e_value_null(Guide(GUIDE), sirna_seqs(), rng=1, n_random=60, workers=1)
    )


def test_e_value_null_rejects_n_random_below_1():
    with pytest.raises(InputError, match=exactly("n_random 0 is less than 1; expected 1 or more")):
        e_value_null(Guide(GUIDE), sirna_seqs(), rng=1, n_random=0)


def test_e_value_null_rejects_an_n_random_that_is_not_an_integer():
    message = "n_random 1.5 has type float; expected an integer"
    with pytest.raises(InputError, match=exactly(message)):
        e_value_null(Guide(GUIDE), sirna_seqs(), rng=1, n_random=1.5)


def test_e_value_null_rejects_options_of_another_type():
    message = "options has type dict; expected a SearchOptions"
    with pytest.raises(InputError, match=exactly(message)):
        e_value_null(Guide(GUIDE), sirna_seqs(), rng=1, options={"max_mismatch": 1})


def test_e_value_null_checks_rng_and_workers_with_the_shared_messages():
    with pytest.raises(InputError, match=r"^rng 1.5 has type float; expected a numpy\.random"):
        e_value_null(Guide(GUIDE), sirna_seqs(), rng=1.5)
    message = "workers 0 is less than 1; expected 1 or more"
    with pytest.raises(InputError, match=exactly(message)):
        e_value_null(Guide(GUIDE), sirna_seqs(), rng=1, workers=0)


def test_e_value_null_rejects_a_guide_of_another_type_before_any_draw():
    generator = np.random.default_rng(1)
    before = generator.bit_generator.state
    message = "guide has type str; expected a Guide"
    with pytest.raises(InputError, match=exactly(message)):
        e_value_null("UUAUAGAG", sirna_seqs(), rng=generator)
    assert generator.bit_generator.state == before


def test_e_value_null_rejects_sequences_of_another_type_before_any_draw():
    generator = np.random.default_rng(1)
    before = generator.bit_generator.state
    message = "sequences has type list; expected a SequenceSet"
    with pytest.raises(InputError, match=exactly(message)):
        e_value_null(Guide(GUIDE), [("a", "ACGT")], rng=generator)
    assert generator.bit_generator.state == before


def test_e_values_reject_a_nan_count_D047():
    """D-047: R gives the e-value 0 to the observed count NA (the notes of evalue.tsv)."""
    de = make_de(["g3", "g4"], [1.0, np.nan])
    message = "column 'mer7m8' of de holds NaN for gene 'g4'; expected a count of 0 or more"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(de, make_null(NOTE_NULL), column="mer7m8")


def test_e_values_reject_a_missing_value_of_a_nullable_integer_column_D047():
    """D-047: a missing count is not a count. R gives the e-value 0 to it."""
    de = make_de(["a", "b"], pd.array([1, None], dtype="Int64"))
    message = "column 'mer7m8' of de holds NaN for gene 'b'; expected a count of 0 or more"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(de, make_null(NOTE_NULL), column="mer7m8")


def test_e_values_reject_a_negative_count_D047():
    """D-047: R gives the e-value 0 to the observed count -1, with the warning NaNs produced."""
    de = make_de(["a", "b", "a"], [1, -1, 3])
    message = "column 'mer7m8' of de holds -1 for gene 'b'; expected a count of 0 or more"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(de, make_null(NOTE_NULL), column="mer7m8")


def test_e_values_name_the_first_bad_count_of_a_float_column_D047():
    de = make_de(["a", "b", "a"], [1.0, -2.5, np.nan])
    message = "column 'mer7m8' of de holds -2.5 for gene 'b'; expected a count of 0 or more"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(de, make_null(NOTE_NULL), column="mer7m8")


@pytest.mark.parametrize(("value", "text"), [(np.inf, "inf"), (-np.inf, "-inf")])
def test_e_values_reject_an_infinite_count_D047(value, text):
    """D-047: R gives the e-value 0 to the observed count Inf, and to -Inf with a warning."""
    de = make_de(["a", "b"], [1.0, value])
    message = f"column 'mer7m8' of de holds {text} for gene 'b'; expected a count of 0 or more"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(de, make_null(NOTE_NULL), column="mer7m8")


@pytest.mark.parametrize(
    ("value", "text"), [(np.nan, "NaN"), (-1.0, "-1"), (np.inf, "inf"), (-np.inf, "-inf")]
)
def test_e_values_reject_an_invalid_count_for_a_gene_that_is_not_in_the_null_D047(value, text):
    """D-047: R drops the row `zz`, which is not in the null, and gives no message."""
    de = make_de(["a", "zz"], [1.0, value])
    message = f"column 'mer7m8' of de holds {text} for gene 'zz'; expected a count of 0 or more"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(de, make_null(NOTE_NULL), column="mer7m8")


@pytest.mark.parametrize(("value", "text"), [(np.inf, "inf"), (-np.inf, "-inf")])
def test_e_values_reject_an_infinite_count_in_the_null_D047(value, text):
    null = make_null([("a", [1, 2, 3, 0]), ("b", [0, value, 0, 0])])
    message = f"null holds {text} for gene 'b'; expected counts of 0 or more"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(make_de(["a"], [2]), null, column="mer7m8")


def test_e_values_reject_a_repeated_gene_of_the_null_D047():
    """D-047: R gives the e-value 1.5 to a = 2 for a null with a second row for a (5, 5, 5, 5)."""
    null = make_null([("a", [1, 2, 3, 0]), ("b", [0, 1, 0, 0]), ("a", [5, 5, 5, 5])])
    message = "null has 2 rows for gene 'a'; expected one row for each gene"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(make_de(["a"], [2]), null, column="mer7m8")


def test_e_values_reject_a_nan_count_in_the_null_D047():
    null = make_null([("a", [1, 2, 3, 0]), ("b", [0, np.nan, 0, 0])])
    message = "null holds NaN for gene 'b'; expected counts of 0 or more"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(make_de(["a"], [2]), null, column="mer7m8")


def test_e_values_reject_a_negative_count_in_the_null_D047():
    null = make_null([("a", [1, 2, 3, 0]), ("b", [0, 1, -3, 0])])
    message = "null holds -3 for gene 'b'; expected counts of 0 or more"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(make_de(["a"], [2]), null, column="mer7m8")


def test_e_values_find_a_bad_count_of_the_null_after_the_first_block_D047():
    genes = [f"g{i}" for i in range(5000)]
    counts = np.ones((5000, 3), dtype=np.int64)
    counts[4500, 1] = -1
    null = pd.DataFrame(counts, index=pd.Index(genes, name="gene_id"))
    message = "null holds -1 for gene 'g4500'; expected counts of 0 or more"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(make_de(["g0"], [1]), null, column="mer7m8")


def test_each_row_of_a_gene_gets_its_own_e_value_D047():
    """D-047: R gives 0.5 to both rows a = 1 and a = 3. One by one they get 0.75 and 0.25."""
    de = make_de(["a", "a"], [1, 3])
    result = e_values(de, make_null(NOTE_NULL), column="mer7m8")
    assert result["e_value"].tolist() == [0.75, 0.25]


def test_a_count_of_0_gives_the_e_value_1():
    result = e_values(make_de(["a", "b"], [0, 0]), make_null(NOTE_NULL), column="mer7m8")
    assert result["e_value"].tolist() == [1.0, 1.0]


def test_e_values_drop_the_rows_whose_gene_is_not_in_the_null():
    """R drops the row of z too (the notes of evalue.tsv)."""
    de = make_de(["z", "a", "y", "b"], [1, 2, 3, 1], column="mer7m8")
    result = e_values(de, make_null(NOTE_NULL), column="mer7m8")
    assert result["gene_id"].tolist() == ["a", "b"]
    assert result["e_value"].tolist() == [0.5, 0.25]
    assert result.index.tolist() == [0, 1]


def test_e_values_reject_a_table_without_a_gene_of_the_null():
    """R stops with `replacement has 1 row, data has 0` (the notes of evalue.tsv)."""
    message = (
        "de has no gene of the null; expected the gene IDs of the null in the column gene_id, "
        "for example 'a'"
    )
    with pytest.raises(InputError, match=exactly(message)):
        e_values(make_de(["z"], [1]), make_null(NOTE_NULL), column="mer7m8")


def test_e_values_reject_a_table_without_rows():
    with pytest.raises(InputError, match=r"^de has no gene of the null; expected"):
        e_values(make_de([], []), make_null(NOTE_NULL), column="mer7m8")


def test_e_values_reject_a_null_without_rows():
    null = pd.DataFrame(np.zeros((0, 3), dtype=np.int32), index=pd.Index([], name="gene_id"))
    message = "null has no row; expected one row for each gene"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(make_de(["a"], [1]), null, column="mer7m8")


def test_e_values_replace_a_column_e_value_that_is_present():
    de = make_de(["a", "b"], [1, 1]).assign(e_value=["old", "old"])
    result = e_values(de, make_null(NOTE_NULL), column="mer7m8")
    assert list(result.columns) == ["gene_id", "mer7m8", "e_value"]
    assert result["e_value"].tolist() == [0.75, 0.25]
    assert result["e_value"].dtype == np.float64


def test_e_values_keep_the_other_columns_and_do_not_change_the_arguments():
    de = make_de(["b", "a"], [1, 2]).assign(note=["x", "y"])
    null = make_null(NOTE_NULL)
    de_before, null_before = de.copy(), null.copy()
    result = e_values(de, null, column="mer7m8")
    assert list(result.columns) == ["gene_id", "mer7m8", "note", "e_value"]
    assert result["note"].tolist() == ["x", "y"]
    pd.testing.assert_frame_equal(de, de_before)
    pd.testing.assert_frame_equal(null, null_before)


def test_e_values_agree_with_a_count_of_the_null_over_blocks_of_rows():
    """9000 rows are two blocks of 4096 rows and a rest."""
    generator = np.random.default_rng(4)
    genes = [f"g{i}" for i in range(300)]
    null = pd.DataFrame(
        generator.integers(0, 6, size=(300, 30)).astype(np.int32),
        index=pd.Index(genes, name="gene_id"),
    )
    picked = generator.integers(0, 300, size=9000)
    de = make_de([genes[i] for i in picked], generator.integers(0, 7, size=9000))
    result = e_values(de, null, column="mer7m8")
    expected = [
        float((null.loc[gene].to_numpy() >= count).sum()) / 30
        for gene, count in zip(de["gene_id"], de["mer7m8"], strict=True)
    ]
    assert result["e_value"].tolist() == expected


def test_e_values_accept_a_float_null_and_a_float_count():
    null = make_null({"a": [1.0, 2.5, 3.0, 0.0]})
    result = e_values(make_de(["a", "a"], [2.0, 2.5]), null, column="mer7m8")
    assert result["e_value"].tolist() == [0.5, 0.5]


def test_e_values_of_an_e_value_null_give_a_fraction_of_the_shuffles():
    sequences = multi_sequence_set()
    guide = Guide(GUIDE)
    options = SearchOptions(max_mismatch=1)
    de = annotate_seed_counts(
        pd.DataFrame({"gene_id": sorted(set(sequences.gene_ids))}),
        sequences,
        guide.target(Seed.MER7M8),
        options=options,
    )
    null = e_value_null(guide, sequences, rng=1, n_random=20, options=options)
    result = e_values(de, null, column="mer7m8")
    assert len(result) == len(de)
    assert (result["e_value"] * 20).round().tolist() == (result["e_value"] * 20).tolist()
    assert result["e_value"].between(0, 1).all()


@pytest.mark.parametrize(("value", "type_name"), NOT_A_DATA_FRAME)
def test_e_values_reject_a_de_that_is_not_a_data_frame(value, type_name):
    message = f"de has type {type_name}; expected a pandas.DataFrame"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(value, make_null(NOTE_NULL), column="mer7m8")


@pytest.mark.parametrize(("value", "type_name"), NOT_A_DATA_FRAME)
def test_e_values_reject_a_null_that_is_not_a_data_frame(value, type_name):
    message = f"null has type {type_name}; expected a pandas.DataFrame with one row for each gene"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(make_de(["a"], [1]), value, column="mer7m8")


def test_e_values_reject_a_de_without_the_column():
    message = "DE table has missing columns: ['mer7m8']; expected ['gene_id', 'mer7m8']"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(pd.DataFrame({"gene_id": ["a"]}), make_null(NOTE_NULL), column="mer7m8")


def test_e_values_reject_a_de_with_a_repeated_column():
    de = pd.DataFrame([["a", 1, 2]], columns=["gene_id", "mer7m8", "mer7m8"])
    message = "DE table has 2 columns 'mer7m8'; expected one column for each name"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(de, make_null(NOTE_NULL), column="mer7m8")


def test_e_values_reject_a_column_name_that_is_not_text():
    message = "column 1 has type int; expected a str"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(make_de(["a"], [1]), make_null(NOTE_NULL), column=1)


def test_e_values_reject_gene_ids_that_are_not_text():
    de = pd.DataFrame({"gene_id": [1, 2], "mer7m8": [1, 1]})
    with pytest.raises(InputError, match=r"^'gene_id' has dtype int64; expected text values"):
        e_values(de, make_null(NOTE_NULL), column="mer7m8")


def test_e_values_reject_a_count_column_that_is_not_numeric():
    de = pd.DataFrame({"gene_id": ["a"], "mer7m8": ["1"]})
    message = f"column 'mer7m8' has dtype {TEXT_DTYPE}; expected a numeric dtype"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(de, make_null(NOTE_NULL), column="mer7m8")


def test_e_values_reject_a_null_with_a_gene_index_that_is_not_text():
    null = pd.DataFrame([[1, 2], [3, 4]], index=[10, 20])
    message = "null index has dtype int64; expected text gene IDs"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(make_de(["a"], [1]), null, column="mer7m8")


def test_e_values_reject_a_null_with_a_missing_gene_id():
    null = pd.DataFrame([[1, 2], [3, 4]], index=pd.Index(["a", None], dtype=object))
    message = "null index has a missing gene ID; expected text gene IDs"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(make_de(["a"], [1]), null, column="mer7m8")


def test_e_values_reject_a_null_without_columns():
    null = pd.DataFrame(index=pd.Index(["a", "b"], name="gene_id"))
    message = "null has no column; expected one column for each random guide"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(make_de(["a"], [1]), null, column="mer7m8")


def test_e_values_reject_a_null_that_repeats_a_column_label():
    null = pd.DataFrame([[1, 2], [3, 4]], index=["a", "b"], columns=[0, 0])
    message = "null repeats the column label 0; expected each label once"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(make_de(["a"], [1]), null, column="mer7m8")


def test_e_values_reject_a_null_with_a_column_that_is_not_numeric():
    null = pd.DataFrame({0: [1, 2], 1: ["x", "y"]}, index=pd.Index(["a", "b"], name="gene_id"))
    message = f"column 1 has dtype {TEXT_DTYPE}; expected a numeric dtype"
    with pytest.raises(InputError, match=exactly(message)):
        e_values(make_de(["a"], [1]), null, column="mer7m8")
