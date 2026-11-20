import math
import re
import sys
import warnings
from functools import partial

import numpy as np
import pandas as pd
import pytest
from helpers import (
    GUIDE,
    STATISTIC_TOLERANCE,
    dts_p_value_error,
    ecdf_golden_row,
    golden_rows,
    sirna_de,
    sirna_seqs,
)

from seedmatch import (
    Alternative,
    EcdfTest,
    EcdfTestResult,
    Guide,
    Seed,
    SeedMatchWarning,
    WassersteinAreas,
    annotate_seed_counts,
    assign_groups,
    ecdf_test,
    filter_de,
    wasserstein_areas,
)
from seedmatch.errors import InputError
from seedmatch.stats import PValueMethod, dts_test, ks_test, wilcoxon_test

R_ROWS = [
    pytest.param(row, id=f"{row.input}-{row.test}-{row.alternative}")
    for row in golden_rows("stats/ecdf.tsv")
]


def inputs(name):
    """Return the DE table, the target genes, and the background genes of an input of ecdf.tsv."""
    de = sirna_de()
    if name == "readme":
        # The table of the first README example: the default filter, then the seed counts.
        de = annotate_seed_counts(
            filter_de(de.drop(columns="mer7m8")), sirna_seqs(), Guide(GUIDE).target(Seed.MER7M8)
        )
    seed = de.loc[de["mer7m8"] >= 1, "gene_id"]
    background = de.loc[de["mer7m8"] == 0, "gene_id"]
    return (de, background, seed) if name == "test_data_swapped" else (de, seed, background)


def run(name, **options):
    de, target, background = inputs(name)
    return ecdf_test(de, target=target, background=background, **options)


def fold_changes(de, genes):
    return de.loc[de["gene_id"].isin(genes), "log2FoldChange"].to_numpy()


@pytest.mark.parametrize("row", R_ROWS)
def test_ecdf_test_matches_r(row):
    test = EcdfTest(row.test)
    alternative = None if row.alternative == "NA" else Alternative(row.alternative)
    options = {"rng": 1} if test is EcdfTest.DTS else {}
    result = run(row.input, test=test, alternative=alternative, **options)
    assert result.test is test
    assert result.alternative is alternative
    assert (result.n_target, result.n_background) == (int(row.n_target), int(row.n_background))
    if test is EcdfTest.WASSERSTEIN:
        assert result.statistic == float(row.statistic)
        assert isinstance(result.areas, WassersteinAreas)
        assert result.statistic == result.areas.negative_censored
        assert (result.p_value, result.method) == (None, None)
        return
    assert result.areas is None
    assert math.isclose(result.statistic, float(row.statistic), **STATISTIC_TOLERANCE)
    assert result.p_value is not None
    if test is EcdfTest.DTS:
        p_r = float(row.p_value)
        assert abs(result.p_value - p_r) <= dts_p_value_error(p_r, 2000)
    else:
        assert math.isclose(result.p_value, float(row.p_value), rel_tol=1e-9)


@pytest.mark.parametrize("test", [EcdfTest.KS, EcdfTest.WILCOXON])
def test_ecdf_test_uses_greater_when_the_alternative_is_none(test):
    default = run("test_data", test=test)
    assert default.alternative is Alternative.GREATER
    assert default == run("test_data", test=test, alternative=Alternative.GREATER)


def test_ecdf_test_uses_the_ks_test_by_default():
    assert run("test_data") == run("test_data", test=EcdfTest.KS)


def test_ecdf_test_reports_the_method_of_the_test():
    methods = {
        test: run("test_data", test=test, **({"rng": 1} if test is EcdfTest.DTS else {})).method
        for test in EcdfTest
    }
    assert methods == {
        EcdfTest.KS: PValueMethod.EXACT,
        EcdfTest.WILCOXON: PValueMethod.ASYMPTOTIC,
        EcdfTest.WASSERSTEIN: None,
        EcdfTest.DTS: PValueMethod.PERMUTATION,
    }


def test_ecdf_test_gives_the_samples_of_the_rows_to_the_tests():
    de, target, background = inputs("test_data")
    x, y = fold_changes(de, set(target)), fold_changes(de, set(background))
    dts = run("test_data", test=EcdfTest.DTS, rng=7, n_boot=50)
    assert dts.p_value == dts_test(x, y, rng=7, n_boot=50).p_value
    wass = run("test_data", test=EcdfTest.WASSERSTEIN, n_interp=100)
    assert wass.areas == wasserstein_areas(x, y, n_interp=100)


def test_ecdf_test_draws_the_permutations_from_the_generator_that_the_caller_gives():
    generator = np.random.default_rng(7)
    before = generator.bit_generator.state
    run("test_data", test=EcdfTest.DTS, rng=generator, n_boot=10)
    assert generator.bit_generator.state != before


@pytest.mark.parametrize("test", list(EcdfTest))
def test_ecdf_test_drops_rows_without_a_fold_change_D036(test):
    de, target, background = inputs("test_data")
    # R drops these rows for KS and Wilcoxon, and it stops for DTS and Wass (D-036).
    missing = pd.DataFrame({"gene_id": [f"missing{i}" for i in range(5)], "log2FoldChange": np.nan})
    options = {"rng": 1} if test is EcdfTest.DTS else {}
    with_missing = ecdf_test(
        pd.concat([de, missing], ignore_index=True),
        target=[*target, "missing0", "missing1"],
        background=[*background, "missing2", "missing3", "missing4"],
        test=test,
        **options,
    )
    assert with_missing == run("test_data", test=test, **options)
    assert (with_missing.n_target, with_missing.n_background) == (18, 161)
    if test in (EcdfTest.KS, EcdfTest.WILCOXON):
        row = ecdf_golden_row("test_data", test.value, "greater")
        assert math.isclose(with_missing.statistic, float(row.statistic), **STATISTIC_TOLERANCE)
        assert math.isclose(with_missing.p_value, float(row.p_value), rel_tol=1e-9)


# The direct functions and the names of their two samples in a message.
DIRECT_TESTS = {
    "ks_test": (ks_test, "x", "y"),
    "wilcoxon_test": (wilcoxon_test, "x", "y"),
    "dts_test": (partial(dts_test, rng=1, n_boot=10), "x", "y"),
    "wasserstein_areas": (wasserstein_areas, "seed", "non_seed"),
}


@pytest.mark.parametrize("name", DIRECT_TESTS)
def test_the_direct_tests_raise_for_nan_where_ecdf_test_drops_the_row_D036(name):
    call, first, second = DIRECT_TESTS[name]
    values = [float(k) for k in range(12)]
    with pytest.raises(InputError, match=rf"^{first}\[3\] is NaN"):
        call([*values[:3], np.nan, *values[3:]], values)
    with pytest.raises(InputError, match=rf"^{second}\[0\] is NaN"):
        call(values, [np.nan, *values])


@pytest.mark.parametrize("test", [EcdfTest.DTS, EcdfTest.WASSERSTEIN])
@pytest.mark.parametrize("alternative", list(Alternative))
def test_ecdf_test_rejects_an_alternative_for_dts_and_wasserstein_D035(test, alternative):
    options = {"rng": 1} if test is EcdfTest.DTS else {}
    message = (
        f"alternative Alternative.{alternative.name} does not apply to EcdfTest.{test.name}; "
        "expected None, because the test has no direction"
    )
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        run("test_data", test=test, alternative=alternative, **options)


def test_ecdf_test_names_an_alternative_that_is_not_a_member_for_dts_by_its_repr():
    message = "alternative 'less' does not apply to EcdfTest.DTS; expected None, because"
    with pytest.raises(InputError, match=f"^{re.escape(message)}"):
        run("test_data", test=EcdfTest.DTS, alternative="less", rng=1)


@pytest.mark.parametrize("alternative", ["greater", "less", 1])
def test_ecdf_test_rejects_an_alternative_that_is_not_a_member(alternative):
    with pytest.raises(InputError, match=r"^alternative .* is not valid; expected one of"):
        run("test_data", alternative=alternative)


@pytest.mark.parametrize("name", ["target", "background"])
def test_ecdf_test_rejects_a_bare_string(name):
    de, target, background = inputs("test_data")
    genes = {"target": target, "background": background}
    genes[name] = "ENSRNOG00000016275"
    message = (
        f"{name} 'ENSRNOG00000016275' is a str; expected an iterable of gene IDs, "
        "for example ('ENSRNOG00000016275',)"
    )
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        ecdf_test(de, **genes)


@pytest.mark.parametrize("name", ["target", "background"])
@pytest.mark.parametrize("genes", [[], set(), ["not_a_gene"], pd.Series(dtype=object)])
def test_ecdf_test_rejects_a_set_without_rows(name, genes):
    de, target, background = inputs("test_data")
    sets = {"target": target, "background": background}
    sets[name] = genes
    message = rf"^{name} has 0 rows with a log2FoldChange in the DE table; expected at least 1"
    with pytest.raises(InputError, match=message):
        ecdf_test(de, **sets)


def test_ecdf_test_rejects_a_set_whose_rows_have_no_fold_change():
    de, target, background = inputs("test_data")
    de.loc[de["gene_id"].isin(set(target)), "log2FoldChange"] = np.nan
    with pytest.raises(InputError, match=r"^target has 0 rows with a log2FoldChange"):
        ecdf_test(de, target=target, background=background)


@pytest.mark.parametrize("genes", [[1, 2], ["gene", None], ["gene", float("nan")], [b"gene"]])
def test_ecdf_test_rejects_a_gene_id_that_is_not_text(genes):
    de, _, background = inputs("test_data")
    with pytest.raises(InputError, match=r"^target\[\d\] .* has type \w+; expected a str$"):
        ecdf_test(de, target=genes, background=background)


@pytest.mark.parametrize("genes", [5, None, 1.5])
def test_ecdf_test_rejects_genes_that_are_not_iterable(genes):
    de, _, background = inputs("test_data")
    with pytest.raises(InputError, match=r"^target .* has type \w+; expected an iterable$"):
        ecdf_test(de, target=genes, background=background)


def test_ecdf_test_accepts_each_kind_of_iterable_of_gene_ids():
    de, target, background = inputs("test_data")
    expected = ecdf_test(de, target=list(target), background=list(background))
    kinds = [
        lambda genes: genes,
        set,
        tuple,
        iter,
        lambda genes: np.array(list(genes)),
        pd.Index,
        dict.fromkeys,
        lambda genes: pd.Series(list(genes), dtype="category"),
    ]
    for kind in kinds:
        result = ecdf_test(de, target=kind(target), background=kind(background))
        assert result == expected


def test_ecdf_test_takes_each_row_of_a_gene_and_each_gene_once():
    de = pd.DataFrame(
        {"gene_id": ["a", "b", "b", "c", "d"], "log2FoldChange": [0.1, 0.2, 0.3, 0.4, 0.5]}
    )
    result = ecdf_test(de, target=["b", "b", "a"], background=["c", "d"])
    assert (result.n_target, result.n_background) == (3, 2)


def test_ecdf_test_needs_rng_for_dts():
    with pytest.raises(
        InputError, match=r"^rng is None; expected a numpy\.random\.Generator or an"
    ):
        run("test_data", test=EcdfTest.DTS)


def test_ecdf_test_reports_a_missing_rng_before_it_selects_the_rows():
    de, _, background = inputs("test_data")
    with pytest.raises(InputError, match=r"^rng is None; expected"):
        ecdf_test(de, target=["not_a_gene"], background=background, test=EcdfTest.DTS)


def sets_of(name, select):
    """Return the table and the sets of ecdf_test. `select` changes the set `name`."""
    de, target, background = inputs("test_data")
    sets = {"target": target, "background": background}
    sets[name] = select(sets[name])
    return de, sets


@pytest.mark.parametrize("test", [EcdfTest.DTS, EcdfTest.WASSERSTEIN])
@pytest.mark.parametrize("set_name", ["target", "background"])
@pytest.mark.parametrize("value", [np.inf, -np.inf])
def test_ecdf_test_rejects_an_infinite_fold_change_for_dts_and_wasserstein(test, set_name, value):
    # The set is in the reverse order of the table, so the first row of the table is not the first
    # gene of the set.
    de, sets = sets_of(set_name, lambda genes: list(reversed(list(genes))))
    first, second = de.index[de["gene_id"].isin(set(sets[set_name]))][:2]
    de.loc[first, "log2FoldChange"] = value
    de.loc[second, "log2FoldChange"] = np.inf
    options = {"rng": 1} if test is EcdfTest.DTS else {}
    message = (
        f"{set_name} gene {de.loc[first, 'gene_id']!r} has log2FoldChange {value}; "
        f"expected a finite value for EcdfTest.{test.name}"
    )
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        ecdf_test(de, **sets, test=test, **options)


@pytest.mark.parametrize("test", [EcdfTest.KS, EcdfTest.WILCOXON])
def test_ecdf_test_accepts_an_infinite_fold_change_for_ks_and_wilcoxon(test):
    de, target, background = inputs("test_data")
    de.loc[de.index[de["gene_id"].isin(set(target))][0], "log2FoldChange"] = np.inf
    direct = {EcdfTest.KS: ks_test, EcdfTest.WILCOXON: wilcoxon_test}[test]
    result = ecdf_test(de, target=target, background=background, test=test)
    x, y = fold_changes(de, set(target)), fold_changes(de, set(background))
    assert (result.n_target, result.n_background) == (18, 161)
    assert result.statistic == direct(x, y, alternative=Alternative.GREATER).statistic


@pytest.mark.parametrize("set_name", ["target", "background"])
def test_ecdf_test_rejects_a_wasserstein_set_of_fewer_than_10_values(set_name):
    de, sets = sets_of(set_name, lambda genes: list(genes)[:9])
    with pytest.raises(InputError, match=rf"^{set_name} has 9 values; expected at least 10$"):
        ecdf_test(de, **sets, test=EcdfTest.WASSERSTEIN)
    # The minimum belongs to the Wasserstein areas only.
    assert ecdf_test(de, **sets).test is EcdfTest.KS


def test_ecdf_test_counts_the_values_after_the_drop_for_the_wasserstein_minimum():
    de, target, background = inputs("test_data")
    chosen = list(target)[:12]
    de.loc[de["gene_id"].isin(chosen[:3]), "log2FoldChange"] = np.nan
    with pytest.raises(InputError, match=r"^target has 9 values; expected at least 10$"):
        ecdf_test(de, target=chosen, background=background, test=EcdfTest.WASSERSTEIN)


@pytest.mark.parametrize(("set_name", "count"), [("target", 18), ("background", 161)])
def test_ecdf_test_rejects_a_wasserstein_set_with_one_distinct_value(set_name, count):
    de, sets = sets_of(set_name, lambda genes: genes)
    de.loc[de["gene_id"].isin(set(sets[set_name])), "log2FoldChange"] = 0.5
    message = rf"^{set_name} has {count} values that are all 0\.5; expected at least 2 distinct"
    with pytest.raises(InputError, match=message + " values$"):
        ecdf_test(de, **sets, test=EcdfTest.WASSERSTEIN)


# The values that an option can have for another test, and the one test that the option applies to.
# 2000 and 10000 are the defaults of DTS and Wasserstein, but ecdf_test rejects them for the others.
OPTIONS = {
    "rng": ([1], EcdfTest.DTS),
    "n_boot": ([10, 2000], EcdfTest.DTS),
    "n_interp": ([100, 10000], EcdfTest.WASSERSTEIN),
}


@pytest.mark.parametrize(
    ("test", "name", "value"),
    [
        (test, name, value)
        for test in EcdfTest
        for name, (values, owner) in OPTIONS.items()
        if test is not owner
        for value in values
    ],
)
def test_ecdf_test_rejects_options_of_another_test(test, name, value):
    owner = OPTIONS[name][1]
    options = {name: value, **({"rng": 1} if test is EcdfTest.DTS else {})}
    message = (
        rf"^{name} {value} applies only to EcdfTest\.{owner.name}, not to EcdfTest\.{test.name}; "
        r"expected None$"
    )
    with pytest.raises(InputError, match=message):
        run("test_data", test=test, **options)


def test_ecdf_test_accepts_none_for_an_option_of_another_test():
    assert run("test_data", n_boot=None, n_interp=None, rng=None) == run("test_data")


def test_ecdf_test_uses_the_default_of_the_test_when_n_boot_and_n_interp_are_none():
    dts = run("test_data", test=EcdfTest.DTS, rng=1)
    assert dts == run("test_data", test=EcdfTest.DTS, rng=1, n_boot=2000)
    assert dts == run("test_data", test=EcdfTest.DTS, rng=1, n_boot=None)
    assert dts != run("test_data", test=EcdfTest.DTS, rng=1, n_boot=10)
    wass = run("test_data", test=EcdfTest.WASSERSTEIN)
    assert wass == run("test_data", test=EcdfTest.WASSERSTEIN, n_interp=10000)
    assert wass == run("test_data", test=EcdfTest.WASSERSTEIN, n_interp=None)


@pytest.mark.parametrize("test", ["KS", "Wass", None, 1])
def test_ecdf_test_rejects_a_test_that_is_not_a_member(test):
    message = r"^test .* is not valid; expected one of EcdfTest\.KS, EcdfTest\.WILCOXON"
    with pytest.raises(InputError, match=message):
        run("test_data", test=test)


@pytest.mark.parametrize("column", ["gene_id", "log2FoldChange"])
def test_ecdf_test_rejects_a_missing_column(column):
    de, target, background = inputs("test_data")
    with pytest.raises(InputError, match=rf"^DE table has missing columns: \['{column}'\]"):
        ecdf_test(de.drop(columns=column), target=target, background=background)


def test_ecdf_test_rejects_a_fold_change_column_that_is_not_numeric():
    de, target, background = inputs("test_data")
    de["log2FoldChange"] = de["log2FoldChange"].astype(str)
    message = r"^column 'log2FoldChange' has dtype \w+; expected a numeric dtype"
    with pytest.raises(InputError, match=message):
        ecdf_test(de, target=target, background=background)


@pytest.mark.parametrize("dtype", ["bool", "complex128"])
def test_ecdf_test_rejects_a_bool_and_a_complex_fold_change_column(dtype):
    de, target, background = inputs("test_data")
    de["log2FoldChange"] = de["log2FoldChange"].astype(dtype)
    message = rf"^column 'log2FoldChange' has dtype {dtype}; expected a numeric dtype"
    with pytest.raises(InputError, match=message):
        ecdf_test(de, target=target, background=background)


def test_ecdf_test_rejects_gene_ids_that_are_not_text():
    de, target, background = inputs("test_data")
    de["gene_id"] = np.arange(len(de))
    with pytest.raises(InputError, match=r"^'gene_id' has dtype int64; expected text values"):
        ecdf_test(de, target=target, background=background)


def test_ecdf_test_does_not_change_the_table():
    de, target, background = inputs("test_data")
    before = de.copy()
    ecdf_test(de, target=target, background=background)
    pd.testing.assert_frame_equal(de, before)


def test_ecdf_test_result_is_a_frozen_record_with_plain_values():
    result = run("test_data")
    assert isinstance(result, EcdfTestResult)
    assert {type(result.statistic), type(result.p_value)} == {float}
    assert (type(result.n_target), type(result.n_background)) == (int, int)
    with pytest.raises(AttributeError):
        result.statistic = 0.0


GROUP_TABLE = pd.DataFrame(
    {"gene_id": ["g1", "g2", "g3", "g4", "g2"], "log2FoldChange": [0.1, 0.2, 0.3, 0.4, 0.5]},
    index=[10, 11, 12, 13, 14],
)


def test_assign_groups_gives_the_last_group_and_warns_for_an_overlap():
    with pytest.warns(SeedMatchWarning, match=r"groups 'a' and 'b' share the gene 'g2'") as caught:
        result = assign_groups(GROUP_TABLE, {"a": ["g1", "g2"], "b": ["g2", "g3"]})
    assert len(caught) == 1
    expected = pd.Series(
        pd.Categorical(["a", "b", "b", None, "b"], categories=["a", "b"]),
        index=GROUP_TABLE.index,
        name="group",
    )
    pd.testing.assert_series_equal(result, expected)


def test_assign_groups_names_the_first_shared_gene_of_the_first_pair():
    groups = {"a": ["g3", "g1", "g2"], "b": ["g4"], "c": ["g9", "g2", "g1", "g3"]}
    with pytest.warns(SeedMatchWarning) as caught:
        assign_groups(GROUP_TABLE, groups)
    assert len(caught) == 1
    assert "groups 'a' and 'c' share the gene 'g2'" in str(caught[0].message)


def test_assign_groups_keeps_the_order_of_the_groups():
    groups = {"z": ["g1"], "m": ["not_in_the_table"], "a": ["g3"]}
    result = assign_groups(GROUP_TABLE, groups)
    assert list(result.cat.categories) == ["z", "m", "a"]
    assert result.cat.ordered is False
    assert result.isna().tolist() == [False, True, False, True, True]
    assert result.dropna().tolist() == ["z", "a"]


def test_assign_groups_gives_no_warning_for_a_repeat_inside_one_group_D039():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        result = assign_groups(GROUP_TABLE, {"a": ["g1", "g1"], "b": ["g2"]})
    assert result.dropna().tolist() == ["a", "b", "b"]


def test_assign_groups_gives_no_warning_without_an_overlap():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        result = assign_groups(GROUP_TABLE, {"a": ["g1"], "b": ["g3", "g4"]})
    assert result.dropna().tolist() == ["a", "b", "b"]


def test_assign_groups_gives_each_row_of_a_gene_the_group_of_the_gene():
    result = assign_groups(GROUP_TABLE, {"a": ["g2"]})
    assert result.dropna().index.tolist() == [11, 14]


def test_assign_groups_with_no_group_gives_a_missing_value_for_each_row():
    result = assign_groups(GROUP_TABLE, {})
    assert result.isna().all()
    assert len(result) == len(GROUP_TABLE)
    assert result.name == "group"


def test_assign_groups_rejects_a_bare_string():
    message = "groups['a'] 'g1' is a str; expected an iterable of gene IDs, for example ('g1',)"
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        assign_groups(GROUP_TABLE, {"a": "g1"})


def test_assign_groups_rejects_groups_that_are_not_a_mapping():
    with pytest.raises(InputError, match=r"^groups .* has type list; expected a mapping"):
        assign_groups(GROUP_TABLE, [["g1"], ["g2"]])


def test_assign_groups_rejects_a_missing_column_and_gene_ids_that_are_not_text():
    with pytest.raises(InputError, match=r"^DE table has missing columns: \['gene_id'\]"):
        assign_groups(GROUP_TABLE.drop(columns="gene_id"), {"a": ["g1"]})
    numbers = GROUP_TABLE.assign(gene_id=range(len(GROUP_TABLE)))
    with pytest.raises(InputError, match=r"^'gene_id' has dtype int64; expected text values"):
        assign_groups(numbers, {"a": ["g1"]})


def test_assign_groups_does_not_change_the_table():
    before = GROUP_TABLE.copy()
    assign_groups(GROUP_TABLE, {"a": ["g1"]})
    pd.testing.assert_frame_equal(GROUP_TABLE, before)


def test_assign_groups_warns_at_the_line_of_the_caller():
    with pytest.warns(SeedMatchWarning) as record:
        assign_groups(GROUP_TABLE, {"a": ["g1"], "b": ["g1"]})
    assert record[0].filename == __file__


def test_assign_groups_keeps_a_group_named_na_apart_D040():
    # R fills the text "NA" for a gene in no group, and the level "NA" then holds these genes too.
    result = assign_groups(GROUP_TABLE, {"a": ["g1"], "NA": ["g2"]})
    assert list(result.cat.categories) == ["a", "NA"]
    assert result.isna().tolist() == [False, False, True, True, False]
    assert result.loc[[11, 14]].tolist() == ["NA", "NA"]
    assert result.loc[[12, 13]].isna().all()


@pytest.mark.parametrize(
    ("test", "maximum"),
    [
        # The DTS limit is the largest float divided by 2 N, and N = 18 + 161 values.
        (EcdfTest.DTS, sys.float_info.max / (2 * 179)),
        (EcdfTest.WASSERSTEIN, sys.float_info.max / 2),
    ],
)
def test_ecdf_test_rejects_a_range_that_a_float_cannot_hold_with_the_names_of_the_sets(
    test, maximum
):
    de, target, background = inputs("test_data")
    de.loc[de.index[de["gene_id"].isin(set(target))][0], "log2FoldChange"] = 1e308
    de.loc[de.index[de["gene_id"].isin(set(background))][0], "log2FoldChange"] = -1e308
    options = {"rng": 1} if test is EcdfTest.DTS else {}
    message = (
        f"target and background span from -1e+308 to 1e+308; expected a range of at most {maximum}"
    )
    with pytest.raises(InputError, match=f"^{re.escape(message)}$"):
        ecdf_test(de, target=target, background=background, test=test, **options)
