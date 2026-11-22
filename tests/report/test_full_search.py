import dataclasses

import numpy as np
import pandas as pd
import pytest
from helpers import GUIDE, TEXT_DTYPE, assert_frame_matches_golden, golden, sirna_seqs

import seedmatch
from seedmatch import (
    CustomSeed,
    Guide,
    SearchOptions,
    Seed,
    SequenceSet,
    count_matches,
    find_matches,
)
from seedmatch.errors import InputError
from seedmatch.report import grid
from seedmatch.report.grid import (
    DEFAULT_SEEDS,
    NONSEED_SEEDS,
    FullSearchOptions,
    FullSearchRangesOptions,
    best_match,
    full_search,
    full_search_ranges,
)

FULL_TARGET = Guide(GUIDE).target(Seed.FULL).sequence
# Every canonical seed matches this subject exactly, so each search of a grid has a match.
SUBJECT_WITH_EVERY_SEED = "GGG" + FULL_TARGET + "AAA"

# The order of R/full_search.R:98-121. The tests spell it out, so that a change of the order fails.
RANKING = [
    *[(f"Full.mm{d}", f"Full-D{d}") for d in range(5)],
    *[(f"18mer.mm{d}", f"18mer-D{d}") for d in range(5)],
    *[(f"15mer.mm{d}", f"15mer-D{d}") for d in range(5)],
    ("mer8.mm0", "8mer"),
    ("mer7m8.mm0", "7mer-m8"),
    ("mer7A1.mm0", "7mer-A1"),
    ("mer6.mm0", "6mer"),
]


def seqs_of(*texts, gene_ids=None):
    records = [(f"t{i}", text) for i, text in enumerate(texts)]
    return SequenceSet.from_records(records, gene_ids=gene_ids)


@pytest.mark.parametrize(
    ("variant", "options"),
    [
        ("default", FullSearchOptions()),
        ("wobble", FullSearchOptions(wobble=True)),
        ("indels", FullSearchOptions(indels=True)),
    ],
)
def test_full_search_counts_match_r(variant, options):
    got = full_search(Guide(GUIDE), sirna_seqs(), options)
    assert_frame_matches_golden(got, f"seedmatchr/full_search_counts_{variant}.tsv")


def test_the_indels_golden_differs_from_the_default_golden():
    default = golden("seedmatchr/full_search_counts_default.tsv")
    indels = golden("seedmatchr/full_search_counts_indels.tsv")
    assert not default.equals(indels)


def test_best_match_skips_absent_columns_D012():
    got = full_search(Guide(GUIDE), sirna_seqs(), FullSearchOptions(nonseed_edit_distances=(0, 1)))
    made = {
        *(f"{seed}-D{d}" for seed in ("Full", "18mer", "15mer") for d in (0, 1)),
        *("8mer", "7mer-m8", "7mer-A1", "6mer", "No.matches"),
    }
    assert set(got["best_match"]) <= made
    # A label of a column that the search makes is the highest present label in both searches.
    default = golden("seedmatchr/full_search_counts_default.tsv")
    kept = default["best.match"].isin(made).to_numpy()
    assert {"Full-D1", "6mer", "No.matches"} <= set(default.loc[kept, "best.match"])
    assert got.loc[kept, "best_match"].tolist() == default.loc[kept, "best.match"].tolist()


def test_best_match_of_a_frame_without_count_columns_is_no_matches_D012():
    frame = pd.DataFrame({"seq_id": ["a", "b"], "other": [3, 4]})
    assert best_match(frame).tolist() == ["No.matches", "No.matches"]


def test_best_match_of_a_frame_without_rows_is_empty_D012():
    got = best_match(pd.DataFrame({"Full.mm0": pd.Series(dtype="int64")}))
    assert got.tolist() == []
    assert str(got.dtype) == TEXT_DTYPE


def test_best_match_ranks_the_columns_in_the_order_of_r():
    names = [column for column, _ in RANKING]
    rows = [{name: int(rank >= first) for rank, name in enumerate(names)} for first in range(19)]
    rows.append(dict.fromkeys(names, 0))
    got = best_match(pd.DataFrame(rows))
    assert got.tolist() == [label for _, label in RANKING] + ["No.matches"]


def test_best_match_needs_a_count_of_at_least_one():
    frame = pd.DataFrame({"Full.mm0": [0, 2], "mer6.mm0": [1, 0]})
    assert best_match(frame).tolist() == ["6mer", "Full-D0"]


@pytest.mark.parametrize("dtype", ["Int64", "float64"])
def test_best_match_treats_a_missing_count_as_no_match(dtype):
    frame = pd.DataFrame(
        {
            "Full.mm0": pd.Series([None, 1, None, 0], dtype=dtype),
            "mer6.mm0": pd.Series([1, None, None, 0], dtype=dtype),
        }
    )
    assert best_match(frame).tolist() == ["6mer", "Full-D0", "No.matches", "No.matches"]


def test_best_match_keeps_the_index_and_is_named():
    got = best_match(pd.DataFrame({"mer8.mm0": [1, 0]}, index=["x", "y"]))
    assert got.index.tolist() == ["x", "y"]
    assert got.name == "best_match"


def test_full_search_ranges_match_patched_r():
    got = full_search_ranges(Guide(GUIDE), sirna_seqs())
    assert_frame_matches_golden(got, "seedmatchr/full_search_ranges.tsv")


def test_dataframe_mode_matches_patched_r():
    seqs = sirna_seqs()
    position = seqs.ids.index("ENSRNOG00000016275")
    one = SequenceSet.from_records([(seqs.ids[position], seqs.sequence(position))])
    guide = Guide(GUIDE)
    partial = find_matches(
        guide.target(CustomSeed("Partial", 2, 19)), one, SearchOptions(max_mismatch=3, indels=True)
    )
    mer7m8 = find_matches(
        guide.target(Seed.MER7M8), one, SearchOptions(max_mismatch=1, indels=True)
    )
    got = pd.concat([mer7m8, partial], ignore_index=True)
    assert_frame_matches_golden(got, "seedmatchr/dataframe_mode.tsv")


def test_short_guide_fails_before_any_search():
    with pytest.raises(InputError, match=r"18mer needs a guide of at least 19 nt, got 18"):
        full_search(Guide("UUAUAGAGCAAGAACACU"), sirna_seqs())


@pytest.mark.parametrize("function", ["full_search", "full_search_ranges"])
def test_short_guide_starts_no_search(function, monkeypatch):
    def fail(*_args, **_kwargs):
        message = "a search ran"
        raise AssertionError(message)

    monkeypatch.setattr(grid, "count_matches", fail)
    monkeypatch.setattr(grid, "find_matches", fail)
    with pytest.raises(InputError, match=r"18mer needs a guide of at least 19 nt, got 18"):
        getattr(grid, function)(Guide("UUAUAGAGCAAGAACACU"), sirna_seqs())


def test_the_options_default_to_the_arguments_of_r():
    options = FullSearchOptions()
    assert options.seeds == DEFAULT_SEEDS
    assert [seed.value for seed in DEFAULT_SEEDS] == [
        "Full",
        "18mer",
        "15mer",
        "mer8",
        "mer7m8",
        "mer7A1",
        "mer6",
    ]
    assert {Seed.FULL, Seed.MER18, Seed.MER15} == NONSEED_SEEDS
    assert options.nonseed_edit_distances == (0, 1, 2, 3, 4)
    assert options.seed_edit_distances == (0, 1)
    assert (options.indels, options.wobble) == (False, False)
    assert (options.fixed_pattern, options.fixed_subject) == (False, False)
    ranges = FullSearchRangesOptions()
    assert ranges.indel_options == (True, False)
    assert ranges.edit_distances == (0, 1, 2, 3, 4)
    assert ranges.seeds == DEFAULT_SEEDS
    assert (ranges.fixed_pattern, ranges.fixed_subject) == (False, False)
    assert ranges.max_seed_edit_distance == 1


@pytest.mark.parametrize("make", [FullSearchOptions, FullSearchRangesOptions])
def test_options_are_frozen(make):
    with pytest.raises(dataclasses.FrozenInstanceError):
        make().seeds = (Seed.FULL,)


def test_options_normalize_their_tuples():
    options = FullSearchOptions(nonseed_edit_distances=[0, np.int64(1)], indels=np.True_)
    assert options.nonseed_edit_distances == (0, 1)
    assert type(options.nonseed_edit_distances[1]) is int
    assert options.indels is True
    ranges = FullSearchRangesOptions(edit_distances=[2, np.int64(3)], indel_options=[np.False_])
    assert ranges.edit_distances == (2, 3)
    assert type(ranges.edit_distances[1]) is int
    assert ranges.indel_options == (False,)
    assert type(ranges.indel_options[0]) is bool


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"seeds": ()}, r"seeds is empty"),
        ({"seeds": ("mer7m8",)}, r"seeds holds 'mer7m8' of type str; expected a Seed"),
        ({"seeds": (CustomSeed("x", 2, 8),)}, r"seeds holds .*CustomSeed.*expected a Seed"),
        ({"seeds": Seed.FULL}, r"seeds .*type Seed; expected a tuple or a list$"),
        ({"nonseed_edit_distances": ()}, r"nonseed_edit_distances is empty"),
        ({"seed_edit_distances": ()}, r"seed_edit_distances is empty"),
        (
            {"nonseed_edit_distances": (0, -1)},
            r"^nonseed_edit_distances -1 is negative; expected 0 or more$",
        ),
        (
            {"seed_edit_distances": (True,)},
            r"^seed_edit_distances True has type bool; expected an integer$",
        ),
        (
            {"seed_edit_distances": (1.5,)},
            r"^seed_edit_distances 1.5 has type float; expected an integer$",
        ),
        (
            {"nonseed_edit_distances": 3},
            r"^nonseed_edit_distances 3 has type int; expected a tuple or a list$",
        ),
        ({"indels": "yes"}, r"indels 'yes' has type str; expected a bool"),
        ({"wobble": 1}, r"wobble 1 has type int; expected a bool"),
        ({"fixed_pattern": None}, r"fixed_pattern None has type NoneType; expected a bool"),
        ({"fixed_subject": 0}, r"fixed_subject 0 has type int; expected a bool"),
    ],
)
def test_full_search_options_reject_invalid_input(kwargs, message):
    with pytest.raises(InputError, match=message):
        FullSearchOptions(**kwargs)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"seeds": ()}, r"seeds is empty"),
        ({"seeds": ("Full",)}, r"seeds holds 'Full' of type str; expected a Seed"),
        ({"indel_options": ()}, r"indel_options is empty"),
        ({"indel_options": (1,)}, r"indel_options 1 has type int; expected a bool"),
        ({"edit_distances": ()}, r"edit_distances is empty"),
        ({"edit_distances": (-2,)}, r"^edit_distances -2 is negative; expected 0 or more$"),
        ({"fixed_pattern": "no"}, r"fixed_pattern 'no' has type str; expected a bool"),
        ({"fixed_subject": 1}, r"fixed_subject 1 has type int; expected a bool"),
    ],
)
def test_full_search_ranges_options_reject_invalid_input(kwargs, message):
    with pytest.raises(InputError, match=message):
        FullSearchRangesOptions(**kwargs)


@pytest.mark.parametrize(
    ("make", "kwargs", "message"),
    [
        (
            FullSearchOptions,
            {"seeds": (Seed.FULL, Seed.MER6, Seed.FULL)},
            r"^seeds repeats 'Full'; expected each seed once$",
        ),
        (
            FullSearchOptions,
            {"nonseed_edit_distances": (0, 0)},
            r"^nonseed_edit_distances repeats 0; expected each distance once$",
        ),
        (
            FullSearchOptions,
            {"seed_edit_distances": (0, 1, 0)},
            r"^seed_edit_distances repeats 0; expected each distance once$",
        ),
        (
            FullSearchRangesOptions,
            {"seeds": (Seed.MER6, Seed.MER6)},
            r"^seeds repeats 'mer6'; expected each seed once$",
        ),
        (
            FullSearchRangesOptions,
            {"edit_distances": (1, 1)},
            r"^edit_distances repeats 1; expected each distance once$",
        ),
        (
            FullSearchRangesOptions,
            {"indel_options": (True, True)},
            r"^indel_options repeats True; expected each value once$",
        ),
    ],
)
def test_a_duplicate_value_is_rejected_D023(make, kwargs, message):
    with pytest.raises(InputError, match=message):
        make(**kwargs)


@pytest.mark.parametrize(
    ("value", "message"),
    [
        (-1, r"^max_seed_edit_distance -1 is negative; expected 0 or more$"),
        (1.5, r"^max_seed_edit_distance 1.5 has type float; expected an integer$"),
        (True, r"^max_seed_edit_distance True has type bool; expected an integer$"),
    ],
)
def test_max_seed_edit_distance_must_be_an_integer_of_0_or_more_D023(value, message):
    with pytest.raises(InputError, match=message):
        FullSearchRangesOptions(max_seed_edit_distance=value)


@pytest.mark.parametrize("unordered", [set, frozenset])
@pytest.mark.parametrize(
    ("make", "field", "values"),
    [
        (FullSearchOptions, "seeds", (Seed.FULL, Seed.MER6)),
        (FullSearchOptions, "nonseed_edit_distances", (0, 1)),
        (FullSearchOptions, "seed_edit_distances", (0, 1)),
        (FullSearchRangesOptions, "seeds", (Seed.FULL, Seed.MER6)),
        (FullSearchRangesOptions, "edit_distances", (0, 1)),
        (FullSearchRangesOptions, "indel_options", (True, False)),
    ],
)
def test_an_unordered_collection_is_rejected_because_its_order_gives_the_output_order(
    make, field, values, unordered
):
    message = rf"^{field} .* has type {unordered.__name__}; expected a tuple or a list$"
    with pytest.raises(InputError, match=message):
        make(**{field: unordered(values)})


@pytest.mark.parametrize(
    ("workers", "message"),
    [
        (0, r"workers 0 is less than 1; expected 1 or more"),
        (-2, r"workers -2 is less than 1; expected 1 or more"),
        (True, r"workers True has type bool; expected an integer"),
        (1.5, r"workers 1.5 has type float; expected an integer"),
        ("4", r"workers '4' has type str; expected an integer"),
    ],
)
@pytest.mark.parametrize("function", [full_search, full_search_ranges])
def test_workers_must_be_an_int_of_one_or_more(function, workers, message):
    with pytest.raises(InputError, match=message):
        function(Guide(GUIDE), seqs_of("ACGT"), workers=workers)


def test_workers_can_be_a_numpy_integer():
    one = full_search(Guide(GUIDE), sirna_seqs(), workers=1)
    pd.testing.assert_frame_equal(full_search(Guide(GUIDE), sirna_seqs(), workers=np.int64(3)), one)


def test_full_search_does_not_depend_on_workers():
    one = full_search(Guide(GUIDE), sirna_seqs(), workers=1)
    four = full_search(Guide(GUIDE), sirna_seqs(), workers=4)
    default = full_search(Guide(GUIDE), sirna_seqs())
    pd.testing.assert_frame_equal(one, four)
    pd.testing.assert_frame_equal(one, default)


def test_full_search_ranges_does_not_depend_on_workers():
    one = full_search_ranges(Guide(GUIDE), sirna_seqs(), workers=1)
    four = full_search_ranges(Guide(GUIDE), sirna_seqs(), workers=4)
    default = full_search_ranges(Guide(GUIDE), sirna_seqs())
    pd.testing.assert_frame_equal(one, four)
    pd.testing.assert_frame_equal(one, default)


def test_full_search_columns_follow_the_options():
    options = FullSearchOptions(
        seeds=(Seed.MER7M8, Seed.FULL), nonseed_edit_distances=(2, 0), seed_edit_distances=(1,)
    )
    got = full_search(Guide(GUIDE), seqs_of("ACGT"), options, group="treated")
    assert got.columns.tolist() == [
        "seq_id",
        "gene_id",
        "mer7m8.mm1",
        "Full.mm2",
        "Full.mm0",
        "group",
        "best_match",
    ]
    assert got["group"].tolist() == ["treated"]
    assert got[["mer7m8.mm1", "Full.mm2", "Full.mm0"]].dtypes.tolist() == [np.dtype("int64")] * 3


def test_full_search_default_columns_and_group():
    got = full_search(Guide(GUIDE), seqs_of("ACGT"))
    assert got.columns[:2].tolist() == ["seq_id", "gene_id"]
    assert got.columns[-2:].tolist() == ["group", "best_match"]
    assert got.shape[1] == 2 + 5 * 3 + 2 * 4 + 2
    assert got["group"].tolist() == ["siRNA"]


def test_full_search_keeps_the_sequence_id_and_the_gene_id():
    seqs = seqs_of("CTCTATA", "AAAA", gene_ids=["g1", "g1"])
    got = full_search(Guide(GUIDE), seqs, FullSearchOptions(seeds=(Seed.MER7M8,)))
    assert got[["seq_id", "gene_id", "mer7m8.mm0", "best_match"]].to_numpy().tolist() == [
        ["t0", "g1", 1, "7mer-m8"],
        ["t1", "g1", 0, "No.matches"],
    ]


@pytest.mark.parametrize("fixed_pattern", [True, False])
@pytest.mark.parametrize("fixed_subject", [True, False])
def test_full_search_passes_the_search_flags(fixed_pattern, fixed_subject):
    seqs = seqs_of("CTCTATN", "CTCTRTA")
    options = FullSearchOptions(
        seeds=(Seed.MER7M8,),
        seed_edit_distances=(0, 1),
        fixed_pattern=fixed_pattern,
        fixed_subject=fixed_subject,
    )
    got = full_search(Guide(GUIDE), seqs, options)
    target = Guide(GUIDE).target(Seed.MER7M8)
    for distance in (0, 1):
        search = SearchOptions(
            max_mismatch=distance, fixed_pattern=fixed_pattern, fixed_subject=fixed_subject
        )
        assert got[f"mer7m8.mm{distance}"].tolist() == count_matches(target, seqs, search).tolist()


def test_full_search_of_an_empty_set_has_the_columns_and_no_rows():
    got = full_search(Guide(GUIDE), SequenceSet.from_records([]))
    assert len(got) == 0
    assert got.columns[:3].tolist() == ["seq_id", "gene_id", "Full.mm0"]
    assert got.columns[-2:].tolist() == ["group", "best_match"]


def test_full_search_of_an_empty_set_has_the_dtypes_of_a_set_with_sequences():
    empty = full_search(Guide(GUIDE), SequenceSet.from_records([]))
    full = full_search(Guide(GUIDE), seqs_of(SUBJECT_WITH_EVERY_SEED))
    assert empty.dtypes.to_dict() == full.dtypes.to_dict()
    assert str(empty["seq_id"].dtype) == TEXT_DTYPE
    assert str(empty["gene_id"].dtype) == TEXT_DTYPE


def test_full_search_ranges_order_follows_the_grid_of_r():
    options = FullSearchRangesOptions(
        edit_distances=(0, 1, 2), seeds=(Seed.MER7M8, Seed.FULL), max_seed_edit_distance=1
    )
    got = full_search_ranges(Guide(GUIDE), seqs_of(SUBJECT_WITH_EVERY_SEED), options)
    assert list(dict.fromkeys(got["query"])) == [
        "mer7m8.mm0.indelTRUE",
        "mer7m8.mm0.indelFALSE",
        "Full.mm0.indelTRUE",
        "Full.mm0.indelFALSE",
        "mer7m8.mm1.indelTRUE",
        "mer7m8.mm1.indelFALSE",
        "Full.mm1.indelTRUE",
        "Full.mm1.indelFALSE",
        "Full.mm2.indelTRUE",
        "Full.mm2.indelFALSE",
    ]


def test_full_search_ranges_with_a_match_for_every_search_D022():
    got = full_search_ranges(Guide(GUIDE), seqs_of(SUBJECT_WITH_EVERY_SEED))
    seed_names = {seed.value for seed in DEFAULT_SEEDS}
    expected = {
        f"{seed}.mm{d}.indel{flag}"
        for seed in seed_names
        for d in range(5 if seed in {"Full", "18mer", "15mer"} else 2)
        for flag in ("TRUE", "FALSE")
    }
    assert set(got["query"]) == expected


def test_full_search_ranges_row_fields():
    options = FullSearchRangesOptions(
        indel_options=(False,), edit_distances=(0,), seeds=(Seed.FULL, Seed.MER7M8)
    )
    got = full_search_ranges(Guide(GUIDE), seqs_of(SUBJECT_WITH_EVERY_SEED), options, group="g")
    assert got.columns.tolist() == [
        "seq_id",
        "gene_id",
        "start",
        "end",
        "width",
        "query",
        "target_seq",
        "edit_distance",
        "indels",
        "group",
    ]
    assert got.to_numpy().tolist() == [
        ["t0", "t0", 4, 25, 22, "Full.mm0.indelFALSE", FULL_TARGET, 0, False, "g"],
        ["t0", "t0", 19, 25, 7, "mer7m8.mm0.indelFALSE", "CTCTATA", 0, False, "g"],
    ]
    assert got.dtypes.astype(str).to_dict() == {
        "seq_id": TEXT_DTYPE,
        "gene_id": TEXT_DTYPE,
        "start": "int64",
        "end": "int64",
        "width": "int64",
        "query": TEXT_DTYPE,
        "target_seq": TEXT_DTYPE,
        "edit_distance": "int64",
        "indels": "bool",
        "group": TEXT_DTYPE,
    }


def test_full_search_ranges_uses_the_distance_and_the_flag_of_each_row():
    options = FullSearchRangesOptions(edit_distances=(0, 1), seeds=(Seed.MER18,))
    got = full_search_ranges(Guide(GUIDE), seqs_of(SUBJECT_WITH_EVERY_SEED), options)
    target = Guide(GUIDE).target(Seed.MER18).sequence
    assert set(got["target_seq"]) == {target}
    pairs = got[["query", "edit_distance", "indels"]].drop_duplicates().to_numpy().tolist()
    assert pairs == [
        ["18mer.mm0.indelTRUE", 0, True],
        ["18mer.mm0.indelFALSE", 0, False],
        ["18mer.mm1.indelTRUE", 1, True],
        ["18mer.mm1.indelFALSE", 1, False],
    ]


def test_full_search_ranges_without_a_search_in_the_grid_is_an_empty_frame():
    options = FullSearchRangesOptions(edit_distances=(2, 3), seeds=(Seed.MER8,))
    got = full_search_ranges(Guide(GUIDE), seqs_of(SUBJECT_WITH_EVERY_SEED), options)
    assert len(got) == 0
    assert got.columns.tolist()[-3:] == ["edit_distance", "indels", "group"]
    assert got["edit_distance"].dtype == np.int64
    assert got["indels"].dtype == np.bool_


def test_full_search_ranges_without_any_match_is_an_empty_frame_D022():
    got = full_search_ranges(Guide(GUIDE), seqs_of("A" * 40))
    assert len(got) == 0
    assert got.columns.tolist()[-3:] == ["edit_distance", "indels", "group"]


def test_full_search_ranges_of_an_empty_set_is_an_empty_frame():
    got = full_search_ranges(Guide(GUIDE), SequenceSet.from_records([]))
    assert len(got) == 0
    assert got.columns[:2].tolist() == ["seq_id", "gene_id"]


def test_full_search_ranges_of_an_empty_set_has_the_dtypes_of_a_set_with_matches():
    empty = full_search_ranges(Guide(GUIDE), SequenceSet.from_records([]))
    full = full_search_ranges(Guide(GUIDE), seqs_of(SUBJECT_WITH_EVERY_SEED))
    assert empty.dtypes.to_dict() == full.dtypes.to_dict()


def test_full_search_ranges_agree_with_the_counts_of_full_search():
    seqs = sirna_seqs()
    ranges = full_search_ranges(Guide(GUIDE), seqs)
    counts = full_search(Guide(GUIDE), seqs, FullSearchOptions(indels=False))
    by_query = ranges.groupby(["query", "seq_id"]).size()
    for column in ("Full.mm4", "18mer.mm4", "mer7m8.mm1", "mer6.mm1"):
        found = by_query.get(f"{column}.indelFALSE", pd.Series(dtype="int64"))
        expected = counts.set_index("seq_id")[column]
        assert found.reindex(expected.index, fill_value=0).tolist() == expected.tolist()


def test_the_package_exports_the_full_search_names():
    for module in (seedmatch, seedmatch.report):
        for name in (
            "FullSearchOptions",
            "FullSearchRangesOptions",
            "best_match",
            "full_search",
            "full_search_ranges",
        ):
            assert callable(getattr(module, name))
