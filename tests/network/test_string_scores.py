import itertools
import math
import warnings

import numpy as np
import pandas as pd
import pytest
from annotation_cases import fixture_reference
from helpers import DATA, exactly, golden

import seedmatch
from seedmatch.errors import InputError, SeedMatchWarning
from seedmatch.network import Aggregator, read_string_links, string_gene_network, string_scores

SAMPLE = DATA / "string_links_sample.txt"
FEATURES = "network/string_features.tsv"
COLUMN = "mer7m8"
NAN = math.nan
EXPECTED = {
    "gene_id": [1, 2, 4, 5, 6, 7, 8],
    "seed_count": [2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "neighbors": [2, 2, 1, 2, 2, 2, 3],
    Aggregator.MEAN: [1.0, 1.5, 1.0, 0.5, 0.75, 0.0, 0.5],
    Aggregator.DEGREE_IMPORTANCE: [1.0, 1.5, 0.5, 0.25, 0.375, 0.0, 0.25],
}


def gene(number):
    return f"ENSSYNG{number:011d}"


def sample_network(**options):
    """Return the edges of the STRING sample, mapped to the genes of the fixture EnsDb."""
    links = read_string_links(SAMPLE, **options).links
    with fixture_reference() as reference:
        return string_gene_network(links, reference).edges


def sample_de(**changes):
    de = pd.DataFrame({"gene_id": [gene(1), gene(2), gene(2), gene(99)], COLUMN: [2, 1, 3, 5]})
    return de.assign(**changes)


def counts_of(*rows):
    return pd.DataFrame({"gene_id": [gene(n) for n, _ in rows], COLUMN: [c for _, c in rows]})


def expected_frame(aggregator):
    return pd.DataFrame(
        {
            "gene_id": pd.Series([gene(n) for n in EXPECTED["gene_id"]], dtype=object),
            "seed_count": pd.Series(EXPECTED["seed_count"], dtype=np.float64),
            "neighbors": pd.Series(EXPECTED["neighbors"], dtype=np.int64),
            "score": pd.Series(EXPECTED[aggregator], dtype=np.float64),
        }
    )


def pair_table(*pairs):
    return pd.DataFrame(
        {"gene1": [first for first, _ in pairs], "gene2": [second for _, second in pairs]}
    )


def symmetric(*pairs):
    return pair_table(
        *[row for first, second in pairs for row in ((first, second), (second, first))]
    )


def test_string_scores_give_the_mean_of_the_neighbors():
    got = string_scores(sample_network(percentile=0), sample_de(), column=COLUMN)
    pd.testing.assert_frame_equal(got, expected_frame(Aggregator.MEAN), check_exact=True)
    assert [str(dtype) for dtype in got.dtypes] == ["object", "float64", "int64", "float64"]
    assert isinstance(got.index, pd.RangeIndex)


def test_string_scores_give_the_degree_importance():
    got = string_scores(
        sample_network(percentile=0),
        sample_de(),
        column=COLUMN,
        aggregator=Aggregator.DEGREE_IMPORTANCE,
    )
    pd.testing.assert_frame_equal(
        got, expected_frame(Aggregator.DEGREE_IMPORTANCE), check_exact=True
    )


def test_string_scores_use_the_mean_for_the_default_aggregator():
    network = sample_network(percentile=0)
    default = string_scores(network, sample_de(), column=COLUMN)
    explicit = string_scores(network, sample_de(), column=COLUMN, aggregator=Aggregator.MEAN)
    pd.testing.assert_frame_equal(default, explicit, check_exact=True)


@pytest.mark.parametrize("counts", [[1, 3], [3, 1], [3, 3]], ids=["first low", "first high", "tie"])
def test_string_scores_take_the_highest_count_of_a_gene_D062(counts):
    de = counts_of((2, counts[0]), (2, counts[1]))
    got = string_scores(sample_network(percentile=0), de, column=COLUMN).set_index("gene_id")
    assert got.loc[gene(2), "seed_count"] == 3.0
    r = golden(FEATURES).query("case == 'repeated'").set_index("name")
    assert float(r.loc[gene(2), "feature"]) == 1.0


# The value, its text in the message, and the case, the gene, and the value of the R golden.
BAD_COUNTS = [
    (NAN, "NaN", "na_count", 6, "NA"),
    (-1.0, "-1", "bad_counts", 1, "-1"),
    (math.inf, "inf", "bad_counts", 4, "Inf"),
    (-math.inf, "-inf", None, None, None),
]


@pytest.mark.parametrize(
    "bad", BAD_COUNTS, ids=["NaN", "negative", "infinite", "negative infinite"]
)
@pytest.mark.parametrize("target", [2, 99], ids=["in the network", "outside the network"])
def test_string_scores_reject_a_bad_count_D062(bad, target):
    value, text, r_case, r_gene, r_text = bad
    de = counts_of((1, 2.0), (target, value))
    message = (
        f"column 'mer7m8' of de holds {text} for gene {gene(target)!r}; "
        "expected a finite count of 0 or more"
    )
    with pytest.raises(InputError, match=exactly(message)):
        string_scores(sample_network(percentile=0), de, column=COLUMN)
    if r_case is not None:
        r = golden(FEATURES).query(f"case == '{r_case}'").set_index("name")
        assert r.loc[gene(r_gene), "feature"] == r_text


def test_string_scores_reject_a_nullable_count_that_is_missing_D062():
    de = counts_of((1, 2), (2, 1)).astype({COLUMN: "Int64"})
    de.loc[1, COLUMN] = pd.NA
    message = (
        f"column 'mer7m8' of de holds NaN for gene {gene(2)!r}; "
        "expected a finite count of 0 or more"
    )
    with pytest.raises(InputError, match=exactly(message)):
        string_scores(sample_network(percentile=0), de, column=COLUMN)


def test_string_scores_warn_without_a_de_gene_in_the_network_D062():
    de = counts_of((99, 5))
    message = (
        "no gene of de is a gene of the network; each seed count is 0. "
        "Expected the gene IDs of the network, for example 'ENSSYNG00000000001'"
    )
    with pytest.warns(SeedMatchWarning, match=exactly(message)) as caught:
        got = string_scores(sample_network(percentile=0), de, column=COLUMN)
    assert len(caught) == 1
    assert caught[0].filename == __file__
    assert got["seed_count"].tolist() == [0.0] * 7
    assert got["score"].tolist() == [0.0] * 7
    r = golden(FEATURES).query("case == 'none_in_network'")
    assert (r["feature"].astype(float) == 0).all()


def test_string_scores_warn_for_a_de_table_without_a_row_D062():
    de = counts_of().astype({COLUMN: "int64"})
    with pytest.warns(SeedMatchWarning, match="^no gene of de is a gene of the network"):
        got = string_scores(sample_network(percentile=0), de, column=COLUMN)
    assert got["score"].tolist() == [0.0] * 7


def test_string_scores_name_the_first_gene_of_the_network_in_the_warning_D062():
    network = symmetric(("REG", "TGT"))
    with pytest.warns(SeedMatchWarning, match=r"for example 'REG'$"):
        string_scores(network, counts_of((99, 1)), column=COLUMN)


def test_string_scores_do_not_warn_with_a_de_gene_in_the_network():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        string_scores(sample_network(percentile=0), counts_of((7, 0)), column=COLUMN)


def test_string_scores_ignore_a_de_gene_outside_the_network():
    network = sample_network(percentile=0)
    without = sample_de().iloc[:3]
    pd.testing.assert_frame_equal(
        string_scores(network, sample_de(), column=COLUMN),
        string_scores(network, without, column=COLUMN),
        check_exact=True,
    )


def test_string_scores_accept_the_dtypes_of_a_de_table():
    network = sample_network(percentile=0)
    expected = string_scores(network, sample_de(), column=COLUMN)
    for changes in (
        {COLUMN: sample_de()[COLUMN].astype("float64")},
        {COLUMN: sample_de()[COLUMN].astype("Int64")},
        {COLUMN: sample_de()[COLUMN].astype("uint8")},
        {"gene_id": sample_de()["gene_id"].astype("category")},
        {"gene_id": sample_de()["gene_id"].astype("string")},
    ):
        got = string_scores(network, sample_de(**changes), column=COLUMN)
        pd.testing.assert_frame_equal(got, expected, check_exact=True)


def test_string_scores_use_the_column_that_the_caller_names():
    de = sample_de().rename(columns={COLUMN: "overall_matches"})
    got = string_scores(sample_network(percentile=0), de, column="overall_matches")
    pd.testing.assert_frame_equal(got, expected_frame(Aggregator.MEAN), check_exact=True)


def test_string_scores_keep_the_gene_order_of_the_network():
    # The genes come in the order of the first occurrence in gene1, and then in gene2.
    network = pair_table(("B", "A"), ("A", "B"), ("A", "C"), ("C", "A"), ("B", "D"), ("D", "B"))
    got = string_scores(network, pd.DataFrame({"gene_id": ["A"], COLUMN: [1]}), column=COLUMN)
    assert got["gene_id"].tolist() == ["B", "A", "C", "D"]
    network = pair_table(("A", "B"), ("C", "A"), ("B", "A"), ("A", "C"))
    got = string_scores(network, pd.DataFrame({"gene_id": ["A"], COLUMN: [1]}), column=COLUMN)
    assert got["gene_id"].tolist() == ["A", "C", "B"]


def test_string_scores_accept_a_network_that_does_not_come_from_string():
    network = symmetric(("A", "B"), ("B", "C")).assign(combined_score=[900, 900, 800, 800])
    de = pd.DataFrame({"gene_id": ["A", "C"], COLUMN: [4, 2]})
    got = string_scores(network, de, column=COLUMN)
    assert got["gene_id"].tolist() == ["A", "B", "C"]
    assert got["neighbors"].tolist() == [1, 2, 1]
    assert got["score"].tolist() == [(0 + 4) / 2, ((4 + 2) / 2 + 0) / 2, (0 + 2) / 2]


def random_case(rng):
    genes = [f"G{n}" for n in range(rng.integers(2, 25))]
    pairs = [pair for pair in itertools.combinations(genes, 2) if rng.random() < 0.2] or [
        (genes[0], genes[1])
    ]
    rows = [row for first, second in pairs for row in ((first, second), (second, first))]
    network = pair_table(*[rows[i] for i in rng.permutation(len(rows))])
    members = [*genes, "OUTSIDE"]
    chosen = rng.choice(members, size=rng.integers(0, 3 * len(genes)))
    de = pd.DataFrame({"gene_id": chosen, COLUMN: rng.integers(0, 6, size=len(chosen))})
    return network, de


def loop_scores(network, de, aggregator):
    """Return the scores of the definition, one gene at a time."""
    genes = list(dict.fromkeys([*network["gene1"], *network["gene2"]]))
    neighbors = {name: [] for name in genes}
    for first, second in zip(network["gene1"], network["gene2"], strict=True):
        neighbors[first].append(second)
    seeds = dict.fromkeys(genes, 0.0)
    for name, count in zip(de["gene_id"], de[COLUMN], strict=True):
        if name in seeds:
            seeds[name] = max(seeds[name], float(count))
    weight = {
        name: 1.0 if aggregator is Aggregator.MEAN else 1.0 / len(neighbors[name]) for name in genes
    }
    scores = [
        (sum(weight[n] * seeds[n] for n in neighbors[name]) / len(neighbors[name]) + seeds[name])
        / 2
        for name in genes
    ]
    return genes, [seeds[name] for name in genes], [len(neighbors[name]) for name in genes], scores


@pytest.mark.parametrize("aggregator", list(Aggregator))
def test_string_scores_match_a_loop_over_the_genes(aggregator):
    rng = np.random.default_rng(7)
    for _ in range(40):
        network, de = random_case(rng)
        genes, seeds, degrees, scores = loop_scores(network, de, aggregator)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SeedMatchWarning)
            got = string_scores(network, de, column=COLUMN, aggregator=aggregator)
        assert got["gene_id"].tolist() == genes
        assert got["seed_count"].tolist() == seeds
        assert got["neighbors"].tolist() == degrees
        assert got["score"].tolist() == pytest.approx(scores, rel=1e-12, abs=0)


def test_string_scores_check_each_row_and_not_the_highest_count_of_a_gene_D062():
    de = counts_of((2, 3), (2, -1))
    message = (
        f"column 'mer7m8' of de holds -1 for gene {gene(2)!r}; expected a finite count of 0 or more"
    )
    with pytest.raises(InputError, match=exactly(message)):
        string_scores(sample_network(percentile=0), de, column=COLUMN)


NETWORK_ERRORS = [
    ([1, 2], "network has type list; expected a pandas.DataFrame"),
    (
        symmetric(("G1", "G4")).drop(columns="gene2"),
        "network table has missing columns: ['gene2']; expected ['gene1', 'gene2']",
    ),
    (
        pair_table((5, "G4"), ("G4", "G1")),
        "network table column 'gene1' holds 5 in row 0; expected a str in each row",
    ),
    (
        pair_table(("G1", NAN), ("G4", "G1")),
        "network table column 'gene2' holds nan in row 0; expected a str in each row",
    ),
    (
        pd.DataFrame({"gene1": [], "gene2": []}),
        "the network table has no row; expected at least one gene pair",
    ),
    (
        pair_table(("G1", "G4"), ("G4", "G1"), ("G1", "G1")),
        "network table: the pair ('G1', 'G1') joins a gene to itself; expected two different genes",
    ),
    (
        pair_table(("G1", "G4"), ("G4", "G1"), ("G1", "G4")),
        "network table: the pair ('G1', 'G4') occurs 2 times; expected one row for each pair",
    ),
    (
        pair_table(("G1", "G4"), ("G4", "G1"), ("G1", "G4"), ("G1", "G4")),
        "network table: the pair ('G1', 'G4') occurs 3 times; expected one row for each pair",
    ),
    (
        pair_table(("G1", "G4")),
        (
            "network table: the pair ('G1', 'G4') has no reverse row ('G4', 'G1'); "
            "expected both orders of each pair"
        ),
    ),
    (
        pair_table(("G1", "G4"), ("G4", "G1"), ("G4", "G5")),
        (
            "network table: the pair ('G4', 'G5') has no reverse row ('G5', 'G4'); "
            "expected both orders of each pair"
        ),
    ),
    (
        # The key of ('X', 'W') is below the key of ('X', 'Y'), but ('X', 'Y') has the first row.
        pair_table(("X", "Y"), ("X", "W"), ("W", "X"), ("Y", "X"), ("X", "Y"), ("X", "W")),
        "network table: the pair ('X', 'Y') occurs 2 times; expected one row for each pair",
    ),
    (
        # A gene with itself comes before a repeated pair and a pair without a reverse row.
        pair_table(("A", "B"), ("A", "B"), ("C", "C"), ("D", "E")),
        "network table: the pair ('C', 'C') joins a gene to itself; expected two different genes",
    ),
    (
        pair_table(("A", "B"), ("B", "A"), ("C", "C"), ("D", "D")),
        "network table: the pair ('C', 'C') joins a gene to itself; expected two different genes",
    ),
    (
        # The first row without a reverse row is in row order, and not in the order of the keys.
        pair_table(("A", "B"), ("B", "A"), ("D", "C"), ("B", "C")),
        (
            "network table: the pair ('D', 'C') has no reverse row ('C', 'D'); "
            "expected both orders of each pair"
        ),
    ),
]
NETWORK_IDS = [
    "list",
    "no gene2",
    "number",
    "NaN",
    "no row",
    "self pair",
    "pair 2 times",
    "pair 3 times",
    "no reverse row",
    "no reverse row after a pair",
    "the first repeated pair by row",
    "the self pair first",
    "the first self pair by row",
    "the first missing reverse by row",
]


@pytest.mark.parametrize(("network", "message"), NETWORK_ERRORS, ids=NETWORK_IDS)
def test_string_scores_reject_a_bad_network(network, message):
    de = pd.DataFrame({"gene_id": ["G1"], COLUMN: [1]})
    with pytest.raises(InputError, match=exactly(message)):
        string_scores(network, de, column=COLUMN)


def test_string_scores_reject_a_network_with_a_column_twice():
    network = pd.concat([symmetric(("G1", "G4")), symmetric(("G1", "G4"))["gene1"]], axis=1)
    message = "network table has 2 columns 'gene1'; expected one column for each name"
    with pytest.raises(InputError, match=exactly(message)):
        string_scores(network, counts_of((1, 1)), column=COLUMN)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"column": 5}, "column 5 has type int; expected a str"),
        (
            {"column": "overall_matches"},
            (
                "DE table has missing columns: ['overall_matches']; expected "
                "['gene_id', 'overall_matches']"
            ),
        ),
        ({"de": [1, 2]}, "de has type list; expected a pandas.DataFrame"),
        (
            {"de": sample_de().assign(**{COLUMN: ["a", "b", "c", "d"]})},
            "column 'mer7m8' has dtype object; expected a numeric dtype",
        ),
        (
            {"de": sample_de().assign(gene_id=[1, 2, 3, 4])},
            "'gene_id' has dtype int64; expected text values. Convert the column with astype(str)",
        ),
        (
            {"de": sample_de().assign(**{COLUMN: [True, False, True, False]})},
            "column 'mer7m8' has dtype bool; expected a numeric dtype",
        ),
        (
            {"aggregator": "mean"},
            (
                "aggregator 'mean' is not valid; expected one of Aggregator.MEAN, "
                "Aggregator.DEGREE_IMPORTANCE"
            ),
        ),
    ],
    ids=[
        "number column",
        "no column",
        "list de",
        "text count",
        "number gene",
        "bool count",
        "aggregator text",
    ],
)
def test_string_scores_reject_a_bad_argument(kwargs, message):
    arguments = {"de": sample_de(), "column": COLUMN, **kwargs}
    with pytest.raises(InputError, match=exactly(message)):
        string_scores(sample_network(percentile=0), **arguments)


def test_string_scores_do_not_change_the_tables():
    network = sample_network(percentile=0)
    de = sample_de()
    network_before, de_before = network.copy(), de.copy()
    string_scores(network, de, column=COLUMN)
    pd.testing.assert_frame_equal(network, network_before)
    pd.testing.assert_frame_equal(de, de_before)


def test_string_scores_is_in_the_top_level_package():
    assert seedmatch.string_scores is string_scores
    assert "string_scores" in seedmatch.__all__
