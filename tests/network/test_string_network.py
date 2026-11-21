import dataclasses
import re

import pandas as pd
import pytest
from annotation_cases import ANNOTATION, ENSDB, copy_ensdb, fixture_reference
from helpers import DATA, exactly, golden

import seedmatch
from seedmatch.annotation.reference import Reference
from seedmatch.errors import InputError
from seedmatch.network import GeneNetwork, read_string_links, string_gene_network

SAMPLE = DATA / "string_links_sample.txt"
NODES = "network/string_nodes.tsv"
PPI = "network/string_ppi.tsv"
LINK_COLUMNS = ["protein1", "protein2", "combined_score"]
NO_PREFIX = "has no taxon prefix; expected the form '<taxon>.<protein ID>'"


def gene(number):
    return f"ENSSYNG{number:011d}"


def string_id(number):
    return f"32630.ENSSYNP{number:011d}"


def edges(*rows):
    return [[gene(first), gene(second), score] for first, second, score in rows]


def edges_of(network):
    return network.edges.to_numpy().tolist()


def links_of(*rows):
    return pd.DataFrame(
        [(string_id(first), string_id(second), score) for first, second, score in rows],
        columns=LINK_COLUMNS,
    )


def r_pairs(case):
    table = golden(PPI)
    rows = table[table["case"] == case]
    return list(zip(rows["gene1"], rows["gene2"], strict=True))


def pairs_of(network):
    return list(zip(network.edges["gene1"], network.edges["gene2"], strict=True))


def sample_links(**options):
    return read_string_links(SAMPLE, **options).links


def network_of(**options):
    with fixture_reference() as reference:
        return string_gene_network(sample_links(**options), reference)


P50_EDGES = edges(
    (1, 4, 900),
    (1, 5, 850),
    (2, 6, 750),
    (2, 8, 700),
    (4, 1, 900),
    (5, 1, 850),
    (6, 2, 750),
    (8, 2, 700),
)
P0_EDGES = edges(
    (1, 4, 900),
    (1, 5, 850),
    (2, 6, 750),
    (2, 8, 700),
    (4, 1, 900),
    (5, 1, 850),
    (5, 7, 600),
    (6, 2, 750),
    (6, 8, 500),
    (7, 5, 600),
    (7, 8, 600),
    (8, 2, 700),
    (8, 6, 500),
    (8, 7, 600),
)


def test_string_gene_network_matches_r_with_one_protein_for_each_gene():
    got = network_of()
    assert edges_of(got) == edges((1, 4, 900), (4, 1, 900))
    assert pairs_of(got) == r_pairs("p95")
    assert got.unmapped_proteins == ()


def test_string_gene_network_gives_the_table_of_the_documentation():
    got = network_of(percentile=0)
    assert isinstance(got, GeneNetwork)
    assert list(got.edges.columns) == ["gene1", "gene2", "combined_score"]
    assert [str(dtype) for dtype in got.edges.dtypes] == ["object", "object", "int64"]
    assert isinstance(got.edges.index, pd.RangeIndex)
    assert type(got.unmapped_proteins) is tuple
    with pytest.raises(dataclasses.FrozenInstanceError):
        got.unmapped_proteins = ()


def test_string_gene_network_joins_the_proteins_of_a_gene_D060():
    got = network_of(percentile=0.5)
    assert edges_of(got) == P50_EDGES
    r_set = set(r_pairs("p50"))
    port_set = set(pairs_of(got))
    assert r_set < port_set
    assert port_set - r_set == {
        (gene(1), gene(5)),
        (gene(5), gene(1)),
        (gene(2), gene(8)),
        (gene(8), gene(2)),
    }


def test_string_gene_network_gives_a_pair_the_highest_score_of_its_links_D060():
    # P1 and P2 are proteins of G1. P1-P2 joins one gene and gives no pair. P1-P1 and P6-P6 link a
    # protein with itself. R keeps such a link as a loop of the gene, and the port drops it (D-060).
    links = links_of((1, 1, 950), (1, 6, 700), (2, 6, 900), (6, 1, 500), (1, 2, 800), (6, 6, 950))
    with fixture_reference() as reference:
        got = string_gene_network(links, reference)
    assert edges_of(got) == edges((1, 4, 900), (4, 1, 900))


def test_string_gene_network_drops_the_unmapped_proteins_D061():
    got = network_of(percentile=0)
    assert edges_of(got) == P0_EDGES
    names = golden(NODES)
    r_proteins = [
        name for name in names[names["case"] == "p0"]["name"] if name.startswith("ENSSYNP")
    ]
    assert r_proteins == ["ENSSYNP00000000098", "ENSSYNP00000000099"]
    assert got.unmapped_proteins == tuple(f"32630.{name}" for name in r_proteins)


def test_string_gene_network_gives_no_pair_for_a_link_of_an_unmapped_protein_D061():
    # P98 comes first and has no gene. The first gene of the network is G1, with P1.
    links = links_of((98, 6, 800), (1, 12, 700))
    with fixture_reference() as reference:
        got = string_gene_network(links, reference)
    assert edges_of(got) == edges((1, 8, 700), (8, 1, 700))
    assert got.unmapped_proteins == (string_id(98),)


def test_string_gene_network_lists_a_protein_of_a_transcript_without_a_gene_as_unmapped(tmp_path):
    path = copy_ensdb(tmp_path, "UPDATE tx SET gene_id = NULL WHERE tx_id = 'ENSSYNT00000000012'")
    links = links_of((1, 6, 900), (4, 12, 700))
    with Reference.from_paths(ensdb=path, genome=ANNOTATION / "synthetic.2bit") as reference:
        got = string_gene_network(links, reference)
    assert edges_of(got) == edges((1, 4, 900), (4, 1, 900))
    assert got.unmapped_proteins == (string_id(12),)


def test_string_gene_network_lists_the_unmapped_proteins_in_the_order_of_the_links():
    # The order is P99, P98: the column protein1 comes before the column protein2.
    links = links_of((1, 98, 800), (99, 6, 700), (1, 6, 900))
    with fixture_reference() as reference:
        got = string_gene_network(links, reference)
    assert edges_of(got) == edges((1, 4, 900), (4, 1, 900))
    assert got.unmapped_proteins == (string_id(99), string_id(98))


def test_string_gene_network_keeps_the_text_after_the_first_dot_of_a_protein_ID_D061(tmp_path):
    path = copy_ensdb(
        tmp_path,
        "UPDATE protein SET protein_id = 'ENSSYNP00000000001.1' "
        "WHERE protein_id = 'ENSSYNP00000000001'",
    )
    links = links_of((1, 6, 900))
    links.loc[0, "protein1"] = "32630.ENSSYNP00000000001.1"
    with Reference.from_paths(ensdb=path, genome=ANNOTATION / "synthetic.2bit") as reference:
        got = string_gene_network(links, reference)
    assert edges_of(got) == edges((1, 4, 900), (4, 1, 900))
    assert got.unmapped_proteins == ()


def test_string_gene_network_orders_the_rows_by_the_first_occurrence_of_a_gene():
    # The genes come in the order G8, G1, G2, G4, and so do the rows.
    links = links_of((12, 4, 700), (1, 6, 900))
    with fixture_reference() as reference:
        got = string_gene_network(links, reference)
    assert edges_of(got) == edges((8, 2, 700), (1, 4, 900), (2, 8, 700), (4, 1, 900))


def test_string_gene_network_maps_with_the_protein_table_D013(tmp_path):
    path = copy_ensdb(
        tmp_path,
        "UPDATE protein SET tx_id = 'ENSSYNT00000000012' WHERE protein_id = 'ENSSYNP00000000006'",
    )
    with Reference.from_paths(ensdb=path, genome=ANNOTATION / "synthetic.2bit") as reference:
        got = string_gene_network(sample_links(), reference)
    assert edges_of(got) == edges((1, 8, 900), (8, 1, 900))


def test_string_gene_network_rejects_links_without_a_gene_pair_D061():
    links = links_of((98, 99, 450), (1, 2, 800))
    with fixture_reference() as reference:
        message = (
            f"no link joins two genes of the annotation {reference.annotation.source!r}; "
            "expected the STRING links of the species of the annotation"
        )
        with pytest.raises(InputError, match=exactly(message)):
            string_gene_network(links, reference)


@pytest.mark.parametrize("column", ["protein1", "protein2"])
@pytest.mark.parametrize("bad", ["ENSSYNP00000000001", ".ENSSYNP00000000001", "32630."])
def test_string_gene_network_rejects_a_protein_id_without_a_taxon_prefix_D061(column, bad):
    links = links_of((1, 6, 900))
    links.loc[0, column] = bad
    with fixture_reference() as reference:
        message = f"links table: protein ID {bad!r} {NO_PREFIX}"
        with pytest.raises(InputError, match=exactly(message)):
            string_gene_network(links, reference)


GTF_PAIRS = [("synthetic.gtf", "synthetic.2bit"), ("synthetic.gencode.gtf.gz", "synthetic.fa.gz")]


@pytest.mark.parametrize(("gtf", "genome"), GTF_PAIRS)
def test_string_gene_network_maps_with_the_protein_ids_of_a_gtf_D013(gtf, genome):
    expected = network_of()
    with Reference.from_gtf(ANNOTATION / gtf, ANNOTATION / genome) as reference:
        got = string_gene_network(sample_links(), reference)
    assert edges_of(got) == edges_of(expected) == edges((1, 4, 900), (4, 1, 900))
    assert got.unmapped_proteins == expected.unmapped_proteins


def test_string_gene_network_rejects_a_gtf_reference_without_protein_ids(tmp_path):
    path = tmp_path / "synthetic.gtf"
    text = (ANNOTATION / "synthetic.gtf").read_text()
    path.write_text(re.sub(r' protein_id "[^"]*";', "", text))
    with Reference.from_gtf(path, ANNOTATION / "synthetic.2bit") as gtf:
        message = (
            f"{gtf.annotation.source!r}: the annotation has no table 'protein'; "
            "expected an EnsDb file with this table, "
            "or a GTF with the attribute protein_id on its CDS lines"
        )
        with pytest.raises(InputError, match=exactly(message)):
            string_gene_network(sample_links(), gtf)


def scored_links(scores, dtype="int64"):
    """Return links of P1 and P6 with these scores, in the dtype."""
    return links_of(*((1, 6, 0) for _ in scores)).assign(
        combined_score=pd.Series(scores, dtype=dtype)
    )


def empty_links():
    return pd.DataFrame(
        {
            "protein1": pd.Series([], dtype=object),
            "protein2": pd.Series([], dtype=object),
            "combined_score": pd.Series([], dtype="int64"),
        }
    )


@pytest.mark.parametrize(
    ("links", "message"),
    [
        ([1, 2], "links has type list; expected a pandas.DataFrame"),
        (
            links_of((1, 6, 900)).drop(columns="combined_score"),
            (
                "links table has missing columns: ['combined_score']; "
                "expected ['protein1', 'protein2', 'combined_score']"
            ),
        ),
        (
            links_of((1, 6, 900)).assign(protein1=5),
            "links table column 'protein1' holds 5 in row 0; expected a str in each row",
        ),
        (
            links_of((1, 6, 900)).astype({"combined_score": "float64"}),
            "links table column 'combined_score' has dtype float64; expected an integer dtype",
        ),
        (
            scored_links([900, None], dtype="Int64"),
            (
                "links table column 'combined_score' has missing values; "
                "expected an integer in each row"
            ),
        ),
        (empty_links(), "the links table has no row; expected at least one link"),
        (
            scored_links([900, 800, 700, 1001]),
            (
                "links table column 'combined_score' holds 1001 in row 3; "
                "expected an integer from 0 to 1000"
            ),
        ),
        (
            scored_links([900, -1, 700, 1001]),
            (
                "links table column 'combined_score' holds -1 in row 1; "
                "expected an integer from 0 to 1000"
            ),
        ),
        (
            scored_links([900, 2**63], dtype="uint64"),
            (
                "links table column 'combined_score' holds 9223372036854775808 in row 1; "
                "expected an integer from 0 to 1000"
            ),
        ),
        (
            scored_links([900, 1001], dtype="Int64"),
            (
                "links table column 'combined_score' holds 1001 in row 1; "
                "expected an integer from 0 to 1000"
            ),
        ),
    ],
    ids=[
        "list",
        "no_score",
        "number_protein",
        "float_score",
        "missing_score",
        "no_row",
        "score_above_the_maximum",
        "score_below_zero",
        "uint64_score",
        "nullable_score",
    ],
)
def test_string_gene_network_rejects_bad_links(links, message):
    with fixture_reference() as reference, pytest.raises(InputError, match=exactly(message)):
        string_gene_network(links, reference)


@pytest.mark.parametrize("dtype", ["int64", "uint16", "uint64", "Int64"])
def test_string_gene_network_accepts_the_scores_from_0_to_the_maximum(dtype):
    links = scored_links([0, 1000], dtype=dtype)
    with fixture_reference() as reference:
        got = string_gene_network(links, reference)
    assert edges_of(got) == edges((1, 4, 1000), (4, 1, 1000))


def test_string_gene_network_rejects_a_reference_of_another_type():
    links = sample_links()
    with pytest.raises(InputError, match=exactly("reference has type str; expected a Reference")):
        string_gene_network(links, str(ENSDB))


def test_string_gene_network_and_gene_network_are_in_the_top_level_package():
    assert seedmatch.string_gene_network is string_gene_network
    assert seedmatch.GeneNetwork is GeneNetwork
    assert {"string_gene_network", "GeneNetwork"} <= set(seedmatch.__all__)
