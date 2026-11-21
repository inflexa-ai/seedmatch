import gzip
import math
import os
import warnings

import numpy as np
import pandas as pd
import pytest
from helpers import DATA, assert_frame_matches_golden, exactly, golden, golden_note

import seedmatch
from seedmatch.errors import InputError, SeedMatchWarning
from seedmatch.network import Aggregator, read_trrust, trrust_scores
from seedmatch.network.messages import mean_messages, pair_keys, sender_weights

SAMPLE = DATA / "trrust_sample.tsv"
ROWS = "network/trrust_rows.tsv"
MATCHES = "network/trrust_matches.tsv"
LINE_FORMAT = "the regulator, the target, the mode, and the PubMed IDs"
MODES = "Activation, Repression, or Unknown"
SEEDS = ["TFA", "GC", "TFC"]
GENES = ["TFA", "TFB", "GA", "TFC", "GB", "GC", "GD"]
NAN = math.nan
SCORE_COLUMNS = ["repression_score", "activation_score"]
PREFIX = ["TFA\tGA\tRepression\t1", "TFB\tGB\tActivation\t2", "TFC\tGC\tUnknown\t3"]


def write_lines(tmp_path, lines, name="x.tsv"):
    path = tmp_path / name
    path.write_text("".join(f"{line}\n" for line in lines), encoding="utf-8")
    return path


def failure(path, text):
    return exactly(f"{os.fspath(path)!r}: {text}")


def sample_lines():
    return SAMPLE.read_text().splitlines()


def note_lines(case):
    """Return the lines of the file that the note `case` of trrust_rows.tsv gives."""
    text = golden_note(
        ROWS, rf"^# {case}: the same call on a file with the lines (.*?) (?:stops|gives)"
    )
    return [line.replace("\\t", "\t") for line in text.split(" | ")]


def expected_scores(*, seed, repressors, repression, activators, activation):
    return pd.DataFrame(
        {
            "gene": pd.Series(GENES, dtype=object),
            "seed": pd.Series(seed, dtype=bool),
            "repressors": pd.Series(repressors, dtype=np.int64),
            "repression_score": pd.Series(repression, dtype=np.float64),
            "activators": pd.Series(activators, dtype=np.int64),
            "activation_score": pd.Series(activation, dtype=np.float64),
        }
    )


MEAN_SCORES = expected_scores(
    seed=[True, False, False, True, False, True, False],
    repressors=[1, 0, 1, 0, 1, 0, 1],
    repression=[0.5, NAN, 0.5, NAN, 0.5, NAN, 0.0],
    activators=[0, 0, 2, 1, 1, 0, 0],
    activation=[NAN, NAN, 0.25, 1.0, 0.5, NAN, NAN],
)
DEGREE_SCORES = MEAN_SCORES.assign(
    repression_score=[0.5, NAN, 0.25, NAN, 0.25, NAN, 0.0],
    activation_score=[NAN, NAN, 0.125, 0.75, 0.25, NAN, NAN],
)


def test_read_trrust_matches_r_without_a_header():
    got = read_trrust(SAMPLE)
    assert len(got) == 10
    assert list(got.columns) == ["source", "target", "edge_type", "evidence"]
    assert [str(dtype) for dtype in got.dtypes] == ["object"] * 4
    assert isinstance(got.index, pd.RangeIndex)
    assert_frame_matches_golden(got, ROWS)


def test_read_trrust_keeps_the_first_line_D057():
    got = read_trrust(SAMPLE)
    assert got.iloc[0].tolist() == ["TFA", "GA", "Repression", "1001"]
    assert golden_note(ROWS, r"^# header_true: .* gives (\d+) rows") == "9"
    assert golden_note(ROWS, r"^# header_true: .* before :30 are (.*)$") == (
        "TFA | GA | Repression | X1001"
    )


def test_read_trrust_keeps_quote_hash_and_na_text_D057(tmp_path):
    quote = read_trrust(write_lines(tmp_path, note_lines("quote")))
    assert quote["source"].tolist() == ["TF1", "TF'2", "TF3", "TF4"]
    assert quote["target"].tolist() == ["G1", "G2", "G3", "G'4"]
    hashed = read_trrust(write_lines(tmp_path, note_lines("hash")))
    assert hashed["source"].tolist() == ["TF1", "TF#2", "TF3"]
    na = read_trrust(write_lines(tmp_path, note_lines("na_text")))
    assert na["source"].tolist() == ["TF1", "NA", "TF3"]
    assert na["target"].tolist() == ["G1", "G2", "NA"]
    for got in (quote, hashed, na):
        assert not got.isna().to_numpy().any()


def test_read_trrust_keeps_a_double_quote_character_D057(tmp_path):
    lines = ['"TFA\tGA\tRepression\t1', 'TFB"\tGB\tActivation\t2', 'TFC\tG"C\tUnknown\t3']
    got = read_trrust(write_lines(tmp_path, lines))
    assert got["source"].tolist() == ['"TFA', 'TFB"', "TFC"]
    assert got["target"].tolist() == ["GA", "GB", 'G"C']


@pytest.mark.parametrize(
    ("lines", "message"),
    [
        (
            ["TFA\tGA\tRepression\t1", "TFB\tGB\tactivation\t2"],
            f"line 2: mode 'activation' is not valid; expected {MODES}",
        ),
        (
            ["TFA\tGA\tRepression\t1", "TFB\tGB\t\t2"],
            f"line 2: mode '' is not valid; expected {MODES}",
        ),
        (
            [*PREFIX, "\tGD\tRepression\t4"],
            "line 4: the regulator is empty; expected a gene symbol",
        ),
        (
            [*PREFIX, "TFD\t\tRepression\t4"],
            "line 4: the target is empty; expected a gene symbol",
        ),
        (
            ["TFA\tGA\tRepression\t1", "TFB\tGB\tActivation\t2", "TFC\tGC\tUnknown"],
            f"line 3 has 3 fields; expected 4: {LINE_FORMAT}",
        ),
        (
            ["TFA\tGA\tRepression\t1", "TFB\tGB\tActivation\t2", "TFC\tGC\tUnknown\t3\t4"],
            f"line 3 has 5 fields; expected 4: {LINE_FORMAT}",
        ),
        (
            ["TFA\tGA\tRepression\t1", "TFB\tGB\tActivation\t" + "9" * 1006],
            (
                "line 2: the line has more than 1024 characters; expected at most 1024, "
                "with the line end"
            ),
        ),
        (
            ["TFA\tGA\tRepression\t1", "TFB\tG\x00B\tActivation\t2"],
            "line 2: the line has a NUL character (U+0000); expected text without NUL characters",
        ),
        (
            ["TFA\tGA\tRepression\t1", "\t\t\t"],
            "line 2: the regulator is empty; expected a gene symbol",
        ),
        (
            ["TFA\tGA\tRepression\t1", "\t\tactivation\t"],
            "line 2: the regulator is empty; expected a gene symbol",
        ),
    ],
    ids=[
        "bad mode",
        "empty mode",
        "empty regulator",
        "empty target",
        "3 fields",
        "5 fields",
        "1025 characters",
        "NUL",
        "tabs only",
        "empty regulator before bad mode",
    ],
)
def test_read_trrust_rejects_a_bad_line_D057(tmp_path, lines, message):
    path = write_lines(tmp_path, lines)
    with pytest.raises(InputError, match=failure(path, message)):
        read_trrust(path)


def test_read_trrust_reads_a_line_of_the_limit_length(tmp_path):
    line = "TFB\tGB\tActivation\t" + "9" * 1005
    assert len(line) + 1 == 1024
    got = read_trrust(write_lines(tmp_path, ["TFA\tGA\tRepression\t1", line]))
    assert got["evidence"][1] == "9" * 1005


def test_read_trrust_accepts_an_empty_evidence(tmp_path):
    got = read_trrust(write_lines(tmp_path, ["TFA\tGA\tRepression\t"]))
    assert got.iloc[0].tolist() == ["TFA", "GA", "Repression", ""]


def test_read_trrust_rejects_a_url():
    message = "'https://x': the path is a URL; expected a local path"
    with pytest.raises(InputError, match=exactly(message)):
        read_trrust("https://x")


@pytest.mark.parametrize("text", ["", "\n", "\r\n\n\n"])
def test_read_trrust_rejects_an_empty_file(tmp_path, text):
    path = tmp_path / "x.tsv"
    path.write_text(text)
    with pytest.raises(
        InputError, match=failure(path, "the file has no line; expected TRRUST lines")
    ):
        read_trrust(path)


def test_read_trrust_rejects_text_that_is_not_utf8(tmp_path):
    path = tmp_path / "x.tsv"
    path.write_bytes(b"TFA\tGA\tRepression\t1\nTF\xff\tGB\tActivation\t2\n")
    message = r"x\.tsv': the file is not valid UTF-8 text .*; expected a text or gzip TRRUST file$"
    with pytest.raises(InputError, match=message):
        read_trrust(path)


def test_read_trrust_reads_gzip(tmp_path):
    path = tmp_path / "trrust_rawdata.human.tsv.gz"
    path.write_bytes(gzip.compress(SAMPLE.read_bytes()))
    pd.testing.assert_frame_equal(read_trrust(path), read_trrust(SAMPLE))


def test_read_trrust_accepts_a_str_path():
    pd.testing.assert_frame_equal(read_trrust(os.fspath(SAMPLE)), read_trrust(SAMPLE))


@pytest.mark.parametrize("ending", ["\r\n", "\r"])
def test_read_trrust_reads_each_line_end(tmp_path, ending):
    path = tmp_path / "x.tsv"
    path.write_bytes(SAMPLE.read_bytes().replace(b"\n", ending.encode()))
    pd.testing.assert_frame_equal(read_trrust(path), read_trrust(SAMPLE))


def test_read_trrust_skips_a_byte_order_mark(tmp_path):
    path = tmp_path / "x.tsv"
    path.write_bytes(b"\xef\xbb\xbf" + SAMPLE.read_bytes())
    pd.testing.assert_frame_equal(read_trrust(path), read_trrust(SAMPLE))


def test_read_trrust_skips_a_blank_line(tmp_path):
    lines = sample_lines()
    path = write_lines(tmp_path, ["", lines[0], "", "", *lines[1:], ""])
    pd.testing.assert_frame_equal(read_trrust(path), read_trrust(SAMPLE))
    bad = write_lines(tmp_path, ["", lines[0], "", "", "TFB\tGB\tactivation\t2"], "bad.tsv")
    message = f"line 5: mode 'activation' is not valid; expected {MODES}"
    with pytest.raises(InputError, match=failure(bad, message)):
        read_trrust(bad)


@pytest.mark.parametrize("chunk", [None, 1])
def test_read_trrust_shares_equal_texts(monkeypatch, chunk):
    if chunk is not None:
        monkeypatch.setattr("seedmatch.network.trrust._CHUNK_LINES", chunk)
    got = read_trrust(SAMPLE)
    assert got["source"][0] is got["source"][1]
    assert got["target"][0] is got["target"][5]
    assert got["edge_type"][0] is got["edge_type"][6]
    assert got["source"][4] is got["target"][0]


def test_read_trrust_reads_more_than_one_chunk(tmp_path, monkeypatch):
    expected = read_trrust(SAMPLE)
    monkeypatch.setattr("seedmatch.network.trrust._CHUNK_LINES", 3)
    pd.testing.assert_frame_equal(read_trrust(SAMPLE), expected)
    lines = sample_lines()
    lines[7] = "TFC\tTFC\tactivation\t1008"
    path = write_lines(tmp_path, lines)
    message = f"line 8: mode 'activation' is not valid; expected {MODES}"
    with pytest.raises(InputError, match=failure(path, message)):
        read_trrust(path)


def test_trrust_scores_use_the_regulators_D007():
    got = trrust_scores(read_trrust(SAMPLE), SEEDS)
    pd.testing.assert_frame_equal(got, MEAN_SCORES, check_exact=True)
    assert isinstance(got.index, pd.RangeIndex)
    r = golden(MATCHES).query("case == 'mean'").set_index("ENSG")
    assert float(r.loc["GA", "repression_matches"]) == 1.0
    assert got.set_index("gene").loc["GA", "repression_score"] == 0.5


def test_trrust_scores_give_the_degree_importance_D007():
    got = trrust_scores(read_trrust(SAMPLE), SEEDS, aggregator=Aggregator.DEGREE_IMPORTANCE)
    pd.testing.assert_frame_equal(got, DEGREE_SCORES, check_exact=True)
    r = golden(MATCHES).query("case == 'degree_importance'").set_index("ENSG")
    assert float(r.loc["TFA", "repression_matches"]) == 0.2
    assert float(r.loc["TFA", "activation_matches"]) == 0.04
    assert got.set_index("gene").loc["TFA", "repression_score"] == 0.5


def test_trrust_scores_put_an_unknown_edge_in_no_channel():
    got = trrust_scores(read_trrust(SAMPLE), SEEDS).set_index("gene")
    assert bool(got.loc["GC", "seed"])
    assert got.loc["GC", ["repressors", "activators"]].tolist() == [0, 0]
    assert got.loc["GC", SCORE_COLUMNS].isna().all()


def test_trrust_scores_count_an_unknown_edge_in_the_out_degree_D007():
    network = pd.DataFrame(
        {
            "source": ["R", "R", "R"],
            "target": ["A", "B", "C"],
            "edge_type": ["Activation", "Unknown", "Unknown"],
        }
    )
    got = trrust_scores(network, ["R"], aggregator=Aggregator.DEGREE_IMPORTANCE)
    assert got.set_index("gene").loc["A", "activation_score"] == pytest.approx(1 / 6)


def test_trrust_scores_count_a_triple_one_time():
    got = trrust_scores(read_trrust(SAMPLE), SEEDS).set_index("gene")
    assert got.loc["GB", "activators"] == 1
    network = pd.DataFrame(
        {
            "source": ["R", "R", "R", "R"],
            "target": ["A", "A", "A", "A"],
            "edge_type": ["Activation", "Activation", "Repression", "Repression"],
        }
    )
    got = trrust_scores(network, ["R"]).set_index("gene")
    assert got.loc["A", ["repressors", "activators"]].tolist() == [1, 1]


def test_trrust_scores_keep_the_gene_order_of_r():
    got = trrust_scores(read_trrust(SAMPLE), SEEDS)
    r = golden(MATCHES).query("case == 'mean'")
    assert got["gene"].tolist() == r["ENSG"].tolist()


def test_trrust_scores_accept_a_table_without_evidence():
    network = read_trrust(SAMPLE).drop(columns="evidence")
    pd.testing.assert_frame_equal(trrust_scores(network, SEEDS), MEAN_SCORES, check_exact=True)


def test_trrust_scores_accept_any_iterable_of_seed_genes_and_ignore_a_gene_outside_the_network():
    expected = trrust_scores(read_trrust(SAMPLE), SEEDS)
    got = trrust_scores(read_trrust(SAMPLE), iter([*SEEDS, "NOT_A_GENE", "TFA"]))
    pd.testing.assert_frame_equal(got, expected, check_exact=True)


def test_trrust_scores_give_a_gene_that_is_its_own_regulator_its_own_message():
    network = pd.DataFrame({"source": ["A"], "target": ["A"], "edge_type": ["Repression"]})
    got = trrust_scores(network, ["A"])
    assert got["repressors"].tolist() == [1]
    assert got["repression_score"].tolist() == [1.0]


def test_trrust_scores_do_not_warn_with_a_seed_gene_in_the_network():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        trrust_scores(read_trrust(SAMPLE), ["GD"])


def test_trrust_scores_warn_without_a_seed_gene_in_the_network_D058():
    seeds = ["ENSG00000000001", "ENSG00000000002"]
    message = (
        "no gene of seed_genes is a gene of the TRRUST table; each seed flag is 0. "
        "TRRUST names a gene by its symbol, for example 'TFA'"
    )
    with pytest.warns(SeedMatchWarning, match=exactly(message)) as caught:
        got = trrust_scores(read_trrust(SAMPLE), seeds)
    assert len(caught) == 1
    assert caught[0].filename == __file__
    assert not got["seed"].any()
    scores = got[SCORE_COLUMNS]
    assert ((scores == 0.0) | scores.isna()).to_numpy().all()
    r = golden(MATCHES).query("case == 'mean_no_seed'")
    assert r["ENSG"].tolist() == got["gene"].tolist()
    assert (r[["repression_matches", "activation_matches"]].astype(float) == 0).to_numpy().all()


def test_trrust_scores_name_the_first_gene_of_the_network_in_the_warning_D058():
    network = pd.DataFrame({"source": ["REG"], "target": ["TGT"], "edge_type": ["Activation"]})
    with pytest.warns(SeedMatchWarning, match=r"for example 'REG'$"):
        trrust_scores(network, [])


@pytest.mark.parametrize(
    ("network", "message"),
    [
        ([1, 2], "network has type list; expected a pandas.DataFrame"),
        (
            pd.DataFrame({"source": ["a"], "target": ["b"]}),
            (
                "TRRUST table has missing columns: ['edge_type']; expected "
                "['source', 'target', 'edge_type']"
            ),
        ),
        (
            pd.DataFrame({"source": ["a", "b"], "target": ["c", 5], "edge_type": ["Unknown"] * 2}),
            "TRRUST table column 'target' holds 5 in row 1; expected a str in each row",
        ),
        (
            pd.DataFrame(
                {"source": [NAN, "b"], "target": ["c", "d"], "edge_type": ["Unknown"] * 2}
            ),
            "TRRUST table column 'source' holds nan in row 0; expected a str in each row",
        ),
        (
            pd.DataFrame({"source": ["a"], "target": ["c"], "edge_type": ["activation"]}),
            f"TRRUST table column 'edge_type' holds 'activation' in row 0; expected {MODES}",
        ),
        (
            pd.DataFrame(
                {"source": ["a", "b"], "target": ["c", "d"], "edge_type": ["Unknown", "Other"]}
            ),
            f"TRRUST table column 'edge_type' holds 'Other' in row 1; expected {MODES}",
        ),
        (
            pd.DataFrame({"source": [], "target": [], "edge_type": []}),
            "the TRRUST table has no row; expected at least one edge",
        ),
    ],
    ids=[
        "list",
        "missing column",
        "number",
        "NaN",
        "bad mode",
        "bad mode in row 1",
        "no row",
    ],
)
def test_trrust_scores_reject_a_bad_network(network, message):
    with pytest.raises(InputError, match=exactly(message)):
        trrust_scores(network, SEEDS)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        (
            {"seed_genes": "TFA"},
            "seed_genes 'TFA' is a str; expected a sequence of str, for example ('TFA',)",
        ),
        ({"seed_genes": [1]}, "seed_genes[0] 1 has type int; expected a str"),
        (
            {"aggregator": "mean"},
            (
                "aggregator 'mean' is not valid; expected one of Aggregator.MEAN, "
                "Aggregator.DEGREE_IMPORTANCE"
            ),
        ),
    ],
    ids=["str seed", "number seed", "aggregator text"],
)
def test_trrust_scores_reject_a_bad_argument(kwargs, message):
    arguments = {"seed_genes": SEEDS, **kwargs}
    with pytest.raises(InputError, match=exactly(message)):
        trrust_scores(read_trrust(SAMPLE), **arguments)


def test_sender_weights_give_one_for_the_mean():
    weights = sender_weights(Aggregator.MEAN, np.array([0, 2, 4], dtype=np.int64))
    assert weights.dtype == np.float64
    assert weights.tolist() == [1.0, 1.0, 1.0]


def test_sender_weights_give_the_inverse_of_the_out_degree_and_zero_without_a_target():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        weights = sender_weights(Aggregator.DEGREE_IMPORTANCE, np.array([0, 2, 4], dtype=np.int64))
    assert weights.dtype == np.float64
    assert weights.tolist() == [0.0, 0.5, 0.25]


def test_pair_keys_give_an_int64_key_to_each_pair_also_for_a_narrow_code():
    first = np.array([0, 1, 60000], dtype=np.int32)
    second = np.array([2, 0, 1], dtype=np.int32)
    keys = pair_keys(first, second, 70000)
    assert keys.dtype == np.int64
    assert keys.tolist() == [2, 70000, 60000 * 70000 + 1]


def test_mean_messages_count_the_senders_and_give_nan_without_a_sender():
    receivers = np.array([1, 1, 2], dtype=np.intp)
    senders = np.array([0, 2, 0], dtype=np.intp)
    values = np.array([1.0, 0.0, 1.0])
    weights = np.array([1.0, 1.0, 0.5])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        counts, scores = mean_messages(receivers, senders, values, weights)
    assert counts.dtype == np.int64
    assert counts.tolist() == [0, 2, 1]
    assert scores.dtype == np.float64
    assert math.isnan(scores[0])
    assert scores[1:].tolist() == [((1 * 1 + 0.5 * 1) / 2 + 0) / 2, ((1 * 1) / 1 + 1) / 2]


def test_mean_messages_give_the_size_of_the_genes_without_an_edge():
    none = np.array([], dtype=np.intp)
    counts, scores = mean_messages(none, none, np.zeros(3), np.ones(3))
    assert counts.tolist() == [0, 0, 0]
    assert np.isnan(scores).all()
    assert len(scores) == 3


def test_top_level_package_exports_the_trrust_names():
    assert seedmatch.read_trrust is read_trrust
    assert seedmatch.trrust_scores is trrust_scores
    assert seedmatch.Aggregator is Aggregator
