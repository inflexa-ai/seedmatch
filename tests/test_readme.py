import gzip
import math
import re
import shutil
from pathlib import Path

import pytest
from helpers import DATA, STATISTIC_TOLERANCE, ecdf_golden_row, golden, sirna_de

from seedmatch import Alternative, Benchmark, EcdfTest

README = Path(__file__).parents[1] / "README.md"
PYTHON_BLOCK = re.compile(r"^```python\n(.*?)^```$", re.MULTILINE | re.DOTALL)


def python_examples() -> list[str]:
    return PYTHON_BLOCK.findall(README.read_text(encoding="utf-8"))


@pytest.fixture
def readme_namespace(tmp_path, monkeypatch):
    shutil.copy(DATA / "Test_seqs.fa", tmp_path / "transcripts.fa")
    shutil.copy(DATA / "annotation" / "synthetic.sqlite", tmp_path / "ensdb.sqlite")
    shutil.copy(DATA / "annotation" / "synthetic.2bit", tmp_path / "genome.2bit")
    # The test table has a ninth column, mer7m8. read_deseq2_table accepts only 8 columns (D-021).
    sirna_de().drop(columns="mer7m8").to_csv(
        tmp_path / "deseq2_results.tsv", sep="\t", index=False, na_rep="NA"
    )
    (tmp_path / "miRDB_v6.0_prediction_result.txt.gz").write_bytes(
        gzip.compress((DATA / "mirdb_sample.txt").read_bytes())
    )
    monkeypatch.chdir(tmp_path)
    namespace: dict[str, object] = {}
    for number, code in enumerate(python_examples(), start=1):
        # The code is the README of this repository, not outside input.
        exec(compile(code, f"README.md, Python example {number}", "exec"), namespace)  # noqa: S102
    return namespace


def test_readme_has_a_python_example():
    assert python_examples()


def test_readme_example_gives_the_r_report(readme_namespace):
    summary = readme_namespace["summary"].set_index("category")
    r_report = golden("seedmatchr/report_default.tsv").set_index("Category")
    row = "In silico predictions"
    assert summary.loc[row, "Total"] == float(r_report.loc[row, "Total"])


def test_readme_example_adds_the_seed_count_column(readme_namespace):
    de = readme_namespace["de"]
    assert de["mer7m8"].dtype == "int64"
    assert de["mer7m8"].sum() > 0


def test_readme_annotation_example_gives_the_r_sequences(readme_namespace):
    utr3 = readme_namespace["utr3"]
    expected = golden("annotation/utr3_sequences.tsv")
    assert list(utr3.ids) == expected["name_tx"].tolist()
    assert [utr3.sequence(i) for i in range(len(utr3))] == expected["sequence"].tolist()
    assert len(readme_namespace["utr3_counts"]) == len(utr3)


def test_readme_ecdf_example_gives_the_r_result(readme_namespace):
    result = readme_namespace["ecdf_result"]
    row = ecdf_golden_row("readme", "KS", "greater")
    assert (result.test, result.alternative) == (EcdfTest.KS, Alternative.GREATER)
    assert (result.n_target, result.n_background) == (int(row.n_target), int(row.n_background))
    assert math.isclose(result.statistic, float(row.statistic), **STATISTIC_TOLERANCE)
    assert result.p_value is not None
    assert math.isclose(result.p_value, float(row.p_value), rel_tol=1e-9)


def test_readme_kmer_example_gives_the_r_target_count(readme_namespace):
    target = readme_namespace["target"]
    rows = golden("analyses/kmer_gene_counts.tsv")
    (row,) = rows[
        (rows["case"] == "test_data_k7_s1") & (rows["kmer"] == target.sequence)
    ].itertuples()
    assert readme_namespace["genes_with_kmer"][target.sequence] == int(row.n)


def test_readme_distance_example_runs(readme_namespace):
    distances = readme_namespace["distances"]
    assert list(distances.columns) == ["seq_id", "gene_id", "distance"]
    assert len(distances) >= 1


def test_readme_alignment_example_gives_the_r_text(readme_namespace):
    rows = golden("analyses/alignment.tsv").query("case == 'mer7m8_1'")
    assert readme_namespace["alignment"] == "\n".join(rows.sort_values("line")["text"])


def test_readme_e_value_example_runs(readme_namespace):
    scored = readme_namespace["scored"]
    assert len(scored) == len(readme_namespace["de"])
    assert scored["e_value"].between(0, 1).all()


def test_readme_benchmark_example_runs(readme_namespace):
    benchmark = readme_namespace["benchmark"]
    de = readme_namespace["de"]
    gold_standard = (de["padj"] <= 0.15) & (de["log2FoldChange"] < 0)
    assert isinstance(benchmark, Benchmark)
    assert benchmark.pr.n_positive == gold_standard.sum()
    assert benchmark.pr.n_negative == (~gold_standard).sum()
    assert 0 < benchmark.pr.auc_integral < 1
    assert 0 < benchmark.roc.auc < 1


def test_readme_mirdb_example_gives_the_rat_rows(readme_namespace):
    rat_targets = readme_namespace["rat_targets"]
    table = golden("analyses/mirdb.tsv")
    rows = table[table["species"] == "rno"]
    assert len(rows) > 0
    assert rat_targets["mirdb_id"].tolist() == rows["miRDB.ID"].tolist()
    assert rat_targets["target_refseq_id"].tolist() == rows["target.REFSEQ.ID"].tolist()
    assert rat_targets["target_score"].tolist() == pytest.approx(
        rows["target.score"].astype(float).tolist(), rel=1e-9
    )
    assert set(rat_targets["species"]) == {"rno"}
