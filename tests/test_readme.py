import re
import shutil
from pathlib import Path

import pytest
from helpers import DATA, golden, sirna_de

README = Path(__file__).parents[1] / "README.md"
PYTHON_BLOCK = re.compile(r"^```python\n(.*?)^```$", re.MULTILINE | re.DOTALL)


def python_examples() -> list[str]:
    return PYTHON_BLOCK.findall(README.read_text(encoding="utf-8"))


@pytest.fixture
def readme_namespace(tmp_path, monkeypatch):
    shutil.copy(DATA / "Test_seqs.fa", tmp_path / "transcripts.fa")
    # The test table has a ninth column, mer7m8. read_deseq2_table accepts only 8 columns (D-021).
    sirna_de().drop(columns="mer7m8").to_csv(
        tmp_path / "deseq2_results.tsv", sep="\t", index=False, na_rep="NA"
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
