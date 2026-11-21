import ast
import re
import subprocess
import sys
from pathlib import Path

import pytest

import seedmatch
from seedmatch.errors import InputError, SeedMatchError

DIVERGENCES = Path(__file__).parents[1] / "docs" / "divergences.md"
DESIGN = Path(__file__).parents[1] / "docs" / "design.md"
DIVERGENCE_FIELDS = ("ID", "R location", "R behavior", "New behavior", "Reason", "Test")
TESTS = Path(__file__).parent


def labels_of_tests(source: str) -> list[str]:
    """Return the name and the docstring of each test function in the source of a test file."""
    functions = (
        node
        for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef)
    )
    return [
        label
        for node in functions
        if node.name.startswith("test_")
        for label in (node.name, ast.get_docstring(node) or "")
    ]


def names_id(label: str, entry_id: str) -> bool:
    # A Python name cannot hold a hyphen, so a test can name D-001 as D001 or d001.
    number = entry_id.removeprefix("D-")
    return re.search(rf"(?<![A-Za-z0-9])D-?{number}(?!\d)", label, flags=re.IGNORECASE) is not None


def ids_without_test(ids: list[str], labels: list[str]) -> list[str]:
    return [i for i in ids if not any(names_id(label, i) for label in labels)]


def bad_entry_headings(text: str) -> list[str]:
    """Return each heading after `## Entries` that does not have the form `## D-NNN: <title>`."""
    entries = text.split("\n## Entries\n", 1)[1]
    headings = re.findall(r"^#.*$", entries, flags=re.MULTILINE)
    return [line for line in headings if not re.fullmatch(r"## D-\d{3}: \S.*", line)]


def test_import_loads_no_fetch_and_no_plot_modules():
    code = (
        "import sys, seedmatch; print(sorted(m for m in sys.modules "
        "if m.startswith(('seedmatch.fetch', 'seedmatch.plot', 'matplotlib', 'pysam'))))"
    )
    # The arguments are the interpreter path and a constant script, so no input reaches the command.
    result = subprocess.run(  # noqa: S603
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "[]"


def test_each_heading_of_the_entries_has_the_form_of_an_entry():
    assert bad_entry_headings(DIVERGENCES.read_text(encoding="utf-8")) == []


def test_each_divergence_id_is_unique_and_has_a_test():
    text = DIVERGENCES.read_text(encoding="utf-8")
    ids = re.findall(r"^## (D-\d{3}):", text, flags=re.MULTILINE)
    assert len(ids) == len(set(ids))
    labels = [
        label
        for path in TESTS.rglob("test_*.py")
        for label in labels_of_tests(path.read_text(encoding="utf-8"))
    ]
    assert ids_without_test(ids, labels) == []


def test_each_divergence_row_of_the_design_has_an_entry_with_each_field():
    design = DESIGN.read_text(encoding="utf-8")
    section = design.split("\n## Divergences (defect fixes)\n", 1)[1].split("\n## ", 1)[0]
    # The table has no ID column: data row n has the ID D-nnn of n. [1:] drops the header row, and
    # the separator line "|---|" does not start with "| ".
    rows = [line for line in section.splitlines() if line.startswith("| ")][1:]
    text = DIVERGENCES.read_text(encoding="utf-8")
    entries = re.split(r"^## (?=D-\d{3}:)", text, flags=re.MULTILINE)[1:]
    ids = [entry.split(":", 1)[0] for entry in entries]
    assert ids == [f"D-{n:03d}" for n in range(1, len(rows) + 1)]
    for entry_id, entry in zip(ids, entries, strict=True):
        fields = re.findall(r"^- ([A-Za-z ]+):", entry, flags=re.MULTILINE)
        assert list(dict.fromkeys(fields)) == list(DIVERGENCE_FIELDS), entry_id
        assert f"\n- ID: `{entry_id}`\n" in entry, entry_id


def test_a_divergence_id_counts_only_in_the_name_or_the_docstring_of_a_test():
    source = """
def test_x_D999(): ...
def test_y():
    '''D-998 keeps the fix.'''
def test_d997_lowercase(): ...
def test_z():
    assert "D-996"  # D-995
def helper_D994(): ...
def test_D9930(): ...
"""
    ids = ["D-999", "D-998", "D-997", "D-996", "D-995", "D-994", "D-993"]
    assert ids_without_test(ids, labels_of_tests(source)) == ["D-996", "D-995", "D-994", "D-993"]


def test_bad_entry_headings_finds_each_heading_without_the_form_of_an_entry():
    text = (
        "# Divergences\n\n## Entries\n\n## D-001: Good\n\n## D-02: Short\n\n"
        "## D-003 No colon\n\n# D-004: Level 1\n\n##D-005: No space\n\n### D-006: Level 3\n"
    )
    assert bad_entry_headings(text) == [
        "## D-02: Short",
        "## D-003 No colon",
        "# D-004: Level 1",
        "##D-005: No space",
        "### D-006: Level 3",
    ]


@pytest.mark.parametrize("base", [ValueError, SeedMatchError])
def test_input_error_is_caught_as_value_error_and_as_package_error(base):
    error = InputError("bad value 'X'")
    with pytest.raises(base, match="bad value 'X'"):
        raise error


# The README, docs/r-map.md, and the docstrings of the options name these, so a user imports
# them from seedmatch.
@pytest.mark.parametrize("name", ["InputError", "SeedMatchError", "DEFAULT_SEEDS", "NONSEED_SEEDS"])
def test_a_name_that_the_documentation_uses_is_in_the_top_level_package(name):
    assert name in seedmatch.__all__
    assert hasattr(seedmatch, name)
