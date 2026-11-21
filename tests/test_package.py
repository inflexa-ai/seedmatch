import re
import subprocess
import sys
from pathlib import Path

import pytest

import seedmatch
from seedmatch.errors import InputError, SeedMatchError

DIVERGENCES = Path(__file__).parents[1] / "docs" / "divergences.md"
TESTS = Path(__file__).parent


def ids_without_test(ids: list[str], test_text: str) -> list[str]:
    # A Python name cannot hold a hyphen, so a test can name D-001 as D001.
    return [i for i in ids if i not in test_text and i.replace("-", "") not in test_text]


def test_import_loads_no_fetch_and_no_plot_modules():
    code = (
        "import sys, seedmatch; print(sorted(m for m in sys.modules "
        "if m.startswith(('seedmatch.fetch', 'seedmatch.plot', 'matplotlib', 'pysam'))))"
    )
    # The arguments are the interpreter path and a constant script, so no input reaches the command.
    out = subprocess.run(  # noqa: S603
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    assert out == "[]"


def test_each_divergence_id_is_unique_and_has_a_test():
    ids = re.findall(r"^## (D-\d{3}):", DIVERGENCES.read_text(), flags=re.MULTILINE)
    assert len(ids) == len(set(ids))
    # This file names a probe ID in its own matcher test, so it must not count as a divergence test.
    test_text = "\n".join(p.read_text() for p in TESTS.rglob("test_*.py") if p != Path(__file__))
    assert ids_without_test(ids, test_text) == []


def test_divergence_id_counts_with_or_without_hyphen_in_test_text():
    assert ids_without_test(["D-999"], "def test_x_D999(): ...") == []
    assert ids_without_test(["D-999"], '"""D-999 keeps the fix."""') == []
    assert ids_without_test(["D-999"], "def test_x(): ...") == ["D-999"]


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
