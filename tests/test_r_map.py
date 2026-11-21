import ast
import importlib
import pkgutil
import re
from pathlib import Path

import pytest

import seedmatch
from seedmatch.errors import InputError
from seedmatch.resources import Group, lookup

ROOT = Path(__file__).parents[1]
R_MAP = ROOT / "docs" / "r-map.md"
README = ROOT / "README.md"
NAMESPACE = ROOT / "tests" / "data" / "SeedMatchR.NAMESPACE"
TICK = re.compile(r"`([^`]+)`")
R_NAME = re.compile(r"(?:[A-Za-z][A-Za-z0-9.]*::)?([A-Za-z.][A-Za-z0-9._]*)")
PYTHON_BLOCK = re.compile(r"^```python\n(.*?)^```$", re.MULTILINE | re.DOTALL)
DOTTED = re.compile(r"seedmatch(?:\.\w+)+")


def r_map_rows() -> list[list[str]]:
    lines = [
        line for line in R_MAP.read_text(encoding="utf-8").splitlines() if line.startswith("| ")
    ]
    rows = [[cell.strip() for cell in line.strip().strip("|").split(" | ")] for line in lines[1:]]
    # The second line is the separator row of the table.
    rows = [row for row in rows if row[0] != "---"]
    assert all(len(row) == 3 for row in rows)
    return rows


def r_names(cell: str) -> set[str]:
    return {found.group(1) for tick in TICK.findall(cell) if (found := R_NAME.match(tick))}


def resolve(dotted: str) -> object:
    parts = dotted.split(".")
    for cut in range(len(parts), 0, -1):
        try:
            value: object = importlib.import_module(".".join(parts[:cut]))
        except ModuleNotFoundError:
            continue
        for part in parts[cut:]:
            value = getattr(value, part)
        return value
    raise ModuleNotFoundError(dotted)


def fetch_command_exists(command: str) -> bool:
    words = command.split()[1:]
    if words == ["list"]:
        return True
    if not words or words[0] not in set(Group):
        return False
    # An uppercase word is a placeholder, for example BUILD or DEST.
    keys = [word for word in words[1:] if not word.isupper()]
    try:
        for key in keys[:1]:
            lookup(Group(words[0]), key)
    except InputError:
        return False
    return True


def code_spans(path: Path) -> list[str]:
    # A fenced block breaks the pairing of the backticks over the whole text. No code span of
    # these documents crosses a line, thus each line pairs its own backticks.
    return [
        tick
        for line in path.read_text(encoding="utf-8").splitlines()
        for tick in TICK.findall(line)
    ]


def public_modules() -> list[str]:
    packages = [
        info.name
        for info in pkgutil.iter_modules(seedmatch.__path__)
        if info.ispkg and not info.name.startswith("_")
    ]
    return ["seedmatch", *(f"seedmatch.{name}" for name in packages)]


def readme_imports() -> set[tuple[str, str]]:
    blocks = PYTHON_BLOCK.findall(README.read_text(encoding="utf-8"))
    return {
        (node.module or "", alias.name)
        for block in blocks
        for node in ast.walk(ast.parse(block))
        if isinstance(node, ast.ImportFrom)
        for alias in node.names
    }


def test_each_exported_r_function_has_a_row_of_the_r_map():
    exported = re.findall(r"^export\((.+)\)$", NAMESPACE.read_text(encoding="utf-8"), re.MULTILINE)
    named = set().union(*(r_names(row[0]) for row in r_map_rows()))
    assert exported
    assert sorted(set(exported) - named) == []


def test_each_python_name_of_the_r_map_exists_D014():
    missing = []
    for row in r_map_rows():
        for tick in (tick for cell in row[1:] for tick in TICK.findall(cell)):
            name = tick.split("(")[0]
            if name.startswith("seedmatch-fetch"):
                if not fetch_command_exists(name):
                    missing.append(tick)
            elif DOTTED.fullmatch(name):
                try:
                    resolve(name)
                except (AttributeError, ModuleNotFoundError):
                    missing.append(tick)
    assert missing == []


@pytest.mark.parametrize("module_name", public_modules())
def test_each_public_name_is_in_the_r_map_or_the_readme(module_name):
    module = importlib.import_module(module_name)
    ticks = code_spans(R_MAP) + code_spans(README)
    imports = readme_imports()
    missing = []
    for name in module.__all__:
        top = name in seedmatch.__all__ and getattr(seedmatch, name) is getattr(module, name)
        # A name of the top-level package can stand alone. Another name needs its module path.
        prefix = (
            r"(?:^|[^\w.])(?:seedmatch(?:\.[a-z_]+)*\.)?" if top else re.escape(module_name) + r"\."
        )
        pattern = re.compile(prefix + re.escape(name) + r"(?!\w)")
        modules = {module_name, "seedmatch"} if top else {module_name}
        if not any(pattern.search(tick) for tick in ticks) and not any(
            (path, name) in imports for path in modules
        ):
            missing.append(name)
    assert missing == []
