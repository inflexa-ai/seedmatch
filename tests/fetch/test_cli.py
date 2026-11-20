import hashlib
import importlib
import runpy
import sys
import tomllib
from pathlib import Path

import pytest
from helpers import DATA
from memory_client import MemoryClient, Route

from seedmatch import FeatureType, Reference
from seedmatch.fetch import FetchOpts
from seedmatch.fetch.cli import main
from seedmatch.resources import REGISTRY, Entry, Group, Resource

PYPROJECT = Path(__file__).parents[2] / "pyproject.toml"
ANNOTATION = DATA / "annotation"
URL_PREFIX = "https://example.org/syn/"
SYNTHETIC_NAMES = ("synthetic.2bit", "synthetic.sqlite")


def synthetic_entry(group=Group.REFERENCE, key="syn"):
    resources = []
    routes = {}
    for name in SYNTHETIC_NAMES:
        body = (ANNOTATION / name).read_bytes()
        resource = Resource(URL_PREFIX + name, len(body), hashlib.sha256(body).hexdigest())
        resources.append(resource)
        routes[resource.url] = Route(body=body)
    entry = Entry(group, key, "synthetic test files", tuple(resources))
    return entry, routes


def run(argv, client, registry):
    return main(argv, opts=FetchOpts(open_url=client), registry=registry)


def reported(directory, *states):
    pairs = zip(SYNTHETIC_NAMES, states, strict=True)
    return "".join(f"{directory / name}\t{state}\n" for name, state in pairs)


def downloading(name):
    size = len((ANNOTATION / name).read_bytes())
    return f"seedmatch-fetch: downloading {name} ({size} bytes)\n"


def test_fetch_reference_then_open_it(tmp_path, capsys):
    entry, routes = synthetic_entry()
    client = MemoryClient(routes)
    dest = tmp_path / "refs" / "syn"

    code = run(["reference", "syn", str(dest)], client, [entry])

    out, err = capsys.readouterr()
    assert code == 0
    assert out == reported(dest, "downloaded", "downloaded")
    assert err == downloading("synthetic.2bit") + downloading("synthetic.sqlite")
    assert client.requests == [URL_PREFIX + name for name in SYNTHETIC_NAMES]
    with Reference.from_directory(dest) as ref:
        assert len(ref.sequences(FeatureType.UTR3)) > 0


def test_a_second_run_reports_present(tmp_path, capsys):
    entry, routes = synthetic_entry()
    client = MemoryClient(routes)
    run(["reference", "syn", str(tmp_path)], client, [entry])
    capsys.readouterr()

    code = run(["reference", "syn", str(tmp_path)], client, [entry])

    out, err = capsys.readouterr()
    assert code == 0
    assert out == reported(tmp_path, "present", "present")
    assert err == ""
    assert len(client.requests) == len(SYNTHETIC_NAMES)


def test_a_present_file_and_a_missing_file_get_their_own_state(tmp_path, capsys):
    entry, routes = synthetic_entry()
    client = MemoryClient(routes)
    (tmp_path / "synthetic.2bit").write_bytes((ANNOTATION / "synthetic.2bit").read_bytes())

    code = run(["reference", "syn", str(tmp_path)], client, [entry])

    out, err = capsys.readouterr()
    assert code == 0
    assert out == reported(tmp_path, "present", "downloaded")
    assert err == downloading("synthetic.sqlite")
    assert client.requests == [URL_PREFIX + "synthetic.sqlite"]


@pytest.mark.parametrize(
    ("argv", "group"),
    [
        (["example", "syn"], Group.EXAMPLE),
        (["trrust", "syn"], Group.TRRUST),
        (["string", "syn"], Group.STRING),
        (["mirdb"], Group.MIRDB),
    ],
)
def test_each_group_has_its_subcommand(tmp_path, capsys, argv, group):
    entry, routes = synthetic_entry(group)
    client = MemoryClient(routes)

    code = run([*argv, str(tmp_path)], client, [entry])

    out, err = capsys.readouterr()
    assert code == 0
    assert out == reported(tmp_path, "downloaded", "downloaded")
    assert err == downloading("synthetic.2bit") + downloading("synthetic.sqlite")
    assert client.requests == [URL_PREFIX + name for name in SYNTHETIC_NAMES]


@pytest.mark.parametrize("entries", [0, 2])
def test_mirdb_needs_exactly_one_entry(tmp_path, capsys, entries):
    registry = [synthetic_entry(Group.MIRDB, f"v{number}")[0] for number in range(entries)]
    client = MemoryClient({})

    code = run(["mirdb", str(tmp_path)], client, registry)

    out, err = capsys.readouterr()
    assert code == 1
    assert out == ""
    assert err == f"seedmatch-fetch: the registry has {entries} mirdb entries; expected 1\n"
    assert client.requests == []


def test_list_prints_each_resource(capsys):
    code = main(["list"])

    out, err = capsys.readouterr()
    assert code == 0
    assert err == ""
    lines = out.splitlines()
    assert len(lines) == sum(len(entry.resources) for entry in REGISTRY)
    assert all(len(line.split("\t")) == 5 for line in lines)
    trrust_human = "trrust\thuman\ttrrust_rawdata.human.tsv\t297659\t"
    assert trrust_human + "TRRUST v2 transcription factor targets, human" in lines
    sizes = {name: size for _, _, name, size, _ in (line.split("\t") for line in lines)}
    assert sizes["EnsDb.Rnorvegicus.v109.sqlite"] == "179650560"


def test_an_unknown_key_returns_1_and_names_the_keys(tmp_path, capsys):
    entry, routes = synthetic_entry()
    client = MemoryClient(routes)
    dest = tmp_path / "refs"

    code = run(["reference", "zzz", str(dest)], client, [entry])

    out, err = capsys.readouterr()
    assert code == 1
    assert out == ""
    assert err == (
        "seedmatch-fetch: reference key 'zzz' is not in the registry; expected one of 'syn'\n"
    )
    assert client.requests == []
    assert not dest.exists()


def test_a_fetch_error_returns_1_with_the_message(tmp_path, capsys):
    entry, routes = synthetic_entry()
    failing = URL_PREFIX + "synthetic.sqlite"
    client = MemoryClient({**routes, failing: Route(status=404)})

    code = run(["reference", "syn", str(tmp_path)], client, [entry])

    out, err = capsys.readouterr()
    assert code == 1
    assert out == f"{tmp_path / 'synthetic.2bit'}\tdownloaded\n"
    assert err == (
        downloading("synthetic.2bit")
        + downloading("synthetic.sqlite")
        + f"seedmatch-fetch: {failing!r} answered with status 404; expected 200\n"
    )
    assert sorted(path.name for path in tmp_path.iterdir()) == ["synthetic.2bit"]


def test_an_input_error_returns_1_with_the_message(tmp_path, capsys):
    entry, routes = synthetic_entry()
    client = MemoryClient(routes)
    dest = tmp_path / "file"
    dest.write_bytes(b"")

    code = run(["reference", "syn", str(dest)], client, [entry])

    out, err = capsys.readouterr()
    assert code == 1
    assert out == ""
    assert err == (
        f"seedmatch-fetch: dest {str(dest)!r} is not a directory; "
        "expected a directory or a path that is absent\n"
    )
    assert client.requests == []


def test_a_final_path_that_cannot_be_accessed_returns_1_with_the_message(
    tmp_path, capsys, monkeypatch
):
    entry, routes = synthetic_entry()
    client = MemoryClient(routes)
    final = tmp_path / "synthetic.2bit"
    error = PermissionError(13, "Permission denied")
    original_exists = Path.exists

    def failing_exists(self, *args, **kwargs):
        if self == final:
            raise error
        return original_exists(self, *args, **kwargs)

    monkeypatch.setattr(Path, "exists", failing_exists)

    code = run(["reference", "syn", str(tmp_path)], client, [entry])

    out, err = capsys.readouterr()
    assert code == 1
    assert out == ""
    assert err == (
        f"seedmatch-fetch: {str(final)!r} cannot be accessed: {error!r}; "
        "expected a path that can be accessed\n"
    )
    assert client.requests == []


def test_a_redirect_prints_one_downloading_line(tmp_path, capsys):
    entry, routes = synthetic_entry()
    first = entry.resources[0].url
    moved = URL_PREFIX + "moved/synthetic.2bit"
    client = MemoryClient(
        {**routes, first: Route(status=302, headers={"Location": moved}), moved: routes[first]}
    )

    code = run(["reference", "syn", str(tmp_path)], client, [entry])

    out, err = capsys.readouterr()
    assert code == 0
    assert out == reported(tmp_path, "downloaded", "downloaded")
    assert err == downloading("synthetic.2bit") + downloading("synthetic.sqlite")
    assert client.requests == [first, moved, URL_PREFIX + "synthetic.sqlite"]


GROUP_CHOICES = "{reference,example,mirdb,trrust,string,list}"


@pytest.mark.parametrize(
    ("argv", "usage", "message"),
    [
        (
            [],
            f"usage: seedmatch-fetch [-h] {GROUP_CHOICES} ...",
            "seedmatch-fetch: error: the following arguments are required: command",
        ),
        (
            ["reference", "rnor7"],
            "usage: seedmatch-fetch reference [-h] BUILD DEST",
            "seedmatch-fetch reference: error: the following arguments are required: DEST",
        ),
        (
            ["mirdb"],
            "usage: seedmatch-fetch mirdb [-h] DEST",
            "seedmatch-fetch mirdb: error: the following arguments are required: DEST",
        ),
    ],
)
def test_a_usage_error_exits_with_2(capsys, monkeypatch, argv, usage, message):
    monkeypatch.setenv("COLUMNS", "200")

    with pytest.raises(SystemExit) as caught:
        main(argv)

    out, err = capsys.readouterr()
    assert caught.value.code == 2
    assert out == ""
    assert err == f"{usage}\n{message}\n"


@pytest.mark.parametrize("group", [Group.REFERENCE, Group.EXAMPLE, Group.TRRUST, Group.STRING])
def test_the_help_of_a_key_points_to_the_list_command(capsys, group):
    with pytest.raises(SystemExit) as caught:
        main([group.value, "--help"])

    out, _ = capsys.readouterr()
    assert caught.value.code == 0
    assert f"the {group.value} key; run 'list' to see the keys" in " ".join(out.split())


def test_the_script_is_declared():
    declared = tomllib.loads(PYPROJECT.read_text(encoding="utf-8"))["project"]["scripts"]
    assert declared["seedmatch-fetch"] == "seedmatch.fetch.cli:main"
    module, _, name = declared["seedmatch-fetch"].partition(":")
    assert getattr(importlib.import_module(module), name) is main


def test_python_m_runs_the_command(capsys, monkeypatch):
    main(["list"])
    expected = capsys.readouterr().out
    monkeypatch.setattr(sys, "argv", ["seedmatch-fetch", "list"])

    with pytest.raises(SystemExit) as caught:
        runpy.run_module("seedmatch.fetch", run_name="__main__")

    out, err = capsys.readouterr()
    assert caught.value.code == 0
    assert err == ""
    assert out == expected
    assert len(out.splitlines()) == sum(len(entry.resources) for entry in REGISTRY)
