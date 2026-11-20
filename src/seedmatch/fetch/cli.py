"""The command `seedmatch-fetch`: download the files of the registry into a directory."""

import argparse
import sys
from collections.abc import Sequence
from dataclasses import replace
from pathlib import Path
from typing import Final, assert_never

from seedmatch.errors import FetchError, InputError
from seedmatch.fetch.download import DEFAULT_FETCH_OPTS, FetchOpts, fetch_resource
from seedmatch.fetch.http import HttpResponse
from seedmatch.resources import REGISTRY, Entry, Group, Resource, lookup

_PROG: Final = "seedmatch-fetch"
_LIST: Final = "list"


def main(
    argv: Sequence[str] | None = None,
    *,
    opts: FetchOpts = DEFAULT_FETCH_OPTS,
    registry: Sequence[Entry] = REGISTRY,
) -> int:
    """Run the command and return its exit status: 0 for success and 1 for an error.

    A usage error ends with `SystemExit(2)`, as `argparse` does it.
    """
    args = _parser().parse_args(argv)
    command: str = args.command
    if command == _LIST:
        _list(registry)
        return 0
    group = Group(command)
    key: str | None = args.key
    try:
        entry = lookup(group, _only_key(group, registry) if key is None else key, registry=registry)
        _fetch(entry, args.dest, opts)
    except (InputError, FetchError) as error:
        sys.stderr.write(f"{_PROG}: {error}\n")
        return 1
    return 0


def _key_name(group: Group) -> str | None:
    match group:
        case Group.REFERENCE:
            return "BUILD"
        case Group.EXAMPLE:
            return "NAME"
        case Group.MIRDB:
            return None
        case Group.TRRUST:
            return "ORGANISM"
        case Group.STRING:
            return "TAXON"
        case _:
            assert_never(group)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=_PROG,
        description="Download the files of the seedmatch registry. Run 'list' to see them.",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    for group in Group:
        command = commands.add_parser(group.value, help=f"download the {group.value} files")
        key_name = _key_name(group)
        if key_name is None:
            command.set_defaults(key=None)
        else:
            command.add_argument(
                "key", metavar=key_name, help=f"the {group.value} key; run 'list' to see the keys"
            )
        command.add_argument("dest", metavar="DEST", help="the directory for the files")
    commands.add_parser(_LIST, help="print each file of the registry")
    return parser


def _only_key(group: Group, registry: Sequence[Entry]) -> str:
    keys = [entry.key for entry in registry if entry.group is group]
    if len(keys) != 1:
        message = f"the registry has {len(keys)} {group.value} entries; expected 1"
        raise InputError(message)
    return keys[0]


def _list(registry: Sequence[Entry]) -> None:
    for entry in registry:
        for resource in entry.resources:
            fields = (entry.group.value, entry.key, resource.name, str(resource.size))
            sys.stdout.write("\t".join((*fields, entry.description)) + "\n")


def _fetch(entry: Entry, dest: str, opts: FetchOpts) -> None:
    for resource in entry.resources:
        path, requested = _fetch_one(resource, dest, opts)
        sys.stdout.write(f"{path}\t{'downloaded' if requested else 'present'}\n")


def _fetch_one(resource: Resource, dest: str, opts: FetchOpts) -> tuple[Path, bool]:
    """Fetch the resource, and tell whether it sent a request. A present file sends none."""
    requested = False

    def open_url(url: str) -> HttpResponse:
        nonlocal requested
        if not requested:
            requested = True
            sys.stderr.write(f"{_PROG}: downloading {resource.name} ({resource.size} bytes)\n")
        return opts.open_url(url)

    return fetch_resource(resource, dest, opts=replace(opts, open_url=open_url)), requested
