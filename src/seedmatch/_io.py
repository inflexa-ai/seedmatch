"""Paths that the readers open."""

import os
import re
from pathlib import Path
from typing import Final

from seedmatch.errors import InputError

_URL: Final = re.compile(r"^[A-Za-z][A-Za-z0-9+.-]*://")


def local_path(path: str | os.PathLike[str]) -> Path:
    """Return the path as a `Path`.

    A URL raises `InputError`, because the core modules do not use the network.
    """
    shown = os.fspath(path)
    if _URL.match(shown):
        message = f"{shown!r}: the path is a URL; expected a local path"
        raise InputError(message)
    return Path(shown)
