import os
import re
from pathlib import Path

import pytest

from seedmatch._io import local_path
from seedmatch.errors import InputError


class Wrapped:
    def __init__(self, text):
        self.text = text

    def __fspath__(self):
        return self.text


@pytest.mark.parametrize(
    "url", ["https://example.org/de.tsv", "ftp://example.org/de.tsv", "s3://bucket/de.tsv"]
)
def test_local_path_refuses_a_url(url):
    message = f"{url!r}: the path is a URL; expected a local path"
    with pytest.raises(InputError, match=re.escape(message)):
        local_path(url)


def test_local_path_refuses_a_path_like_object_that_holds_a_url():
    url = "https://example.org/de.tsv"
    message = f"{url!r}: the path is a URL; expected a local path"
    with pytest.raises(InputError, match=re.escape(message)):
        local_path(Wrapped(url))


@pytest.mark.parametrize("text", ["table.tsv", "data/table.tsv", "/mnt/refs/a.2bit", "dir/a://b"])
def test_local_path_gives_a_path_for_a_str(text):
    path = local_path(text)
    assert isinstance(path, Path)
    assert path == Path(text)


def test_local_path_gives_a_path_for_a_path_like_object():
    path = local_path(Wrapped("data/table.tsv"))
    assert isinstance(path, Path)
    assert path == Path("data/table.tsv")
    assert os.fspath(path) == "data/table.tsv"


def test_local_path_gives_a_path_for_a_path():
    assert local_path(Path("data/table.tsv")) == Path("data/table.tsv")
