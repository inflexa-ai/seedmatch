import bz2
import gzip
import lzma
import os
import re
from pathlib import Path

import pytest
from helpers import damaged_compressed_files

from seedmatch._io import local_path, open_text
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


@pytest.mark.parametrize(
    "compress",
    [bytes, gzip.compress, bz2.compress, lzma.compress],
    ids=["none", "gzip", "bz2", "xz"],
)
def test_open_text_reads_the_text_of_a_file_whatever_its_compression_and_suffix(tmp_path, compress):
    # The first bytes show the compression, and each line end stays as the file has it.
    path = tmp_path / "table.txt"
    path.write_bytes(compress("a\tb\r\nc\u00e9\n".encode()))
    with open_text(path) as handle:
        assert handle.read() == "a\tb\r\nc\u00e9\n"


def test_open_text_reads_the_lines_of_a_file_one_at_a_time(tmp_path):
    path = tmp_path / "table.txt.gz"
    path.write_bytes(gzip.compress(b"first\nsecond\n"))
    with open_text(path) as handle:
        assert list(handle) == ["first\n", "second\n"]


def test_open_text_of_a_missing_file_raises_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError), open_text(tmp_path / "absent.txt"):
        pytest.fail("the file does not exist")


DAMAGED = damaged_compressed_files()


@pytest.mark.parametrize("name", DAMAGED)
def test_open_text_of_a_damaged_compressed_file_raises_input_error_with_the_cause(tmp_path, name):
    data, cause = DAMAGED[name]
    path = tmp_path / "table.bin"
    path.write_bytes(data)
    prefix = f"{str(path)!r}: the compressed data is damaged or truncated ("
    with pytest.raises(InputError, match=re.escape(prefix)) as raised, open_text(path) as handle:
        handle.read()
    assert str(raised.value).endswith("); expected a complete gzip, bzip2, or xz file")
    assert type(raised.value.__cause__) is cause


def test_open_text_of_a_directory_raises_input_error(tmp_path):
    message = f"{str(tmp_path)!r}: the path is a directory; expected a file"
    with pytest.raises(InputError, match=re.escape(message)) as raised, open_text(tmp_path):
        pytest.fail("the path is a directory")
    assert isinstance(raised.value.__cause__, IsADirectoryError)


def test_open_text_gives_an_error_of_the_consumer_as_it_is_for_a_compressed_file(tmp_path):
    path = tmp_path / "table.gz"
    path.write_bytes(gzip.compress(b"text\n"))
    with pytest.raises(OSError, match="disk error") as raised, open_text(path):
        raise OSError(5, "disk error")
    assert not isinstance(raised.value, InputError)
