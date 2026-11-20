from collections import Counter
from urllib.parse import urlsplit

import pytest

from seedmatch.errors import InputError
from seedmatch.resources import REGISTRY, STRING_VERSION, Entry, Group, Resource, lookup

REFERENCE_KEYS = ("hg38", "mm39", "mm10", "rnor6", "rnor7", "rnor7.113")
DIGEST = "0" * 64


def resources():
    return [resource for entry in REGISTRY for resource in entry.resources]


def test_each_registry_url_is_https():
    assert resources()
    assert [r.url for r in resources() if urlsplit(r.url).scheme != "https"] == []


def test_the_name_is_the_last_part_of_the_url_path():
    assert Resource("https://example.org/a/b/file.tsv.gz?x=1/y", 1, DIGEST).name == "file.tsv.gz"
    twobit, ensdb = lookup(Group.REFERENCE, "rnor7").resources
    assert twobit.name == "Rattus_norvegicus.mRatBN7.2.dna_sm.toplevel.2bit"
    assert ensdb.name == "EnsDb.Rnorvegicus.v109.sqlite"


def test_each_group_and_key_occurs_one_time():
    counts = Counter((entry.group, entry.key) for entry in REGISTRY)
    assert [pair for pair, count in counts.items() if count > 1] == []


def test_the_names_in_one_entry_are_different():
    for entry in REGISTRY:
        names = [resource.name for resource in entry.resources]
        assert len(names) == len(set(names)), entry


def test_the_reference_keys_are_the_builds_of_r():
    entries = [entry for entry in REGISTRY if entry.group is Group.REFERENCE]
    assert tuple(entry.key for entry in entries) == REFERENCE_KEYS
    for entry in entries:
        twobit, ensdb = entry.resources
        assert twobit.name.endswith(".2bit"), entry
        assert ensdb.name.endswith(".sqlite"), entry


def test_the_string_urls_name_the_pinned_version_D041():
    assert STRING_VERSION == "12.5"
    for entry in REGISTRY:
        if entry.group is Group.STRING:
            (resource,) = entry.resources
            assert resource.name == f"{entry.key}.protein.links.v{STRING_VERSION}.txt.gz"
            assert f"/protein.links.v{STRING_VERSION}/" in resource.url


def test_lookup_gives_the_entry_of_the_key():
    entry = lookup(Group.TRRUST, "human")
    assert (entry.group, entry.key) == (Group.TRRUST, "human")
    assert [resource.name for resource in entry.resources] == ["trrust_rawdata.human.tsv"]


def test_lookup_reads_the_registry_that_the_caller_gives():
    own = Entry(
        Group.TRRUST, "rat", "a test entry", (Resource("https://example.org/a", 1, DIGEST),)
    )
    assert lookup(Group.TRRUST, "rat", registry=(own,)) is own
    message = "trrust key 'human' is not in the registry; expected one of 'rat'"
    with pytest.raises(InputError) as raised:
        lookup(Group.TRRUST, "human", registry=(own,))
    assert str(raised.value) == message


def test_lookup_names_the_permitted_keys():
    with pytest.raises(InputError) as raised:
        lookup(Group.TRRUST, "rat")
    assert str(raised.value) == (
        "trrust key 'rat' is not in the registry; expected one of 'human', 'mouse'"
    )


def test_lookup_names_the_permitted_keys_in_registry_order():
    with pytest.raises(InputError) as raised:
        lookup(Group.REFERENCE, "hg19")
    assert str(raised.value) == (
        "reference key 'hg19' is not in the registry; "
        "expected one of 'hg38', 'mm39', 'mm10', 'rnor6', 'rnor7', 'rnor7.113'"
    )


def test_lookup_names_a_group_that_has_no_entry():
    with pytest.raises(InputError) as raised:
        lookup(Group.TRRUST, "human", registry=())
    assert str(raised.value) == "the registry has no trrust entry"


def test_lookup_rejects_text_for_the_group():
    with pytest.raises(InputError) as raised:
        lookup("trrust", "human")
    assert str(raised.value) == (
        "group 'trrust' is not valid; expected one of "
        "Group.REFERENCE, Group.EXAMPLE, Group.MIRDB, Group.TRRUST, Group.STRING"
    )


def test_lookup_rejects_a_key_that_is_not_text():
    with pytest.raises(InputError) as raised:
        lookup(Group.STRING, 9606)
    assert str(raised.value) == "key 9606 has type int; expected a str"


def test_rnor7_and_rnor7_113_share_the_2bit_resource():
    rnor7 = lookup(Group.REFERENCE, "rnor7")
    rnor7_113 = lookup(Group.REFERENCE, "rnor7.113")
    assert rnor7.resources[0] == rnor7_113.resources[0]
    assert rnor7.resources[1] != rnor7_113.resources[1]


@pytest.mark.parametrize(
    ("url", "name"),
    [
        ("https://example.org/", ""),
        ("https://example.org/data/", ""),
        ("https://example.org/data/.", "."),
        ("https://example.org/a/..", ".."),
    ],
)
def test_a_resource_refuses_a_url_without_a_usable_file_name(url, name):
    with pytest.raises(InputError) as raised:
        Resource(url, 1, DIGEST)
    assert str(raised.value) == (
        f"{url!r} file name {name!r} is not valid; expected a name that is not empty, '.', or '..'"
    )


def test_a_resource_refuses_a_size_below_1():
    with pytest.raises(InputError) as raised:
        Resource("https://example.org/a", 0, DIGEST)
    assert str(raised.value) == (
        "'https://example.org/a' size 0 is less than 1; expected 1 or more"
    )


@pytest.mark.parametrize("size", [True, 1.0, "1"])
def test_a_resource_refuses_a_size_that_is_not_an_int(size):
    with pytest.raises(InputError) as raised:
        Resource("https://example.org/a", size, DIGEST)
    assert str(raised.value) == (
        f"'https://example.org/a' size {size!r} has type {type(size).__name__}; expected an integer"
    )


@pytest.mark.parametrize("sha256", ["", "0" * 63, "0" * 65, "A" * 64, "g" * 64, "0" * 64 + "\n"])
def test_a_resource_refuses_a_digest_that_is_not_64_lowercase_hexadecimal_digits(sha256):
    with pytest.raises(InputError) as raised:
        Resource("https://example.org/a", 1, sha256)
    assert str(raised.value) == (
        f"'https://example.org/a' sha256 {sha256!r} is not valid; "
        "expected 64 lowercase hexadecimal digits"
    )
