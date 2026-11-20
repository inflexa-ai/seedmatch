"""The registry of the upstream files that seedmatch downloads.

Each `Resource` pins the URL, the size, and the SHA-256 digest of one file. An `Entry` groups the
files of one reference build, one example data set, or one network source.
"""

import re
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from typing import Final
from urllib.parse import urlsplit

from seedmatch._checks import as_limit, as_member, as_text, check_at_least
from seedmatch.errors import InputError

STRING_VERSION: Final = "12.5"

_ANNOTATION_HUB: Final = "https://mghp.osn.xsede.org/bir190004-bucket01/AnnotationHub/"
_STRING_DOWNLOADS: Final = "https://stringdb-downloads.org/download/"
_GEO_SIRNA: Final = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184929/suppl/"
_TRRUST: Final = "https://www.grnpedia.org/trrust/data/"
_SHA256: Final = re.compile(r"[0-9a-f]{64}")
_NOT_A_FILE_NAME: Final = frozenset({"", ".", ".."})


class Group(StrEnum):
    """The kind of an entry. The value is the name of the group in the `seedmatch-fetch` command."""

    REFERENCE = "reference"
    EXAMPLE = "example"
    MIRDB = "mirdb"
    TRRUST = "trrust"
    STRING = "string"


@dataclass(frozen=True)
class Resource:
    """One file of an upstream source, with the size and the SHA-256 digest that it must have."""

    url: str
    size: int
    sha256: str

    def __post_init__(self) -> None:
        """Reject a file name, a size, or a digest that a download cannot use."""
        if self.name in _NOT_A_FILE_NAME:
            message = (
                f"{self.url!r} file name {self.name!r} is not valid; "
                "expected a name that is not empty, '.', or '..'"
            )
            raise InputError(message)
        size_name = f"{self.url!r} size"
        check_at_least(size_name, as_limit(size_name, self.size), 1)
        if not _SHA256.fullmatch(self.sha256):
            message = (
                f"{self.url!r} sha256 {self.sha256!r} is not valid; "
                "expected 64 lowercase hexadecimal digits"
            )
            raise InputError(message)

    @property
    def name(self) -> str:
        """The last part of the URL path. A download stores the file under this name."""
        return urlsplit(self.url).path.rpartition("/")[2]


@dataclass(frozen=True)
class Entry:
    """The files of one `key` of one `group`. `description` names the source of the files."""

    group: Group
    key: str
    description: str
    resources: tuple[Resource, ...]


def _hub(path: str, size: int, sha256: str) -> Resource:
    return Resource(_ANNOTATION_HUB + path, size, sha256)


def _string(taxon: str, size: int, sha256: str) -> Resource:
    name = f"{taxon}.protein.links.v{STRING_VERSION}.txt.gz"
    return Resource(f"{_STRING_DOWNLOADS}protein.links.v{STRING_VERSION}/{name}", size, sha256)


_RNOR7_2BIT: Final = _hub(
    "ensembl/release-107/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna_sm.toplevel.2bit",
    661989722,
    "73eacccd66a8843cced837d0cb0a272b42b7214b4a7c62a3c48963509ea914f3",
)

REGISTRY: Final[tuple[Entry, ...]] = (
    Entry(
        Group.REFERENCE,
        "hg38",
        "AH106283 (2bit) and AH109606 (EnsDb), Ensembl 109",
        (
            _hub(
                "ensembl/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.2bit",
                774952244,
                "6464ece374014a850de7da50bfaab5083b73197bdeb1393731065b6ab7560be6",
            ),
            _hub(
                "AHEnsDbs/v109/EnsDb.Hsapiens.v109.sqlite",
                445898752,
                "a3d9a352b0bc07ac6f27e89df31308ff8ead36380608644cc5939de3212b8253",
            ),
        ),
    ),
    Entry(
        Group.REFERENCE,
        "mm39",
        "AH106440 (2bit) and AH109655 (EnsDb), Ensembl 109",
        (
            _hub(
                "ensembl/release-107/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.2bit",
                682060023,
                "eac802b955d7885f0eb620ba5d49f119397bf5c630593e20d16f0f5bf3d4e8e9",
            ),
            _hub(
                "AHEnsDbs/v109/EnsDb.Mmusculus.v109.sqlite",
                269914112,
                "6d2efcbfd5d3cb14ba2037d67583e1d4027b1eba81d489a9a4938f948107a7c2",
            ),
        ),
    ),
    Entry(
        Group.REFERENCE,
        "mm10",
        "AH88475 (2bit) and AH89211 (EnsDb), Ensembl 102",
        (
            _hub(
                "ensembl/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.2bit",
                682724912,
                "b12e356f498d74a8af43e374d3240ab058d1d6b9f76baefb277d5edb469f488a",
            ),
            _hub(
                "AHEnsDbs/v102/EnsDb.Mmusculus.v102.sqlite",
                265134080,
                "b8b915e96131fadf494d3f7a79e9e7b31a54aa002cc60e4ba2278797a2ce0c08",
            ),
        ),
    ),
    Entry(
        Group.REFERENCE,
        "rnor6",
        "AH93578 (2bit) and AH95846 (EnsDb), Ensembl 104",
        (
            _hub(
                "ensembl/release-104/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna_sm.toplevel.2bit",
                718742859,
                "e5eb2cd9a17d8f52c33d07bb423085addbc74a3c5ec3b1b2a8f2c6eb0ef6a4ee",
            ),
            _hub(
                "AHEnsDbs/v104/EnsDb.Rnorvegicus.v104.sqlite",
                117067776,
                "174312139350c7cb86af5fd82c9758c98c4fe43d6f75ca21797f2c154fa1bdae",
            ),
        ),
    ),
    Entry(
        Group.REFERENCE,
        "rnor7",
        "AH106786 (2bit) and AH109732 (EnsDb), Ensembl 109",
        (
            _RNOR7_2BIT,
            _hub(
                "AHEnsDbs/v109/EnsDb.Rnorvegicus.v109.sqlite",
                179650560,
                "2453c59e1054157798df5a4114650911d046b5f7795af91a6a2171614ba0289b",
            ),
        ),
    ),
    Entry(
        Group.REFERENCE,
        "rnor7.113",
        "AH106786 (2bit) and AH119437 (EnsDb), Ensembl 113",
        (
            _RNOR7_2BIT,
            _hub(
                "AHEnsDbs/v113/EnsDb.Rnorvegicus.v113.sqlite",
                182464512,
                "ed63640b3d6d4cd7d26f1410a0746c341a17f7d8028df03db8dc772e77faff89",
            ),
        ),
    ),
    Entry(
        Group.EXAMPLE,
        "sirna",
        "Schlegel 2022 siRNA DE tables of Ttr, GEO series GSE184929",
        (
            Resource(
                f"{_GEO_SIRNA}GSE184929_Fig3-TTR-2_groupD1_30mkg.vs.Mock.DESeq2.txt.gz",
                1124255,
                "c418c8aebd2138b892d91a0ca0a99952b083dd71f1043533d2c35045b2551142",
            ),
            Resource(
                f"{_GEO_SIRNA}GSE184929_Fig3-TTR-2_groupD4_30mkg.vs.Mock.DESeq2.txt.gz",
                1085045,
                "a2e39fb9215aa647ffd8b6e7e4e491fffccf438e8f706e826e24c88313624696",
            ),
            Resource(
                f"{_GEO_SIRNA}GSE184929_Fig3-TTR-1_groupD1_10mkg.vs.Mock.DESeq2.txt.gz",
                1094976,
                "97598daf5670d94875149cca0a3afebb7309758e41c99c35bb10274498186248",
            ),
        ),
    ),
    Entry(
        Group.MIRDB,
        "6.0",
        "miRDB v6.0 predicted miRNA targets of all species",
        (
            Resource(
                "https://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz",
                59483982,
                "446636d2182212f44ff91623f1429f7922f9575b82226962a6213d592a7e109d",
            ),
        ),
    ),
    Entry(
        Group.TRRUST,
        "human",
        "TRRUST v2 transcription factor targets, human",
        (
            Resource(
                f"{_TRRUST}trrust_rawdata.human.tsv",
                297659,
                "9b909319ccc8e36588b5a1bd3640e0dfba03936c6fc154cc4304d742bec11bed",
            ),
        ),
    ),
    Entry(
        Group.TRRUST,
        "mouse",
        "TRRUST v2 transcription factor targets, mouse",
        (
            Resource(
                f"{_TRRUST}trrust_rawdata.mouse.tsv",
                220353,
                "4847c960c7c2ce0173b056184137b1cc71e8d61d1b1edc825247504ad882dd55",
            ),
        ),
    ),
    Entry(
        Group.STRING,
        "9606",
        f"STRING v{STRING_VERSION} protein links, human",
        (
            _string(
                "9606",
                62813591,
                "3d6894ef3a0a8578b70421b7d6bac159c41ed2f6c1b5da21ab7fabb4ac21cd25",
            ),
        ),
    ),
    Entry(
        Group.STRING,
        "10090",
        f"STRING v{STRING_VERSION} protein links, mouse",
        (
            _string(
                "10090",
                75044046,
                "b3318ac1cbd0fe0e81a3d6b0b432352f72733b7cf36d2437423a5cbf2a98e9d0",
            ),
        ),
    ),
    Entry(
        Group.STRING,
        "10116",
        f"STRING v{STRING_VERSION} protein links, rat",
        (
            _string(
                "10116",
                75736490,
                "4d53bee2bc10daa046aae67532152c7389ba2374dacb18cf2b2cd3886b253387",
            ),
        ),
    ),
)


def lookup(group: Group, key: str, *, registry: Sequence[Entry] = REGISTRY) -> Entry:
    """Return the entry of the group and the key.

    Raise `InputError` for a group that is not a `Group` member, a key that is not text, and a key
    that the group does not have. The last message gives the keys of the group in registry order.
    """
    group = as_member("group", group, Group)
    key = as_text("key", key)
    keys: list[str] = []
    for entry in registry:
        if entry.group is group:
            if entry.key == key:
                return entry
            keys.append(entry.key)
    if not keys:
        message = f"the registry has no {group.value} entry"
        raise InputError(message)
    expected = ", ".join(repr(known) for known in keys)
    message = f"{group.value} key {key!r} is not in the registry; expected one of {expected}"
    raise InputError(message)
