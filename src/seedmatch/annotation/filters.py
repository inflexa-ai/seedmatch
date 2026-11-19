"""The transcript filter and the options of the feature extraction."""

from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from enum import Enum, StrEnum, auto
from typing import Final

from seedmatch._checks import as_flag, as_limit, shown
from seedmatch.errors import InputError

_MIN_SUPPORT_LEVEL: Final = 1
_MAX_SUPPORT_LEVEL: Final = 5
_MIN_ENTREZ_ID: Final = 1
_CHROMOSOMES: Final = [*(str(number) for number in range(1, 23)), "X", "Y"]
STANDARD_CHROMOSOMES: Final = frozenset(_CHROMOSOMES) | {f"chr{name}" for name in _CHROMOSOMES}


class FeatureType(StrEnum):
    """The feature of a transcript. The values are the names of `feature.type` in R."""

    UTR3 = "3UTR"
    UTR5 = "5UTR"
    EXONS = "exons"
    CDS = "cds"


class PerGene(Enum):
    """The transcripts of a gene that give a result."""

    EACH_TRANSCRIPT = auto()
    LONGEST = auto()
    UNION = auto()


def _as_text(name: str, item: object) -> str:
    if not isinstance(item, str):
        message = f"{name} {shown(item)} has type {type(item).__name__}; expected a str"
        raise InputError(message)
    return item


def _as_support_level(value: object) -> int:
    level = as_limit("support_level", value)
    if not _MIN_SUPPORT_LEVEL <= level <= _MAX_SUPPORT_LEVEL:
        message = (
            f"support_level {level} is out of range; "
            f"expected an integer from {_MIN_SUPPORT_LEVEL} to {_MAX_SUPPORT_LEVEL}"
        )
        raise InputError(message)
    return level


def _as_entrez_id(name: str, item: object) -> int:
    entrez_id = as_limit(name, item)
    if entrez_id < _MIN_ENTREZ_ID:
        message = (
            f"{name} {entrez_id} is below {_MIN_ENTREZ_ID}; "
            f"expected an Entrez ID of {_MIN_ENTREZ_ID} or more"
        )
        raise InputError(message)
    return entrez_id


def _as_set[T](
    name: str,
    value: Iterable[object] | None,
    as_item: Callable[[str, object], T],
    item_type: str,
    example: str | None = None,
) -> tuple[T, ...] | None:
    """Return the values as a tuple of checked items. `None` means that the field is not set.

    The example for a bare `str` is the `str` in a tuple, unless the caller gives another one.
    """
    if value is None:
        return None
    if isinstance(value, str):
        message = (
            f"{name} {shown(value)} is a str; expected a sequence of {item_type}, "
            f"for example ({example or repr(value)},)"
        )
        raise InputError(message)
    try:
        items = tuple(value)
    except TypeError:
        message = f"{name} {shown(value)} has type {type(value).__name__}; expected a sequence"
        raise InputError(message) from None
    if not items:
        message = f"{name} is empty; expected at least one value, or None for no filter"
        raise InputError(message)
    return tuple(as_item(f"{name}[{position}]", item) for position, item in enumerate(items))


@dataclass(frozen=True)
class TranscriptFilter:
    """The transcripts to keep. The defaults are the defaults of `load_annotations`.

    Each field that is set adds one condition, and a transcript must pass each condition. A field
    with a set of values keeps the transcripts that have one of the values. `None` leaves the field
    unset, and an empty sequence raises `InputError`, because it keeps no transcript.
    """

    standard_chromosomes: bool = True
    require_tx_id: bool = True
    protein_coding: bool = True
    canonical: bool = True
    support_level: int | None = None
    gene_ids: Sequence[str] | None = None
    tx_ids: Sequence[str] | None = None
    symbols: Sequence[str] | None = None
    tx_biotypes: Sequence[str] | None = None
    gene_biotypes: Sequence[str] | None = None
    entrez_ids: Sequence[int] | None = None

    def __post_init__(self) -> None:
        """Make sure that the fields are valid, and store each set of values as a tuple."""
        for name in ("standard_chromosomes", "require_tx_id", "protein_coding", "canonical"):
            object.__setattr__(self, name, as_flag(name, getattr(self, name)))
        if self.support_level is not None:
            object.__setattr__(self, "support_level", _as_support_level(self.support_level))
        for name in ("gene_ids", "tx_ids", "symbols", "tx_biotypes", "gene_biotypes"):
            object.__setattr__(self, name, _as_set(name, getattr(self, name), _as_text, "str"))
        entrez_ids = _as_set("entrez_ids", self.entrez_ids, _as_entrez_id, "int", "1006")
        object.__setattr__(self, "entrez_ids", entrez_ids)
