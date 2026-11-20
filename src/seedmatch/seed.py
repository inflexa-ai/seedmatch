"""The guide strand, its seeds, and the target sequences that a search looks for."""

import re
from dataclasses import dataclass
from enum import StrEnum
from typing import Final, assert_never

from seedmatch._checks import as_flag, as_limit, shown
from seedmatch.errors import InputError
from seedmatch.sequence import reverse_complement

_NOT_RNA_LETTER: Final = re.compile(r"[^ACGUacgu]")
_MIN_GUIDE_LENGTH: Final = 8
_CANONICAL_START: Final = 2


class Seed(StrEnum):
    """A canonical seed. The value is the seed name in SeedMatchR."""

    MER7M8 = "mer7m8"
    MER8 = "mer8"
    MER6 = "mer6"
    MER7A1 = "mer7A1"
    FULL = "Full"
    MER18 = "18mer"
    MER15 = "15mer"


@dataclass(frozen=True)
class CustomSeed:
    """A window of the guide from `start` to `stop`, 1-based and inclusive."""

    name: str
    start: int
    stop: int

    def __post_init__(self) -> None:
        """Store `start` and `stop` as plain ints."""
        # The dataclass is frozen, so the constructor is the one writer of the normalized values.
        object.__setattr__(self, "start", as_limit("start", self.start))
        object.__setattr__(self, "stop", as_limit("stop", self.stop))


SeedSpec = Seed | CustomSeed


@dataclass(frozen=True)
class Target:
    """The pattern that a search looks for.

    `sequence` is DNA. `seed_rna` is the reverse complement of `sequence` as RNA. Thus it holds the
    U of the A suffix of `mer8` and `mer7A1`, and a Y for each wobble position.
    """

    name: str
    sequence: str
    seed_rna: str


def _canonical_window(seed: Seed) -> tuple[int | None, str]:
    """Return the last guide position of the window (`None`: the guide end) and the suffix."""
    match seed:
        case Seed.FULL:
            stop, suffix = None, ""
        case Seed.MER18:
            stop, suffix = 19, ""
        case Seed.MER15:
            stop, suffix = 16, ""
        case Seed.MER7M8:
            stop, suffix = 8, ""
        case Seed.MER6:
            stop, suffix = 7, ""
        case Seed.MER8:
            stop, suffix = 8, "A"
        case Seed.MER7A1:
            stop, suffix = 7, "A"
        case _:
            assert_never(seed)
    return stop, suffix


def seed_window(seed: SeedSpec, guide_length: int) -> tuple[int, int]:
    """Return the first and the last guide position, 1-based, that the target of `seed` pairs with.

    The `A` at the end of the target of `Seed.MER8` and `Seed.MER7A1` pairs with g1. The function
    does not check that the seed fits the guide: `Guide.target` does.
    """
    if isinstance(seed, CustomSeed):
        return seed.start, seed.stop
    last, suffix = _canonical_window(seed)
    return (1 if suffix else _CANONICAL_START), (guide_length if last is None else last)


def _check_custom_window(seed: CustomSeed, length: int) -> None:
    prefix = f"custom seed {seed.name!r}"
    if seed.start < 1:
        message = (
            f"{prefix}: start {seed.start} is less than 1; expected a position from 1 to {length}"
        )
        raise InputError(message)
    if seed.stop > length:
        message = (
            f"{prefix}: stop {seed.stop} is beyond the guide length {length}; "
            f"expected a position from {seed.start} to {length}"
        )
        raise InputError(message)
    if seed.start > seed.stop:
        message = f"{prefix}: start {seed.start} is after stop {seed.stop}; expected start <= stop"
        raise InputError(message)


@dataclass(frozen=True)
class Guide:
    """The guide strand as RNA, 5' to 3'. The letters are A, C, G, and U, in either case."""

    sequence: str

    def __post_init__(self) -> None:
        """Reject a letter that is not A, C, G, or U, and store the guide in uppercase."""
        bad = _NOT_RNA_LETTER.search(self.sequence)
        if bad is not None:
            letter, position = bad.group(), bad.start() + 1
            hint = "use 'U', not 'T', for uracil; " if letter in "Tt" else ""
            message = (
                f"invalid guide letter {letter!r} at position {position}; {hint}"
                "expected A, C, G, or U in either case"
            )
            raise InputError(message)
        # The dataclass is frozen, so this is the one write of the field, before any reader sees the
        # instance. The check above leaves ASCII letters only, so upper() keeps their positions.
        object.__setattr__(self, "sequence", self.sequence.upper())

    def target(self, seed: SeedSpec, *, wobble: bool = False) -> Target:
        """Return the DNA pattern that finds the window of the guide that `seed` names.

        With `wobble`, each U of the guide becomes Y before the window is cut. Thus the target has R
        at each such position, and only a search with `fixed_pattern=False` matches an A or a G
        there.
        """
        if not isinstance(seed, Seed | CustomSeed):
            message = (
                f"seed {shown(seed)} has type {type(seed).__name__}; "
                "expected a Seed or a CustomSeed"
            )
            raise InputError(message)
        wobble = as_flag("wobble", wobble)
        length = len(self.sequence)
        if isinstance(seed, CustomSeed):
            _check_custom_window(seed, length)
            name, start, stop, suffix = seed.name, seed.start, seed.stop, ""
        else:
            last, suffix = _canonical_window(seed)
            minimum = _MIN_GUIDE_LENGTH if last is None else max(_MIN_GUIDE_LENGTH, last)
            if length < minimum:
                message = f"{seed.value} needs a guide of at least {minimum} nt, got {length}"
                raise InputError(message)
            name, start, stop = seed.value, _CANONICAL_START, length if last is None else last
        guide = self.sequence.replace("U", "Y") if wobble else self.sequence
        window = guide[start - 1 : stop].replace("U", "T")
        sequence = reverse_complement(window) + suffix
        return Target(name, sequence, reverse_complement(sequence).replace("T", "U"))
