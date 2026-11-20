"""Open a genome file of either format."""

import os
from typing import Final

from seedmatch._io import local_path, open_binary
from seedmatch.annotation.fasta_genome import GZIP_MAGIC, FastaGenome
from seedmatch.annotation.genome import Genome
from seedmatch.annotation.twobit import SIGNATURES, TwoBitGenome
from seedmatch.errors import InputError

_FASTA_START: Final = b">"
_HEAD_SIZE: Final = 4


def open_genome(path: str | os.PathLike[str]) -> Genome:
    """Open the genome file, and select the reader from the first bytes of the file.

    The 2bit signature in either byte order gives a `TwoBitGenome`. The gzip signature or `>` gives
    a `FastaGenome`. Other bytes raise `InputError`.
    """
    local = local_path(path)
    with open_binary(local) as handle:
        head = handle.read(_HEAD_SIZE)
    if not head:
        message = f"{os.fspath(local)!r}: the file is empty; expected a 2bit or FASTA file"
        raise InputError(message)
    if head in SIGNATURES.values():
        return TwoBitGenome(local)
    if head.startswith((GZIP_MAGIC, _FASTA_START)):
        return FastaGenome(local)
    signatures = " or ".join(signature.hex() for signature in SIGNATURES.values())
    message = (
        f"{os.fspath(local)!r}: not a 2bit or FASTA file; the first bytes are {head.hex()}, "
        f"expected {signatures} (2bit), {GZIP_MAGIC.hex()} (gzip), or {_FASTA_START.hex()} (FASTA)"
    )
    raise InputError(message)
