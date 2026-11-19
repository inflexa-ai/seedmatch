"""Sequences: the IUPAC alphabet as masks, `SequenceSet`, and FASTA files."""

from seedmatch.sequence.alphabet import reverse_complement
from seedmatch.sequence.fasta import read_fasta, write_fasta
from seedmatch.sequence.sequence_set import SequenceSet

__all__ = ["SequenceSet", "read_fasta", "reverse_complement", "write_fasta"]
