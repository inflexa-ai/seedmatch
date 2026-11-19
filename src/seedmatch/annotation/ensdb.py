"""The EnsDb reader: the transcripts of an ensembldb SQLite file, selected in SQL."""

import json
import os
import sqlite3
from array import array
from dataclasses import dataclass
from typing import Final, Self

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from seedmatch._io import local_path
from seedmatch.annotation.filters import STANDARD_CHROMOSOMES, TranscriptFilter
from seedmatch.errors import InputError

_REQUIRED_TABLES: Final = ("gene", "tx", "tx2exon", "exon")
# A result becomes a frame chunk by chunk, so Python never holds the tuples of all rows at once.
_CHUNK_ROWS: Final = 50_000
_STANDARD_CHROMOSOMES_JSON: Final = json.dumps(sorted(STANDARD_CHROMOSOMES))
_IN_SET: Final = "IN (SELECT value FROM json_each(?))"
_ENTREZ_IDS: Final = (
    "tx.gene_id IN (SELECT gene_id FROM entrezgene "
    "WHERE entrezid IN (SELECT value FROM json_each(?)))"
)
_COLUMNS: Final = [
    "tx_id",
    "gene_id",
    "gene_name",
    "contig",
    "start",
    "end",
    "strand",
    "tx_biotype",
    "gene_biotype",
    "cds_start",
    "cds_end",
    "support_level",
    "canonical",
]
_DTYPES: Final = {
    "start": "int64",
    "end": "int64",
    "cds_start": "Int64",
    "cds_end": "Int64",
    "support_level": "Int64",
    "canonical": "boolean",
}
_STRANDS: Final = {1: "+", -1: "-"}


@dataclass(frozen=True, eq=False)
class ExonRows:
    """The exons of the transcripts of a filter, one entry of each column for each exon.

    The rows come in the order of `tx_id`, then of `exon_idx`, thus the rows of a transcript are
    adjacent and in transcript order. `seq_name` and `strand` (1 or -1) are the columns of the gene.
    `cds_start` and `cds_end` include the stop codon. They are 0 for a transcript without a CDS,
    because a coordinate starts at 1. `tx_row` is the row order of the `tx` table.
    """

    tx_id: NDArray[np.object_]
    gene_id: NDArray[np.object_]
    seq_name: NDArray[np.object_]
    strand: NDArray[np.int8]
    cds_start: NDArray[np.int64]
    cds_end: NDArray[np.int64]
    exon_start: NDArray[np.int64]
    exon_end: NDArray[np.int64]
    tx_seq_start: NDArray[np.int64]
    tx_row: NDArray[np.int64]


class AnnotationDb:
    """The tables of an EnsDb in an SQLite database. A filter is a condition of the SQL query.

    The connection can be any `sqlite3.Connection` with the tables `gene`, `tx`, `tx2exon`, and
    `exon`. The object owns the connection, and `close` closes it. `source` names the database in
    the messages.
    """

    def __init__(self, connection: sqlite3.Connection, source: str) -> None:
        """Make sure that the database has the tables of an EnsDb, and read its optional columns."""
        self._connection = connection
        self.source = source
        self._tables = frozenset(
            name
            for (name,) in connection.execute("SELECT name FROM sqlite_master WHERE type = 'table'")
        )
        for table in _REQUIRED_TABLES:
            if table not in self._tables:
                message = (
                    f"{source!r}: the file has no table {table!r}; "
                    f"expected the EnsDb tables {', '.join(_REQUIRED_TABLES)}"
                )
                raise InputError(message)
        self._tx_columns = frozenset(
            name for _, name, *_ in connection.execute("PRAGMA table_info(tx)")
        )

    @classmethod
    def open(cls, path: str | os.PathLike[str]) -> Self:
        """Open the file read-only and immutable, so that the file can be in a read-only store.

        A URL and a file that is not an EnsDb raise `InputError`, and a missing file raises
        `FileNotFoundError`.
        """
        local = local_path(path)
        source = os.fspath(local)
        # With mode=ro, SQLite reports a missing file only as "unable to open database file".
        if not local.exists():
            message = f"{source!r}: no such file; expected an EnsDb SQLite file"
            raise FileNotFoundError(message)
        try:
            connection = sqlite3.connect(
                f"{local.resolve().as_uri()}?mode=ro&immutable=1", uri=True
            )
            try:
                return cls(connection, source)
            except Exception:
                connection.close()
                raise
        except sqlite3.DatabaseError as error:
            message = f"{source!r}: SQLite cannot read the file ({error}); expected an EnsDb file"
            raise InputError(message) from error

    def contig_lengths(self) -> dict[str, int]:
        """Map the name of each contig to its length, or give `{}` without a `chromosome` table."""
        if "chromosome" not in self._tables:
            return {}
        rows = self._connection.execute("SELECT seq_name, seq_length FROM chromosome")
        return dict(rows)

    def contig_names(self) -> list[str]:
        """Return the contigs of the genes, each one once, in sorted order."""
        rows = self._connection.execute("SELECT DISTINCT seq_name FROM gene ORDER BY seq_name")
        return [name for (name,) in rows]

    def where(self, transcript_filter: TranscriptFilter) -> tuple[str, list[object]]:
        """Return the SQL condition of the filter and its parameters, in the order of the `?`.

        The condition names the tables `tx` and `gene`, so it fits a query that joins them. It is
        `1` for a filter without a condition. Each value is a parameter, and each set of values is
        one JSON array that `json_each` reads, so the number of values has no limit of SQLite.
        """
        conditions: list[str] = []
        parameters: list[object] = []
        if transcript_filter.standard_chromosomes:
            conditions.append(f"gene.seq_name {_IN_SET}")
            parameters.append(_STANDARD_CHROMOSOMES_JSON)
        if transcript_filter.require_tx_id:
            conditions.append("tx.tx_id != ''")
        if transcript_filter.protein_coding:
            conditions.append(
                "gene.gene_biotype = 'protein_coding' AND tx.tx_biotype = 'protein_coding'"
            )
        if transcript_filter.canonical:
            self._require_tx_column("tx_is_canonical", "canonical=False")
            conditions.append("tx.tx_is_canonical = 1")
        if transcript_filter.support_level is not None:
            self._require_tx_column("tx_support_level", "support_level=None")
            conditions.append("tx.tx_support_level = ?")
            parameters.append(transcript_filter.support_level)
        for column, values in (
            ("tx.gene_id", transcript_filter.gene_ids),
            ("tx.tx_id", transcript_filter.tx_ids),
            ("gene.gene_name", transcript_filter.symbols),
            ("tx.tx_biotype", transcript_filter.tx_biotypes),
            ("gene.gene_biotype", transcript_filter.gene_biotypes),
        ):
            if values is not None:
                conditions.append(f"{column} {_IN_SET}")
                parameters.append(json.dumps(values))
        if transcript_filter.entrez_ids is not None:
            if "entrezgene" not in self._tables:
                message = (
                    f"{self.source!r}: the file has no table 'entrezgene'; "
                    "expected a file with this table, or entrez_ids=None"
                )
                raise InputError(message)
            conditions.append(_ENTREZ_IDS)
            parameters.append(json.dumps(transcript_filter.entrez_ids))
        return " AND ".join(f"({condition})" for condition in conditions) or "1", parameters

    def transcripts(self, transcript_filter: TranscriptFilter) -> pd.DataFrame:
        """Return the transcripts that pass the filter, one row for each transcript.

        The columns are `tx_id`, `gene_id`, `gene_name`, `contig`, `start`, `end`, `strand`,
        `tx_biotype`, `gene_biotype`, `cds_start`, `cds_end`, `support_level`, and `canonical`. A
        column that the file lacks (`support_level`, `canonical`) has missing values, and so does
        a transcript without a CDS. The rows come in the order of `contig`, `start`, and then the
        row order of the tables.
        """
        condition, parameters = self.where(transcript_filter)
        support = "tx.tx_support_level" if "tx_support_level" in self._tx_columns else "NULL"
        canonical = "(tx.tx_is_canonical = 1)" if "tx_is_canonical" in self._tx_columns else "NULL"
        # ensembldb sorts on the contig and the start only. The row orders make the order of a tie
        # the same for each filter (D-030). The text holds only fixed names, and each value of the
        # filter is a parameter.
        sql = (
            "SELECT tx.tx_id, tx.gene_id, gene.gene_name, gene.seq_name, tx.tx_seq_start, "  # noqa: S608
            "tx.tx_seq_end, gene.seq_strand, tx.tx_biotype, gene.gene_biotype, "
            f"tx.tx_cds_seq_start, tx.tx_cds_seq_end, {support}, {canonical} "
            "FROM tx JOIN gene ON tx.gene_id = gene.gene_id "
            f"WHERE {condition} "
            "ORDER BY gene.seq_name, tx.tx_seq_start, gene.rowid, tx.rowid"
        )
        cursor = self._connection.execute(sql, parameters)
        frames = []
        while rows := cursor.fetchmany(_CHUNK_ROWS):
            frames.append(self._frame(rows))
        return pd.concat(frames, ignore_index=True) if frames else self._frame([])

    def exon_rows(self, transcript_filter: TranscriptFilter, *, coding_only: bool) -> ExonRows:
        """Return the exons of the transcripts that pass the filter, as `ExonRows`.

        With `coding_only`, the rows hold only the transcripts with a CDS.
        """
        condition, parameters = self.where(transcript_filter)
        if coding_only:
            condition += " AND tx.tx_cds_seq_start IS NOT NULL"
        # The text holds only fixed names, and each value of the filter is a parameter.
        sql = (
            "SELECT tx.tx_id, tx.gene_id, gene.seq_name, gene.seq_strand, "  # noqa: S608
            "IFNULL(tx.tx_cds_seq_start, 0), IFNULL(tx.tx_cds_seq_end, 0), exon.exon_seq_start, "
            "exon.exon_seq_end, tx.tx_seq_start, tx.rowid "
            "FROM tx JOIN gene ON tx.gene_id = gene.gene_id "
            "JOIN tx2exon ON tx2exon.tx_id = tx.tx_id "
            "JOIN exon ON exon.exon_id = tx2exon.exon_id "
            f"WHERE {condition} "
            "ORDER BY tx.tx_id, tx2exon.exon_idx"
        )
        cursor = self._connection.execute(sql, parameters)
        texts: tuple[list[str], ...] = ([], [], [])
        numbers = tuple(array("q") for _ in range(7))
        # The rows of a transcript repeat its text values. One str object for each value keeps the
        # memory near the size of the integer columns.
        shared: dict[str, str] = {}
        while rows := cursor.fetchmany(_CHUNK_ROWS):
            tx_ids, gene_ids, seq_names, strands, *others = zip(*rows, strict=True)
            if not set(strands) <= _STRANDS.keys():
                position = next(i for i, strand in enumerate(strands) if strand not in _STRANDS)
                raise self._strand_error(tx_ids[position], strands[position])
            for text, values in zip(texts, (tx_ids, gene_ids, seq_names), strict=True):
                text.extend(shared.setdefault(value, value) for value in values)
            for number, values in zip(numbers, (strands, *others), strict=True):
                number.extend(values)
        tx_id, gene_id, seq_name = (np.array(text, dtype=object) for text in texts)
        strand, cds_start, cds_end, exon_start, exon_end, tx_seq_start, tx_row = (
            np.array(number, dtype=np.int64) for number in numbers
        )
        return ExonRows(
            tx_id=tx_id,
            gene_id=gene_id,
            seq_name=seq_name,
            strand=strand.astype(np.int8),
            cds_start=cds_start,
            cds_end=cds_end,
            exon_start=exon_start,
            exon_end=exon_end,
            tx_seq_start=tx_seq_start,
            tx_row=tx_row,
        )

    def close(self) -> None:
        """Close the connection."""
        self._connection.close()

    def _require_tx_column(self, column: str, instead: str) -> None:
        if column not in self._tx_columns:
            message = (
                f"{self.source!r}: the table 'tx' has no column {column!r}; "
                f"expected a file with this column, or {instead}"
            )
            raise InputError(message)

    def _frame(self, rows: list[tuple[object, ...]]) -> pd.DataFrame:
        frame = pd.DataFrame.from_records(rows, columns=_COLUMNS)
        strand = frame["strand"].map(_STRANDS)
        unknown = strand.isna()
        if unknown.any():
            position = int(unknown.to_numpy().argmax())
            raise self._strand_error(frame["tx_id"].iloc[position], frame["strand"].iloc[position])
        return frame.assign(strand=strand).astype(_DTYPES)

    def _strand_error(self, tx_id: object, strand: object) -> InputError:
        message = f"{self.source!r}: transcript {tx_id!r} has strand {strand}; expected 1 or -1"
        return InputError(message)
