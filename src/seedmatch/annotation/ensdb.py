"""The EnsDb reader: the transcripts of an ensembldb SQLite file, selected in SQL."""

import json
import operator
import os
import sqlite3
from array import array
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from types import NoneType
from typing import Final, Self, cast

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from seedmatch._frames import with_text_dtype
from seedmatch._io import COORDINATE_MAX, limit_text, local_path
from seedmatch.annotation.filters import STANDARD_CHROMOSOMES, TranscriptFilter
from seedmatch.errors import InputError

_REQUIRED_COLUMNS: Final = {
    "gene": ("gene_id", "gene_name", "gene_biotype", "seq_name", "seq_strand"),
    "tx": (
        "tx_id",
        "tx_biotype",
        "tx_seq_start",
        "tx_seq_end",
        "tx_cds_seq_start",
        "tx_cds_seq_end",
        "gene_id",
    ),
    "tx2exon": ("tx_id", "exon_id", "exon_idx"),
    "exon": ("exon_id", "exon_seq_start", "exon_seq_end"),
}
_OPTIONAL_COLUMNS: Final = {
    "chromosome": ("seq_name", "seq_length"),
    "entrezgene": ("gene_id", "entrezid"),
    "protein": ("tx_id", "protein_id"),
}
_COORDINATE: Final = (1, COORDINATE_MAX)
# The valid values of each integer column that the checks read: the lowest and the highest value,
# or `None` for each integer, and if NULL is valid. `exon_rows` reads a NULL CDS start or end as 0,
# so 0 is valid for the CDS columns.
_INTEGER_COLUMNS: Final[dict[str, tuple[tuple[int, int] | None, bool]]] = {
    "tx.tx_seq_start": (_COORDINATE, False),
    "tx.tx_seq_end": (_COORDINATE, False),
    "tx.tx_cds_seq_start": ((0, COORDINATE_MAX), True),
    "tx.tx_cds_seq_end": ((0, COORDINATE_MAX), True),
    "tx.tx_support_level": (None, True),
    "exon.exon_seq_start": (_COORDINATE, False),
    "exon.exon_seq_end": (_COORDINATE, False),
}
# The position of each integer column in a row of the query of `exon_rows` and of `transcripts`.
_EXON_ROW_INTEGERS: Final = {
    4: "tx.tx_cds_seq_start",
    5: "tx.tx_cds_seq_end",
    6: "exon.exon_seq_start",
    7: "exon.exon_seq_end",
    8: "tx.tx_seq_start",
}
_TRANSCRIPT_INTEGERS: Final = {
    4: "tx.tx_seq_start",
    5: "tx.tx_seq_end",
    9: "tx.tx_cds_seq_start",
    10: "tx.tx_cds_seq_end",
    11: "tx.tx_support_level",
}
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
_TEXT_COLUMNS: Final = tuple(name for name in _COLUMNS if name not in _DTYPES)
STRAND_SYMBOLS: Final = {1: "+", -1: "-"}


@dataclass(frozen=True, eq=False)
class ExonRows:
    """The exons of the transcripts of a filter, one entry of each column for each exon.

    The rows come in the order of `tx_id`, then of `exon_idx`, thus the rows of a transcript are
    adjacent and in transcript order. `seq_name` and `strand` (1 or -1) are the columns of the gene.
    `cds_start` and `cds_end` include the stop codon. They are 0 for a transcript without a CDS,
    because a coordinate starts at 1. Each coordinate is at most `COORDINATE_MAX`. `tx_row` is the
    row order of the `tx` table.
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
        """Make sure that the database has the EnsDb tables and the columns that the queries read.

        An optional table (`chromosome`, `entrezgene`, `protein`) that the file has must have the
        columns that the queries read. Then read which optional columns `tx` has.
        """
        self._connection = connection
        self.source = source
        self._tables = frozenset(
            name
            for (name,) in connection.execute("SELECT name FROM sqlite_master WHERE type = 'table'")
        )
        for table in _REQUIRED_COLUMNS:
            if table not in self._tables:
                message = (
                    f"{source!r}: the file has no table {table!r}; "
                    f"expected the EnsDb tables {', '.join(_REQUIRED_COLUMNS)}"
                )
                raise InputError(message)
        checked = _REQUIRED_COLUMNS | {
            table: required
            for table, required in _OPTIONAL_COLUMNS.items()
            if table in self._tables
        }
        columns = {
            table: frozenset(
                name
                for (name,) in connection.execute("SELECT name FROM pragma_table_info(?)", (table,))
            )
            for table in checked
        }
        for table, required in checked.items():
            for column in required:
                if column not in columns[table]:
                    raise _missing_column(
                        source, table, column, f"the EnsDb columns {', '.join(required)}"
                    )
        self._tx_columns = columns["tx"]

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
            "exon.exon_seq_end, tx.tx_seq_start, tx.rowid, exon.exon_id "
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
            columns = tuple(zip(*rows, strict=True))
            self._check_exon_columns(columns)
            (
                tx_ids,
                gene_ids,
                seq_names,
                strands,
                cds_starts,
                cds_ends,
                starts,
                ends,
                tx_seq_starts,
                tx_rows,
                _exon_ids,
            ) = columns
            for text, values in zip(texts, (tx_ids, gene_ids, seq_names), strict=True):
                text.extend(shared.setdefault(value, value) for value in values)
            for number, values in zip(
                numbers,
                (strands, cds_starts, cds_ends, starts, ends, tx_seq_starts, tx_rows),
                strict=True,
            ):
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

    def protein_genes(self, protein_ids: Sequence[str]) -> dict[str, str]:
        """Map each protein ID of the table `protein` to the gene of its transcript.

        A protein ID that the table lacks is not in the result. A protein of a transcript without a
        gene (`tx.gene_id` is NULL) is not in the result either. An annotation without a `protein`
        table raises `InputError`. An example is an annotation of `load_gtf` without `protein_id`
        attributes. A protein with more than one gene raises `InputError` too. A protein of two
        transcripts of one gene is valid. The IDs go to SQLite as one JSON array, so their number
        has no limit of SQLite. The query uses the index `protein_protein_id_idx` of the `protein`
        table, and the memory is O(P) for P IDs.
        """
        if "protein" not in self._tables:
            message = (
                f"{self.source!r}: the annotation has no table 'protein'; "
                "expected an EnsDb file with this table, "
                "or a GTF with the attribute protein_id on its CDS lines"
            )
            raise InputError(message)
        # The text holds only fixed names, and the IDs are one parameter.
        sql = (
            "SELECT protein.protein_id, tx.gene_id FROM protein "  # noqa: S608
            "JOIN tx ON tx.tx_id = protein.tx_id "
            f"WHERE protein.protein_id {_IN_SET} AND tx.gene_id IS NOT NULL"
        )
        genes: dict[str, str] = {}
        for protein_id, gene_id in self._connection.execute(sql, (json.dumps(list(protein_ids)),)):
            known = genes.setdefault(protein_id, gene_id)
            if known != gene_id:
                first, second = sorted((known, gene_id))
                message = (
                    f"{self.source!r}: protein {protein_id!r} has the genes {first!r} and "
                    f"{second!r}; expected one gene for each protein"
                )
                raise InputError(message)
        return genes

    def close(self) -> None:
        """Close the connection."""
        self._connection.close()

    def _require_tx_column(self, column: str, instead: str) -> None:
        if column not in self._tx_columns:
            raise _missing_column(
                self.source, "tx", column, f"a file with this column, or {instead}"
            )

    def _frame(self, rows: list[tuple[object, ...]]) -> pd.DataFrame:
        if rows:
            columns = tuple(zip(*rows, strict=True))
            tx_ids = columns[0]

            def place(index: int) -> str:
                return f"transcript {tx_ids[index]!r}"

            self._check_integers(columns, _TRANSCRIPT_INTEGERS, place)
            self._check_starts_not_after_ends(columns[4], columns[5], place)
        frame = pd.DataFrame.from_records(rows, columns=_COLUMNS)
        strand = frame["strand"].map(STRAND_SYMBOLS)
        unknown = strand.isna()
        if unknown.any():
            position = int(unknown.to_numpy().argmax())
            raise self._strand_error(frame["tx_id"].iloc[position], frame["strand"].iloc[position])
        return with_text_dtype(frame.assign(strand=strand).astype(_DTYPES), _TEXT_COLUMNS)

    def _check_exon_columns(self, columns: Sequence[tuple[object, ...]]) -> None:
        """Raise `InputError` for a number of a chunk of `exon_rows` that is not valid.

        The strand must be 1 or -1, each integer column must pass `_check_integers`, and the exon
        start must not be after its end.
        """
        tx_ids, strands, starts, ends, exon_ids = (columns[i] for i in (0, 3, 6, 7, 10))
        if not (set(map(type, strands)) <= {int} and set(strands) <= STRAND_SYMBOLS.keys()):
            index = next(
                i
                for i, strand in enumerate(strands)
                if type(strand) is not int or strand not in STRAND_SYMBOLS
            )
            raise self._strand_error(tx_ids[index], strands[index])

        def place(index: int) -> str:
            return f"exon {exon_ids[index]!r} of transcript {tx_ids[index]!r}"

        self._check_integers(columns, _EXON_ROW_INTEGERS, place)
        self._check_starts_not_after_ends(starts, ends, place)

    def _check_integers(
        self,
        columns: Sequence[tuple[object, ...]],
        positions: Mapping[int, str],
        place: Callable[[int], str],
    ) -> None:
        """Raise `InputError` for the first value of an integer column that is not valid.

        `positions` maps the position of each integer column to its name in `_INTEGER_COLUMNS`.
        `place` names the row of an index for the message.
        """
        for position, column in positions.items():
            values = columns[position]
            index = _first_bad_integer(values, column)
            if index is not None:
                limits, nullable = _INTEGER_COLUMNS[column]
                accepted = "NULL or an integer" if nullable else "an integer"
                if limits is not None:
                    accepted += f" from {limits[0]} to {limit_text(limits[1])}"
                message = (
                    f"{self.source!r}: {place(index)} has {column} {values[index]!r}; "
                    f"expected {accepted}"
                )
                raise InputError(message)

    def _check_starts_not_after_ends(
        self,
        starts: tuple[object, ...],
        ends: tuple[object, ...],
        place: Callable[[int], str],
    ) -> None:
        """Raise `InputError` for the first start after its end. `place` names the row."""
        # `_check_integers` runs first, and it makes sure that each start and end is an int.
        first, last = cast("tuple[int, ...]", starts), cast("tuple[int, ...]", ends)
        if not all(map(operator.le, first, last)):
            index = next(
                i for i, (start, end) in enumerate(zip(first, last, strict=True)) if start > end
            )
            message = (
                f"{self.source!r}: {place(index)} has start {first[index]} and end "
                f"{last[index]}; expected integers with start <= end"
            )
            raise InputError(message)

    def _strand_error(self, tx_id: object, strand: object) -> InputError:
        message = f"{self.source!r}: transcript {tx_id!r} has strand {strand}; expected 1 or -1"
        return InputError(message)


def _first_bad_integer(values: tuple[object, ...], column: str) -> int | None:
    """Return the index of the first value that is not valid for the column, or `None`.

    A valid value is an integer inside the limits of the column in `_INTEGER_COLUMNS`, or NULL if
    the column permits NULL.
    """
    limits, nullable = _INTEGER_COLUMNS[column]
    lowest, highest = limits or (None, None)
    kinds = set(map(type, values))
    if kinds <= ({int, NoneType} if nullable else {int}):
        if lowest is None or highest is None:
            return None
        present = [value for value in values if value is not None] if NoneType in kinds else values
        # The types of the present values are int only, so min and max compare integers.
        numbers = cast("Sequence[int]", present)
        if not numbers or (lowest <= min(numbers) and max(numbers) <= highest):
            return None
    return next(
        i
        for i, value in enumerate(values)
        if not (value is None and nullable)
        and not (
            type(value) is int and (lowest is None or highest is None or lowest <= value <= highest)
        )
    )


def _missing_column(source: str, table: str, column: str, expected: str) -> InputError:
    message = f"{source!r}: the table {table!r} has no column {column!r}; expected {expected}"
    return InputError(message)
