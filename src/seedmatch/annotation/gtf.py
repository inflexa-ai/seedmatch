"""The GTF reader: an Ensembl or GENCODE GTF becomes an in-memory database with the EnsDb tables."""

import os
import re
import sqlite3
from collections.abc import Iterable
from typing import Final, NamedTuple, NoReturn

from seedmatch._io import local_path, open_text
from seedmatch.annotation.ensdb import AnnotationDb
from seedmatch.errors import InputError

_FIELDS: Final = 9
_BATCH_ROWS: Final = 10_000
# One pair of a key and a value. A quoted value can hold a semicolon, so the quotes come first.
_ATTRIBUTE: Final = re.compile(r'([^\s;"]+)\s+(?:"([^"]*)"|([^\s;"]*))')
_SUPPORT_LEVEL: Final = re.compile(r"\d+")
_NOT_ASSIGNED: Final = re.compile(r"NA\b")
_CANONICAL_TAG: Final = "Ensembl_canonical"
_STRANDS: Final = {"+": 1, "-": -1}
_SCHEMA: Final = """
CREATE TABLE gene (
  gene_id TEXT, gene_name TEXT, gene_biotype TEXT, seq_name TEXT, seq_strand INTEGER
);
CREATE TABLE tx (
  tx_id TEXT, tx_biotype TEXT, tx_seq_start INTEGER, tx_seq_end INTEGER,
  tx_cds_seq_start INTEGER, tx_cds_seq_end INTEGER, gene_id TEXT,
  tx_support_level INTEGER, tx_is_canonical INTEGER
);
CREATE TABLE exon (exon_id TEXT, exon_seq_start INTEGER, exon_seq_end INTEGER);
CREATE TABLE tx2exon (tx_id TEXT, exon_id TEXT, exon_idx INTEGER);
CREATE TABLE gtf_exon (
  tx_id TEXT, exon_id TEXT, exon_seq_start INTEGER, exon_seq_end INTEGER, line INTEGER
);
"""
# The statements run after the last line, so that the inserts do not maintain an index.
_UNIQUE_IDS: Final = """
CREATE UNIQUE INDEX gene_gene_id_idx ON gene (gene_id);
CREATE UNIQUE INDEX tx_tx_id_idx ON tx (tx_id);
CREATE INDEX tx_gene_id_idx ON tx (gene_id);
"""
_TX2EXON: Final = """
INSERT INTO tx2exon
SELECT staged.tx_id, staged.exon_id, ROW_NUMBER() OVER (
    PARTITION BY staged.tx_id
    ORDER BY CASE gene.seq_strand WHEN 1 THEN staged.exon_seq_start ELSE -staged.exon_seq_start END,
             staged.line
)
FROM gtf_exon AS staged
JOIN tx ON tx.tx_id = staged.tx_id
JOIN gene ON gene.gene_id = tx.gene_id
ORDER BY staged.line
"""


class _Gene(NamedTuple):
    name: str
    biotype: str | None
    seq_name: str
    strand: int


class _Transcript(NamedTuple):
    biotype: str | None
    start: int
    end: int
    gene_id: str
    support_level: int | None
    canonical: bool


def load_gtf(path: str | os.PathLike[str]) -> AnnotationDb:
    """Read an Ensembl or GENCODE GTF into an in-memory database with the tables of an EnsDb.

    The file is UTF-8 text, or gzip, bzip2, or xz of UTF-8 text. The reader uses only the lines
    `transcript`, `exon`, `CDS`, and `stop_codon`, and it skips the lines that start with `#`. The
    rows of the `gene` and `tx` tables come in the line order of the first transcript line of
    each gene and of each transcript. A `gene` row holds the name and the biotype of the first
    transcript line of its gene. The CDS of a transcript goes from the first base to the last base
    of its `CDS` and `stop_codon` lines. The columns `tx_is_canonical` and `tx_support_level` are
    in the database only if some transcript has the tag `Ensembl_canonical` or the attribute
    `transcript_support_level`. The database has no `chromosome` table and no `entrezgene` table.
    The reader writes no file.

    A URL, a line that is not valid, a repeated transcript, and an exon without a transcript line
    raise `InputError`, and a missing file raises `FileNotFoundError`.
    """
    local = local_path(path)
    source = os.fspath(local)
    if not local.exists():
        message = f"{source!r}: no such file; expected a GTF file"
        raise FileNotFoundError(message)
    connection = sqlite3.connect(":memory:")
    try:
        # The sorts of the queries below stay in memory, so that the reader makes no temporary file.
        connection.execute("PRAGMA temp_store = MEMORY")
        connection.executescript(_SCHEMA)
        reader = _Reader(connection, source)
        try:
            with open_text(local) as lines:
                reader.read(lines)
        except UnicodeDecodeError as exc:
            message = (
                f"{source!r}: the file is not valid UTF-8 text ({exc.reason}); "
                "expected a text or gzip GTF file"
            )
            raise InputError(message) from exc
        reader.finish()
        return AnnotationDb(connection, source)
    except Exception:
        connection.close()
        raise


class _Reader:
    """The state of one pass over the lines of a GTF."""

    def __init__(self, connection: sqlite3.Connection, source: str) -> None:
        self._connection = connection
        self._source = source
        self._genes: dict[str, _Gene] = {}
        self._transcripts: dict[str, _Transcript] = {}
        self._cds: dict[str, tuple[int, int]] = {}
        self._exons: list[tuple[str, str, int, int, int]] = []
        self._has_support_level = False
        self._has_canonical = False

    def read(self, lines: Iterable[str]) -> None:
        """Parse each line once. The exons go to the database in batches of `_BATCH_ROWS`."""
        for number, line in enumerate(lines, start=1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if fields == [""]:
                continue
            if len(fields) != _FIELDS:
                self._fail(
                    number, f"the line has {len(fields)} tab-separated fields; expected {_FIELDS}"
                )
            match fields[2]:
                case "transcript":
                    self._transcript(number, fields)
                case "exon":
                    self._exon(number, fields)
                case "CDS" | "stop_codon":
                    self._coding(number, fields)
        self._flush()

    def finish(self) -> None:
        """Write the genes, the transcripts, and the exon tables, then the indexes."""
        connection = self._connection
        connection.executemany(
            "INSERT INTO gene VALUES (?, ?, ?, ?, ?)",
            ((gene_id, *gene) for gene_id, gene in self._genes.items()),
        )
        connection.executemany(
            "INSERT INTO tx VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (
                (
                    tx_id,
                    tx.biotype,
                    tx.start,
                    tx.end,
                    *self._cds.get(tx_id, (None, None)),
                    tx.gene_id,
                    tx.support_level,
                    tx.canonical,
                )
                for tx_id, tx in self._transcripts.items()
            ),
        )
        connection.executescript(_UNIQUE_IDS)
        self._write_exons()
        for present, column in (
            (self._has_support_level, "tx_support_level"),
            (self._has_canonical, "tx_is_canonical"),
        ):
            if not present:
                connection.execute(f"ALTER TABLE tx DROP COLUMN {column}")
        connection.commit()

    def _transcript(self, number: int, fields: list[str]) -> None:
        start, end = self._span(number, fields)
        attributes = _attributes(fields[8])
        tx_id = self._required(number, attributes, "transcript_id", "transcript")
        gene_id = self._required(number, attributes, "gene_id", "transcript")
        strand = _STRANDS.get(fields[6])
        if strand is None:
            self._fail(number, f"strand {fields[6]!r}; expected + or -")
        if tx_id in self._transcripts:
            self._fail(
                number,
                f"transcript {tx_id!r} occurs twice; "
                "expected one transcript line for each transcript_id",
            )
        gene = self._genes.get(gene_id)
        if gene is None:
            self._genes[gene_id] = _Gene(
                attributes.get("gene_name", ""),
                attributes.get("gene_biotype") or attributes.get("gene_type"),
                fields[0],
                strand,
            )
        elif (gene.seq_name, gene.strand) != (fields[0], strand):
            self._fail(
                number,
                f"transcript {tx_id!r} is on contig {fields[0]!r} and strand {fields[6]!r}, "
                f"but its gene {gene_id!r} is on contig {gene.seq_name!r} and strand "
                f"{'+' if gene.strand == 1 else '-'!r}, from its first transcript line; "
                "expected one contig and one strand for each gene",
            )
        level = attributes.get("transcript_support_level")
        canonical = _CANONICAL_TAG in _tags(fields[8])
        self._has_support_level = self._has_support_level or level is not None
        self._has_canonical = self._has_canonical or canonical
        self._transcripts[tx_id] = _Transcript(
            attributes.get("transcript_biotype") or attributes.get("transcript_type"),
            start,
            end,
            gene_id,
            None if level is None else self._support_level(number, level),
            canonical,
        )

    def _exon(self, number: int, fields: list[str]) -> None:
        start, end = self._span(number, fields)
        attributes = _attributes(fields[8])
        tx_id = self._required(number, attributes, "transcript_id", "exon")
        exon_id = self._required(number, attributes, "exon_id", "exon")
        self._exons.append((tx_id, exon_id, start, end, number))
        if len(self._exons) >= _BATCH_ROWS:
            self._flush()

    def _coding(self, number: int, fields: list[str]) -> None:
        start, end = self._span(number, fields)
        tx_id = self._required(number, _attributes(fields[8]), "transcript_id", fields[2])
        first, last = self._cds.get(tx_id, (start, end))
        self._cds[tx_id] = (min(first, start), max(last, end))

    def _flush(self) -> None:
        self._connection.executemany("INSERT INTO gtf_exon VALUES (?, ?, ?, ?, ?)", self._exons)
        self._exons.clear()

    def _write_exons(self) -> None:
        """Move the exon lines to `exon` and `tx2exon`, and drop the staging table."""
        connection = self._connection
        orphan = connection.execute(
            "SELECT tx_id, line FROM gtf_exon WHERE tx_id NOT IN (SELECT tx_id FROM tx) "
            "ORDER BY line LIMIT 1"
        ).fetchone()
        if orphan is not None:
            tx_id, number = orphan
            self._fail(
                number,
                f"the exon of transcript {tx_id!r} has no transcript line; "
                "expected a line with the feature transcript",
            )
        conflict = connection.execute(
            "SELECT exon_id FROM gtf_exon GROUP BY exon_id "
            "HAVING COUNT(DISTINCT exon_seq_start) > 1 OR COUNT(DISTINCT exon_seq_end) > 1 "
            "ORDER BY MIN(line) LIMIT 1"
        ).fetchone()
        if conflict is not None:
            self._fail_conflict(conflict[0])
        connection.execute(
            "INSERT INTO exon SELECT exon_id, MIN(exon_seq_start), MIN(exon_seq_end) "
            "FROM gtf_exon GROUP BY exon_id ORDER BY MIN(line)"
        )
        connection.execute("CREATE UNIQUE INDEX exon_exon_id_idx ON exon (exon_id)")
        connection.execute(_TX2EXON)
        connection.execute("CREATE INDEX tx2exon_tx_id_idx ON tx2exon (tx_id)")
        connection.execute("DROP TABLE gtf_exon")

    def _fail_conflict(self, exon_id: str) -> NoReturn:
        rows = self._connection.execute(
            "SELECT exon_seq_start, exon_seq_end, line FROM gtf_exon WHERE exon_id = ? "
            "ORDER BY line",
            (exon_id,),
        ).fetchall()
        first = rows[0]
        other = next(row for row in rows if row[:2] != first[:2])
        self._fail_at(
            f"exon {exon_id!r} has the coordinates {first[0]}-{first[1]} at line {first[2]} "
            f"and {other[0]}-{other[1]} at line {other[2]}; "
            "expected the same coordinates for each line of an exon_id"
        )

    def _span(self, number: int, fields: list[str]) -> tuple[int, int]:
        start, end = fields[3], fields[4]
        if not (start.isdecimal() and end.isdecimal() and 0 < int(start) <= int(end)):
            self._fail(
                number, f"start {start!r} and end {end!r}; expected integers with 1 <= start <= end"
            )
        return int(start), int(end)

    def _required(self, number: int, attributes: dict[str, str], key: str, feature: str) -> str:
        value = attributes.get(key)
        if not value:
            self._fail(
                number, f'the {feature} line has no {key}; expected the attribute {key} "..."'
            )
        return value

    def _support_level(self, number: int, text: str) -> int | None:
        if _NOT_ASSIGNED.match(text):
            return None
        found = _SUPPORT_LEVEL.match(text)
        if found is None:
            self._fail(
                number,
                f"transcript_support_level {text!r}; "
                "expected a text that starts with NA or with an integer",
            )
        return int(found[0])

    def _fail(self, number: int, text: str) -> NoReturn:
        self._fail_at(f"line {number}: {text}")

    def _fail_at(self, text: str) -> NoReturn:
        message = f"{self._source!r}: {text}"
        raise InputError(message)


def _attributes(text: str) -> dict[str, str]:
    """Return the value of each key of the attribute field. A repeated key keeps its last value."""
    return {key: quoted or bare for key, quoted, bare in _ATTRIBUTE.findall(text)}


def _tags(text: str) -> list[str]:
    return [quoted or bare for key, quoted, bare in _ATTRIBUTE.findall(text) if key == "tag"]
