"""Write the synthetic annotation fixture in tests/data/annotation/ from one specification.

The specification below defines each gene, transcript, exon, CDS, N block, and soft-masked run.
The script writes an EnsDb with schema 2.2, 2bit files of version 0 and 1, a FASTA file, an
Ensembl GTF, and a GENCODE GTF. bgzip of htslib and samtools write the BGZF file and its indexes,
thus the FASTA reader of seedmatch reads the files of an independent writer.

Run it from the repository root: uv run python tools/fixtures/annotation.py
"""

import gzip
import shutil
import sqlite3
import struct
import subprocess
import sys
from collections.abc import Callable, Iterator
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

OUT = Path(__file__).resolve().parents[2] / "tests" / "data" / "annotation"
SEED = 20261119
LINE_BASES = 60

type Span = tuple[int, int]

# 1-based inclusive coordinates, as in the EnsDb and the GTF.
# No length is a multiple of 4, thus the last byte of each 2bit record is partial.
CONTIGS = {"1": 150_002, "X": 5003, "SCAF1.1": 2001}
SCAFFOLDS = frozenset({"SCAF1.1"})
GENOME_BUILD = "SYN1"
N_BLOCKS: dict[str, tuple[Span, ...]] = {"1": ((1, 100), (22790, 22809), (64320, 64339))}
SOFT_MASKED: dict[str, tuple[Span, ...]] = {"1": ((22400, 22440), (25350, 25380), (66100, 66130))}


@dataclass(frozen=True)
class Transcript:
    number: int
    biotype: str
    # In transcript order, thus in descending order on the minus strand.
    exons: tuple[Span, ...]
    # The EnsDb columns tx_cds_seq_start and tx_cds_seq_end: the CDS with the stop codon.
    cds: Span | None = None
    support_level: int | None = None
    canonical: bool = False

    @property
    def tx_id(self) -> str:
        return f"ENSSYNT{self.number:011d}"

    @property
    def protein_id(self) -> str:
        return f"ENSSYNP{self.number:011d}"

    @property
    def span(self) -> Span:
        return min(s for s, _ in self.exons), max(e for _, e in self.exons)


@dataclass(frozen=True)
class Gene:
    number: int
    name: str
    biotype: str
    contig: str
    strand: int
    # In the row order of the tx table.
    transcripts: tuple[Transcript, ...]
    entrez: tuple[int, ...] = ()

    @property
    def gene_id(self) -> str:
        return f"ENSSYNG{self.number:011d}"

    @property
    def span(self) -> Span:
        spans = [tx.span for tx in self.transcripts]
        return min(s for s, _ in spans), max(e for _, e in spans)


PC = "protein_coding"

# In the row order of the gene table. The ID numbers give the row order of the exon, tx2exon, and
# entrezgene tables.
GENES = (
    Gene(7, "Syn7", PC, "1", -1, (
        Transcript(10, PC, ((67000, 67300), (64000, 64600)), (64451, 67149), 2, canonical=True),
        Transcript(11, "nonsense_mediated_decay", ((67600, 67900), (67000, 67300)), (67151, 67749)),
    )),
    Gene(1, "Syn1", PC, "1", 1, (
        Transcript(2, PC, ((64000, 64300), (65000, 65200), (66000, 66700)), (64151, 65200)),
        Transcript(1, PC, ((64000, 64300), (65000, 65400), (66000, 66500)), (64151, 65200), 1,
                   canonical=True),
    ), entrez=(1001,)),
    Gene(2, "Syn2", PC, "1", -1, (
        Transcript(3, PC, ((29000, 29200), (27000, 27400), (25000, 25600), (22000, 22800)),
                   (25599, 27201), canonical=True),
        Transcript(4, PC, ((29000, 29200), (26800, 27400), (25300, 25700), (22801, 23000),
                           (21500, 21998)), (26800, 27201)),
    ), entrez=(1002,)),
    Gene(3, "Syn3", "lncRNA", "1", 1, (
        Transcript(5, "lncRNA", ((40001, 40300), (41001, 41500)), canonical=True),
    )),
    Gene(4, "Syn4", PC, "SCAF1.1", 1, (
        Transcript(6, PC, ((101, 600), (900, 1800)), (301, 1199), canonical=True),
    )),
    Gene(5, "Syn5", PC, "X", 1, (
        Transcript(7, PC, ((1001, 1500), (2001, 2600)), (1201, 2600), 1, canonical=True),
        Transcript(8, PC, ((1001, 1500), (2001, 3400)), (1201, 2300), 2),
    ), entrez=(1005,)),
    # An Ensembl EnsDb stores an absent gene name as an empty text, and the GTF omits gene_name.
    Gene(6, "", PC, "1", 1, (
        Transcript(9, PC, ((50001, 50600),), (50101, 50595), canonical=True),
    ), entrez=(1006, 1007)),
    Gene(8, "Syn8", PC, "1", 1, (
        Transcript(12, PC, ((80001, 80400), (81001, 81600)), (80200, 81600), canonical=True),
    )),
)  # fmt: skip

TX_ROW_ORDER = (2, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
T09_UTR3_BASES = 5

METADATA = (
    ("Db type", "EnsDb"),
    ("Type of Gene ID", "Ensembl Gene ID"),
    ("Supporting package", "ensembldb"),
    ("Db created by", "seedmatch tools/fixtures/annotation.py"),
    ("script_version", "1"),
    ("Creation time", "Thu Nov 19 00:00:00 2026"),
    ("ensembl_version", "109"),
    ("ensembl_host", "localhost"),
    ("Organism", "Synthetic organism"),
    ("taxonomy_id", "32630"),
    ("genome_build", GENOME_BUILD),
    ("DBSCHEMAVERSION", "2.2"),
    ("common_name", "synthetic"),
    ("species", "synthetic_organism"),
)

# The tables and indexes of schema 2.2, copied from EnsDb.Rnorvegicus.v109.sqlite. The file has the
# protein table, but not the uniprot and protein_domain tables. Thus hasProteinData of ensembldb is
# FALSE, and its queries do not change.
SCHEMA = """
CREATE TABLE `metadata` (
  `name` TEXT,
  `value` TEXT
);
CREATE TABLE `chromosome` (
  `seq_name` TEXT,
  `seq_length` INTEGER,
  `is_circular` INTEGER
);
CREATE TABLE `gene` (
  `gene_id` TEXT,
  `gene_name` TEXT,
  `gene_biotype` TEXT,
  `gene_seq_start` INTEGER,
  `gene_seq_end` INTEGER,
  `seq_name` TEXT,
  `seq_strand` INTEGER,
  `seq_coord_system` TEXT,
  `description` TEXT,
  `gene_id_version` TEXT,
  `canonical_transcript` TEXT
);
CREATE TABLE `entrezgene` (
  `gene_id` TEXT,
  `entrezid` INTEGER
);
CREATE TABLE `tx` (
  `tx_id` TEXT,
  `tx_biotype` TEXT,
  `tx_seq_start` INTEGER,
  `tx_seq_end` INTEGER,
  `tx_cds_seq_start` INTEGER,
  `tx_cds_seq_end` INTEGER,
  `gene_id` TEXT,
  `tx_support_level` INTEGER,
  `tx_id_version` TEXT,
  `gc_content` REAL,
  `tx_external_name` TEXT,
  `tx_is_canonical` INTEGER
);
CREATE TABLE `exon` (
  `exon_id` TEXT,
  `exon_seq_start` INTEGER,
  `exon_seq_end` INTEGER
);
CREATE TABLE `tx2exon` (
  `tx_id` TEXT,
  `exon_id` TEXT,
  `exon_idx` INTEGER
);
CREATE TABLE `protein` (`tx_id` TEXT, `protein_id` TEXT, `protein_sequence` TEXT);
"""
INDEXES = """
CREATE UNIQUE INDEX chromosome_seq_name_idx on chromosome (seq_name);
CREATE UNIQUE INDEX gene_gene_id_idx on gene (gene_id);
CREATE INDEX gene_gene_name_idx on gene (gene_name);
CREATE INDEX gene_seq_name_idx on gene (seq_name);
CREATE UNIQUE INDEX tx_tx_id_idx on tx (tx_id);
CREATE INDEX tx_gene_id_idx on tx (gene_id);
CREATE UNIQUE INDEX exon_exon_id_idx on exon (exon_id);
CREATE INDEX tx2exon_tx_id_idx on tx2exon (tx_id);
CREATE INDEX tx2exon_exon_id_idx on tx2exon (exon_id);
CREATE INDEX entrezgene_gene_id_idx on entrezgene (gene_id);
CREATE INDEX entrezgene_entrezid_idx on entrezgene (entrezid);
CREATE INDEX tx2exon_exon_idx_idx on tx2exon (exon_idx);
CREATE INDEX protein_tx_id_idx on protein (tx_id);
CREATE INDEX protein_protein_id_idx on protein (protein_id);
"""

TWOBIT_MAGIC = 0x1A412743
TWOBIT_CODES = {"T": 0, "C": 1, "A": 2, "G": 3}
CODON = 3
GENCODE_CONTIGS = {"1": "chr1", "X": "chrX"}
OUTPUTS = (
    "synthetic.sqlite",
    "synthetic.2bit",
    "synthetic.v1.2bit",
    "synthetic.fa",
    "synthetic.fa.fai",
    "synthetic.fa.gz",
    "synthetic.fa.gz.fai",
    "synthetic.fa.gz.gzi",
    "synthetic.gtf",
    "synthetic.gencode.gtf.gz",
)

type Part = tuple[int, Span]  # the 1-based exon number and the part of that exon


def gene_of(number: int) -> Gene:
    return next(g for g in GENES if any(tx.number == number for tx in g.transcripts))


def tx_of(number: int) -> Transcript:
    return next(tx for tx in gene_of(number).transcripts if tx.number == number)


def transcripts_by_id() -> list[tuple[Gene, Transcript]]:
    pairs = [(gene, tx) for gene in GENES for tx in gene.transcripts]
    return sorted(pairs, key=lambda pair: pair[1].number)


def clip(tx: Transcript, lo: int, hi: int) -> list[Part]:
    parts = [(k, (max(s, lo), min(e, hi))) for k, (s, e) in enumerate(tx.exons, 1)]
    return [(k, (s, e)) for k, (s, e) in parts if s <= e]


def utr3(gene: Gene, tx: Transcript) -> list[Part]:
    if tx.cds is None:
        return []
    start, end = tx.cds
    return clip(tx, end + 1, sys.maxsize) if gene.strand == 1 else clip(tx, 1, start - 1)


def utr5(gene: Gene, tx: Transcript) -> list[Part]:
    if tx.cds is None:
        return []
    start, end = tx.cds
    return clip(tx, 1, start - 1) if gene.strand == 1 else clip(tx, end + 1, sys.maxsize)


def cds_parts(tx: Transcript) -> list[Part]:
    return [] if tx.cds is None else clip(tx, *tx.cds)


def cds_bases(gene: Gene, tx: Transcript) -> list[tuple[int, int]]:
    """The exon number and the position of each CDS base, in transcript order."""
    bases = []
    for k, (s, e) in cds_parts(tx):
        positions = range(s, e + 1) if gene.strand == 1 else range(e, s - 1, -1)
        bases.extend((k, p) for p in positions)
    return bases


def join_bases(bases: list[tuple[int, int]]) -> list[Part]:
    exons = sorted({k for k, _ in bases})
    return [
        (k, (min(p for j, p in bases if j == k), max(p for j, p in bases if j == k))) for k in exons
    ]


def start_codon(gene: Gene, tx: Transcript) -> list[Part]:
    return join_bases(cds_bases(gene, tx)[:CODON])


def stop_codon(gene: Gene, tx: Transcript) -> list[Part]:
    return join_bases(cds_bases(gene, tx)[-CODON:])


def cds_without_stop(gene: Gene, tx: Transcript) -> list[Part]:
    return join_bases(cds_bases(gene, tx)[:-CODON])


def width(parts: list[Part]) -> int:
    return sum(e - s + 1 for _, (s, e) in parts)


def inside(inner: Span, outer: Span) -> bool:
    return outer[0] <= inner[0] and inner[1] <= outer[1]


def overlaps(a: Span, b: Span) -> bool:
    return a[0] <= b[1] and b[0] <= a[1]


def gap(a: Span, b: Span) -> int:
    """The number of bases between two ranges that do not overlap."""
    return max(a[0], b[0]) - min(a[1], b[1]) - 1


def structure_checks() -> list[tuple[str, bool]]:
    checks = []
    for gene, tx in transcripts_by_id():
        exons = tx.exons
        ordered = all(
            (a[1] < b[0]) if gene.strand == 1 else (b[1] < a[0]) for a, b in pairwise(exons)
        )
        checks += [
            (f"{tx.tx_id}: the exons are in transcript order and do not overlap", ordered),
            (
                f"{tx.tx_id}: each exon is inside contig {gene.contig}",
                all(1 <= s <= e <= CONTIGS[gene.contig] for s, e in exons),
            ),
        ]
        if tx.cds is not None:
            checks += [
                (
                    f"{tx.tx_id}: the CDS starts and ends inside an exon",
                    all(any(inside((p, p), x) for x in exons) for p in tx.cds),
                ),
                (
                    f"{tx.tx_id}: the CDS length is a multiple of 3",
                    width(cds_parts(tx)) % CODON == 0,
                ),
            ]
    for gene in GENES:
        canonical = [tx for tx in gene.transcripts if tx.canonical]
        checks.append((f"{gene.gene_id} has one canonical transcript", len(canonical) == 1))
    numbers = sorted(tx.number for g in GENES for tx in g.transcripts)
    checks += [
        ("TX_ROW_ORDER holds each transcript once", sorted(TX_ROW_ORDER) == numbers),
        (
            "the transcripts of each gene are in the order of TX_ROW_ORDER",
            all(
                [tx.number for tx in g.transcripts]
                == [n for n in TX_ROW_ORDER if n in {tx.number for tx in g.transcripts}]
                for g in GENES
            ),
        ),
        ("bases 1 to 100 of contig 1 are N", (1, 100) in N_BLOCKS["1"]),
    ]
    return checks


def exon_numbers(parts: list[Part]) -> list[int]:
    return [k for k, _ in parts]


def case_checks() -> list[tuple[str, bool]]:
    g01, g02, g07 = gene_of(1), gene_of(3), gene_of(10)
    t01, t02, t03, t04, t05, t10 = (tx_of(n) for n in (1, 2, 3, 4, 5, 10))
    t01_utr3, t02_utr3 = utr3(g01, t01), utr3(g01, t02)
    t03_utr3 = [span for _, span in utr3(g02, t03)]
    t04_utr3 = [span for _, span in utr3(g02, t04)]
    t04_gaps = [
        min((gap(a, b) for b in t03_utr3 if not overlaps(a, b)), default=None) for a in t04_utr3
    ]
    t01_stop = t01.cds[1] if t01.cds else 0
    t02_stop = t02.cds[1] if t02.cds else 0
    rows = [tx_of(n).tx_id for n in TX_ROW_ORDER]
    utr3_of = {n: width(utr3(gene_of(n), tx_of(n))) for n in TX_ROW_ORDER}
    utr5_of = {n: width(utr5(gene_of(n), tx_of(n))) for n in TX_ROW_ORDER}
    return [
        ("T01: the CDS ends inside exon 2", t01.exons[1][0] <= t01_stop < t01.exons[1][1]),
        ("T01: the 3' UTR is across exons 2 and 3", exon_numbers(t01_utr3) == [2, 3]),
        ("T02: the same tx_seq_start as T01", t02.span[0] == t01.span[0]),
        ("T02: the CDS ends at the last base of exon 2", t02_stop == t02.exons[1][1]),
        ("T02: the 3' UTR length of T01", width(t02_utr3) == width(t01_utr3)),
        ("T02: a lower tx row than T01", rows.index(t02.tx_id) < rows.index(t01.tx_id)),
        ("T03: 4 exons", exon_numbers(clip(t03, 1, sys.maxsize)) == [1, 2, 3, 4]),
        ("T03: the 5' UTR is across 2 exons", exon_numbers(utr5(g02, t03)) == [1, 2]),
        ("T03: the 3' UTR is across 2 exons", exon_numbers(utr3(g02, t03)) == [3, 4]),
        ("T03: the stop codon is across an exon boundary", len(stop_codon(g02, t03)) > 1),
        (
            "T03: an N block is across an exon boundary",
            any(
                any(overlaps(block, x) for x in t03.exons)
                and not any(inside(block, x) for x in t03.exons)
                for block in N_BLOCKS.get(g02.contig, ())
            ),
        ),
        (
            "T04: a 3' UTR range overlaps a T03 range",
            any(overlaps(a, b) for a in t04_utr3 for b in t03_utr3),
        ),
        ("T04: a 3' UTR range touches a T03 range (gap 0)", 0 in t04_gaps),
        ("T04: a 3' UTR range is 1 base from a T03 range (gap 1)", 1 in t04_gaps),
        (
            "T05: a lncRNA with 2 exons and no CDS",
            gene_of(5).biotype == t05.biotype == "lncRNA"
            and exon_numbers(clip(t05, 1, sys.maxsize)) == [1, 2]
            and t05.cds is None,
        ),
        ("T06: a 3' UTR on SCAF1.1", gene_of(6).contig == "SCAF1.1" and utr3_of[6] > 0),
        ("T07: the CDS ends at the transcript end", utr3_of[7] == 0),
        ("T08: a 3' UTR", utr3_of[8] > 0),
        ("T09: a 3' UTR of 5 bases", utr3_of[9] == T09_UTR3_BASES),
        ("T10: the same tx_seq_start as T01", t10.span[0] == t01.span[0]),
        (
            "T10: the 3' UTR holds an N block",
            any(
                inside(block, span)
                for block in N_BLOCKS.get(g07.contig, ())
                for _, span in utr3(g07, t10)
            ),
        ),
        ("T10: a higher tx row than T01", rows.index(t10.tx_id) > rows.index(t01.tx_id)),
        ("G07: a lower gene row than G01", GENES.index(g07) < GENES.index(g01)),
        ("T11: nonsense_mediated_decay", tx_of(11).biotype == "nonsense_mediated_decay"),
        (
            "TSL: T01 1, T07 1, T08 2, T10 2, each other transcript NULL",
            {n: tx_of(n).support_level for n in TX_ROW_ORDER}
            == {**dict.fromkeys(TX_ROW_ORDER), 1: 1, 7: 1, 8: 2, 10: 2},
        ),
        (
            "canonical: T01, T03, T05, T06, T07, T09, T10, T12",
            {n for n in TX_ROW_ORDER if tx_of(n).canonical} == {1, 3, 5, 6, 7, 9, 10, 12},
        ),
        (
            "gene names: Syn<n>, and no name for G06",
            {g.number: g.name for g in GENES} == {n: f"Syn{n}" for n in range(1, 9)} | {6: ""},
        ),
        (
            "Entrez IDs: G01 1001, G02 1002, G05 1005, G06 1006 and 1007",
            {g.number: g.entrez for g in GENES if g.entrez}
            == {1: (1001,), 2: (1002,), 5: (1005,), 6: (1006, 1007)},
        ),
        (
            "gene contigs, strands, and biotypes",
            {g.number: (g.contig, g.strand, g.biotype) for g in GENES}
            == dict.fromkeys(range(1, 9), ("1", 1, PC))
            | {2: ("1", -1, PC), 3: ("1", 1, "lncRNA"), 4: ("SCAF1.1", 1, PC)}
            | {5: ("X", 1, PC), 7: ("1", -1, PC)},
        ),
        ("T12: a 5' UTR and no 3' UTR", utr5_of[12] > 0 and utr3_of[12] == 0),
        (
            "each 3' UTR of G01 and G02 holds a soft-masked run",
            all(
                any(
                    inside(run, span)
                    for run in SOFT_MASKED.get(g.contig, ())
                    for _, span in utr3(g, tx)
                )
                for g in (g01, g02)
                for tx in g.transcripts
            ),
        ),
    ]


def contig_sequences() -> dict[str, NDArray[np.uint8]]:
    """The uppercase letters of each contig, with N in each N block."""
    rng = np.random.default_rng(SEED)
    letters = np.frombuffer(b"ACGT", dtype=np.uint8)
    genome = {}
    for name, length in CONTIGS.items():
        seq = letters[rng.integers(0, len(letters), size=length)]
        for s, e in N_BLOCKS.get(name, ()):
            seq[s - 1 : e] = ord("N")
        genome[name] = seq
    return genome


def blocks(spans: tuple[Span, ...]) -> bytes:
    starts = [s - 1 for s, _ in spans]
    sizes = [e - s + 1 for s, e in spans]
    return struct.pack(f"<I{len(spans)}I{len(spans)}I", len(spans), *starts, *sizes)


def twobit_record(name: str, seq: NDArray[np.uint8]) -> bytes:
    # The packed DNA holds T (code 0) at each N position, as faToTwoBit writes it.
    codes = np.zeros(256, dtype=np.uint8)
    for letter, code in TWOBIT_CODES.items():
        codes[ord(letter)] = code
    padded = np.zeros(-(-len(seq) // 4) * 4, dtype=np.uint8)
    padded[: len(seq)] = codes[seq]
    quads = padded.reshape(-1, 4)
    packed = (quads[:, 0] << 6) | (quads[:, 1] << 4) | (quads[:, 2] << 2) | quads[:, 3]
    header = struct.pack("<I", len(seq)) + blocks(N_BLOCKS.get(name, ()))
    header += blocks(SOFT_MASKED.get(name, ())) + struct.pack("<I", 0)
    return header + packed.astype(np.uint8).tobytes()


def twobit(genome: dict[str, NDArray[np.uint8]], version: int) -> bytes:
    """The UCSC 2bit file, little-endian. Version 1 has 64-bit offsets in the index."""
    offset_format = "<I" if version == 0 else "<Q"
    records = [twobit_record(name, seq) for name, seq in genome.items()]
    index_size = sum(1 + len(name) + struct.calcsize(offset_format) for name in genome)
    offset = struct.calcsize("<4I") + index_size
    index = b""
    for name, record in zip(genome, records, strict=True):
        index += bytes([len(name)]) + name.encode() + struct.pack(offset_format, offset)
        offset += len(record)
    return struct.pack("<4I", TWOBIT_MAGIC, version, len(genome), 0) + index + b"".join(records)


def fasta(genome: dict[str, NDArray[np.uint8]]) -> bytes:
    out = []
    for name, seq in genome.items():
        letters = seq.copy()
        for s, e in SOFT_MASKED.get(name, ()):
            letters[s - 1 : e] = np.frombuffer(letters[s - 1 : e].tobytes().lower(), np.uint8)
        text = letters.tobytes()
        lines = [text[i : i + LINE_BASES] for i in range(0, len(text), LINE_BASES)]
        out.append(b">" + name.encode() + b"\n" + b"".join(line + b"\n" for line in lines))
    return b"".join(out)


def run_htslib(bgzip: str, samtools: str) -> None:
    # The executables come from shutil.which and the arguments are fixed names.
    with (OUT / "synthetic.fa.gz").open("wb") as out:
        subprocess.run(  # noqa: S603
            [bgzip, "--index", "--index-name", "synthetic.fa.gz.gzi", "--stdout", "synthetic.fa"],
            cwd=OUT,
            stdout=out,
            check=True,
        )
    for name in ("synthetic.fa.gz", "synthetic.fa"):
        subprocess.run([samtools, "faidx", name], cwd=OUT, check=True)  # noqa: S603


def block_boundary_check() -> tuple[str, bool]:
    """Read the first BGZF block end from the .gzi and the base offsets from the .fai."""
    gzi = (OUT / "synthetic.fa.gz.gzi").read_bytes()
    (count,) = struct.unpack_from("<Q", gzi)
    boundary = struct.unpack_from("<QQ", gzi, 8)[1] if count else 0
    rows = [line.split("\t") for line in (OUT / "synthetic.fa.gz.fai").read_text().splitlines()]
    fai = {fields[0]: [int(x) for x in fields[1:]] for fields in rows}
    gene = gene_of(10)
    _, offset, line_bases, line_width = fai[gene.contig][:4]

    def byte(position: int) -> int:
        i = position - 1
        return offset + (i // line_bases) * line_width + i % line_bases

    crosses = any(byte(s) < boundary <= byte(e) for _, (s, e) in utr3(gene, tx_of(10)))
    return f"T10: the 3' UTR is across the first BGZF block end (byte {boundary})", crosses


def require(checks: list[tuple[str, bool]]) -> None:
    failures = [message for message, ok in checks if not ok]
    if failures:
        sys.exit("The fixture breaks these cases:\n" + "\n".join(failures))


def exon_ids() -> dict[tuple[str, int, Span], str]:
    """Transcripts of a gene share the ID of an exon with the same coordinates, as in Ensembl."""
    ids: dict[tuple[str, int, Span], str] = {}
    for gene, tx in transcripts_by_id():
        for exon in tx.exons:
            ids.setdefault((gene.contig, gene.strand, exon), f"ENSSYNE{len(ids) + 1:011d}")
    return ids


def write_ensdb(path: Path) -> None:
    ids = exon_ids()
    gene_rows = []
    for g in GENES:
        canonical = next(tx.tx_id for tx in g.transcripts if tx.canonical)
        system = "scaffold" if g.contig in SCAFFOLDS else "chromosome"
        description, version = f"Synthetic gene {g.number}", f"{g.gene_id}.1"
        gene_rows.append(
            (
                g.gene_id,
                g.name,
                g.biotype,
                *g.span,
                g.contig,
                g.strand,
                system,
                description,
                version,
                canonical,
            )
        )
    tx_rows = []
    for number in TX_ROW_ORDER:
        g, tx = gene_of(number), tx_of(number)
        cds = tx.cds or (None, None)
        # tx_id_version, gc_content, and tx_external_name are NULL.
        tx_rows.append(
            (
                tx.tx_id,
                tx.biotype,
                *tx.span,
                *cds,
                g.gene_id,
                tx.support_level,
                None,
                None,
                None,
                int(tx.canonical),
            )
        )
    pairs = transcripts_by_id()
    exon_rows = {ids[g.contig, g.strand, x]: x for g, tx in pairs for x in tx.exons}
    tx2exon_rows = [
        (tx.tx_id, ids[g.contig, g.strand, x], k)
        for g, tx in pairs
        for k, x in enumerate(tx.exons, 1)
    ]
    entrez_rows = [(g.gene_id, e) for g in sorted(GENES, key=lambda g: g.number) for e in g.entrez]
    protein_rows = [
        (tx.tx_id, tx.protein_id, None)
        for number in TX_ROW_ORDER
        if (tx := tx_of(number)).cds is not None
    ]
    with sqlite3.connect(path) as db:
        db.executescript(SCHEMA)
        db.executemany("INSERT INTO metadata VALUES (?, ?)", METADATA)
        db.executemany("INSERT INTO chromosome VALUES (?, ?, 0)", CONTIGS.items())
        db.executemany("INSERT INTO gene VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", gene_rows)
        db.executemany("INSERT INTO entrezgene VALUES (?, ?)", entrez_rows)
        db.executemany("INSERT INTO tx VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", tx_rows)
        db.executemany("INSERT INTO exon VALUES (?, ?, ?)", [(i, *x) for i, x in exon_rows.items()])
        db.executemany("INSERT INTO tx2exon VALUES (?, ?, ?)", tx2exon_rows)
        db.executemany("INSERT INTO protein VALUES (?, ?, ?)", protein_rows)
        db.executescript(INDEXES)
    db.close()


@dataclass(frozen=True)
class GtfStyle:
    gencode: bool
    source: str
    utr5: str
    utr3: str

    def contig(self, name: str) -> str:
        return GENCODE_CONTIGS.get(name, name) if self.gencode else name


ENSEMBL = GtfStyle(gencode=False, source="ensembl", utr5="five_prime_utr", utr3="three_prime_utr")
GENCODE = GtfStyle(gencode=True, source="ENSEMBL", utr5="UTR", utr3="UTR")


def attributes(
    style: GtfStyle, gene: Gene, tx: Transcript | None, exon: tuple[int, str] | None, feature: str
) -> str:
    """Write the attributes in the order of the Ensembl or GENCODE files."""
    q = '"{}"'.format
    fields = [("gene_id", q(gene.gene_id))]
    if tx is not None:
        fields.append(("transcript_id", q(tx.tx_id)))
    name = [("gene_name", q(gene.name))] if gene.name else []
    tsl = "NA" if tx is None or tx.support_level is None else str(tx.support_level)
    protein = [("protein_id", q(tx.protein_id))] if tx is not None and feature == "CDS" else []
    if style.gencode:
        fields += [("gene_type", q(gene.biotype)), *name]
        if tx is not None:
            fields.append(("transcript_type", q(tx.biotype)))
        if exon is not None:
            fields += [("exon_number", str(exon[0])), ("exon_id", q(exon[1]))]
        fields += [("level", "2"), *protein]
        if tx is not None:
            fields += [("transcript_support_level", q(tsl)), ("tag", q("basic"))]
            fields += [("tag", q("Ensembl_canonical"))] if tx.canonical else []
    else:
        if exon is not None and feature not in {style.utr5, style.utr3}:
            fields.append(("exon_number", q(exon[0])))
        fields += [*name, ("gene_biotype", q(gene.biotype))]
        if tx is not None:
            fields += [("transcript_biotype", q(tx.biotype)), *protein]
            fields += [("tag", q("Ensembl_canonical"))] if tx.canonical else []
            if tsl != "NA":
                tsl = f"{tsl} (assigned to previous version 2)"
            fields.append(("transcript_support_level", q(tsl)))
        if feature == "exon" and exon is not None:
            fields.append(("exon_id", q(exon[1])))
    return " ".join(f"{key} {value};" for key, value in fields)


def frames(parts: list[Part]) -> Iterator[tuple[Part, int]]:
    """Give the GTF frame of each part: the bases before the first codon start in the part."""
    done = 0
    for part in parts:
        yield part, (CODON - done % CODON) % CODON
        done += part[1][1] - part[1][0] + 1


CODING_FEATURES: tuple[tuple[str, Callable[[Gene, Transcript], list[Part]]], ...] = (
    ("CDS", cds_without_stop),
    ("start_codon", start_codon),
    ("stop_codon", stop_codon),
)


def transcript_features(
    style: GtfStyle, gene: Gene, tx: Transcript
) -> Iterator[tuple[str, Span, str, tuple[int, str]]]:
    """Give the lines after the transcript line in the order of Ensembl: the exons with their
    CDS and codon parts, then the UTRs."""
    ids = [exon_ids()[gene.contig, gene.strand, exon] for exon in tx.exons]
    pieces: dict[int, list[tuple[str, Span, str]]] = {}
    for feature, parts_of in CODING_FEATURES:
        for (k, span), frame in frames(parts_of(gene, tx)):
            pieces.setdefault(k, []).append((feature, span, str(frame)))
    for k, exon in enumerate(tx.exons, 1):
        yield "exon", exon, ".", (k, ids[k - 1])
        for feature, span, frame in pieces.get(k, []):
            yield feature, span, frame, (k, ids[k - 1])
    for feature, parts in ((style.utr5, utr5(gene, tx)), (style.utr3, utr3(gene, tx))):
        for k, span in parts:
            yield feature, span, ".", (k, ids[k - 1])


def gtf_text(style: GtfStyle) -> str:
    if style.gencode:
        lines = ["##description: synthetic annotation of seedmatch, GENCODE style", "##format: gtf"]
    else:
        lines = [f"#!genome-build {GENOME_BUILD}", f"#!genome-version {GENOME_BUILD}"]

    def line(gene: Gene, feature: str, span: Span, frame: str, attrs: str) -> None:
        strand = "+" if gene.strand == 1 else "-"
        fields = [style.contig(gene.contig), style.source, feature, *map(str, span), ".", strand]
        lines.append("\t".join([*fields, frame, attrs]))

    for gene in GENES:
        line(gene, "gene", gene.span, ".", attributes(style, gene, None, None, "gene"))
        for tx in gene.transcripts:
            line(gene, "transcript", tx.span, ".", attributes(style, gene, tx, None, "transcript"))
            for feature, span, frame, exon in transcript_features(style, gene, tx):
                line(gene, feature, span, frame, attributes(style, gene, tx, exon, feature))
    return "\n".join(lines) + "\n"


def main() -> None:
    tools = {name: shutil.which(name) for name in ("bgzip", "samtools")}
    missing = [name for name, path in tools.items() if path is None]
    if missing:
        sys.exit(f"{' and '.join(missing)} not found on PATH; install htslib and samtools")
    require(structure_checks() + case_checks())
    OUT.mkdir(parents=True, exist_ok=True)
    for name in OUTPUTS:
        (OUT / name).unlink(missing_ok=True)
    genome = contig_sequences()
    (OUT / "synthetic.2bit").write_bytes(twobit(genome, version=0))
    (OUT / "synthetic.v1.2bit").write_bytes(twobit(genome, version=1))
    (OUT / "synthetic.fa").write_bytes(fasta(genome))
    run_htslib(str(tools["bgzip"]), str(tools["samtools"]))
    require([block_boundary_check()])
    write_ensdb(OUT / "synthetic.sqlite")
    (OUT / "synthetic.gtf").write_text(gtf_text(ENSEMBL))
    (OUT / "synthetic.gencode.gtf.gz").write_bytes(
        gzip.compress(gtf_text(GENCODE).encode(), mtime=0)
    )
    print(f"wrote {len(OUTPUTS)} files to {OUT}")


if __name__ == "__main__":
    main()
