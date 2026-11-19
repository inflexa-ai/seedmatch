# Annotation: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port the SeedMatchR annotation path (`load_annotations` and its helpers) to `seedmatch`: the EnsDb and 2bit path with output identical to R, and the GTF and FASTA path of the inflexa catalog.

**Architecture:** `seedmatch.annotation` reads an EnsDb with `sqlite3` and puts each filter in the SQL `WHERE` clause. The GTF reader writes the same tables into an in-memory SQLite database. Thus the two paths share the filters, the feature derivation, and the order rules. The genome readers (2bit, plain FASTA, BGZF FASTA) decode only the requested ranges into the 4-bit masks of `SequenceSet`. `Reference` joins one annotation and one genome, and gives transcripts, feature ranges, and sequences.

**Tech Stack:** Python 3.12, numpy, pandas, the standard library (`sqlite3`, `mmap`, `zlib`, `struct`). R 4.5.3 with ensembldb 2.34.0, GenomicFeatures 1.62.0, rtracklayer 1.70.1, Rsamtools 2.26.0, and SeedMatchR 2.0.0 for the goldens. `bgzip` of htslib 1.24 and `samtools` 1.24 for the fixture only.

**Spec:** `docs/design.md`, sections "Annotation" and "Tests > Annotation".

**Later plans (not in this plan):**

- statistics and the ECDF test
- resources and fetch (the registry gets the rnor7 sizes and digests of Task 8)
- the `hg38` parity digests
- `docs/performance.md`, the other analyses, and the plots

## Global Constraints

- No new runtime dependency. The genome readers use only the standard library and numpy. There is no `fasta` extra.
- The filters use the SQLite JSON function `json_each`, thus SQLite 3.38 or later is necessary.
- No core module imports a network module or `seedmatch.fetch` (ruff `TID251`). No function writes a file unless it has an explicit output path parameter. No reader writes next to its input file.
- Coordinates are 1-based and inclusive, in the API and in the goldens.
- Errors derive from `seedmatch.errors.SeedMatchError`. Invalid input raises `InputError`, and each message names the bad value and the expected value. A missing file raises `FileNotFoundError`.
- Each change of behavior from R gets an entry `## D-NNN: <title>` in `docs/divergences.md` and a test whose name contains the ID. `D-002` and `D-009` are the IDs of their rows in the design table. `D-027`, `D-028`, `D-029`, and `D-030` are the new rows of this plan.
- A rename of a column or a change of a value type is not a divergence. `docs/r-map.md` records it.
- The R source of truth is `../SeedMatchR` at commit `2fc0c4f`. The R scripts run with `tools/golden/run.sh`. The tests never run R, `bgzip`, or `samtools`, and they never use the network.
- The rnor7 files are in `~/Development/refs/seedmatch/rnor7/`:
  - `EnsDb.Rnorvegicus.v109.sqlite`: 179650560 bytes, SHA-256 `2453c59e1054157798df5a4114650911d046b5f7795af91a6a2171614ba0289b`
  - `Rattus_norvegicus.mRatBN7.2.dna_sm.toplevel.2bit`: 661989722 bytes, SHA-256 `73eacccd66a8843cced837d0cb0a272b42b7214b4a7c62a3c48963509ea914f3`
  - `Rattus_norvegicus.mRatBN7.2.109.gtf.gz`: 18029537 bytes, SHA-256 `ff7d5ecd69ed7df5e641ec90383e5de9e10e9d82316a0570c71fa80a98e8edca` (Task 8 Step 1)
- Documents use ASD-STE100. Code comments obey the admission test of the project rules. The word "seam" does not occur.
- Git: an implementer does no `git add` and no commit. Each task ends with a checkpoint. The controller commits only with the consent of the user.

## Review Focus

1. An EnsDb and a genome in a read-only directory, as `/mnt/refs` in the sandbox: each reader opens the file and reads it, and it makes no file next to it. Owner: Task 4 (EnsDb) and Task 2 (2bit).
2. A genome of a different build from the EnsDb, for example a rat 2bit with an EnsDb of the same contig names but other lengths: `InputError` when `Reference` opens, with the contig and the two lengths. Owner: Task 4.
3. A GENCODE GTF (names `chr1`) with an Ensembl FASTA (names `1`): the names agree and the sequences are correct. A GTF whose contigs agree with no genome contig: `InputError`. Owner: Task 7.
4. A set filter with a bare `str`, for example `TranscriptFilter(gene_ids="ENSG...")`: `InputError`, not a filter on the letters of the ID. Owner: Task 4.
5. A filter with more IDs than the SQLite parameter limit (32766), for example 40000 gene IDs from a DE table: the query runs and keeps the correct transcripts. Owner: Task 4.

---

### Task 1: Synthetic fixture and R goldens

**Files:**
- Make: `tools/fixtures/annotation.py`
- Make: `tests/data/annotation/` (written by the script): `synthetic.sqlite`, `synthetic.2bit`, `synthetic.v1.2bit`, `synthetic.fa`, `synthetic.fa.fai`, `synthetic.fa.gz`, `synthetic.fa.gz.fai`, `synthetic.fa.gz.gzi`, `synthetic.gtf`, `synthetic.gencode.gtf.gz`
- Make: `tools/golden/R/annotation.R`, `tools/golden/patches/get_longest_utr.R`, `tools/golden/patches/reduce_features.R`
- Make: `tests/golden/annotation/*.tsv` (written by the R script)
- Modify: `tools/golden/R/check_env.R`, `tools/golden/README.md`, `tests/data/README.md`

**Interfaces:**
- Produces: the fixture files above, and the goldens below. Later tasks read them through `tests/annotation_cases.py` (Task 4).
- Produces: the IDs `ENSSYNG00000000001` to `ENSSYNG00000000008` (G01 to G08) and `ENSSYNT00000000001` to `ENSSYNT00000000012` (T01 to T12).

The script writes each file from one specification. A fixed `numpy.random.default_rng` seed makes the contig sequences. Run it with `uv run python tools/fixtures/annotation.py`. It needs `bgzip` and `samtools` on `PATH`, and it stops with a message if they are absent.

Contigs: `1` (150002 bases), `X` (5003 bases), `SCAF1.1` (2001 bases). No length is a multiple of 4, thus the last byte of each 2bit record is partial. Bases 1 to 100 of contig `1` are `N`. Each 3' UTR of G01 and G02 has a soft-masked (lowercase) run.

| Gene | Contig, strand | Gene biotype, name | Transcripts and cases |
|---|---|---|---|
| G01 | `1`, + | protein_coding, `Syn1` | T01 (canonical, TSL 1): 3 exons. The CDS ends inside exon 2, thus the 3' UTR is across exons 2 and 3. T02 (TSL NULL): the same `tx_seq_start` as T01. Its CDS ends at the last base of exon 2. Its 3' UTR length is equal to the length of T01. T02 has a lower row in `tx` than T01. |
| G02 | `1`, − | protein_coding, `Syn2` | T03 (canonical): 4 exons. The 5' UTR and the 3' UTR are each across 2 exons. The stop codon is across an exon boundary. An N block is across an exon boundary. T04: one 3' UTR range overlaps a range of T03, one range touches a range of T03 (gap 0), and one range is 1 base from a range of T03 (gap 1). |
| G03 | `1`, + | lncRNA, `Syn3` | T05 (canonical, lncRNA): 2 exons, no CDS. |
| G04 | `SCAF1.1`, + | protein_coding, `Syn4` | T06 (canonical): a 3' UTR on a non-standard contig. |
| G05 | `X`, + | protein_coding, `Syn5` | T07 (canonical, TSL 1): the CDS ends at the transcript end, thus it has no 3' UTR. T08 (TSL 2): a 3' UTR. |
| G06 | `1`, + | protein_coding, no name | T09 (canonical): a 3' UTR of 5 bases. Entrez IDs 1006 and 1007. |
| G07 | `1`, − | protein_coding, `Syn7` | T10 (canonical, TSL 2): the same `tx_seq_start` as T01. Its 3' UTR is across the end of the first BGZF block of `synthetic.fa.gz`, and it holds an N block. T11 (nonsense_mediated_decay). The `gene` row of G07 comes before the row of G01, but T10 comes after T01 in `tx`. |
| G08 | `1`, + | protein_coding, `Syn8` | T12 (canonical): a 5' UTR and no 3' UTR. |

Entrez IDs: G01 1001, G02 1002, G05 1005, G06 1006 and 1007.

- [ ] **Step 1: Write the EnsDb.** Use the schema 2.2 tables and indexes of the rnor7 EnsDb (`sqlite3 ~/Development/refs/seedmatch/rnor7/EnsDb.Rnorvegicus.v109.sqlite .schema`): `metadata`, `chromosome`, `gene`, `entrezgene`, `tx`, `exon`, and `tx2exon`, with the indexes of these tables. Omit `protein`, `uniprot`, and `protein_domain`. The `tx` CDS columns include the stop codon. The `metadata` table has the names of the rnor7 file, with `Organism` = `Synthetic organism`, `genome_build` = `SYN1`, and `DBSCHEMAVERSION` = `2.2`. Insert the rows in the order of the table above, except for the G07 and T02 rows that the table moves.

- [ ] **Step 2: Write the genome files.** `synthetic.2bit` is UCSC 2bit version 0, little-endian, with N blocks and mask blocks. `synthetic.v1.2bit` is version 1 (64-bit offsets) with the same records. `synthetic.fa` has 60 bases on each line and the soft-masked letters in lowercase. Then run `bgzip --index --index-name synthetic.fa.gz.gzi --stdout synthetic.fa > synthetic.fa.gz`, `samtools faidx synthetic.fa.gz`, and `samtools faidx synthetic.fa`.

- [ ] **Step 3: Write the GTF files.** `synthetic.gtf` has the Ensembl style: `gene`, `transcript`, `exon`, `CDS` (without the stop codon), `start_codon`, `stop_codon`, `five_prime_utr`, and `three_prime_utr` lines. The attributes are `gene_id`, `transcript_id`, `gene_name` (absent for G06), `gene_biotype`, `transcript_biotype`, `tag "Ensembl_canonical"`, `transcript_support_level "1 (assigned to previous version 2)"` or `"NA"`, `exon_number`, and `exon_id`. The genes and transcripts come in the row order of the EnsDb. `synthetic.gencode.gtf.gz` is the same annotation in the GENCODE style: contig names `chr1` and `chrX` (`SCAF1.1` stays), `gene_type` and `transcript_type`, `tag "basic"; tag "Ensembl_canonical";`, `transcript_support_level "1"`, and the unquoted `level 2;`. Write it with `gzip` and `mtime=0`.

- [ ] **Step 4: Write the patches.** `get_longest_utr.R` defines `.get_longest_utr` with one change: before `group_by`, it keeps only the rows of `tx_lens` whose `tx_id` is in `names(utr.grlist)`. `reduce_features.R` defines `.reduce_features` with one change: after `reduce`, it reverses the ranges of each element on the minus strand (`S4Vectors::revElements`).

- [ ] **Step 5: Write `annotation.R`.** It loads the fixture with `ensembldb::EnsDb` and `rtracklayer::TwoBitFile`. It puts a replacement of `load_species_anno_db` into the SeedMatchR namespace with `assignInNamespace`. The replacement returns `list(gtf = ensembldb::transcripts(txdb, filter = filter), dna = dna, txdb = txdb)`. For each config below, it calls `load_annotations("rnor7", ...)` with `return_gene_name = FALSE` and with `TRUE`. A patched config puts the patch into the namespace and restores the original after the config. Each golden has the notes `ensembldb: <version>`, `GenomicFeatures: <version>`, `rtracklayer: <version>`, and `Rsamtools: <version>`.

| Config | Feature | R arguments | Notes |
|---|---|---|---|
| `utr3` | `3UTR` | the defaults | |
| `utr5` | `5UTR` | the defaults | |
| `exons` | `exons` | the defaults | |
| `cds` | `cds` | the defaults | |
| `utr3_all` | `3UTR` | `canonical = FALSE` | |
| `utr3_min0` | `3UTR` | `canonical = FALSE, min.feature.width = 0` | |
| `exons_all` | `exons` | `canonical = FALSE` | |
| `exons_everything` | `exons` | `canonical = FALSE, protein.coding = FALSE, standard.chroms = FALSE` | |
| `utr3_longest` | `3UTR` | `canonical = FALSE, longest.utr = TRUE` | |
| `utr3_longest_noncoding` | `3UTR` | `canonical = FALSE, protein.coding = FALSE, longest.utr = TRUE` | patched `.get_longest_utr`. A note gives the error of the unpatched function. |
| `utr3_union_r` | `3UTR` | `canonical = FALSE, reduce.features = TRUE` | unpatched, R only |
| `utr3_union` | `3UTR` | as `utr3_union_r` | patched `.reduce_features` |
| `utr3_union_min0` | `3UTR` | as `utr3_union`, `min.feature.width = 0` | patched `.reduce_features` |
| `exons_longest_r` | `exons` | `canonical = FALSE, longest.utr = TRUE` | R only |
| `tsl1` | `3UTR` | `canonical = FALSE, transcript.support = 1` | |
| `gene_ids` | `3UTR` | `canonical = FALSE, gene_id = c(G01, G02)` | |
| `tx_ids` | `3UTR` | `canonical = FALSE, tx_id = c(T02, T04)` | |
| `symbols` | `3UTR` | `canonical = FALSE, symbol = c("Syn1", "Syn5")` | |
| `entrez` | `3UTR` | `canonical = FALSE, entrez_id = c(1006, 1007, 1001)` | |
| `nmd` | `exons` | `canonical = FALSE, protein.coding = FALSE, add.filter = AnnotationFilterList(TxBiotypeFilter("nonsense_mediated_decay"))` | |
| `lnc` | `exons` | `canonical = FALSE, protein.coding = FALSE, add.filter = AnnotationFilterList(GeneBiotypeFilter("lncRNA"))` | |

For each config, the script writes these goldens in `tests/golden/annotation/`:

- `<config>_features.tsv`: `name`, `rank`, `seqnames`, `start`, `end`, `strand`. One row for each range of `$features` with a width above 0, in the order of the list. `rank` counts only these ranges.
- `<config>_sequences.tsv`: `name_tx` (the name with `return_gene_name = FALSE`), `name_gene` (the name with `TRUE`), `width`, `sequence`.
- `transcripts_<config>.tsv`: `tx_id`, `gene_id`, `seqnames`, `start`, `end`, `strand`, `tx_biotype`, `tx_cds_seq_start`, `tx_cds_seq_end`, `tx_support_level`, `tx_is_canonical`, from `as.data.frame(ensembldb::transcripts(txdb, filter = afl))` with the filter of the config, in the R order.

It also writes `genome_ranges.tsv`: `label`, `contig`, `start`, `end`, `twobit` (`getSeq` on `synthetic.2bit`), `twobit_v1` (on `synthetic.v1.2bit`, or `unsupported` if rtracklayer cannot read version 1), and `bgzf` (`Rsamtools::scanFa` on `synthetic.fa.gz`). The ranges are:

- all of contig `X` and of contig `SCAF1.1`
- bases 1 to 120 of contig `1`
- a range that starts at each of the 4 positions in a 2bit byte
- the range of T10 across the block boundary
- a soft-masked run and a range of width 1
- the last 7 bases of each contig

- [ ] **Step 6: Run the script.** Run `tools/golden/run.sh tools/golden/R/annotation.R`. Expected: the goldens exist. In `genome_ranges.tsv`, `twobit` and `bgzf` are equal in each row, and each soft-masked base is uppercase. In `transcripts_exons_everything.tsv`, read the order of T01, T02, and T10 (the tie of the same start). If T10 does not come first, the order of the `tx` rows decides a tie before the order of the `gene` rows. Then change the order rule in `docs/design.md` and in Task 4, and record the ruling.

- [ ] **Step 7: Document.** Add the annotation packages to `check_env.R`. In `tools/golden/README.md`, add the commands of the fixture script and of `annotation.R`, and the purpose of each patch. In `tests/data/README.md`, tell that `annotation/` is synthetic and give the script. Checkpoint.

### Task 2: Path helper and 2bit genome

**Files:**
- Make: `src/seedmatch/_io.py`, `src/seedmatch/annotation/__init__.py`, `src/seedmatch/annotation/genome.py`, `src/seedmatch/annotation/twobit.py`
- Modify: `src/seedmatch/de.py` (use `local_path`), `src/seedmatch/sequence/alphabet.py`
- Test: `tests/annotation/test_twobit.py`, `tests/test_io.py`, `tests/sequence/test_alphabet.py`

**Interfaces:**
- Produces: `seedmatch._io.local_path(path: str | os.PathLike[str]) -> Path`. It raises `InputError` `f"{shown!r}: the path is a URL; expected a local path"` with the `_URL` pattern of `de.py`. `read_deseq2_table` uses it.
- Produces: `seedmatch.sequence.alphabet.reverse_complement_codes(codes: NDArray[np.uint8]) -> NDArray[np.uint8]`: the complement mask of each code, in reverse order.
- Produces: `seedmatch.annotation.genome.Genome` (abstract base class) with `path: Path`, the property `contigs -> Mapping[str, int]` (name to length), `codes(contig: str, start: int, end: int) -> NDArray[np.uint8]`, the abstract `_codes(contig, start, end)` that a subclass implements, `sequence(contig, start, end) -> str` (`decode_dna` of `codes`), `resolve(name: str) -> str | None`, `close()`, `__enter__`, and `__exit__`.
- Produces: `seedmatch.annotation.twobit.TwoBitGenome(path: str | os.PathLike[str])`, a `Genome`.

`codes` takes 1-based inclusive coordinates. It checks the contig and the range, then calls `_codes`. Thus each reader gives the same errors. `end == start - 1` gives an empty array, for `0 <= end <= length`. A duplicate contig name in a genome file raises `InputError`. An unknown contig raises `InputError` that names the contig and the first contigs of the genome. A range outside `1..length` raises `InputError` with the range and the length. The codes do not mark the soft-masked bases, thus `sequence` gives uppercase text.

`resolve` gives the genome name of an annotation contig: the same name, then the name without `chr`, then the name with `chr`. `chrM` and `MT` are the same contig. It gives `None` if no name exists.

- [ ] **Step 1: Write the failing tests**

```python
@pytest.mark.parametrize("name", ["synthetic.2bit", "synthetic.v1.2bit"])
@pytest.mark.parametrize("row", golden_rows("annotation/genome_ranges.tsv"), ids=label)
def test_twobit_ranges_match_r(name, row):
    with TwoBitGenome(ANNOTATION / name) as genome:
        assert genome.sequence(row.contig, int(row.start), int(row.end)) == row.twobit

def test_contigs_give_the_lengths():
    with TwoBitGenome(ANNOTATION / "synthetic.2bit") as genome:
        assert dict(genome.contigs) == {"1": 150002, "X": 5003, "SCAF1.1": 2001}

def test_big_endian_file_gives_the_same_sequences(tmp_path): ...
def test_range_after_the_contig_end_raises():  # match=r"5004.*length 5003"
def test_unknown_contig_raises_and_names_the_contigs():  # match=r"'Y'.*1, SCAF1.1, X"
def test_file_that_is_not_2bit_raises(tmp_path):  # match="not a 2bit file"
def test_unknown_version_raises(tmp_path):  # match=r"version 2; expected 0 or 1"
def test_truncated_file_raises(tmp_path):  # match="ends before"
def test_reads_from_a_read_only_directory(tmp_path): ...  # chmod 0o555 directory, 0o444 file
def test_resolve_maps_chr_names_and_mt(): ...  # "chr1" -> "1", "MT" <-> "chrM", "Z" -> None
```

`test_big_endian_file_gives_the_same_sequences` rewrites each `uint32` of the header, the index, and the record headers of `synthetic.2bit` in big-endian order, with a helper in the test. In `tests/test_io.py`, move the URL test of `read_deseq2_table` to `local_path`, and keep one test that `read_deseq2_table` still rejects a URL. In `test_alphabet.py`, assert that `reverse_complement_codes(encode_dna(s))` equals `encode_dna(reverse_complement(s))` for each IUPAC letter.

- [ ] **Step 2: Run `uv run pytest tests/annotation/test_twobit.py tests/test_io.py -v`.** Expected: FAIL with `ModuleNotFoundError`.

- [ ] **Step 3: Write `_io.py`, `genome.py`, and `twobit.py`.** The reader maps the file with `mmap` (`ACCESS_READ`). At open, it reads the header and the index. It reads the N blocks of a contig at the first request of the contig. The magic `0x1A412743` in the native order or in the swapped order sets the byte order. A range decodes only the bytes of the packed DNA that hold the range, through a 256 × 4 table. The 2-bit codes `T`, `C`, `A`, `G` (0 to 3, the first base in the high bits) give the masks 8, 2, 1, 4. Then the N blocks that overlap the range set the mask 15. The reader ignores the mask blocks.

- [ ] **Step 4: Run the tests.** Expected: PASS. Run `uv run ruff check` and `uv run mypy`. Checkpoint.

### Task 3: FASTA genome, plain and BGZF

**Files:**
- Make: `src/seedmatch/annotation/fasta_genome.py`
- Make: `src/seedmatch/annotation/genome_file.py` (`open_genome`). It imports each reader, and each reader imports `genome.py`. Thus `genome.py` cannot hold `open_genome` without an import cycle.
- Test: `tests/annotation/test_fasta_genome.py`

**Interfaces:**
- Consumes: `Genome` (Task 2).
- Produces: `seedmatch.annotation.fasta_genome.FastaGenome(path: str | os.PathLike[str])`, a `Genome`.
- Produces: `seedmatch.annotation.genome_file.open_genome(path: str | os.PathLike[str]) -> Genome`. The first bytes select the reader: the 2bit magic in each order gives `TwoBitGenome`, and the gzip magic or `>` gives `FastaGenome`. Other bytes raise `InputError`.

The reader reads `<path>.fai` (name, length, offset, line bases, line width). A BGZF file also needs `<path>.gzi`. Base `i` (0-based) is at the uncompressed byte `offset + (i // line_bases) * line_width + i % line_bases`. The reader reads the bytes of the range and removes `\n` and `\r`. For BGZF, the `.gzi` pairs (compressed offset, uncompressed offset) and the implicit pair (0, 0) give the first block of the range (`bisect`). Each block is a gzip member with the extra subfield `BC`. `zlib.decompress(data, wbits=-15)` decompresses its data. The reader keeps the last decompressed block.

- [ ] **Step 1: Write the failing tests**

```python
@pytest.mark.parametrize("name", ["synthetic.fa", "synthetic.fa.gz"])
@pytest.mark.parametrize("row", golden_rows("annotation/genome_ranges.tsv"), ids=label)
def test_fasta_ranges_match_r(name, row):
    with FastaGenome(ANNOTATION / name) as genome:
        assert genome.sequence(row.contig, int(row.start), int(row.end)) == row.bgzf

def test_fasta_and_twobit_give_the_same_contigs(): ...  # each contig, from 1 to its length
def test_open_genome_selects_the_reader_from_the_first_bytes(tmp_path): ...  # names without suffixes
def test_missing_fai_raises_and_names_the_command(tmp_path):  # match=r"synthetic\.fa\.fai.*samtools faidx"
def test_missing_gzi_raises(tmp_path):  # match=r"\.gzi"
def test_gzip_that_is_not_bgzf_raises(tmp_path):  # match=r"not BGZF.*bgzip"
def test_crlf_line_ends_give_the_same_sequence(tmp_path): ...  # the test writes the file and its .fai
def test_fai_with_too_few_columns_raises(tmp_path):  # match=r"line 1.*5 columns"
def test_t10_range_crosses_a_block_boundary(): ...  # the .gzi has a block start inside the bytes of the range
```

- [ ] **Step 2: Run `uv run pytest tests/annotation/test_fasta_genome.py -v`.** Expected: FAIL with `ImportError`.

- [ ] **Step 3: Write `fasta_genome.py` and `open_genome`.**

- [ ] **Step 4: Run the tests.** Expected: PASS. Checkpoint.

### Task 4: Filters, EnsDb reader, and `Reference.transcripts`

**Files:**
- Make: `src/seedmatch/annotation/filters.py`, `src/seedmatch/annotation/ensdb.py`, `src/seedmatch/annotation/reference.py`
- Make: `tests/annotation_cases.py`
- Modify: `pyproject.toml` (`extend-immutable-calls`)
- Test: `tests/annotation/test_filters.py`, `tests/annotation/test_ensdb.py`

**Interfaces:**
- Consumes: `open_genome`, `Genome.contigs`, `Genome.resolve` (Tasks 2 and 3).
- Produces: in `filters.py`:
  - `FeatureType(StrEnum)`: `UTR3 = "3UTR"`, `UTR5 = "5UTR"`, `EXONS = "exons"`, `CDS = "cds"`
  - `PerGene(Enum)`: `EACH_TRANSCRIPT`, `LONGEST`, `UNION`
  - `STANDARD_CHROMOSOMES: Final[frozenset[str]]`: `1` to `22`, `X`, `Y`, and each with the prefix `chr`
  - `TranscriptFilter` (frozen dataclass): `standard_chromosomes: bool = True`, `require_tx_id: bool = True`, `protein_coding: bool = True`, `canonical: bool = True`, `support_level: int | None = None`, `gene_ids`, `tx_ids`, `symbols`, `tx_biotypes`, `gene_biotypes: Sequence[str] | None = None`, `entrez_ids: Sequence[int] | None = None`. `__post_init__` stores each set field as a `tuple`.
- Produces: in `ensdb.py`, `AnnotationDb(connection: sqlite3.Connection, source: str)` with:
  - `AnnotationDb.open(path: str | os.PathLike[str]) -> AnnotationDb` (classmethod)
  - `contig_lengths() -> dict[str, int]` (the `chromosome` table, or `{}` without it)
  - `contig_names() -> list[str]` (the distinct `gene.seq_name` values)
  - `transcripts(transcript_filter: TranscriptFilter) -> pd.DataFrame`
  - `where(transcript_filter: TranscriptFilter) -> tuple[str, list[object]]`: the SQL condition and its parameters
  - `close()`
- Produces: `Reference(annotation: AnnotationDb, genome: Genome)` in `reference.py`, with `Reference.from_paths(*, ensdb, genome) -> Reference`, `transcripts(transcript_filter: TranscriptFilter = TranscriptFilter()) -> pd.DataFrame`, `close()`, `__enter__`, and `__exit__`.
- Produces: `tests/annotation_cases.py` with `ANNOTATION = DATA / "annotation"`, `Case(feature, transcript_filter, per_gene, min_width)` (a `NamedTuple`), `CASES: dict[str, Case]` (each config of Task 1 except the configs that end with `_r`), and `fixture_reference() -> Reference` (`synthetic.sqlite` with `synthetic.2bit`).

`AnnotationDb.open` connects to `f"{path.resolve().as_uri()}?mode=ro&immutable=1"` with `uri=True`. A file that is not SQLite raises `InputError` with the path. A file without one of the tables `gene`, `tx`, `tx2exon`, and `exon` raises `InputError` with the absent table.

`where` gives one SQL condition for each field that is set, joined with `AND`. Each value is a bound parameter. Each set is one JSON parameter: `<column> IN (SELECT value FROM json_each(?))`.

| Field | Condition |
|---|---|
| `standard_chromosomes` | `gene.seq_name` in `STANDARD_CHROMOSOMES` |
| `require_tx_id` | `tx.tx_id != ''` |
| `protein_coding` | `gene.gene_biotype = 'protein_coding' AND tx.tx_biotype = 'protein_coding'` |
| `canonical` | `tx.tx_is_canonical = 1` |
| `support_level` | `tx.tx_support_level = ?` |
| `gene_ids`, `tx_ids`, `symbols` | `tx.gene_id`, `tx.tx_id`, `gene.gene_name` in the set |
| `tx_biotypes`, `gene_biotypes` | `tx.tx_biotype`, `gene.gene_biotype` in the set |
| `entrez_ids` | `tx.gene_id IN (SELECT gene_id FROM entrezgene WHERE entrezid IN (SELECT value FROM json_each(?)))` |

A field that needs an absent column or table (`tx_is_canonical`, `tx_support_level`, `entrezgene`) raises `InputError` with the source and the column. The SQL text holds only fixed names, thus each `S608` ignore gets that reason.

`transcripts` gives the columns `tx_id`, `gene_id`, `gene_name`, `contig`, `start`, `end`, `strand` (`"+"` or `"-"`), `tx_biotype`, `gene_biotype`, `cds_start`, `cds_end`, `support_level` (`Int64`), and `canonical` (`boolean`). An absent optional column gives missing values. The order is `ORDER BY gene.seq_name, tx.tx_seq_start, gene.rowid, tx.rowid`.

`Reference.__init__` makes sure that the annotation and the genome agree. If no name of `contig_names()` resolves in the genome, it raises `InputError` with some names of each. For each contig of `contig_lengths()` that resolves, a different genome length raises `InputError` with the contig and the two lengths.

- [ ] **Step 1: Write the failing tests**

```python
@pytest.mark.parametrize("config", CASES)
def test_transcripts_match_r(config):  # exons_everything: the tie T10, T02, T01 (D-030)
    expected = golden(f"annotation/transcripts_{config}.tsv")
    with fixture_reference() as ref:
        got = ref.transcripts(CASES[config].transcript_filter)
    assert got["tx_id"].tolist() == expected["tx_id"].tolist()
    # also gene_id, contig, start, end, strand, tx_biotype, cds_start, cds_end, support_level, canonical
    # R gives NA for a missing value

def test_set_field_rejects_a_str():  # TranscriptFilter(gene_ids="ENSSYNG00000000001"), match=r"gene_ids.*str"
def test_set_field_rejects_an_empty_sequence():  # match="gene_ids is empty"
@pytest.mark.parametrize("level", [0, 6, True, 1.5, "1"])
def test_support_level_must_be_an_integer_from_1_to_5(level): ...
def test_entrez_ids_reject_text():  # entrez_ids=("1006",)
def test_flag_rejects_a_value_that_is_not_bool():  # canonical="yes"
def test_filter_with_more_ids_than_the_parameter_limit_keeps_the_correct_transcripts(): ...
def test_canonical_on_a_file_without_the_column_raises(tmp_path):  # ALTER TABLE tx DROP COLUMN, match="tx_is_canonical"
def test_entrez_on_a_file_without_the_table_raises(tmp_path):  # match="entrezgene"
def test_file_that_is_not_sqlite_raises(tmp_path): ...
def test_missing_ensdb_raises_file_not_found(tmp_path): ...  # sqlite3 says "unable to open"
def test_sqlite_file_without_the_tx_table_raises(tmp_path):  # match="table 'tx'"
def test_reference_opens_from_a_read_only_directory(tmp_path): ...  # and makes no new file there
def test_genome_with_another_length_raises(tmp_path):  # chromosome X length 4000, match=r"'X'.*4000.*5003"
def test_genome_without_an_agreeing_contig_raises(tmp_path): ...
def test_tie_order_does_not_depend_on_the_filter_D030(): ...
```

`D030`: R gives T10, T02, T01 in `transcripts_exons_all.tsv` (with the standard-chromosome filter), but T02, T01, T10 in `transcripts_exons_everything.tsv` (without it). Python gives T10, T02, T01 with each filter. For `exons_everything`, `test_transcripts_match_r` expects the R rows with this one tie in the Python order. Add the `D-030` entry to `docs/divergences.md`.

The parameter-limit test gives 40000 IDs: the 8 fixture gene IDs and 39992 IDs that do not exist. The result has each transcript of the fixture genes.

- [ ] **Step 2: Run `uv run pytest tests/annotation -v`.** Expected: the new tests FAIL with `ImportError`.

- [ ] **Step 3: Write `filters.py`, `ensdb.py`, and `reference.py`.** The flags use `as_flag`. `support_level` accepts an `int` from 1 to 5, not a `bool`. Add `TranscriptFilter` at its 3 import paths to `extend-immutable-calls`.

- [ ] **Step 4: Run the tests.** Expected: PASS. Checkpoint.

### Task 5: Feature ranges and the per-gene modes

**Files:**
- Make: `src/seedmatch/annotation/features.py`
- Modify: `src/seedmatch/annotation/ensdb.py` (add `exon_rows`), `src/seedmatch/annotation/reference.py` (add `features`)
- Modify: `docs/divergences.md` (`D-027`, `D-028`, `D-029`)
- Test: `tests/annotation/test_features.py`

**Interfaces:**
- Consumes: `AnnotationDb`, `TranscriptFilter`, `FeatureType`, `PerGene`, `CASES` (Task 4).
- Produces: `Features` (frozen dataclass, `eq=False`) in `features.py`: `ids: tuple[str, ...]`, `gene_ids: tuple[str, ...]`, `contigs: tuple[str, ...]`, `strands: NDArray[np.int8]` (1 or −1), `starts` and `ends: NDArray[np.int64]` (one entry for each range), `offsets: NDArray[np.int64]` (feature `i` has the ranges `offsets[i]:offsets[i + 1]`, in transcript order). Methods: `__len__`, `widths() -> NDArray[np.int64]`, `to_frame() -> pd.DataFrame` with the columns `id`, `gene_id`, `contig`, `strand`, `rank`, `start`, `end`.
- Produces: `AnnotationDb.exon_rows(transcript_filter: TranscriptFilter, *, coding_only: bool) -> ExonRows`. `ExonRows` holds one numpy array or tuple for each column: `tx_id`, `gene_id`, `seq_name`, `strand`, `cds_start`, `cds_end`, `exon_start`, `exon_end`, `tx_seq_start`, `tx_row`. The rows are in `ORDER BY tx.tx_id, tx2exon.exon_idx`. `coding_only` adds `tx.tx_cds_seq_start IS NOT NULL`.
- Produces: `derive_features(rows: ExonRows, feature: FeatureType, per_gene: PerGene) -> Features`.
- Produces: `Reference.features(feature: FeatureType, transcript_filter: TranscriptFilter = TranscriptFilter(), *, per_gene: PerGene = PerGene.EACH_TRANSCRIPT) -> Features`.

The part of an exon `[s, e]` for each feature, with the CDS `[cs, ce]` (it includes the stop codon). A part with `lo > hi` is not a range, but its transcript stays in the result.

| Feature | Plus strand | Minus strand |
|---|---|---|
| `UTR3` | `[max(s, ce + 1), e]` | `[s, min(e, cs - 1)]` |
| `UTR5` | `[s, min(e, cs - 1)]` | `[max(s, ce + 1), e]` |
| `CDS` | `[max(s, cs), min(e, ce)]` | the same |
| `EXONS` | `[s, e]` | the same |

`UTR3`, `UTR5`, and `CDS` read the rows with `coding_only=True`.

- `EACH_TRANSCRIPT`: one feature for each transcript, in the `tx_id` order of the rows.
- `LONGEST`: the genes in sorted `gene_id` order. For each gene, the first transcript with the maximum sum of widths, in the order `(tx_seq_start, tx_row)`. The gene of a transcript has one contig, thus this is the `Reference.transcripts` order (`D-028`).
- `UNION`: the genes in sorted `gene_id` order. Sort the ranges of all transcripts of the gene by start. Merge a range into the previous range if its start is not more than the previous end + 1. Give the merged ranges in descending order on the minus strand (`D-027`). The ID of the feature is the gene ID.

- [ ] **Step 1: Write the failing tests**

```python
@pytest.mark.parametrize("config", CASES)
def test_features_match_r(config):
    case = CASES[config]
    with fixture_reference() as ref:
        got = ref.features(case.feature, case.transcript_filter, per_gene=case.per_gene).to_frame()
    expected = golden(f"annotation/{config}_features.tsv")
    assert got[["id", "rank", "contig", "start", "end", "strand"]].astype(str).values.tolist() == (
        expected[["name", "rank", "seqnames", "start", "end", "strand"]].values.tolist()
    )

def test_union_gives_the_minus_strand_ranges_in_transcript_order_D027(): ...
def test_longest_skips_a_gene_without_a_coding_transcript_D028(): ...
def test_per_gene_applies_to_the_exons_D029(): ...
def test_longest_tie_goes_to_the_lower_tx_row(): ...  # G01 gives T02
def test_union_merges_touching_ranges_but_not_a_gap_of_one_base(): ...  # G02
def test_transcript_without_a_3utr_has_no_range(): ...  # T07 and T12 are features with 0 ranges
```

- `D027`: the G02 ranges of `utr3_union_r_features.tsv` ascend. The Python G02 ranges are the same ranges in descending order.
- `D028`: with `protein_coding=False`, the result has no G03 and has each coding gene of the patched golden.
- `D029`: `exons_longest_r_features.tsv` is equal to `exons_all_features.tsv`. The Python result has one feature for each gene, and each feature is the longest exon set of its gene.

- [ ] **Step 2: Run `uv run pytest tests/annotation/test_features.py -v`.** Expected: FAIL with `ImportError`.

- [ ] **Step 3: Write `exon_rows`, `features.py`, and `Reference.features`.** The cut uses numpy on all rows. Add the `D-027`, `D-028`, and `D-029` entries to `docs/divergences.md`. The R evidence of each entry is the rnor7 probe (the design table) and the fixture golden.

- [ ] **Step 4: Run the tests.** Expected: PASS. Checkpoint.

### Task 6: Sequences, `from_directory`, the exports, and the documents

**Files:**
- Make: `src/seedmatch/annotation/extract.py`
- Modify: `src/seedmatch/annotation/reference.py`, `src/seedmatch/annotation/__init__.py`, `src/seedmatch/__init__.py`
- Modify: `tests/helpers.py` (`run_core_pipeline`), `README.md`, `tests/test_readme.py`, `docs/r-map.md`, `docs/divergences.md` (`D-002`, `D-009`)
- Test: `tests/annotation/test_reference.py`

**Interfaces:**
- Consumes: `Features`, `Reference.features` (Task 5), `Genome.resolve`, `reverse_complement_codes` (Task 2).
- Produces: `extract_sequences(features: Features, genome: Genome) -> SequenceSet` in `extract.py`. The `SequenceSet` has `ids = features.ids` and `gene_ids = features.gene_ids`.
- Produces: `Reference.sequences(feature: FeatureType, transcript_filter: TranscriptFilter = TranscriptFilter(), *, per_gene: PerGene = PerGene.EACH_TRANSCRIPT, min_width: int = 8) -> SequenceSet` and `Reference.from_directory(path: str | os.PathLike[str]) -> Reference`.
- Produces: the top-level exports `Reference`, `TranscriptFilter`, `FeatureType`, `PerGene`, and `Features`. `seedmatch.annotation` also exports `TwoBitGenome`, `FastaGenome`, `Genome`, and `open_genome`.

`extract_sequences` joins the codes of the ranges of each feature in list order. On the minus strand, it uses `reverse_complement_codes` of each range. A contig that `resolve` cannot find raises `InputError` with the contig and the genome path. `sequences` removes each feature with `widths() < min_width` before the extraction. `min_width` uses `as_limit`.

`from_directory` finds exactly one `*.sqlite` and one genome (`*.2bit`, or `*.fa.gz` with `.fai` and `.gzi`). An absent file or more than one candidate raises `InputError` that names the directory and the candidates.

- [ ] **Step 1: Write the failing tests**

```python
@pytest.mark.parametrize("config", CASES)
def test_sequences_match_r(config):
    case = CASES[config]
    with fixture_reference() as ref:
        got = ref.sequences(case.feature, case.transcript_filter,
                            per_gene=case.per_gene, min_width=case.min_width)
    expected = golden(f"annotation/{config}_sequences.tsv")
    assert list(got.ids) == expected["name_tx"].tolist()
    assert list(got.gene_ids) == expected["name_gene"].tolist()
    assert [got.sequence(i) for i in range(len(got))] == expected["sequence"].tolist()

def test_union_sequences_keep_the_gene_id_D002(): ...  # utr3_union_r gives NA for each name_gene
def test_sequence_ids_are_the_transcript_ids_D009(): ...
def test_from_directory_finds_the_ensdb_and_the_genome(tmp_path): ...
def test_from_directory_with_two_genomes_raises():  # ANNOTATION itself, match="more than one genome"
def test_from_directory_without_an_ensdb_raises(tmp_path): ...
def test_min_width_rejects_a_negative_value(): ...
def test_contig_absent_from_the_genome_raises(tmp_path): ...  # a copy of the EnsDb, G04 moved to "Y"
```

In `test_sequences_match_r`, the configs `utr3_union` and `utr3_union_min0` take `name_tx` as the expected gene IDs, because R gives `NA` to each `name_gene` (`D-002`).

- [ ] **Step 2: Run `uv run pytest tests/annotation/test_reference.py -v`.** Expected: FAIL with `AttributeError` for `sequences`.

- [ ] **Step 3: Write `extract.py`, `Reference.sequences`, `Reference.from_directory`, and the exports.**

- [ ] **Step 4: Update the documents and the core pipeline.**
  - `run_core_pipeline` also makes the 3' UTR sequences of the fixture and counts a target in them.
  - `README.md`: an example block with `Reference.from_paths` and `ref.sequences(FeatureType.UTR3)`. Remove EnsDb, 2bit, and FASTA from "not in this version", and add them to the scope. Task 7 adds the GTF path, because each task is a separate commit. `test_readme.py` copies `synthetic.sqlite` and `synthetic.2bit` to the names that the example uses.
  - `docs/r-map.md`: rows for `load_annotations`, `load_species_anno_db`, `build_annotation_filter`, `get_3UTR_granges`, `.get_longest_utr`, `.reduce_features`, `ensembldb::transcripts`, `extractTranscriptSeqs`, and `rtracklayer::TwoBitFile`. Record that a transcript without the feature has no range, but a range of width 0 in R. Record that `PerGene` replaces the two flags, and that `Reference.transcripts` replaces `$gtf` for all transcripts of the filter.
  - `docs/divergences.md`: `D-002` and `D-009`.

- [ ] **Step 5: Run `uv run pytest tests/annotation tests/test_readme.py tests/test_offline.py tests/test_package.py -v`.** Expected: PASS. Checkpoint.

### Task 7: GTF path

**Files:**
- Make: `src/seedmatch/annotation/gtf.py`
- Modify: `src/seedmatch/_io.py` (move `_decompressed` and `_open_text` from `de.py` as `open_text`), `src/seedmatch/de.py`, `src/seedmatch/annotation/reference.py`, `README.md` (the GTF path in the scope and in "not in this version")
- Test: `tests/annotation/test_gtf.py`

**Interfaces:**
- Consumes: `AnnotationDb(connection, source)`, `Reference`, `CASES` (Tasks 4 to 6).
- Produces: `seedmatch._io.open_text(path: Path) -> AbstractContextManager[TextIO]`, with the behavior of the old `_open_text`.
- Produces: `load_gtf(path: str | os.PathLike[str]) -> AnnotationDb` in `gtf.py`, and `Reference.from_gtf(gtf: str | os.PathLike[str], genome: str | os.PathLike[str]) -> Reference`.

`load_gtf` streams the lines through `open_text`. It skips the lines that start with `#`. A line without 9 tab-separated fields raises `InputError` with the line number. It uses only the `transcript`, `exon`, `CDS`, and `stop_codon` lines:

- A `transcript` line gives the `tx` row and, at the first transcript of a gene, the `gene` row. The attribute parser accepts quoted and unquoted values and a repeated `tag`.
- An `exon` line goes into a staging table at once, in batches of 10000 rows, thus the memory does not hold the exon lines.
- `CDS` and `stop_codon` lines update the minimum start and the maximum end of the CDS of their transcript.
- After the last line, `tx2exon.exon_idx` comes from `ROW_NUMBER() OVER (PARTITION BY tx_id ORDER BY <start ascending on +, descending on −>)`.
- An exon of a transcript without a `transcript` line raises `InputError` with the transcript ID.

The database has no `chromosome` table and no `entrezgene` table, thus `entrez_ids` raises `InputError`. `Reference.from_gtf` uses the check of `Reference.__init__` without lengths.

The GTF has the genes in the row order of the EnsDb, and it groups the transcripts by gene. The order rule of `D-030` uses the `gene` rows first, thus the GTF gives the tie of T10, T02, and T01 in the EnsDb order.

- [ ] **Step 1: Write the failing tests**

```python
GTF_CASES = {k: v for k, v in CASES.items() if v.transcript_filter.entrez_ids is None}
PAIRS = [("synthetic.gtf", "synthetic.2bit"), ("synthetic.gencode.gtf.gz", "synthetic.fa.gz")]

@pytest.mark.parametrize(("gtf", "genome"), PAIRS)
@pytest.mark.parametrize("config", GTF_CASES)
def test_gtf_gives_the_sequences_of_the_ensdb(gtf, genome, config): ...

@pytest.mark.parametrize(("gtf", "genome"), PAIRS)
def test_gtf_gives_the_transcripts_of_the_ensdb(gtf, genome): ...  # all columns except contig for GENCODE
def test_entrez_filter_on_a_gtf_raises():  # match="entrezgene"
def test_gtf_without_an_agreeing_contig_raises(tmp_path):  # contigs renamed to chrZ
def test_line_without_9_fields_raises_with_the_line_number(tmp_path):  # match="line 3"
def test_exon_without_a_transcript_line_raises(tmp_path):  # match="ENSSYNT00000000005"
def test_support_level_is_the_leading_integer(): ...  # "1 (assigned to previous version 2)" gives 1
```

- [ ] **Step 2: Run `uv run pytest tests/annotation/test_gtf.py -v`.** Expected: FAIL with `ImportError`.

- [ ] **Step 3: Move `open_text` to `_io.py`.** Run `uv run pytest tests/test_de.py`. Expected: PASS without a change to the tests.

- [ ] **Step 4: Write `gtf.py` and `Reference.from_gtf`.**

- [ ] **Step 5: Run the tests.** Expected: PASS. Checkpoint.

### Task 8: rnor7 parity and the real GTF

**Files:**
- Make: `tools/golden/R/annotation_rnor7.R`, `tests/annotation/test_refs.py`
- Make: `tests/golden/annotation/rnor7/<config>.tsv.gz` (written by the R script)
- Modify: `tools/golden/R/common.R` (`write_golden` writes gzip for a `.gz` path), `tools/golden/README.md`, `CONTRIBUTING.md`

**Interfaces:**
- Consumes: `Reference.from_directory`, `Reference.from_gtf`, `Reference.sequences` (Tasks 6 and 7).

The R script reads `file.path(Sys.getenv("SEEDMATCH_TEST_REFS"), "rnor7")` and stops if the variable is not set. It uses the replacement of `load_species_anno_db` of Task 1 and the patch of `.reduce_features`. Configs: `utr3`, `utr5`, `exons`, `cds`, `utr3_longest`, and `utr3_union`, with the R arguments of Task 1. Each file has the columns `name_tx`, `name_gene`, `width`, and `digest`: the first 16 hexadecimal digits of `openssl::sha256` of the sequence.

- [x] **Step 1: Download the Ensembl GTF.** (Done by the controller; the size and the digest are in the Global Constraints.) Run `curl -fsSL -o ~/Development/refs/seedmatch/rnor7/Rattus_norvegicus.mRatBN7.2.109.gtf.gz https://ftp.ensembl.org/pub/release-109/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.109.gtf.gz`. Record its size and SHA-256 in the report.

- [ ] **Step 2: Run the R script.** Run `SEEDMATCH_TEST_REFS=~/Development/refs/seedmatch tools/golden/run.sh tools/golden/R/annotation_rnor7.R`. Expected: 6 files. The `utr3` file has fewer than 23036 rows, because `min.feature.width` removes the 4168 sequences of width 0 and the sequences of width 1 to 7.

- [ ] **Step 3: Write the tests.** Each test has the marker `refs`. A fixture skips with the reason `SEEDMATCH_TEST_REFS is not set` or `<file> is absent`.

```python
def test_rnor7_files_are_the_pinned_files(rnor7): ...  # size and SHA-256 of the Global Constraints

@pytest.mark.parametrize("config", RNOR7_CONFIGS)
def test_rnor7_sequences_match_r(rnor7, config):
    expected = golden(f"annotation/rnor7/{config}.tsv.gz")
    # ids, gene_ids, widths, and digest16 of each sequence equal the columns, in order

@pytest.mark.parametrize("feature", list(FeatureType))
def test_rnor7_gtf_gives_the_sequences_of_the_ensdb(rnor7_gtf, feature): ...  # canonical=False
def test_rnor7_gtf_union_gives_the_sequences_of_the_ensdb(rnor7_gtf): ...
def test_rnor7_gtf_longest_gives_the_widths_of_the_ensdb(rnor7_gtf): ...
```

The GTF line order is not the EnsDb row order. Thus a `LONGEST` tie can select a different transcript of the same width, and the third consistency test compares only the widths.

- [ ] **Step 4: Run `SEEDMATCH_TEST_REFS=~/Development/refs/seedmatch uv run pytest -m refs -v`.** Expected: PASS. If the two paths disagree for a transcript, find the cause before you change a test. A data difference between the GTF and the EnsDb goes into the test as an explicit list with the reason. Record the wall time and the peak memory of `load_gtf` on the rat GTF in the report.

- [ ] **Step 5: Document.** In `tools/golden/README.md`, add the command of `annotation_rnor7.R` and the expected directory layout. In `CONTRIBUTING.md`, tell how to run the `refs` tests. Checkpoint.

### Task 9: Verification of the plan

- [ ] **Step 1: Run the checks one time.** Run `uv run ruff check`, `uv run ruff format --check`, `uv run mypy`, and `NUMBA_DISABLE_JIT=1 uv run pytest --cov`. Expected: each passes, and the coverage is not below 95%. Give a reason for each line of `seedmatch.annotation` without coverage.

- [ ] **Step 2: Compare the documents with the code.** Each public name of `seedmatch.annotation` has a row in `docs/r-map.md` or is a reader class that the README names. Each `D-NNN` of this plan has an entry and a test. `docs/design.md` agrees with the code, also with the ruling of Task 1 Step 6, if there is one. Checkpoint.
