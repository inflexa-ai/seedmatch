# Open defects: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the defects of `tmp/claude/SESSION.md` that the section "Findings of the planner probes" gives. Each public function must give a correct result or raise `InputError`, or `FileNotFoundError` for an absent file. A raw pandas, numpy, or matplotlib error or warning breaks this contract.

**Architecture:** Each task changes one module group and its tests. Tasks 3 and 4 change the line helpers of `seedmatch._lines` and the three table readers, in this order. The other tasks are independent. The plan adds no divergence ID. Task 7 changes the entry `D-013`.

**Tech Stack:** Python 3.12, numpy 2.5.2, pandas 2.3.3, matplotlib 3.11.1, SQLite through `sqlite3`, pytest. R 4.5.3 with SeedMatchR 2.0.0 and ggplot2 4.0.2 for the evidence probes. htslib `bgzip` and `samtools` for the fixture script of Task 7.

**Spec:** `docs/design.md`, sections "Rules", "Annotation", "Divergences (defect fixes)", and "Performance and bounds". `docs/divergences.md`, entry `D-013`.

## Findings of the planner probes

Each probe script and its output is in `tmp/claude/defects-plan/` of the planner. Each item of `tmp/claude/SESSION.md` gets one task.

### Item 1: a BGZF file with empty blocks at its start (Task 1)

Probe `probe_bgzf_empty_first.py`. `bgzip -r` of htslib writes no `.gzi` entry for an empty block, and no entry for the first block with data. `samtools faidx` reads each file below. `FastaGenome` reads the block at compressed offset 0 as the first data block:

- One or two empty blocks before many data blocks: `the block at compressed offset 0 has 0 bytes, but the .gzi gives 65273`.
- One or three empty blocks before one data block (the `.gzi` has no entry): `contig 's': the range 1..10 has 0 bases in the file; expected 10 bases, as the .fai says`. `tmp/claude/SESSION.md` does not give this second form.
- An empty block between two data blocks: valid, because the `.gzi` entry points to the data block after it.

Decision: at open, the reader finds the first block with data, from offset 0 to the first `.gzi` entry or to the end of the file.

### Item 2: feature widths that wrap int64 (Task 2)

Probe `probe_feature_widths.py`. The description of `tmp/claude/SESSION.md` is not correct in three points:

- The wrap starts at 3 ranges of `2^62 - 1` bases, not at 4. 2 ranges give `2^63 - 2`, which int64 holds.
- The offsets stay correct. The widths are wrong: 3 ranges give `-4611686018427387907`, 4 give `-4`, and 5 give `4611686018427387899`, a positive wrong value.
- The public symptom is a silent drop. On a copy of the fixture EnsDb with 3 exons `1..2^62 - 1` for `ENSSYNT00000000001`, `Reference.sequences` gives 8 sequences with `EACH_TRANSCRIPT` and 6 with `LONGEST`, and the transcript is absent. `at_least` drops the negative width. `LONGEST` can also select the wrong transcript.

The probe also found a second wrap in `_merge` (`src/seedmatch/annotation/features.py:251`). The band `gene * (max(ends) + 2)` overflows for 2 genes with an end near `2^62`. Then `UNION` keeps a range inside an earlier range: the gene with the ranges `(C - 100, C)` and `(C - 50, C - 40)` for `C = 2^62 - 1` keeps both ranges.

Decision: the sum of the range widths of one feature must be at most `2^62 - 1` (`COORDINATE_MAX`). `Features` and `LONGEST` raise `InputError` for a feature above it. `_merge` uses the ranks of the ends in its band. R `IRanges` holds 32-bit integers, thus R cannot hold these coordinates, and the plan adds no divergence.

### Item 3: the line number of a bad table line (Task 3)

Probe `probe_table_lines.py`. Each reader reads a file with a bad value on line 3 and a structural error on line 4. At the default chunk size, the message names line 4. At the chunk size 1, it names line 3. A structural error is a bad field count, a line above the length limit, or a NUL character.

A NUL character on line 3 and a structural error on line 4 also give line 4, because `data_bytes` looks for NUL only after `field_lines` accepts the chunk. `tmp/claude/SESSION.md` does not give this second form.

Decision: `field_lines` returns the data lines before the first bad line, and the error of that line. The reader parses these data lines first, and then it raises the error.

### Item 4: the pandas read of lone CR and U+FEFF (Task 4)

Probes `probe_table_lines.py`, `probe_pandas_fuzz.py`, `probe_cr_wrong_values.py`, `probe_cr_silent.py`, and `probe_string_feff.py`. pandas 2.3.3 reads the bytes of `data_bytes` in another way than the lines of `field_lines`:

- A lone CR line end before two lines that start with a space (a tab for STRING) gives the raw `pandas.errors.ParserError: Error tokenizing data. C error: Buffer overflow caught - possible malformed input file.` One such line is not sufficient. `tmp/claude/SESSION.md` does not give this condition.
- On 30,000 random chunks for each separator, the lone CR gives 673 to 700 raw `ParserError` and 11 to 27 chunks with wrong values. One chunk of 3 lines gives 131,073 rows. The readers raise in these cases, and no silent wrong table came from a reader.
- pandas removes a U+FEFF at the start of each input. Thus a U+FEFF at the start of the first data line of a chunk goes away. For the default chunk, `read_string_links` of the lines `protein1 protein2 combined_score`, `﻿9.A 9.B 600`, `9.C 9.D 500` gives the protein `9.A`.

Decision: `data_bytes` writes LF line ends and puts one blank line before the data. On 20,000 random files for each separator, this gives the values of `str.split` in each file. The current code differs in 1,372 to 1,465 files. A read of the miRDB 6.0 file (6,831,595 rows) takes 6.51 to 6.58 seconds before the change and 6.68 to 6.76 seconds after it, with Tasks 3 and 4 together. The peak RSS is 736 to 792 MiB before and 755 to 793 MiB after, thus the change is in the spread of the runs.

### Item 5: the Wasserstein overflow warning (Task 5)

Probe `probe_wasserstein_limits.py` and R probe `wasserstein_limits.R`. `_interpolate` (`src/seedmatch/stats/wasserstein.py:218-225`) divides by a segment width. A grid point outside the knots uses the first or the last segment, and a denormal width gives an infinite value. Then the two masks replace each such value with 0 or 1. Evidence:

- 2,276 random sample pairs, at the scales from 5e-324 to 4e307, give 436,278 values that are not finite. Each is outside `[first knot, last knot)`, and no value inside it is above 1.
- The unrounded areas and the frames of `wasserstein_difference` are bit-equal to a variant that never divides outside the knots, for each of the 2,276 pairs.
- R gives the areas of the port: `(4.992, 4.992, 0)` for a seed of `(0:9) * 5e-324` against `0:11`, and `(5.903, 0, 5.903)` against `-(0:11)`. For `(0:9) * 1e-300` against `(0:11) * 1e300`, R gives `4.9916666666666682e+300`, and the port gives `4.991666666666667e+300`.

Decision: the infinite value is exact in its use and the result is correct. Thus `_interpolate` divides inside `np.errstate(over="ignore")`, and `check_wasserstein_range` does not change.

### Item 6: a plot size that collapses the layout (Task 6)

Probes `probe_plot_sizes.py`, `probe_plot_layout.py`, `probe_plot_positions.py`, `probe_plot_cap.py`, and `probe_plot_hold.py`. matplotlib 3.11.1 warns `constrained_layout not applied because axes sizes collapsed to zero` when the figure draws. The smallest size without the warning, PNG at 100 dpi, with the test data of `tests/helpers.py`:

| Function | Default size | Smallest width (height 20) | Smallest height (width 40) |
|---|---|---|---|
| `plot_ranges` | (7, 5) | 0.684 | 0.730 |
| `plot_ecdf` | (5, 5) | 1.981 | 1.472 |
| `plot_kmer_histogram` | (7, 5) | 0.799 | 0.755 |
| `plot_wasserstein` | (15, 5) | 4.339 | 1.006 |
| `plot_benchmark` | (12, 5) | 2.301 | 0.947 |

One bound for each function is not possible:

- The bound moves with the dpi of the save. For `plot_benchmark`, the smallest width is 1.996 at 50 dpi and 2.308 at 300 dpi.
- The bound moves with the texts of the arguments. A group name of 10, 30, or 60 characters gives the smallest `plot_ecdf` width 1.981, 3.547, or 5.924. Thus the default size (5, 5) warns for a group name of 60 characters.
- The bound moves with the tick labels of the data. Thresholds near 9e305 move the smallest `plot_benchmark` width from 2.301 to 3.003.
- A layout at the figure dpi does not predict a save at another dpi. Near each bound, the two disagree at 8 to 13 of 78 sizes.

The current code also gives raw errors at the save. A side of 0.001 inches gives `ValueError: cannot write empty image`. A side of 10^6 inches gives `ValueError: Image size of 100000000x500 pixels is too large`. A side of 1e300 inches gives a `TypeError`.

Decision: each plot function applies the constrained layout one time before it returns. A collapse raises `InputError`. Then the figure keeps the positions of its Axes, thus a later save at each dpi does not lay out again. Near each bound, 30 of 63 sizes raise `InputError`. The other 33 save without a warning at 50, 72, 100, 150, 300, and 600 dpi, and to PDF and SVG.

`as_size` also rejects a side of 50 inches or more, as `ggsave` of ggplot2 4.0.2 does (`ggplot2:::plot_dim`, `any(dim >= 50)`, probe `ggsave_limit.R`). This bound keeps the layout canvas small. Costs:

- At 100 dpi, the PNG bytes of a default figure do not change. At 72 and 300 dpi, a position changes by at most 5.4e-3 of the figure.
- The build and one save take a little more time, before and after the change:
  - `plot_ranges` with 10^5 rows: 253.5 ms and 253.8 ms.
  - `plot_ecdf`: 48.6 ms and 52.7 ms.
  - `plot_wasserstein`: 100.4 ms and 104.7 ms.
- At 49.9 × 49.9 inches, the layout takes 40 to 56 ms. The peak memory grows by about 95 MiB for each such figure that the caller keeps.

### Item 7: the protein map of a GTF reference (Task 7)

Probes `probe_gtf_protein.py` and `probe_fixture_protein.py`. In `Homo_sapiens.GRCh38.109.gtf.gz`, only the `CDS` lines have the attribute `protein_id` (885,638 lines). R has no GTF path: `get_stringPPI_to_stringGGI` maps the proteins with `biomaRt::getBM` (`../SeedMatchR/R/get_stringPPI_to_stringGGI.R:47-50`). Thus the GTF map extends `D-013`, and the plan adds no ID. A `protein` table from the first `protein_id` of the `CDS` lines of each transcript gives these results:

- Human 109: 111,317 rows, one for each transcript with a CDS, with no conflict. Each row is a row of the EnsDb 109 `protein` table, and each protein has the gene of the EnsDb.
- The EnsDb has 12,072 rows more: 10,449 proteins of transcripts on the patch scaffolds (`CHR_*`) and 1,623 proteins of LRG transcripts. The primary-assembly GTF has neither.
- Rat 109: 45,936 rows, equal to the EnsDb 109 table. The default links of the STRING 12.5 rat file give the same network with each reference. `string_gene_network` gives 84,654 rows, 7,979 genes, and 628 unmapped proteins.
- `load_gtf` of the human GTF: 29.2 to 29.5 seconds and a peak RSS of 941 to 942 MiB before, and 29.8 to 29.9 seconds and 964 to 967 MiB after.
- The fixture script with a `protein_id` on each `CDS` line writes each other fixture file with the same bytes. The two GTF fixtures then give the 11 rows of the fixture EnsDb table, and the edges of the EnsDb reference.

## Global Constraints

- No new runtime dependency.
- Each public function gives a correct result or raises `InputError` (`FileNotFoundError` for an absent file). A raw pandas, numpy, or matplotlib error or warning is a defect.
- The message templates of this plan are normative. Each error test matches the full message or a stable part of it.
- The tests never run R and never use the network. The tests with the marker `refs` read the directory of `SEEDMATCH_TEST_REFS`.
- The plan adds no divergence ID. Task 7 changes the fields "New behavior" and "Test" of `D-013`.
- The R source of truth is `../SeedMatchR` at commit `2fc0c4f`. The R probes ran with `tools/golden/run.sh`.
- Documents use ASD-STE100. Code comments obey the admission test of the project rules. Documents, names, and comments obey the banned-word rule of the project.
- Run only the test files of a task in the task. Task 8 runs the full suite one time.
- The planner ran the test code of Tasks 1 to 6 on the current code. Each Step 2 gives the observed failures.
- Git: an implementer does no `git add` and no commit. Each task ends with a checkpoint. The controller commits only with the consent of the user.
- Each time comes from an Apple M2 Max with 12 CPUs and 32 GB, with other processes on the machine. Thus a time can be shorter on a quiet machine.

## Review Focus

1. A table file with a mix of line ends, for example LF and CR in one file. The random probes used one line end for each chunk, and no task test mixes them. Owner: Task 4.
2. A caller that sets `matplotlib.rcParams["figure.dpi"]` to a value other than 100. The layout canvas then has another size, thus the smallest valid size moves. No test changes `rcParams`. Owner: Task 6.
3. A caller that changes a returned figure, for example with `set_size_inches` or a new title. The figure keeps its positions and does not lay out again. The docstrings state this, but no test covers it. Owner: Task 6.
4. A GENCODE GTF of a real release. Its `protein_id` can keep a version suffix, as its other IDs do. Then no STRING ID matches, and `string_gene_network` raises `no link joins two genes of the annotation`. The machine has no GENCODE GTF. Owner: Task 7.
5. The human STRING links with a GTF reference. The machine has no human STRING file. The GTF map lacks the 12,072 patch and LRG proteins of the EnsDb, thus a link of such a protein can map with the EnsDb only. Owner: Task 7.

---

### Task 1: BGZF files with empty blocks at the start

**Files:**
- Modify: `src/seedmatch/annotation/fasta_genome.py`, `tests/annotation/test_fasta_genome.py`, `docs/design.md`

**Interfaces:**
- Produces: `FastaGenome._block_end(self, offset: int, size: int) -> tuple[int, int]` and `FastaGenome._first_data_block(self, limit: int, size: int) -> int`.
- Consumes: the test helpers `bgzf_block`, `write_genome`, `read_gzi`, `gzi_bytes`, and `fasta_bytes` of `tests/annotation/test_fasta_genome.py`.

- [ ] **Step 1: Write the failing tests.** Add to `tests/annotation/test_fasta_genome.py`, after `test_gzi_of_a_file_with_one_block_has_no_entries`:

```python
EMPTY_BLOCK = bgzf_block(b"")


@pytest.mark.parametrize("empty_blocks", [1, 2])
@pytest.mark.parametrize("block_size", [50, 1000], ids=["many blocks", "one block"])
def test_bgzf_file_that_starts_with_empty_blocks_is_valid(tmp_path, empty_blocks, block_size):
    sequences = {"a": "ACGT" * 30, "b": "TTGCA" * 20}
    path = write_genome(tmp_path, sequences, block_size=block_size)
    gzi = tmp_path / "genome.fa.gzi"
    lead = EMPTY_BLOCK * empty_blocks
    # The empty blocks move each later block, and they get no entry, as `bgzip -r` gives.
    pairs = [(compressed + len(lead), uncompressed) for compressed, uncompressed in read_gzi(gzi)]
    path.write_bytes(lead + path.read_bytes())
    gzi.write_bytes(gzi_bytes(pairs))
    assert len(EMPTY_BLOCK) == 28
    with FastaGenome(path) as genome:
        found = {name: genome.sequence(name, 1, len(bases)) for name, bases in sequences.items()}
    assert found == sequences


def test_bgzf_file_cut_inside_the_block_after_an_empty_block_raises(tmp_path):
    path = write_genome(tmp_path, {"a": "ACGT" * 10}, block_size=1000)
    path.write_bytes((EMPTY_BLOCK + path.read_bytes())[: len(EMPTY_BLOCK) + 30])
    message = (
        f"{os.fspath(path)!r}: the BGZF block at compressed offset 28 is damaged: "
        "the file ends inside it; expected an intact file"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        FastaGenome(path)
```

- [ ] **Step 2: Run the tests and see them fail.** Run `uv run pytest tests/annotation/test_fasta_genome.py -k "empty_block"`. Expected:
  - `many blocks`: `InputError` with `the block at compressed offset 0 has 0 bytes, but the .gzi gives 50`.
  - `one block`: `InputError` with `contig 'a': the range 1..120 has 0 bases in the file; expected 120 bases, as the .fai says`.
  - The cut file: `Failed: DID NOT RAISE <class 'seedmatch.errors.InputError'>`.

- [ ] **Step 3: Find the first block with data at open.** In `src/seedmatch/annotation/fasta_genome.py`:
  - Add `_block_end(self, offset: int, size: int) -> tuple[int, int]`. It returns the end of the block at `offset` and the ISIZE field of the block. It reads the header with `_read_header`, and then the last 4 bytes of the block. A block that ends after `size` raises `self._damaged(offset, "the file ends inside it")`.
  - `_check_uncompressed_offsets(self, size: int)` uses `_block_end` for each block. Remove its inline ISIZE read and the comment on the 4 bytes.
  - Add `_first_data_block(self, limit: int, size: int) -> int`. It walks the blocks from offset 0 with `_block_end`, and it returns the offset of the first block with an ISIZE above 0 before `limit`. If no block before `limit` has data, it returns 0, and the checks of `_read_gzi` give their present errors.
  - In `_read_gzi`, the first compressed offset is `self._first_data_block(offsets[0] if offsets else size, size)`. Replace the comment `# The first block has no entry in the file.` with `# bgzip -r writes no entry for the first block with data, and empty blocks can come before it.` Pass `size` to `_check_uncompressed_offsets`.
  - The walk reads one header and one ISIZE for each empty block at the start. `_check_uncompressed_offsets` already reads each header before the last entry, thus the open stays O(blocks).

- [ ] **Step 4: Run the tests and see them pass.** Run `uv run pytest tests/annotation/test_fasta_genome.py`. Expected: each test passes.

- [ ] **Step 5: Document.** In `docs/design.md`, section "GTF and FASTA", add this sentence to the bullet of the FASTA genome reader: "A BGZF file can start with empty blocks, for example after `cat` of an empty BGZF file and another file. The `.gzi` has no entry for the first block with data, thus the reader finds it at open." Checkpoint.

### Task 2: The width limit of a feature, and the union near the coordinate limit

**Files:**
- Modify: `src/seedmatch/annotation/features.py`, `tests/annotation/test_features.py`, `docs/design.md`

**Interfaces:**
- Produces: `range_totals(widths: NDArray[np.int64], offsets: NDArray[np.int64], ids: Sequence[str] | NDArray[np.object_]) -> NDArray[np.int64]` in `seedmatch.annotation.features`.
- Produces: the message `feature '<id>' has more than 2^62 - 1 bases in its ranges; expected at most 2^62 - 1`.
- Consumes: `COORDINATE_MAX` and `limit_text` of `seedmatch._io`. The test helpers `feature_parts`, `hand_rows`, `ranges_of`, `copy_ensdb`, `ALL`, and `ANNOTATION`.

- [ ] **Step 1: Write the failing tests.** Add to `tests/annotation/test_features.py`. Import `ANNOTATION` from `annotation_cases` and `Reference` from `seedmatch.annotation.reference`:

```python
C = 2**62 - 1
TOO_WIDE = "has more than 2^62 - 1 bases in its ranges; expected at most 2^62 - 1"


@pytest.mark.parametrize("count", [2, 3, 4, 5])
def test_features_reject_a_feature_with_more_than_2_62_minus_1_bases(count):
    parts = feature_parts(
        starts=np.array([10] + [1] * count, dtype=np.int64),
        ends=np.array([15] + [C] * count, dtype=np.int64),
        offsets=np.array([0, 1, 1 + count], dtype=np.int64),
    )
    with pytest.raises(InputError, match=re.escape(f"feature 'B' {TOO_WIDE}")):
        Features(**parts)


def test_features_accept_a_feature_of_exactly_2_62_minus_1_bases():
    parts = feature_parts(
        starts=np.array([10, 1, 1], dtype=np.int64),
        ends=np.array([15, 2**61, 2**61 - 1], dtype=np.int64),
    )
    assert Features(**parts).widths().tolist() == [6, C]


def test_features_reject_a_feature_of_2_62_bases():
    parts = feature_parts(
        starts=np.array([10, 1, 1], dtype=np.int64),
        ends=np.array([15, 2**61, 2**61], dtype=np.int64),
    )
    with pytest.raises(InputError, match=re.escape(f"feature 'B' {TOO_WIDE}")):
        Features(**parts)


WIDE_ROWS = hand_rows(
    ("T1", "G1", 1, None, ((1, C), (1, C), (1, C))),
    ("T2", "G1", 1, None, ((1, 100),)),
)


@pytest.mark.parametrize("per_gene", [PerGene.EACH_TRANSCRIPT, PerGene.LONGEST])
def test_derive_features_reject_a_transcript_with_more_than_2_62_minus_1_bases(per_gene):
    with pytest.raises(InputError, match=re.escape(f"feature 'T1' {TOO_WIDE}")):
        derive_features(WIDE_ROWS, FeatureType.EXONS, per_gene)


def test_union_of_exons_with_more_than_2_62_minus_1_bases_is_one_range():
    got = derive_features(WIDE_ROWS, FeatureType.EXONS, PerGene.UNION)
    assert ranges_of(got.to_frame(), "G1") == [(1, C)]


def test_union_merges_a_range_inside_a_range_that_ends_at_the_coordinate_limit():
    rows = hand_rows(
        ("T1", "GA", 1, None, ((1, 10), (C - 10, C))),
        ("T2", "GB", 1, None, ((C - 100, C),)),
        ("T3", "GB", 1, None, ((C - 50, C - 40),)),
    )
    frame = derive_features(rows, FeatureType.EXONS, PerGene.UNION).to_frame()
    assert ranges_of(frame, "GA") == [(1, 10), (C - 10, C)]
    assert ranges_of(frame, "GB") == [(C - 100, C)]


# ENSSYNT00000000001 has 3 exons. New exons keep the exons that other transcripts share.
WIDE_T01 = (
    *(f"INSERT INTO exon VALUES ('BIG{k}', 1, {C})" for k in (1, 2, 3)),
    "UPDATE tx2exon SET exon_id = 'BIG' || exon_idx WHERE tx_id = 'ENSSYNT00000000001'",
    f"UPDATE tx SET tx_seq_start = 1, tx_seq_end = {C} WHERE tx_id = 'ENSSYNT00000000001'",
)


@pytest.mark.parametrize("per_gene", [PerGene.EACH_TRANSCRIPT, PerGene.LONGEST])
def test_sequences_reject_a_transcript_with_more_than_2_62_minus_1_bases(tmp_path, per_gene):
    path = copy_ensdb(tmp_path, *WIDE_T01)
    with Reference.from_paths(ensdb=path, genome=ANNOTATION / "synthetic.2bit") as reference:
        with pytest.raises(InputError, match=re.escape(f"feature '{T01}' {TOO_WIDE}")):
            reference.sequences(FeatureType.EXONS, ALL, per_gene=per_gene)
```

- [ ] **Step 2: Run the tests and see them fail.** Run `uv run pytest tests/annotation/test_features.py -k "2_62 or coordinate_limit"`. Expected:
  - Each `reject` test: `Failed: DID NOT RAISE <class 'seedmatch.errors.InputError'>`.
  - `test_union_merges_a_range_inside_a_range_that_ends_at_the_coordinate_limit`: the assertion shows `[(4611686018427387803, 4611686018427387903), (4611686018427387853, 4611686018427387863)]` for `GB`.
  - The `accept` test and the `UNION` test of one gene pass before the change. They keep the bound inclusive and the union correct.
  - Make sure that the `copy_ensdb` statements run. A failure of `sqlite3` is not the reason under test.

- [ ] **Step 3: Bound the widths, and rank the ends of the union.** In `src/seedmatch/annotation/features.py`:
  - Add `range_totals`. It returns the sum of `widths[offsets[i]:offsets[i + 1]]` for each group i. A group above `COORDINATE_MAX` raises `InputError` with the message of the Interfaces. The time is O(n) for n widths, with one cumulative sum, and no Python loop over the groups.
  - Give `range_totals` this comment: `# A width is at most 2^63 - 1. Thus the partial sums of a group are exact up to the first one above COORDINATE_MAX, and that one is exact or negative after the wrap.`
  - `Features.__post_init__` calls `range_totals(self.ends - self.starts + 1, self.offsets, self.ids)` after the check of `starts > ends`. Then `widths()` is exact, because each total fits in int64.
  - The `LONGEST` branch of `derive_features` takes the transcript widths from `range_totals(row_widths, bounds, rows.tx_id[bounds[:-1]])`.
  - In `_merge`, replace the ends in the band with their dense ranks from `np.unique(ends, return_inverse=True)`. The band is `gene * len(levels)`, and `reach` maps the running maximum of the ranks back to `levels`. Replace the comment above the band with `# Each gene gets its own band of end ranks, above the bands of the genes before it. A rank and a band are below the square of the range count, thus they fit in int64.`
  - The extra sort of `np.unique` is O(n log n), the order of the `lexsort` of `_merge`.

- [ ] **Step 4: Run the tests and see them pass.** Run `uv run pytest tests/annotation/test_features.py tests/annotation/test_reference.py`. Expected: each test passes.

- [ ] **Step 5: Document.**
  - `docs/design.md`, section "EnsDb and 2bit", after the bullet of `min_width`: "The ranges of a feature hold at most 2^62 - 1 bases. `Features` and `PerGene.LONGEST` raise `InputError` for a feature with more bases, because the sums use int64."
  - The docstring of `Features`: "The ranges of a feature hold at most 2^62 - 1 bases, and a feature with more bases raises `InputError`."
  - Checkpoint.

### Task 3: The first bad line of a table file

**Files:**
- Modify: `src/seedmatch/_lines.py`, `src/seedmatch/mirdb.py`, `src/seedmatch/network/trrust.py`, `src/seedmatch/network/stringdb.py`
- Modify: `tests/test_lines.py`, `tests/test_mirdb.py`, `tests/network/test_trrust.py`, `tests/network/test_string_links.py`, `docs/design.md`

**Interfaces:**
- Produces: `class DataLines(NamedTuple)` in `seedmatch._lines`, with the fields `is_data: NDArray[np.bool_]` and `error: InputError | None`.
- Produces: `field_lines(chunk, first, source, *, separator, fields, described) -> DataLines`.
- Produces: `data_bytes(chunk: list[str], is_data: NDArray[np.bool_]) -> bytes`, and `parse_data_lines(chunk, is_data, *, separator, names, dtype=str) -> pd.DataFrame`. The arguments `first` and `source` go away.
- Consumes: the messages of `field_lines` and `data_bytes` stay as they are.

- [ ] **Step 1: Write the failing reader tests.** Add one test to each reader test file. Each test uses its file helpers `write_lines` and `failure`:

```python
# tests/test_mirdb.py
SCORE_ABC = "line 3: score 'abc' is not a finite number; expected a decimal number"
NUL_3 = "line 3: the line has a NUL character (U+0000); expected text without NUL characters"


@pytest.mark.parametrize("chunk", [None, 1, 2], ids=["default", "1", "2"])
@pytest.mark.parametrize(
    ("line3", "line4", "message"),
    [
        ("hsa-c\tNM_3\tabc", "hsa-d\tNM_4", SCORE_ABC),
        ("hsa-c\tNM_3\tabc", "hsa-d\tNM_4\t4\x000", SCORE_ABC),
        ("hsa-c\tNM_3\tabc", "hsa-d\tNM_4\t" + "4" * 30, SCORE_ABC),
        ("hsa-c\tNM_3\t3\x000", "hsa-d\tNM_4", NUL_3),
        ("hsa-c\tNM_3\t3\x000", "hsa-d\tNM_4\t" + "4" * 30, NUL_3),
    ],
    ids=["value, fields", "value, NUL", "value, long", "NUL, fields", "NUL, long"],
)
def test_read_mirdb_names_the_first_bad_line_for_each_chunk_size(
    tmp_path, monkeypatch, chunk, line3, line4, message
):
    if chunk is not None:
        monkeypatch.setattr("seedmatch.mirdb._CHUNK_LINES", chunk)
    monkeypatch.setattr("seedmatch._lines.LINE_LIMIT", 25)
    path = write_lines(tmp_path, ["hsa-a\tNM_1\t50", "hsa-b\tNM_2\t60", line3, line4])
    with pytest.raises(InputError, match=failure(path, message)):
        read_mirdb(path)
```

  `tests/network/test_trrust.py` gets `test_read_trrust_names_the_first_bad_line_for_each_chunk_size` with the same parameters and these values:
  - Lines 1 and 2: `"TFA\tGA\tRepression\t1"` and `"TFB\tGB\tActivation\t2"`. `_CHUNK_LINES` of `seedmatch.network.trrust`, and `LINE_LIMIT` 25.
  - The bad value of line 3: `"TFC\tGC\tactivation\t3"`, with the message `f"line 3: mode 'activation' is not valid; expected {MODES}"`. The NUL of line 3: `"TFC\tGC\tUnknown\t3\x00"`.
  - Line 4: `"TFD\tGD\tRepression"` (3 fields), `"TFD\tGD\tRepression\t4\x00"` (NUL), and `"TFD\tGD\tRepression\t" + "4" * 20` (long).

  `tests/network/test_string_links.py` gets `test_read_string_links_name_the_first_bad_line_for_each_chunk_size` with the same parameters and these values:
  - Lines 1 and 2: `HEADER` and `"9.A 9.B 500"`. `_CHUNK_LINES` of `seedmatch.network.stringdb`, and `LINE_LIMIT` 40, because the header has 33 characters with its line end.
  - The bad value of line 3: `"9.C 9.D abc"`, with the message `f"line 3: score 'abc' {NOT_A_SCORE}"`. The NUL of line 3: `"9.C 9.D 6\x0000"`.
  - Line 4: `"9.E 9.F"` (2 fields), `"9.E 9.F 7\x0000"` (NUL), and `"9.E 9.F " + "7" * 40` (long).

- [ ] **Step 2: Run the reader tests and see them fail.** Run `uv run pytest tests/test_mirdb.py tests/network/test_trrust.py tests/network/test_string_links.py -k first_bad_line`. Expected:
  - Chunk 1: each case passes.
  - Default chunk and chunk 2: each case fails with the `AssertionError` of `match`. The actual message names line 4, for example `line 4 has 2 fields; expected 3: the miRDB ID, the RefSeq ID, and the score`.

- [ ] **Step 3: Change the unit tests of `tests/test_lines.py` to the new interface.**
  - `data_flags` and `tab_lines` return the `DataLines`. Each test that read the array reads `.is_data`.
  - Each test that expected `field_lines` to raise reads `.error`. It asserts `str(lines.error) == message`, and `.is_data` is `False` from the bad line on. For example, `data_flags(["a b 1\n", "a b\n"], first=7)` gives `[True, False]` and `'f': line 8 has 2 fields; expected 3: x, y, and z`.
  - Add `test_field_lines_mark_no_data_line_from_the_first_bad_line_on`: `data_flags(["a b 1\n", "a b\n", "c d 2\n", "\n"])` gives `is_data == [True, False, False, False]` and the error of line 2.
  - Add `test_field_lines_check_the_length_then_the_fields_then_the_nul_of_a_line`. With `LINE_LIMIT` 10, `"a\0" + "b" * 10 + "\n"` gives the length message of line 1. With the default limit, `"a\0 b\n"` gives `line 1 has 2 fields`.
  - Move `test_data_bytes_name_the_line_of_a_nul_character` to `field_lines`, as `test_field_lines_name_the_line_of_a_nul_character`, with the same chunks, `first=10`, and the numbers 12 and 13.
  - `test_data_bytes_join_the_data_lines_in_utf8` calls `data_bytes(chunk, is_data)`. Task 4 changes its expected bytes.
  - The `parse_data_lines` tests drop the arguments `1, "f"` and `11, "f"`. Remove `test_parse_data_lines_name_the_line_of_a_nul_character`, because `field_lines` now owns the NUL check.
  - `test_field_lines_reject_fewer_than_two_fields` stays: a bad call still raises `ValueError`.
  - Run `uv run pytest tests/test_lines.py`. Expected: the changed tests fail, for example with `AttributeError: 'numpy.ndarray' object has no attribute 'error'` and `TypeError: data_bytes() missing 2 required positional arguments: 'first' and 'source'`.

- [ ] **Step 4: Return the first bad line, and raise it after the data lines.**
  - `src/seedmatch/_lines.py`: add `DataLines` with the docstring "The data lines of a chunk before its first bad line, and the error of that line, or `None`."
  - `field_lines` returns `DataLines`. On each line, the checks run in this order: the length, the field count of a line that is not blank, and a NUL character. The first bad line gives `error`, with the present messages. `is_data` is `False` for that line and for each line after it. The NUL check is one pass over the chunk, as the separator count is. A Python loop runs only over the lines that a check flags.
  - `data_bytes(chunk, is_data)` joins the data lines and returns the UTF-8 bytes. It no longer looks for NUL. `parse_data_lines(chunk, is_data, *, separator, names, dtype=str)` calls it.
  - Docstring of `field_lines`: "Return the data lines of the chunk and the error of its first bad line. A bad line has more than `LINE_LIMIT` characters, another field count, or a NUL character. The caller parses the data lines first, and then raises the error. Thus a bad value on an earlier line raises first."
  - `src/seedmatch/mirdb.py` and `src/seedmatch/network/trrust.py`: the loop calls `rows.add(chunk, lines.is_data, first, source, ...)`, and then raises `lines.error` if it is not `None`.
  - `src/seedmatch/network/stringdb.py`: `_Rows.add` keeps `_check_header` first. Move the body after `is_data[0] = False` into `_Rows._keep(self, chunk, is_data, first, source)`, which keeps the early returns. `add` calls `_keep`, and then raises `lines.error` if it is not `None`. `_read_chunk` loses the argument `first`. Remove its clause `except InputError: raise`, because `parse_data_lines` no longer raises `InputError`.

- [ ] **Step 5: Run the tests and see them pass.** Run `uv run pytest tests/test_lines.py tests/test_mirdb.py tests/network/test_trrust.py tests/network/test_string_links.py`. Expected: each test passes. `test_read_mirdb_numbers_a_bad_line_across_chunks` and the other present message tests pass without a change.

- [ ] **Step 6: Document.**
  - `docs/design.md`, section "Rules", a new bullet: "A reader of a table file (`read_mirdb`, `read_trrust`, `read_string_links`) names the first bad line of the file, for each chunk size."
  - The docstring of each of the three readers, after "The line number in a message counts the blank lines.": "The message names the first bad line of the file."
  - Checkpoint.

### Task 4: The line ends and U+FEFF of the table readers

**Files:**
- Modify: `src/seedmatch/_lines.py`, `tests/test_lines.py`, `tests/test_mirdb.py`, `tests/network/test_trrust.py`, `tests/network/test_string_links.py`, `src/seedmatch/mirdb.py`, `src/seedmatch/network/trrust.py`, `src/seedmatch/network/stringdb.py`, `docs/design.md`

**Interfaces:**
- Consumes: `data_bytes(chunk, is_data) -> bytes` of Task 3.
- Produces: the bytes of `data_bytes` start with `\n`, and each line end is `\n`.

- [ ] **Step 1: Write the failing tests.** In `tests/test_lines.py`, change the expected bytes of `test_data_bytes_join_the_data_lines_in_utf8` to `"\na b 1\né b 2\n".encode()`, and add:

```python
@pytest.mark.parametrize("ending", ["\n", "\r\n", "\r"])
def test_data_bytes_give_lf_line_ends_after_a_blank_line(ending):
    chunk = [f"a b 1{ending}", ending, f"﻿c d 2{ending}", "e f 3"]
    is_data = data_flags(chunk).is_data
    assert data_bytes(chunk, is_data) == "\na b 1\n﻿c d 2\ne f 3".encode()
```

  Add these tests to `tests/network/test_trrust.py`:

```python
ENDINGS = pytest.mark.parametrize("ending", ["\n", "\r\n", "\r"], ids=["LF", "CRLF", "CR"])
SPACE_LINES = ["TFA\tGA\tRepression\t1", " TFB\tGB\tActivation\t2", " TFC\tGC\tUnknown\t3"]


def write_ends(tmp_path, lines, ending):
    path = tmp_path / "x.tsv"
    path.write_bytes("".join(f"{line}{ending}" for line in lines).encode())
    return path


@ENDINGS
def test_read_trrust_reads_each_line_end_before_lines_that_start_with_a_space(tmp_path, ending):
    got = read_trrust(write_ends(tmp_path, SPACE_LINES, ending))
    assert got["source"].tolist() == ["TFA", " TFB", " TFC"]


@ENDINGS
def test_read_trrust_names_a_bad_line_after_lines_that_start_with_a_space(tmp_path, ending):
    path = write_ends(tmp_path, [*SPACE_LINES, "TFD\tGD\tactivation\t4"], ending)
    message = f"line 4: mode 'activation' is not valid; expected {MODES}"
    with pytest.raises(InputError, match=failure(path, message)):
        read_trrust(path)


@pytest.mark.parametrize(
    ("text", "chunk", "sources"),
    [
        ("TFA\tGA\tRepression\t1\n﻿TFB\tGB\tActivation\t2\n", 1, ["TFA", "﻿TFB"]),
        ("﻿﻿TFA\tGA\tRepression\t1\nTFB\tGB\tActivation\t2\n", None, ["﻿TFA", "TFB"]),
        ("\n﻿TFA\tGA\tRepression\t1\nTFB\tGB\tActivation\t2\n", None, ["﻿TFA", "TFB"]),
    ],
    ids=["first line of a chunk", "after the byte order mark", "after a blank line"],
)
def test_read_trrust_keeps_a_u_feff_after_the_start_of_the_file(
    tmp_path, monkeypatch, text, chunk, sources
):
    if chunk is not None:
        monkeypatch.setattr("seedmatch.network.trrust._CHUNK_LINES", chunk)
    path = tmp_path / "x.tsv"
    path.write_text(text, encoding="utf-8")
    assert read_trrust(path)["source"].tolist() == sources
```

  `tests/test_mirdb.py` gets the same tests and helpers, named with `read_mirdb`, and these values:
  - The lines `["hsa-a\tNM_1\t50", " hsa-b\tNM_2\t60", " hsa-c\tNM_3\t70"]` give `mirdb_id == ["hsa-a", " hsa-b", " hsa-c"]`.
  - The bad line `"hsa-d\tNM_4\tabc"` gives `"line 4: score 'abc' is not a finite number; expected a decimal number"`.
  - The U+FEFF cases use `"hsa-a\tNM_1\t50\n"` and `"hsa-b\tNM_2\t60\n"` in the three forms of TRRUST, with `_CHUNK_LINES` of `seedmatch.mirdb`. They give `["hsa-a", "﻿hsa-b"]`, `["﻿hsa-a", "hsa-b"]`, and `["﻿hsa-a", "hsa-b"]` in `mirdb_id`.

  `tests/network/test_string_links.py` gets the same tests and helpers, named with `read_string_links` and with `a_tab` in place of `a_space`. A read that gives a table uses `score_threshold=0, percentile=0`. The values:
  - The lines `[HEADER, "9.A 9.B 500", "\t9.C 9.D 600", "\t9.E 9.F 700"]` give `protein1 == ["\t9.C", "\t9.E"]`. The cutoff of `percentile=0` is the lowest score, 500, and it drops the first link.
  - The bad line `"9.G 9.H abc"` gives `f"line 5: score 'abc' {NOT_A_SCORE}"`.
  - The text `f"{HEADER}\n9.A 9.B 500\n﻿9.C 9.D 600\n"` with the chunk 2 gives `protein1 == ["﻿9.C"]`.
  - The text `f"{HEADER}\n﻿9.A 9.B 600\n9.C 9.D 500\n"` with the default chunk gives `protein1 == ["﻿9.A"]`. The probe `probe_string_feff.py` gives `9.A` before the change.

- [ ] **Step 2: Run the tests and see them fail.** Run `uv run pytest tests/test_lines.py tests/test_mirdb.py tests/network/test_trrust.py tests/network/test_string_links.py -k "line_end or u_feff or blank_line or start_with"`. Expected:
  - The `CR` cases of each reader: `pandas.errors.ParserError: Error tokenizing data. C error: Buffer overflow caught - possible malformed input file.`
  - The U+FEFF cases: an `AssertionError`, for example `['TFA', 'TFB'] == ['TFA', '﻿TFB']`.
  - The `data_bytes` tests: the bytes differ in the first `\n` and in the CR line ends.

- [ ] **Step 3: Write LF line ends after a blank line.** In `data_bytes` of `src/seedmatch/_lines.py`:
  - The text starts with `\n`. pandas skips the blank line, and a U+FEFF at the start of the first data line stays text.
  - If the text holds `\r`, `\r\n` and then `\r` become `\n`. `open_text` reads with `newline=""`, thus a `\r` can only be a line end.
  - Docstring: "Return the data lines of the chunk as UTF-8 bytes, after one blank line, and with LF line ends. pandas 2.3.3 removes a U+FEFF at the start of its input, and it misreads a lone CR before a line that starts with a space."
  - Keep the comment on the bytes and the `StringIO`.

- [ ] **Step 4: Run the tests and see them pass.** Run `uv run pytest tests/test_lines.py tests/test_mirdb.py tests/network/test_trrust.py tests/network/test_string_links.py`. Expected: each test passes.

- [ ] **Step 5: Document.**
  - `docs/design.md`, section "Rules", add to the bullet of Task 3: "A line end is LF, CRLF, or CR. The reader removes a byte order mark only at the start of the file, and a U+FEFF at another place is text."
  - The docstring of each of the three readers, after the sentence of Task 3: "A line end is LF, CRLF, or CR. A U+FEFF after the start of the file is text."
  - Checkpoint.

### Task 5: The Wasserstein ECDF of a denormal spread

**Files:**
- Modify: `src/seedmatch/stats/wasserstein.py`, `tests/stats/test_wasserstein.py`, `tests/plot/test_plot_wasserstein.py`

**Interfaces:**
- Consumes: `wasserstein_areas`, `wasserstein_difference`, and `plot_wasserstein`, without a change of signature.

- [ ] **Step 1: Write the failing tests.** Add to `tests/stats/test_wasserstein.py`:

```python
DENORMAL_SEED = np.arange(10) * 5e-324
# The areas that wass_dist of SeedMatchR 2.0.0 gives under R 4.5.3.
DENORMAL_CASES = [
    pytest.param(DENORMAL_SEED, np.arange(12.0), (4.992, 4.992, 0.0), id="above the seed"),
    pytest.param(DENORMAL_SEED, -np.arange(12.0), (5.903, 0.0, 5.903), id="below the seed"),
]


@pytest.mark.parametrize(("seed", "non_seed", "expected"), DENORMAL_CASES)
def test_wasserstein_gives_the_r_areas_without_a_warning_for_a_denormal_spread(
    seed, non_seed, expected
):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        areas = wasserstein_areas(seed, non_seed, n_interp=100)
    assert (areas.absolute, areas.negative_censored, areas.positive_censored) == expected


def test_wasserstein_gives_the_r_area_of_a_tiny_seed_against_a_huge_non_seed():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        areas = wasserstein_areas(np.arange(10) * 1e-300, np.arange(12) * 1e300, n_interp=100)
    assert areas.absolute == pytest.approx(4.9916666666666682e300, rel=1e-12)
    assert areas.positive_censored == 0.0


@pytest.mark.parametrize(("seed", "non_seed", "expected"), DENORMAL_CASES)
def test_wasserstein_difference_gives_no_warning_for_a_denormal_spread(seed, non_seed, expected):
    del expected
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        got = wasserstein_difference(seed, non_seed, n_interp=100)
    assert np.isfinite(got.to_numpy(dtype=np.float64)).all()
```

  In `tests/plot/test_plot_wasserstein.py`, remove the mark `@pytest.mark.filterwarnings("ignore:overflow encountered in divide:RuntimeWarning")` and the comment above it from `test_plot_wasserstein_rejects_a_sample_whose_bandwidth_is_zero`. Put the `pytest.raises` block inside `warnings.catch_warnings()` with `warnings.simplefilter("error")`.

- [ ] **Step 2: Run the tests and see them fail.** Run `uv run pytest tests/stats/test_wasserstein.py tests/plot/test_plot_wasserstein.py -k "denormal or huge_non_seed or bandwidth_is_zero"`. Expected: each case fails with `RuntimeWarning: overflow encountered in divide`, also the `seed` and the `non_seed` cases of the plot test. The planner ran these tests on the current code with this result. The areas are correct with the warnings off.

- [ ] **Step 3: Divide inside `np.errstate`.** In `_interpolate` of `src/seedmatch/stats/wasserstein.py`, put only the line `values /= curve.widths[segment]` inside `with np.errstate(over="ignore"):`. Put this comment above it: `# Outside the knots, a denormal width can give an infinite value. The two masks below replace each such value with 0 or 1, thus no infinity reaches a result.` Do not change `check_wasserstein_range`.

- [ ] **Step 4: Run the tests and see them pass.** Run `uv run pytest tests/stats/test_wasserstein.py tests/plot/test_plot_wasserstein.py`. Expected: each test passes. The R golden tests stay bit-equal, because the change touches no value.

- [ ] **Step 5: Document.** The code comment of Step 3 holds the fact. `docs/design.md` does not change, because the result does not change. Checkpoint.

### Task 6: The layout of the plot figures

**Files:**
- Modify: `src/seedmatch/plot/_style.py`, `src/seedmatch/plot/ranges.py`, `src/seedmatch/plot/ecdf.py`, `src/seedmatch/plot/kmers.py`, `src/seedmatch/plot/wasserstein.py`, `src/seedmatch/plot/benchmark.py`
- Modify: `tests/plot/test_plot_package.py`, `tests/plot/test_plot_ranges.py`, `tests/plot/test_plot_ecdf.py`, `tests/plot/test_plot_kmers.py`, `tests/plot/test_plot_wasserstein.py`, `tests/plot/test_plot_benchmark.py`, `docs/design.md`

**Interfaces:**
- Produces: `apply_layout(figure: Figure) -> Figure` in `seedmatch.plot._style`.
- Produces: the message `size (0.3, 0.3) is too small for the texts of the figure; expected a larger width or height`.
- Produces: the message of `as_size`, `size (0, 5) is not valid; expected (width, height) in inches, two finite numbers above 0 and below 50`.

- [ ] **Step 1: Write the failing tests.** In `tests/plot/test_plot_package.py`:
  - In `test_as_size_rejects_a_bad_size`, add the sizes `(50, 5)`, `(5, 50.0)`, and `(10**6, 5)`, and end the message with `two finite numbers above 0 and below 50`.
  - Add `test_as_size_accepts_a_side_just_below_50`: `as_size((49.99, 0.01)) == (49.99, 0.01)`.
  - Add the tests below. Import `io`, `warnings`, `pandas as pd`, the helpers `sirna_de`, `sirna_seqs`, and `wasserstein_samples`, the functions `benchmark_predictions`, `ecdf_test`, `kmer_counts`, and `kmer_gene_counts` of `seedmatch`, and the five plot functions.

```python
def sized_calls():
    """Return a call of each plot function that takes `size`, on the test data."""
    de = sirna_de()
    groups = {
        "Background": de.loc[de["mer7m8"] == 0, "gene_id"],
        "mer7m8": de.loc[de["mer7m8"] >= 1, "gene_id"],
    }
    result = ecdf_test(de, target=groups["mer7m8"], background=groups["Background"])
    seed, non_seed = wasserstein_samples("test_data")
    counts = kmer_gene_counts(kmer_counts(sirna_seqs()))
    ranges = pd.DataFrame({"start": [1, 5, 9], "end": [3, 9, 20]})
    return {
        "plot_ranges": lambda **options: plot_ranges(ranges, **options),
        "plot_ecdf": lambda **options: plot_ecdf(de, groups, result, **options),
        "plot_kmer_histogram": lambda **options: plot_kmer_histogram(counts, "CTCTATA", **options),
        "plot_wasserstein": lambda **options: plot_wasserstein(
            seed, non_seed, n_interp=100, **options
        ),
        "plot_benchmark": lambda **options: plot_benchmark(benchmark_predictions(de), **options),
    }


SIZED = list(sized_calls())


def too_small(size):
    return f"size {size!r} is too small for the texts of the figure; expected a larger width or height"


@pytest.mark.parametrize("name", SIZED)
def test_each_plot_function_rejects_a_size_too_small_for_its_texts(name):
    with pytest.raises(InputError, match=exactly(too_small((0.3, 0.3)))):
        sized_calls()[name](size=(0.3, 0.3))


def test_plot_ecdf_rejects_a_legend_too_wide_for_the_default_size():
    de = sirna_de()
    groups = {"g" * 60: de.loc[de["mer7m8"] >= 1, "gene_id"]}
    groups["Background"] = de.loc[de["mer7m8"] == 0, "gene_id"]
    result = ecdf_test(de, target=groups["g" * 60], background=groups["Background"])
    with pytest.raises(InputError, match=exactly(too_small((5.0, 5.0)))):
        plot_ecdf(de, groups, result)


# The color bars of plot_benchmark take their width at the first draw, thus the test draws one
# time before it reads the positions.
@pytest.mark.parametrize("name", SIZED)
def test_each_plot_function_keeps_its_layout_for_each_save(name):
    figure = sized_calls()[name]()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        figure.savefig(io.BytesIO(), format="png", dpi=100)
        positions = [axes.get_position().bounds for axes in figure.axes]
        for file_format, dpi in (("png", 50), ("png", 300), ("pdf", 72), ("svg", 72)):
            figure.savefig(io.BytesIO(), format=file_format, dpi=dpi)
            assert [axes.get_position().bounds for axes in figure.axes] == positions
```

  In `tests/plot/test_plot_ranges.py`, `tests/plot/test_plot_ecdf.py`, `tests/plot/test_plot_kmers.py`, `tests/plot/test_plot_wasserstein.py`, and `tests/plot/test_plot_benchmark.py`, end each size message with `two finite numbers above 0 and below 50`.

- [ ] **Step 2: Run the tests and see them fail.** Run `uv run pytest tests/plot/`. Expected:
  - `test_each_plot_function_rejects_a_size_too_small_for_its_texts` and the ECDF legend test: `Failed: DID NOT RAISE <class 'seedmatch.errors.InputError'>`.
  - `test_each_plot_function_keeps_its_layout_for_each_save`: an `AssertionError` after the save at 50 dpi, because the layout runs again at each save.
  - The size message tests: the `AssertionError` of `match`, and `DID NOT RAISE` for the sides of 50 or more.

- [ ] **Step 3: Bound the size, and apply the layout one time.** In `src/seedmatch/plot/_style.py`:
  - `as_size` rejects a side of 50 or more, with the message of the Interfaces. Add `_MAX_SIDE: Final = 50.0` with this comment: `# ggsave of ggplot2 4.0.2 rejects a side of 50 inches or more. The bound also keeps the canvas of apply_layout below about 95 MiB at 100 dpi.`
  - Add `apply_layout(figure: Figure) -> Figure`. It runs `figure.get_layout_engine().execute(figure)` one time, inside `warnings.catch_warnings()` with `warnings.filterwarnings("error", message="constrained_layout not applied", category=UserWarning)`. Another warning stays a warning.
  - The `UserWarning` of a collapse becomes `InputError` with the message of the Interfaces, `raise ... from None`. The size of the message is `tuple(figure.get_size_inches().tolist())`.
  - After the layout, `figure.set_layout_engine("none")` keeps the positions for each later draw. Docstring: "Apply the constrained layout of the figure one time, and keep the positions of the Axes for each later draw. Thus a save at each dpi gives no layout warning. A size too small for the texts raises `InputError`."
  - `plot_ranges`, `plot_ecdf`, `plot_kmer_histogram`, `plot_wasserstein`, and `plot_benchmark` end with `return apply_layout(figure)`. `plot_alignment` takes no `size`, and it does not change.

- [ ] **Step 4: Run the tests and see them pass.** Run `uv run pytest tests/plot/`. Expected: each test passes. `test_new_figure_has_the_size_and_the_constrained_layout` stays as it is, because `new_figure` does not change.

- [ ] **Step 5: Document.**
  - The `size` sentence of each of the five docstrings: "`size` is the width and the height in inches, each above 0 and below 50, as `ggsave` requires. A size too small for the texts, the legends, and the labels raises `InputError`. The figure keeps its layout, thus a later change of the figure does not move the Axes. To lay out again, call `figure.set_layout_engine("constrained")`."
  - `docs/design.md`, section "Rules", a new bullet: "A plot function applies its layout one time before it returns. A size too small for the texts raises `InputError`, and a later save at each dpi keeps the positions."
  - `docs/design.md`, section "Performance and bounds", a new bullet with the costs of Item 6 of "Findings of the planner probes": the build and one save of each measured function, and the 40 to 56 ms and 95 MiB at 49.9 × 49.9 inches. State that the layout canvas holds 4 bytes for each pixel at the figure dpi.
  - Checkpoint.

### Task 7: The protein table of a GTF reference

**Files:**
- Modify: `src/seedmatch/annotation/gtf.py`, `src/seedmatch/annotation/ensdb.py`, `src/seedmatch/network/stringdb.py` (docstring)
- Modify: `tools/fixtures/annotation.py`, `tests/data/annotation/synthetic.gtf`, `tests/data/annotation/synthetic.gencode.gtf.gz` (written by the script)
- Modify: `tests/annotation/test_gtf.py`, `tests/annotation/test_ensdb.py`, `tests/annotation/test_refs.py`, `tests/network/test_string_network.py`
- Modify: `docs/design.md`, `docs/divergences.md` (`D-013`), `tests/data/README.md`

**Interfaces:**
- Produces: a table `protein (tx_id TEXT, protein_id TEXT)` with the indexes `protein_tx_id_idx` and `protein_protein_id_idx` in the database of `load_gtf`, if a `CDS` line has a `protein_id`.
- Produces: the message `'<source>': line <n>: the CDS line of transcript 'T1' has the protein_id 'P2', but an earlier CDS line has 'P1'; expected one protein_id for each transcript`.
- Produces: the message of `AnnotationDb.protein_genes` without a table, `'<source>': the annotation has no table 'protein'; expected an EnsDb file with this table, or a GTF with the attribute protein_id on its CDS lines`.
- Consumes: `AnnotationDb.protein_genes` and `_OPTIONAL_COLUMNS["protein"]` of `src/seedmatch/annotation/ensdb.py`, without a change of the query.

- [ ] **Step 1: Write the failing tests.** In `tests/annotation/test_gtf.py`, with its helpers `line`, `write_gtf`, `TX1`, and `EXON1`:

```python
NO_PROTEIN_TABLE = (
    "the annotation has no table 'protein'; expected an EnsDb file with this table, "
    "or a GTF with the attribute protein_id on its CDS lines"
)


def coding(protein_id=None):
    protein = f' protein_id "{protein_id}";' if protein_id else ""
    return f'gene_id "G1"; transcript_id "T1";{protein}'


def test_load_gtf_maps_the_protein_id_of_the_cds_lines_to_the_gene(tmp_path):
    path = write_gtf(
        tmp_path,
        line("transcript", 1, 100, TX1),
        line("exon", 1, 100, EXON1),
        line("CDS", 10, 40, coding("P1")),
        line("CDS", 50, 60, coding("P1")),
        line("stop_codon", 61, 63, coding("P9")),
    )
    annotation = load_gtf(path)
    try:
        assert annotation.protein_genes(["P1", "P9"]) == {"P1": "G1"}
    finally:
        annotation.close()


def test_load_gtf_rejects_two_protein_ids_for_one_transcript(tmp_path):
    path = write_gtf(
        tmp_path,
        line("transcript", 1, 100, TX1),
        line("CDS", 10, 40, coding("P1")),
        line("CDS", 50, 60, coding("P2")),
    )
    message = (
        f"{str(path)!r}: line 3: the CDS line of transcript 'T1' has the protein_id 'P2', "
        "but an earlier CDS line has 'P1'; expected one protein_id for each transcript"
    )
    with pytest.raises(InputError, match=re.escape(message)):
        load_gtf(path)


def test_load_gtf_without_a_protein_id_has_no_protein_table(tmp_path):
    path = write_gtf(tmp_path, line("transcript", 1, 100, TX1), line("CDS", 10, 40, coding()))
    annotation = load_gtf(path)
    try:
        with pytest.raises(InputError, match=re.escape(f"{str(path)!r}: {NO_PROTEIN_TABLE}")):
            annotation.protein_genes(["P1"])
    finally:
        annotation.close()
```

  In `tests/annotation/test_ensdb.py`, end the message of `test_protein_genes_reject_an_annotation_without_a_protein_table` with `, or a GTF with the attribute protein_id on its CDS lines`.

  In `tests/network/test_string_network.py`, replace `test_string_gene_network_rejects_a_gtf_reference` with these tests. Import `re`:

```python
GTF_PAIRS = [("synthetic.gtf", "synthetic.2bit"), ("synthetic.gencode.gtf.gz", "synthetic.fa.gz")]


@pytest.mark.parametrize(("gtf", "genome"), GTF_PAIRS)
def test_string_gene_network_maps_with_the_protein_ids_of_a_gtf_D013(gtf, genome):
    expected = network_of()
    with Reference.from_gtf(ANNOTATION / gtf, ANNOTATION / genome) as reference:
        got = string_gene_network(sample_links(), reference)
    assert edges_of(got) == edges_of(expected) == edges((1, 4, 900), (4, 1, 900))
    assert got.unmapped_proteins == expected.unmapped_proteins


def test_string_gene_network_rejects_a_gtf_reference_without_protein_ids(tmp_path):
    path = tmp_path / "synthetic.gtf"
    text = (ANNOTATION / "synthetic.gtf").read_text()
    path.write_text(re.sub(r' protein_id "[^"]*";', "", text))
    with Reference.from_gtf(path, ANNOTATION / "synthetic.2bit") as gtf:
        message = (
            f"{gtf.annotation.source!r}: the annotation has no table 'protein'; "
            "expected an EnsDb file with this table, "
            "or a GTF with the attribute protein_id on its CDS lines"
        )
        with pytest.raises(InputError, match=exactly(message)):
            string_gene_network(sample_links(), gtf)
```

  In `tests/annotation/test_refs.py`, import `sqlite3`, and add:

```python
# The hg38 EnsDb also has the proteins of the patch scaffolds (CHR_*) and of the LRG transcripts,
# which the primary-assembly GTF lacks.
GTF_PROTEINS = {"rnor7": 45936, "hg38": 111317}


def test_gtf_gives_the_protein_genes_of_the_ensdb(refs_dir, build, ensdb_ref, gtf_ref):
    _, ensdb = lookup(Group.REFERENCE, build).resources
    uri = f"{(refs_dir / build / ensdb.name).resolve().as_uri()}?mode=ro&immutable=1"
    with sqlite3.connect(uri, uri=True) as connection:
        ids = [protein_id for (protein_id,) in connection.execute("SELECT protein_id FROM protein")]
    connection.close()
    expected = ensdb_ref.annotation.protein_genes(ids)
    got = gtf_ref.annotation.protein_genes(ids)
    assert len(got) == GTF_PROTEINS[build]
    assert got.items() <= expected.items()
```

- [ ] **Step 2: Run the tests and see them fail.** Run `uv run pytest tests/annotation/test_gtf.py tests/annotation/test_ensdb.py tests/network/test_string_network.py -k "protein"`. Expected:
  - The map test and the `D013` GTF tests: `InputError` with `the annotation has no table 'protein'; expected an EnsDb file with this table`.
  - The test of two protein IDs: `Failed: DID NOT RAISE <class 'seedmatch.errors.InputError'>`.
  - The message tests: the `AssertionError` of `match`, because the message lacks the GTF part.

- [ ] **Step 3: Put a `protein_id` on each `CDS` line of the fixture GTFs.**
  - In `attributes` of `tools/fixtures/annotation.py`, a `CDS` line of a transcript gets `("protein_id", q(tx.protein_id))`. In the Ensembl style, it comes after `transcript_biotype`. In the GENCODE style, it comes after `level`. The Ensembl 109 and GENCODE files use this order.
  - Run `uv run python tools/fixtures/annotation.py`, then `git status --short tests/data/annotation/`. The script uses `bgzip` and `samtools` from the `PATH`.
  - Expected: only `synthetic.gtf` and `synthetic.gencode.gtf.gz` change. The planner probe wrote each other file with the same bytes. If another file changes, stop and report it to the controller.

- [ ] **Step 4: Read the `protein_id` of the `CDS` lines.** In `src/seedmatch/annotation/gtf.py`:
  - `_Reader.__init__` adds `self._proteins: dict[str, str] = {}`, the protein of each transcript in the order of its first `CDS` line with a `protein_id`.
  - `_coding` parses the attributes one time. For a `CDS` line with a `protein_id` that is not empty, it keeps the ID for the transcript. Another `protein_id` for the same transcript raises the message of the Interfaces through `self._fail`. A `stop_codon` line gives no protein.
  - `finish` writes the table `protein` only if `self._proteins` is not empty. It inserts the rows, and then it makes the two indexes, as `_UNIQUE_IDS` does for the other tables.
  - The docstring of `load_gtf`: "The database has a `protein` table only if a `CDS` line has the attribute `protein_id`. A transcript gets the `protein_id` of its `CDS` lines, and two different IDs raise `InputError`."
  - In `src/seedmatch/annotation/ensdb.py`, `protein_genes` gives the message of the Interfaces. Its docstring names an annotation of `load_gtf` without `protein_id` attributes in place of "one from `load_gtf`".
  - In `src/seedmatch/network/stringdb.py`, the docstring of `string_gene_network` states that the map uses the `protein` table of the EnsDb or of the GTF.

- [ ] **Step 5: Run the tests and see them pass.** Run `uv run pytest tests/annotation/test_gtf.py tests/annotation/test_ensdb.py tests/network/test_string_network.py`. Then run `SEEDMATCH_TEST_REFS=~/Development/refs/seedmatch uv run pytest tests/annotation/test_refs.py -k protein_genes`. Expected: each test passes, and the refs test runs for `rnor7` and `hg38`.

- [ ] **Step 6: Document.**
  - `docs/divergences.md`, `D-013`, "New behavior": "`string_gene_network` maps the protein IDs with the tables `protein` and `tx` of the annotation of the `Reference`. For a GTF reference, `load_gtf` makes the `protein` table from the attribute `protein_id` of the `CDS` lines."
  - `D-013`, "Test": add `test_string_gene_network_maps_with_the_protein_ids_of_a_gtf_D013` in `tests/network/test_string_network.py`.
  - `docs/design.md`, section "GTF and FASTA", a new bullet: "The `protein` table comes from the attribute `protein_id` of the `CDS` lines. On Ensembl 109, each GTF protein has the gene of the EnsDb. The human EnsDb also has 12,072 proteins of patch scaffolds and LRG transcripts."
  - `docs/design.md`, the `protein_genes` sentence of the section "EnsDb and 2bit" and the row `get_stringPPI_to_stringGGI` of the divergence table: the map uses the `protein` table of the EnsDb or of the GTF.
  - `docs/design.md`, the `load_gtf` bullet of "Performance and bounds": "The `protein` table adds about 25 MiB and 0.5 seconds on the human GTF, for 111,317 rows."
  - `tests/data/README.md`, the bullet of the two GTF fixtures: "Each `CDS` line has the attribute `protein_id` of its transcript, the protein of the table `protein`."
  - Checkpoint.

### Task 8: Verification of the plan

- [ ] **Step 1: Run the checks one time.** Run `uv run ruff check`, `uv run ruff format --check`, and `uv run mypy`. Then run `uv run pytest --cov` one time, without `SEEDMATCH_TEST_REFS`. Expected: each passes, the refs tests skip, and the coverage stays at 95% or more. For each new branch without coverage, give the reason in the report.

- [ ] **Step 2: Compare the documents with the code.**
  - `git status --short` shows the files of the tasks only. In `tests/data/annotation/`, only the two GTF fixtures change.
  - Each message of the Interfaces of a task occurs in its source file and in a test.
  - `docs/divergences.md` has no new ID, and `D-013` names the new test.
  - The controller removes the items of this plan from `tmp/claude/SESSION.md`.
  - Checkpoint.
