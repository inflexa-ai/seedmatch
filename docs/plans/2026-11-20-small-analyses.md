# Small analyses: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port the small analyses of SeedMatchR to `seedmatch`. The output is identical to R, except for the recorded divergences. The R functions are `sequence_kmer_counts`, the counts of `kmer_hist`, `plotAlignment` with `output.type = "data.frame"`, `match_distances`, `generate_e_value`, `benchmark_preds`, and the miRDB table of `get_example_data`. The port also writes the PR and ROC curves of PRROC 1.4.

**Architecture:** Each analysis has its own module: `seedmatch.kmers`, `seedmatch.alignment`, `seedmatch.distances`, `seedmatch.evalue`, and `seedmatch.mirdb`. The PR and ROC curves are in `seedmatch.stats.curves`, and the benchmark is in `seedmatch.stats.benchmark`.

A function takes the types of the core: `SequenceSet`, `Guide`, `SeedSpec`, `SearchOptions`, a DE table, or a `find_matches` table. It returns a `DataFrame`, a `str`, or a frozen dataclass.

The benchmark goes in `seedmatch.stats`, because the design puts the PR and ROC curves there (`docs/design.md`, "Packages"). `benchmark_predictions` selects two classes from a DE table and calls the curves. In the same way, `ecdf_test` selects two samples from a DE table and calls a test. The miRDB reader goes in its own module, `seedmatch.mirdb`, because the miRDB table is not a DE table. The module holds the schema of one upstream resource, as `seedmatch.de` holds the schema of the GEO tables.

**Tech Stack:** Python 3.12, numpy, pandas, scipy 1.16.3 (`scipy.sparse`). R 4.5.3 with Biostrings 2.78.0, GenomicRanges 1.62.1, and PRROC 1.4 for the goldens.

**Spec:** `docs/design.md`, sections "Parity contract", "Architecture", "Divergences (defect fixes)", and "Performance and bounds".

**Later plans (not in this plan):**

- The k-mer histogram, the `ggplot` output of `plotAlignment`, and the PR and ROC plots (delivery step 9).
- The wobble and bulge search, the compensatory search, and the network scores.
- In R, `run_wobble_bulges_pipeline` makes the column `overall_matches` that `generate_e_value` reads. This plan gives the e-value of each count column of a DE table.

## Global Constraints

- Runtime dependencies: `numpy`, `pandas`, `numba`, and `scipy`. No new dependency. `seedmatch.kmers` uses `scipy.sparse.csr_array` of the present dependency `scipy`.
- No core module imports a network module or `seedmatch.fetch` (ruff `TID251`). No function writes a file or draws a plot.
- Errors derive from `seedmatch.errors.SeedMatchError`. Invalid input raises `InputError`, and each message names the bad value and the expected value.
- The message templates of this plan are normative. Each error test matches a stable part of its template. A template that names a helper of `seedmatch._checks` is the message of that helper.
- A float result agrees with R within a relative tolerance of 1e-9. Integers and strings are equal. An e-value is equal, because R and Python divide the same two integers.
- A PR or ROC value also agrees within an absolute tolerance of 1e-12 (`STATISTIC_TOLERANCE` of `tests/helpers.py`), because a value near 0 has no relative scale.
- Each change of behavior from R gets an entry `## D-NNN: <title>` in `docs/divergences.md`. It also gets a row in the table of `docs/design.md`, and a test whose name contains the ID.
- This plan writes the entries of the design rows `D-006`, `D-008`, `D-010`, and `D-011`, and of the new IDs `D-042` to `D-049`. Put each entry in ID order.
- A rename of an argument or a value, or a change of a type, is not a divergence. `docs/r-map.md` records it.
- The R source of truth is `../SeedMatchR` at commit `2fc0c4f`. The R scripts run with `tools/golden/run.sh`. The tests never run R and never use the network.
- Licenses (`docs/design.md`, "Licenses"): PRROC 1.4 has the GPL-3 license. Write the curves from Keilwagen, Grosse, and Grau (2014, PLOS ONE 9(3)) and from Davis and Goadrich (2006, ICML).
- The planner read `R/PRROC.R` of PRROC 1.4 to settle the details that the papers leave open. Task 6 states them as rules. Do not translate the PRROC source line by line.
- Write the k-mer count from the documentation of `Biostrings::oligonucleotideFrequency` and from the rules of Task 2. SeedMatchR has the MIT license, thus a golden patch can copy a SeedMatchR function.
- Each function obeys "Performance and bounds" of the design. Its docstring gives the time and the memory, and its task adds a line to that section of `docs/design.md`.
- Documents use ASD-STE100. Code comments obey the admission test of the project rules. The word "seam" does not occur.
- Git: an implementer does no `git add` and no commit. Each task ends with a checkpoint. The controller commits.

## Shared interfaces

A task makes each of these names, and a later task uses it:

```python
# seedmatch/seed.py (Task 3)
def seed_window(seed: SeedSpec, guide_length: int) -> tuple[int, int]: ...

# seedmatch/_checks.py (Task 4): keyword-only names for the messages. The defaults keep each caller.
def as_data_frame(value: object, *, name: str = "de") -> pd.DataFrame: ...
def check_columns(
    table: object, required: Sequence[str], *, name: str = "de", label: str = "DE table"
) -> None: ...

# seedmatch/_checks.py (Task 5): moved from seedmatch.stats.dts._as_generator, same messages.
def as_generator(rng: object) -> np.random.Generator: ...

# seedmatch/_pool.py (Task 5): moved from seedmatch.report.grid (_pool_width and _run).
def pool_width(workers: int | None) -> int: ...
def map_in_pool[S, T](work: Callable[[S], T], items: Sequence[S], width: int) -> Iterator[T]: ...
```

`seed_window` returns the first and the last guide position (1-based) that the target of `seed` pairs with. The `A` at the end of the target of `MER8` and `MER7A1` takes the position g1. The function uses `_canonical_window` of the same module. Let N be the length of the guide:

- `MER7M8` gives (2, 8), `MER8` (1, 8), `MER6` (2, 7), and `MER7A1` (1, 7).
- `FULL` gives (2, N), `MER18` (2, 19), and `MER15` (2, 16).
- `CustomSeed(name, start, stop)` gives (start, stop).

`map_in_pool` gives the result of each item in the order of the items. It yields a result while the pool works on the next items, thus a caller can write each result into a preallocated array. `report.grid` keeps its behavior with `list(map_in_pool(...))`. The callers of `as_data_frame` and `check_columns` give the table by position, thus the new parameter names change no caller.

## Review Focus

1. The k-mer counts of a real set, for example `ref.sequences(FeatureType.UTR3)` of `rnor7` with k = 7. Give the time and the peak memory of `kmer_counts` and `kmer_gene_counts`. `to_frame()` on this set is dense, about 3 GB, and its docstring must say so. Owner: Task 2.
2. A seed match near the start of a 3' UTR, and a mismatch match of `find_matches` that starts before 1. For a guide of 23 bases, the window of a `mer7m8` match starts 15 bases before the match. `alignment_text` shows spaces outside the sequence (`D-046`), and R stops. Owner: Task 3.
3. `e_value_null` with the default 1000 shuffles on the canonical 3' UTRs of `rnor7`, with `workers=None`. Give the wall time and the memory of the null. The README example uses a small `n_random`. Owner: Task 5.
4. `benchmark_predictions` on a DE table without `filter_de`. A row with a NaN `padj` or a NaN `log2FoldChange` is a negative, as in R. A count column that `annotate_seed_counts` makes with `shared_genes=False` gives 0 to a gene without a sequence. Owner: Task 7.
5. The real miRDB v6.0 file (59,483,982 bytes of gzip) with `species="rno"` and with `species=None`. Give the time, the peak memory, and the species codes of the file. No test reads this file. Owner: Task 8.

---

### Task 1: R goldens of the small analyses, and PRROC

**Files:**
- Modify: `tools/golden/install.R`, `tools/golden/R/check_env.R`, `tools/golden/R/common.R`, `tools/golden/R/annotation_common.R`, `tools/golden/README.md`, `tests/data/README.md`
- Make: `tools/golden/patches/plotAlignment.R`, `tools/golden/patches/benchmark_preds.R`, `tools/golden/R/analyses.R`
- Make: `tests/data/kmer_cases.tsv`, `tests/data/mirdb_sample.txt`
- Make: the files of `tests/golden/analyses/` (written by the script)

**Interfaces:**
- Produces: the goldens below. Each later task reads them with `tests/helpers.py::golden("analyses/<file>")`.

- [ ] **Step 1 (controller): Pin PRROC 1.4 in `install.R`.** PRROC is a Suggests of SeedMatchR, thus `install.R` does not install it now, and `tools/golden/.rlib` does not have it. `benchmark_preds` calls it. Add a block after the install of the imports:
  - If `lib` has no PRROC 1.4, download `https://cran.r-project.org/src/contrib/PRROC_1.4.tar.gz` to `tempdir()`. If this URL fails, use `https://cran.r-project.org/src/contrib/Archive/PRROC/PRROC_1.4.tar.gz`.
  - Stop unless `tools::md5sum` of the file is `229140aa270cdad059daf08783b083bb`. The SHA-256 is `9dbb4a7a9becb5a8f281a05fe263c4c96b187e459f1ed96b27381e838db493b6`. The planner computed the two digests from the CRAN file on 2026-11-20.
  - Install the file into `lib` with `repos = NULL, type = "source"`. PRROC depends on `rlang`, and the system library has `rlang` 1.1.7.
  - `plotAlignment` calls only SeedMatchR, Biostrings, and XVector. A planner probe ran it without `pwalign`, thus no other package is necessary.
  - `check_env.R` stops unless the PRROC version is exactly `1.4`. A prefix check with `"1.4."` cannot match, because the version has no third part.
  - Run `Rscript tools/golden/install.R`, then `tools/golden/run.sh tools/golden/R/check_env.R`. Expected: the last line starts with `golden environment OK` and names `PRROC 1.4`.
- [ ] **Step 2: Move `with_seedmatchr_patch`.** Move the function from `annotation_common.R` to `common.R` without a change. Each annotation script sources `common.R` before `annotation_common.R`, thus it keeps the function. In `tools/golden/README.md`, section "Write a golden file", move the bullet of the function to the list of `common.R`.
- [ ] **Step 3: Write the patches.** Each patch is a copy of one SeedMatchR function with one fix. Its first comment names the divergence and the R line.
  - `plotAlignment.R` (`D-045`): the guide target is `.get_custom_sequence(guide.seq, 1, guide.len)`, not `get_seed(guide.seq, "Full", ...)` (`R/plotAlignment.R:32`). Thus the guide string has the full guide, g1 to gN.
  - `benchmark_preds.R` (`D-006`): the filter compares the column `padj` with the argument (`padj <= !!padj`). `roc.curve` gets `weights.class0 = res$GS.label`, as `pr.curve` does.
- [ ] **Step 4: Write the fixtures.** Both fixtures are synthetic. `tests/data/kmer_cases.tsv` has a header line and these rows, with a tab between the fields. The sequence of `s4` is empty:

```text
seq_id	gene_id	sequence
s1	g1	ACGTNACGTA
s2	g1	acgtRacgt
s3	g2	AC
s4	g3	
s5	g4	AAAAAAAAAA
s6	g2	GATTACAGATTACA
s7	g5	TTTTMKTTTT
s8	g6	CGCGCGCGCGCGCGCGCGCG
```

  `tests/data/mirdb_sample.txt` has the format of the miRDB file: no header, a tab between the fields, and the fields ID, RefSeq ID, and score. The IDs are not copies of miRDB data:

```text
hsa-let-7a-5p	NM_000101	96.4871830265
hsa-let-7a-5p	NM_000102	50
hsa-miR-21-5p	NM_000103	54.527
hsa-miR-21-5p	XM_000104	99.7897051377461
mmu-miR-122-5p	NM_000201	88.1
mmu-miR-122-5p	NM_000202	62
rno-miR-192-5p	NM_000301	91.25
rno-miR-192-5p	XM_000302	77
rno-let-7a-2-3p	NM_000303	58.06
cfa-miR-1185	XM_000401	59.3438099752
gga-miR-1a-3p	NM_000501	71
lin4	NM_000601	55
```

- [ ] **Step 5: Write `tools/golden/R/analyses.R`.** The script sources `common.R` and writes each golden with `write_golden` into `tests/golden/analyses/`. The notes of a golden that uses PRROC give its version. Each section that draws random numbers starts with its own `set.seed`, thus a new case moves no other section. These are the sections:
  - **k-mers.** Read `kmer_cases.tsv` with `read.delim(..., colClasses = "character", na.strings = character(0))`, thus the empty sequence stays `""`. `kmer_counts.tsv` has the columns `case`, `seq_id`, `gene_id`, `kmer`, and `count`. It has one row for each count above 0 of `sequence_kmer_counts`, in the order of its rows, and then of its columns. The cases:
    - `test_data_first10_k7_s1`: the first 10 sequences of `Test_seqs.fa`, `width = 7`, `step = 1`. This is the call of the k-mer vignette. `tx2gene` maps each name to itself.
    - `test_data_k3_s2`: each sequence of `Test_seqs.fa`, `width = 3`, `step = 2`.
    - `cases_k1_s1`, `cases_k2_s1`, `cases_k3_s2`, and `cases_k4_s3`: the sequences of `kmer_cases.tsv`, with its gene IDs in `tx2gene`.
  - `kmer_frame.tsv`: the data frame of `sequence_kmer_counts` for `cases_k2_s1`, as R gives it. Its columns are `tx_id`, `gene_id`, and one column for each 2-mer.
  - `kmer_gene_counts.tsv` has the columns `case`, `kmer`, and `n`. `n` is `colSums` of the presence matrix of `R/kmer_hist.R:19-28`, for each k-mer in column order. The cases:
    - `test_data_k7_s1`: each sequence of `Test_seqs.fa`, with k = 7.
    - `cases_k2_s1`: the sequences of `kmer_cases.tsv`, with k = 2.
  - Notes of the k-mer section:
    - `n` of the `mer7m8` target `CTCTATA` for `test_data_k7_s1`.
    - The value of `kmer_hist(<cases>, "ACG", width = 3)`, with `pdf(NULL)` as the device, thus the script writes no `Rplots.pdf`. The planner probe gives `integer(0)`.
    - The error of `sequence_kmer_counts` for a set with a repeated name. The planner probe gives `duplicate 'row.names' are not allowed`.
    - `dim` of `oligonucleotideFrequency(DNAStringSet("ACGTACGTACGTAC"), width = 13)`. The planner probe gives 67108864 columns.
  - **Alignment.** `alignment.tsv` has the columns `case`, `seq_id`, `start`, `end`, `seed`, `start_pos`, `stop_pos`, `left_pad`, `right_pad`, `line` (1 to 3), and `text`. `text` is a value of the column `alignment` of the patched `plotAlignment`, with each space. The cases use `GUIDE` of `common.R` and the sequences `seqs` of `Test_seqs.fa`. In a case ID, `i` is the 1-based row of the match, and `seed` holds the R seed name:
    - `mer7m8_<i>`: each row of `vmatchPattern2(<mer7m8 target>, seqs)`, with the positions 2 and 8 and the pads 0 and 0.
    - `mer7m8_pad_<i>`: the first 3 of these rows, with `left.pad = 3` and `right.pad = 2`.
    - `mer8_<i>`: each match of the `mer8` target, with the positions 1 and 8.
    - `mer6_<i>`: the first 5 matches of the `mer6` target, with the positions 2 and 7.
    - `full_mm2_<i>`: each match of `Biostrings::vmatchPattern(<Full target>, seqs, max.mismatch = 2)` whose window is in its sequence, with the positions 2 and 23. These rows show mismatches. If this search gives no such match, use the smallest `max.mismatch` that gives one, and give the value in the notes.
    - `custom_<i>`: each match of `get_seed(GUIDE, "custom", start.pos = 3, stop.pos = 12)$Target.Seq`, with the positions 3 and 12. The `seed` column holds `custom`.
  - Notes of the alignment section: the R error for each of these inputs.
    - `mer7m8_1` with the unpatched `plotAlignment`. The planner probe gives `nchar(s1) == nchar(s2) is not TRUE`.
    - A row of `SeedMatchR(res.format = "data.frame")` with the unpatched function. The planner probe gives `subscript out of bounds`.
    - A window that starts before 1, with the patched function. The planner probe gives `solving row 1: the supplied start/end lead to a negative width`.
    - A window that ends after its sequence, with the patched function. The planner probe gives `solving row 1: 'allow.nonnarrowing' is FALSE and the supplied end (147) is > refwidth`.
    - A match with an indel, with the patched function. The planner probe gives `nchar(s1) == nchar(s2) is not TRUE`.
  - **Distances.** Each case splits its ranges by name with `factor(names, levels = unique(names))`, thus the list keeps the order of the first range of each name. The files:
    - `distances_ranges.tsv` has the columns `case`, `seq_id`, `start`, `end`, and `width`. These are the ranges that the script gives to `match_distances`, in list order.
    - `distances.tsv` has the columns `case`, `anchor` (`fivePrimeEnd`, `threePrimeEnd`, or `center`), `seq_id`, and `distance`. It has one row for each value of the R list, in list order and in the order of each vector. `NA` stays `NA`.
  - The cases of the distances:
    - `mer6`: the matches of `vmatchPattern2(<mer6 target>, seqs)`. Most names have one match, and some names have more than one.
    - `mer7m8_indel`: the matches of `vmatchPattern2(<mer7m8 target>, seqs, max.mismatch = 1, with.indels = TRUE)`. Their widths differ, thus `center` has half values.
    - `synthetic`: the names `a` (starts 10, 30, 35, widths 7, 7, 8), `b` (start 5, width 7), and `c` (starts 100, 50, width 7).
  - A note of the distances section: `match_distances` of a list element without a range. The planner probe gives `NA NA`.
  - **E-value.** Start with `set.seed(20261122)`. The data:
    - The genes are `g01` to `g40`. The observed counts are `sample(0:4, 40, replace = TRUE)`, and then `g01` gets 0.
    - The null is a 40 × 25 matrix of `sample(0:5, 1000, replace = TRUE)`. Its rows have the reverse order of the genes, thus a match by row position fails. `symbol` holds the gene IDs.
    - Write `evalue_observed.tsv` (`gene_id`, `overall_matches`) and `evalue_null.tsv` (`gene_id`, `shuffle` from 1 to 25, `count`).
    - Write `evalue.tsv` (`gene_id`, `overall_matches`, `e_value`) in the row order of the R result.
  - Notes of the e-value section: the R results for the inputs of the planner probe. Let `null` be `data.frame(symbol = c("a", "b"), r1 = c(1, 0), r2 = c(2, 1), r3 = c(3, 0), r4 = c(0, 0))`.
    - `generate_e_value(data.frame(symbol = "a", overall_matches = 1))` stops with `object 'random_siRNA_results_overall' not found`.
    - The observed counts `c(a = NA)` and `c(a = -1)` give the e-value 0. The second gives the warning `NaNs produced`.
    - The observed rows `a = 1` and `a = 3` both get 0.5. Their e-values are 0.75 and 0.25.
    - A null with two rows for `a` (the second row `c(5, 5, 5, 5)`) gives the e-value 1.5 for `a = 2`.
    - An observed table without a gene of the null stops with `replacement has 1 row, data has 0`.
  - **Curves.** The R calls are `PRROC::pr.curve(scores.class0 = score, weights.class0 = label, curve = TRUE)` and `PRROC::roc.curve(scores.class0 = score, weights.class0 = label, curve = TRUE)`. The files:
    - `curve_samples.tsv` has the columns `case`, `score`, and `label` (0 or 1).
    - `pr.tsv` has `case`, `auc_integral`, and `auc_davis_goadrich`. `pr_curve.tsv` has `case`, `row`, `recall`, `precision`, and `threshold`.
    - `roc.tsv` has `case` and `auc`. `roc_curve.tsv` has `case`, `row`, `fpr`, `tpr`, and `threshold`.
  - The cases of the curves start with `set.seed(20261123)`:
    - `ties_small`: 40 scores from `sample(0:3, 40, replace = TRUE)`, and labels from `rbinom(40, 1, 0.4)`.
    - `continuous`: 60 scores from `rnorm(60)`, and labels from `rbinom(60, 1, 0.5)`.
    - `large`: 400 scores from `sample(0:5, 400, replace = TRUE)`, and the label 1 with the probability `0.3 + 0.1 * score`. `stopifnot` that it has at least 100 positives, thus its curve uses the recall step 1/P.
    - `one_value`: the score 1 for each of 10 points, and the labels `rep(c(0, 1), 5)`.
    - `perfect`: the scores `5:1`, and the labels `c(1, 1, 0, 0, 0)`.
    - `first_negative`: the scores `c(9, 9, 5, 5, 1, 1)`, and the labels `c(0, 0, 1, 0, 1, 1)`. The top score has only negatives.
    - `inverse`: the scores `1:6`, and the labels `c(1, 1, 1, 0, 0, 0)`.
    - `stopifnot` that each case has at least one positive and one negative.
  - **Benchmark.** Run the patched `benchmark_preds` with `with_seedmatchr_patch` on `load_test_data()$res` and `preds_col = "mer7m8"`. The script counts the labels with the filter of the patch.
    - `benchmark.tsv` has the columns `case`, `column`, `padj`, `l2fc`, `n_positive`, `n_negative`, `auc_integral`, `auc_davis_goadrich`, and `roc_auc`.
    - The cases are `benchmark_1` (`padj = 0.15`, `l2fc = 0`), `benchmark_2` (`1e-10`, `0`), and `benchmark_3` (`1e-10`, `-1`).
    - Append the scores, the labels, and the curves of each case to `curve_samples.tsv`, `pr_curve.tsv`, and `roc_curve.tsv`, with the same case ID.
  - Notes of the benchmark section:
    - The unpatched `benchmark_preds(res, padj = 1e-10)`, with its PR `auc.integral`, its ROC `auc`, and the number of genes that its filter keeps. The planner probe gives 63 genes for the unpatched filter, and 24 for the patched filter.
    - The values that the patched function returns for a NaN score in one row, and for `padj = 0` (no positive).
  - **miRDB.** `mirdb.tsv` is the table of `R/get_data.R:43-52` for `tests/data/mirdb_sample.txt`, with the R column names. The notes give the R result for these lines in a file of `tempdir()`:
    - A score `abc`. The planner probe gives a text column.
    - An empty score. The planner probe gives `NA`.
    - The ID `NA`. The planner probe gives the value `NA` for the ID, and the species `NA`.
    - A line with 2 fields, and a field with `#`. The planner probe gives `line 2 did not have 3 elements` for both.
    - A field with `'`. The planner did not do a probe of this line.
- [ ] **Step 6: Run the script two times.** Expected: the second run of `tools/golden/run.sh tools/golden/R/analyses.R` changes no byte. For `mer7m8_1`, `alignment.tsv` has these lines of 37 characters, with their spaces:

```text
     3'-UUUUGUCACAAGAACGAGAUAUU-5'   
         ||||||||||||||||||||||      
  5'-...GAAACAGUGUUCUUGCUCUAUAA...-3'
```

  Also expected:
  - `benchmark.tsv` gives `n_positive` 63 for `benchmark_1` and 24 for `benchmark_2`.
  - In `mirdb.tsv`, the row of `lin4` has the species `lin4` and an empty `miRNA.name`.
- [ ] **Step 7: Document the goldens.** Add a section "Golden files of the small analyses" to `tools/golden/README.md`: the command, the PRROC pin, the patches, the cases, and the columns. Add the two fixtures to `tests/data/README.md`, and state that they are synthetic. Checkpoint.

### Task 2: The k-mer counts

**Files:**
- Make: `src/seedmatch/kmers.py`, `tests/test_kmers.py`
- Modify: `src/seedmatch/__init__.py`, `tests/helpers.py` (`run_core_pipeline`), `tests/test_readme.py`, `README.md`, `docs/r-map.md`, `docs/divergences.md` (`D-010`, `D-042`, `D-043`, `D-044`), `docs/design.md`

**Interfaces:**
- Consumes: `SequenceSet`, `as_limit` and `check_at_least` of `seedmatch._checks`, and the k-mer goldens of Task 1.
- Produces:

```python
MAX_K: Final = 12

@dataclass(frozen=True, eq=False)
class KmerCounts:
    k: int
    step: int
    ids: tuple[str, ...]
    gene_ids: tuple[str, ...]
    matrix: scipy.sparse.csr_array  # shape (len(ids), 4**k), dtype int64

    def kmers(self) -> tuple[str, ...]: ...
    def to_frame(self) -> pd.DataFrame: ...

def kmer_counts(sequences: SequenceSet, *, k: int = 7, step: int = 1) -> KmerCounts: ...
def kmer_gene_counts(counts: KmerCounts) -> pd.Series: ...
```

Inputs of `kmer_counts`. The defaults are the R defaults `width = 7` and `step = 1`:
- `sequences` is a `SequenceSet`. Another type raises `InputError`: `sequences has type list; expected a SequenceSet`.
- `k` goes through `as_limit`, thus a bool raises. A value outside 1 to `MAX_K` raises `k 13 is out of range; expected an integer from 1 to 12` (`D-043`).
- `step` goes through `as_limit` and `check_at_least(..., 1)`: `step 0 is less than 1; expected 1 or more`.

The counts obey `oligonucleotideFrequency` of Biostrings 2.78.0, as `sequence_kmer_counts` calls it:
- The windows of a sequence of length L start at 1, 1 + step, 1 + 2 · step, and so on. The last start is at most L − k + 1. A sequence shorter than k has no window, and its row holds only 0.
- A window counts only if each of its letters is A, C, G, or T. A window with another IUPAC letter, for example `N` or `R`, does not count.
- A planner probe shows the rule. For `ACGTNACGTA` and k = 2, Biostrings gives `AC` 2, `CG` 2, `GT` 2, and `TA` 1. With `step = 2`, it gives `AC`, `GT`, `CG`, and `TA` 1 each.
- Column j is the k-mer whose base-4 digits give j, with A = 0, C = 1, G = 2, and T = 3. The first letter is the most significant digit. This is the column order of Biostrings: `AA`, `AC`, `AG`, `AT`, `CA`, and so on.
- Row i is sequence i of the set, also if its ID occurs again (`D-044`). `ids` and `gene_ids` are the fields of the set.
- `matrix` is a `scipy.sparse.csr_array` of shape (n, 4^k) and dtype `int64`. It has the canonical format: sorted indices, no duplicate entry, and no stored 0.

The algorithm keeps the work memory bounded:
- A table of 16 entries maps each mask to a digit: A to 0, C to 1, G to 2, T to 3, and each other mask to 4.
- The function reads the start positions of the buffer in blocks of `_BLOCK_STARTS = 2**20` positions. A block can split a sequence. Block [a, b) reads the digits of `codes[a : b + k − 1]`, cut at the end of the buffer.
- k vectorized passes of `code = 4 * code + digit` give the code of each window. A cumulative count of the digit 4 marks each window with another letter.
- `np.searchsorted(offsets, start, side="right") - 1` gives the sequence of each start.
- A start counts if each of these conditions is true:
  - Its window ends in its sequence.
  - `start - offset` is a multiple of `step`.
  - Each letter of the window is A, C, G, or T.
- `np.unique(row * 4**k + code, return_counts=True)` counts the windows of the block.
- The counts of the blocks make the `csr_array` through the COO form. The COO form sums the entries of a sequence that two blocks split.
- The work memory of a block is about 40 bytes for each start, thus about 40 MB. The matrix holds at most one entry for each window.

Rules of `KmerCounts`:
- `kmers()` returns the label of each column, in column order. It makes 4^k strings: 16,384 for k = 7 and 16,777,216 for k = 12.
- `to_frame()` returns the layout of `sequence_kmer_counts` with `simplify.as = "matrix"`. The columns are `seq_id`, `gene_id`, and then one `int64` column for each k-mer, in column order. The frame has one row for each sequence and a `RangeIndex`.
- `to_frame()` is dense: it holds n · 4^k · 8 bytes, for example 3.0 GB for 23,000 sequences and k = 7. The docstring gives this cost, and it names `matrix` and `kmer_gene_counts` as the small forms.

Rules of `kmer_gene_counts` (`kmer_hist` without the plot, `R/kmer_hist.R:19-30`):
- `counts` is a `KmerCounts`. Another type raises `InputError`: `counts has type DataFrame; expected a KmerCounts`.
- The result is an `int64` Series named `genes`. Its index is named `kmer` and holds `counts.kmers()`. Each k-mer has a value, also a k-mer that no gene has, because the histogram of a later plan uses each value.
- The value of a k-mer is the number of distinct gene IDs with at least one sequence whose count is 1 or more (`D-042`).
- R gives the target count as `gene.counts[target.seq]`. The port gives it as `result[target.sequence]`, and the docstring shows this.
- The work is the sparse product of a gene indicator matrix with the presence matrix. Factorize the gene IDs one time. The memory is O(the entries of the matrix + 4^k).

- [ ] **Step 1: Write the failing tests** in `tests/test_kmers.py`:
  - `test_kmer_counts_match_r(case)`: parametrized over the cases of `kmer_counts.tsv`. Build the set from `Test_seqs.fa` or from `kmer_cases.tsv` with `SequenceSet.from_records`. Assert that the entries above 0 of `matrix`, as (row, k-mer, count) in row order and then in column order, equal the golden rows.
  - `test_to_frame_has_the_r_layout`: `cases_k2_s1` against `kmer_frame.tsv`, with `tx_id` as `seq_id`. Assert the column names and their order, the zero columns, and the `int64` dtype.
  - `test_kmer_counts_equal_an_oracle(sequences, k, step)`: a hypothesis test. The sequences have IUPAC letters, runs of `N`, empty sequences, and lowercase letters. k is from 1 to 4, and step is from 1 to 4.
  - In the oracle test, the oracle counts the window texts at each start in Python. Monkeypatch `seedmatch.kmers._BLOCK_STARTS` to 7, thus a block splits a sequence.
  - `test_kmer_counts_skip_a_window_with_an_ambiguous_letter`: `ACGTNACGTA` with k = 2 gives the counts of the planner probe.
  - `test_kmer_gene_counts_match_r_for_one_sequence_per_gene`: `test_data_k7_s1` against `kmer_gene_counts.tsv`, each of the 16,384 rows in order. The set gives each sequence its own ID as the gene ID.
  - `test_kmer_gene_counts_give_the_counts_of_kmer_hist_D010`: the value of `CTCTATA` equals the note of the golden. The docstring gives the R value `integer(0)` of `kmer_hist`.
  - `test_kmer_gene_counts_count_genes_not_sequences_D042`: for `cases_k2_s1`, the port value of each k-mer is the number of distinct genes in the rows of `kmer_counts.tsv`. The R value of at least one k-mer is greater, because `s1` and `s2` of `g1` both have it.
  - `test_kmer_counts_reject_k_above_12_D043`, `test_kmer_counts_reject_k_below_1`, `test_kmer_counts_reject_a_bool_k`, and `test_kmer_counts_reject_step_below_1`: `InputError` with a `match` on the template.
  - `test_kmer_counts_keep_each_sequence_with_a_repeated_id_D044`: two sequences with the ID `t1` give two rows with their own counts.
  - `test_kmer_counts_of_an_empty_set`: the shape is (0, 4^k), and `to_frame()` has the columns and no row.
- [ ] **Step 2: Run the tests.** Run `uv run pytest tests/test_kmers.py -q`. Expected: FAIL with `ModuleNotFoundError: No module named 'seedmatch.kmers'`.
- [ ] **Step 3: Write `kmers.py`.** Export `KmerCounts`, `kmer_counts`, and `kmer_gene_counts` from `seedmatch`. Add a `kmer_counts` call with k = 3 to `run_core_pipeline` in `tests/helpers.py`.
- [ ] **Step 4: Run the tests.** Expected: PASS. Also run `uv run ruff check`, `uv run ruff format --check`, and `uv run mypy`.
- [ ] **Step 5: Write the README part and the documents.**
  - `README.md`: add the section "## Other analyses" after the `Reference` example. Its first Python block calls `kmer_gene_counts(kmer_counts(sequences))` and prints the value of `target.sequence`.
  - `tests/test_readme.py`: add `test_readme_kmer_example_gives_the_r_target_count`.
  - `README.md`: add the k-mer counts to the first list of "Scope". Add a bullet for `D-042` to "Differences".
  - `docs/r-map.md`, row `sequence_kmer_counts` → `seedmatch.kmer_counts` and `KmerCounts.to_frame`:
    - `seqs` becomes `sequences`, and `width` becomes `k`.
    - `tx2gene` has no counterpart, because each sequence of a `SequenceSet` has its gene ID. R gives `NA` to a name that `tx2gene` lacks, and the IDs `1`, `2`, and so on to a set without names.
    - `out.path` and `path.basename` have no counterpart, because the core writes no file. Use `to_frame().to_csv(path, sep="\t")`.
    - `to_frame()` has `seq_id` for `tx_id`, and `int64` for the integers of R.
    - `simplify.as = "list"` and `"collapsed"` have no counterpart. Use a row of `matrix`, or `matrix.sum(axis=0)`.
    - R stops with `'width' must be an integer >= 1` and `'step' must be a positive integer`. Name `D-043` and `D-044`.
  - `docs/r-map.md`, row `kmer_hist` → `seedmatch.kmer_gene_counts`:
    - `seqs`, `width`, and `step` become the arguments of `kmer_counts`.
    - `target.seq` has no argument, because `result[target.sequence]` gives the target count.
    - `title` and the plot are in a later plan. Name `D-010` and `D-042`.
  - `D-010`: "kmer_gene_counts returns the counts of kmer_hist."
    - R location: `kmer_hist` (`R/kmer_hist.R:32-36`).
    - R behavior: `plot(hist(...))` returns `NULL`, and `NULL + abline(...)` returns `integer(0)`. Evidence: the note of `kmer_gene_counts.tsv`.
    - New behavior: `kmer_gene_counts` returns the gene count of each k-mer. The plot module of a later plan draws the histogram.
    - Reason: R returns no data, thus a caller cannot read the counts or the target count.
  - `D-042`: "kmer_gene_counts counts genes, not sequences."
    - R location: `kmer_hist` (`R/kmer_hist.R:24-28`).
    - R behavior: `colSums` of the presence matrix counts the sequences. The axis label is `Number of genes with k-mer` (`:32`), and the title is `K-mer gene frequencies`. Evidence: `cases_k2_s1` of `kmer_gene_counts.tsv`.
    - New behavior: the count of the distinct gene IDs. A set with one sequence for each gene gives the R values.
    - Reason: the label of R names genes. A gene-level result uses all the sequences of a gene (`D-003`).
  - `D-043`: "kmer_counts rejects a k above 12."
    - R location: `sequence_kmer_counts` calls `oligonucleotideFrequency` (`R/kmer_hist.R:63-66`).
    - R behavior: R accepts `width = 13`. It makes a dense integer matrix with 4^13 = 67,108,864 columns for each sequence. Evidence: the note of the k-mer section.
    - New behavior: `InputError` for a `k` above 12, and a sparse matrix.
    - Reason: the labels of `kmers()`, `to_frame()`, and `kmer_gene_counts` grow as 4^k. 4^12 is 16,777,216 labels, about 1 GB.
  - `D-044`: "kmer_counts gives a row for each sequence, also for a repeated ID."
    - R location: `sequence_kmer_counts` (`R/kmer_hist.R:71`).
    - R behavior: `rownames(olig.freq) = names(seqs)` stops with `duplicate 'row.names' are not allowed`. Evidence: the note of the k-mer section.
    - New behavior: one row for each sequence, in the order of the set.
    - Reason: a row of the matrix is a sequence, thus a repeated ID does not make a count ambiguous.
  - `docs/design.md`: the table rows of `D-042`, `D-043`, and `D-044`. Add to "Performance and bounds":
    - `kmer_counts` works in blocks of 2^20 window starts, and its matrix has at most one entry for each window.
    - `to_frame` is dense and holds n · 4^k · 8 bytes.
  - Checkpoint.

### Task 3: The alignment text

**Files:**
- Make: `src/seedmatch/alignment.py`, `tests/test_alignment.py`
- Modify: `src/seedmatch/seed.py` (`seed_window`), `tests/test_seed.py`, `src/seedmatch/__init__.py`, `tests/helpers.py` (`run_core_pipeline`), `tests/test_readme.py`, `README.md`, `docs/r-map.md`, `docs/divergences.md` (`D-045`, `D-046`), `docs/design.md`

**Interfaces:**
- Consumes: `Guide`, `Seed`, `CustomSeed`, `SeedSpec`, `SequenceSet`, `as_limit` and `check_not_negative` of `seedmatch._checks`, and `alignment.tsv` of Task 1.
- Produces: `seed_window` of "Shared interfaces", and:

```python
def alignment_text(
    match: pd.Series | Mapping[str, object],
    guide: Guide,
    sequences: SequenceSet,
    *,
    seed: SeedSpec = Seed.MER7M8,
    left_pad: int = 0,
    right_pad: int = 0,
) -> str: ...
```

Inputs of `alignment_text` (`plotAlignment` with the fix of `D-045`):
- `match` is a row of `find_matches`: a `Series` or a mapping with the keys `seq_id`, `start`, and `end`. Other keys can be present.
- An absent key raises `InputError`: `match has no 'start'; expected a row of find_matches with seq_id, start, and end`.
- `start` and `end` go through `as_limit`. `start` can be below 1, as in a mismatch match.
- Exactly one sequence of `sequences` has the ID `seq_id`. No sequence raises `seq_id 'g9' is not an ID of the sequences; expected the ID of one sequence`. More than one raises `seq_id 'g1' names 2 sequences; expected the ID of one sequence`.
- `guide.target(seed)` makes sure that the seed fits the guide. Let N be the length of the guide, and (a, b) = `seed_window(seed, N)`.
- The width `end − start + 1` must be `b − a + 1`. Another width raises `match 'g1' 73-80 has width 8, but the window of mer7m8 has width 7; expected a match without an insertion or a deletion`. R stops with `nchar(s1) == nchar(s2) is not TRUE`.
- `left_pad` and `right_pad` go through `as_limit` and `check_not_negative`: `left_pad -1 is negative; expected 0 or more`.

The text obeys `.matchStrings` and `.create_alignment_df`:
- The window goes from `start − (N − b) − left_pad` to `end + (a − 1) + right_pad`. Its width is W = N + `left_pad` + `right_pad`.
- The window text has the letter of the sequence at each position from 1 to the length of the sequence. It has a space at each other position (`D-046`).
- Let `r` be the guide in reverse order (3' to 5'). The three lines are:
  1. `" " * (5 + left_pad) + "3'-" + r + "-5'" + " " * (3 + right_pad)`
  2. `" " * 8 + pairs + " " * 6`
  3. `"  5'-..." + window.replace("T", "U") + "...-3'"`
- `pairs[i]` is `|` if `left_pad <= i < left_pad + N`, and the window letter at i is the DNA complement of `r[i - left_pad]`. The complement of A, C, G, and U is T, G, C, and A.
- Each other position of `pairs` is a space. Thus an IUPAC letter, a space, and a G:U pair give a space, as `.matchStrings` does.
- The result is the three lines joined with `"\n"`, without a final line end. Each line has W + 14 characters, as the left fill of `.create_alignment_df` makes them. The row names of R are not in the text.
- The time is O(W), and the scan of the IDs is O(n). The docstring states that a loop over m matches costs O(m · n).

- [ ] **Step 1: Write the failing tests** in `tests/test_alignment.py` and `tests/test_seed.py`:
  - `test_alignment_text_matches_r(case)`: parametrized over the cases of `alignment.tsv`. Map `seed` to a `Seed`, or to `CustomSeed("custom", start_pos, stop_pos)`. Build the match from `seq_id`, `start`, and `end`. Assert that the three lines equal the golden texts, with each space.
  - `test_alignment_text_shows_the_full_guide_D045`: the guide line of `mer7m8_1` has the 23 letters of the guide in reverse order. The docstring gives the error of the unpatched R from the golden note.
  - `test_alignment_text_puts_spaces_outside_the_sequence_D046`: a synthetic set holds the `mer7m8` target at position 3 and at the end. The target line and the pair line have a space for each position outside the sequence, and each line has W + 14 characters.
  - `test_alignment_text_accepts_a_row_of_find_matches`: `find_matches(target, seqs).iloc[0]` and the same values in a `dict` give the text of `mer7m8_1`.
  - `test_alignment_text_accepts_a_mismatch_match_before_the_start_D046`: a row of `find_matches` with `start` below 1 gives a text.
  - `test_alignment_text_rejects_an_indel_match`, `test_alignment_text_rejects_an_unknown_seq_id`, `test_alignment_text_rejects_a_repeated_seq_id`, `test_alignment_text_rejects_a_negative_pad`, and `test_alignment_text_rejects_a_match_without_start`: `InputError` with a `match` on the template.
  - In `tests/test_seed.py`: `test_seed_window_gives_g1_to_the_a_of_mer8_and_mer7a1`. Assert `seed_window(Seed.MER8, 23) == (1, 8)` and `seed_window(Seed.MER7A1, 23) == (1, 7)`. The `mer8_<i>` goldens pass only with (1, 8).
- [ ] **Step 2: Run the tests.** Run `uv run pytest tests/test_alignment.py tests/test_seed.py -q`. Expected: FAIL with `ModuleNotFoundError: No module named 'seedmatch.alignment'`.
- [ ] **Step 3: Write `seed_window` and `alignment.py`.** Export `alignment_text` from `seedmatch`. Add an `alignment_text` call to `run_core_pipeline`.
- [ ] **Step 4: Run the tests.** Expected: PASS, and ruff and mypy pass.
- [ ] **Step 5: Write the README part and the documents.**
  - `README.md`: in "Other analyses", a Python block prints `alignment_text(find_matches(target, sequences).iloc[0], guide, sequences)`. Add the alignment text to "Scope", and a bullet for `D-045` to "Differences".
  - `tests/test_readme.py`: add `test_readme_alignment_example_gives_the_r_text`. The text equals `mer7m8_1`.
  - `docs/r-map.md`, row `plotAlignment` with `output.type = "data.frame"` → `seedmatch.alignment_text`:
    - `match.info` becomes `match`, a row of `find_matches`. R reads the column `names`, and the port reads `seq_id`.
    - `guide.seq` becomes a `Guide`, and `seqs` becomes `sequences`.
    - `seq.start.pos` and `seq.stop.pos` become `seed`, a `Seed` or a `CustomSeed(name, start, stop)`. `seed_window` gives the positions.
    - The `A` of `MER8` and `MER7A1` takes g1, thus an R call for `mer8` gives `seq.start.pos = 1`.
    - `left.pad` and `right.pad` become `left_pad` and `right_pad`. `aln.name` has no counterpart.
    - The function returns the three values of the column `alignment` as one `str`, joined with line ends. `output.type = "ggplot"` is in a later plan.
    - Name `D-045` and `D-046`, and the `InputError` for a match with an indel.
  - `docs/r-map.md`, row `.matchStrings`, `.make_space`, and `.create_alignment_df` → `seedmatch.alignment_text`. These helpers make the lines. The `gsub` calls of `.matchStrings` (`R/plotAlignment.R:89-90`) change no line, because the line holds only `|` and spaces.
  - `D-045`: "alignment_text shows the full guide."
    - R location: `plotAlignment` (`R/plotAlignment.R:28-57`), `get_seed` (`R/get_seed.R:121`), and `.matchStrings` (`R/plotAlignment.R:84`).
    - R behavior: `Full` is a canonical seed, thus `get_seed(guide.seq, "Full", start.pos = 1, stop.pos = guide.len)` ignores the positions. Its target is g2 to gN, with N − 1 letters, but the window has N letters.
    - Thus each call with the default positions stops with `nchar(s1) == nchar(s2) is not TRUE`. With `seq.start.pos = 1`, the call works, but the guide line omits g1.
    - A row of `SeedMatchR(res.format = "data.frame")` has no column `names`, thus R stops with `subscript out of bounds`. Evidence: the notes of `alignment.tsv`.
    - New behavior: the guide line holds g1 to gN, and `match` gives the sequence by `seq_id`.
    - Reason: the R documentation calls `guide.seq` "The full siRNA guide sequence". `tests/testthat/helper.R:27-45` of SeedMatchR uses strings of 23 letters for the guide of 23 letters. The function must work for its documented input.
  - `D-046`: "alignment_text shows a window past a sequence end."
    - R location: `plotAlignment` (`R/plotAlignment.R:43-45`, `XVector::subseq`).
    - R behavior: a window that starts before 1 stops with `solving row 1: the supplied start/end lead to a negative width`.
    - A window that ends after the sequence stops with `solving row 1: 'allow.nonnarrowing' is FALSE and the supplied end (147) is > refwidth`. Evidence: the notes of `alignment.tsv`.
    - New behavior: a space in the target line and in the pair line for each position outside the sequence.
    - Reason: a seed match in the first N − b bases of a 3' UTR has no alignment in R. `find_matches` also gives mismatch matches that start before 1 (search rule 2).
  - `docs/design.md`: the table rows of `D-045` and `D-046`. Add `alignment_text` to "Performance and bounds". Checkpoint.

### Task 4: The match distances

**Files:**
- Make: `src/seedmatch/distances.py`, `tests/test_distances.py`
- Modify: `src/seedmatch/_checks.py` (`as_data_frame`, `check_columns`), `tests/test_checks.py`, `src/seedmatch/__init__.py`, `tests/helpers.py` (`run_core_pipeline`), `tests/test_readme.py`, `README.md`, `docs/r-map.md`, `docs/divergences.md` (`D-011`), `docs/design.md`

**Interfaces:**
- Consumes: `as_member` of `seedmatch._checks`, and `distances_ranges.tsv` and `distances.tsv` of Task 1.
- Produces: the new keywords of `as_data_frame` and `check_columns` ("Shared interfaces"). The messages become `f"{name} has type {type}; expected a pandas.DataFrame"` and `f"{label} has missing columns: {missing}; expected {required}"`. The defaults give the present messages. Also:

```python
class MatchAnchor(StrEnum):
    FIVE_PRIME_END = "fivePrimeEnd"
    THREE_PRIME_END = "threePrimeEnd"
    CENTER = "center"

def match_distances(
    matches: pd.DataFrame, *, anchor: MatchAnchor = MatchAnchor.FIVE_PRIME_END
) -> pd.DataFrame: ...
```

Inputs of `match_distances` (`R/match_distances.R`):
- `matches` is a `DataFrame` with the columns `seq_id`, `gene_id`, `start`, `end`, and `width`, for example the result of `find_matches`. Other columns can be present.
- An absent column raises `matches table has missing columns: ['width']; expected ['seq_id', 'gene_id', 'start', 'end', 'width']`.
- `start`, `end`, and `width` have an integer dtype, and bool is not an integer: `column 'start' has dtype float64; expected an integer dtype`.
- `anchor` goes through `as_member`. The default is the R default `fivePrimeEnd`.

Rules of the result:
- The coordinate of a match is `end` for `FIVE_PRIME_END`, `start` for `THREE_PRIME_END`, and `start + width / 2` for `CENTER`. These are the coordinates of `.get_fivePrimeEnd`, `.get_threePrimeEnd`, and `.get_centers`.
- The rows with one `seq_id` are one sequence. The sequences come in the order of their first row. The matches of a sequence keep the order of the rows.
- With the coordinates c_1 to c_n, the pairs come in the R order: (1, 2), (1, 3), ..., (1, n), (2, 3), and so on. The distance of a pair is |c_i − c_j|.
- A sequence with n matches gives n(n − 1)/2 rows. A sequence with one match gives no row (`D-011`).
- The result has the columns `seq_id`, `gene_id` (of the first row of the sequence), and `distance`, and a `RangeIndex`. `distance` is `float64`, because `CENTER` gives half values.
- An empty `matches` gives an empty result with these columns.
- Count the pairs of each sequence first, and allocate the result one time. Fill the pairs of match i with one vectorized difference. Thus the work memory of a step is O(n).
- No bound below the size of the result is possible: 10,000 matches in one sequence give 49,995,000 rows. The docstring states this.

- [ ] **Step 1: Write the failing tests** in `tests/test_distances.py`, and add to `tests/test_checks.py`:
  - `test_match_distances_match_r(case, anchor)`: for each case of `distances_ranges.tsv` and each anchor, build `matches` from the golden ranges. Compare the rows with `distances.tsv`, without the sequences that have one match.
  - `test_a_sequence_with_one_match_gives_no_distance_D011`: the golden has `NA` and 0 for each such sequence of `mer6`, and the port has no row.
  - `test_center_is_the_start_plus_half_the_width`: the `synthetic` case gives 20, 25.5, and 5.5 for `a`.
  - `test_each_pair_occurs_once(starts)`: a hypothesis test against an oracle with `itertools.combinations`, in the R order.
  - `test_rows_keep_the_order_of_the_sequences`: the result of `find_matches` with `Seed.MER6` keeps the order of the set.
  - `test_match_distances_rejects_a_missing_column`, `test_match_distances_rejects_a_float_start`, and `test_match_distances_rejects_an_anchor_as_text`: `InputError` with a `match` on the template.
  - `test_an_empty_table_gives_an_empty_result`.
  - In `tests/test_checks.py`: `test_check_columns_names_the_table`, with `label="matches table"`.
- [ ] **Step 2: Run the tests.** Run `uv run pytest tests/test_distances.py tests/test_checks.py -q`. Expected: FAIL with `ModuleNotFoundError: No module named 'seedmatch.distances'`.
- [ ] **Step 3: Write `distances.py` and extend `_checks.py`.** Export `MatchAnchor` and `match_distances` from `seedmatch`. Add a `match_distances` call to `run_core_pipeline`.
- [ ] **Step 4: Run the tests.** Expected: PASS, and ruff and mypy pass. Also run `uv run pytest tests/report tests/stats tests/test_de.py -q`, because they use the messages of `check_columns`.
- [ ] **Step 5: Write the README part and the documents.**
  - `README.md`: in "Other analyses", a Python block gives `match_distances(find_matches(guide.target(Seed.MER6), sequences))`. Add the distances to "Scope".
  - `tests/test_readme.py`: add `test_readme_distance_example_runs`. The result has the three columns and at least one row.
  - `docs/r-map.md`, row `match_distances` (with `.get_centers`, `.get_fivePrimeEnd`, and `.get_threePrimeEnd`) → `seedmatch.match_distances`:
    - `gr_list` becomes `matches`, a `find_matches` table. The rows of one `seq_id` replace an element of the list.
    - `feature.pos` becomes `anchor`, a `MatchAnchor`.
    - The function returns a `DataFrame` with `seq_id`, `gene_id`, and `distance`, not a list of vectors.
    - A list without names stops in R with `attempt to select less than one element in OneIndex`. Name `D-011`.
  - `D-011`: "A sequence with one match gives no distance."
    - R location: `match_distances` (`R/match_distances.R:33-44`).
    - R behavior: for one match, `1:(length(feature_coordinates) - 1)` is `c(1, 0)`, thus the distances are `c(NA, 0)`. An element without a range gives `c(NA, NA)`.
    - For a name that occurs again, the list keeps the distances of the last element (`:44`). Evidence: `distances.tsv`, its note, and a planner probe.
    - New behavior: a sequence with n matches gives each of the n(n − 1)/2 pair distances, and one match gives no row. The rows of one `seq_id` are one sequence.
    - Reason: a sequence with one match has no pair. `NA` and the self-distance 0 change each summary of the distances.
  - `docs/design.md`: the row of `D-011` says "sequence", not "gene", because a coordinate belongs to one sequence. Add the bound of `match_distances` to "Performance and bounds". Checkpoint.

### Task 5: The e-value and its null

**Files:**
- Make: `src/seedmatch/evalue.py`, `src/seedmatch/_pool.py`, `tests/test_evalue.py`
- Modify: `src/seedmatch/_checks.py` (`as_generator`), `src/seedmatch/stats/dts.py`, `src/seedmatch/report/grid.py`, `src/seedmatch/__init__.py`, `tests/helpers.py` (`run_core_pipeline`), `tests/test_readme.py`, `README.md`, `docs/r-map.md`, `docs/divergences.md` (`D-008`, `D-047`), `docs/design.md`

**Interfaces:**
- Consumes: `Guide`, `SeedSpec`, `SearchOptions`, `count_matches`, `SequenceSet`, `check_columns`, `check_gene_ids_are_text`, `check_numeric_columns`, and the e-value goldens of Task 1.
- Produces: `as_generator`, `pool_width`, and `map_in_pool` of "Shared interfaces". `dts.py` and `grid.py` use them, and their behavior does not change. Also:

```python
def e_value_null(
    guide: Guide,
    sequences: SequenceSet,
    *,
    rng: np.random.Generator | int,
    n_random: int = 1000,
    seed: SeedSpec = Seed.MER7M8,
    options: SearchOptions = SearchOptions(indels=True),
    workers: int | None = None,
) -> pd.DataFrame: ...

def e_values(de: pd.DataFrame, null: pd.DataFrame, *, column: str) -> pd.DataFrame: ...
```

Inputs of `e_value_null` (the null of `D-008`):
- `rng` goes through `as_generator`, and `workers` goes through `pool_width`.
- `n_random` goes through `as_limit` and `check_at_least(..., 1)`: `n_random 0 is less than 1; expected 1 or more`.
- `options` is a `SearchOptions`: `options has type dict; expected a SearchOptions`.
- `guide.target(seed)` makes sure that the seed fits the guide, before any draw.

Rules of the null:
- Shuffle j, for j = 0 to `n_random` − 1, takes the letters of the guide in the order `generator.permutation(len(guide.sequence))`. Thus a shuffle keeps the length and the letters of the guide.
- The function draws each shuffle in this order before the first search. Thus the null does not depend on `workers`.
- The count of shuffle j for a gene is the maximum of `count_matches(Guide(shuffle_j).target(seed), sequences, options)` over the sequences of the gene. This is the gene rule of `annotate_seed_counts` (`D-003`), and its default `options` are the default of `annotate_seed_counts`.
- The result is a `DataFrame` with an index named `gene_id`. The index holds the distinct gene IDs of `sequences`, in the order of their first sequence.
- The result has one `int32` column for each shuffle, with the labels 0 to `n_random` − 1. A count of matches in one gene fits in `int32`.
- Factorize the gene IDs one time, before the searches. `map_in_pool` runs the searches, and each result goes into a preallocated array.
- The time is `n_random` searches of the set on `workers` threads. The memory is genes × `n_random` × 4 bytes, for example 80 MB for 20,000 genes and 1000 shuffles.

Inputs of `e_values` (`generate_e_value`, `R/generate_e_value.R`):
- `de` has the columns `gene_id` (text) and `column` (numeric).
- Each value of `column` is a count of 0 or more (`D-047`): `column 'mer7m8' of de holds NaN for gene 'g4'; expected a count of 0 or more`. A negative count gets the same text, with its value for `NaN`.
- `null` is a `DataFrame`: `null has type dict; expected a pandas.DataFrame with one row for each gene`.
- The index of `null` holds text gene IDs, and each ID occurs once (`D-047`): `null has 2 rows for gene 'a'; expected one row for each gene`.
- `null` has at least one column, and each value is a count of 0 or more: `null holds NaN for gene 'a'; expected counts of 0 or more`.

Rules of the e-values:
- The result has the rows of `de` whose `gene_id` is in the index of `null`, in the order of `de`, with a new index. Genes match by `gene_id` (`D-008`).
- A `de` without such a row raises `de has no gene of the null; expected the gene IDs of the null in the column gene_id, for example 'g01'`. R stops with `replacement has 1 row, data has 0`.
- The result is a copy of these rows with the `float64` column `e_value`. A column `e_value` that is present gets the new values.
- The e-value of a row is a fraction. Its numerator is the number of null columns whose count for the gene is at least the count of the row. Its denominator is the number of null columns.
- Each row has its own e-value, also for a gene with two rows (`D-047`). A count of 0 gives the e-value 1, as in R.
- Compare in blocks of 4096 rows, with `null.index.get_indexer` one time before the blocks. Thus the work memory is O(4096 × `n_random`), and no copy of the null is made.

- [ ] **Step 1: Write the failing tests** in `tests/test_evalue.py`:
  - `test_e_values_match_r`: build `null` from `evalue_null.tsv` and `de` from `evalue_observed.tsv`. `e_values(de, null, column="overall_matches")` gives the e-values of `evalue.tsv`, with exact equality, in the R row order.
  - `test_e_values_match_genes_by_gene_id_D008`: the null rows in another order give the same result.
  - `test_e_value_null_counts_each_shuffle_like_annotate_seed_counts_D008`: `n_random=5` and `rng=7`. The oracle draws the shuffles with the rule above. For each shuffle, `annotate_seed_counts(..., shared_genes=False)` on a table of each gene gives the column of the null.
  - `test_each_shuffle_keeps_the_letters_of_the_guide_D008`: each shuffle of the private shuffle function has the sorted letters of the guide.
  - `test_e_value_null_is_the_same_for_the_same_seed` and `test_e_value_null_does_not_depend_on_workers`, with `workers=1` and `workers=4`.
  - `test_e_values_reject_a_nan_count_D047`, `test_e_values_reject_a_negative_count_D047`, and `test_e_values_reject_a_repeated_gene_of_the_null_D047`: `InputError` with a `match` on the template. The docstring of each test gives the R value of the golden note.
  - `test_each_row_of_a_gene_gets_its_own_e_value_D047`: the rows `a = 1` and `a = 3` of the note get 0.75 and 0.25. R gives 0.5 to both.
  - `test_e_values_reject_a_table_without_a_gene_of_the_null`, `test_e_value_null_rejects_n_random_below_1`, and `test_e_value_null_rejects_options_of_another_type`.
- [ ] **Step 2: Run the tests.** Run `uv run pytest tests/test_evalue.py -q`. Expected: FAIL with `ModuleNotFoundError: No module named 'seedmatch.evalue'`.
- [ ] **Step 3: Move the helpers, and write `evalue.py`.**
  - Move `_as_generator` to `seedmatch._checks.as_generator`, and `_pool_width` and `_run` to `seedmatch/_pool.py`.
  - Run `uv run pytest tests/stats/test_dts.py tests/report/test_full_search.py -q`. Expected: PASS without a change of a test.
  - Write `evalue.py`. Export `e_value_null` and `e_values` from `seedmatch`. Add an `e_value_null` call with `n_random=2` to `run_core_pipeline`.
- [ ] **Step 4: Run the tests.** Expected: PASS, and ruff and mypy pass.
- [ ] **Step 5: Write the README part and the documents.**
  - `README.md`: in "Other analyses", a Python block gives `e_values(de, e_value_null(guide, sequences, rng=1, n_random=100), column="mer7m8")`. Add the e-value to "Scope".
  - `tests/test_readme.py`: add `test_readme_e_value_example_runs`. Each e-value is from 0 to 1.
  - `docs/r-map.md`, row `generate_e_value` → `seedmatch.e_values` and `seedmatch.e_value_null`:
    - `result_output` becomes `de`, and `random_seq_results` becomes `null`, a required argument.
    - The column `overall_matches` becomes the argument `column`.
    - The null is a `DataFrame` with the gene IDs as its index and one column for each random guide. The R null has a `symbol` column.
    - R compares the logarithms of the counts, which gives the same result as the counts. Name `D-008` and `D-047`.
  - `D-008`: "e_values has no default null, and genes match by gene_id." Write the entry from the row of `docs/design.md`.
    - R location: `generate_e_value` (`R/generate_e_value.R:12-24`).
    - R behavior: the default `random_siRNA_results_overall` does not exist. Thus a call without a null stops with `object 'random_siRNA_results_overall' not found`. SeedMatchR has no function that makes the null.
    - Genes match by `symbol` (`:14`, `:18-20`). The e-value is the fraction of the null columns with `log(null) >= log(observed)`, and `log(0)` is `-Inf`. Evidence: the notes of `evalue.tsv`.
    - New behavior: the null is a required argument. `e_value_null` makes it from `n_random` shuffles of the guide (default 1000), with the random generator of the caller. Genes match by `gene_id`.
    - Reason: without a null, R stops for each call. A symbol can name two genes, and a DE table can have no symbol.
  - `D-047`: "e_values rejects a count that is not valid, and gives each row its own e-value."
    - R location: `generate_e_value` (`R/generate_e_value.R:16-21`).
    - R behavior: a count of `NA` gives the e-value 0. A negative count gives 0 and the warning `NaNs produced`.
    - A gene with two rows in the null gives an e-value above 1, for example 1.5.
    - Two rows of one gene in `result_output` both get the value of a recycled comparison, for example 0.5 for 0.75 and 0.25. Evidence: the notes of `evalue.tsv`.
    - New behavior: `InputError` for a NaN or negative count, and for a gene with two rows in the null. Each row gets its own e-value.
    - Reason: the e-value 0 is the strongest result, thus a count that is not known must not give it. A fraction above 1 has no meaning.
  - `docs/design.md`: the table row of `D-047`. Add the bounds of `e_value_null` and `e_values` to "Performance and bounds". Checkpoint.

### Task 6: The PR and ROC curves

**Files:**
- Make: `src/seedmatch/stats/curves.py`, `tests/stats/test_curves.py`
- Modify: `src/seedmatch/stats/__init__.py`, `src/seedmatch/__init__.py`, `README.md` ("Scope"), `docs/r-map.md`, `docs/design.md`

**Interfaces:**
- Consumes: `as_sample` of `seedmatch.stats.types`, and the curve goldens of Task 1.
- Produces:

```python
@dataclass(frozen=True)
class PrCurve:
    auc_integral: float
    auc_davis_goadrich: float
    curve: pd.DataFrame  # recall, precision, threshold (float64)
    n_positive: int
    n_negative: int

@dataclass(frozen=True)
class RocCurve:
    auc: float
    curve: pd.DataFrame  # fpr, tpr, threshold (float64)
    n_positive: int
    n_negative: int

def pr_curve(scores: ArrayLike, labels: ArrayLike) -> PrCurve: ...
def roc_curve(scores: ArrayLike, labels: ArrayLike) -> RocCurve: ...
```

Inputs. These rules replace the call `pr.curve(scores.class0 = scores, weights.class0 = labels, curve = TRUE)` of PRROC 1.4, and the same call of `roc.curve`:
- `scores` goes through `as_sample("scores", scores, finite=False)`. Thus NaN raises `scores[2] is NaN; expected a number`, and an infinite score is valid.
- `labels` is 1-D, with the length of `scores`: `scores has 5 values and labels has 4; expected the same number`.
- Each label is a bool, or a number equal to 0 or 1: `labels[3] is 2; expected 0, 1, False, or True`. A label of 1 is a positive. PRROC calls the positives class 0.
- `pr_curve` raises `labels have no positive; expected at least one label 1 or True`.
- `roc_curve` raises the same error, and also `labels have no negative; expected at least one label 0 or False`.

Definitions:
- Let s_1 > s_2 > ... > s_K be the distinct scores. TP_k and FP_k count the positives and the negatives with a score of at least s_k, and TP_0 = FP_0 = 0.
- Let P = TP_K and N = FP_K. The recall is r_k = TP_k / P, and the precision is p_k = TP_k / (TP_k + FP_k).
- If TP_k > TP_{k−1}, then h_k = (FP_k − FP_{k−1}) / (TP_k − TP_{k−1}). Between the points k − 1 and k, FP grows on a line: FP = FP_{k−1} + h_k · (TP − TP_{k−1}).

`auc_integral` (Keilwagen, Grosse, and Grau 2014):
- On step k, the precision at the recall r is r / (a_k · r + b_k), with a_k = 1 + h_k and b_k = (FP_{k−1} − h_k · TP_{k−1}) / P.
- The area of step k is the integral of this precision from r_{k−1} to r_k: v_k = ((r_k − r_{k−1}) − (b_k / a_k) · (ln(a_k · r_k + b_k) − ln(a_k · r_{k−1} + b_k))) / a_k.
- If b_k = 0, then v_k = (r_k − r_{k−1}) / a_k. If TP_k = TP_{k−1}, then v_k = 0.
- `auc_integral` is the sum of v_k.

`auc_davis_goadrich` (Davis and Goadrich 2006):
- Step 1 adds (TP_1 / P) · p_1. PRROC uses the precision of point 1 for the full width of step 1.
- A step k ≥ 2 with TP_k − TP_{k−1} ≥ 2 and p_k ≠ p_{k−1} gets a point at each integer TP. Let Δ = TP_k − TP_{k−1} and q(x) = (TP_{k−1} + x) / (TP_{k−1} + x + FP_{k−1} + h_k · x).
- Such a step adds (1 / P) · (q(x) + q(x + 1)) / 2 for x = 0 to Δ − 1.
- Each other step adds the trapezoid ((TP_k − TP_{k−1}) / P) · (p_{k−1} + p_k) / 2.
- `auc_davis_goadrich` is the sum. The labels are hard, thus this area is always defined. PRROC gives `NA` for soft labels, which the port does not take.

The PR curve (`curve`, with the columns `recall`, `precision`, and `threshold`):
- The base rows are (r_k, p_k, s_k) for k = 1 to K.
- The recall step is σ = min(1, P / 100) / P. It comes from the PRROC default `minStepSize = min(1, P / 100)`.
- Each step k ≥ 2 with (TP_k − TP_{k−1}) / P > σ and p_k ≠ p_{k−1} gets intermediate rows. Let d = r_k − r_{k−1} and m = `int(d / σ + 1e-10)`.
- For j = 0 to m, let x = min(j · σ, d). The row is (r_{k−1} + x, (r_{k−1} + x) / (r_{k−1} + x + FP_{k−1} / P + h_k · x), s_k). This is `seq(0, d, by = σ)` of R 4.5.3.
- The intermediate rows come after the base rows, in the order of k and then of j.
- Sort the rows by `recall` from the highest to the lowest value, then by `threshold` from the lowest to the highest value. The sort is stable, thus equal keys keep their order. A planner probe shows that `order(..., decreasing = TRUE)` of R is stable.
- Then add a first row with the recall 1 and the precision and the threshold of the first sorted row. Add a last row with the recall 0 and the precision and the threshold of the last sorted row.
- Compute each value with the IEEE operations of these formulas, in their order. Then the recalls that R gives as equal stay equal, and the stable sort gives the R order.
- The curve has at most 2K + 2 + max(100, P) rows, thus O(n).

The ROC curve:
- The points are (0, 0) and (FP_k / N, TP_k / P) for k = 1 to K.
- `auc` is the sum over k of (fpr_k − fpr_{k−1}) · (tpr_k + tpr_{k−1}) / 2.
- `curve` has the columns `fpr`, `tpr`, and `threshold`, with one row for each point.
- The threshold column is the threshold column of PRROC. The point (0, 0) has s_1, the point k has s_{k+1} for k < K, and the point K has s_K.
- Thus the threshold of a row is the score of the points that the next row adds. The docstring states this rule.

The time of each curve is O(n log n) for one sort of the scores.

- [ ] **Step 1: Write the failing tests** in `tests/stats/test_curves.py`:
  - `test_pr_curve_matches_r(case)` and `test_roc_curve_matches_r(case)`: parametrized over the cases of `curve_samples.tsv` whose ID does not start with `benchmark_`. Assert the areas, the number of rows, and each value of each row, with the tolerances of the Global Constraints.
  - `test_roc_auc_is_the_probability_that_a_positive_scores_higher(samples)`: a hypothesis test with ties. The oracle is (the pairs with a higher positive score + 0.5 · the pairs with equal scores) / (P · N).
  - `test_pr_auc_integral_equals_a_numeric_integral(samples)`: a hypothesis test. The oracle integrates the precision r / (a_k · r + b_k) of each step with a midpoint rule of 10^4 points. The tolerance is 1e-6.
  - `test_pr_curve_goes_from_recall_1_to_0(samples)`: a hypothesis test of the order. The recall does not increase, the first recall is 1, and the last recall is 0.
  - `test_large_uses_the_recall_step_1_over_p`: the case `large` has at least 100 positives, and its intermediate recalls differ by 1/P.
  - `test_curves_reject_a_nan_score`, `test_curves_reject_a_label_of_2`, `test_curves_reject_labels_of_another_length`, `test_pr_curve_rejects_labels_without_a_positive`, and `test_roc_curve_rejects_labels_without_a_negative`: `InputError` with a `match` on the template.
- [ ] **Step 2: Run the tests.** Run `uv run pytest tests/stats/test_curves.py -q`. Expected: FAIL with `ModuleNotFoundError: No module named 'seedmatch.stats.curves'`.
- [ ] **Step 3: Write `curves.py`.** Export `PrCurve`, `RocCurve`, `pr_curve`, and `roc_curve` from `seedmatch.stats` and from `seedmatch`.
- [ ] **Step 4: Run the tests.** Expected: PASS, and ruff and mypy pass.
- [ ] **Step 5: Write the documents.**
  - `README.md`: remove "the PR and ROC curves" from the second list of "Scope", and add them to the first list.
  - `docs/r-map.md`, row `PRROC::pr.curve` (with `weights.class0` and `curve = TRUE`) → `seedmatch.stats.pr_curve`:
    - `scores.class0` becomes `scores`, and `weights.class0` becomes `labels`, with the values 0 and 1.
    - Soft labels, `scores.class1`, `weights.class1`, `sorted`, `minStepSize`, `max.compute`, `min.compute`, `rand.compute`, and `dg.compute` have no counterpart.
    - `auc.integral` and `auc.davis.goadrich` become `auc_integral` and `auc_davis_goadrich`.
    - The matrix `curve` becomes a `DataFrame` with the columns `recall`, `precision`, and `threshold`.
  - `docs/r-map.md`, row `PRROC::roc.curve` (with the same arguments) → `seedmatch.stats.roc_curve`. `auc` keeps its name. The matrix `curve` becomes a `DataFrame` with the columns `fpr`, `tpr`, and `threshold`. The threshold column keeps the rule of PRROC.
  - `docs/design.md`: add a section "PR and ROC curves" to "Statistics rules", with the rules of this task. Add the bound of the curves to "Performance and bounds". Checkpoint.

### Task 7: The benchmark

**Files:**
- Make: `src/seedmatch/stats/benchmark.py`, `tests/stats/test_benchmark.py`
- Modify: `src/seedmatch/stats/__init__.py`, `src/seedmatch/__init__.py`, `tests/test_readme.py`, `README.md`, `docs/r-map.md`, `docs/divergences.md` (`D-006`, `D-048`), `docs/design.md`

**Interfaces:**
- Consumes: `pr_curve`, `roc_curve`, `PrCurve`, and `RocCurve` (Task 6). `as_cutoff`, `as_text`, `check_columns`, `check_gene_ids_are_text`, and `check_numeric_columns` of `seedmatch._checks`. `benchmark.tsv` and the benchmark curves of Task 1.
- Produces:

```python
@dataclass(frozen=True)
class Benchmark:
    pr: PrCurve
    roc: RocCurve

def benchmark_predictions(
    de: pd.DataFrame,
    *,
    column: str = "mer7m8",
    padj_cutoff: float = 0.15,
    l2fc_cutoff: float = 0.0,
) -> Benchmark: ...
```

Inputs of `benchmark_predictions` (`benchmark_preds` with the fix of `D-006`):
- `de` has the columns `gene_id` (text), `padj`, `log2FoldChange`, and `column` (numeric). `column` goes through `as_text`.
- The cutoffs go through `as_cutoff`, thus a bool, text, and NaN raise `InputError` (`D-048`). For NaN, the message is `padj_cutoff nan is not a number; expected an int or a float, not NaN`.
- The scores are the values of `column` as `float64`. A NaN score raises `InputError` (`D-048`): `column 'mer7m8' of de holds NaN for gene 'g4'; expected a score for each row`.

Rules of the result:
- The gold standard is the set G of the gene IDs of the rows with `padj <= padj_cutoff` and `log2FoldChange < l2fc_cutoff`. A NaN fails each comparison, as `dplyr::filter` drops the row.
- A row is a positive if its `gene_id` is in G, as `%in%` does. Thus each row of a gene in G is a positive.
- No positive raises `no gene has padj <= 1e-10 and log2FoldChange < 0.0; expected at least one positive gene` (`D-048`).
- No negative raises `each row is a positive with padj <= 1.0 and log2FoldChange < 10.0; expected at least one negative row` (`D-048`).
- The result is `Benchmark(pr=pr_curve(scores, labels), roc=roc_curve(scores, labels))`. The ROC curve takes the labels as the PR curve does (`D-006`).

- [ ] **Step 1: Write the failing tests** in `tests/stats/test_benchmark.py`:
  - `test_benchmark_matches_r(case)`: each row of `benchmark.tsv` on `tests/helpers.py::sirna_de()`. Assert `n_positive`, `n_negative`, the areas, and each row of the curves of the same case.
  - `test_the_padj_cutoff_applies_D006`: `padj_cutoff=1e-10` gives 24 positives. The golden note gives 63 for the unpatched R filter.
  - `test_the_roc_curve_uses_the_labels_D006`: the port ROC area equals the `roc_auc` of the golden. The unpatched R area of the note differs.
  - `test_benchmark_rejects_a_nan_score_D048`, `test_benchmark_rejects_a_nan_cutoff_D048`, `test_benchmark_rejects_a_gold_standard_without_a_positive_D048`, and `test_benchmark_rejects_a_gold_standard_without_a_negative_D048`: `InputError` with a `match` on the template. The docstring of each test gives the R value of the golden note.
  - `test_each_row_of_a_gene_in_the_gold_standard_is_a_positive`: a gene with two rows, of which one passes the cutoffs.
  - `test_a_row_without_padj_is_a_negative`.
  - `test_benchmark_rejects_a_missing_column` and `test_benchmark_rejects_a_text_column`.
- [ ] **Step 2: Run the tests.** Run `uv run pytest tests/stats/test_benchmark.py -q`. Expected: FAIL with `ModuleNotFoundError: No module named 'seedmatch.stats.benchmark'`.
- [ ] **Step 3: Write `benchmark.py`.** Export `Benchmark` and `benchmark_predictions` from `seedmatch.stats` and from `seedmatch`.
- [ ] **Step 4: Run the tests.** Expected: PASS, and ruff and mypy pass.
- [ ] **Step 5: Write the README part and the documents.**
  - `README.md`: in "Other analyses", a Python block prints `benchmark_predictions(de).pr.auc_integral` and `.roc.auc`. Add the benchmark to "Scope", and a bullet for `D-006` to "Differences".
  - `tests/test_readme.py`: add `test_readme_benchmark_example_runs`.
  - `docs/r-map.md`, row `benchmark_preds` → `seedmatch.benchmark_predictions`:
    - `res` becomes `de`, and `preds_col` becomes `column`.
    - `padj` and `l2fc` become `padj_cutoff` and `l2fc_cutoff`, keyword-only, with the R defaults.
    - The function returns a `Benchmark` with the fields `pr` and `roc`, not a list with `PRC` and `ROC`. Name `D-006` and `D-048`.
  - `D-006`: "The benchmark applies the padj cutoff, and the ROC curve uses the labels." Write the entry from the row of `docs/design.md`.
    - R location: `benchmark_preds` (`R/benchmark.R:22-38`).
    - R behavior: in `dplyr::filter(padj <= padj, ...)`, both names are the column `padj` of the data mask. Thus each row with a `padj` passes.
    - On the test data with `padj = 1e-10`, R keeps 63 genes, and the cutoff keeps 24.
    - `roc.curve(scores.class0 = preds, scores.class1 = res$GS.label)` takes the labels as the scores of the negative class. Thus R compares the seed counts with the numbers 0 and 1. Evidence: the notes of `benchmark.tsv`.
    - New behavior: the cutoffs apply, and `roc_curve` takes the labels as `pr_curve` does.
    - Reason: a cutoff without an effect gives a wrong gold standard. The R ROC area does not measure the predictions.
  - `D-048`: "The benchmark rejects a NaN score, a cutoff that is not a number, and a gold standard with one class."
    - R location: `benchmark_preds` (`R/benchmark.R:22-38`) and PRROC 1.4.
    - R behavior: give the R values of the notes of `benchmark.tsv` for a NaN score and for no positive. A NaN `l2fc` keeps no gene, because each comparison with NaN fails.
    - New behavior: `InputError` for a NaN score, and for a bool, text, or NaN cutoff. Also `InputError` for a gold standard without a positive or without a negative row.
    - Reason: an area from a broken input looks like a result, and no message tells the caller.
  - `docs/design.md`: the table row of `D-048`. In "Packages", the row of `seedmatch.stats` also names the benchmark. Checkpoint.

### Task 8: The miRDB reader

**Files:**
- Make: `src/seedmatch/mirdb.py`, `tests/test_mirdb.py`
- Modify: `src/seedmatch/__init__.py`, `tests/test_readme.py`, `README.md`, `docs/r-map.md`, `docs/divergences.md` (`D-049`), `docs/design.md`

**Interfaces:**
- Consumes: `local_path` and `open_text` of `seedmatch._io`, `tests/data/mirdb_sample.txt`, and `mirdb.tsv` of Task 1.
- Produces:

```python
def read_mirdb(path: str | os.PathLike[str], *, species: str | None = None) -> pd.DataFrame: ...
```

Inputs of `read_mirdb` (`R/get_data.R:43-52`):
- The path is local. A URL raises `InputError`: `'https://x': the path is a URL; expected a local path` (`local_path`).
- The file is UTF-8 text, plain or compressed, and the first bytes show the compression (`open_text`). This is the file that `seedmatch-fetch mirdb` downloads, `miRDB_v6.0_prediction_result.txt.gz`.
- `species` is a `str` or `None`: `species 5 has type int; expected a str or None`.

Rules of the lines:
- Each line has 3 fields with a tab between them: the miRDB ID, the RefSeq ID of the target, and the score. The file has no header.
- A blank line is skipped, as `read.table` does. An empty file raises `'x.txt': the file has no line; expected miRDB lines`.
- A line with another number of fields raises `'x.txt': line 2 has 2 fields; expected 3: the miRDB ID, the RefSeq ID, and the score`. R stops with `line 2 did not have 3 elements`.
- The line number in a message is the number of the line in the file, with the blank lines. The reader counts the lines, because a parser that skips blank lines numbers the rows in a different way.
- The score is a finite decimal number. Other text raises `'x.txt': line 2: score 'abc' is not a finite number; expected a decimal number` (`D-049`). This includes an empty score, `NA`, `NaN`, and `Inf`.
- The ID fields stay as they are. No quote character and no comment character has a special meaning, and the text `NA` is text (`D-049`).
- `species` is the text of the ID before the first `-`, and `mirna_name` is the text after it. An ID without `-` gives the full ID as `species` and an empty `mirna_name`, as `stringr::str_split_fixed(id, "-", 2)` does.

Rules of the result:
- The columns are `mirdb_id`, `target_refseq_id`, `target_score` (`float64`), `species`, and `mirna_name`, in this order. The text columns have the dtype `object`. The rows keep the order of the file, with a `RangeIndex`.
- With `species`, the result keeps only the rows of this species code.
- A code without a row raises `'x.txt': no row has the species 'rat'; expected one of cfa, gga, hsa, lin4, mmu, rno`. The codes of the message are the sorted codes of the file.
- Read the file in chunks of 2^18 lines, without a quote character (`csv.QUOTE_NONE` for pandas), and with the score as text. Drop the rows of other species in each chunk.
- Make `species` and `mirna_name` from the distinct IDs, and then map them to the rows. Thus equal values share one object.
- The result holds about 40 bytes for each row, and the distinct strings. The docstring states this.

- [ ] **Step 1: Write the failing tests** in `tests/test_mirdb.py`:
  - `test_read_mirdb_matches_r`: `read_mirdb(DATA / "mirdb_sample.txt")` against `mirdb.tsv`, with the R names mapped to the port names.
  - `test_read_mirdb_reads_gzip`: the sample in gzip in `tmp_path` gives the same table.
  - `test_read_mirdb_keeps_only_the_species`: `species="rno"` gives the rows of `rno`, in file order.
  - `test_read_mirdb_rejects_an_unknown_species`: the message lists the codes of the sample.
  - `test_read_mirdb_rejects_a_score_that_is_not_a_number_D049`: parametrized over `abc`, an empty score, `NA`, `NaN`, and `Inf`. The message names the line.
  - `test_read_mirdb_keeps_a_hash_a_quote_and_the_text_na_D049`: an ID with `#`, an ID with `'`, and the ID `NA` stay text.
  - `test_read_mirdb_rejects_a_line_with_2_fields` and `test_read_mirdb_rejects_a_line_with_4_fields`: the message names the line, also after a blank line.
  - `test_read_mirdb_rejects_a_url` and `test_read_mirdb_rejects_an_empty_file`.
  - `test_equal_text_values_share_one_object`: the `species` values of two rows of one species are the same object.
- [ ] **Step 2: Run the tests.** Run `uv run pytest tests/test_mirdb.py -q`. Expected: FAIL with `ModuleNotFoundError: No module named 'seedmatch.mirdb'`.
- [ ] **Step 3: Write `mirdb.py`.** Export `read_mirdb` from `seedmatch`.
- [ ] **Step 4: Run the tests.** Expected: PASS, and ruff and mypy pass.
- [ ] **Step 5: Write the README part and the documents.**
  - `README.md`: in "Other analyses", a Python block gives `read_mirdb("miRDB_v6.0_prediction_result.txt.gz", species="rno")`.
  - `tests/test_readme.py`: the fixture writes the sample in gzip under this name. Add `test_readme_mirdb_example_gives_the_rat_rows`.
  - `README.md`, "Resources": write that `read_mirdb` reads the miRDB file. Also write that no function of this version reads the TRRUST and STRING files.
  - `README.md`, "Scope": add the reader to the first list, and remove "the other small analyses" from the second list.
  - `docs/r-map.md`: in the row `get_example_data`, replace "This version has no reader of the miRDB table." with a reference to `read_mirdb`.
  - `docs/r-map.md`, a new row: `get_example_data` with `example.type = "mirna"` (the table) and `load_example_data("mirna")` → `seedmatch.read_mirdb`:
    - R downloads the file, reads it with `read.table`, and saves `mirdb.Rdata`. The port reads the file of `seedmatch-fetch mirdb`.
    - `miRDB.ID`, `target.REFSEQ.ID`, `target.score`, `species`, and `miRNA.name` become `mirdb_id`, `target_refseq_id`, `target_score`, `species`, and `mirna_name`.
    - The keyword `species` is new. Name `D-049`.
  - `D-049`: "read_mirdb rejects a score that is not a number, and keeps the ID text as it is."
    - R location: `get_example_data` with `example.type = "mirna"` (`R/get_data.R:43-46`, `read.table` with the defaults).
    - R behavior: a score `abc` makes the score column text. Then the vignette comparison `mirdb.rat[,3] >= 90` compares text.
    - An empty score gives `NA`, and the ID text `NA` becomes `NA`, the R value for no data.
    - A `#` in a field starts a comment, and R stops with `line 2 did not have 3 elements`. Evidence: the notes of `mirdb.tsv`.
    - New behavior: `InputError` for a score that is not a finite decimal number. The ID fields keep `#`, `'`, and `NA` as text.
    - Reason: a text score gives wrong comparisons and no message. The miRDB format has no quotes and no comments.
  - `docs/design.md`: the table row of `D-049`. In "Packages", add the row of `seedmatch.mirdb`: the reader of the miRDB table. Add the bound of `read_mirdb` to "Performance and bounds". Checkpoint.

### Task 9: Verification of the plan

- [ ] **Step 1: Run the checks one time.** Run `uv run ruff check`, `uv run ruff format --check`, `uv run mypy`, and `NUMBA_DISABLE_JIT=1 uv run pytest --cov`. Expected: each passes, and the coverage is not below 95%. Give a reason for each line of the new modules without coverage.
- [ ] **Step 2: Compare the documents with the code.**
  - Each public name of the new modules has a row in `docs/r-map.md`, or is a type of a function that has a row.
  - Each `D-NNN` of this plan has an entry, a table row, and a test.
  - Each README example runs in `tests/test_readme.py`.
  - `tests/test_package.py::test_import_loads_no_fetch_and_no_plot_modules` passes.
  - Checkpoint.
