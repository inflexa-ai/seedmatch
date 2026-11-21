# Documentation and performance: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Do the last delivery step of `docs/design.md` ("Delivery order", step 10): the README, the map from R names, the divergences, and the performance record. The plan adds `docs/performance.md` with its measurement scripts. It also corrects each defect that a comparison of the documents with the code found.

**Architecture:** `tools/perf/` holds the scripts that measure `seed_match_report` and `SeedMatchReport` on one FASTA file of the canonical 3' UTRs of `rnor7`. New tests keep the documents true. `tests/test_r_map.py` compares `docs/r-map.md` and `README.md` with the R exports and with the public names of seedmatch. A new test in `tests/test_package.py` compares `docs/divergences.md` with the table of `docs/design.md`. This plan changes no file in `src/`.

**Tech Stack:** Python 3.12.13, pytest, numba 0.67.0. R 4.5.3 with SeedMatchR 2.0.0, Biostrings 2.78.0, dplyr 1.2.1, and gt 1.4.0 for the R measurement.

**Spec:** `docs/design.md`, sections "Goals", "Divergences (defect fixes)", "Tests > Documentation examples", "Performance and bounds", and "Delivery order".

**Later plans (not in this plan):** none.

## The measurement of the planner

The planner measured the two calls on this machine. The scripts and their output are in `tmp/claude/plan-docs/`: `tools/perf/` (the scripts of Task 2), `r_run_1.out`, `r_run_2.out`, `py_fasta_run1.out`, `py_fasta_run2.out`, `py_run_2.out`, `py_fasta_workers1.out`, and `py_run_2_workers1.out`.

The calls:

- Python: `seed_match_report(de, sequences, guide)`, with the default `ReportOptions()` and the default `workers=None`.
- R: `SeedMatchR::SeedMatchReport(res, seqs, guide.seq)`, with the default arguments of `R/SeedMatchReport.R:46-55`. The row `SeedMatchReport` of `docs/r-map.md` maps each argument, and the defaults agree.
- The sequences are the canonical 3' UTRs of `rnor7`: 18,614 sequences and 17,913,962 bases, with one sequence for each gene.
- Both calls read one FASTA file, with the gene ID as the record ID. R names the sequences of `load_annotations` with their gene IDs (`R/load_annotations.R:105-107`). `SeedMatchReport` joins the sequences with the DE table by these names (`R/SeedMatchReport.R:67-70`).
- The DE table is `tests/data/DESEQ2_results.tsv`, the test data that `tools/golden/R/report.R` reads. The example data of `get_example_data("sirna")` is not on this machine. The table changes only the join after the search.
- The two reports are equal in each cell. The row `In silico predictions` has the total 11,477 in each language.

The machine: Apple M2 Max, 12 CPUs (8 performance cores, 4 efficiency cores), 32 GiB, macOS 26.5.2. Python 3.12.13, seedmatch 0.1.0 at commit `6d17d56`, numba 0.67.0, llvmlite 0.49.0, numpy 2.5.2, pandas 2.3.3. R 4.5.3, SeedMatchR 2.0.0 (`2fc0c4f`), Biostrings 2.78.0, dplyr 1.2.1, gt 1.4.0.

Other processes: the machine was not quiet. Two `pytest` processes of an earlier session (PIDs 7858 and 16083) used 100% of one CPU each during each run. They ran `pytest tests/test_kmers.py tests/plot/test_plot_kmers.py tests/plot/test_plot_package.py -q -x` for 4.5 hours. The load average was 2.8 to 4.7.

Each measure has two runs. The planner measured `workers=1` only with a filled numba cache.

| Measure | Python, `workers=None` | Python, `workers=1` | R |
|---|---|---|---|
| First call, empty numba cache | 2.86 s, 2.88 s | not measured | 9.59 s, 9.63 s |
| First call, filled numba cache | 1.93 s, 1.94 s | 8.69 s, 8.73 s | no cache |
| Later calls | 1.74 to 2.06 s | 8.62 to 8.64 s | 9.32 to 9.35 s |
| Resident memory before the first call | 144 MiB | 145 MiB | 296 MiB, 313 MiB |
| Peak resident memory of the process | 244 MiB (empty cache), 205 MiB and 206 MiB (filled cache) | 203 MiB, 205 MiB | 432 MiB, 443 MiB |
| Peak above the memory before the first call | 100 MiB (empty cache), 61 MiB and 62 MiB (filled cache) | 58 MiB, 60 MiB | 136 MiB, 130 MiB |

- After the JIT compilation, R takes 4.5 to 5.4 times the Python time. With `workers=1`, Python takes 0.92 to 0.93 times the R time. Thus the target of the design holds, and the planner made no profile.
- The JIT compilation adds about 0.9 seconds to the first call. `import seedmatch` takes 0.42 to 0.63 seconds, and it compiles no kernel.
- The R call also makes the `gt` table (`R/SeedMatchReport.R:135`). In Python, `report_html` makes the table in less than 1 ms.

## The comparison of the documents with the code

The scripts and their output are in `tmp/claude/plan-docs/`. Task 3 and Task 4 fix each defect.

`check_r_map.py ../SeedMatchR` (output `check_r_map.out`) gives these lines for the R functions and the Python names:

```text
exported R functions without a row: ['.get_counts_column', '.get_match_df']
internal R functions that the map names without a row: none
internal R functions that the map does not name: ['.get_algo', '.get_alignment', '.get_match_seq', '.normalize_fixed_arg', '.run_scanning', '.valid_algos']
Python names that do not resolve:
  0 problems
```

`check_public_names.py` (output `check_public_names.out`) lists the public names that neither `docs/r-map.md` nor `README.md` names. A name that `seedmatch` does not re-export counts only with its module path, for example `seedmatch.resources.Group`. Without this rule, the R column `Group` of the report hides `seedmatch.resources.Group`.

```text
seedmatch: ['DEFAULT_SEEDS', 'NONSEED_SEEDS', 'SeedMatchError', 'write_fasta']
seedmatch.annotation: ['FastaGenome', 'Genome']
seedmatch.fetch: ['DEFAULT_FETCH_OPTS', 'FetchError', 'FetchOpts', 'HttpResponse', 'fetch_resource', 'open_url']
seedmatch.report: ['DEFAULT_SEEDS', 'NONSEED_SEEDS']
seedmatch.resources: ['REGISTRY', 'STRING_VERSION', 'Entry', 'Group', 'Resource', 'lookup']
seedmatch.search: ['MAX_PATTERN_LENGTH', 'match_ranges']
seedmatch.sequence: ['write_fasta']
seedmatch.stats: ['TwoSampleResult']
```

The other packages give `none`.

`check_divergences.py` (output `check_divergences.out`): the row order of the design table gives the IDs `D-001` to `D-070`. Row 14 is the row `Documentation`. Each entry has each field, but `D-014` has no entry.

The defects:

1. `docs/r-map.md` has no row for the exported `.get_counts_column` (`R/SeedMatchR.R:183`) and `.get_match_df` (`R/SeedMatchR.R:252`). Fix: Task 3, Step 3.
2. `docs/r-map.md:3` maps "each SeedMatchR function", but no row names `.run_scanning` (`R/full_search.R:191`), `.get_match_seq` (`:225`), or `.get_alignment` (`:231`). No row names `.normalize_fixed_arg`, `.get_algo`, or `.valid_algos` (`R/vmatchPattern2.R:66`, `:104`, `:119`). Fix: Task 3, Step 3.
3. `docs/r-map.md:16` maps `load_example_data` only for `"mirna"`. The default `"sirna"` (`R/get_data.R:96-99`) has no row. Fix: Task 3, Step 3.
4. The public names of the second list above have no row and no README text. `docs/r-map.md:32` names `FastaGenome` and `Genome` without their module. `docs/r-map.md:40` names `TwoSampleResult` without its module. Fix: Task 3, Steps 4 and 5.
5. `README.md:78` gives p = 1.04e-4, 0.9996, and 3.9e-4 "for this table". These are the values of the raw table with its column `mer7m8`, but the `de` of the example above gives 0.0011, 0.9986, and 0.0015. The evidence is the probe `readme_claims.out` and the golden rows `readme` of `tests/golden/stats/ecdf.tsv`. Fix: Task 3, Step 5.
6. `README.md:374` "This version has the core path:" is stale. The list below it now holds each part of the port. Fix: Task 3, Step 5.
7. The comment of `tests/test_package.py:56` states that the README names `SeedMatchError`, `DEFAULT_SEEDS`, and `NONSEED_SEEDS`. The README names none of them. Fix: Task 3, Step 6.
8. `docs/divergences.md` has no entry `D-014` (design row 14, `docs/design.md:464`), and no test name contains `D014`. The entries go from `D-013` (`:179`) to `D-015` (`:189`). Fix: Task 4.
9. The R behavior of `D-052` to `D-055` (`docs/divergences.md:655`, `:666`, `:677`, `:689`) names golden cases, but not their golden file `tests/golden/noncanonical/compensatory.tsv`. Fix: Task 4, Step 4.
10. `docs/design.md:592` states that a test runs each code example in `README.md` and in `docs/`. `tests/test_readme.py` runs only the Python blocks of `README.md`. The one Python block of `docs/` is the sketch of `docs/design.md:125-146`, with placeholder paths. Fix: Task 4, Step 5.

`README.md:9` ("Each change of behavior from R has an entry") becomes true with fix 8. `README.md:96` ("the Python name of each R function") becomes true with fixes 1 to 3.

The search found no stale "later plan" or "not in this plan" text in `docs/design.md` or `README.md`.

The design sketch (`docs/design.md:125-146`) agrees with the code. The probe `design_example.py` runs it with the fixture paths. It uses `Test_seqs.fa` as `utr3` after the `Reference` call. Each name and each keyword of the sketch exists (`design_example.out`).

`tests/test_readme.py` passes on the current tree.

The probe `simulate_after.py` applies the texts of Task 2, Task 3, and Task 4 to copies of the documents in `tmp/claude/plan-docs/after/`. The planned tests pass on the copies (`simulate_after.out`). `check_r_map.py` gives `none` for each line on the copies (`check_r_map_after.out`).

## Global Constraints

- No new runtime dependency. No change in `src/`.
- The tests never run R, and they never use the network.
- The R source of truth is `../SeedMatchR` at commit `2fc0c4f`. The R scripts run with `tools/golden/run.sh`.
- Target (`docs/design.md`, "Performance and bounds"): "after the JIT compilation, Python is not slower than R by more than a factor of 2." If a measurement does not obey the target, stop and report it to the controller.
- Each change of behavior from R has an entry `## D-NNN: <title>` in `docs/divergences.md` and a row in the table of `docs/design.md`. The name of a test contains the ID.
- This plan writes the entry `D-014` and makes no new ID.
- A rename of an argument or a value, or a change of a type, is not a divergence. `docs/r-map.md` records it.
- Documents use ASD-STE100: sentences of 25 words or fewer, no semicolons, no `-ing` verb forms. Code comments obey the admission test of the project rules. The banned words of the project rules do not occur.
- Do not state the size of a collection that another source owns. A measured input, a limit, and a result are permitted.
- Git: an implementer does no `git add` and no commit. Each task ends with a checkpoint. The controller commits only with the consent of the user.

## Review Focus

1. A cold run that reuses an old numba cache. The cache of `src/seedmatch/search/__pycache__` or of the user cache hides the JIT compilation. The cold command sets `NUMBA_CACHE_DIR` to a new empty directory, thus the first call takes about 0.9 seconds more than with the filled cache. Owner: Task 2.
2. Two calls on different input. The FASTA IDs must be gene IDs, because R joins by the names and Python by `gene_ids`. Both scripts print `In silico predictions, Total: 11477`. Owner: Task 2.
3. A cell of `docs/r-map.md` with ` | ` in it. The split then gives a row of more than 3 cells, and a test can miss a name. `r_map_rows` asserts 3 cells for each row. Owner: Task 3.
4. A public name that occurs only as an R column name, for example `Group`. A name that `seedmatch` does not re-export counts only with its module path. Owner: Task 3.
5. A new row of the design table without an entry, or an entry without a field. The new test of `tests/test_package.py` fails for both. An entry with a repeated field passes, for example the `R behavior` fields of `D-064`. Owner: Task 4.

---

### Task 1: The plan file and the R exports

**Files:**
- Make: `docs/plans/2026-11-21-docs-and-performance.md` (this plan, written by the planner)
- Make: `tests/data/SeedMatchR.NAMESPACE`
- Modify: `tests/data/README.md`

**Interfaces:**
- Produces: `tests/data/SeedMatchR.NAMESPACE`, a byte copy of `../SeedMatchR/NAMESPACE`. Task 3 reads its `export(...)` lines.

This plan makes no R golden. The R probe of `D-014` is `tmp/claude/plan-docs/d014_probe.R`, and its output is `d014_probe.out`.

- [ ] **Step 1: Copy the NAMESPACE file.** Run `cp ../SeedMatchR/NAMESPACE tests/data/SeedMatchR.NAMESPACE`. Expected: `git -C ../SeedMatchR rev-parse --short HEAD` gives `2fc0c4f`. Then `cmp ../SeedMatchR/NAMESPACE tests/data/SeedMatchR.NAMESPACE` gives no output.
- [ ] **Step 2: Document the copy.** In `tests/data/README.md`, after the first paragraph, add this paragraph:

```markdown
`SeedMatchR.NAMESPACE` is a copy of the file `NAMESPACE` of SeedMatchR at commit `2fc0c4f`. The bytes are the same as in the source. `tests/test_r_map.py` reads its `export(...)` lines.
```

- [ ] **Step 3: Checkpoint.**

### Task 2: The performance record

**Files:**
- Make: `tools/perf/utr3_fasta.py`, `tools/perf/report.py`, `tools/perf/report.R`
- Make: `docs/performance.md`
- Modify: `README.md` (section "Development")

**Interfaces:**
- Consumes: `Reference.from_directory`, `Reference.sequences`, `write_fasta`, `read_fasta`, `seed_match_report`, `report_html`, and `load_test_data` of `tools/golden/R/common.R`.
- Produces the commands below. No other task uses them.
  - `uv run python tools/perf/utr3_fasta.py REFERENCE OUT`. `REFERENCE` is the directory of the EnsDb and the 2bit file. The script writes the canonical 3' UTRs to `OUT`, with the gene ID as the record ID. It prints `<OUT>: <n> sequences, <bases> bases`. If a gene has more than one sequence, the script stops with a message.
  - `uv run python tools/perf/report.py FASTA [--calls N] [--workers W]`. The default of `N` is 3, and the default of `W` is `None`. The script prints the import time and the peak before the first call. Then it prints the time and the peak after each call, the time of `report_html`, the total line, and the versions.
  - `PERF_FASTA=FASTA [PERF_CALLS=N] tools/golden/run.sh tools/perf/report.R`. The script prints the resident memory before the first call. Then it prints the time and the resident memory after each call, the total line, and the versions.
  - The total line is `In silico predictions, Total: <n>`.

- [ ] **Step 1: Copy the scripts.** Copy the three files of `tmp/claude/plan-docs/tools/perf/` to `tools/perf/` without a change. The planner ran each file. `ruff check` and `ruff format --check` pass on the two Python files. `report.py` reads the DE table from the relative path `tests/data/DESEQ2_results.tsv`, thus each command runs from the repository root.
- [ ] **Step 2: Lint.** Run `uv run ruff check tools/perf` and `uv run ruff format --check tools/perf`. Expected: no finding.
- [ ] **Step 3: Measure.** Run these commands from the repository root, one at a time. Do not run other work of this plan at the same time. `SEEDMATCH_TEST_REFS` is the directory of the `refs` tests, for example `~/Development/refs/seedmatch`.

```bash
uv run python tools/perf/utr3_fasta.py "$SEEDMATCH_TEST_REFS/rnor7" tmp/perf/rnor7_utr3.fa
shasum -a 256 tmp/perf/rnor7_utr3.fa
rm -rf tmp/perf/numba && mkdir -p tmp/perf/numba
NUMBA_CACHE_DIR="$PWD/tmp/perf/numba" uv run python tools/perf/report.py tmp/perf/rnor7_utr3.fa
NUMBA_CACHE_DIR="$PWD/tmp/perf/numba" uv run python tools/perf/report.py tmp/perf/rnor7_utr3.fa
NUMBA_CACHE_DIR="$PWD/tmp/perf/numba" uv run python tools/perf/report.py tmp/perf/rnor7_utr3.fa --workers 1
PERF_FASTA=tmp/perf/rnor7_utr3.fa /usr/bin/time -l tools/golden/run.sh tools/perf/report.R
```

  Expected:
  - The FASTA file has 18614 sequences and 17913962 bases. Its SHA-256 is `ac679b985d96394f943fcefcec3b94e91c92c12b5683b7b0fae063389f738edb`.
  - Each run prints `In silico predictions, Total: 11477`.
  - The times and the memory are near the table of "The measurement of the planner". The later Python calls with `workers=None` take less than 2 times the later R calls.
  - If a later Python call takes 2 times the R call or more, stop. Profile the slow call with `python -m cProfile -s cumtime`, and report the profile to the controller.
  - If a time differs from the table by more than 25%, put the new values in `docs/performance.md`. Also give the load average of `uptime` and the busy processes of `ps -Ao pcpu,comm -r`.
- [ ] **Step 4: Write `docs/performance.md`.** Write the three parts below in this order. Put the second part in a fenced `bash` block. If Step 3 changed a value, use the new value, and keep the two runs of each measure.

  The first part:

```markdown
# Performance of seed_match_report

This document records the wall time and the peak memory of `seed_match_report` on the canonical 3' UTRs of `rnor7`, for Python and for R. `docs/design.md` ("Performance and bounds") sets the target: after the JIT compilation, Python is not slower than R by more than a factor of 2.

## The calls

- Python: `seed_match_report(de, sequences, guide)`, with the default `ReportOptions()` and the default `workers=None`. `workers=None` gives a thread pool with one thread for each CPU.
- R: `SeedMatchR::SeedMatchReport(res, seqs, guide.seq)`, with the default arguments (`R/SeedMatchReport.R:46-55`). R uses one thread.
- The defaults of the two calls agree. The row `SeedMatchReport` of [`docs/r-map.md`](r-map.md) maps each argument.
- The guide is `UUAUAGAGCAAGAACACUGUUUU`, the guide of the R goldens.
- The R call also makes the `gt` table (`R/SeedMatchReport.R:135`). In Python, `report_html` makes the table in a separate call, and this call takes less than 1 ms.

## The input

- The sequences are the canonical 3' UTRs of `rnor7`. They come from the EnsDb 109 and the 2bit file of `seedmatch-fetch reference rnor7`, with the default `TranscriptFilter` and `min_width=8`.
- The set has 18,614 sequences and 17,913,962 bases, and each gene has one sequence.
- `tools/perf/utr3_fasta.py` writes the sequences to one FASTA file, with the gene ID as the record ID. R reads the file with `Biostrings::readDNAStringSet`, and Python reads it with `read_fasta`. Thus the two calls search the same sequences.
- R names each sequence of `load_annotations` with its gene ID (`R/load_annotations.R:105-107`), as the file does. The `refs` tests of `tests/annotation/test_refs.py` show that `Reference.sequences` gives the sequences of `load_annotations`.
- The DE table is `tests/data/DESEQ2_results.tsv`, the test data of SeedMatchR that `tools/golden/R/report.R` reads. The table changes only the join after the search.
- The two calls give the same report. For example, the row `In silico predictions` has the total 11,477 in each language.

## The machine

- Apple M2 Max with 12 CPUs (8 performance cores and 4 efficiency cores), 32 GiB of memory, and macOS 26.5.2.
- Python 3.12.13, seedmatch 0.1.0 (commit `6d17d56`), numba 0.67.0, llvmlite 0.49.0, numpy 2.5.2, and pandas 2.3.3.
- R 4.5.3, SeedMatchR 2.0.0 (commit `2fc0c4f`), Biostrings 2.78.0, dplyr 1.2.1, and gt 1.4.0.
- Two other processes used one CPU each during each run. The load average was 2.8 to 4.7.

## The commands

Run the commands from the repository root. `SEEDMATCH_TEST_REFS` is the directory of the `refs` tests. On Linux, use `/usr/bin/time -v` for the peak memory of the R process.
```

  The second part, the commands:

```bash
# The canonical 3' UTRs of rnor7, with the gene ID as the record ID.
uv run python tools/perf/utr3_fasta.py "$SEEDMATCH_TEST_REFS/rnor7" tmp/perf/rnor7_utr3.fa

# Python: the first run compiles the kernels into an empty cache. The second run reads the cache.
rm -rf tmp/perf/numba && mkdir -p tmp/perf/numba
NUMBA_CACHE_DIR="$PWD/tmp/perf/numba" uv run python tools/perf/report.py tmp/perf/rnor7_utr3.fa
NUMBA_CACHE_DIR="$PWD/tmp/perf/numba" uv run python tools/perf/report.py tmp/perf/rnor7_utr3.fa
NUMBA_CACHE_DIR="$PWD/tmp/perf/numba" uv run python tools/perf/report.py tmp/perf/rnor7_utr3.fa --workers 1

# R: three calls in one session.
PERF_FASTA=tmp/perf/rnor7_utr3.fa /usr/bin/time -l tools/golden/run.sh tools/perf/report.R
```

  The third part:

```markdown
Each script makes three calls in one process. The Python script gives the peak resident memory of its process after each call. The R script gives the resident memory after each call, and `/usr/bin/time` gives the peak of the R process.

## The results

Each cell gives the values of two runs. `workers=1` has only runs with a filled numba cache.

| Measure | Python, `workers=None` | Python, `workers=1` | R |
|---|---|---|---|
| First call, empty numba cache | 2.86 s, 2.88 s | not measured | 9.59 s, 9.63 s |
| First call, filled numba cache | 1.93 s, 1.94 s | 8.69 s, 8.73 s | no cache |
| Later calls | 1.74 to 2.06 s | 8.62 to 8.64 s | 9.32 to 9.35 s |
| Resident memory before the first call | 144 MiB | 145 MiB | 296 MiB, 313 MiB |
| Peak resident memory of the process | 244 MiB (empty cache), 205 MiB and 206 MiB (filled cache) | 203 MiB, 205 MiB | 432 MiB, 443 MiB |
| Peak above the memory before the first call | 100 MiB (empty cache), 61 MiB and 62 MiB (filled cache) | 58 MiB, 60 MiB | 136 MiB, 130 MiB |

- After the JIT compilation, R takes 4.5 to 5.4 times the Python time. Thus the target holds.
- With `workers=1`, Python takes 0.92 to 0.93 times the R time. Thus the target also holds on one thread.
- The JIT compilation adds about 0.9 seconds to the first call. `import seedmatch` takes 0.42 to 0.63 seconds, and it compiles no kernel.
- The memory above the start is the work memory of the call and its result. The empty cache adds the memory of the compiler.
```

- [ ] **Step 5: Link the record.** In `README.md`, section "Development", after the line that names `tools/golden/README.md` (`README.md:440`), add this paragraph:

```markdown
[`docs/performance.md`](docs/performance.md) gives the time and the memory of `seed_match_report` on `rnor7`, and the commands that measure them.
```

- [ ] **Step 6: Run the README tests.** Run `uv run pytest tests/test_readme.py -q`. Expected: PASS.
- [ ] **Step 7: Checkpoint.**

### Task 3: The map from R names and the public names

**Files:**
- Make: `tests/test_r_map.py`
- Modify: `docs/r-map.md`, `README.md`, `tests/test_package.py:56` (a comment)

**Interfaces:**
- Consumes: `tests/data/SeedMatchR.NAMESPACE` (Task 1).
- Produces: `test_each_python_name_of_the_r_map_exists_D014` in `tests/test_r_map.py`. The entry `D-014` of Task 4 names it.

- [ ] **Step 1: Write the tests.** Copy `tmp/claude/plan-docs/test_r_map_scratch.py` to `tests/test_r_map.py` with these changes:
  - Set `ROOT = Path(__file__).parents[1]` and `NAMESPACE = ROOT / "tests" / "data" / "SeedMatchR.NAMESPACE"`.
  - Remove the comment about the planner copy.

  The file has these helpers and tests:
  - `r_map_rows()` gives the cells of each row of the table, without the header and the separator. It asserts that each row has 3 cells.
  - `r_names(cell)` gives the R name of each code span of a cell, without a `package::` prefix and without the arguments. For example, `load_example_data("mirna")` gives `load_example_data`.
  - `code_spans(path)` gives the code spans of each line of a document. A fenced block breaks the pairing of the backticks over the whole text. Thus a search over the whole README misses names.
  - `test_each_exported_r_function_has_a_row_of_the_r_map`: each name of an `export(...)` line of the NAMESPACE copy is an R name of the first column of a row.
  - `test_each_python_name_of_the_r_map_exists_D014`: each dotted `seedmatch.` code span of the second and the third column resolves with `importlib` and `getattr`. Each `seedmatch-fetch` command names a `Group` member. Its first key that is not uppercase must pass `lookup`. An uppercase word is a placeholder, for example `BUILD` or `DEST`.
  - `test_each_public_name_is_in_the_r_map_or_the_readme(module_name)`: the parameters are `seedmatch` and each package of `pkgutil.iter_modules(seedmatch.__path__)` whose name does not start with `_`. Each name of `__all__` occurs in a code span of `docs/r-map.md` or `README.md`, or in an `import` of a Python block of `README.md`.
  - A name of `seedmatch.__all__` can stand alone. Another name counts only after its module path, for example `seedmatch.resources.lookup`.
- [ ] **Step 2: Run the tests.** Run `uv run pytest tests/test_r_map.py -q`. Expected: 9 failed and 4 passed. The failures:
  - `test_each_exported_r_function_has_a_row_of_the_r_map`: `['.get_counts_column', '.get_match_df']`.
  - `test_each_public_name_is_in_the_r_map_or_the_readme`: the lists of "The comparison of the documents with the code", for `seedmatch`, `seedmatch.annotation`, `seedmatch.fetch`, `seedmatch.report`, `seedmatch.resources`, `seedmatch.search`, `seedmatch.sequence`, and `seedmatch.stats`.
  - `test_each_python_name_of_the_r_map_exists_D014` passes now. It keeps the names that exist.
- [ ] **Step 3: Add the rows of the R functions.** Add each row to `docs/r-map.md` at the given place.
  - After the row `vmatchPattern2`:

```markdown
| `.normalize_fixed_arg`, `.get_algo`, `.valid_algos` | `seedmatch.SearchOptions`, `seedmatch.search.MAX_PATTERN_LENGTH` | These helpers make sure that the arguments of `vmatchPattern2` are valid (`R/vmatchPattern2.R:29-33`). `.normalize_fixed_arg` takes `fixed` as one or two flags or as names, and `SearchOptions` has the two flags `fixed_pattern` and `fixed_subject`. For an empty pattern or a pattern of more than 20000 letters, `.valid_algos` stops (`:131-134`), and `count_matches` and `find_matches` raise `InputError`. For `min.mismatch` above 0 with indels and `max.mismatch` above 0, `.valid_algos` stops (`:135-137`). `SearchOptions` raises `InputError` for `min_mismatch` above 0 with `indels` on. `vmatchPattern2` gets an algorithm from `.get_algo`, but it gives `"indels"` to the C code (`:43`, `D-001`). |
```

  - After the row `get_example_data` with `example.type = "mirna"` (the table):

```markdown
| `load_example_data` with `example.type = "sirna"` | `seedmatch.read_deseq2_table` | `"sirna"` is the default. R loads the DESeq2 tables that `get_example_data("sirna")` saves in `Schlegel_2022.Rdata`, and it returns them as a list (`R/get_data.R:96-99`). `seedmatch-fetch example sirna DEST` downloads the files, and `read_deseq2_table` reads each file. The port saves no `.Rdata` file. |
```

  - After the row `SeedMatchR` with `res.format = "DESEQ2"`:

```markdown
| `.get_counts_column` | `seedmatch.annotate_seed_counts` | `SeedMatchR` with `res.format = "DESEQ2"` calls this helper (`R/SeedMatchR.R:95-104`), thus the row of that mode gives the arguments. A direct call has the defaults `with.indels = F` and `fixed = F` (`R/SeedMatchR.R:190-191`). `options=SearchOptions(fixed_pattern=False, fixed_subject=False)` gives these defaults. The default of `annotate_seed_counts` is the default of `SeedMatchR`. `seed.name` and `col.name` become `column`. The count of a gene is the maximum over its sequences (`D-003`), and the counts replace a column with the same name (`D-015`). |
```

  - After the row `full_search_ranges`:

```markdown
| `.run_scanning` | `seedmatch.full_search_ranges` | `full_search` calls this helper for each search (`R/full_search.R:170-179`). It names a search `<seed>.mm<distance>.indel<flag>`, and it adds the columns `query`, `TargetSeq`, `seqnames`, `EditDistance`, and `AllowIndels` (`:191-223`). The row of `full_search_ranges` gives the Python columns. |
```

  - After the row `.get_best_match`:

```markdown
| `.get_match_seq`, `.get_alignment` | no counterpart | No SeedMatchR function calls them. The calls of `get_match_seq` are comments (`R/SeedMatchR.R:276`, `R/full_search.R:217`). `.get_alignment` reads the column `MatchString`, and no code makes this column (`R/full_search.R:231-238`). `alignment_text` gives the alignment of a match. |
```

  - After the row `SeedMatchR` with `res.format = "data.frame"`:

```markdown
| `.get_match_df` | `seedmatch.find_matches` | `SeedMatchR` with `res.format = "data.frame"` or `"granges"` calls this helper (`R/SeedMatchR.R:110-117`), thus the row of that mode gives the arguments. A direct call has the defaults `with.indels = TRUE` and `fixed = TRUE` (`R/SeedMatchR.R:257-258`), thus `SearchOptions(indels=True)`. `seed.name` and `sirna.name` have no counterpart. R uses them only for the columns `sirna`, `query`, `TargetSeq`, `EditDistance`, and `AllowIndels` (`:272-279`), and the port does not make these columns. The ranges obey `indels` (`D-001`). |
```

- [ ] **Step 4: Name the public names in the map.** Change these rows of `docs/r-map.md`:
  - Row `Biostrings::vmatchPattern`: at the end of the third cell, add "`seedmatch.search.match_ranges` gives the same matches as three `int64` arrays: the index of the sequence in the set, the start, and the end."
  - Row `full_search` with `res_format = "counts"`: after "`seqs.to.search` becomes `seeds` (a tuple of `Seed`).", add "The default is `DEFAULT_SEEDS`, the seeds of R in the order of R. `nonseed_edit_distances` applies to the seeds of `NONSEED_SEEDS`: `Seed.FULL`, `Seed.MER18`, and `Seed.MER15`, as in `R/full_search.R:44`."
  - Row `rtracklayer::TwoBitFile`: replace "`Genome.codes(contig, start, end)` gives" with "`seedmatch.annotation.Genome` is the type of each genome reader. `Genome.codes(contig, start, end)` gives". Replace "`FastaGenome` reads" with "`seedmatch.annotation.FastaGenome` reads".
  - Row `stats::ks.test`: replace "The function returns a `TwoSampleResult` with" with "The function returns a `seedmatch.stats.TwoSampleResult` with".
- [ ] **Step 5: Correct the README.** Change these parts of `README.md`:
  - The third bullet of "A part of an agentic toolchain for bioinformatics" (`README.md:19`): at the end, add "Each error of seedmatch is a `SeedMatchError`, and `InputError` is also a `ValueError`."
  - After the line "If a gene has more than one sequence, ..." (`README.md:60`), add a paragraph: "`write_fasta(sequences, path)` writes a set to a FASTA file, with the ID of each sequence in its header line."
  - Replace the paragraph of `README.md:78` with these two paragraphs:

```markdown
The default alternative `GREATER` has another meaning for `EcdfTest.KS` and for `EcdfTest.WILCOXON`. With KS, a small p-value shows that the target values are lower. With WILCOXON, a small p-value shows that the target values are higher. To find lower target values with WILCOXON, set `alternative=Alternative.LESS`. Import `Alternative` from `seedmatch`.

In the test table of this repository (`tests/data/DESEQ2_results.tsv`), the genes with a seed match have lower values. The table has a column `mer7m8`. With this column and without `filter_de`, KS with `GREATER` gives p = 1.04e-4, but WILCOXON with `GREATER` gives p = 0.9996. WILCOXON with `Alternative.LESS` gives the small p-value 3.9e-4. For the `de` of the example above, the three p-values are 0.0011, 0.9986, and 0.0015, as in R.
```

  - Section "Resources", after the paragraph that starts "`python -m seedmatch.fetch` runs the same command." (`README.md:347`), add these two paragraphs:

```markdown
Python code can use the registry without the command. `seedmatch.resources.REGISTRY` holds each `seedmatch.resources.Entry`, and an entry holds its files as `seedmatch.resources.Resource` objects. `seedmatch.resources.lookup(group, key)` returns the entry of a `seedmatch.resources.Group` member and a key. `seedmatch.resources.STRING_VERSION` is the pinned STRING version.

`seedmatch.fetch.fetch_resource(resource, dest)` downloads one resource into the directory `dest`, and it returns the path of the file. A failed download raises `seedmatch.fetch.FetchError`. The keyword `opts` takes a `seedmatch.fetch.FetchOpts`, and its field `open_url` is the HTTP client. A client takes a URL and returns a `seedmatch.fetch.HttpResponse`. `seedmatch.fetch.DEFAULT_FETCH_OPTS` uses `seedmatch.fetch.open_url`, a client of `urllib.request`.
```

  - Section "Scope" (`README.md:374`): replace "This version has the core path:" with "This version has these parts:".
- [ ] **Step 6: Correct the comment of `tests/test_package.py:56`.** Replace "The README and the docstrings of the options name these" with "The README, docs/r-map.md, and the docstrings of the options name these".
- [ ] **Step 7: Run the tests.** Run `uv run pytest tests/test_r_map.py tests/test_readme.py tests/test_package.py -q`, and expect PASS. Then run `uv run ruff check` and `uv run ruff format --check`, and expect no finding.
- [ ] **Step 8: Run the full check of the map.** Run `uv run python tmp/claude/plan-docs/check_r_map.py ../SeedMatchR`. Expected:
  - `exported R functions without a row: none`
  - `internal R functions that the map names without a row: none`
  - `internal R functions that the map does not name: []`
  - `0 problems` for the Python names
  - `none` for each package of the public names
- [ ] **Step 9: Checkpoint.**

### Task 4: The divergence record

**Files:**
- Modify: `tests/test_package.py`, `docs/divergences.md`, `docs/design.md:592`

**Interfaces:**
- Consumes: `test_each_python_name_of_the_r_map_exists_D014` (Task 3).

- [ ] **Step 1: Write the failing test.** Add to `tests/test_package.py`:

```python
DESIGN = Path(__file__).parents[1] / "docs" / "design.md"
DIVERGENCE_FIELDS = ("ID", "R location", "R behavior", "New behavior", "Reason", "Test")


def test_each_divergence_row_of_the_design_has_an_entry_with_each_field():
    design = DESIGN.read_text(encoding="utf-8")
    section = design.split("\n## Divergences (defect fixes)\n", 1)[1].split("\n## ", 1)[0]
    # The table has no ID column: row n has the ID D-00n. The separator line "|---|" is not a row.
    rows = [line for line in section.splitlines() if line.startswith("| ")][1:]
    text = DIVERGENCES.read_text(encoding="utf-8")
    entries = re.split(r"^## (?=D-\d{3}:)", text, flags=re.MULTILINE)[1:]
    ids = [entry.split(":", 1)[0] for entry in entries]
    assert ids == [f"D-{n:03d}" for n in range(1, len(rows) + 1)]
    for entry_id, entry in zip(ids, entries, strict=True):
        fields = re.findall(r"^- ([A-Za-z ]+):", entry, flags=re.MULTILINE)
        assert list(dict.fromkeys(fields)) == list(DIVERGENCE_FIELDS), entry_id
        assert f"\n- ID: `{entry_id}`\n" in entry, entry_id
```

- [ ] **Step 2: Run the test.** Run `uv run pytest tests/test_package.py -q`, and expect FAIL in `test_each_divergence_row_of_the_design_has_an_entry_with_each_field`, with "At index 13 diff: 'D-015' != 'D-014'". The planner ran this test from `tmp/claude/plan-docs/test_divergence_rows_scratch.py`.
- [ ] **Step 3: Write the entry `D-014`.** Put it between `D-013` and `D-015` of `docs/divergences.md`:

```markdown
## D-014: The documentation uses only functions that exist

- ID: `D-014`
- R location: the vignettes `vignettes/articles/SeedMatchR_miRNA_analysis.Rmd` (`:64`, `:76`, `:178`) and `vignettes/articles/SeedMatchR_build_PPI_network.Rmd` (`:30`).
- R behavior: The miRNA vignette calls `filter_deseq` (`:64`), `get_feature_seqs` (`:76`), and `de_fc_ecdf` (`:178`). The network vignette calls `get_network_matches` without `organism` (`:30`). An R probe with R 4.5.3 and SeedMatchR 2.0.0 gives the evidence:
  - `exists(name, envir = asNamespace("SeedMatchR"), inherits = FALSE)` is `FALSE` for `filter_deseq`, `de_fc_ecdf`, and `get_feature_seqs`. It is `TRUE` for `filter_res` and `deseq_fc_ecdf`.
  - `get_network_matches(c("A"), aggregator = "mean")` stops with `argument "organism" is missing, with no default`.
- New behavior: The examples of `README.md` use only the public names of seedmatch, and a test runs each Python example. Each Python name of `docs/r-map.md` exists.
- Reason: An agent copies the examples of the documentation. An example that calls an absent function stops with an error, and the agent cannot use it.
- Test: `test_each_python_name_of_the_r_map_exists_D014` in `tests/test_r_map.py`. The tests of `tests/test_readme.py` run each Python example of `README.md`.
```

- [ ] **Step 4: Name the golden file of the compensatory cases.** In `D-052`, `D-053`, `D-054`, and `D-055`, add this item after the last item of the R behavior. It goes before the line `- New behavior:`.

```markdown
  - Evidence: the cases and the notes of `tests/golden/noncanonical/compensatory.tsv`.
```

  The item obeys the form of `docs/divergences.md:123`. The planner found each case and each note of these entries in this file.
- [ ] **Step 5: Correct the claim of the design.** Replace `docs/design.md:592` ("A test runs each code example in `README.md` and in `docs/`.") with "A test runs each Python example in `README.md`. The example of "Public API example" is a sketch with placeholder paths, and no test runs it."
- [ ] **Step 6: Run the tests.** Run `uv run pytest tests/test_package.py tests/test_r_map.py -q`, and expect PASS. `test_each_divergence_id_is_unique_and_has_a_test` finds `D014` in the name of the test of Task 3. Then run `uv run ruff check`, `uv run ruff format --check`, and `uv run mypy`, and expect no finding.
- [ ] **Step 7: Run the check of the record.** Run `uv run python tmp/claude/plan-docs/check_divergences.py`, and expect `design IDs without an entry: none`. The lines "the R behavior names no golden, probe, or R line" for `D-052` to `D-055` do not occur. The lines with repeated fields, for example of `D-064`, stay, and they are permitted.
- [ ] **Step 8: Checkpoint.** The controller runs the full suite one time, after the last task.
