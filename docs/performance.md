# Performance of seed_match_report

This document records the wall time and the peak memory of `seed_match_report` on the canonical 3' UTRs of `rnor7`, for Python and for R. `docs/design.md` ("Performance and bounds") sets the target: Python is not slower than R by more than a factor of 2.

## The calls

- Python: `seed_match_report(de, sequences, guide)`, with the default `ReportOptions()` and the default `workers=None`. `workers=None` gives a thread pool with one thread for each CPU.
- R: `SeedMatchR::SeedMatchReport(res, seqs, guide.seq)`, with the default arguments (`R/SeedMatchReport.R:46-55`). R uses one thread.
- The defaults of the two calls agree. The row `SeedMatchReport` of [`src/seedmatch/r-map.md`](../src/seedmatch/r-map.md) maps each argument.
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
- Python 3.12.13, seedmatch 0.1.0 (commit `07648b7`), Rust 1.100.0, PyO3 0.30.1, rust-numpy 0.30.0, numpy 2.5.2, and pandas 2.3.3.
- R 4.5.3, SeedMatchR 2.0.0 (commit `2fc0c4f`), Biostrings 2.78.0, dplyr 1.2.1, and gt 1.4.0.
- Other processes used the CPUs during the Python runs. The 15-minute load average after the runs was 9.7.

## The commands

Run the commands from the repository root. The R command reads the R library `tools/golden/.rlib`. Make this library first with `tools/golden/run.sh tools/golden/install.R`. `SEEDMATCH_TEST_REFS` is the directory of the `refs` tests. On Linux, use `/usr/bin/time -v` for the peak memory of the R process.

```bash
# The canonical 3' UTRs of rnor7, with the gene ID as the record ID.
uv run python tools/perf/utr3_fasta.py "$SEEDMATCH_TEST_REFS/rnor7" tmp/perf/rnor7_utr3.fa

# Python: --hashes gives the SHA-256 of the report and of the matches.
uv run python tools/perf/report.py tmp/perf/rnor7_utr3.fa --hashes
uv run python tools/perf/report.py tmp/perf/rnor7_utr3.fa --workers 1 --hashes

# Python: the time of each search of the full_search grids, and the SHA-256 of its ranges.
uv run python tools/perf/grid.py tmp/perf/rnor7_utr3.fa

# R: three calls in one session.
PERF_FASTA=tmp/perf/rnor7_utr3.fa /usr/bin/time -l tools/golden/run.sh tools/perf/report.R
```

Each report script makes three calls in one process. The Python script gives the peak resident memory of its process after each call. The R script gives the resident memory after each call, and `/usr/bin/time` gives the peak of the R process. Thus the row "Resident memory before the first call" is a peak for Python and the current value for R.

## The results

Each cell gives the values of two runs.

| Measure | Python, `workers=None` | Python, `workers=1` | R |
|---|---|---|---|
| First call | 0.19 s, 0.20 s | 1.19 s, 1.19 s | 9.59 s, 9.63 s |
| Later calls | 0.18 to 0.19 s | 1.19 to 1.22 s | 9.32 to 9.35 s |
| Resident memory before the first call | 108 MiB, 109 MiB | 107 MiB, 108 MiB | 296 MiB, 313 MiB |
| Peak resident memory of the process | 137 MiB, 140 MiB | 135 MiB, 136 MiB | 432 MiB, 443 MiB |
| Peak above the memory before the first call | 28 MiB, 32 MiB | 27 MiB, 29 MiB | 136 MiB, 130 MiB |

- R takes 47 to 54 times the Python time. Thus the target holds.
- With `workers=1`, R takes 7.6 to 8.1 times the Python time. Thus the target also holds on one thread.
- The first call and a later call take the same time, because maturin compiles the search kernels when it builds the wheel. `import seedmatch` takes 0.33 to 0.45 seconds.
- The row "Peak above the memory before the first call" also holds a part that does not grow with the input. In R, the first call loads the SeedMatchR namespace, and the load adds about 111 MiB.

## Before the change to the Rust kernels

Before commit `13f5cf2`, numba compiled the search kernels at run time. These runs used the same machine, the same input, and seedmatch commit `5e0783c` with numba 0.67.0 and llvmlite 0.49.0. The numba cache was full. The report and the matches have the same SHA-256 before and after the change.

| Measure | Python, `workers=None` | Python, `workers=1` |
|---|---|---|
| First call | 2.03 s, 2.03 s | 8.87 s, 8.88 s |
| Later calls | 1.87 to 2.26 s | 8.70 to 8.86 s |
| Resident memory before the first call | 143 MiB, 143 MiB | 142 MiB, 144 MiB |

An empty numba cache added about 0.9 seconds to the first call.

## The searches of full_search

`tools/perf/grid.py` times each search of `FullSearchOptions()` and of `FullSearchOptions(indels=True)`. Each time is the best of 3 calls of `count_matches` on one thread, after one call to warm up. The SHA-256 of the ranges of each search is the same for numba and for the Rust kernels.

The grid without indels takes 3.495 seconds with numba and 1.430 seconds with the Rust kernels. The grid with indels takes 9.835 seconds with numba and 1.440 seconds with the Rust kernels. The table gives each search of the grid with indels. A search with `mm0` is an exact search, thus it does not use the indel kernel.

| Search | numba | Rust | Factor |
|---|---|---|---|
| `Full.mm0` | 0.080 s | 0.050 s | 1.6 |
| `Full.mm1` | 0.189 s | 0.068 s | 2.8 |
| `Full.mm2` | 0.313 s | 0.067 s | 4.6 |
| `Full.mm3` | 0.498 s | 0.067 s | 7.4 |
| `Full.mm4` | 0.825 s | 0.067 s | 12.2 |
| `18mer.mm0` | 0.087 s | 0.049 s | 1.8 |
| `18mer.mm1` | 0.302 s | 0.067 s | 4.5 |
| `18mer.mm2` | 0.541 s | 0.067 s | 8.1 |
| `18mer.mm3` | 0.946 s | 0.067 s | 14.1 |
| `18mer.mm4` | 1.501 s | 0.067 s | 22.3 |
| `15mer.mm0` | 0.095 s | 0.049 s | 1.9 |
| `15mer.mm1` | 0.308 s | 0.067 s | 4.6 |
| `15mer.mm2` | 0.509 s | 0.067 s | 7.6 |
| `15mer.mm3` | 0.823 s | 0.067 s | 12.2 |
| `15mer.mm4` | 1.287 s | 0.077 s | 16.8 |
| `mer8.mm0` | 0.084 s | 0.048 s | 1.7 |
| `mer8.mm1` | 0.289 s | 0.068 s | 4.3 |
| `mer7m8.mm0` | 0.085 s | 0.049 s | 1.8 |
| `mer7m8.mm1` | 0.284 s | 0.069 s | 4.1 |
| `mer7A1.mm0` | 0.093 s | 0.048 s | 1.9 |
| `mer7A1.mm1` | 0.303 s | 0.070 s | 4.3 |
| `mer6.mm0` | 0.093 s | 0.048 s | 1.9 |
| `mer6.mm1` | 0.300 s | 0.076 s | 3.9 |

The time of the Rust kernels changes little with the edit distance. The filter of the indel search reads each letter of a sequence one time. Then the exact DP runs only at the positions that can start a match.
