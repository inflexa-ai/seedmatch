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
- Other processes used 2 to 4 CPUs during each run, and the load average was 2.8 to 8.7.
- In a second set of runs on the same day, no time was more than 17% above the value of the table.

## The commands

Run the commands from the repository root. The R command reads the R library `tools/golden/.rlib`. Make this library first with `tools/golden/run.sh tools/golden/install.R`. `SEEDMATCH_TEST_REFS` is the directory of the `refs` tests. On Linux, use `/usr/bin/time -v` for the peak memory of the R process.

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

Each script makes three calls in one process. The Python script gives the peak resident memory of its process after each call. The R script gives the resident memory after each call, and `/usr/bin/time` gives the peak of the R process. Thus the row "Resident memory before the first call" is a peak for Python and the current value for R.

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
- The row "Peak above the memory before the first call" also holds a part that does not grow with the input. A probe on `tests/data/Test_seqs.fa` gives this part:
  - In R, the first call loads the SeedMatchR namespace, and the load adds about 111 MiB.
  - In Python, the first call adds about 33 MiB, and the second call adds none.
- The empty cache adds the memory of the compiler.
