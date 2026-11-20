# seedmatch: design

## Purpose

`seedmatch` is a Python 3.12 port of the R package SeedMatchR 2.0.0. SeedMatchR finds the seed matches of an siRNA guide strand in transcript sequences. It also does a test of whether the genes with seed matches have a shift in differential expression (DE) results.

Agents in the inflexa sandbox write Python code that imports `seedmatch`. The sandbox has no network egress. The reference store (`/mnt/refs`) and the package store (`/mnt/libs`) are read-only in the sandbox.

The source of truth for behavior is the SeedMatchR repository at commit `2fc0c4f` in the sibling directory `../SeedMatchR`. Each `R/...` reference in this document is a path in that repository.

## Goals

- Make a better version of SeedMatchR in Python, not only a copy.
- Port each SeedMatchR function to Python.
- Give output that is identical to SeedMatchR on the core paths, except for each recorded divergence. Refer to "Parity contract".
- Improve the design: typed options, no hidden state, and computation that is separate from I/O and from plots.
- Read each external resource from a path that the caller gives.
- Give a separate fetch module and command that download each resource from its upstream.
- Fix each defect of the R code, also on the core paths.
- Improve a behavior when a better behavior is possible.
- Record each change of behavior in `docs/divergences.md`.
- Do a test of each behavior with tests that can find a real defect: R goldens, oracles, and invariants.

## Non-goals

- A compatibility layer with the R function names. The documentation gives a map from each R name to its Python name.
- GRanges and IRanges objects. Results are pandas DataFrames.
- Plots that are identical to the ggplot2 plots at the pixel level.
- Publication to PyPI, an entry in the inflexa package manifest, an inflexa skill pack, and entries in the inflexa reference catalog. These are later work in other repositories.

## Anti-goals

- Network access in a core module, also as a fallback.
- A write next to a resource file. The stores are read-only, thus such a write fails with `EACCES` in the sandbox.
- A divergence from R that `docs/divergences.md` does not record.
- A test that only restates the code. Examples: a test that a dataclass keeps a field, a test of a default value, or a test that a mock received a call.

## Constraints

- Python 3.12. `pyproject.toml` sets `requires-python = ">=3.12"`, and `.python-version` is `3.12`.
- uv manages the project, with the `uv_build` backend and the `src/` layout.
- ruff does the lint and the format.
- The runtime dependencies come only from the inflexa package-store lock (`../inflexa/images/package-store/lock.amd64.json`):
  - required: `numpy`, `pandas`, `scipy`, `numba`
  - extra `plot`: `matplotlib`
- The genome readers use only the standard library (`mmap`, `zlib`, `sqlite3`) and numpy. Thus a BGZF FASTA genome needs no extra.
- A new runtime dependency needs an entry in the inflexa manifest. Thus, prevent it.
- The wheel is pure Python. The inflexa build images have no Rust toolchain. Numba compiles the search kernels at run time.
- The kernels use `@njit(cache=True)`. In the sandbox, the package directory is read-only. Thus numba writes its cache to `NUMBA_CACHE_DIR` or to its user cache. An inflexa warm script can fill this cache later.

## Parity contract

### Meaning of "identical"

- The rows are the same, in the same order.
- Integer and string values are equal.
- A float value that R rounds (for example, a Wasserstein area) is equal.
- Another float value agrees with R within a relative tolerance of 1e-9.

### Core paths

These paths must be identical to R. A recorded divergence replaces the R behavior.

| R function | Python replacement |
|---|---|
| `get_seed`, `.get_canonical_seed`, `.get_custom_sequence` | `Guide.target` |
| `Biostrings::vcountPattern`, `vmatchPattern`, `vmatchPattern2`, as SeedMatchR uses them | `count_matches`, `find_matches` |
| `SeedMatchR`, each `res.format` | `annotate_seed_counts`, `find_matches` |
| `full_search`, `full_search_ranges`, `.get_best_match` | `full_search`, `full_search_ranges`, `best_match` |
| `SeedMatchReport` (`$report` and `$matches`) | `seed_match_report` |
| `filter_res` | `filter_de` |
| `ecdf_stat_test`, `deseq_fc_ecdf` (statistics and groups) | `ecdf_test`, `assign_groups` |
| `wass_dist` (the areas) | `wasserstein_areas` |
| `load_annotations`, `get_3UTR_granges`, `.get_longest_utr`, `.reduce_features`, `build_annotation_filter` | `Reference.features`, `Reference.sequences`, `TranscriptFilter` |
| `sequence_kmer_counts` | `kmer_counts` |
| `plotAlignment` with `output.type = "data.frame"` | `alignment_text` |
| `download_parse_file`, `load_example_data` | `read_deseq2_table` |

### Not identical by design

- The DTS p-value comes from permutations. The statistic is identical. The p-value agrees with R only within random error. The caller gives the random generator.
- Plots use matplotlib.
- The report table uses a pandas `Styler` in HTML, not `gt`.
- Each entry in `docs/divergences.md` changes the R behavior.

## Architecture

### Packages

| Module | Responsibility |
|---|---|
| `seedmatch.sequence` | IUPAC alphabet as 4-bit masks. `SequenceSet`: names, one `uint8` buffer, and offsets. FASTA read and write. |
| `seedmatch.seed` | `Guide`, `Seed` (enum of the canonical seeds, plus a custom window), and the target derivation. |
| `seedmatch.search` | Numba kernels and the search API. |
| `seedmatch.annotation` | EnsDb reader, GTF reader, 2bit and FASTA genome readers (plain and BGZF), filters, feature derivation, sequence extraction, `Reference`. |
| `seedmatch.de` | DE table schema, the reader for the GEO DESeq2 files, the filter. |
| `seedmatch.stats` | KS and Wilcoxon with the R rules, DTS, Wasserstein areas, ECDF test, PR and ROC curves. |
| `seedmatch.report` | Full search, best match, report. |
| `seedmatch.noncanonical` | Wobble and bulge search, 3' compensatory search. |
| `seedmatch.network` | TRRUST message passing and the STRING gene network. |
| `seedmatch.kmers`, `seedmatch.alignment`, `seedmatch.distances`, `seedmatch.evalue` | The smaller analyses. |
| `seedmatch.plot` | Figures from results. Extra `plot`. |
| `seedmatch.resources` | The registry of the upstream resources: identifiers, URLs, file names, digests. No I/O. |
| `seedmatch.fetch` | Downloads, and the command `seedmatch-fetch`. |

### Rules

- Only `seedmatch.fetch` can import a network module. The ruff rule `TID251` (banned API) bans these modules in all other code: `urllib.request`, `http.client`, `socket`, `ssl`, `urllib3`, `requests`, `httpx`. A per-file ignore lets `src/seedmatch/fetch/**` use them.
- The same rule bans `seedmatch.fetch` in each module outside `seedmatch.fetch`.
- `import seedmatch` does not import matplotlib or `seedmatch.fetch`.
- A computation function returns data. A plot function takes a result and returns a `matplotlib.figure.Figure`. No function writes a file unless it has an explicit output path parameter.
- A reader reads only the part of a file that the call needs:
  - The 2bit reader uses `mmap` and decodes only the requested ranges.
  - The EnsDb reader puts the filter in the SQL `WHERE` clause.
- Options are frozen dataclasses. Modes are enums. No function takes a string mode flag.
- All errors derive from one base class, `SeedMatchError`. Each message names the bad value and the expected value.
- A function validates untrusted input at the boundary only: the guide alphabet and length, the seed window, and the DE table columns.

### Public API example

The final names are part of the approval of this document.

```python
import seedmatch as sm

ref = sm.Reference.from_paths(
    ensdb="/mnt/refs/managed/<id>/<version>/EnsDb.Rnorvegicus.v109.sqlite",
    genome="/mnt/refs/managed/<id>/<version>/Rattus_norvegicus.mRatBN7.2.dna_sm.toplevel.2bit",
)
utr3 = ref.sequences(sm.FeatureType.UTR3, sm.TranscriptFilter(canonical=True))

guide = sm.Guide("UUAUAGAGCAAGAACACUGUUUU")
target = guide.target(sm.Seed.MER7M8)

de = sm.filter_de(sm.read_deseq2_table(path), fdr_cutoff=1.0, fc_cutoff=0.0, basemean_cutoff=30)
de = sm.annotate_seed_counts(de, utr3, target, column="mer7m8")

seed_genes = de.loc[de["mer7m8"] >= 1, "gene_id"]
background = de.loc[de["mer7m8"] == 0, "gene_id"]
result = sm.ecdf_test(de, target=seed_genes, background=background,
                      test=sm.EcdfTest.KS, alternative=sm.Alternative.GREATER)

report = sm.seed_match_report(de, utr3, guide)
```

## Search rules

Biostrings 2.78.0 shows these rules. Each rule is normative for the port.

1. Matches overlap. `countPattern("AA", "AAAA")` gives 3.
2. If `max_mismatch > 0` and indels are off, the search obeys these rules:
   - A match can extend past either end of a sequence. Each position outside the sequence counts as a mismatch. Thus, the reported coordinates can be less than 1 or greater than the sequence length.
   - Let m be the pattern length and n be the sequence length. The search tries only the starts from 2 − m to n. Thus, if the sequence is not empty, at least one pattern letter is in the sequence.
   - This bound changes the result only if m ≤ `max_mismatch`. For example, `matchPattern("AC", DNAString("GG"), max.mismatch = 3)` gives 0–1, 1–2, and 2–3. An empty sequence still gets the starts from 2 − m to 0. For example, `matchPattern("AC", DNAString(""), max.mismatch = 2)` gives 0–1.
   - Biostrings sets this bound in `match_naive_inexact` (`src/match_pattern.c`).
3. The `fixed` argument is a pair (pattern, subject). If a side is not fixed, its IUPAC codes are ambiguous. For example, with `fixed = FALSE`, an `N` in the subject matches each base. An R script makes the complete table of letter compatibility for each `fixed` pair. The port uses this table.
4. If `with_indels` is on and `0 < max_mismatch < m`, the search does these steps. Biostrings does them in `_match_pattern_indels` (`src/match_pattern_indels.c`).
   1. Each sequence position j can give one candidate that starts at j. Find the first pattern letter that matches the sequence letter at j. Let i0 be the number of pattern letters before it. If no pattern letter matches, or if i0 > `max_mismatch`, position j gives no candidate.
   2. Align the pattern letters after this letter to a substring that starts at j + 1. The substring can have each width, also 0. Find the least edit distance d, and the shortest width with this distance.
   3. The cost of the candidate is i0 + d. The candidate starts at j and ends at the end of this substring. If the cost is greater than `max_mismatch`, discard the candidate.
   4. Examine the candidates in the order of their starts, and hold one candidate. If the end of the next candidate is not after the held end, compare the two costs. If the cost of the next candidate is not greater, the next candidate replaces the held candidate. Else, discard the next candidate. If the next candidate ends after the held candidate, report the held candidate, and then hold the next candidate. At the end, report the held candidate.

   Step 1 uses only the first matching pattern letter. A later matching letter cannot give a lower cost, because the alignment of the first letter can delete the letters between them. Thus, if two letters give the same cost, the first letter sets the width.

   The matches are in ascending order of start. No match contains a different match, thus the ends are also in ascending order.

   If `with_indels` is on and `max_mismatch ≥ m`, the search gives the result of rule 2 with `min_mismatch = 0`. Biostrings selects this path in `_match_pattern_XString` (`src/match_pattern.c`).

   The documentation of `matchPattern` defines the result as the "best local matches", but Biostrings does not obey this definition. For example, a reported range can contain a shorter substring with the same edit distance. The port obeys Biostrings, not the documentation.
5. If `max_mismatch = 0`, the search ignores `with_indels`.
6. `min_mismatch` applies only when indels are off.
7. Coordinates are 1-based and inclusive, as in R.

SeedMatchR gets its ranges from `vmatchPattern2`, which always uses the indel algorithm (`R/vmatchPattern2.R:43`). Thus, for a pattern longer than `max.mismatch`, R ignores `with.indels` and `min.mismatch` for ranges, and its ranges disagree with its counts. An R probe shows this: with `with.indels = FALSE`, `vmatchPattern2` returns a deletion match that `vmatchPattern` does not find. If the pattern is not longer than `max.mismatch`, the C code sends the indel algorithm to the mismatch search, and this search applies `min.mismatch` (`_match_pattern_XString`, `src/match_pattern.c`). The port obeys rules 5 and 6 for counts and for ranges. This is a divergence.

The default `SeedMatchReport` uses indels with the edit distances 0 to 4 for `Full`, `18mer`, and `15mer` (`R/SeedMatchReport.R:50`, `R/full_search.R:34`). Thus rules 2 and 4 change the default report. `SeedMatchReport` passes `fixed = TRUE` (`R/SeedMatchReport.R:55`, `:64`), but a direct call of `full_search` uses `fixed = FALSE` (`R/full_search.R:33`). Thus rule 3 changes only the direct call.

### Licenses

- `seedmatch` uses the MIT License, the same license as SeedMatchR.
- SeedMatchR has the copyright of Eli Lilly and Company. The port is a work from SeedMatchR, thus `LICENSE` gives this copyright next to the copyright of `seedmatch`.
- Biostrings uses Artistic-2.0. R `stats` uses GPL. PRROC and twosamples are R packages with their own licenses.
- The port writes each algorithm from its documentation or its published paper. The R output then makes sure that the result is the same.
- The builder can read the R or C source to settle a behavior that the documentation leaves open. The builder must not translate GPL source line by line.

## Statistics rules

### KS test

- R uses the exact p-value only if `n_x * n_y < 10000`. scipy `method="auto"` uses a different rule. Thus the port uses the R rule.
- For the exact p-value with ties, R 4.3 and later use the algorithm of Schröer and Trenkler (1995). The port computes the same upper tail with the recursion of Viehmann (2021). A small p-value keeps its relative precision.
- The statistic is d / (n_x · n_y). The integer d is the largest gap `i · n_y − j · n_x` at a checkpoint. `LESS` uses the negative gap, and the two-sided test uses the absolute gap. Here i and j count the x values and the y values up to the checkpoint. Thus d is exact, and D never exceeds 1. A checkpoint is the end of a run of equal pooled values.
- The one-sided asymptotic p-value is exp(-2 · en · D²), with en = n_x · n_y / (n_x + n_y). scipy adds the Hodges correction, thus the port does not use the scipy value.
- The two-sided asymptotic p-value is the limit Kolmogorov distribution, `scipy.special.kolmogorov(sqrt(en) · D)`.
- The asymptotic formulas fix numeric defects of R, which give a wrong p-value or a p-value of 0 (`D-032`).
- R warns if ties make the asymptotic p-value approximate. The port gives no warning, because `method` shows it (`D-037`).
- On the test data, the scipy exact p-value agrees with R: D = 0.51035, p = 0.0001037731.

### Wilcoxon rank-sum test

- R uses the exact p-value if `n_x < 50`, `n_y < 50`, and there are no ties. Otherwise, R uses the normal approximation with the tie correction and the continuity correction. The port uses the R rule.
- The statistic W is the rank sum of x minus n_x(n_x + 1)/2. A tie gets the average rank.
- The exact p-value comes from the count of the subsets of ranks for each W. The port counts with a table that has one row for each subset size and one column for each rank sum. It adds the ranks to the table one at a time. The cost is bounded, because each size is below 50.
- The counts pass 2^63, thus the table holds float64 values. Each step adds non-negative numbers, thus the relative error stays small.
- `LESS` gives P(W ≤ W_obs), and `GREATER` gives P(W ≥ W_obs). Each tail is a sum of counts, not 1 minus the other tail. Thus a small p-value keeps its relative precision.
- `TWO_SIDED` gives min(1, 2 · p). Here p is P(W ≥ W_obs) if W_obs > n_x · n_y / 2, and else P(W ≤ W_obs).
- The normal approximation starts with `z = W - n_x * n_y / 2`.
- The standard deviation is `SIGMA = sqrt((n_x * n_y / 12) * ((N + 1) - sum(t^3 - t) / (N * (N - 1))))`. Here N = n_x + n_y, and t is the size of a group of ties.
- The continuity correction is `0.5 * sign(z)` for `TWO_SIDED`, 0.5 for `GREATER`, and -0.5 for `LESS`. The corrected value is `z' = (z - correction) / SIGMA`.
- Let Phi be the normal CDF. `LESS` gives Phi(z'), and `GREATER` gives Phi(-z'). `TWO_SIDED` gives 2 * min(Phi(z'), Phi(-z')).
- The port computes Phi(x) as `exp(log_ndtr(x))`, with `scipy.special.log_ndtr`. `ndtr(x)` gives 0 below about -37.5, and `1 - ndtr(x)` gives 0 for a small tail.
- Thus a small p-value keeps its subnormal float down to 5e-324, as `pnorm` of R does.
- If each value is equal, SIGMA is 0. R gives NaN for `TWO_SIDED`. The port gives the p-value 1 for each alternative (`D-033`).
- R warns if ties make the p-value approximate. The port gives no warning, because `method` shows it (`D-037`).
- On the test data, `scipy.stats.mannwhitneyu(method="asymptotic", use_continuity=True)` agrees with R: W = 2150, p = 0.9996.

### DTS test

- The statistic is a sum over the gaps between the sorted pooled values.
- Each term is |E − F| / sqrt(2 · G · (1 − G) / n), multiplied by the width of the gap.
- E and F are the ECDFs of the two samples. G is the pooled ECDF. n is the pooled size. Each ECDF is at the left value of the gap.
- A gap of width 0 adds nothing to the sum. Thus only the end of a run of equal values counts.
- A Python test of this formula gives 74.91157, the same as the R snapshot.
- The port computes E − F from integer counts: (c · n − k · n_x) / (n_x · n_y). Here k values lie at or below the gap, and c of them come from the first sample.
- A permutation chooses n_x of the pooled sorted positions for the first sample. This is the method of twosamples.
- The p-value is the fraction of `n_boot` permutation statistics that are not less than the observed statistic (default 2000, as in twosamples). If the fraction is 0, the p-value is 1 / (2 · `n_boot`).
- The observed statistic and each permutation statistic come from the same code, and each row adds its terms in the same order. Thus a permutation with the same counts as the observed arrangement gives the same value, and it counts.
- The caller gives a `numpy.random.Generator` or a non-negative int seed.
- The port draws the permutations in blocks of 64, and numpy does the work of each block. Thus the memory is O(64 · n) and does not grow with `n_boot`. The time is O(`n_boot` · n).
- An infinite value and an empty sample raise `InputError`. R gives the statistic `Inf`, or a statistic that does not compare two samples, and the p-value 1 (`D-034`).

### Wasserstein areas

The port copies the algorithm of `R/wass_dist.R:61-81`:

- Evaluate each ECDF at its own data points. The value at a point is the count of values at or below it, divided by the sample size.
- Merge equal values to one point, as `stats::approx` does with `ties = mean`. The ECDF has the same value at equal points, thus the mean is that value. R warns if a sample has equal values. The port gives no warning (`D-037`).
- Interpolate linearly between the merged points with the expression of `approx`: y_i + (y_j − y_i) · ((v − x_i) / (x_j − x_i)). The value is 0 below the first point and 1 above the last point (`yleft = 0`, `yright = 1`).
- Use a grid of `n_interp` points from the minimum to the maximum of the pooled values, as `seq(from, to, length.out)` of R does. Point k is min + k · ((max − min) / (`n_interp` − 1)). The first point is the minimum and the last point is the maximum.
- Let diff be the seed ECDF minus the non-seed ECDF. Integrate max(diff, 0), max(−diff, 0), and |diff| with the trapezoid rule of `caTools::trapz`. The areas are `negative_censored`, `positive_censored`, and `absolute`.
- Round each area to 3 decimals with the R rounding rule.

`numpy.interp` does not merge equal x values, and it uses another expression for the line. Thus the port does not call it.

The R rounding rule (`r_round` in `seedmatch._rmath`) applies in R 4.0.0 and later. Python `round` is different:

- Let p = 10^d. The two candidates are xd = floor(x · p) / p and xu = ceil(x · p) / p, in float64.
- The result is xu if xu − x < x − xd. If the two differences are equal, the result is the candidate whose integer numerator is even. The numerator is floor(x · p) for xd, and ceil(x · p) for xu. Thus xu wins if floor(x · p) is odd. Otherwise, the result is xd.
- Python `round` uses the exact binary value of x. For x = 0.0025 and d = 3, both float differences are 0.0005. Thus R gives 0.002, and Python gives 0.003.
- NaN, an infinite value, and 0 stay the same. R rounds a negative x as −x, and then puts the sign back.
- x stays the same if log10(2) · (0.5 + logb(x)) + d is greater than 15. Here logb(x) is the binary exponent of x. R uses this estimate, not log10(x).
- d = 0 gives the nearest integer, and a tie goes to the even integer.
- `r_round` accepts d from 0 to 22, because 10^d is exact in float64 up to 10^22. A different d raises `InputError`.

These inputs raise `InputError`. The message gives the sample and the bad value:

- A sample with fewer than 10 values. R stops with the same limit (`R/wass_dist.R:51-56`).
- A sample with only one distinct value. R stops in `stats::approx` with the message `need at least two non-NA values to interpolate`.
- NaN or an infinite value. R stops in `seq` with the message `'from' must be a finite number` for NaN and −Inf. For +Inf, R stops with the message `'to' must be a finite number`.
- An `n_interp` below 2 (`D-038`).
- A pooled range above 8.988465674311579e307, half of the largest float (`D-038`). The sum of the products of the trapezoid rule can be twice the range, and a larger sum overflows.

Let n be the size of a sample. The time is O(n log n + `n_interp` log n), because `np.unique` sorts each sample and `np.searchsorted` searches for each grid point. The port processes the grid in chunks of 2^16 trapezoids, and numpy does the work of each chunk. Thus the memory does not grow with `n_interp`. For two samples of 20000 values, `n_interp` = 10^7 takes 0.4 seconds and 6.9 MB.

The unrounded areas of R depend on the platform. R 4.5.3 on arm64 macOS uses a fused multiply-add in `approx` and in the product `%*%` of `trapz`. R on x86-64 uses no fused multiply-add. Thus the last bits of an unrounded area differ between the two platforms. The port uses the plain IEEE expression of `approx`. It adds the products of `trapz` in order, with `numpy.cumsum`, because `numpy.sum` adds in another order. A chunk adds the sum of the earlier chunks to its first product. Thus the order is the same as for one sum over the grid. On the goldens, the unrounded areas differ from the arm64 values by a relative 5e-16 at most. A rounded area can differ from R only if its unrounded area lies within about 1e-15 of a rounding boundary. The nearest golden area is 1.3e-5 from a boundary. Thus each rounded golden area is equal.

## Annotation

### EnsDb and 2bit (the parity path)

- The EnsDb is an SQLite file with the ensembldb schema. The reader uses the standard `sqlite3` module. It opens the file read-only and immutable (`file:...?mode=ro&immutable=1`, URI open). Thus the reader can open a file in a read-only store.
- The reader needs the tables `gene`, `tx`, `tx2exon`, and `exon`. A filter that needs an optional column (`tx_is_canonical`, `tx_support_level`) or the `entrezgene` table raises `InputError` if the file does not have it. The message gives the file and the column.
- The EnsDb `chromosome` table gives the length of each contig. If the genome gives a different length for a contig, `Reference` raises `InputError`. Thus a genome of a different build causes an error.
- `TranscriptFilter` has these fields. Their semantics come from `build_annotation_filter` (`R/load_annotations.R:242-339`) and from the SQL that ensembldb 2.34.0 runs:
  - `standard_chromosomes`: `seq_name` is one of `1` to `22`, `X`, and `Y`, for each species, or one of these names with the prefix `chr`
  - `require_tx_id`: `tx_id != ''`
  - `protein_coding`: the gene biotype and the transcript biotype are `protein_coding`
  - `canonical`: `tx_is_canonical = 1`
  - `support_level`: `tx_support_level` is equal to the value, from 1 to 5, as `TxSupportLevelFilter` does
  - `gene_ids`, `tx_ids`, `symbols` (`gene_name`), `entrez_ids`: the value is one of the given values
  - `tx_biotypes`, `gene_biotypes` (these replace the free `add.filter`)
- The defaults are the defaults of `load_annotations`. Thus `canonical` is true, but the default of `build_annotation_filter` is false.
- A field with a set of values takes a sequence. A `str` raises `InputError`, because it is one ID and not a set. An empty sequence raises `InputError`, because it keeps no transcript.
- The feature derivation copies `ensembldb`. The EnsDb CDS includes the stop codon. The rnor7 probe shows this: the 3' UTR starts at `tx_cds_seq_end + 1`.
  - On the plus strand, the 3' UTR is the part of each exon after `tx_cds_seq_end`. The 5' UTR is the part before `tx_cds_seq_start`. On the minus strand, the two parts change places.
  - The CDS is the part of each exon from `tx_cds_seq_start` to `tx_cds_seq_end`.
  - The 5' UTR, 3' UTR, and CDS results have each transcript with a CDS. The exon result has each transcript.
  - R gives a range of width 0 to a transcript without the feature, for example to 4168 of the 23036 canonical rat 3' UTRs. The port gives such a transcript no range. The sequence is empty in the two cases, thus `min_width` removes it.
- The order of each result is the R order, because it sets the row order of the searches. The SQL that ensembldb runs shows the order:
  - Each transcript: the transcripts in the order of `tx_id`. The ranges of a transcript come in the order of `exon_idx`, thus in transcript order.
  - `Reference.transcripts`: the order of `seq_name`, then of `tx_seq_start`. ensembldb sorts in SQL on these two columns only, thus a tie keeps the scan order of SQLite. The port sorts on 4 keys explicitly: `seq_name`, `tx_seq_start`, the row order of the `gene` table, and the row order of the `tx` table (`D-030`).
  - The scan order of SQLite depends on the query plan. With the filter on the standard chromosomes, SQLite reads the `gene` table first. Without it, SQLite reads the `tx` table first. An EnsDb from ensembldb has the `tx` rows of each gene in the order of the `gene` rows, thus the two scan orders agree. The rnor7 probe shows this for each tie. The fixture has a tie of two genes whose row orders disagree.
- The sequence extraction copies `GenomicFeatures::extractTranscriptSeqs`. It joins the ranges in the order of the list. On the minus strand, it uses the reverse complement of each range.
- `PerGene` selects the transcripts of a gene:
  - `EACH_TRANSCRIPT` is the default. The result has one sequence for each transcript.
  - `LONGEST` copies `.get_longest_utr`. The genes come in the order of `gene_id`. For each gene, the result has the first transcript with the maximum length, in the order of `Reference.transcripts`. Only a transcript in the feature result is a candidate (`D-028`).
  - `UNION` copies `.reduce_features`. The genes come in the order of `gene_id`. The union merges the ranges that overlap or touch. The ranges of the union come in transcript order (`D-027`).
  - R applies `longest.utr` and `reduce.features` only to the 3' UTR. The port applies `PerGene` to each feature type (`D-029`).
- `min_width` removes each sequence that is shorter than the width. The default is 8, the default of `min.feature.width`.
- Each sequence keeps its transcript ID and its gene ID as two fields of `SequenceSet`. With `UNION`, the ID of a sequence is its gene ID. R replaces the names with gene IDs, thus two sequences can have the same name. With `reduce.features`, R gives the name `NA` to each sequence (`R/load_annotations.R:103-107`).
- A gene-level result uses an explicit rule. A seed count for a gene is the maximum over its sequences. The best match of a gene is the best over its sequences. R uses the first sequence for the counts (`R/SeedMatchR.R:207`), but each sequence for the report (`R/SeedMatchReport.R:85-99`). These are divergences.
- The rat 2bit file is soft-masked. Biostrings changes lowercase letters to uppercase, thus the reader also does this.
- The 2bit reader supports the UCSC format, version 0 and version 1 (64-bit offsets), in each byte order. It applies the N blocks.

### GTF and FASTA (the catalog path)

- `Reference.from_gtf(gtf, genome)` reads an Ensembl or GENCODE GTF, plain or compressed. The genome is a 2bit file or a FASTA file.
- The GTF reader streams the file. It writes the transcripts and the exons into an in-memory SQLite database with the EnsDb tables. Thus the GTF path uses the same filters, feature derivation, and order rules as the EnsDb path. The memory is linear in the number of exon lines.
- The row order of the in-memory tables is the line order of the GTF.
- An Ensembl GTF CDS line excludes the stop codon, but an EnsDb CDS includes it. Thus the CDS of a transcript goes from the first base to the last base of its `CDS` and `stop_codon` lines.
- The GTF reader raises `InputError` for a transcript with lines that do not agree (two exons that overlap, an exon with a contig or a strand different from the transcript line, or a `CDS` or `stop_codon` line without a transcript line).
- The canonical flag comes from the `Ensembl_canonical` tag. The support level is the integer at the start of `transcript_support_level`. If the value starts with `NA`, the transcript has no support level. If the value starts with neither an integer nor `NA`, the reader raises `InputError`. The biotypes come from `gene_biotype` and `transcript_biotype`, or from the GENCODE attributes `gene_type` and `transcript_type`.
- A GTF has no Entrez IDs, thus a filter on `entrez_ids` raises `InputError`.
- The IDs stay as the GTF gives them. Thus a GENCODE ID keeps its version suffix.
- The FASTA genome reader reads a plain FASTA file with its `.fai` index, or a BGZF file with its `.fai` and `.gzi` indexes. It uses only `zlib` and the file functions of Python. It decompresses only the blocks of the requested range. A gzip file that is not BGZF raises `InputError`.
- GENCODE contig names have the prefix `chr`, but the Ensembl FASTA names do not. For each contig, the extraction uses the first genome name that exists: the same name, the name without `chr`, or the name with `chr`. `chrM` and `MT` are the same contig. If no contig name of the annotation agrees with a genome contig, `Reference` raises `InputError`.
- A GTF has no contig lengths. Thus, on this path, `Reference` does not compare the contig lengths with the genome.
- This path makes no parity claim, because the catalog has a different Ensembl release from SeedMatchR. A consistency test compares this path with the EnsDb path on one Ensembl release.

### Reference objects

- `Reference.from_paths(ensdb=..., genome=...)` takes explicit paths. The genome is a 2bit file or a FASTA file. The first bytes of the file show the format.
- `Reference.from_directory(path)` finds exactly one EnsDb (`*.sqlite`) and one genome (`*.2bit`, or `*.fa.gz` with its indexes). It stops with an error if a file is absent or if there is more than one candidate.
- This fits the inflexa catalog layout `managed/<id>/<version>/`, which keeps the upstream file names.
- `Reference.transcripts(transcript_filter)` gives a `DataFrame` of the transcripts, as `ensembldb::transcripts` does. `Reference.features(feature, transcript_filter, per_gene=...)` gives the ranges of each feature. `Reference.sequences(feature, transcript_filter, per_gene=..., min_width=8)` gives a `SequenceSet`.
- A `Reference` holds an open database and an open genome. `close()` and the `with` statement release them.

## Resources and fetch

### Registry

`seedmatch.resources` holds the upstream data that the fetch module downloads. It does no I/O.

| Build | 2bit AnnotationHub ID | EnsDb AnnotationHub ID |
|---|---|---|
| `hg38` | AH106283 | AH109606 |
| `mm39` | AH106440 | AH109655 |
| `mm10` | AH88475 | AH89211 |
| `rnor6` | AH93578 | AH95846 |
| `rnor7` | AH106786 | AH109732 |
| `rnor7.113` | AH106786 | AH119437 |

- The download URL of an AnnotationHub resource is `location_prefix + rdatapath`, from `https://annotationhub.bioconductor.org/ahid/<id>`. For example, AH109732 is `https://mghp.osn.xsede.org/bir190004-bucket01/AnnotationHub/AHEnsDbs/v109/EnsDb.Rnorvegicus.v109.sqlite`.
- The registry pins the URL, the file name, the size, and the SHA-256 digest. The builder records the size and the digest after the first download.
- The other resources come from the URLs in the R code:
  - the DESeq2 tables of GEO GSE184929 (`R/get_data.R:12-21`)
  - the miRDB v6.0 predictions (`R/get_data.R:28`)
  - the TRRUST human and mouse tables (`R/get_network_matches.R:13-14`)
  - the STRING protein links for a species, with a pinned STRING version

### Fetch behavior

- The fetch module uses `urllib.request` from the standard library. It adds no dependency.
- It streams each file in chunks to `<dest>/<name>.part`. Then it makes sure that the size and the SHA-256 digest agree with the registry. Then it renames the file to its final name.
- If the final file is present and its digest agrees, it does not download the file again.
- It refuses a URL that is not `https`, and a redirect to a URL that is not `https`.
- A failed or interrupted download leaves no final file.
- The command `seedmatch-fetch` has a subcommand for each resource group. Examples:
  - `seedmatch-fetch reference rnor7 DEST`
  - `seedmatch-fetch example sirna DEST`
  - `seedmatch-fetch trrust human DEST`
  - `seedmatch-fetch string 10116 DEST`
  - `seedmatch-fetch mirdb DEST`
  - `seedmatch-fetch list`

## Divergences (defect fixes)

`docs/divergences.md` is the record of each change of behavior from R: each defect fix and each improvement. It is part of the repository and of the documentation. Each entry has these fields:

- a stable ID, for example `D-001`
- the R function and the R location
- the R behavior, with the evidence (an R golden or an R probe)
- the new behavior
- the reason
- the test that keeps the new behavior

A change of behavior that has no entry is a defect of the port. The builder adds an entry for each new divergence that the work finds. The table below gives the known entries. The IDs `D-001` to `D-014` follow the order of the first rows. Each row after the `Documentation` row has the next ID, from `D-015`.

| R function | R defect | Corrected behavior |
|---|---|---|
| `vmatchPattern2` (used by `SeedMatchR`, `full_search`, `full_search_ranges`) | Always uses the indel algorithm (`R/vmatchPattern2.R:43`). For a pattern longer than `max.mismatch`, it ignores `with.indels` and `min.mismatch`, thus its ranges disagree with the counts. | Ranges obey the same rules as counts. `full_search_ranges` gives no duplicate rows for the two indel options. |
| `load_annotations` with `reduce.features` | Each sequence gets the name `NA` (`R/load_annotations.R:103-107`). | Each sequence keeps its gene ID. |
| `.get_counts_column`, `SeedMatchReport` | With more than one sequence for a gene, the counts use the first sequence (`R/SeedMatchR.R:207`), but the report counts each sequence (`R/SeedMatchReport.R:85-99`). | Gene-level results use the maximum count and the best match over the sequences of the gene. The report counts genes. |
| `get_compensatory_sequences` | Always returns 0 rows: the empty-result return runs before the search (`R/get_compensatory_sequences.R:51-60`). | The search runs, and the function returns the found regions. |
| `run_wobble_bulges_pipeline` | Stops with `object 'bulge_allowance' not found` (`R/run_wobble_bulges_pipeline.R:137`). It swaps two arguments (`:12`, `:115`), uses only sets of exactly 4 U positions (`:63`), and sums the counts of variants that overlap. | A guide U can pair with a target G: one IUPAC pattern covers this at each U position. A bulge is 1 to `bulges_lim` extra target bases between two adjacent seed positions. No bulge occurs next to the `nt` seed positions from g2 (default g2 to g4), as the R documentation states. The result counts each distinct start position once for each transcript. |
| `benchmark_preds` | `padj <= padj` ignores the threshold (`R/benchmark.R:9`). `roc.curve` gets the labels as `scores.class1` (`:22-24`). | The threshold applies. The ROC curve uses the labels as weights, as the PR curve does. |
| `message_passing` | `mean(a, b)` uses `b` as `trim` (`R/get_network_matches.R:69-70`). `length(neighbors != 0)` (`:46`). `features_act` reads `features_rep` (`:70`, `:95`). | TRRUST edges go from a regulator to a target. The score of a gene comes from its regulators. The type of each edge selects the repression channel or the activation channel. An `Unknown` edge is in no channel. The score is (mean of the seed flags of the regulators + own seed flag) / 2. `degree_importance` divides each regulator flag by the out-degree of the regulator. A STRING network has no edge type, thus it gets one channel with the mean of the neighbors. |
| `generate_e_value` | The default argument `random_siRNA_results_overall` does not exist (`R/generate_e_value.R:12`). | The null matrix is a required argument. A helper makes it from `n_random` shuffles of the guide (default 1000), which keep the length and the base composition. The caller gives the random generator. The e-value of a gene is the fraction of shuffles with a site count that is not less than the observed count. Genes match by `gene_id`, not by `symbol`. |
| `load_annotations` | `hg38-old` has no ID (`R/load_annotations.R:42`, `:152`). `tx.names` is always `NULL` (`:116`). | There is no `hg38-old`. The result gives the transcript IDs. |
| `kmer_hist` | Returns `NULL + NULL` (`R/kmer_hist.R:32-36`). | Returns the counts. The plot module makes the histogram. |
| `match_distances` | A gene with one match gives `c(NA, 0)` (`R/match_distances.R:33-34`). | A gene with one match gives no distances. A gene with n matches gives each of the n(n − 1)/2 pair distances. |
| `.get_best_match` | Stops if the search did not make each expected column, for example with other edit distances (`R/full_search.R:96-125`). | The best match uses only the columns that the search made. The default output does not change. |
| `get_stringPPI_to_stringGGI` | Maps the STRING protein IDs through BioMart on the network. The human query uses the GRCh37 archive (`R/get_stringPPI_to_stringGGI.R:21-28`). | Maps the protein IDs with the `protein` table of the EnsDb of the same build. |
| Documentation | The vignettes call functions that do not exist: `filter_deseq`, `de_fc_ecdf`, `get_feature_seqs`. The network vignette omits `organism`. | The Python documentation uses only functions that exist. A test runs each documentation example. |
| `.get_counts_column` (used by `SeedMatchR`) | If the DE table has a column with the seed name, the result has two columns with this name (`R/SeedMatchR.R:214-217`). Then `res$mer7m8` reads the old column. | The counts replace the values of the column with the same name. A new column goes at the end. |
| `Biostrings::readDNAStringSet` (the sequences of SeedMatchR) | The name of a sequence is the complete header line after `>`. Then a header with a description does not match a `gene_id`. | The ID of a FASTA record ends at the first white space. |
| `Biostrings::DNAString`, `Biostrings::readDNAStringSet` | The gap letter `-`, the mask letter `+`, and `.` are kept in a sequence. | These letters raise `InputError` with the letter and its position. |
| `get_seed` with `allow_wobbles` | The wobble replaces only an uppercase `U` (`R/get_seed.R:106-108`, `:160-162`). A lowercase guide gets no wobble, and R gives no message. | The guide is changed to uppercase, thus a lowercase guide gets the wobble. |
| `get_seed` | The guide can hold each IUPAC letter (`R/get_seed.R:68`, `:111`, `:165`). | A guide holds only `A`, `C`, `G`, and `U`. The option `wobble` gives an ambiguous base. |
| `get_seed` with a custom window | A window with the start after the stop gives an empty target (`R/get_seed.R:157-168`). | This window raises `InputError`. |
| `download_parse_file` | The 8 names label the wrong columns of a table with more columns, or with an extra field in each row. A text value makes a number column text (`R/get_data.R:74-76`). R gives no message. | `read_deseq2_table` raises `InputError` for these tables. |
| `full_search_ranges` | Stops if every search finds a match, or if no search finds a match (`R/full_search.R:170-183`). | The result does not depend on the searches that find a match. |
| `full_search`, `full_search_ranges` | `full_search` stops for a duplicate seed or distance. `full_search_ranges` gives duplicate rows, and it takes any `max.seed.edit.dist` (`R/full_search.R:156-162`). | The options raise `InputError` for a duplicate and for a `max_seed_edit_distance` that is not an integer of 0 or more. |
| `SeedMatchReport` | A NaN cutoff gives no off-target. A bool cutoff, a text cutoff, or a text `log2FoldChange` or `padj` column gives a wrong count (`R/SeedMatchReport.R:94-99`). R gives no message. | The report raises `InputError` for these inputs. |
| `filter_res` | A NaN cutoff keeps no row, and a bool or text cutoff gives a wrong result (`R/filter_res.R:64-70`). R gives no message. A `NULL` cutoff stops with a dplyr error. | `filter_de` raises `InputError` for these cutoffs. |
| `SeedMatchReport` with `allow_wobbles` | The default `fixed = T` matches the wobble letter `R` only to `R`, thus the wobble finds fewer matches (`R/SeedMatchReport.R:51,55`). The R documentation tells the user to set `fixed` to `FALSE` with wobbles (`R/SeedMatchR.R:37`). | `ReportOptions` and `FullSearchOptions` raise `InputError` for `wobble=True` with `fixed_pattern=True`. With a fixed pattern, `annotate_seed_counts` raises `InputError` for a target with an ambiguity letter. |
| `.reduce_features` (used by `load_annotations` with `reduce.features`) | `GenomicRanges::reduce` gives the ranges of a gene in the order of the position (`R/load_annotations.R:210`). On the minus strand, `extractTranscriptSeqs` then joins the ranges in the wrong order. The rnor7 probe shows a wrong sequence. | The ranges of the union come in transcript order. |
| `.get_longest_utr` | A gene without a coding transcript stops the function with `subscript contains invalid names` (`R/load_annotations.R:408-415`), for example with `protein.coding = FALSE`. | Only a transcript with the feature is a candidate. A gene without a candidate has no sequence. |
| `load_annotations` with `longest.utr` or `reduce.features` | R applies the two options only to the 3' UTR. For the other feature types, it ignores them without a message (`R/load_annotations.R:82-91`). | `PerGene` applies to each feature type. |
| `ensembldb::transcripts` (used by `load_annotations` for `$gtf`, and by `.get_longest_utr` through `transcriptLengths`) | The SQL sorts on `seq_name` and `tx_seq_start` only. A tie of two genes then gets the scan order of the SQLite query plan, and the plan depends on the filter. The fixture golden gives two orders for one tie. | The order of a tie is the row order of the `gene` table, then of the `tx` table, for each filter. |
| `build_annotation_filter` (used by `load_annotations`) | A `transcript.support` of 0 or less, or an `entrez_id` of 0 or less, selects no transcript and gives no message (`R/load_annotations.R:298-330`). | `TranscriptFilter` raises `InputError` for a support level below 1 and for an Entrez ID below 1. |
| `stats::ks.test` (used by `ecdf_stat_test`) | The asymptotic p-value has numeric defects. The two-sided series keeps one term for `sqrt(en) * D` below 1 (`ks.c:100-111`). The one-sided value `1 - (-expm1(...))` gives 0 for a p-value below about 1e-16. D from a cumulative sum can be above 1, and then the p-value is 0. | The statistic comes from integer counts. The p-value is `scipy.special.kolmogorov(sqrt(en) * D)` for two-sided and `exp(-2 * en * D^2)` for one-sided. |
| `stats::wilcox.test` (used by `ecdf_stat_test`) | If each value is equal, the variance is 0 (`R/deseq_fc_ecdf.R:274-276`). The two-sided z is 0/0, and the p-value is NaN. The one-sided p-values are 1. | The p-value is 1 for each alternative. |
| `twosamples::dts_test` (used by `ecdf_stat_test`) | The function does not check its samples. An infinite value gives the statistic `Inf` and the p-value 1. An empty sample gives a statistic and the p-value 1 (`R/deseq_fc_ecdf.R:259-264`). | `dts_test` raises `InputError` for an infinite value and for an empty sample. |
| `stats::ks.test`, `stats::wilcox.test` (used by `ecdf_stat_test`), `stats::approx` (used by `wass_dist`) | Each test warns if ties make the p-value approximate (`R/deseq_fc_ecdf.R:255-257`, `:274-276`). `approx` warns `collapsing to unique 'x' values` for a sample with equal values (`R/wass_dist.R:70-71`). | No warning. The field `method` of a test result is `EXACT` or `ASYMPTOTIC`. The areas use the merged values, as in R. |
| `wass_dist` (the areas) | `n.interp` of 1 gives `NA` areas. A value of 0 or less stops with a message that does not name the cause (`R/wass_dist.R:69`). A pooled range above half of the largest float can overflow the sum. R then gives `Inf` or a meaningless area, and no message. | `wasserstein_areas` raises `InputError` for an `n_interp` below 2 and for a range above half of the largest float. |

## Tests

### Principles

- Each test asserts a behavior that a caller depends on. Each test can fail for a plausible defect.
- An error test asserts the exception type and a part of the message.
- A bug fix starts with a test that fails.
- The tests never run R. The core tests never use the network.

### Golden files from R

- R scripts in `tools/golden/` make the golden files in `tests/golden/`. The golden files are committed.
- The scripts use a project-local R library (`tools/golden/.rlib`, ignored by git). They install SeedMatchR from `../SeedMatchR` at commit `2fc0c4f`, with its Bioconductor dependencies.
- The scripts load each EnsDb and 2bit file from a local path. They replace `load_species_anno_db` in the SeedMatchR namespace. Thus `load_annotations` runs unchanged, without AnnotationHub.
- Each golden file records the R version and the package versions.
- Tables are TSV files. Scalar results are JSON files.
- For a defect fix, a patch in `tools/golden/patches/` applies the fix to a copy of the R function, if the fix is unambiguous. The golden of the patched function then does a test of the rest of the ported logic.
- The test of a divergence names its ID from `docs/divergences.md`. A test asserts that each ID in the record has at least one test that names it.

### Search

- Compatibility table: R gives the match result for each pair of IUPAC letters, for each `fixed` pair. A test compares the complete table.
- Random corpus: a seeded generator makes patterns and subjects. They contain IUPAC codes, runs of `N`, repeats, and matches near the ends. R gives the match ranges for each case and each option set: mismatches 0 to 4, indels on and off, each `fixed` pair, and `min_mismatch`. A test compares the ranges exactly.
- Oracle: a slow Python version of each rule, written from the definitions. Hypothesis property tests compare the numba kernels with the oracle on generated input. The R corpus also does a test of the oracle.
- Invariants: the count equals the number of ranges. An independent edit-distance function gives each reported range a distance that is not more than the limit.

### Statistics

- R goldens cover a grid: sample sizes on the two sides of each branch threshold, data with and without ties, and each alternative.
- DTS: the statistic equals the R value. With a fixed generator, the p-value equals the fraction of the permutation statistics that are not less than the observed statistic, with the twosamples formula. The p-value agrees with R within the binomial error of the permutations.

### Annotation

- Synthetic fixture: a Python script in `tools/fixtures/` writes the fixture from one specification: an EnsDb with schema 2.2, a 2bit file, a FASTA file, an Ensembl GTF, and a GENCODE GTF. `bgzip` of htslib and `samtools faidx` write the BGZF file and its indexes. Thus the FASTA reader reads files from an independent writer. R reads the EnsDb, the 2bit file, and the BGZF file to make the goldens. Thus R also makes sure that the 2bit writer of the script is correct.
- `ensembldb::ensDbFromGtf` does not make the fixture, because it writes schema 1.0 and it downloads the contig lengths.
- The fixture has these cases:
  - a minus-strand transcript with a 3' UTR across more than one exon
  - a CDS that ends at an exon boundary
  - a coding transcript without a 3' UTR
  - a non-coding transcript
  - a transcript on a non-standard contig
  - a gene with more than one transcript, with two transcripts that have the same 3' UTR length and the same start
  - 3' UTRs of one gene that overlap or touch
  - soft-masked bases and N blocks
  - a range across a BGZF block boundary
- Real references: tests with the marker `refs` read the directory in `SEEDMATCH_TEST_REFS`. They compare the feature sequences of `rnor7` with R digests: the names, the length, and the first 16 hexadecimal digits of the SHA-256 digest of each sequence. The digest files are gzip TSV files. A later plan adds `hg38`. If the variable is not set, the tests skip and give the reason.
- Consistency: on the fixture, each GTF gives the same sequences as the EnsDb. On `rnor7`, the Ensembl release 109 GTF gives the same sequences as the EnsDb of release 109.

### Core without network

- `pytest-socket` disables sockets for the test session.
- A test imports each core module and runs the core pipeline on the fixtures. It makes sure that no network module is in `sys.modules` after the run.

### Fetch

- `seedmatch.fetch` takes a `FetchOpts` value with the HTTP client (`DEFAULT_FETCH_OPTS` uses `urllib.request`). The tests replace the HTTP client with one that serves fixture bytes from memory.
- Cases:
  - a digest mismatch leaves no final file
  - an interrupted transfer leaves no final file
  - a present file with the correct digest causes no download
  - a URL or a redirect that is not `https` causes an error

### Documentation examples

A test runs each code example in `README.md` and in `docs/`.

### Coverage

`pytest-cov` measures branch coverage. The floor is 95% for `src/seedmatch`, except `seedmatch.plot`. Each line without coverage needs a reason in the review.

## Performance and bounds

- `SequenceSet` keeps all bases in one buffer. Memory is linear in the total length.
- The 2bit reader decodes by range, and the EnsDb reader filters in SQL. Thus the request, not the file, sets the size of a read.
- A mismatch search is O(N · m) at worst, with an early stop when the mismatch count passes the limit. An indel search is O(N · m) dynamic programming.
- The kernels use `nogil=True`. `full_search` runs its patterns on a thread pool of fixed width. The width is a parameter, with the CPU count as the default.
- `docs/performance.md` records the wall time and the peak memory of `seed_match_report` on the canonical 3' UTRs of `rnor7`, for Python and for R.
- Target: after the JIT compilation, Python is not slower than R by more than a factor of 2.

## Tooling

- uv: `uv sync` and `uv run`. The `dev` group holds pytest, pytest-cov, hypothesis, pytest-socket, mypy, pandas-stubs, and ruff.
- ruff: `select = ["ALL"]`. Each ignore has a reason in `pyproject.toml`, for example `COM812` and `ISC001`, which conflict with the formatter. `line-length = 100`, as in the sibling project `tcga`. `target-version = "py312"`.
- mypy with `--strict`. An override applies only to the modules with numba decorators.
- `README.md` gives the commands: `uv run ruff check`, `uv run ruff format --check`, `uv run mypy`, `uv run pytest`.

## Delivery order

1. Scaffold: the uv project, the tool configuration, the empty package, and the golden R environment. `git init` waits for the consent of the user.
2. `sequence` and `seed`.
3. `search`: kernels, oracle, corpus.
4. `report`: the SeedMatchR modes, full search, best match, report.
5. `de`, `stats`, and the ECDF test.
6. `annotation`: 2bit, EnsDb, filters, features, extraction, `Reference`, the GTF path, and the real-reference parity.
7. `resources`, `fetch`, and the command.
8. The other analyses: k-mers, alignment, distances, non-canonical searches, network, benchmark, e-value.
9. `plot` and the report HTML.
10. Documentation: README, the map from R names, divergences, performance.

Each step ends with ruff, mypy, and the tests of that step.

## Decisions

### Decided by the user

- The name is `seedmatch`, in the sibling directory `../seedmatch`.
- uv, ruff, and Python 3.12.
- The goal is a better version of SeedMatchR, not only a copy. Improvements are welcome.
- Output identical to R on the core paths, except for the recorded divergences.
- Each SeedMatchR function gets a port. Each defect gets a fix, also on the core paths.
- `docs/divergences.md` records each change of behavior, and it is part of the repository.
- The core reads resources from given paths and never downloads. A separate module downloads.
- The users are agents in the inflexa sandbox.
- The license is MIT, the same license as SeedMatchR.
- The semantics of the wobble and bulge search, of the message passing, and of the e-value null model in "Divergences". The designer proposed them, and the user did not change them.

### Calls of the designer

- Numba kernels. Runtime dependencies only from the inflexa store lock, with the `plot` extra. The genome readers use only the standard library and numpy.
- R coordinates: 1-based and inclusive.
- No layer with the R function names.
- `SequenceSet` keeps the transcript ID and the gene ID. Gene-level results use the maximum count and the best match.
- The EnsDb `protein` table maps the STRING protein IDs.
- ruff `ALL`, mypy `--strict`, and the 95% coverage floor.
- TSV and JSON golden files.
- The thread pool for the patterns, and the performance target.

### Left to the builder

- The internal names below the public API.
- The kernel algorithms, within the search rules.
- The size and the seeds of the random corpus.
- The layout of the golden files.
- An entry in `docs/divergences.md` for each new divergence that the work finds.