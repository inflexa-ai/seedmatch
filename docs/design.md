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
- The genome readers use only the standard library (`mmap`, `zlib`, `sqlite3`) and numpy. Thus the port reads a BGZF FASTA genome without an extra.
- A new runtime dependency must have an entry in the inflexa manifest. Thus, do not add one.
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
| `wass_dist` (the areas and the data of `diff.plot`) | `wasserstein_areas`, `wasserstein_difference` |
| `load_annotations`, `get_3UTR_granges`, `.get_longest_utr`, `.reduce_features`, `build_annotation_filter` | `Reference.features`, `Reference.sequences`, `TranscriptFilter` |
| `sequence_kmer_counts` | `kmer_counts` |
| `plotAlignment` with `output.type = "data.frame"` | `alignment_text` |
| `download_parse_file`, `load_example_data` | `read_deseq2_table` |

### Not identical by design

- The DTS p-value comes from permutations. The statistic is identical. The p-value agrees with R only within random error. The caller gives the random generator.
- Plots use matplotlib.
- The report table is the HTML of `report_html`, not `gt`. pandas 2.3.3 gives a `Styler` only with jinja2 3.1.2 or later, and the store lock has jinja2 3.0.3.
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
| `seedmatch.mirdb` | The reader of the miRDB table of predicted miRNA targets. |
| `seedmatch.stats` | KS and Wilcoxon with the R rules, DTS, Wasserstein areas, ECDF test, PR and ROC curves, and the benchmark of seed counts. |
| `seedmatch.report` | Full search, best match, report, and the HTML of the report table. |
| `seedmatch.noncanonical` | Wobble and bulge search, 3' compensatory search. |
| `seedmatch.network` | TRRUST message passing and the STRING gene network. |
| `seedmatch.kmers`, `seedmatch.alignment`, `seedmatch.distances`, `seedmatch.evalue` | The smaller analyses. |
| `seedmatch.plot` | Figures from results: the ECDF, the Wasserstein areas, the k-mer histogram, the alignment, the ranges, and the PR and ROC curves. Extra `plot`. |
| `seedmatch.resources` | The registry of the upstream resources: identifiers, URLs, file names, digests. No I/O. |
| `seedmatch.fetch` | Downloads, and the command `seedmatch-fetch`. |

### Rules

- Only `seedmatch.fetch` can import a network module. The ruff rule `TID251` (banned API) bans these modules in all other code: `urllib.request`, `http.client`, `socket`, `ssl`, `urllib3`, `requests`, `httpx`. A per-file ignore lets `src/seedmatch/fetch/**` and `tests/fetch/**` use them.
- The same rule bans `seedmatch.fetch` in each module outside `seedmatch.fetch`.
- `import seedmatch` does not import matplotlib or `seedmatch.fetch`.
- A computation function returns data. A plot function takes a result and returns a `matplotlib.figure.Figure`. No function writes a file unless it has an explicit output path parameter.
- A plot function uses the object API of matplotlib. It does not import `matplotlib.pyplot`, and it does not change `matplotlib.rcParams`.
- A reader reads only the part of a file that the call uses:
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
- On the test data (`ecdf.tsv`), the port gives D = 0.51035 and p = 0.0001037731, as R does.

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
- On the test data (`ecdf.tsv`), the port gives W = 748 and p = 0.9996 for `GREATER`, as R does.

### DTS test

- The statistic is a sum over the gaps between the sorted pooled values.
- Each term is |E − F| / sqrt(2 · G · (1 − G) / n), multiplied by the width of the gap.
- E and F are the ECDFs of the two samples. G is the pooled ECDF. n is the pooled size. Each ECDF is at the left value of the gap.
- A gap of width 0 adds nothing to the sum. Thus only the end of a run of equal values counts.
- On the test data (`ecdf.tsv`), the statistic is 74.91157, as in R.
- The port computes E − F from integer counts: (c · n − k · n_x) / (n_x · n_y). Here k values lie at or below the gap, and c of them come from the first sample.
- A permutation chooses n_x of the pooled sorted positions for the first sample. This is the method of twosamples.
- The p-value is the fraction of `n_boot` permutation statistics that are not less than the observed statistic (default 2000, as in twosamples). If the fraction is 0, the p-value is 1 / (2 · `n_boot`).
- The observed statistic and each permutation statistic come from the same code, and each row adds its terms in the same order. Thus a permutation with the same counts as the observed arrangement gives the same value, and it counts.
- The caller gives a `numpy.random.Generator` or a non-negative int seed.
- The port draws the permutations in blocks, and numpy does the work of each block. A block has `max(1, min(64, 2^20 // n))` rows, thus at most `max(2^20, n)` elements. Thus the memory does not grow with `n_boot`. The peak is 27 MB for 20000 values, 35 MB for 200000 values, and 80 MB for 10^6 values. The time is O(`n_boot` · n).
- An infinite value and an empty sample raise `InputError`. R gives the statistic `Inf`, or a statistic that does not compare two samples, and the p-value 1 (`D-034`).
- A pooled range above the largest float divided by 2 n raises `InputError` (`D-034`). Each term is at most the width of its gap times n, thus the statistic is at most the range times n. The port multiplies the width by a factor of at most n / (n − 1). It does not multiply by n^1.5 first. Thus no product overflows below the limit. R overflows to the statistic `Inf` and the p-value 1.

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

### ECDF test

- `ecdf_test` uses `GREATER` for KS and for WILCOXON if `alternative` is `None`, as R does. The direction of `GREATER` is not the same for KS and for WILCOXON.
- For KS, `GREATER` means that the target ECDF is above the background ECDF. Thus the target values are lower.
- For WILCOXON, `GREATER` means that the target values are higher.
- On the test data (`ecdf.tsv`), the target genes have lower values. KS `greater` gives p = 1.04e-4, and WILCOXON `greater` gives p = 0.9996. WILCOXON `less` gives p = 3.9e-4.
- The port keeps the default of R. The docstring of `ecdf_test` and the README give the direction.

### PR and ROC curves

`pr_curve` and `roc_curve` replace the calls `pr.curve(scores.class0 = scores, weights.class0 = labels, curve = TRUE)` and `roc.curve(...)` of PRROC 1.4. PRROC has the GPL-3 license. Thus the port writes the curves from Keilwagen, Grosse, and Grau (2014, PLOS ONE 9(3)) and from Davis and Goadrich (2006, ICML). The R code of PRROC settles the details that the papers leave open. The port does not translate it.

Inputs:

- `scores` goes through `as_sample` with `finite=False`. NaN raises `InputError`, and an infinite score is valid.
- `labels` is 1-D and has the length of `scores`. Each label is a bool, or a number equal to 0 or 1. A label of 1 or `True` is a positive. PRROC calls the positives class 0.
- Labels without a positive raise `InputError`. `roc_curve` also raises for labels without a negative. `pr_curve` accepts labels without a negative, and then the precision is 1 at each point.

Definitions:

- Let s_1 > s_2 > ... > s_K be the distinct scores. TP_k and FP_k count the positives and the negatives with a score of at least s_k. TP_0 = FP_0 = 0.
- The scores -0.0 and 0.0 are one score. The threshold of a run of equal scores is its last score in the order of the input, as in PRROC.
- Let P = TP_K and N = FP_K. The recall is r_k = TP_k / P, and the precision is p_k = TP_k / (TP_k + FP_k).
- If TP_k > TP_{k−1}, then h_k = (FP_k − FP_{k−1}) / (TP_k − TP_{k−1}). Between the points k − 1 and k, FP grows on a line: FP = FP_{k−1} + h_k · (TP − TP_{k−1}).

`auc_integral` (Keilwagen, Grosse, and Grau 2014):

- On step k, the precision at the recall r is r / (a_k · r + b_k), with a_k = 1 + h_k and b_k = (FP_{k−1} − h_k · TP_{k−1}) / P.
- The area of step k is the integral of this precision from r_{k−1} to r_k. Let c_{k−1} = (TP_{k−1} + FP_{k−1}) / P, which is a_k · r_{k−1} + b_k. Then v_k = ((r_k − r_{k−1}) − (b_k / a_k) · log1p(a_k · (r_k − r_{k−1}) / c_{k−1})) / a_k.
- If b_k = 0, then v_k = (r_k − r_{k−1}) / a_k. If TP_k = TP_{k−1}, then v_k = 0.
- `auc_integral` is the sum of v_k.
- PRROC writes the logarithm as ln(a_k · r_k + b_k) − ln(a_k · r_{k−1} + b_k). This difference loses digits if the negatives far outnumber the positives (`D-050`). The port uses `log1p`, which gives the same value in exact arithmetic. It loses fewer digits, but a smaller error remains, because it still subtracts a term from the width of the step. For 11 positives in 200000 scores, the relative error is 2.4e-13 for the port and 4.2e-8 for PRROC. For 1 positive in 1000000 scores, it is 3.2e-10 for the port and 3.7e-4 for PRROC.

`auc_davis_goadrich` (Davis and Goadrich 2006):

- Step 1 adds (TP_1 / P) · p_1. PRROC uses the precision of point 1 for the full width of step 1.
- A step k ≥ 2 with TP_k − TP_{k−1} ≥ 2 and p_k ≠ p_{k−1} gets a point at each integer TP. Let Δ = TP_k − TP_{k−1} and q(x) = (TP_{k−1} + x) / (TP_{k−1} + x + FP_{k−1} + h_k · x). Such a step adds (1 / P) · (q(x) + q(x + 1)) / 2 for x = 0 to Δ − 1.
- Each other step adds the trapezoid ((TP_k − TP_{k−1}) / P) · (p_{k−1} + p_k) / 2.
- `auc_davis_goadrich` is the sum. The labels are hard, thus this area is always defined. PRROC gives `NA` for soft labels, which the port does not take.

The PR curve (`curve`, with the columns `recall`, `precision`, and `threshold`):

- The base rows are (r_k, p_k, s_k) for k = 1 to K.
- The recall step is σ = min(1, P / 100) / P. It comes from the PRROC default `minStepSize = min(1, P / 100)`.
- Each step k ≥ 2 with (TP_k − TP_{k−1}) / P > σ and p_k ≠ p_{k−1} gets intermediate rows. Let d = r_k − r_{k−1} and m = `int(d / σ + 1e-10)`. For j = 0 to m, let x = min(j · σ, d). The row is (r_{k−1} + x, (r_{k−1} + x) / (r_{k−1} + x + FP_{k−1} / P + h_k · x), s_k). This is `seq(0, d, by = σ)` of R 4.5.3.
- The intermediate rows come after the base rows, in the order of k and then of j.
- The sort puts the rows in the order of `recall`, from the highest to the lowest value. Rows with equal `recall` come in the order of `threshold`, from the lowest to the highest value. The sort is stable, as `order(..., decreasing = TRUE)` of R is.
- A first row with the recall 1 and a last row with the recall 0 follow. Each takes the precision and the threshold of the sorted row next to it.
- The code computes each value with the IEEE operations of these formulas, in their order. Thus the recalls that R gives as equal stay equal, and the stable sort gives the R order.
- The curve has at most 2K + 2 + max(100, P) rows.

The ROC curve:

- The points are (0, 0) and (FP_k / N, TP_k / P) for k = 1 to K. `auc` is the sum over k of (fpr_k − fpr_{k−1}) · (tpr_k + tpr_{k−1}) / 2.
- `curve` has the columns `fpr`, `tpr`, and `threshold`, with one row for each point. The threshold column keeps the rule of PRROC. The point (0, 0) has s_1, the point k has s_{k+1} for k < K, and the point K has s_K. Thus the threshold of a row is the score of the points that the next row adds.

## Annotation

### EnsDb and 2bit (the parity path)

- The EnsDb is an SQLite file with the ensembldb schema. The reader uses the standard `sqlite3` module. It opens the file read-only and immutable (`file:...?mode=ro&immutable=1`, URI open). Thus the reader can open a file in a read-only store.
- The file must have the tables `gene`, `tx`, `tx2exon`, and `exon`. A filter that uses an optional column (`tx_is_canonical`, `tx_support_level`) or the `entrezgene` table raises `InputError` if the file does not have it. The message gives the file and the column. `protein_genes` uses the `protein` table, and a file without it raises `InputError`.
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
- The registry pins the URL, the file name, the size, and the SHA-256 digest. The values come from a streamed download of each file.
- The other resources come from the URLs in the R code:
  - the DESeq2 tables of GEO GSE184929 (`R/get_data.R:23`, `:27`, and `:31`)
  - the miRDB v6.0 predictions (`R/get_data.R:39`)
  - the TRRUST human and mouse tables (`R/get_network_matches.R:21-22`)
  - the STRING protein links for a species (`STRINGdb::STRINGdb$new` in `R/get_stringPPI_to_stringGGI.R:15`), with the pinned STRING version 12.5 (`D-041`)

### Fetch behavior

- The fetch module uses `urllib.request` from the standard library. It adds no dependency.
- The fetch module uses `fcntl.flock` for its lock, thus it runs only on POSIX systems. The CI runs on Linux.
- It streams each file in chunks to `<dest>/<name>.part`. The option `chunk_size` (default 1 MiB) limits the memory. The read stops at the registry size plus 1 byte. Thus a server that sends more bytes, or that sends no `Content-Length` header, cannot make the read longer. Then it makes sure that the size and the SHA-256 digest agree with the registry. Then it renames the file to its final name.
- If the response has a `Content-Length` header that is not the registry size, it raises `FetchError` before it reads the body. The message names the URL, the header value, and the registry size. A header that is not a decimal number raises `FetchError` too.
- It locks the `.part` file with `fcntl.flock` before the request. The run keeps the lock until the file gets its final name or the run deletes the file.
- If another run holds the lock, this run raises `FetchError` that names the `.part` path. This run sends no request. It does not delete the file of the other run.
- After it gets the lock, it does a second check of the final path. The rules for a present file apply to a file that the other run made.
- If the run cannot make, open, or lock the `.part` file, it raises `FetchError` before the request. The message names the `.part` path. A file system without `flock` support thus stops the download.
- If the final file is present with the size and the digest of the registry, it sends no request.
- If the final file is present with another size or digest, it raises `FetchError`. The message names the path, the two sizes, and the two digests. The file stays, and the message tells the user to remove the file before the next download.
- It refuses a URL that is not `https`.
- It goes to the target of a redirect itself, a maximum of `max_redirects` times (default 5). The default client (`open_url`) does not go to the target of a redirect. The module makes sure that the target is `https` before it sends a request to the target. A redirect to `http` raises `FetchError`, and the `http` URL gets no request.
- A failed download leaves no final file. An interrupted transfer raises `FetchError` that names the URL and the count of received bytes.
- An error or Ctrl-C removes the `.part` file while the run holds the lock. If the run cannot delete the `.part` file, the file stays, and the next run replaces it. If a signal kills the process, the next run also replaces the `.part` file.
- A `.part` path that is not a regular file raises `FetchError` before the request.
- The command `seedmatch-fetch` has a subcommand for each resource group. Examples:
  - `seedmatch-fetch reference rnor7 DEST`
  - `seedmatch-fetch example sirna DEST`
  - `seedmatch-fetch trrust human DEST`
  - `seedmatch-fetch string 10116 DEST`
  - `seedmatch-fetch mirdb DEST`
  - `seedmatch-fetch list`
- The command prints the path of each file with the word `present` or `downloaded`. An `InputError` or a `FetchError` gives exit status 1 and the message on standard error. Ctrl-C prints the message `seedmatch-fetch: interrupted`. Then the process ends by SIGINT, thus a shell gives the status 130 and stops a loop. The command flushes standard output after each path line.

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
| `get_compensatory_sequences` | Always returns 0 rows: the empty-result return runs before the search (`R/get_compensatory_sequences.R:57-77`). Without it, the function stops at the first site (`:99`). | The search runs, and the function returns the found regions. |
| `run_wobble_bulges_pipeline` | Stops with `object 'bulge_allowance' not found` (`R/run_wobble_bulges_pipeline.R:167`, `:192`). It swaps two arguments (`:27`, `:160`), uses only sets of exactly 4 U positions (`:82-93`), and sums the counts of variants that overlap (`:35-55`). | A guide U can pair with a target G: one IUPAC pattern covers this at each U position. A bulge is 1 to `bulges_lim` extra target bases between two adjacent seed positions. No bulge occurs next to the `nt` seed positions from g2 (default g2 to g4), as the R documentation states. The result counts each distinct start position once for each transcript. |
| `benchmark_preds` | `padj <= padj` ignores the threshold (`R/benchmark.R:23`). `roc.curve` gets the labels as `scores.class1` (`:36-38`). | The threshold applies. The ROC curve uses the labels as weights, as the PR curve does. |
| `message_passing` | `igraph::neighbors` has the default mode `out`, thus a gene gets the messages of its targets (`R/get_network_matches.R:59`). `length(neighbors != 0)` is the number of neighbors (`:61`). `mean(a, b)` uses `b` as `trim`, and `features_act` reads `features_rep` after its change (`:84-85`, `:110`). | TRRUST edges go from a regulator to a target. The score of a gene comes from its regulators. The type of each edge selects the repression channel or the activation channel. An `Unknown` edge is in no channel. The score is (mean of the seed flags of the regulators + own seed flag) / 2. `degree_importance` divides each regulator flag by the out-degree of the regulator. |
| `generate_e_value` | The default argument `random_siRNA_results_overall` does not exist (`R/generate_e_value.R:12`). | The null matrix is a required argument. A helper makes it from `n_random` shuffles of the guide (default 1000), which keep the length and the base composition. The caller gives the random generator. The e-value of a gene is the fraction of shuffles with a site count that is not less than the observed count. Genes match by `gene_id`, not by `symbol`. |
| `load_annotations` | `hg38-old` has no ID (`R/load_annotations.R:42`, `:152`). `tx.names` is always `NULL` (`:116`). | There is no `hg38-old`. The result gives the transcript IDs. |
| `kmer_hist` | Returns `NULL + NULL` (`R/kmer_hist.R:32-36`). | Returns the counts. `plot_kmer_histogram` draws the histogram. |
| `match_distances` | A sequence with one match gives `c(NA, 0)` (`R/match_distances.R:33-34`). | A sequence with one match gives no distances. A sequence with n matches gives each of the n(n − 1)/2 pair distances. |
| `.get_best_match` | Stops if the search did not make each expected column, for example with other edit distances (`R/full_search.R:96-125`). | The best match uses only the columns that the search made. The default output does not change. |
| `get_stringPPI_to_stringGGI` | Maps the STRING protein IDs through BioMart on the network. The human query uses the GRCh37 archive (`R/get_stringPPI_to_stringGGI.R:32-60`). | Maps the protein IDs with the `protein` table of the EnsDb of the same build. |
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
| `stats::ks.test` (used by `ecdf_stat_test`) | The asymptotic p-value has numeric defects. The two-sided series keeps one term for `sqrt(en) * D` below 1 (`ks.c:100-111`). The one-sided value `1 - (-expm1(...))` has an absolute error of about 1e-16. Its relative error is above 1e-9 for a p-value below about 1e-7, and the value is 0 below about 1e-16. D from a cumulative sum can be above 1, and then the p-value is 0. | The statistic comes from integer counts. The p-value is `scipy.special.kolmogorov(sqrt(en) * D)` for two-sided and `exp(-2 * en * D^2)` for one-sided. |
| `stats::wilcox.test` (used by `ecdf_stat_test`) | If each value is equal, the variance is 0 (`R/deseq_fc_ecdf.R:274-276`). The two-sided z is 0/0, and the p-value is NaN. The one-sided p-values are 1. | The p-value is 1 for each alternative. |
| `twosamples::dts_test` (used by `ecdf_stat_test`) | The function does not check its samples. An infinite value gives the statistic `Inf` and the p-value 1. An empty sample gives a statistic and the p-value 1. Finite values near ±1e308 overflow to the statistic `Inf` and the p-value 1 (`R/deseq_fc_ecdf.R:259-264`). | `dts_test` raises `InputError` for an infinite value, for an empty sample, and for a pooled range above the largest float divided by 2 n. |
| `ecdf_stat_test` with `twosamples::dts_test` or `wass_dist` | R passes `alternative` only to `ks.test` and `wilcox.test`. For `DTS` and `Wass`, R ignores it and gives no message (`R/deseq_fc_ecdf.R:255-277`). | `ecdf_test` raises `InputError` for an `alternative` with `DTS` or `WASSERSTEIN`. |
| `ecdf_stat_test` (the rows without a fold change) | `ks.test` and `wilcox.test` drop the NA values. `dts_test` stops with `missing value where TRUE/FALSE needed`, and `wass_dist` stops with `'from' must be a finite number` (`R/deseq_fc_ecdf.R:247-277`). | `ecdf_test` drops the rows with a NaN `log2FoldChange` for each test. The direct functions raise `InputError` for NaN. |
| `stats::ks.test`, `stats::wilcox.test` (used by `ecdf_stat_test`), `stats::approx` (used by `wass_dist`) | Each test warns if ties make the p-value approximate (`R/deseq_fc_ecdf.R:255-257`, `:274-276`). `approx` warns `collapsing to unique 'x' values` for a sample with equal values (`R/wass_dist.R:70-71`). | No warning. The field `method` of a test result is `EXACT` or `ASYMPTOTIC`. The areas use the merged values, as in R. |
| `wass_dist` (the areas) | `n.interp` of 1 gives `NA` areas. A value of 0 or less stops with a message that does not name the cause (`R/wass_dist.R:69`). A pooled range above half of the largest float can overflow the sum. R then gives `Inf` or a meaningless area, and no message. | `wasserstein_areas` raises `InputError` for an `n_interp` below 2 and for a range above half of the largest float. |
| `.check_gene_list_overlap` (used by `deseq_fc_ecdf`) | R warns `Gene lists overlap!` for `any(duplicated(unlist(gene_lists)))`, thus also for a gene that occurs two times in one list (`R/deseq_fc_ecdf.R:165-169`). | `assign_groups` warns with `SeedMatchWarning` only if two groups share a gene. The message names the two groups and one gene. |
| `.assign_groups` (used by `deseq_fc_ecdf`) | The text `"NA"` marks a gene in no group, and the `factor` call makes it a level if a group has the name `NA`. Then this level also holds each gene in no group (`R/deseq_fc_ecdf.R:109`, `:171-182`). | `assign_groups` gives no value to a gene in no group. A group named `NA` holds only its own genes. |
| `get_stringPPI_to_stringGGI` (the STRING download) | `STRINGdb$new` gets no `version`, thus STRINGdb 2.22.0 asks the STRING server for the current version at each run (`R/get_stringPPI_to_stringGGI.R:15`, `R/rstring.R:145-151` of STRINGdb 2.22.0). The network depends on the date of the run. | The registry pins STRING 12.5, with the size and the SHA-256 digest of each file. |
| `kmer_hist` (the gene count) | `colSums` of the presence matrix counts the sequences with a k-mer, but the label says genes (`R/kmer_hist.R:24-28`, `:32`). | `kmer_gene_counts` counts the distinct gene IDs. |
| `sequence_kmer_counts` (the width) | `oligonucleotideFrequency` accepts a width of 13 and more. It makes a dense matrix with 4^13 columns for each sequence (`R/kmer_hist.R:63-66`). | `kmer_counts` raises `InputError` for a `k` above 12, and its matrix is sparse. |
| `sequence_kmer_counts` (the names) | A repeated sequence name stops the function with `duplicate 'row.names' are not allowed` (`R/kmer_hist.R:71`). | `kmer_counts` gives one row for each sequence, in the order of the set. |
| `plotAlignment` | The guide target comes from `get_seed(..., "Full", ...)`, which starts at g2 (`R/plotAlignment.R:32`, `R/get_seed.R:121`). The window has one base more, thus each call with the default positions stops. A row of `SeedMatchR` with `res.format = "data.frame"` has no column `names`. | `alignment_text` shows g1 to gN, and it finds the sequence with `seq_id`. |
| `plotAlignment` (a window past a sequence end) | `XVector::subseq` stops for a window that starts before 1 or ends after the sequence (`R/plotAlignment.R:43-45`). | `alignment_text` shows a space for each position outside the sequence. |
| `generate_e_value` (the counts) | A count of `NA`, a negative count, `Inf`, and `-Inf` in `result_output` give the e-value 0 (`R/generate_e_value.R:16-21`). In the null, R does not count a cell with `NA`, a negative number, or `-Inf`, but the cell stays in the denominator. A gene with two rows in the null gives an e-value above 1. Two rows of one gene in `result_output` both get the value of a recycled comparison. | `e_values` raises `InputError` for a count that is NaN, infinite, or negative, and for a gene with two rows in the null. Each row gets its own e-value. |
| `benchmark_preds` (the inputs) | A NaN score stops with `NAs are not allowed in subscripted assignments`. A NaN cutoff keeps no gene, a text cutoff compares the values as text, and a bool cutoff is the number 1. A gold standard without a positive gives NaN areas. One without a negative gives the PR areas 1 and 1, and the ROC area NaN (`R/benchmark.R:22-38`). | `benchmark_predictions` raises `InputError` for a NaN score and for a cutoff that is a bool, text, or NaN. It also raises `InputError` for a gold standard without a positive row or without a negative row. |
| `get_example_data` with `example.type = "mirna"` (the miRDB table) | `read.table` with the defaults gives a text column for a score that is not a number. It gives `NA` for an empty score and for the ID text `NA`. A `#` starts a comment, and R stops. A quote character in a field loses rows, and R gives no message (`R/get_data.R:43-46`). | `read_mirdb` raises `InputError` for a score that is not a finite decimal number. The ID fields keep `#`, quote characters, and `NA` as text, and each line is a row. It also raises `InputError` for a line with a NUL character. It does the same for a line of more than 1024 characters, which bounds the memory of a chunk. R reads a long line, and it cuts a score at a NUL. |
| `pr.curve` (the integral area) | `compute.pr` of PRROC 1.4 gives the area of a step as the difference of two logarithms (lines 71 and 72 of the printed function). The difference loses digits if the negatives far outnumber the positives. The relative error is 4.2e-8 for 11 positives in 200000 scores. | `pr_curve` computes the same area with `log1p`. The relative error is 2.4e-13 on the same sample, and 3.2e-10 for 1 positive in 1000000 scores (PRROC: 3.7e-4). |
| `ecdf_stat_test`, `.assign_groups`, `.get_counts_column` with `shared_genes = FALSE` (a repeated column) | `res$log2FoldChange` and `res$gene_id` take the first column of that name, and R gives no message (`R/deseq_fc_ecdf.R:176`, `:248`, `:252`, `R/SeedMatchR.R:207`). The result depends on the column order. `.get_counts_column` appends the counts after a column with the name of the count (`R/SeedMatchR.R:214-217`). | `ecdf_test`, `assign_groups`, and `annotate_seed_counts` raise `InputError` for a repeated column that they use. `annotate_seed_counts` also raises it for a repeated column with the name of the count column. |
| `get_compensatory_sequences` (the sequences) | Selects the sequences with a seed match by name, thus a repeated name gives the first sequence two times (`R/get_compensatory_sequences.R:36-41`). The merge by name pairs each seed row with each region row of the name (`:90`). | Each sequence of the set is searched, and a pair holds a seed match and a region of the same sequence. |
| `get_compensatory_sequences` (the filter) | Two `dplyr` filters keep the longest region of each transcript and region end, and then of each transcript and region start. The groups ignore the seed match (`R/get_compensatory_sequences.R:123-125`). | A pair is dropped only for a longer region of the same seed match that contains its region. |
| `get_compensatory_sequences` (the columns) | `relative_position_seq` gives each start of the region in the guide as a text with commas (`R/get_compensatory_sequences.R:3-8`, `:127-128`). `mismatch_in_seed` is 1 if the seed text differs from the target (`:106-108`). | `guide_start` and `guide_end` are the integers of the window of the row. `seed_mismatches` is the number of mismatches. |
| `get_compensatory_sequences` (the lengths) | A minimum length above the maximum searches the lengths from the minimum down to the maximum. A minimum above the guide positions searches the same region again. A `limit` of 2 or less gives no row (`R/get_compensatory_sequences.R:46-54`, `:92`). | `CompensatoryOptions` and `compensatory_sites` raise `InputError` for these values, and for a longest region above the pattern limit of the search. |
| `get_compensatory_sequences` (the region start) | The region is always g8 to the end of the guide (`R/get_compensatory_sequences.R:44`). | The default region starts after the seed window. `first_position=8` gives the R region. |
| `get_network_matches` (the TRRUST file) | `read.table(..., header = TRUE)` reads line 1 as the column names, but the TRRUST file has no header. R loses the first edge (`R/get_network_matches.R:25`, `:27`). A `'` or a `#` stops R, and the text `NA` becomes the value `NA`. | `read_trrust` keeps each line as an edge. The fields keep quote characters, `#`, and `NA` as text. A bad mode, an empty gene, another number of fields, a NUL character, and a line of more than 1024 characters raise `InputError`. |
| `get_network_matches` (the seed genes) | If no gene of `ENSG_Genelist` is a vertex, each value is 0, and R gives no message (`R/get_network_matches.R:36-37`). The vignette gives Ensembl gene IDs, but TRRUST names each gene by its symbol (`vignettes/articles/SeedMatchR_build_PPI_network.Rmd:28-30`). | `trrust_scores` warns with `SeedMatchWarning` if no seed gene is in the network. |
| `get_stringPPI_to_stringGGI` (the STRING links file) | STRINGdb reads the file with `read.table(..., sep = " ", header = TRUE, fill = TRUE)`. It keeps the rows with `combined_score >= score_threshold` (`R/get_stringPPI_to_stringGGI.R:15-16`, `R/rstring.R:1064-1067` of STRINGdb 2.22.0). A text score makes the comparison use text. A `#`, a `'`, and a line with another number of fields give rows of `NA` values or move the columns. R gives no message. | In `read_string_links`, line 1 is the header and each other line has 3 fields. It raises `InputError` for a score that is not an integer from 0 to 1000 and for an empty protein ID. It does the same for a NUL character and for a line of more than 1024 characters. `#` and quote characters are text. |
| `get_stringPPI_to_stringGGI` (the genes with more than one protein) | `adj_matrix[!duplicated(newnames), !duplicated(newnames)]` keeps only the first protein of each gene. R drops the links of the other proteins (`R/get_stringPPI_to_stringGGI.R:62-64`). | `string_gene_network` joins the links of each protein of a gene. A gene pair gets the highest score of its links. |
| `get_stringPPI_to_stringGGI` (the proteins without a gene) | A protein that BioMart does not map keeps its protein ID as its name. Thus the result has gene IDs and protein IDs (`R/get_stringPPI_to_stringGGI.R:43-60`). R drops the text after a second `.` of the ID. | `string_gene_network` has only gene pairs. `GeneNetwork.unmapped_proteins` lists each protein without a gene. The protein ID is the text after the first `.`. An ID without a taxon prefix, and links that join no two genes, raise `InputError`. |
| `create_graph` (the seed counts) | `f[gene] = res$overall_matches[which(res$gene_id == gene)]` keeps the first count of a gene with more than one row and warns `number of items to replace is not a multiple of replacement length`. A count of `NA`, a negative count, and `Inf` stay in the graph. If no gene of `res` is in the network, each feature is 0, and R gives no message (`R/create_graph.R:13-23`). | `string_scores` takes the highest count of a gene. It raises `InputError` for a count that is NaN, infinite, or negative. If no gene of `de` is in the network, it warns with `SeedMatchWarning`. |
| `.plot_ecdf` (the p-value label) | The branches compare `p < .01 && p >= .001`, `p > .01`, and `p < .001`. No branch takes p = 0.01, thus R stops with `object 'pvalue_label' not found`. A NaN or NA p-value stops R with `missing value where TRUE/FALSE needed`. A DTS p-value of 20 of 2000 permutations is exactly 0.01 (`R/deseq_fc_ecdf.R:312-324`). | `plot_ecdf` gives the label `0.01` for p = 0.01, and `NaN` for NaN. The other labels are the labels of R. |
| `deseq_fc_ecdf`, `.plot_ecdf` (the lines and the legend) | A row of `res` in no gene list gets the group `NA`, and `stat_ecdf` draws a gray line for it (`R/deseq_fc_ecdf.R:100-121`). The count of a legend label is the length of its gene list, thus an absent gene and a repeated gene add to it. The counts are in the order of `gene_lists`, not in the order of `factor_order` (`:112`, `:121`). | `plot_ecdf` draws one line for each key of `groups`, in the order of the keys. The count of a label is the number of finite `log2FoldChange` values of its line. |
| `wass_dist` (the ECDF plot) | `stat_ecdf(geom = "line")` adds a point at `-Inf` with 0 and a point at `Inf` with 1, and `geom_line` joins them to the values. Thus the line goes from 0 at the left edge of the panel to the height of the first value. Its slope changes with the width of the panel. The areas use 0 below the first value (`R/wass_dist.R:70-71`, `:89-93`). | `plot_wasserstein` draws the interpolated ECDF of the areas. The line is 0 up to the first value, goes up to the height of the first value, and is 1 after the last value. |
| `kmer_hist` (the histogram) | `plot(graphics::hist(...))` draws the histogram two times, and `abline` draws the target line on page 2 only. `breaks = 50` gives the breaks of `pretty`, and the first bin also holds its lower break. For other counts, `pretty` gives more than 50 bins, or bins of a fractional width. A target that is not a k-mer gives `NA`, and R draws no line (`R/kmer_hist.R:30-34`). | `kmer_histogram` gives bins of w integer counts each, with at most `max_bins` bins. `plot_kmer_histogram` draws one Axes with the bars, the title, and the target line. A target that is not a k-mer raises `InputError`. |
| `plotAlignment` (the figure) | The figure names the guide line `miRNA`, but the data.frame names it `siRNA ` (`R/plotAlignment.R:63`, `:135-137`). | `plot_alignment` names the guide line `siRNA `, as the data.frame does. |
| `plotRangesDF`, `plotIRanges` (the bars) | `plotRangesDF` draws only the first 6 rows, because of `head()`, and it ignores `xlim`, `main`, and `col`. The y value of a row is its row name. Thus a subset keeps the old positions of its rows, and text row names give NA. The bar of a range of one base has the width 0, and `plotIRanges` skips such a range with a warning (`R/plotIRanges.R:3-21`). | `plot_ranges` draws one bar for each row, at the position of the row in the table. A bar goes from `start - 0.5` to `end + 0.5`, and each option has an effect. |
| `plot.PRROC` (the curves of `benchmark_preds`) | For one distinct score s, the color range goes from 0.9 s to s + 0.09 s, and for s = 0 each color is NA. For s <= 0, `image` stops R. An infinite score gives the color NA to its threshold, and R stops (`R/benchmark.R:32-40`). | `plot_benchmark` takes the color range from the finite thresholds, and an infinite threshold gets the color of the nearest end of the scale. Without two distinct finite thresholds, each segment gets the middle color, and the color bar goes from c - 0.5 to c + 0.5. |
| `.make_table` (the labels and the ID of the report table) | The labels of the windows are fixed texts: `Full Guide Strand (g2:g23)`, `18-mer (g2:g19)`, and `15-mer (g2:g19)` (`R/SeedMatchReport.R:158`, `:162`, `:166`). The window of `15mer` is g2 to g16 (`R/get_seed.R:121-123`), and a guide of another length than 23 bases has another `Full` window. `gt::gt(..., id = "mygt")` gives each table the element ID `mygt` (`:145`). | `report_html` gives the labels from `seed_window` for the guide. The table has the class `seedmatch-report` and no ID. |

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
  - a second run into one `dest` while the first run downloads raises `FetchError` and keeps the `.part` file of the first run
  - a `Content-Length` header that is not the registry size causes an error before the read
  - a URL or a redirect that is not `https` causes an error

### Plots

The tests of `seedmatch.plot` read the data of the artists: the points of each line, the bars, the texts, and the colors. They do not compare pixels. `tools/golden/R/plots.R` gives the layer data of ggplot2 and the cells of `gt`.

### Documentation examples

A test runs each code example in `README.md` and in `docs/`.

### Coverage

`pytest-cov` measures branch coverage. The floor is 95% for `src/seedmatch`, except `seedmatch.plot`. The review must give a reason for each line without coverage.

## Performance and bounds

- `SequenceSet` keeps all bases in one buffer. Memory is linear in the total length.
- The 2bit reader decodes by range, and the EnsDb reader filters in SQL. Thus the request, not the file, sets the size of a read.
- A mismatch search is O(N · m) at worst, with an early stop when the mismatch count passes the limit. An indel search is O(N · m) dynamic programming.
- The kernels use `nogil=True`. `full_search` runs its patterns on a thread pool of fixed width. The width is a parameter, with the CPU count as the default.
- `kmer_counts` works in blocks of 2^20 window starts, and its matrix has at most one entry for each window.
- `KmerCounts.to_frame` is dense and holds n · 4^k · 8 bytes.
- `kmer_gene_counts` multiplies a sparse gene indicator matrix with the matrix of `kmer_counts`. Its memory is O(the entries of the matrix + 4^k), because each k-mer has a label and a value. On the 18,614 canonical 3' UTRs of `rnor7` with k = 7, the call takes 0.11 seconds, and its peak memory is 360 MB. The product has an entry for each pair of a gene and a k-mer. For k = 12, the call takes about 3 seconds, and its peak memory is near 2 GB, because of the 16,777,216 labels. Machine: Apple M2 Max.
- `alignment_text` takes O(W) time and memory for a window of W bases. Each pad is at most the length L of the sequence, thus W is at most N + 2L for a guide of N bases. It takes O(n) to find the sequence ID among n sequences, thus a loop over m matches costs O(m · n).
- `match_distances` gives n(n − 1)/2 rows for a sequence with n matches, thus the result sets the memory. 10,000 matches in one sequence give 49,995,000 rows of 24 bytes. The function allocates the result one time. The work memory, without the result, is O(m) for m matches. The time is O(m log m + P) for P pairs.
- `e_value_null` and `wobble_bulge_null` call `shuffle_null`, which owns the draw, the pool, and the result. `e_value_null` runs `n_random` searches of the set on `workers` threads. Its memory is the result, genes × `n_random` × 4 bytes, and O(`workers` × (genes + sequences)) more, for example 80 MB for 20,000 genes and 1000 shuffles. `map_in_pool` takes the shuffles one at a time from a generator, and the pool holds at most 2 × `workers` + 1 of them. The result exists before the first draw, thus a table that the system cannot allocate raises `MemoryError` before any search.
- `e_values` compares blocks of 4096 rows of `de` with `null.index.get_indexer`, which it calls one time. The work memory is O(4096 × `n_random`), and the function makes no copy of the null. The time is O((R + G) × `n_random`) for R rows of `de` and G rows of the null. The check of the null reads each of its rows.
- `pr_curve` and `roc_curve` sort the scores one time. The time is O(n log n) and the memory is O(n) for n scores. The curve has at most 2K + 2 + max(100, P) rows for K distinct scores. For 1,000,000 distinct scores, `pr_curve` takes 0.18 seconds and 161 MB, and `roc_curve` takes 0.11 seconds and 104 MB.
- `benchmark_predictions` finds the gold standard with one `isin` over the gene IDs, and then it calls `pr_curve` and `roc_curve` one time each. The time is O(n log n) and the memory is O(n) for n rows. For 1,000,000 rows with distinct scores, it takes 0.33 seconds, and the peak memory above the start is 220 MB.
- `read_mirdb` reads 2^18 lines at a time, and it rejects a line of more than 1024 characters (`D-049`). Thus the work memory is one chunk and 16 bytes for each kept row. The 16 bytes are a 4-byte code for the ID, a 4-byte code for the RefSeq ID, and the 8-byte score. The time is O(n) for n lines. The result holds 40 bytes for each row (a pointer for each text column, and the score) and the distinct strings. The miRDB 6.0 file of `seedmatch-fetch mirdb` has 59,483,982 bytes of gzip. It has 6,831,595 lines, 7,086 miRNA IDs, and 147,959 RefSeq IDs. The longest line has 51 characters with its line end. For all species, the read takes 6.2 seconds, and the peak memory above the start is about 670 MB. This resident size varies by up to 100 MB between runs. The result is 273 MB and 8.4 MB of strings. The step that makes the table needs 1.2 times the result in allocated memory (`tracemalloc`). With `species="rno"` (200,619 rows), the read takes 5.8 seconds, and the peak is about 155 MB. Machine: Apple M2 Max, Python 3.12.
- `wobble_bulge_matches` runs 1 + `max_bulge` · max(0, 6 − `bulge_free`) searches of the set on `workers` threads, 13 for the defaults. `max_bulge` is at most 19993, thus no pattern is longer than the 20000 letters that a search accepts. `wobble_bulge_patterns` holds all patterns, thus its memory is O(`max_bulge`² · max(0, 6 − `bulge_free`)) letters, about 1.2 billion for `max_bulge=19993` and `bulge_free=0`. `wobble_bulge_matches` keeps the sequence index, the start, and the pattern number of each hit until it sorts all hits, thus the memory is O(the hits of all patterns). A sequence of `N` is the worst case, because each pattern matches at each position. Each hit takes about 70 bytes, and there is at most one hit for each pattern and base. With the defaults, 4,000,000 `N` take 3.6 GB. On the 18,614 canonical 3' UTRs of `rnor7`, the defaults give 66,374 rows. One search takes 0.16 seconds, and the call takes 0.40 seconds with `workers=None` and 2.1 seconds with `workers=1`. The peak memory above the start is about 20 MB. Machine: Apple M2 Max with 12 CPUs.
- `wobble_bulge_null` runs the pool over the shuffles. The patterns of one shuffle run one after the other in the thread of that shuffle, with no second pool. The time is the sum, over the shuffles, of the time of `wobble_bulge_matches` for the shuffle on one thread, divided by `workers`. A shuffle can have more sites than the guide, and then it takes more time. The memory is the result, genes × `n_random` × 4 bytes, and O(`workers` × (genes + sequences)) more for the pool. Each thread that runs holds the hits of one shuffle, about 70 bytes for each hit. A sequence of `N` is the worst case, and then the bound is `workers` times the bound of `wobble_bulge_matches`. On the 18,614 canonical 3' UTRs of `rnor7`, the guide takes 2.3 seconds on one thread and has 66,374 rows. A shuffle takes 2.1 to 3.6 seconds and has 11,799 to 176,430 rows. The 20 shuffles with the defaults take 7.5 to 8.4 seconds with `workers=None` and 57 to 63 seconds with `workers=1`. The peak memory above the start is 40 to 60 MB with `workers=None` and 23 to 29 MB with `workers=1`. The frame holds 2.6 MB: the table of counts, 18,614 × 20 × 4 bytes (1.4 MB), and the gene IDs of the index. The default 1000 shuffles take about 6 minutes at this speed. Machine: Apple M2 Max with 12 CPUs, shared with a coverage run.
- `compensatory_sites` runs `max_seed_mismatch` + 1 searches for the seed matches, and then one search for each window of the guide, on `workers` threads. A window is a part of the guide from `first_position`, with a length from `min_length` to `max_length`. A guide of 23 bases with the defaults has 55 windows, thus the call does 57 searches. `first_position=8` gives 66 windows. The number of windows grows with the square of the guide length, and `max_length` bounds it. The function raises `InputError` before the first search if the longest region is above the pattern limit of the search. A set without a seed match skips the window searches. The function pairs each region result when `map_in_pool` gives it. Two `np.searchsorted` calls on the sorted seed keys give the pairs of a region. The upper bound is the end of its sequence. The time is O((S + R) log S + P log P) for S seed matches, R region matches, and P pairs. A window has at most one match for each end, thus P is at most S · W · `max_distance` for W windows. With the defaults, this is 5,390 pairs for each seed match. The memory is O(S + P) and the results that wait in the pool. The function holds each pair until it drops the contained regions (`D-053`). On the 18,614 canonical 3' UTRs of `rnor7`, the defaults give 8,273 rows. The call takes 1.4 seconds with `workers=None` and 10 seconds with `workers=1`. The peak memory above the start is 12 to 17 MB. P grows if `min_length` is short, if `max_region_mismatch` or `max_seed_mismatch` is high, or if the sequence has low complexity. In these cases the peak is from 150 to 390 bytes for each pair, and a large result adds to it. The guide of 23 `U` on a sequence of 2,000 `A` gives 10,437,735 pairs and 189,243 rows. The peak is 1.6 GB, and the call takes 1.8 seconds with `workers=1`. The default guide on the `rnor7` UTRs with `min_length=3` and `max_region_mismatch=1` gives 7,851,633 pairs and 1,900,173 rows. The peak is 2.2 GB, and the call takes 6.9 seconds with `workers=None`. Machine: Apple M2 Max with 12 CPUs, shared with a coverage run.
- `read_trrust` reads 2^16 lines at a time, and it rejects a line of more than 1024 characters (`D-057`). Thus the work memory is one chunk and 16 bytes for each kept row. The 16 bytes are four 4-byte codes. The four columns share one vocabulary, thus the result holds one string for each distinct text and 32 bytes for each row. The time is O(n) for n lines. The TRRUST human file of `seedmatch-fetch trrust` has 297,659 bytes and 9,396 lines. The longest line has 273 characters with its line end. The read takes 0.008 seconds, and the peak memory above the start is about 8 MB. A synthetic file of 2,000,000 lines with distinct PubMed IDs takes 2.7 seconds and about 440 MB. Machine: Apple M2 Max, Python 3.12.
- `trrust_scores` factorizes the genes with one `pd.factorize` call, and it sorts the edge keys with `np.unique`. The time is O(R log R) for R rows, and the memory is O(R + G) for G genes. `mean_messages` uses `np.bincount`, thus a channel takes O(E + G) time for E edges. The TRRUST human file has 9,396 rows and 2,862 genes. The call takes 0.003 seconds, and the peak memory above the start is about 1 MB. On the synthetic file of 2,000,000 rows and 100,000 genes, it takes 0.63 seconds and about 150 MB. Machine: Apple M2 Max, Python 3.12.
- `read_string_links` reads 2^16 lines at a time, and it rejects a line of more than 1024 characters (`D-059`). Thus the work memory is one chunk and 10 bytes for each kept row. The 10 bytes are two 4-byte codes and a 2-byte score, and the two protein columns share one vocabulary. The time is O(n) for n lines. `r_quantile` partitions a float64 copy of the kept scores with `numpy.ndarray.partition`, thus it takes O(k) time and 8 bytes for each of k kept rows. The result holds 24 bytes for each link and the distinct strings. The STRING rat file of `seedmatch-fetch string 10116` has 761,435,565 bytes of text and 14,100,659 lines. The longest line has 54 characters with its line end. With the defaults, the read takes 8.6 seconds, and the peak memory above the start is about 95 MB. The result has 98,920 links and 8,809 proteins. With `score_threshold=0` and `percentile=0`, the result has 13,894,314 links, the read takes 10.1 seconds, and the peak is about 660 MB. Machine: Apple M2 Max, Python 3.12.
- `string_gene_network` factorizes the protein IDs of the links with one `pd.factorize` call. It asks the EnsDb for the genes of the distinct proteins with one query on the index `protein_protein_id_idx` (`D-013`). The IDs go to SQLite as one JSON array, thus their number has no limit of SQLite. The function sorts the keys `gene1 * G + gene2` for G genes one time, and `np.maximum.reduceat` takes the highest score of each pair (`D-060`). The time is O(E log E) for E links, and the memory is O(E). The default links of the STRING rat file have 98,920 links and 8,809 proteins. With the EnsDb 109 of `rnor7`, the result has 84,654 rows and 7,979 genes. 628 proteins have no gene, and 96 genes have more than one protein (`D-061`, `D-060`). The first call takes 0.2 seconds, and a later call takes 0.04 seconds. The peak memory above the start is about 25 MB. With `score_threshold=0` and `percentile=0`, the links have 13,894,314 rows. The call takes 3.6 seconds, and the result has 11,097,016 rows. The peak memory above the start is about 1.8 GB, or 130 bytes for each link. Machine: Apple M2 Max, Python 3.12.
- `string_scores` factorizes the genes of the network with one `pd.factorize` call. The check of the pairs sorts the keys `gene1 * G + gene2` for G genes one time, and it finds the repeated pairs and the reverse rows. `np.maximum.at` takes the highest count of each gene of `de`, and `mean_messages` uses `np.bincount`. The time is O(E log E + D) for E rows of the network and D rows of `de`, and the memory is O(E + G + D). The default network of the STRING rat file has 84,654 rows and 7,979 genes. With a `de` of 20,000 rows, the call takes 0.015 seconds, and the peak memory above the start is about 5 MB. With `score_threshold=0` and `percentile=0`, the network has 11,097,016 rows and 18,223 genes. The call takes 1.8 seconds, and the peak memory above the start is about 630 MB, or 57 bytes for each row. Machine: Apple M2 Max, Python 3.12.
- `plot_ecdf` sorts the finite values of each group with `numpy.unique`. The time is O(n log n) for n rows. The work memory is O(n), and the figure holds two floats for each distinct value. For 22,000 rows in two groups, the call takes 0.014 seconds, and 0.06 seconds with `Figure.savefig` to PNG at 100 dpi. The peak memory above the start of the two is 3.3 MiB. Machine: Apple M2 Max, Python 3.12.
- `wasserstein_difference` sorts each sample with `numpy.unique`, and `np.searchsorted` reads the grid. The time is O(n log n + `n_interp` log n). The result holds 32 bytes for each grid point. The function does not read the grid in chunks, thus the peak memory is about twice the result. For 20,000 and 2,000 values and `n_interp` = 10^6, the call takes 0.04 seconds, and the peak is 61 MiB. Machine: Apple M2 Max, Python 3.12.
- `plot_wasserstein` sums the Gaussian kernels of each sample at 512 grid points, in blocks of 1024 values. The work memory is 512 · 1024 · 8 bytes, 4 MiB, and the time is O(512 · n). The other Axes take O(`n_interp` + n log n). The figure holds 512 points for each density and `n_interp` points for the difference. For 20,000 and 2,000 values and `n_interp` = 10000, the call takes 0.053 seconds, and 0.14 seconds with `Figure.savefig` to PNG at 100 dpi. The peak memory of `tracemalloc` is 5.4 MiB. Machine: Apple M2 Max, Python 3.12.
- `kmer_histogram` takes the lowest and the highest count, and one `np.bincount` counts the bins. The time is O(n) for n k-mers, and the memory is O(n + `max_bins`). For the 16,384 counts of k = 7, the call takes 0.1 ms. For the 4^12 counts of k = 12, the call takes 0.05 seconds. The peak above the input is 128 MiB, one temporary array of n values. Machine: Apple M2 Max, Python 3.12.
- `plot_kmer_histogram` calls `kmer_histogram`, and one lookup in the index finds the target. The first lookup in an index builds its hash table, which takes O(n) time and about 32 bytes for each of n labels. For the 4^10 labels of k = 10, the first call takes 0.06 seconds and 33 MiB above the input. A later call takes 0.009 seconds. The figure holds one bar for each bin. Machine: Apple M2 Max, Python 3.12.
- `plot_alignment` puts the lines of the text, each after its name, in one `Figure.text`. The time and the memory are O(W) for a label whose longest line has W characters. The label is the text with the names. The figure is 0.1 inch wide for each character, thus `Figure.savefig` at 100 dpi draws a canvas of 3,600 bytes for each character. For the label of the first `mer7m8` match of the test data with `left_pad=20` and `right_pad=20`, W is 84. The call takes 0.1 ms, and 13 ms with `Figure.savefig` to PNG at 100 dpi. For W = 2,000, `Figure.savefig` takes 0.27 seconds. Machine: Apple M2 Max, Python 3.12.
- `plot_ranges` puts the bars of all rows in one `PolyCollection`, with 4 vertices for each row. The time and the memory are O(n) for n rows. The function has no row limit, because a limit drops rows without a message, as `head()` of R does. For 10^5 ranges, the call takes 0.065 seconds, 0.65 µs for each row. With `Figure.savefig` to PNG at 100 dpi, it takes 0.24 seconds, 2.4 µs for each row. The peak of `tracemalloc` is 40.7 MiB, 427 bytes for each row. Machine: Apple M2 Max, Python 3.12.
- `plot_benchmark` puts the segments of each curve in one `LineCollection`, with a segment from each row to the next row. The time and the memory are O(R) for R rows of the two curves. The PR curve has at most 2K + 2 + max(100, P) rows for K distinct scores and P positives. The function has no row limit, because a limit drops rows without a message. For 10^5 distinct scores (200,003 rows), the call takes 0.32 seconds, 1.6 µs for each row. With `Figure.savefig` to PNG at 100 dpi, it takes 0.71 seconds, 3.5 µs for each row. The peak of `tracemalloc` is 69.1 MiB, 362 bytes for each row, and 71.4 MiB, 374 bytes for each row, with `savefig`. For 10^6 scores of 8 values with 30% positives (262,411 rows), the call takes 0.41 seconds, and 0.89 seconds with `savefig`. The peak is 94.7 MiB, 379 bytes for each row, and 98.3 MiB, 393 bytes for each row, with `savefig`. Machine: Apple M2 Max, Python 3.12.
- `report_html` writes the table in one pass over the rows of `report.report`. A dict of lists groups the row positions by group, and each row takes one f-string for each cell. The time and the memory are O(R) for R rows. For 10^5 rows in 25,000 groups, the call takes 1.0 seconds, 10 µs for each row. The peak of `tracemalloc` is 119 MiB, and the result is 46 MiB. The report of the test data has 4 rows, and the call takes 0.3 ms. Machine: Apple M2 Max, Python 3.12.
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