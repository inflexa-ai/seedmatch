# Statistics: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port the statistics of SeedMatchR (`ecdf_stat_test`, `wass_dist`, and the groups of `deseq_fc_ecdf`) to `seedmatch.stats`, with output identical to R except for the recorded divergences.

**Architecture:** `seedmatch.stats` has one module for each test: KS, Wilcoxon rank-sum, DTS, and the Wasserstein areas. Each test takes two samples as arrays and returns a frozen result. `ecdf_test` selects the `log2FoldChange` values of two gene sets from a DE table, and it calls one test. `assign_groups` gives the group of each row of a DE table. The plot of `deseq_fc_ecdf` is in a later plan.

**Tech Stack:** Python 3.12, numpy, pandas, scipy 1.16.3 (`scipy.special` only). R 4.5.3 with twosamples 2.0.1 and caTools 1.18.3 for the goldens.

**Spec:** `docs/design.md`, sections "Parity contract", "Statistics rules", and "Divergences".

**Later plans (not in this plan):** the ECDF plot and the p-value label of `.plot_ecdf`, and `benchmark_preds` (PR and ROC curves) with the other analyses: closed. `docs/plans/2026-11-21-plot-and-report.md` ports the plots.

## Global Constraints

- Runtime dependencies: `numpy`, `pandas`, `numba`, and `scipy`. Task 1 adds `scipy>=1.14` to `dependencies` and `scipy==1.16.3` to `constraint-dependencies`. No other new dependency. The code uses only `scipy.special` (`kolmogorov`, `ndtr`).
- No core module imports a network module or `seedmatch.fetch` (ruff `TID251`). No function writes a file.
- Errors derive from `seedmatch.errors.SeedMatchError`. Invalid input raises `InputError`, and each message names the bad value and the expected value.
- A float result agrees with R within a relative tolerance of 1e-9. A KS statistic also agrees within an absolute tolerance of 1e-12, because R adds 1/m and -1/n in a cumulative sum. A Wasserstein area that R rounds is equal.
- Each change of behavior from R gets an entry `## D-NNN: <title>` in `docs/divergences.md`. It also gets a row in the table of `docs/design.md`, and a test whose name contains the ID. This plan adds `D-032` to `D-036`.
- A rename of an argument or a value, or a change of a type, is not a divergence. `docs/r-map.md` records it.
- The R source of truth is `../SeedMatchR` at commit `2fc0c4f`. The R scripts run with `tools/golden/run.sh`. The tests never run R and never use the network.
- Licenses (`docs/design.md`, "Licenses"): write each algorithm from its definition or its paper. The R output then makes sure that the result is the same. Do not translate GPL source line by line. The controller read `ks.c` and `fround.c` of R 4.5.3, and this plan states the behavior that they show.
- Documents use ASD-STE100. Code comments obey the admission test of the project rules. The word "seam" does not occur.
- Git: an implementer does no `git add` and no commit. Each task ends with a checkpoint. The controller commits.

## Shared interfaces

`src/seedmatch/stats/types.py` (Task 2) defines these names. Each later task uses them:

```python
class Alternative(StrEnum):
    TWO_SIDED = "two.sided"
    LESS = "less"
    GREATER = "greater"

class PValueMethod(StrEnum):
    EXACT = "exact"
    ASYMPTOTIC = "asymptotic"
    PERMUTATION = "permutation"

@dataclass(frozen=True)
class TwoSampleResult:
    statistic: float
    p_value: float
    method: PValueMethod
    alternative: Alternative | None  # None for DTS
    n_x: int
    n_y: int

def as_sample(name: str, values: ArrayLike, *, finite: bool) -> npt.NDArray[np.float64]: ...
```

`as_sample` returns a new 1-D float64 array. It raises `InputError` for an empty sample, for a value that is not a number, and for NaN. With `finite=True`, it also raises for an infinite value. Each message gives `name`, and for a bad value its index.

## Review Focus

1. The default call of R on a DE table without a filter: some `log2FoldChange` values are NaN. Each test drops these rows, and the result gives the counts after the drop (`D-036`). Owner: Task 6.
2. A target set that is a pandas Series of gene IDs, or a set, or a bare `str`. The Series and the set work. A bare `str` raises `InputError`, because it is one ID and not a set. Owner: Task 6.
3. A target set that has no row in the DE table: `InputError` that names the set and gives the count 0. R stops with `not enough 'x' data`. Owner: Task 6.
4. A large DE table (20000 genes) with DTS: the run takes seconds, not minutes, and the memory of the permutations is bounded. Owner: Task 4.
5. A Wasserstein area that is near a rounding boundary, for example 0.0025: the R rule of `round` gives 0.002, and Python `round` gives 0.003. Owner: Task 5.

---

### Task 1: R goldens of the statistics, and the scipy dependency

**Files:**
- Make: `tools/golden/R/stats.R`
- Make: `tests/golden/stats/samples.tsv`, `ks.tsv`, `wilcoxon.tsv`, `dts.tsv`, `wasserstein.tsv`, `round.tsv`, `ecdf.tsv` (written by the script)
- Modify: `tools/golden/README.md`, `pyproject.toml`, `uv.lock`

**Interfaces:**
- Produces: the goldens below. Each later task reads them with `tests/helpers.py::golden`.

- [ ] **Step 1 (controller): Add scipy.** Add `"scipy>=1.14"` to `dependencies` and `"scipy==1.16.3"` to `constraint-dependencies` in `pyproject.toml`. Run `uv lock` and `uv sync`. Expected: `uv run python -c "import scipy; print(scipy.__version__)"` prints `1.16.3`.

- [ ] **Step 2: Write `tools/golden/R/stats.R`.** The script sources `common.R` and writes each golden with `write_golden`. The notes of each golden give the versions of twosamples and caTools. Each case has an ID. `samples.tsv` has the columns `case`, `sample` (`x` or `y`), and `value` (17 significant digits), one row for each value. The script makes the samples with `set.seed(20261120)`. A case with ties rounds its values to 1 decimal. These cases are necessary:
  - KS: the sizes (1, 12), (5, 7), (30, 40), (99, 101), (100, 100), (120, 150), and (400, 300), each with and without ties. 99 × 101 = 9999 gives the exact p-value, and 100 × 100 = 10000 gives the asymptotic p-value. Also a case where each value is equal, and a case with a large shift (D near 1) for each size above 100.
  - KS two-sided asymptotic: one case with sqrt(en)·D between 0.9 and 1 (the cut series of R), and one case with a one-sided p-value below 1e-17 (R gives 0).
  - Wilcoxon: the sizes (1, 5), (10, 12), (49, 49), (49, 50), and (60, 80), each with and without ties. Also a case where each value is equal.
  - DTS: the sizes (5, 7), (30, 40), and (200, 300), with and without ties.
  - Wasserstein: the sizes (10, 10), (10, 50), (40, 300), with and without ties, each with `n.interp` 2, 100, and 10000. Also a seed sample with only one distinct value. R stops, and the notes record the message.
- [ ] **Step 3: Write the result goldens.**
  - `ks.tsv` and `wilcoxon.tsv`: `case`, `alternative`, `statistic`, `p_value`, `method` (`exact` or `asymptotic`, from the `method` text of R), for each alternative.
  - `dts.tsv`: `case`, `statistic`, `p_value`, `n_boot` (2000). The p-value comes from `twosamples::dts_test` after `set.seed(1)`.
  - `wasserstein.tsv`: `case`, `n_interp`, `abs_auc`, `neg_cens_auc`, `pos_cens_auc`.
  - `round.tsv`: `x`, `digits`, `r`. The rows are `round(x, 3)` for each `x` in `c(0.0005, 0.0015, 0.0025, 0.0035, 0.1235, 2.675, 1.0005, -0.0025, -0.1235)`, for each `x` in `(0:4000)/4000 + 0.0005`, and for 200 values of `runif(200, 0, 3)`. Also `round(x, 2)` and `round(x, 0)` for the first list.
  - `ecdf.tsv`: `input`, `test`, `alternative`, `statistic`, `p_value`, `n_target`, `n_background`. `input` is `test_data` or `readme`. For `test_data`, the table is `load_test_data()$res`, with the target `res$gene_id[res$mer7m8 >= 1]` and the background `res$gene_id[res$mer7m8 == 0]`, as in `tests/testthat/helper.R` of SeedMatchR. For `readme`, the table comes from the R steps that the first README example ports (`filter_res` with the defaults of `filter_de`, then the `mer7m8` counts of `SeedMatchR`). `tools/golden/R/de_and_counts.R` has these steps. Write a row for KS and Wilcoxon with each alternative, and for Wass and DTS (`set.seed(1)`) with the default alternative.
- [ ] **Step 4: Run the script two times.** Run `tools/golden/run.sh tools/golden/R/stats.R` two times. Expected: the second run changes no byte. In `ecdf.tsv`, the `test_data` KS `greater` row has `statistic` 0.51035… and `p_value` 0.0001037731…, and the Wilcoxon `less` row with the swapped sets has `W` 2150. These are the R snapshots of SeedMatchR. Write the row with the swapped sets too.
- [ ] **Step 5: Document the goldens.** Add a section "Golden files of the statistics" to `tools/golden/README.md`: the command, the cases, and the columns. Checkpoint.

### Task 2: The shared types and the KS test

**Files:**
- Make: `src/seedmatch/stats/__init__.py`, `src/seedmatch/stats/types.py`, `src/seedmatch/stats/ks.py`
- Make: `tests/stats/test_ks.py`, `tests/stats/test_types.py`
- Modify: `docs/divergences.md` (`D-032`), `docs/design.md` (the KS rules and the table row)

**Interfaces:**
- Consumes: `ks.tsv` and `samples.tsv` of Task 1.
- Produces: the names of "Shared interfaces" above, and `ks_test(x: ArrayLike, y: ArrayLike, *, alternative: Alternative = Alternative.TWO_SIDED) -> TwoSampleResult`.

Rules of `ks_test`. The default alternative is the default of `stats::ks.test`.
- The samples go through `as_sample(..., finite=False)`. Let m = n_x and n = n_y.
- The statistic comes from the integer lattice. Sort the pooled values. A checkpoint is each index k from 1 to m + n with one of these conditions: k = m + n, or the value at k is less than the value at k + 1. At a checkpoint, let i and j be the counts of x and y values up to k. Then d = max(i·n − j·m) for `GREATER`, max(j·m − i·n) for `LESS`, and max(|i·n − j·m|) for `TWO_SIDED`, over the checkpoints. The statistic is d / (m·n). A d below 1 gives the p-value 1.
- The exact p-value applies if m·n < 10000, also with ties. It is P(D ≥ d) under the permutation distribution: all C(m + n, m) choices of the x positions in the pooled order have the same probability. A path is a choice. It fails at a checkpoint (i, j) if the bound of the alternative reaches d there. Compute the upper tail directly with the recursion of Viehmann (2021): C(i, j) = 1 at a point that fails, else (i·C(i−1, j) + j·C(i, j−1)) / (i + j), with C(0, 0) = 0. Thus a small p-value keeps its relative precision. For `LESS`, swap the roles of x and y. The work is O(m·n) and m·n < 10000.
- The asymptotic p-value applies if m·n ≥ 10000. Let en = m·n / (m + n) and D = d / (m·n). `TWO_SIDED` gives `scipy.special.kolmogorov(sqrt(en) * D)`. `GREATER` and `LESS` give `exp(-2 * en * D**2)`.
- The p-value is clamped to [0, 1].

- [ ] **Step 1: Write the failing tests.**
  - `test_ks_matches_r(case, alternative)`: parametrized over the rows of `ks.tsv` that are `exact`, and over the `asymptotic` rows that are not the two `D-032` cases. Assert the statistic and the p-value with the tolerances of the Global Constraints, and assert the method.
  - `test_exact_p_value_equals_the_count_of_the_permutations(samples)`: a hypothesis test with m + n ≤ 12 and values from a list of 4 numbers, thus with ties. The oracle lists each choice of the x positions, and it computes d for each choice. The p-value is the fraction with a d not below the observed d. Assert equality within 1e-12 for each alternative.
  - `test_ks_asymptotic_two_sided_uses_the_limit_distribution_D032` and `test_ks_asymptotic_one_sided_keeps_a_small_p_value_D032`: for the two R cases, assert the port value against `scipy.special.kolmogorov` and `exp(...)` within 1e-12. Also assert that the R value differs: by more than 1e-7 relative in the first case, and R gives 0 in the second case.
  - `test_ks_rejects_an_empty_sample` and `test_ks_rejects_nan`: `InputError`, `match` the sample name.
  - `test_ks_accepts_an_infinite_value`: the R golden case or a direct assertion that `inf` counts as the greatest value.
- [ ] **Step 2: Run the tests.** Run `uv run pytest tests/stats -q`. Expected: FAIL with `ModuleNotFoundError: No module named 'seedmatch.stats'`.
- [ ] **Step 3: Write `types.py` and `ks.py`.** Export the names of "Shared interfaces" and `ks_test` from `seedmatch.stats`.
- [ ] **Step 4: Run the tests.** Expected: PASS. Also `uv run ruff check`, `uv run ruff format --check`, `uv run mypy`.
- [ ] **Step 5: Write `D-032`.** "The asymptotic KS p-values are accurate." R location: `stats::psmirnov` through `ks.test` (`R/deseq_fc_ecdf.R` calls `stats::ks.test`). R behavior: for x < 1, the two-sided limit series of R 4.5.3 keeps one term (`src/library/stats/src/ks.c:100-111`). The one-sided value is `1 - (-expm1(-2 * en * D^2))`, which gives 0 for a p-value below about 1e-16. Evidence: the two R cases of `ks.tsv`. Update "KS test" in `docs/design.md`, and add the row to the divergence table. Checkpoint.

### Task 3: The Wilcoxon rank-sum test

**Files:**
- Make: `src/seedmatch/stats/wilcoxon.py`, `tests/stats/test_wilcoxon.py`
- Modify: `src/seedmatch/stats/__init__.py`, `docs/divergences.md` (`D-033`), `docs/design.md`

**Interfaces:**
- Consumes: `as_sample`, `Alternative`, `PValueMethod`, `TwoSampleResult` (Task 2), and `wilcoxon.tsv` (Task 1).
- Produces: `wilcoxon_test(x: ArrayLike, y: ArrayLike, *, alternative: Alternative = Alternative.TWO_SIDED) -> TwoSampleResult`.

Rules:
- The samples go through `as_sample(..., finite=False)`. The ranks are the average ranks of the pooled values. W is the rank sum of x minus n_x(n_x + 1)/2.
- The exact p-value applies if n_x < 50, n_y < 50, and no value ties. Let P(W ≤ w) be the distribution of the rank sum. `LESS` gives P(W ≤ W_obs). `GREATER` gives P(W ≥ W_obs). `TWO_SIDED` gives min(1, 2·p), where p = P(W ≥ W_obs) if W_obs > n_x·n_y/2, else P(W ≤ W_obs). The counts of the distribution can pass 2^63. Use float64 counts, or Python integers, and divide by C(n_x + n_y, n_x).
- Otherwise, the normal approximation applies, with the tie correction and the continuity correction of `docs/design.md`. If each pooled value is equal, the variance is 0, and the p-value is 1 for each alternative (`D-033`).

- [ ] **Step 1: Write the failing tests.**
  - `test_wilcoxon_matches_r(case, alternative)`: each row of `wilcoxon.tsv` except the case where each value is equal.
  - `test_exact_p_value_equals_the_count_of_the_subsets(samples)`: a hypothesis test without ties and with n_x + n_y ≤ 12. The oracle lists each subset of the ranks.
  - `test_wilcoxon_with_all_values_equal_gives_1_D033`: the R golden row has NaN for `two.sided`, and the port gives 1 for each alternative.
  - The rejection tests of Task 2, for `wilcoxon_test`.
- [ ] **Step 2: Run the tests.** Expected: FAIL with an `ImportError` for `wilcoxon_test`.
- [ ] **Step 3: Write `wilcoxon.py`.** Use `scipy.special.ndtr` for the normal distribution.
- [ ] **Step 4: Run the tests.** Expected: PASS, and ruff and mypy pass.
- [ ] **Step 5: Write `D-033`.** "The Wilcoxon test gives the p-value 1 when each value is equal." R gives NaN for `two.sided`, because z is 0/0. Update "Wilcoxon rank-sum test" in `docs/design.md`, and add the table row. Checkpoint.

### Task 4: The DTS test

**Files:**
- Make: `src/seedmatch/stats/dts.py`, `tests/stats/test_dts.py`
- Modify: `src/seedmatch/stats/__init__.py`, `docs/divergences.md` (`D-034`), `docs/design.md`

**Interfaces:**
- Consumes: the Task 2 types and `dts.tsv` (Task 1).
- Produces: `dts_test(x: ArrayLike, y: ArrayLike, *, rng: np.random.Generator | int, n_boot: int = 2000) -> TwoSampleResult` with `method=PValueMethod.PERMUTATION` and `alternative=None`.

Rules:
- The samples go through `as_sample(..., finite=True)` (`D-034`). An `int` seed becomes `np.random.default_rng(seed)`. `n_boot` is an `int` of 1 or more. A bool raises `InputError`.
- The statistic is the formula of `docs/design.md`, "DTS test", with p = 1.
- Each permutation draws a random choice of n_x positions of the pooled sorted values for x. The p-value is the fraction of the `n_boot` permutation statistics that are not less than the observed statistic. If the fraction is 0, the p-value is 1 / (2·`n_boot`), as twosamples does.
- Draw the permutations in blocks of a fixed size, for example 64. Thus the memory is O(64·(n_x + n_y)), and numpy does the work of each block.

- [ ] **Step 1: Write the failing tests.**
  - `test_dts_statistic_matches_r(case)`: each row of `dts.tsv`, the statistic only.
  - `test_dts_statistic_equals_the_definition(samples)`: a hypothesis test against an oracle that evaluates the three ECDFs with `np.searchsorted` at the left value of each gap.
  - `test_dts_statistic_is_symmetric(samples)`: swapping x and y gives the same statistic.
  - `test_dts_p_value_agrees_with_r_within_the_permutation_error(case)`: with `rng=1`, |p − p_R| ≤ 4·sqrt(2·p_R·(1 − p_R)/2000) + 2/2000.
  - `test_dts_p_value_is_the_same_for_the_same_seed` and `test_dts_p_value_is_never_0`.
  - `test_dts_rejects_an_infinite_value_D034` and `test_dts_rejects_an_empty_sample_D034`: `InputError` with a `match`.
  - `test_dts_runs_on_a_large_table`: n_x = 2000 and n_y = 18000 with `n_boot=2000` finishes in less than 30 s on CI. Give the time in the report.
- [ ] **Step 2: Run the tests.** Expected: FAIL with an `ImportError` for `dts_test`.
- [ ] **Step 3: Write `dts.py`.**
- [ ] **Step 4: Run the tests.** Expected: PASS, and ruff and mypy pass.
- [ ] **Step 5: Write `D-034`.** "The DTS test rejects an infinite value and an empty sample." Evidence (controller R probe): `twosamples::dts_test(c(a, Inf), b)` gives the statistic `Inf` with p 1, and `dts_test(numeric(0), b)` gives a statistic with p 1. Update "DTS test" in `docs/design.md`, and add the table row. Checkpoint.

### Task 5: The Wasserstein areas and the R rounding rule

**Files:**
- Make: `src/seedmatch/_rmath.py`, `src/seedmatch/stats/wasserstein.py`
- Make: `tests/test_rmath.py`, `tests/stats/test_wasserstein.py`
- Modify: `src/seedmatch/stats/__init__.py`, `docs/design.md`

**Interfaces:**
- Consumes: `as_sample` (Task 2), `wasserstein.tsv` and `round.tsv` (Task 1).
- Produces: `r_round(x: float, digits: int) -> float` in `seedmatch._rmath`. The later plot plan uses it for the labels.
- Produces: `WassersteinAreas` (frozen: `absolute: float`, `negative_censored: float`, `positive_censored: float`) and `wasserstein_areas(seed: ArrayLike, non_seed: ArrayLike, *, n_interp: int = 10000) -> WassersteinAreas`.

The R rule of `round(x, d)` for d > 0 applies in R 4.0.0 and later. `?Round` says that R "measures" and rounds to the nearest even:
- NaN, an infinite value, and 0 stay the same. A negative x is rounded as -x, and the sign is put back.
- If log10(x) + d is greater than 15, x stays the same.
- Let p = 10^d, xd = floor(x·p)/p, and xu = ceil(x·p)/p, in float64. The result is xu if xu − x < x − xd, or if the two differences are equal and floor(x·p) is odd. Otherwise, the result is xd.
- d = 0 rounds to the nearest integer, and a tie goes to the even integer.

Rules of `wasserstein_areas` (the algorithm of `R/wass_dist.R:61-81`):
- Each sample goes through `as_sample(..., finite=True)`. A sample with fewer than 10 values raises `InputError`, as in R. A sample with fewer than 2 distinct values raises `InputError`. R stops in `stats::approx`. `n_interp` is an `int` of 2 or more.
- The ECDF of each sample is evaluated at its own values. Merge the tied values to one point (`approx`, `ties = mean`).
- The grid is `n_interp` points from the minimum to the maximum of the pooled values, as `seq(from, to, length.out)`: the first point is the minimum, the last point is the maximum, and point k is min + k·((max − min)/(n_interp − 1)).
- Interpolate each ECDF on the grid with the expression y_i + (y_j − y_i)·((v − x_i)/(x_j − x_i)), with 0 below the first point and 1 above the last point.
- diff = seed ECDF − non-seed ECDF. The three areas integrate max(diff, 0), max(−diff, 0), and |diff| with the trapezoid rule, and `r_round` rounds each to 3 decimals.

- [ ] **Step 1: Write the failing tests.**
  - `test_r_round_matches_r(row)`: each row of `round.tsv`, with exact equality.
  - `test_r_round_differs_from_python_round`: `r_round(0.0025, 3) == 0.002` and `round(0.0025, 3) == 0.003`. Also `r_round(0.1235, 3) == 0.124`.
  - `test_wasserstein_matches_r(case, n_interp)`: each row of `wasserstein.tsv`, with exact equality of the three areas.
  - `test_the_censored_areas_swap_with_the_samples`: swapping seed and non-seed swaps `negative_censored` and `positive_censored`, and `absolute` does not change.
  - `test_wasserstein_rejects_fewer_than_10_values`, `test_wasserstein_rejects_one_distinct_value`, and `test_wasserstein_rejects_n_interp_below_2`: `InputError` with a `match`.
- [ ] **Step 2: Run the tests.** Expected: FAIL with an `ImportError`.
- [ ] **Step 3: Write `_rmath.py` and `wasserstein.py`.**
- [ ] **Step 4: Run the tests.** Expected: PASS, and ruff and mypy pass.
- [ ] **Step 5: Update "Wasserstein areas" in `docs/design.md`.** Give the rounding rule and the two `InputError` cases. Checkpoint.

### Task 6: `ecdf_test`, `assign_groups`, and the documents

**Files:**
- Make: `src/seedmatch/stats/ecdf.py`, `tests/stats/test_ecdf.py`
- Modify: `src/seedmatch/stats/__init__.py`, `src/seedmatch/__init__.py`, `src/seedmatch/errors.py`, `tests/test_readme.py`, `README.md`, `docs/r-map.md`, `docs/divergences.md` (`D-035`, `D-036`), `docs/design.md`

**Interfaces:**
- Consumes: `ks_test`, `wilcoxon_test`, `dts_test`, `wasserstein_areas`, and the Task 2 types. `check_columns` and `check_numeric_columns` in `seedmatch._checks`. `ecdf.tsv` (Task 1).
- Produces:

```python
class EcdfTest(StrEnum):
    KS = "KS"
    WILCOXON = "Wilcoxon"
    WASSERSTEIN = "Wass"
    DTS = "DTS"

@dataclass(frozen=True)
class EcdfTestResult:
    test: EcdfTest
    statistic: float
    p_value: float | None        # None for WASSERSTEIN
    method: PValueMethod | None  # None for WASSERSTEIN
    alternative: Alternative | None
    n_target: int
    n_background: int
    areas: WassersteinAreas | None  # only for WASSERSTEIN

def ecdf_test(
    de: pd.DataFrame,
    *,
    target: Iterable[str],
    background: Iterable[str],
    test: EcdfTest = EcdfTest.KS,
    alternative: Alternative | None = None,
    rng: np.random.Generator | int | None = None,
    n_boot: int = 2000,
    n_interp: int = 10000,
) -> EcdfTestResult: ...

def assign_groups(de: pd.DataFrame, groups: Mapping[str, Iterable[str]]) -> pd.Series: ...

class SeedMatchWarning(UserWarning): ...  # in seedmatch.errors
```

Rules of `ecdf_test` (`ecdf_stat_test`, `R/deseq_fc_ecdf.R`):
- `de` has the columns `gene_id` and `log2FoldChange`, and `log2FoldChange` is numeric.
- `target` and `background` are each an iterable of gene IDs. A bare `str` raises `InputError`. The target values are the `log2FoldChange` of each row whose `gene_id` is in `target`, in the row order of `de`. A gene ID with more than one row gives each of its rows, as `%in%` does. The same applies to `background`.
- Each test drops a row whose `log2FoldChange` is NaN (`D-036`). `n_target` and `n_background` are the counts after the drop. A count of 0 raises `InputError` that names the set.
- `alternative=None` gives `GREATER`, the default of R, for KS and Wilcoxon. An `alternative` with DTS or Wasserstein raises `InputError` (`D-035`). The result has `alternative=None` for these two tests.
- DTS needs `rng`. `None` raises `InputError` that names `rng`.
- `rng` and `n_boot` apply only to DTS, and `n_interp` applies only to Wasserstein. A value other than the default for another test raises `InputError`.
- For Wasserstein, the seed sample is the target and the non-seed sample is the background. `statistic` is `areas.negative_censored`, as in R.

Rules of `assign_groups` (`.assign_groups` and the `factor` of `deseq_fc_ecdf`):
- The result is a Series named `group`, with the index of `de`. Its dtype is categorical, and the categories are the keys of `groups` in their order.
- A row whose `gene_id` is in no group has no value. If a gene is in more than one group, the last group wins, as in the loop of R.
- If two groups share a gene, `assign_groups` warns with `SeedMatchWarning`. The message names the two groups and one shared gene. R warns with `Gene lists overlap!`.

- [ ] **Step 1: Write the failing tests.**
  - `test_ecdf_test_matches_r(row)`: each row of `ecdf.tsv`. For DTS, compare only the statistic, and the p-value within the error of Task 4.
  - `test_ecdf_test_drops_rows_without_a_fold_change_D036`: a table with NaN values. KS gives the R result of the table without these rows, and DTS and Wasserstein give a result (R stops).
  - `test_ecdf_test_rejects_an_alternative_for_dts_and_wasserstein_D035`.
  - `test_ecdf_test_rejects_a_bare_string`, `test_ecdf_test_rejects_a_set_without_rows`, `test_ecdf_test_needs_rng_for_dts`, `test_ecdf_test_rejects_options_of_another_test`, `test_ecdf_test_rejects_a_missing_column`.
  - `test_assign_groups_gives_the_last_group_and_warns_for_an_overlap`, and `test_assign_groups_keeps_the_order_of_the_groups`.
  - `test_readme_ecdf_example_gives_the_r_result` in `tests/test_readme.py`: the README example gives the `readme` KS `greater` row of `ecdf.tsv`.
- [ ] **Step 2: Run the tests.** Expected: FAIL with an `ImportError` for `ecdf_test`.
- [ ] **Step 3: Write `ecdf.py`.** Export `EcdfTest`, `EcdfTestResult`, `ecdf_test`, `assign_groups`, `Alternative`, `WassersteinAreas`, and `wasserstein_areas` from `seedmatch`, and `SeedMatchWarning` from `seedmatch` and `seedmatch.errors`.
- [ ] **Step 4: Write the README example.** After the report example, add a Python block. It compares the `mer7m8` genes with the genes without a match, with `ecdf_test(de, target=..., background=...)`. Then it prints the statistic and the p-value. Add a "Differences" bullet for `D-035` and `D-036`, and add the statistics to "Scope".
- [ ] **Step 5: Run the tests.** Run `uv run pytest tests/stats tests/test_readme.py -q`. Expected: PASS, and ruff and mypy pass.
- [ ] **Step 6: Write the documents.**
  - `docs/r-map.md`: rows for `ecdf_stat_test`, `deseq_fc_ecdf` (the statistics and the groups, without the plot), `.check_gene_list_overlap`, `.assign_groups`, and `wass_dist`. Give the renames: `stats_test` becomes `test`, `Wilcoxen` becomes `EcdfTest.WILCOXON`, `abs.auc`, `neg.cens.auc`, and `pos.cens.auc` become the fields of `WassersteinAreas`, and `n.interp` becomes `n_interp`.
  - `D-035`: "ecdf_test rejects an alternative for DTS and Wasserstein." R passes `alternative` only to `ks.test` and `wilcox.test` (`R/deseq_fc_ecdf.R`), thus R ignores it for the other tests without a message.
  - `D-036`: "ecdf_test drops the rows without a fold change for each test." R drops them for KS and Wilcoxon only. Evidence (controller R probe): `dts_test` stops with `missing value where TRUE/FALSE needed`, and `wass_dist` stops with `'from' must be a finite number`.
  - `docs/design.md`: the two table rows, and the public API example uses the final names. Checkpoint.

### Task 7: Verification of the plan

- [ ] **Step 1: Run the checks one time.** Run `uv run ruff check`, `uv run ruff format --check`, `uv run mypy`, and `NUMBA_DISABLE_JIT=1 uv run pytest --cov`. Expected: each passes, and the coverage is not below 95%. Give a reason for each line of `seedmatch.stats` without coverage.
- [ ] **Step 2: Compare the documents with the code.** Each public name of `seedmatch.stats` has a row in `docs/r-map.md` or is a type of a function that has a row. Each `D-NNN` of this plan has an entry, a table row, and a test. Checkpoint.
