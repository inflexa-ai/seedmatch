# Non-canonical searches: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port the non-canonical searches of SeedMatchR to `seedmatch.noncanonical`: the wobble and bulge search of `run_wobble_bulges_pipeline`, and the 3' compensatory search of `get_compensatory_sequences`. In R, the first function stops for each input, and the second function returns 0 rows for each input. The port gives the corrected semantics of `D-004` and `D-005`.

**Architecture:** `seedmatch.noncanonical` is a package with a module for each search. `wobble.py` makes the patterns of the wobble and bulge search, finds their matches, and counts the matches of each gene. `compensatory.py` pairs each seed match with the 3' compensatory regions of the same sequence. Each module calls the kernels of the port through `seedmatch.search`.

Task 2 adds the shared helpers. A new function `match_ranges` gives the raw ranges with the index of each sequence, and `find_matches` uses it. The gene counts of `annotate_seed_counts` move to helpers in `report/annotate.py`, and the wobble counts use the same helpers.

No new kernel is necessary. Each wobble and bulge pattern is an IUPAC pattern of fixed length, and each compensatory region is a pattern of fixed length. The mismatch search of `match_ranges` finds each of them.

**Tech Stack:** Python 3.12, numpy, pandas. R 4.5.3 with Biostrings 2.78.0, dplyr 1.2.1, and gtools 3.9.5 for the goldens.

**Spec:** `docs/design.md`, sections "Parity contract", "Search rules", "Divergences (defect fixes)" (the rows `get_compensatory_sequences` and `run_wobble_bulges_pipeline`), and "Performance and bounds".

**Later plans (not in this plan):**

- The network scores (`seedmatch.network`).
- The plots and the report HTML.
- A null of `e_value_null` for the wobble and bulge counts. `e_value_null` counts the canonical seed matches of each shuffle, thus `e_values` with the column `overall_matches` compares two different counts.

## Global Constraints

- Runtime dependencies: `numpy`, `pandas`, `numba`, and `scipy`. No new dependency.
- No core module imports a network module or `seedmatch.fetch` (ruff `TID251`). No function writes a file or draws a plot.
- Errors derive from `seedmatch.errors.SeedMatchError`. Invalid input raises `InputError`, and each message names the bad value and the expected value.
- The message templates of this plan are normative. Each error test matches a stable part of its template. A template that names a helper of `seedmatch._checks` is the message of that helper.
- The R source of truth is `../SeedMatchR` at commit `2fc0c4f`. The R scripts run with `tools/golden/run.sh`. The tests never run R and never use the network.
- Each change of behavior from R gets an entry `## D-NNN: <title>` in `docs/divergences.md` and a row in the table of `docs/design.md`. It also gets a test whose name or docstring contains the ID.
- This plan writes the entries of the design rows `D-004` and `D-005`, and of the new IDs `D-052` to `D-056`. Put each entry in ID order.
- A rename of an argument, a column, or a value, or a change of a type, is not a divergence. `docs/r-map.md` records it.
- The search obeys "Search rules" of `docs/design.md`. Coordinates are 1-based and inclusive.
- SeedMatchR has the MIT license, thus a golden patch can copy a SeedMatchR function.
- Each function obeys "Performance and bounds" of the design. Its docstring gives the time and the memory, and its task adds a line to that section of `docs/design.md`.
- Documents use ASD-STE100. Code comments obey the admission test of the project rules. The word "seam" does not occur.
- Git: an implementer does no `git add` and no commit. Each task ends with a checkpoint. The controller commits.

## Shared interfaces

Task 2 makes each of these names, and Task 3 uses them:

```python
# seedmatch/search/api.py. find_matches builds its frame from the same arrays.
def match_ranges(
    pattern: str | Target, sequences: SequenceSet, options: SearchOptions = SearchOptions()
) -> tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.int64]]: ...

# seedmatch/_checks.py
def check_type(name: str, value: object, expected: type) -> None: ...

# seedmatch/report/annotate.py
def check_count_column(de: object, name: str) -> pd.DataFrame: ...
def add_gene_counts(
    de: pd.DataFrame,
    sequences: SequenceSet,
    counts: NDArray[np.int64],
    *,
    name: str,
    shared_genes: bool,
) -> pd.DataFrame: ...
```

- `match_ranges` returns the index of the sequence in the set, the start, and the end of each match. The arrays have the row order of `find_matches`. `find_matches` returns the same rows, thus a sequence ID that occurs more than once does not hide which sequence has the match.
- `check_type` raises `f"{name} has type {type(value).__name__}; expected a {expected.__name__}"` unless `isinstance(value, expected)`.
- `check_count_column` does the checks of `annotate_seed_counts` on the DE table, in their present order: the type, the column `gene_id`, the name `gene_id`, a repeated column (`D-051`), and text gene IDs. It returns the table. The messages do not change.
- `add_gene_counts` gives each gene the maximum of `counts` over its sequences (`D-003`). It then does the `shared_genes` rule and the column rule of `annotate_seed_counts` (`D-015`). `counts` has one value for each sequence of the set.

## Review Focus

1. The canonical 3' UTRs of `rnor7` (18,614 sequences, 17,913,962 bases). One exact scan of this set takes about 0.15 seconds on the planner machine. The default wobble search makes 13 scans, and the default compensatory search makes 57 scans. Give the wall time and the peak memory of each function with `workers=None` and `workers=1`. Task 2 and Task 3 own this item.
2. A run of `N` in a sequence. The wobble search uses `fixed = FALSE` of R for the pattern and the subject, thus a subject `N` matches each pattern letter, and 11 `N` give 5 rows. The compensatory search uses `fixed = TRUE` of R, thus a subject `N` is a mismatch. Owners: Task 2 (`test_an_n_run_matches_each_pattern`) and Task 3 (`test_an_n_in_the_seed_site_is_a_mismatch`).
3. A region at the end of one sequence and a seed match at the start of the next sequence. The pair search uses positions in the buffer of the set, and a window must stop at the end of its sequence. Owner: Task 3 (`test_a_window_stops_at_the_end_of_its_sequence`).
4. A guide that is too short for the region, for example a guide of 13 bases with the default options. The region g9 to g13 has 5 positions, which is less than `min_length` 6. Owner: Task 3 (`test_a_short_guide_has_no_room_for_the_region_D055`).
5. An empty set, an empty sequence, and a sequence shorter than each pattern. Each function gives an empty frame with the columns and the dtypes of the result. Owners: Task 2 and Task 3.

---

### Task 1: R goldens of the non-canonical searches, and this plan

**Files:**
- Add: `docs/plans/2026-11-20-noncanonical.md` (this plan)
- Make: `tools/golden/patches/get_compensatory_sequences.R`, `tools/golden/R/noncanonical.R`
- Make: `tests/data/compensatory_cases.tsv`, `tests/data/wobble_bulge_cases.tsv`
- Make: `tests/golden/noncanonical/compensatory.tsv` and `tests/golden/noncanonical/wobble_bulge_r.tsv` (written by the script)
- Modify: `tools/golden/README.md`, `tests/data/README.md`

**Interfaces:**
- Produces: the goldens below. Tasks 2 and 3 read them with `tests/helpers.py::golden("noncanonical/<file>")` and `golden_note`.

- [ ] **Step 1: Write the patch.** `tools/golden/patches/get_compensatory_sequences.R` is a copy of `get_compensatory_sequences` (`R/get_compensatory_sequences.R:26-133`) with these changes. Its first comment names `D-004` and the changed R lines. The planner copy is `tmp/claude/plan-noncanonical/get_compensatory_sequences.R`, and `diff` against the original shows only these changes:
  - Remove the test `if (nrow(compensatory_regions_captured) == 0)`, its early return, the `else {`, and the brace that closes the `else` (`:67-77`, `:130`). Keep the body of the `else` branch (`:78-129`). The frame of `:57-66` is always empty, thus the test is always true.
  - At `:99`, write `Biostrings::reverseComplement(Biostrings::DNAString(seed))`. `gsub` at `:34` makes `seed` a character string, and `reverseComplement` has no method for it. A planner probe without this change stops at the first match with `unable to find an inherited method for function 'reverseComplement' for signature 'x = "character"'`.
- [ ] **Step 2: Write the fixtures.** Each fixture is synthetic, and a tab separates the fields. `tests/data/compensatory_cases.tsv` has these lines. The column `limit` is the R argument. `NA` in `max_length` is the R default.

```text
case	guide	min_length	max_length	limit	seed_mismatch	seq_id	sequence
seed_past_end	UUAUAGAGCAAGAACACUGUUUU	6	NA	100	1	t1	GGGGGTTCTTGCGGGGGGGGGGCTCTAT
no_site	UUAUAGAGCAAGAACACUGUUUU	6	NA	100	1	t1	GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
repeated_name_r	UUAUAGAGCAAGAACACUGUUUU	6	NA	100	1	g	GGGGGTTCTTGCGGGGGGGGGGCTCTATAGGGGG
repeated_name_r	UUAUAGAGCAAGAACACUGUUUU	6	NA	100	1	g	GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTCTTGCGGGGGGGGGGGGGGGGGGGGCTCTATAGGGGG
two_seed_matches_r	UUAUAGAGCAAGAACACUGUUUU	6	NA	100	1	t1	GGGGGGTTCTTGCTCTATAGGGGGGGGGGCTCTATAGGGGG
repeated_region_r	UUAUAGAGCAAGCAAGCAAGUUU	6	NA	100	1	t1	GGGGGAGCTTGCAGGGGGGGGGGCTCTATAGGGGG
region_at_g2_r	UUAUAGAUAUAGACACUGUUUUA	6	6	100	1	t1	GGGGGCTCTATACGGGGGGGGGGATCTATAGGGGG
seed_mm2_r	UUAUAGAGCAAGAACACUGUUUU	6	NA	100	2	t1	GGGGGTTCTTGCGGGGGGGGGGCTGTAGAGGGGG
min_above_max_r	UUAUAGAGCAAGAACACUGUUUU	6	3	100	1	t1	GGGGGTTCTTGCGGGGGGGGGGCTCTATAGGGGG
```

  `tests/data/wobble_bulge_cases.tsv` has these lines:

```text
case	guide	sequence
perfect	UUAUAGAGCAAGAACACUGUUUU	GGGGGGGGCTCTATAGGGGGGGG
wobble_g4	UUAUAGAGCAAGAACACUGUUUU	GGGGGGGGCTCTGTAGGGGGGGG
mismatch_t4	UUAUAGAGCAAGAACACUGUUUU	GGGGGGGGCTCCATAGGGGGGGG
bulge_after_g5	UUAUAGAGCAAGAACACUGUUUU	GGGGGGGGCTCATATAGGGGGGG
bulge_after_g4	UUAUAGAGCAAGAACACUGUUUU	GGGGGGGGCTCTGATAGGGGGGG
bulge_4_after_g7	UUAUAGAGCAAGAACACUGUUUU	GGGGGGGGCAAAATCTATAGGGG
bulge_5_after_g7	UUAUAGAGCAAGAACACUGUUUU	GGGGGGGGCAAAAATCTATAGGG
one_start_two_patterns	UUAUAGAGCAAGAACACUGUUUU	GGGGGGGGCTCTATATATAGGGG
five_wobbles	AUUUUUAGCAAGAACACAGAUAA	AAAACTGGGGGAAAA
```

- [ ] **Step 3: Write `tools/golden/R/noncanonical.R`.** The script sources `common.R` and writes each golden with `write_golden` into `tests/golden/noncanonical/`. The planner draft is `tmp/claude/plan-noncanonical/noncanonical_draft.R`, and it writes the same files into the scratch directory. Read each fixture with `read.delim(..., colClasses = "character", na.strings = character(0))`, thus the text `NA` stays text. The sections:
  - **Compensatory, test data.** Run the patched function with `with_seedmatchr_patch` on `load_test_data()$seqs`, `GUIDE`, and `seed = get_seed(GUIDE, <seed>)$Target.Seq`. The cases are in the table below. Keep the 15 R columns of `:116-128` in their order, after a first column `case`.
  - **Compensatory, synthetic cases.** For each case of `compensatory_cases.tsv`, make a `DNAStringSet` with the names `seq_id`. Run the patched function with the arguments of the case and the seed `mer7m8` of its guide. If R stops, add the note `<case>: the patched function stops with: <message>`.
  - Write all rows to `compensatory.tsv`. The notes give the dplyr version, the patch, and these R results:
    - `original_vignette`: the number of rows of the unpatched function for the vignette case, and its message. Expected: 0 rows and `No compensatory regions found.`
    - `line_99`: the error of `Biostrings::reverseComplement(gsub("U", "T", get_seed(GUIDE, "mer7m8")$Target.Seq))`.
    - `repeated_name_r`: the number of times that `seqs[c("g", "g")]` holds the first sequence. Expected: 2.
    - `relative_position_seq_g6`: `SeedMatchR:::relative_position_seq("TATAGA", "UUAUAGAUAUAGACACUGUUUUA")`. Expected: start `2,8` and end `7,13`.
    - `min_17`: the substrings that the loop of `:52-54` makes for `compensatory_region_length_min = 17`. Expected: the 16 letters of g8 to g23, three times. The loop `17:16` gives the starts 1 and 0 for the length 17, and the start 1 for the length 16.
    - `gc_content`: `.gc_content("GCAT")` and `.gc_content("")`. Expected: 50 and `NaN`.
  - **Wobble and bulge.** The notes of `wobble_bulge_r.tsv` give these R results. Wrap each call of the pipeline in `utils::capture.output`, because it prints a line also with `silent = TRUE` (`:22`).
    - `pipeline_error` and `pipeline_error_bulges_lim_1`: the error of `run_wobble_bulges_pipeline(GUIDE, res = d$res, features = d$seqs, silent = TRUE)`, with the default `bulges_lim` and with 1. Expected: `object 'bulge_allowance' not found`.
    - `wobble_seeds_guide`, `wobble_seeds_4_u`, `wobble_seeds_5_u_in_seed`: the `seed_combinations` of `.get_all_combination_w_b` for `GUIDE`, for `AUAUAGAGCAAGAACACAGAUAU`, and for `AUUUUUAGCAAGAACACAGAUAA`.
    - `wobble_seeds_2_u` and `wobble_seeds_0_u`: the error of `.get_all_combination_w_b` for `AUAGAGAGCAAGAACACAGAUAA` and for `CAGAGAGCAAGAACACAGACAAC`.
    - `formals_get_target_bulged_guides`: the formals of `.get_target_bulged_guides`.
    - `check_bulges_nt`: `check_bulges_nt(RNAStringSet(c("CNUCUAUA", "CUCUANUA", "CUCNUAUA")), 3)`. Expected: `CUCUANUA` and `CUCNUAUA`.
  - **Wobble and bulge, the R count rule.** R stops before it counts. To show the count rule of `:35-55`, the script sets `bulge_allowance` in the global environment before each call of `.sequence_with_tb`. It removes the variable after the call, and it gives the arguments in the documented order. For each case of `wobble_bulge_cases.tsv`, the columns `r_rule_bulges_lim_1` and `r_rule_bulges_lim_4` hold the sum over the wobble seeds and the strings of `check_bulges_nt(..., 3)` of `vcountPattern(..., fixed = FALSE)`. Write `wobble_bulge_r.tsv` with the columns `case`, `guide`, `sequence`, `r_rule_bulges_lim_1`, and `r_rule_bulges_lim_4`.

| case | R arguments |
|---|---|
| `vignette` | `seed` `mer7m8`, and the defaults: `compensatory_region_length_min = 6`, `compensatory_region_length_max = NA`, `limit = 100`, `mismatch_ = 1`, `mismatch_compensatory_region_ = 0` |
| `seed_mm0` | `mer7m8`, `mismatch_ = 0` |
| `min8_max10` | `mer7m8`, `compensatory_region_length_min = 8`, `compensatory_region_length_max = 10` |
| `limit40` | `mer7m8`, `limit = 40` |
| `region_mm1` | `mer7m8`, `compensatory_region_length_min = 9`, `mismatch_compensatory_region_ = 1` |
| `mer6` | `seed` `mer6`, and the defaults |

- [ ] **Step 4: Run the script two times.** Run `tools/golden/run.sh tools/golden/R/noncanonical.R` two times. Expected: the second run changes no byte. The planner draft gives these results.
  - The row counts of `compensatory.tsv` are `vignette` 104, `seed_mm0` 5, `min8_max10` 10, `limit40` 39, `region_mm1` 43, `mer6` 454, `seed_past_end` 1, `repeated_name_r` 4, `two_seed_matches_r` 1, `repeated_region_r` 2, `region_at_g2_r` 1, `seed_mm2_r` 1, and `min_above_max_r` 3.
  - `no_site` has a note and no row.
  - In `wobble_bulge_r.tsv`, the table below gives the values of `r_rule_bulges_lim_1` and `r_rule_bulges_lim_4`.

| case | `r_rule_bulges_lim_1` | `r_rule_bulges_lim_4` |
|---|---|---|
| `perfect` | 12 | 24 |
| `wobble_g4` | 8 | 20 |
| `mismatch_t4` | 4 | 12 |
| `bulge_after_g5`, `bulge_after_g4`, `bulge_4_after_g7`, `bulge_5_after_g7` | 0 | 0 |
| `one_start_two_patterns` | 12 | 24 |
| `five_wobbles` | 3 | 25 |

- [ ] **Step 5: Document the goldens.**
  - Add a section "Golden files of the non-canonical searches" to `tools/golden/README.md`. Give the command, the patch and its changes, the cases, the notes, and the R count rule.
  - Add each fixture of Step 2 to `tests/data/README.md`, and state that it is synthetic. Give the purpose of each case. A case that ends with `_r` shows the R behavior that a divergence changes.
  - Checkpoint.

### Task 2: The wobble and bulge search

**Files:**
- Make: `src/seedmatch/noncanonical/__init__.py`, `src/seedmatch/noncanonical/wobble.py`, `tests/noncanonical/test_wobble.py`
- Modify: `src/seedmatch/search/api.py`, `src/seedmatch/search/__init__.py`, `src/seedmatch/_checks.py`, `src/seedmatch/report/annotate.py`, `src/seedmatch/__init__.py`, `tests/search/test_api.py`, `tests/test_checks.py`, `tests/helpers.py` (`run_core_pipeline`), `tests/test_readme.py`, `README.md`, `docs/r-map.md`, `docs/divergences.md` (`D-005`), `docs/design.md`

**Interfaces:**
- Consumes: `Guide`, `Seed`, `SearchOptions`, `SequenceSet`, `as_limit`, `check_not_negative`, `pool_width`, `map_in_pool`, and the goldens of Task 1.
- Produces: the names of "Shared interfaces", and:

```python
@dataclass(frozen=True)
class WobbleBulgeOptions:
    max_bulge: int = 4
    bulge_free: int = 3

def wobble_bulge_patterns(
    guide: Guide, options: WobbleBulgeOptions = WobbleBulgeOptions()
) -> pd.DataFrame: ...

def wobble_bulge_matches(
    guide: Guide,
    sequences: SequenceSet,
    options: WobbleBulgeOptions = WobbleBulgeOptions(),
    *,
    workers: int | None = None,
) -> pd.DataFrame: ...

def annotate_wobble_bulge_counts(
    de: pd.DataFrame,
    sequences: SequenceSet,
    guide: Guide,
    *,
    column: str = "overall_matches",
    options: WobbleBulgeOptions = WobbleBulgeOptions(),
    shared_genes: bool = True,
    workers: int | None = None,
) -> pd.DataFrame: ...
```

Rules of the options (`bulges_lim` and `nt` of R):
- `max_bulge` is the largest number of extra target bases in one bulge. It goes through `as_limit` and `check_not_negative`: `max_bulge -1 is negative; expected 0 or more`. The value 0 gives the search without a bulge.
- `bulge_free` is the number of seed positions from g2 that no bulge is next to. It goes through `as_limit`. A value outside 0 to 7 raises `bulge_free 8 is out of range; expected an integer from 0 to 7`.
- A bool raises the message of `as_limit`, for example `max_bulge True has type bool; expected an integer`.

Rules of the patterns (the row `run_wobble_bulges_pipeline` of the design):
- The seed is `Seed.MER7M8`, as in R (`R/run_wobble_bulges_pipeline.R:113`). Let t be `guide.target(Seed.MER7M8, wobble=True).sequence`. Thus t has `R` at each target position that pairs with a guide `U`, and one pattern covers each set of G:U pairs.
- Target position k (1 to 7) pairs with the guide position g(9 − k). Thus t1 pairs with g8, and t7 pairs with g2.
- A bulge is b extra target bases between the guide positions g_i and g_(i+1). Its pattern is `t[: 8 - i] + "N" * b + t[8 - i :]`. `bulge_after` holds i, and `bulge_size` holds b.
- The values of i are `bulge_free` + 2 to 7. Thus no bulge is next to g2 to g(`bulge_free` + 1). The default 3 gives the gaps g5|g6, g6|g7, and g7|g8. The values of b are 1 to `max_bulge`. A site has one bulge at most.
- The first pattern is t, with `bulge_after` 0 and `bulge_size` 0. Then the bulged patterns come in the order of b, and then of i.
- `wobble_bulge_patterns` returns the columns `pattern` (text), `bulge_after`, and `bulge_size` (`int64`), with one row for each pattern in this order and a `RangeIndex`. The number of rows is 1 + `max_bulge` · max(0, 6 − `bulge_free`).

Rules of the matches:
- Each pattern goes to `match_ranges` with `SearchOptions(fixed_pattern=False, fixed_subject=False)`. This is `fixed = FALSE` of R (`R/run_wobble_bulges_pipeline.R:39`, `:42`). Thus a subject `N` matches each pattern letter.
- The searches run on `map_in_pool` with `pool_width(workers)`.
- A sequence has a row for each distinct start of a match of a pattern. The row takes the first pattern in pattern order that matches at this start. Its `end` and `width` come from this pattern.
- The columns are `seq_id`, `gene_id`, `start`, `end`, `width`, `bulge_after`, and `bulge_size`. The integer columns are `int64`. The rows come in the order of the set, and then of `start`, with a `RangeIndex`.
- Sort the triples (sequence index, start, pattern number) of all hits with `np.lexsort`, and keep the first triple of each pair (sequence index, start). The memory is O(the hits of all patterns).
- `annotate_wobble_bulge_counts` counts the rows of each sequence with `np.bincount` over the sequence index. Then it calls `add_gene_counts`. The checks of `check_count_column` run before the searches. Its default column is `overall_matches`, the column of R (`:51-55`).
- The time is 1 + `max_bulge` · max(0, 6 − `bulge_free`) searches of the set, 13 for the defaults, on `workers` threads.
- `guide` goes through `check_type(..., Guide)`, `sequences` through `check_type(..., SequenceSet)`, and `options` through `check_type(..., WobbleBulgeOptions)`.

The expected rows of the cases, with the default options. A row is (`start`, `end`, `bulge_after`, `bulge_size`). The planner computed them with a regex oracle of these rules:

| case | rows |
|---|---|
| `perfect` | (9, 15, 0, 0) |
| `wobble_g4` | (9, 15, 0, 0) |
| `mismatch_t4` | none |
| `bulge_after_g5` | (9, 16, 5, 1) |
| `bulge_after_g4` | none. With `bulge_free=2`: (9, 16, 4, 1) |
| `bulge_4_after_g7` | (9, 19, 7, 4) |
| `bulge_5_after_g7` | none. With `max_bulge=5`: (9, 20, 7, 5) |
| `one_start_two_patterns` | (9, 15, 0, 0). Three patterns match at 9. |
| `five_wobbles` (its own guide) | (5, 11, 0, 0) |

- [ ] **Step 1: Write the failing tests** in `tests/noncanonical/test_wobble.py`, `tests/search/test_api.py`, and `tests/test_checks.py`. Read the cases with `golden_rows("noncanonical/wobble_bulge_r.tsv")`.
  - `test_match_ranges_give_the_rows_of_find_matches(options)`: parametrized over `SearchOptions()`, `SearchOptions(max_mismatch=1)`, and `SearchOptions(max_mismatch=1, indels=True)` on `sirna_seqs()` with the `mer7m8` target. `seq_ids[index]`, `start`, and `end` equal the columns of `find_matches`.
  - `test_check_type_names_the_argument_and_the_type`: `check_type("guide", "UUAU", Guide)` raises `guide has type str; expected a Guide`.
  - `test_wobble_bulge_patterns_of_the_default_guide`: the rows for `GUIDE` are these triples, in this order: `("CTCTRTR", 0, 0)`, `("CTCNTRTR", 5, 1)`, `("CTNCTRTR", 6, 1)`, `("CNTCTRTR", 7, 1)`, `("CTCNNTRTR", 5, 2)`, `("CTNNCTRTR", 6, 2)`, `("CNNTCTRTR", 7, 2)`, `("CTCNNNTRTR", 5, 3)`, `("CTNNNCTRTR", 6, 3)`, `("CNNNTCTRTR", 7, 3)`, `("CTCNNNNTRTR", 5, 4)`, `("CTNNNNCTRTR", 6, 4)`, `("CNNNNTCTRTR", 7, 4)`.
  - `test_wobble_bulge_patterns_obey_the_options(max_bulge, bulge_free)`: parametrized over (0, 3), (1, 0), (2, 5), (2, 6), and (3, 7). The number of rows is 1 + `max_bulge` · max(0, 6 − `bulge_free`). The set of `bulge_after` values above 0 is `set(range(bulge_free + 2, 8))` if `max_bulge` is above 0.
  - `test_wobble_bulge_matches_of_the_cases_D005(case)`: parametrized over the rows of the golden, with the guide of the row. The rows equal the table above. The docstring gives the note `pipeline_error`: R stops for each call.
  - `test_one_site_counts_once_D005`: `annotate_wobble_bulge_counts(pd.DataFrame({"gene_id": ["perfect"]}), SequenceSet.from_records([("perfect", <sequence>)]), Guide(GUIDE))` gives `overall_matches` 1. The golden gives 12 and 24 for the R count rule.
  - `test_a_mismatch_is_not_a_bulge_D005`: `mismatch_t4` gives no row. The R rule gives 4 and 12, because its "bulge" replaces a seed letter with `N` (`R/run_wobble_bulges_pipeline.R:197`).
  - `test_an_extra_target_base_is_a_bulge_D005`: `bulge_after_g5` and `bulge_4_after_g7` give one row each. The R rule gives 0.
  - `test_no_bulge_is_next_to_g2_to_g4_D005`: `bulge_after_g4` gives no row with the defaults, and (9, 16, 4, 1) with `bulge_free=2`. The docstring gives the note `check_bulges_nt`. R drops an `N` after t1, between g8 and g7. R keeps an `N` after t3 or t5, between g6 and g5 or between g4 and g3.
  - `test_a_bulge_above_max_bulge_is_not_a_site`: `bulge_5_after_g7` gives no row with the defaults, and (9, 20, 7, 5) with `max_bulge=5`.
  - `test_each_set_of_wobbles_is_in_one_pattern_D005`: the first pattern of `CAGAGAGCAAGAACACAGACAAC` (no U) is `GCTCTCT`, and of `AUAGAGAGCAAGAACACAGAUAA` (2 U) is `CTCTCTR`. `five_wobbles` gives (5, 11, 0, 0). The docstring gives the notes `wobble_seeds_0_u`, `wobble_seeds_2_u`, and `wobble_seeds_5_u_in_seed`: R stops for fewer than 4 U, and no R seed has `R` at the 5 U positions.
  - `test_each_distinct_start_counts_once_D005`: the rows of `ENSRNOG00000016275` in `sirna_seqs()` are (69, 79, 6, 4) and (73, 79, 0, 0). `one_start_two_patterns` gives one row.
  - `test_max_bulge_0_gives_the_matches_of_the_wobble_target`: with `max_bulge=0`, the columns `seq_id` to `width` equal `find_matches` of `Guide(GUIDE).target(Seed.MER7M8, wobble=True)` with `SearchOptions(fixed_pattern=False, fixed_subject=False)`. The fixture gives 91 rows. `bulge_after` and `bulge_size` are 0 in each row.
  - `test_each_canonical_match_is_a_row_without_a_bulge`: each of the 20 rows of `find_matches(Guide(GUIDE).target(Seed.MER7M8), sirna_seqs())` has a row with the same `seq_id`, `start`, and `end`, and `bulge_size` 0.
  - `test_the_fixture_gives_the_counts_of_the_oracle`: the defaults on `sirna_seqs()` give 757 rows, in 145 sequences. The counts of `bulge_size` 0 to 4 are 91, 189, 162, 159, and 156.
  - `test_wobble_bulge_matches_equal_an_oracle(guide, records, max_bulge, bulge_free)`: a hypothesis test with `deadline=None`. The guide has 8 to 23 letters of `ACGU`. The records have 0 to 4 sequences of 0 to 40 letters of `ACGTNRY`, and half of them hold the target of the guide. `max_bulge` is 0 to 3, and `bulge_free` is 0 to 7.
  - In the oracle test, the oracle builds the patterns from the rules above. A pattern matches at start s if it ends in the sequence and each of its letters matches the subject letter. Two letters match if their masks in `seedmatch.sequence.alphabet.MASKS` share a bit. The oracle gives, for each start in the union, the first pattern in pattern order.
  - `test_an_n_run_matches_each_pattern`: the sequence `NNNNNNNNNNN` gives the starts 1 to 5, each with `bulge_size` 0.
  - `test_empty_and_short_sequences_give_an_empty_frame`: an empty set, and a set of `("e", "")` and `("s", "CTCTAT")`. The result has the columns and the dtypes of "Rules of the matches", and no row.
  - `test_annotate_wobble_bulge_counts_give_the_gene_maximum_D003`: `SequenceSet.from_records([("a", <sequence of perfect>), ("b", "GGCTCTATAGGCTCTATAGG")], gene_ids=["g", "g"])` and a DE table with the gene `g` give the count 2.
  - `test_annotate_wobble_bulge_counts_equal_the_rows_of_each_gene`: on `sirna_de()` and `sirna_seqs()`, the value of each gene is the number of its rows in `wobble_bulge_matches`.
  - `test_wobble_bulge_matches_do_not_depend_on_workers`: `workers=1` and `workers=4` give equal frames.
  - `test_wobble_bulge_options_reject_a_negative_max_bulge`, `test_wobble_bulge_options_reject_bulge_free_above_7`, `test_wobble_bulge_options_reject_a_bool`, `test_wobble_bulge_matches_reject_a_guide_of_another_type`, `test_wobble_bulge_matches_reject_sequences_of_another_type`, and `test_wobble_bulge_matches_reject_options_of_another_type`: `InputError` with a `match` on the template.
- [ ] **Step 2: Run the tests.** Run `uv run pytest tests/noncanonical/test_wobble.py tests/search/test_api.py tests/test_checks.py -q`. Expected: FAIL with `ModuleNotFoundError: No module named 'seedmatch.noncanonical'` and `ImportError` for `match_ranges` and `check_type`.
- [ ] **Step 3: Write the shared helpers.**
  - Add `match_ranges`, and build `find_matches` from it. Export `match_ranges` from `seedmatch.search`, not from `seedmatch`.
  - Add `check_type` to `_checks.py`. Then make the inline type checks of `kmers.py` (`sequences` and `counts`) and of `evalue.py` (`options`) call it. Their messages do not change.
  - Move the checks and the gene counts of `annotate_seed_counts` to `check_count_column` and `add_gene_counts`.
  - Run `uv run pytest tests/report tests/search -q`. Expected: PASS without a change of a test.
- [ ] **Step 4: Write `wobble.py`.** Export `WobbleBulgeOptions`, `wobble_bulge_patterns`, `wobble_bulge_matches`, and `annotate_wobble_bulge_counts` from `seedmatch.noncanonical` and from `seedmatch`. `annotate_wobble_bulge_counts` has 7 parameters, thus it carries `# noqa: PLR0913` with the comment of `annotate_seed_counts`. Add a `wobble_bulge_matches` call to `run_core_pipeline`.
- [ ] **Step 5: Run the tests.** Expected: PASS. Also run `uv run ruff check`, `uv run ruff format --check`, and `uv run mypy`.
- [ ] **Step 6: Measure.** Do Review Focus 1 for `wobble_bulge_matches`. Give the numbers in the report of the task and in "Performance and bounds".
- [ ] **Step 7: Write the README part and the documents.**
  - `README.md`, "Other analyses": a paragraph and a Python block after the miRDB block. The block calls `wobble_matches = wobble_bulge_matches(guide, sequences)` and `wobble_de = annotate_wobble_bulge_counts(de, sequences, guide)`, and prints the head of each. The paragraph gives the pattern rule, `max_bulge`, `bulge_free`, the distinct-start rule, and `D-005`.
  - `tests/test_readme.py`: add `test_readme_wobble_example_gives_the_count_of_the_rows`. `wobble_matches` has 757 rows, and the `overall_matches` of each gene of `wobble_de` is the number of its rows.
  - `README.md`, "Scope": add the wobble and bulge search to the first list, and remove it from the second list. Add a bullet for `D-005` to "Differences".
  - `docs/r-map.md`, row `run_wobble_bulges_pipeline` → `seedmatch.wobble_bulge_matches`, `seedmatch.annotate_wobble_bulge_counts`:
    - `guide.seq` becomes a `Guide`, `features` becomes `sequences`, and `res` becomes `de` of `annotate_wobble_bulge_counts`. R ignores `res`.
    - `bulges_lim` and `nt` become `max_bulge` and `bulge_free` of `WobbleBulgeOptions`, with the R defaults 4 and 3. `silent` has no counterpart, because the port prints nothing.
    - R gives one column `SeedMatches_<j>` for each wobble seed, and `overall_matches`. The port has one pattern for each set of wobbles, thus it gives only `overall_matches`.
    - `transcript` becomes `seq_id` and `gene_id`. Name `D-005`.
  - `docs/r-map.md`, row `.get_all_combination_w_b`, `.find_factorial_internal`, `.get_target_bulged_guides`, `.sequence_with_tb`, `check_bulges_nt` → `seedmatch.wobble_bulge_patterns`. One IUPAC pattern replaces the wobble seeds. The bulged patterns insert bases, and `bulge_free` replaces `check_bulges_nt`. Name `D-005`.
  - `D-005`: "The wobble and bulge search runs, with one pattern for the wobbles and extra bases for a bulge." Write the entry from the row of `docs/design.md` and from these facts:
    - R location: `run_wobble_bulges_pipeline` (`R/run_wobble_bulges_pipeline.R:18-62`), `.get_all_combination_w_b` (`:75-129`), `.get_target_bulged_guides` (`:160-173`), `.sequence_with_tb` (`:185-198`), and `check_bulges_nt` (`:208-225`).
    - R behavior: `.get_target_bulged_guides` sets `bulge_allowance` in its own frame (`:167`). `.sequence_with_tb` reads it at `:192` and `:197`, and its environment is the namespace. Thus each call stops with `object 'bulge_allowance' not found` (note `pipeline_error`).
    - The call at `:27` gives `bulges_lim` to the parameter `guide_seq_length` and the guide length to `bulges_lim` (`:160`, note `formals_get_target_bulged_guides`).
    - `.get_all_combination_w_b` sets exactly 4 U positions to `Y` (`:82-93`). A guide with fewer than 4 U stops (`:101`, notes `wobble_seeds_2_u` and `wobble_seeds_0_u`). No seed has `R` at 5 positions (note `wobble_seeds_5_u_in_seed`).
    - `.sequence_with_tb` replaces seed letters with `N` (`:197`). Thus a "bulge" of R matches a mismatch, and no R pattern matches an extra base (the cases `mismatch_t4` and `bulge_after_g5`).
    - `check_bulges_nt` looks for an `N` in the first `nt` letters of the target (`:213`). They pair with g8, g7, and g6, not with g2 to g4 (note `check_bulges_nt`). The R documentation of `nt` names the start of the siRNA (`:11`) and the bulges (`:203`).
    - The port counts each distinct start. Two patterns with different starts and one end thus count two times, for example `ENSRNOG00000016275` at 69 and 73.
    - The counts add each bulged string of each wobble seed (`:35-55`). With the R count rule, one site gives 12 for `bulges_lim = 1` and 24 for 4 (case `perfect`).
    - New behavior: the rules of this task.
    - Reason: R gives no result. The R count rule counts one site more than one time, counts mismatches as bulges, and protects the wrong end of the seed.
  - `docs/design.md`: in the row of `run_wobble_bulges_pipeline`, replace `:137` with `:167` and `:192`. Replace `:12` and `:115` with `:27` and `:160`, and `:63` with `:82-93`. Add `:35-55` to the sum of the counts. Add the bound of `wobble_bulge_matches` to "Performance and bounds". Checkpoint.

### Task 3: The 3' compensatory search

**Files:**
- Make: `src/seedmatch/noncanonical/compensatory.py`, `tests/noncanonical/test_compensatory.py`
- Modify: `src/seedmatch/noncanonical/__init__.py`, `src/seedmatch/__init__.py`, `tests/helpers.py` (`run_core_pipeline`), `tests/test_readme.py`, `README.md`, `docs/r-map.md`, `docs/divergences.md` (`D-004`, `D-052` to `D-056`), `docs/design.md`

**Interfaces:**
- Consumes: `match_ranges`, `check_type` (Task 2), `Guide`, `Seed`, `SeedSpec`, `seed_window`, `SearchOptions`, `SequenceSet`, `reverse_complement`, `as_limit`, `check_at_least`, `check_not_negative`, `pool_width`, `map_in_pool`, and the goldens of Task 1.
- Produces:

```python
@dataclass(frozen=True)
class CompensatoryOptions:
    min_length: int = 6
    max_length: int | None = None
    first_position: int | None = None
    max_distance: int = 98
    max_seed_mismatch: int = 1
    max_region_mismatch: int = 0

def compensatory_sites(
    guide: Guide,
    sequences: SequenceSet,
    *,
    seed: SeedSpec = Seed.MER7M8,
    options: CompensatoryOptions = CompensatoryOptions(),
    workers: int | None = None,
) -> pd.DataFrame: ...
```

Rules of the options (the arguments of `get_compensatory_sequences`):
- `min_length` and `max_length` are `compensatory_region_length_min` and `compensatory_region_length_max`. `max_length=None` is the R default `NA`: each position of the region.
- `first_position` is the first guide position of a region. `None` is the guide position after the seed window: `seed_window(seed, len(guide.sequence))[1] + 1`. R always uses g8 (`R/get_compensatory_sequences.R:44`, `D-056`). Thus the default is g9 for `MER7M8` and g8 for `MER6`.
- `max_distance` is the largest number of bases between a region and its seed match. R `limit = L` gives `max_distance = L - 2`, because R keeps `seed_start < end + limit` (`:92`). Thus the default 98 is the R default 100.
- `max_seed_mismatch` is `mismatch_`, and `max_region_mismatch` is `mismatch_compensatory_region_`.
- Each field goes through `as_limit`. `min_length`, `first_position`, and `max_distance` go through `check_at_least(..., 1)`. The mismatch fields go through `check_not_negative`.
- `max_length` below `min_length` raises `max_length 3 is less than min_length 6; expected None or an integer of 6 or more` (`D-055`).

Rules of the call:
- `guide`, `sequences`, and `options` go through `check_type`. `guide.target(seed)` makes the target t and checks the seed. Let N be the length of the guide, and F the first position.
- F above N raises `the region starts at g24, but the guide ends at g23; expected a seed or a first_position that leaves a region`.
- Let A = N − F + 1, the positions of the region. `max_length` above A raises `max_length 20 is more than the 15 positions from g9 to g23; expected at most 15`. With `max_length=None`, `min_length` above A raises the same text with `min_length` (`D-055`).

Rules of the sites:
- The seed matches are the matches of t with `SearchOptions(max_mismatch=k, min_mismatch=k)`, for each k from 0 to `max_seed_mismatch`. `seed_mismatches` holds k (`D-054`). This is `vmatchPattern` with `fixed = TRUE` of R (`:78`). A match can reach past an end, and a position outside counts as a mismatch (search rule 2).
- The regions are the guide windows g_a to g_b with F ≤ a, b ≤ N, and `min_length` ≤ b − a + 1 ≤ `max_length`. The pattern of a window is `reverse_complement` of its DNA text, `region_target`. Its matches come from `SearchOptions(max_mismatch=max_region_mismatch)`. R makes the same windows (`:44-54`), but always with a = 8 or more.
- A pair is a seed match and a region match of the same sequence, by its index in the set (`D-052`). The distance is `seed_start − region_end − 1`, and it is 1 to `max_distance`.
- A pair is dropped if a pair of the same seed match has a region that contains its region and is longer (`D-053`). Pairs with the same region and different guide windows stay.
- `gc_content` is (the count of `G` and `C` in `region_target`) / its length × 100, as `.gc_content` of R (`:172-186`).
- The columns are `seq_id`, `gene_id`, `seed_start`, `seed_end`, `seed_mismatches`, `region_start`, `region_end`, `region_length`, `guide_start`, `guide_end`, `region_target`, `distance`, and `gc_content`. The text columns are `object`, `gc_content` is `float64`, and the other columns are `int64`.
- The rows come in the order of the set, then of `seed_start`, `region_start`, `region_end`, and `guide_start`, with a `RangeIndex`. A call without a site gives an empty frame with these columns and dtypes.

The algorithm keeps the cost bounded:
- The seed searches and the region searches run on `map_in_pool` with `pool_width(workers)`. The time is `max_seed_mismatch` + 1 + W searches of the set. W = Σ (A − l + 1) for l from `min_length` to `max_length`. For a guide of 23 bases with the defaults, W is 55, thus 57 searches. With `first_position=8`, W is 66.
- Drop the seed matches with a start below 1, because a region ends at 1 or later, thus a seed start of a pair is 3 or more. The key of a seed match is `offsets[index] + start - 1`. Sort the keys one time.
- For a region match in sequence i with the end e, the seed keys from `offsets[i] + e + 1` to `min(offsets[i] + e + max_distance, offsets[i + 1] - 1)` are the pairs. Two calls of `np.searchsorted` give the range. The upper bound stops the window at the end of the sequence.
- Pair each region result when `map_in_pool` gives it, thus the memory holds the seed matches, the pairs, and the results that wait in the pool.
- Sort the pairs of each seed match by `region_start`, and then by `region_end` from the highest to the lowest value. Drop a pair if an earlier pair of its seed match has another interval and a `region_end` that is not less than its `region_end`.
- The time of the pairs is O((S + R) log S + P) for S seed matches, R region matches, and P pairs. The memory is O(S + P) and the results in the pool.

The R golden cases (Task 1). The test maps the R columns with this table:

| R column | port column |
|---|---|
| `TX_ID` | `seq_id` |
| `Seed_End`, `Seed_Start` | `seed_start`, `seed_end` |
| `mismatch_in_seed` | `seed_mismatches` |
| `Compensatory_End`, `Compensatory_Start` | `region_start`, `region_end` |
| `CompensatoryLength` | `region_length` |
| `relative_position_guide_seq_start`, `relative_position_guide_seq_end` | `guide_start`, `guide_end` |
| `Compensatory_Region_in_Transcript` | `region_target` |
| `Distance` | `distance` |
| `GC_Content` | `gc_content` |

| case | port arguments |
|---|---|
| `vignette` | `CompensatoryOptions(first_position=8)` |
| `seed_mm0` | `CompensatoryOptions(first_position=8, max_seed_mismatch=0)` |
| `min8_max10` | `CompensatoryOptions(first_position=8, min_length=8, max_length=10)` |
| `limit40` | `CompensatoryOptions(first_position=8, max_distance=38)` |
| `region_mm1` | `CompensatoryOptions(first_position=8, min_length=9, max_region_mismatch=1)` |
| `mer6` | `seed=Seed.MER6`, and the default options |

The expected port rows of the synthetic cases, with `first_position=8`. A row is (`seed_start`, `seed_end`, `seed_mismatches`, `region_start`, `region_end`, `guide_start`, `guide_end`, `region_target`, `distance`). The planner computed them with a Python version of these rules:

| case | rows |
|---|---|
| `seed_past_end` | (23, 29, 1, 5, 12, 8, 15, `GTTCTTGC`, 10). This is the R row. |
| `no_site` | none |
| `repeated_name_r` | sequence 0: (23, 29, 0, 5, 12, 8, 15, `GTTCTTGC`, 10). Sequence 1: (58, 64, 0, 30, 37, 8, 15, `GTTCTTGC`, 20). |
| `two_seed_matches_r` | (13, 19, 0, 6, 11, 10, 15, `GTTCTT`, 1) and (30, 36, 0, 6, 13, 8, 15, `GTTCTTGC`, 16) |
| `repeated_region_r` | (24, 30, 0, 7, 12, 8, 13, `GCTTGC`, 11) and (24, 30, 0, 7, 12, 12, 17, `GCTTGC`, 11) |
| `region_at_g2_r` | (24, 30, 0, 7, 12, 8, 13, `TCTATA`, 11) |
| `seed_mm2_r` (`max_seed_mismatch=2`) | (23, 29, 2, 5, 12, 8, 15, `GTTCTTGC`, 10) |
| `min_above_max_r` | `InputError` from `CompensatoryOptions(min_length=6, max_length=3)` |

- [ ] **Step 1: Write the failing tests** in `tests/noncanonical/test_compensatory.py`. A helper reads the cases of `tests/data/compensatory_cases.tsv` into a `SequenceSet` and `CompensatoryOptions(max_distance=limit - 2, ...)`.
  - `test_compensatory_sites_match_r_D004(case)`: parametrized over the cases of the test data. Map the R rows with the table above, and sort them by the index of `seq_id` in `sirna_seqs()`, `seed_start`, `region_start`, `region_end`, and `guide_start`. The port rows equal the R rows. `gc_content` agrees within a relative tolerance of 1e-9, and the other values are equal. The docstring gives the notes `original_vignette` and `line_99`.
  - `test_the_rows_come_in_the_order_of_the_set`: the vignette rows come in the order of the index of the sequence, `seed_start`, `region_start`, `region_end`, and `guide_start`.
  - `test_a_set_without_a_site_gives_an_empty_frame`: `no_site` gives the columns and the dtypes of "Rules of the sites", and no row. The docstring gives its note: the patched R stops, and the original R returns 0 rows.
  - `test_a_seed_match_past_the_end_counts_the_outside_position`: `seed_past_end` gives the R row of the golden.
  - `test_each_sequence_is_searched_for_a_repeated_id_D052`: `repeated_name_r` gives the rows of the table, in the order of the set. The golden has 4 rows, each of the first sequence.
  - `test_a_region_of_one_seed_match_does_not_drop_a_region_of_another_D053`: `two_seed_matches_r` gives the rows of the table. The golden has one row, at `Seed_End` 30.
  - `test_a_region_inside_a_longer_region_of_its_seed_match_is_dropped`: with the default options, `ENSRNOG00000016275` has one row: (73, 79, 0, 59, 71, 10, 22, `AAACAGTGTTCTT`, 1), with `region_length` 13.
  - `test_each_row_gives_the_guide_positions_of_its_region_D054`: `repeated_region_r` and `region_at_g2_r` give their rows. The golden has the text `8,12` and `13,17` in the two equal rows of `repeated_region_r`, and `2,8` and `7,13` for `region_at_g2_r`.
  - `test_seed_mismatches_is_the_number_of_mismatches_D054`: `seed_mm2_r` gives `seed_mismatches` 2. The golden gives `mismatch_in_seed` 1.
  - `test_the_region_starts_after_the_seed_D056`: with the default options, the vignette arguments give 96 rows in 51 genes, and no `guide_start` is below 9. With `first_position=8`, they give the 104 rows of the golden. `repeated_region_r` gives only the row with `guide_start` 12.
  - `test_compensatory_options_reject_lengths_that_select_no_region_D055`: `CompensatoryOptions(min_length=6, max_length=3)`, `min_length=17` with `GUIDE`, `max_length=20` with `GUIDE`, and `max_distance=0` raise `InputError` with a `match` on the template. The docstring gives the 3 rows of `min_above_max_r`, the note `min_17`, and the R `limit = 2`, which gives no row.
  - `test_a_short_guide_has_no_room_for_the_region_D055`: `Guide("UUAUAGAGCAAGA")` with the defaults raises `min_length 6 is more than the 5 positions from g9 to g13; expected at most 5`.
  - `test_a_seed_that_ends_at_the_guide_end_leaves_no_region`: `seed=Seed.FULL` raises `the region starts at g24, but the guide ends at g23`.
  - `test_compensatory_options_reject_a_first_position_below_1` and `test_compensatory_sites_reject_a_first_position_after_the_guide`: `first_position=0` and `first_position=24`.
  - `test_compensatory_sites_equal_an_oracle(guide, records, options)`: a hypothesis test with `deadline=None`. The guide has 10 to 16 letters of `ACGU`. The records have 1 to 3 sequences of 0 to 50 letters of `ACGTN`, and they hold copies of the target and of `region_target` texts. `min_length` is 2 to 4, `max_length` is `None` or `min_length` to `min_length` + 3, `max_distance` is 1 to 20, `max_seed_mismatch` is 0 to 2, and `max_region_mismatch` is 0 to 1.
  - In the oracle test, the oracle tries each start from 2 − m to L of each pattern of length m. It counts a position outside the sequence as a mismatch, and a letter that is not equal as a mismatch. It then pairs, filters, and drops the pairs with the rules above, one pair at a time.
  - `test_a_window_stops_at_the_end_of_its_sequence`: the set `("a", "GGGGGTTCTTGC")`, `("b", "CTCTATAGGGG")` gives no row. A window from the region of `a` that runs into `b` would pair it with the seed match of `b`.
  - `test_an_n_in_the_seed_site_is_a_mismatch`: `GGGGGTTCTTGCGGGGGGGGGGCTCNATAGGGGG` gives one row with `seed_mismatches` 1.
  - `test_compensatory_sites_do_not_depend_on_workers`: `workers=1` and `workers=4` give equal frames.
  - `test_compensatory_sites_reject_a_guide_of_another_type`, `test_compensatory_sites_reject_sequences_of_another_type`, `test_compensatory_sites_reject_options_of_another_type`, and `test_compensatory_options_reject_a_negative_mismatch`: `InputError` with a `match` on the template.
- [ ] **Step 2: Run the tests.** Run `uv run pytest tests/noncanonical/test_compensatory.py -q`. Expected: FAIL with `ImportError` for `compensatory_sites`.
- [ ] **Step 3: Write `compensatory.py`.** Export `CompensatoryOptions` and `compensatory_sites` from `seedmatch.noncanonical` and from `seedmatch`. Add a `compensatory_sites` call to `run_core_pipeline`.
- [ ] **Step 4: Run the tests.** Expected: PASS, and ruff and mypy pass.
- [ ] **Step 5: Measure.** Do Review Focus 1 for `compensatory_sites`. Give the numbers in the report of the task and in "Performance and bounds".
- [ ] **Step 6: Write the README part and the documents.**
  - `README.md`, "Other analyses": a paragraph and a Python block after the wobble block. The block calls `sites = compensatory_sites(guide, sequences)` and prints the head. Then it takes `site_genes = sites["gene_id"].unique()` and `other_genes`, the genes of `de` without a site. It calls `ecdf_test(de, target=site_genes, background=other_genes)`, as the R vignette does. The paragraph gives the window rule, `first_position`, `max_distance`, and the R `limit`.
  - `tests/test_readme.py`: add `test_readme_compensatory_example_gives_the_sites`. `sites` has 96 rows in 51 genes, and the p-value of the ECDF test is from 0 to 1.
  - `README.md`, "Scope": add the 3' compensatory search to the first list, and remove it from the second list. Add a bullet for `D-004` to "Differences".
  - `docs/r-map.md`, row `get_compensatory_sequences` → `seedmatch.compensatory_sites`:
    - `guide.seq` becomes a `Guide`, and `sequences` keeps its name. `seed` becomes a `SeedSpec`, and the default is `Seed.MER7M8`, the seed of the R vignette. R takes any target text.
    - The other arguments become the fields of `CompensatoryOptions`. `limit = L` becomes `max_distance = L - 2`.
    - The R columns map to the port columns as in the table of this task. R gives the transcript coordinates of the seed and of the region in swapped order: `Seed_Start` is the end of the seed match.
    - `Seed_in_Transcript` and `Seed_Region_in_Guide` have no column, because each row holds `guide.target(seed).sequence` and its reverse complement. `Compensatory_Region_in_Guide` has no column: it is the DNA text of `guide.sequence[guide_start - 1 : guide_end]`.
    - The gene list of the R vignette, `unique(comp_matches$TX_ID)`, becomes `sites["gene_id"].unique()`. Name `D-004` and `D-052` to `D-056`.
  - `docs/r-map.md`, row `relative_position_seq`, `.gc_content` → the columns `guide_start`, `guide_end`, and `gc_content` of `compensatory_sites` (`D-054`).
  - `docs/r-map.md`, row `.update_sequence_at` → no counterpart. No SeedMatchR function calls it.
  - `D-004`: "compensatory_sites runs the search and returns the sites."
    - R location: `get_compensatory_sequences` (`R/get_compensatory_sequences.R:57-77`, `:34`, `:99`).
    - R behavior: the frame of `:57-66` has no row, thus the test of `:67` is true. The function returns 0 rows with the message `No compensatory regions found.` for each input (note `original_vignette`).
    - Without the early return, `:99` stops at the first site. `gsub` at `:34` makes the seed a character string, and `reverseComplement` has no method for it (note `line_99`).
    - New behavior: the search runs. Each row is a pair of a seed match and a region. The golden of the patched function gives the R rows of each case of the test data.
    - Reason: R gives no result for each input.
  - `D-052`: "compensatory_sites searches each sequence, also for a repeated ID."
    - R location: `:36-41`, `:90`, and `:106-108`.
    - R behavior: R selects the sequences with a seed match by name (`sequences[names]`). For a name that occurs two times, this gives the first sequence two times (note `repeated_name_r`). The merge by `gene_name` then pairs each seed row with each region row of the name. The case `repeated_name_r` gives 4 rows of the first sequence, and no row of the second.
    - New behavior: each sequence of the set is searched, and a pair holds a seed match and a region of the same sequence.
    - Reason: a `SequenceSet` can hold two sequences with one ID, for example with `PerGene.EACH_TRANSCRIPT` and gene IDs as names.
  - `D-053`: "A region is dropped only for a longer region of the same seed match."
    - R location: `:123-125`.
    - R behavior: the `dplyr` filters keep the longest region of each transcript and region end, and then of each transcript and region start. The groups ignore the seed match. Thus a longer region of one seed match drops the region of another seed match. The case `two_seed_matches_r` gives one row, and the seed match at 13 has no row.
    - New behavior: the rule of this task.
    - Reason: each seed match is a site of its own. On the R goldens without this case, the R rule and the port rule give the same rows.
  - `D-054`: "Each row gives the guide positions of its region and the mismatches of its seed match."
    - R location: `relative_position_seq` (`:3-8`), `:127-128`, and `:106-108`.
    - R behavior: `relative_position_seq` searches the full guide and joins each start with a comma. Thus a region that occurs two times in the guide gives two equal rows with `8,12`. A region at g8 that also occurs at g2 gives `2,8`.
    - `mismatch_in_seed` is 1 if the text of the seed match differs from the target. Thus a match with 2 mismatches gives 1 (case `seed_mm2_r`).
    - New behavior: `guide_start` and `guide_end` are the integers of the window of the row. `seed_mismatches` is the number of mismatches, with each position outside the sequence. For `max_seed_mismatch` of 0 or 1, it equals the R value.
    - Reason: a text of positions is not a number, and the position of g2 is not in the region. The flag hides the number of mismatches.
  - `D-055`: "CompensatoryOptions rejects lengths and distances that select no region."
    - R location: `:46-54` and `:92`.
    - R behavior: `compensatory_region_length_min = 6` with `compensatory_region_length_max = 3` gives the lengths 6, 5, 4, and 3, and R returns 3 rows of length 6 (case `min_above_max_r`). A minimum above the 16 positions searches g8 to g23 three times (note `min_17`). `limit = 2` gives no row.
    - New behavior: `InputError` for `max_length` below `min_length`, a length above the positions of the region, and a `max_distance` below 1.
    - Reason: R gives a result for other lengths than the arguments, and no message.
  - `D-056`: "The compensatory region starts after the seed."
    - R location: `:44`.
    - R behavior: the region is always g8 to gN. With the `mer7m8` seed, g8 pairs with the seed match and with the region of one row. In the golden `vignette`, 10 of the 104 rows have a region from g8.
    - New behavior: the default region starts after the seed window. `first_position=8` gives the R region.
    - Reason: one guide base pairs with one target base.
  - `docs/design.md`: in the row of `get_compensatory_sequences`, replace `R/get_compensatory_sequences.R:51-60` with `:57-77`, and add the stop at `:99`. Add the table rows of `D-052` to `D-056`. Add the bound of `compensatory_sites` to "Performance and bounds". Checkpoint.

### Task 4: Verification of the plan

- [ ] **Step 1: Run the checks one time.** The controller runs `uv run ruff check`, `uv run ruff format --check`, `uv run mypy`, and `PYTHONPATH=tmp/claude/covfix NUMBA_DISABLE_JIT=1 uv run pytest -p cov_serial --cov`. The plugin `cov_serial` runs the pool of `seedmatch._pool` in the caller thread, except in `tests/test_pool.py`. Without it, the coverage tracer stops the pool threads for hours. Expected: each passes, and the coverage is not below 95%. Give a reason for each line of the new modules without coverage.
- [ ] **Step 2: Compare the documents with the code.**
  - Each public name of `seedmatch.noncanonical` has a row in `docs/r-map.md`, or is a type of a function that has a row.
  - Each `D-NNN` of this plan has an entry, a table row, and a test.
  - Each README example runs in `tests/test_readme.py`.
  - `tests/test_package.py::test_import_loads_no_fetch_and_no_plot_modules` passes.
  - Checkpoint.
