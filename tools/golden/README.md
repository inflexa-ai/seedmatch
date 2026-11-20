# Golden files from R

The scripts in this directory run SeedMatchR and write the golden files in `tests/golden/`. The Python tests read the golden files and compare them with the results of seedmatch.

The tests never run R. Run the commands below only to make or to update a golden file.

## Commands

Run each command from the repository root.

Install the R packages into `tools/golden/.rlib`. The script stops if `../SeedMatchR` is not at commit `2fc0c4f` or if it has a change to a tracked file:

```bash
Rscript tools/golden/install.R
```

Make sure that the environment is correct. The last line of the output is `golden environment OK` and the versions:

```bash
tools/golden/run.sh tools/golden/R/check_env.R
```

To make golden files, replace `check_env.R` with the name of another script in `tools/golden/R/`.

The search corpora have a Python step. `make_corpus.py` writes the cases of each corpus with a fixed seed. Then `search_corpus.R` writes the Biostrings result for each case:

```bash
uv run python tools/golden/make_corpus.py
tools/golden/run.sh tools/golden/R/search_corpus.R
```

The script `full_search.R` makes the golden files of `full_search`, `full_search_ranges`, and the `data.frame` mode of `SeedMatchR`:

```bash
tools/golden/run.sh tools/golden/R/full_search.R
```

The script `report.R` makes the golden files of `SeedMatchReport`. The report uses the counts of `full_search`, thus the script does not source the patch of `vmatchPattern2`:

```bash
tools/golden/run.sh tools/golden/R/report.R
```

The annotation golden files have a Python step. `tools/fixtures/annotation.py` writes the synthetic fixture in `tests/data/annotation/` from one specification. The script runs `bgzip` of htslib and `samtools`, thus they must be on `PATH`. Then `annotation.R` reads the fixture and writes the golden files in `tests/golden/annotation/`:

```bash
uv run python tools/fixtures/annotation.py
tools/golden/run.sh tools/golden/R/annotation.R
```

If the versions of numpy, SQLite, zlib, htslib, and samtools do not change, the fixture script writes the same bytes each time. `annotation.R` replaces `load_species_anno_db` in the SeedMatchR namespace. The replacement gives the local EnsDb and 2bit file, not the AnnotationHub files. Thus `load_annotations` runs unchanged.

For each config of `ANNOTATION_CONFIGS`, `annotation.R` writes the feature ranges, the sequences, and the transcripts of the filter. The script also writes `genome_ranges.tsv`. For each range, this file holds the sequence from the 2bit reader of rtracklayer and from the BGZF reader of Rsamtools. rtracklayer 1.70.1 cannot read the 2bit file of version 1, thus the column `twobit_v1` holds `unsupported`.

The Python tests do not compare the configs whose names end with `_r`. These golden files show the R behavior that a divergence changes.

## Golden files of rnor7

The script `annotation_rnor7.R` writes the golden files of the rat reference `rnor7`, Ensembl release 109 of mRatBN7.2. The script reads the directory `rnor7` in the directory of the variable `SEEDMATCH_TEST_REFS`. If the variable is not set, the script stops:

```bash
SEEDMATCH_TEST_REFS=~/Development/refs/seedmatch tools/golden/run.sh tools/golden/R/annotation_rnor7.R
```

The directory of the variable has this layout:

```text
$SEEDMATCH_TEST_REFS/
└── rnor7/
    ├── EnsDb.Rnorvegicus.v109.sqlite
    ├── Rattus_norvegicus.mRatBN7.2.109.gtf.gz
    └── Rattus_norvegicus.mRatBN7.2.dna_sm.toplevel.2bit
```

- The EnsDb file is the AnnotationHub resource `AH109732`, and the 2bit file is the AnnotationHub resource `AH106786`. SeedMatchR uses these resources for `rnor7`.
- The GTF file is `https://ftp.ensembl.org/pub/release-109/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.109.gtf.gz`. The script does not read it, but the Python tests compare the GTF with the EnsDb.
- `tests/annotation/test_refs.py` gives the size and the SHA-256 digest of each file.

The script uses the configs `utr3`, `utr5`, `exons`, `cds`, `utr3_longest`, and `utr3_union` of `ANNOTATION_CONFIGS`. For each config, it writes `tests/golden/annotation/rnor7/<config>.tsv.gz` with the columns `name_tx`, `name_gene`, `width`, and `digest`. The sequences of the features of `rnor7` are too large for a golden file. Thus the column `digest` holds the first 16 hexadecimal digits of the SHA-256 digest of each sequence, from `openssl::sha256`.

The script uses the R package `openssl`. SeedMatchR does not import it directly. If the package is not installed, the script stops at the start and names the package.

## Golden files of the statistics

The script `stats.R` writes the golden files in `tests/golden/stats/`. They give the R results of the KS test, the Wilcoxon rank-sum test, DTS, the Wasserstein areas, the `round` function, and `ecdf_stat_test`:

```bash
tools/golden/run.sh tools/golden/R/stats.R
```

`twosamples` and `caTools` are imports of SeedMatchR, thus `install.R` makes them available. The comment lines of each file give their versions.

The script makes the samples with `set.seed(20261120)`. One random stream feeds each case in the order of the script. A new case that draws moves each later case, thus a new run changes the goldens of the later sections. Only a case without a random draw keeps the other cases. A case at the end of the Wass section moves no other case. The random values of `round.tsv` have their own `set.seed(20261120)`, thus no case moves them. The second run of the script writes the same bytes as the first run. The script stops if a value of the R snapshots of `ecdf_stat_test` changes (`D^+` 0.51035, `W` 2150, and the Wass statistic 0).

The ID of a case is `<test>_<m>x<n>_<kind>`. The test is `ks`, `wilcoxon`, `dts`, or `wass`. The sizes of the samples `x` and `y` are `m` and `n`. The kind is one of these:

- `plain`: no value occurs two times.
- `ties`: the values of the `plain` case, rounded to 1 decimal. At least one value occurs in both samples.
- `shift`: the sample `y` is 6 standard deviations above the sample `x`, thus D is near 1. Only KS has this kind.

The sizes of the cases:

- KS: (1, 12), (5, 7), (30, 40), (99, 101), (100, 100), (120, 150), and (400, 300). The sizes 99 × 101 give the exact p-value, and the sizes 100 × 100 give the asymptotic p-value. The kind `shift` exists for the sizes 99 × 101 and larger.
- Wilcoxon: (1, 5), (10, 12), (49, 49), (49, 50), (50, 49), and (60, 80). Only a case without ties and with both sizes below 50 gives the exact p-value.
- DTS: (5, 7), (30, 40), and (200, 300).
- Wass: (10, 10), (10, 50), and (40, 300), each with `n.interp` 2, 100, and 10000.

Some cases have a name that tells their purpose:

- `ks_all_equal` and `wilcoxon_all_equal`: each value is equal.
- `ks_inf` and `wilcoxon_inf`: the sample `x` has the value `Inf`.
- `ks_asym_cut_series`: a two-sided asymptotic case with `sqrt(en) * D` between 0.9 and 1. R uses the cut series here.
- `ks_asym_small_p`: a one-sided asymptotic case with a p-value below 1e-17. R gives 0.
- `wass_seed_one_value` and `wass_non_seed_one_value`: R stops, because a sample has one distinct value. These cases have no row in `wasserstein.tsv`. The comment lines of the file give the message of R.

The files:

- `samples.tsv`: the columns `case`, `sample` (`x` or `y`), and `value`. Each value has one row, in the order of the sample. The text `Inf` is an infinite value, and Python `float()` reads it. Each other file refers to the `case` of this file.
- `ks.tsv` and `wilcoxon.tsv`: the columns `case`, `alternative` (`two.sided`, `less`, or `greater`), `statistic`, `p_value`, and `method`. `method` is `exact` or `asymptotic`. The script maps the `method` text of R: `Exact` or `Asymptotic` for KS, and `exact` or `continuity correction` for Wilcoxon.
- `dts.tsv`: the columns `case`, `statistic`, `p_value`, and `n_boot` (2000). The script calls `set.seed(1)` before each `dts_test`.
- `wasserstein.tsv`: the columns `case`, `n_interp`, `abs_auc`, `neg_cens_auc`, and `pos_cens_auc`. The sample `x` is the seed sample, and the sample `y` is the non-seed sample. R rounds each area to 3 decimals.
- `round.tsv`: the columns `x`, `digits`, and `r`. `r` is `round(x, digits)`. The values of `x` are the boundary values of the rule, the grid `(0:4000)/4000 + 0.0005`, and 200 random values.
- `ecdf.tsv`: the columns `input`, `test`, `alternative`, `statistic`, `p_value`, `n_target`, and `n_background`.

Notes on `ks.tsv` and `wilcoxon.tsv`:

- For KS, `greater` gives `D^+` and `less` gives `D^-`.
- R finds D with a cumulative sum, thus D can differ from the exact value by about 1e-15.
- For Wilcoxon, `NA` in `p_value` means that R gives NaN, because each value is equal.

The input of `ecdf.tsv` is one of these:

- `test_data`: the table `tests/data/DESEQ2_results.tsv`. The target is the genes with `mer7m8 >= 1`, and the background is the genes with `mer7m8 == 0`.
- `test_data_swapped`: the same table, with the target and the background in the other order. It has only the rows of the R snapshots: Wilcoxon `less` and Wass.
- `readme`: the table that the first example of the README ports. It comes from `filter_res` with the defaults and the `mer7m8` counts of `SeedMatchR`. The target and the background are the same sets.

The columns of `ecdf.tsv` obey these rules:

- `test` is `KS`, `Wilcoxon`, `Wass`, or `DTS`. R spells the argument `Wilcoxen`.
- `Wass` and `DTS` ignore the alternative in R. Their rows have `NA` in `alternative`.
- The `statistic` of `Wass` is `neg_cens_auc`, and its `p_value` is `NA`.
- The `p_value` of `DTS` comes from `set.seed(1)` before the call.
- `n_target` and `n_background` count the fold changes that are not `NA`.

## Golden files of the small analyses

The script `analyses.R` writes the golden files in `tests/golden/analyses/`. They give the R results of these functions and files:

- the k-mer counts of `sequence_kmer_counts` and `kmer_hist`
- the alignment text of `plotAlignment`
- the distances of `match_distances`
- the e-value of `generate_e_value`
- the PR and ROC curves of PRROC
- `benchmark_preds`
- the miRDB table of `get_example_data`

Run the script from the repository root:

```bash
tools/golden/run.sh tools/golden/R/analyses.R
```

The script reads `Test_seqs.fa` and `DESEQ2_results.tsv` from `tests/data/`. It also reads the synthetic fixtures `kmer_cases.tsv` and `mirdb_sample.txt` from the same directory.

The second run of the script writes the same bytes as the first run. A section that draws random numbers starts with its own `set.seed`. The e-value section uses 20261122, and the curve section uses 20261123. Thus a new case in one section moves no case of another section. A new random case in the curve section moves the random cases after it.

### PRROC

`benchmark_preds` calls PRROC, which is a Suggests of SeedMatchR. The golden files use PRROC 1.4. `install.R` installs the CRAN source file `PRROC_1.4.tar.gz` into `tools/golden/.rlib`. It stops unless the MD5 digest of the file is `229140aa270cdad059daf08783b083bb`. `check_env.R` stops unless the PRROC version is exactly 1.4. The comment lines of each file with a PRROC result give the version.

### Patches

`analyses.R` puts these patches into the SeedMatchR namespace with `with_seedmatchr_patch`. A patch is a copy of one SeedMatchR function. Its first comment gives the divergence, the R line, and the fix:

- `tools/golden/patches/plotAlignment.R` (`D-045`): `plotAlignment` takes the guide target from `get_seed(guide.seq, "Full", ...)` (`R/plotAlignment.R:32`). This target starts at g2, thus the guide string is one base shorter than the window, and R stops. The patch takes the target from `.get_custom_sequence(guide.seq, 1, guide.len)`. As a result, the guide string has the full guide, g1 to gN.
- `tools/golden/patches/benchmark_preds.R` (`D-006`): the filter `padj <= padj` of `benchmark_preds` (`R/benchmark.R:23`) compares the column `padj` with itself, thus the argument `padj` has no effect. The patch compares the column with the argument. It also gives `roc.curve` the labels as `weights.class0`, as `pr.curve` gets them.

The comment lines of the golden files give the results of the unpatched functions where they differ.

### K-mers

- `kmer_counts.tsv`: the columns `case`, `seq_id`, `gene_id`, `kmer`, and `count`. It has one row for each count above 0 of `sequence_kmer_counts`, in the order of its rows and then of its columns.
- `kmer_frame.tsv`: the data frame of `sequence_kmer_counts` for `cases_k2_s1`, as R gives it. Its columns are `tx_id`, `gene_id`, and one column for each 2-mer.
- `kmer_gene_counts.tsv`: the columns `case`, `kmer`, and `n`. The value `n` is the number of sequences that have the k-mer. It is `colSums` of the presence matrix of `R/kmer_hist.R:19-28`. The file has one row for each k-mer, in column order.

The ID of a case is `<set>_k<width>_s<step>`. The set is one of these:

- `test_data_first10`: the first 10 sequences of `Test_seqs.fa`. With the width 7 and the step 1, this is the call of the k-mer vignette.
- `test_data`: each sequence of `Test_seqs.fa`.
- `cases`: the sequences of `kmer_cases.tsv`.

For `Test_seqs.fa`, `tx2gene` maps each name to itself. For `kmer_cases.tsv`, it maps `seq_id` to `gene_id`.

R counts a window only if each base is A, C, G, or T. Biostrings reads a lowercase base as an uppercase base. A sequence that is shorter than the width has no row in `kmer_counts.tsv`. The empty sequence `s4` has no row either.

### Alignment

`alignment.tsv` has the columns `case`, `seq_id`, `start`, `end`, `seed`, `start_pos`, `stop_pos`, `left_pad`, `right_pad`, `line`, and `text`.

- `seq_id`, `start`, and `end` give the match. `seed` holds the R seed name.
- `start_pos` and `stop_pos` are the arguments `seq.start.pos` and `seq.stop.pos`.
- A case has 3 rows, for `line` 1 to 3: the siRNA, the match marks, and the target.
- `text` is a value of the column `alignment` of the patched `plotAlignment`. Each text starts with spaces, and the three texts of a case have the same length.

R takes the window of the sequence from `start - (N - stop_pos) - left_pad` to `end + (start_pos - 1) + right_pad`. N is the length of the guide. For `mer7m8` and a guide of 23 bases, the window starts 15 bases before the match.

The cases use `GUIDE` of `common.R` and the sequences of `Test_seqs.fa`. In a case ID, `i` is the 1-based row of the match:

- `mer7m8_<i>`: each match of `vmatchPattern2` for the `mer7m8` target. The positions are 2 and 8, and the pads are 0.
- `mer7m8_pad_<i>`: the first 3 of these matches, with `left.pad` 3 and `right.pad` 2.
- `mer8_<i>`: each match of the `mer8` target. The positions are 1 and 8.
- `mer6_<i>`: the first 5 matches of the `mer6` target. The positions are 2 and 7.
- `full_mm2_<i>`: each match of `Biostrings::vmatchPattern` for the `Full` target with `max.mismatch` 2, if its window is in its sequence. The positions are 2 and 23.
- `custom_<i>`: each match of the custom seed from position 3 to 12. The `seed` column holds `custom`.
- `edge_first_1` and `edge_last_1`: a `mer7m8` window that starts at position 1 of the first sequence, and a window that ends at its last base. These positions are not search results.

The comment lines give the R error for these inputs:

- the unpatched `plotAlignment`
- a row of `SeedMatchR(res.format = "data.frame")`
- a window that starts before 1
- a window that ends after its sequence
- a match with an indel

### Distances

- `distances_ranges.tsv`: the columns `case`, `seq_id`, `start`, `end`, and `width`. These are the ranges that the script gives to `match_distances`, in list order.
- `distances.tsv`: the columns `case`, `anchor`, `seq_id`, and `distance`. The `anchor` is `fivePrimeEnd`, `threePrimeEnd`, or `center`. The file has one row for each value of the R list, in list order and in the order of each vector. `NA` stays `NA`.

Each case splits its ranges by name with `factor(names, levels = unique(names))`. Thus the list keeps the order of the first range of each name. A list element with one range gives the values `NA` and 0. An element without a range gives `NA` and `NA`.

The cases:

- `mer6`: the matches of `vmatchPattern2` for the `mer6` target. Some names have one match, and some names have more.
- `mer7m8_indel`: the matches of `vmatchPattern2` for the `mer7m8` target with `max.mismatch` 1 and `with.indels` `TRUE`. Their widths differ, thus `center` has half values.
- `synthetic`: the names `a` (starts 10, 30, 35, widths 7, 7, 8), `b` (start 5, width 7), and `c` (starts 100, 50, width 7).

### E-value

The genes are `g01` to `g40`. The files:

- `evalue_observed.tsv`: the columns `gene_id` and `overall_matches`. The gene `g01` has the count 0.
- `evalue_null.tsv`: the columns `gene_id`, `shuffle` (1 to 25), and `count`. The rows follow the null table of R. This table has the genes in reverse order, thus a match by row position fails.
- `evalue.tsv`: the columns `gene_id`, `overall_matches`, and `e_value`. The rows are the result of `generate_e_value`, in its row order. The e-value is the number of shuffles with a count at or above the observed count, divided by the number of shuffles.

The comment lines of `evalue.tsv` give the R results for these inputs:

- no null table
- an observed count `NA` or -1
- observed rows with the same gene
- a null with two rows for a gene
- a gene that is not in the null
- observed rows where only some genes are in the null

### Curves

The R calls are `PRROC::pr.curve(scores.class0 = score, weights.class0 = label, curve = TRUE)` and `PRROC::roc.curve(scores.class0 = score, weights.class0 = label, curve = TRUE)`. The files:

- `curve_samples.tsv`: the columns `case`, `score`, and `label` (0 or 1).
- `pr.tsv`: the columns `case`, `auc_integral`, and `auc_davis_goadrich`.
- `pr_curve.tsv`: the columns `case`, `row`, `recall`, `precision`, and `threshold`.
- `roc.tsv`: the columns `case` and `auc`.
- `roc_curve.tsv`: the columns `case`, `row`, `fpr`, `tpr`, and `threshold`.

The value `row` is the 1-based row of the matrix `curve` of PRROC. The benchmark cases also have rows in `curve_samples.tsv`, `pr_curve.tsv`, and `roc_curve.tsv`. Their areas are in `benchmark.tsv`, not in `pr.tsv` and `roc.tsv`.

The curve cases start with `set.seed(20261123)`. Each case has at least one positive and one negative:

- `ties_small`: 40 scores from `sample(0:3, 40, replace = TRUE)`, and labels from `rbinom(40, 1, 0.4)`.
- `continuous`: 60 scores from `rnorm(60)`, and labels from `rbinom(60, 1, 0.5)`.
- `large`: 400 scores from `sample(0:5, 400, replace = TRUE)`, and the label 1 with the probability `0.3 + 0.1 * score`. It has at least 100 positives, thus its curve uses the recall step 1/P.
- `one_value`: the score 1 for each of 10 points, and the labels `rep(c(0, 1), 5)`.
- `perfect`: the scores `5:1`, and the labels `c(1, 1, 0, 0, 0)`.
- `first_negative`: the scores `c(9, 9, 5, 5, 1, 1)`, and the labels `c(0, 0, 1, 0, 1, 1)`. The top score has only negatives.
- `inverse`: the scores `1:6`, and the labels `c(1, 1, 1, 0, 0, 0)`.

### Benchmark

`benchmark.tsv` has the columns `case`, `column`, `padj`, `l2fc`, `n_positive`, `n_negative`, `auc_integral`, `auc_davis_goadrich`, and `roc_auc`. The script runs the patched `benchmark_preds` on `load_test_data()$res` with `preds_col = "mer7m8"`. The cases:

- `benchmark_1`: `padj` 0.15 and `l2fc` 0.
- `benchmark_2`: `padj` 1e-10 and `l2fc` 0.
- `benchmark_3`: `padj` 1e-10 and `l2fc` -1.
- `benchmark_na`: the filter of `benchmark_2`, on a copy of `res` where `padj` is `NA` for two genes that pass the filter. R makes these rows negatives.

`n_positive` and `n_negative` count the labels of the filter of the patch. The comment lines give the result of the unpatched `benchmark_preds(res, padj = 1e-10)`. They also give the result of the patched function for a `NaN` score, for `padj = 0` (no positive), and for `padj = 1` with `l2fc = Inf` (no negative).

### miRDB

`mirdb.tsv` is the table of `R/get_data.R:43-52` for `tests/data/mirdb_sample.txt`. Its columns are `miRDB.ID`, `target.REFSEQ.ID`, `target.score`, `species`, and `miRNA.name`.

R splits the ID at the first hyphen. The ID `rno-let-7a-2-3p` has the species `rno` and the name `let-7a-2-3p`. The ID `lin4` has no hyphen, thus its species is `lin4` and its `miRNA.name` is empty.

The comment lines give the R result for these lines of a file:

- a text score
- an empty score
- the ID `NA`
- a line with 2 fields
- a field with `#`
- a field with `'`

## Golden files of the non-canonical searches

The script `noncanonical.R` writes the golden files in `tests/golden/noncanonical/`. They give the R results of `get_compensatory_sequences` and of the wobble and bulge pipeline:

```bash
tools/golden/run.sh tools/golden/R/noncanonical.R
```

The script reads `Test_seqs.fa` and `DESEQ2_results.tsv` from `tests/data/`. It also reads the synthetic fixtures `compensatory_cases.tsv` and `wobble_bulge_cases.tsv` from the same directory. The script draws no random numbers, thus the second run writes the same bytes as the first run.

The script stops if the R result of a note changes from the value that this section gives. The notes `wobble_seeds_*` and `formals_get_target_bulged_guides` have no such check.

### Patch of get_compensatory_sequences

`get_compensatory_sequences` of SeedMatchR always returns an empty data frame (`D-004`). The file `tools/golden/patches/get_compensatory_sequences.R` is a copy of the function with these changes. `noncanonical.R` puts it into the SeedMatchR namespace with `with_seedmatchr_patch`.

- The patch removes the test `if (nrow(compensatory_regions_captured) == 0)`, its early return, the `else {`, and the brace that closes the `else` (`R/get_compensatory_sequences.R:67-77` and `:130`). The frame of `:57-66` is always empty, thus the test is always true. The patch keeps the body of the `else` (`:78-129`).
- At `:99`, the patch writes `Biostrings::reverseComplement(Biostrings::DNAString(seed))`. The `gsub` at `:34` makes `seed` a character string, and `reverseComplement` has no method for it. Without this change, the function stops at the first match.

### Compensatory search

`compensatory.tsv` has the column `case`, then each R column of the result in R order, from `Seed_in_Transcript` to `relative_position_guide_seq_end` (`R/get_compensatory_sequences.R:116-128`). In each case, the seed is `get_seed(guide, <seed>)$Target.Seq`.

The test-data cases use `GUIDE` and the sequences of `Test_seqs.fa`. An argument that a case does not name has the R default: `compensatory_region_length_min = 6`, `compensatory_region_length_max = NA`, `limit = 100`, `mismatch_ = 1`, and `mismatch_compensatory_region_ = 0`.

- `vignette`: the seed `mer7m8` and the defaults.
- `seed_mm0`: `mer7m8` and `mismatch_ = 0`.
- `min8_max10`: `mer7m8`, `compensatory_region_length_min = 8`, and `compensatory_region_length_max = 10`.
- `limit40`: `mer7m8` and `limit = 40`.
- `region_mm1`: `mer7m8`, `compensatory_region_length_min = 9`, and `mismatch_compensatory_region_ = 1`.
- `mer6`: the seed `mer6` and the defaults.

The synthetic cases come from `tests/data/compensatory_cases.tsv`. Each case ID of this file is a case of the golden file. The seed is `mer7m8` of the guide of the case, and the fields of the file give the other arguments:

- `min_length` is `compensatory_region_length_min`, and `max_length` is `compensatory_region_length_max`. The text `NA` is the R default.
- `limit` is `limit`, and `seed_mismatch` is `mismatch_`.
- `seq_id` gives the names of the `DNAStringSet`, and `sequence` gives its sequences.

A case where the patched function stops has no row. The file has a note `<case>: the patched function stops with: <message>` for it. The case `no_site` has such a note, because no seed match has a region.

The comment lines of `compensatory.tsv` give the version of `dplyr`, the patch, and these R results:

- `original_vignette`: the unpatched function for the `vignette` case gives 0 rows and the message `No compensatory regions found.`
- `line_99`: the error of `Biostrings::reverseComplement(gsub("U", "T", get_seed(GUIDE, "mer7m8")$Target.Seq))`. This is the error that the second change of the patch prevents.
- `repeated_name_r`: the number of times that `seqs[c("g", "g")]` holds the first sequence. The value is 2. R selects the sequences with a seed match by name.
- `relative_position_seq_g6`: `SeedMatchR:::relative_position_seq("TATAGA", "UUAUAGAUAUAGACACUGUUUUA")`, with the guide of `region_at_g2_r`, gives the start `2,8` and the end `7,13`.
- `min_17`: the substrings that the loop of `R/get_compensatory_sequences.R:52-54` makes for `compensatory_region_length_min = 17`. These are the 16 bases of g8 to g23, three times. The maximum is 16, thus the loop runs the lengths 17 and 16. The length 17 gives the starts 1 and 0, and the length 16 gives the start 1.
- `gc_content`: `.gc_content("GCAT")` gives 50, and `.gc_content("")` gives `NaN`.

### Wobble and bulge pipeline

`wobble_bulge_r.tsv` has the columns `case`, `guide`, and `sequence` of `tests/data/wobble_bulge_cases.tsv`. It also has the columns `r_rule_bulges_lim_1` and `r_rule_bulges_lim_4`. They hold the count of the R rule for `bulges_lim` 1 and 4, with `nt = 3`.

R stops before it counts. The pipeline gives `.get_target_bulged_guides` its arguments in the wrong order (`R/run_wobble_bulges_pipeline.R:27`). Also, `.sequence_with_tb` reads `bulge_allowance` from the global environment, where the variable does not exist.

To show the count rule of `:35-55`, the script does the loop of `.get_target_bulged_guides` itself. It gives the length of the guide and `bulges_lim` in the documented order. Before each call of `.sequence_with_tb`, it sets `bulge_allowance` in the global environment. After the call, it removes the variable.

The count of a case is a sum of `vcountPattern(..., fixed = FALSE)` for the sequence of the case. The sum goes over the wobble seeds of `.get_all_combination_w_b` and over the strings that `check_bulges_nt(..., 3)` keeps.

A bulged string of R has `N` in place of seed letters (`.sequence_with_tb`). Thus R counts a mismatch and not an extra target base. Each wobble seed adds its own count, thus one site counts more than one time.

The comment lines of `wobble_bulge_r.tsv` give the version of `gtools` and these R results:

- `pipeline_error` and `pipeline_error_bulges_lim_1`: the error of `run_wobble_bulges_pipeline(GUIDE, res = d$res, features = d$seqs, silent = TRUE)`, with the default `bulges_lim` and with 1. Each gives `object 'bulge_allowance' not found`. The script wraps each call in `utils::capture.output`, because the pipeline prints a line also with `silent = TRUE` (`:22`).
- `wobble_seeds_guide`, `wobble_seeds_4_u`, and `wobble_seeds_5_u_in_seed`: the `seed_combinations` of `.get_all_combination_w_b` for `GUIDE`, for `AUAUAGAGCAAGAACACAGAUAU`, and for `AUUUUUAGCAAGAACACAGAUAA`.
- `wobble_seeds_2_u` and `wobble_seeds_0_u`: the error of `.get_all_combination_w_b` for `AUAGAGAGCAAGAACACAGAUAA` and for `CAGAGAGCAAGAACACAGACAAC`.
- `formals_get_target_bulged_guides`: the names of the formals of `.get_target_bulged_guides`, in their order.
- `check_bulges_nt`: `check_bulges_nt(RNAStringSet(c("CNUCUAUA", "CUCUANUA", "CUCNUAUA")), 3)` keeps `CUCUANUA` and `CUCNUAUA`.

## Patch of vmatchPattern2

`vmatchPattern2` of SeedMatchR ignores `with.indels` and `min.mismatch` (`docs/divergences.md`, `D-001`). The file `tools/golden/patches/vmatchPattern2.R` defines a patched function with the same arguments and the same data frame. `full_search.R` puts it into the SeedMatchR namespace with `assignInNamespace`. As a result, the ranges in the golden files obey the same rules as the counts. A script that needs the original function does not source the patch.

## Patches of the annotation functions

`run_config` applies a patch only to the configs that name the patch. After each of these configs, it replaces the patched function with the original function.

- `tools/golden/patches/get_longest_utr.R`: `.get_longest_utr` of SeedMatchR selects the longest transcript of each gene from all the transcripts of the filter. If a gene has no coding transcript, the function stops with `subscript contains invalid names` (`D-028`). The patched function selects only from the transcripts of the feature result. A note in the golden file gives the error of the original function.
- `tools/golden/patches/reduce_features.R`: `.reduce_features` of SeedMatchR gives the ranges of a gene in the order of the position. On the minus strand, `extractTranscriptSeqs` then joins the ranges in the wrong order (`D-027`). The patched function reverses the ranges of each gene on the minus strand. Thus the ranges come in transcript order.

## Write a golden file

A script in `tools/golden/R/` starts with `source("tools/golden/R/common.R")`. Then it calls these functions:

- `load_test_data()` returns the list `res`, `seqs`, and `guide` from `tests/data/`.
- `write_golden(df, path, script, notes)` writes the data frame `df` to `tests/golden/<path>`. The argument `script` is the file name of the script, for example `check_env.R`. Each text in the optional argument `notes` becomes one more comment line. If `path` ends with `.gz`, the file is gzip. The gzip header has no time stamp and no file name, thus each run writes the same bytes.
- `with_seedmatchr_patch(patch, code)` puts each function of the file `patch` into the SeedMatchR namespace. Then it evaluates `code` and replaces the patched functions with the original functions.
- `stop_on_nan(df)` stops the script if a double column of `df` holds `NaN`, because `write_golden` writes a `NaN` as text.
- `one_line(text)` replaces each run of white space with one space, thus a note stays one comment line.
- `error_of(expr)` returns the message of the error of `expr`. For an error that dplyr chains, it returns the message of the first error in the chain. The script stops if `expr` gives no error.

An annotation script also sources `tools/golden/R/annotation_common.R`, which has these values and functions:

- `local_anno_db(txdb, dna)` returns the replacement of `load_species_anno_db` for the EnsDb `txdb` and the 2bit file `dna`.
- `ANNOTATION_CONFIGS` holds the configs of `load_annotations`. A config gives the feature, the arguments, and the optional patch.
- `run_config(name, config)` runs `load_annotations` for the config, with the patch of the config. It returns a list. The element `result` holds the result of `load_annotations` for each value of `return_gene_name`. The element `notes` holds the notes for the golden files.
- `sequences_table(result)` returns the names, the widths, and the sequences of the element `result` of the `run_config` list.

A golden file is a TSV file. Its first lines are comments that give the script, the R version, the Biostrings version, and the SeedMatchR version. Read it in Python with `comment="#"` and `dtype=str`.

A double in the file has 17 significant digits. A correctly rounded reader, for example Python `float()`, gives the same double as R.
