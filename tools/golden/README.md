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

An annotation script also sources `tools/golden/R/annotation_common.R`, which has these values and functions:

- `local_anno_db(txdb, dna)` returns the replacement of `load_species_anno_db` for the EnsDb `txdb` and the 2bit file `dna`.
- `with_seedmatchr_patch(patch, code)` puts each function of the file `patch` into the SeedMatchR namespace. Then it evaluates `code` and replaces the patched functions with the original functions.
- `ANNOTATION_CONFIGS` holds the configs of `load_annotations`. A config gives the feature, the arguments, and the optional patch.
- `run_config(name, config)` runs `load_annotations` for the config, with the patch of the config. It returns a list. The element `result` holds the result of `load_annotations` for each value of `return_gene_name`. The element `notes` holds the notes for the golden files.
- `sequences_table(result)` returns the names, the widths, and the sequences of the element `result` of the `run_config` list.

A golden file is a TSV file. Its first lines are comments that give the script, the R version, the Biostrings version, and the SeedMatchR version. Read it in Python with `comment="#"` and `dtype=str`.

A double in the file has 17 significant digits. A correctly rounded reader, for example Python `float()`, gives the same double as R.
