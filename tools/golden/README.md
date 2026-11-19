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

For each config of `annotation.R`, the script writes the feature ranges, the sequences, and the transcripts of the filter. The script also writes `genome_ranges.tsv`. For each range, this file holds the sequence from the 2bit reader of rtracklayer and from the BGZF reader of Rsamtools. rtracklayer 1.70.1 cannot read the 2bit file of version 1, thus the column `twobit_v1` holds `unsupported`.

The Python tests do not compare the configs whose names end with `_r`. These golden files show the R behavior that a divergence changes.

## Patch of vmatchPattern2

`vmatchPattern2` of SeedMatchR ignores `with.indels` and `min.mismatch` (`docs/divergences.md`, `D-001`). The file `tools/golden/patches/vmatchPattern2.R` defines a patched function with the same arguments and the same data frame. `full_search.R` puts it into the SeedMatchR namespace with `assignInNamespace`. As a result, the ranges in the golden files obey the same rules as the counts. A script that needs the original function does not source the patch.

## Patches of the annotation functions

`annotation.R` applies a patch only to the configs that name the patch. After each of these configs, it replaces the patched function with the original function.

- `tools/golden/patches/get_longest_utr.R`: `.get_longest_utr` of SeedMatchR selects the longest transcript of each gene from all the transcripts of the filter. If a gene has no coding transcript, the function stops with `subscript contains invalid names` (`D-028`). The patched function selects only from the transcripts of the feature result. A note in the golden file gives the error of the original function.
- `tools/golden/patches/reduce_features.R`: `.reduce_features` of SeedMatchR gives the ranges of a gene in the order of the position. On the minus strand, `extractTranscriptSeqs` then joins the ranges in the wrong order (`D-027`). The patched function reverses the ranges of each gene on the minus strand. Thus the ranges come in transcript order.

## Write a golden file

A script in `tools/golden/R/` starts with `source("tools/golden/R/common.R")`. Then it calls these functions:

- `load_test_data()` returns the list `res`, `seqs`, and `guide` from `tests/data/`.
- `write_golden(df, path, script, notes)` writes the data frame `df` to `tests/golden/<path>`. The argument `script` is the file name of the script, for example `check_env.R`. Each text in the optional argument `notes` becomes one more comment line.

An annotation script also sources `tools/golden/R/annotation_common.R`, which has these functions:

- `local_anno_db(txdb, dna)` returns the replacement of `load_species_anno_db` for the EnsDb `txdb` and the 2bit file `dna`.
- `with_seedmatchr_patch(patch, code)` puts each function of the file `patch` into the SeedMatchR namespace. Then it evaluates `code` and replaces the patched functions with the original functions.

A golden file is a TSV file. Its first lines are comments that give the script, the R version, the Biostrings version, and the SeedMatchR version. Read it in Python with `comment="#"` and `dtype=str`.

A double in the file has 17 significant digits. A correctly rounded reader, for example Python `float()`, gives the same double as R.
