# seedmatch

[![Test](https://github.com/inflexa-ai/seedmatch/actions/workflows/test.yml/badge.svg)](https://github.com/inflexa-ai/seedmatch/actions/workflows/test.yml)
[![OpenSSF Scorecard](https://api.scorecard.dev/projects/github.com/inflexa-ai/seedmatch/badge)](https://scorecard.dev/viewer/?uri=github.com/inflexa-ai/seedmatch)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

seedmatch is a Python port of the R package [SeedMatchR](https://github.com/tacazares/SeedMatchR), ready for AI agent sandboxes. It finds the seed matches of an siRNA guide strand in transcript sequences. It then joins the matches with the results of a differential expression (DE) analysis. Thus you can see which genes with a seed match go down.

The output is the same as the output of SeedMatchR on the core paths. R golden files in `tests/golden/` keep this true. Each change of behavior from R has an entry in [`docs/divergences.md`](docs/divergences.md).

## A part of an agentic toolchain for bioinformatics

Inflexa, Inc. builds an agentic toolchain for bioinformatics. The toolchain ports existing tools to languages and architectures that an agent sandbox can run. An example of such a sandbox is the sandbox of [Inflexa](https://github.com/inflexa-ai/inflexa). seedmatch is a tool of this toolchain.

A sandbox such as the Inflexa sandbox has no network egress, and an agent corrects its calls from the error messages. Thus seedmatch obeys these rules:

- **No network in the core.** The core modules never use the network. Each function reads its input from a path or an object that you give. The separate module `seedmatch.fetch` downloads the resources (see [Resources](#resources)).
- **Python only.** The package is pure Python on numpy, pandas, numba, and scipy. It does not use R or a system library.
- **Errors that an agent can use.** Each `InputError` names the bad value and the expected value. Thus an agent can correct the call.

## Example

```python
from seedmatch import (
    Guide,
    SearchOptions,
    Seed,
    annotate_seed_counts,
    count_matches,
    filter_de,
    find_matches,
    read_deseq2_table,
    read_fasta,
    seed_match_report,
)

# The ID of each FASTA record is the gene ID of its sequence.
sequences = read_fasta("transcripts.fa")

# The guide strand, 5' to 3'. The target is the DNA site that the seed pairs with.
guide = Guide("UUAUAGAGCAAGAACACUGUUUU")
target = guide.target(Seed.MER7M8)

# The number of matches in each sequence, and the 1-based range of each match.
counts = count_matches(target, sequences)
ranges = find_matches(target, sequences, SearchOptions(max_mismatch=1))

# A DESeq2 table with 8 columns, filtered as filter_res does it in R.
de = filter_de(read_deseq2_table("deseq2_results.tsv"))

# Add the column "mer7m8": the count of the seed matches of each gene.
de = annotate_seed_counts(de, sequences, target)

# The SeedMatchReport of R: the genes of each match class, and the off-targets among them.
report = seed_match_report(de, sequences, guide)
summary = report.report
matches = report.matches
```

If a gene has more than one sequence, make the set with `SequenceSet.from_records(records, gene_ids=...)`.

`ecdf_test` compares the `log2FoldChange` values of two gene sets. This example compares the genes with a seed match to the genes without a match. With the default KS test, a small p-value shows that the genes with a seed match go down.

```python
from seedmatch import ecdf_test

# The genes with a seed match, and the genes without a match.
seed_genes = de.loc[de["mer7m8"] >= 1, "gene_id"]
background = de.loc[de["mer7m8"] == 0, "gene_id"]

# The KS test of R. The default alternative is GREATER, as in R.
ecdf_result = ecdf_test(de, target=seed_genes, background=background)
print(ecdf_result.statistic, ecdf_result.p_value)
```

Import the members with `from seedmatch import EcdfTest`. Then use `test=EcdfTest.WILCOXON`, `EcdfTest.WASSERSTEIN`, or `EcdfTest.DTS` for another test. The DTS test must have `rng`, a `numpy.random.Generator` or an int seed, because its p-value comes from permutations. `assign_groups(de, groups)` gives the group of each row of `de`.

The default alternative `GREATER` has another meaning for `EcdfTest.KS` and for `EcdfTest.WILCOXON`. With KS, a small p-value shows that the target values are lower. With WILCOXON, a small p-value shows that the target values are higher. In the test table of this repository (`tests/data/DESEQ2_results.tsv`), the genes with a seed match have lower values. For this table, KS with `GREATER` gives p = 1.04e-4, but WILCOXON with `GREATER` gives p = 0.9996. WILCOXON with `Alternative.LESS` gives the small p-value 3.9e-4. To find lower target values with WILCOXON, set `alternative=Alternative.LESS`. Import `Alternative` from `seedmatch`.

To get the sequences from an annotation, open a `Reference`. A `Reference` joins an EnsDb file and a genome file of the same build. The genome file is a 2bit file, or a FASTA file with its index files.

```python
from seedmatch import FeatureType, Reference, count_matches

with Reference.from_paths(ensdb="ensdb.sqlite", genome="genome.2bit") as ref:
    # The 3' UTR of each canonical protein-coding transcript.
    # The ID of a sequence is its transcript ID.
    utr3 = ref.sequences(FeatureType.UTR3)

# The number of seed matches in each 3' UTR, with the target of the example above.
utr3_counts = count_matches(target, utr3)
```

`Reference.sequences` also takes a `TranscriptFilter` and a `PerGene` option. The defaults are the defaults of `load_annotations` in R. `PerGene.LONGEST` selects, for each gene, the transcript with the longest feature. For `FeatureType.UTR3`, this is the transcript with the longest 3' UTR. `PerGene.UNION` merges the ranges of the transcripts of a gene. `Reference.from_directory(path)` opens the one `.sqlite` file and the one `.2bit` or `.fa.gz` file of a directory. `Reference.from_gtf(gtf, genome)` reads a GTF file instead of an EnsDb file.

[`docs/r-map.md`](docs/r-map.md) gives the Python name of each R function, and the changes of the arguments.

## Other analyses

`kmer_counts` counts the k-mers of each sequence. Its result holds a sparse matrix with one row for each sequence and one column for each k-mer. `kmer_gene_counts` gives the number of genes that have each k-mer. R draws a histogram of this number in `kmer_hist`. This example reads the number for the target of the example above.

```python
from seedmatch import kmer_counts, kmer_gene_counts

# The 7-mers of the sequences. The default k is 7, the default of R.
kmers = kmer_counts(sequences)

# The number of genes with each 7-mer, as a Series. The index is the 7-mer.
genes_with_kmer = kmer_gene_counts(kmers)

# The number of genes that have the 7-mer of the target.
print(genes_with_kmer[target.sequence])
```

`kmers.matrix` is a `scipy.sparse.csr_array`. `kmers.to_frame()` gives the layout of `sequence_kmer_counts` in R. This frame is dense and holds n · 4^k · 8 bytes for n sequences. The value of `k` is at most 12 (D-043).

`alignment_text` shows how the guide pairs with the sequence at one match. It takes a row of `find_matches` and the seed of the search. This example shows the first match of the example above.

```python
from seedmatch import alignment_text

# The three lines: the guide, the pairs, and the sequence. The default seed is Seed.MER7M8.
alignment = alignment_text(find_matches(target, sequences).iloc[0], guide, sequences)
print(alignment)
```

Use `seed=` for the seed of another search. `left_pad` and `right_pad` show more bases of the sequence on each side, and each pad is at most the length of the sequence. A `|` shows a Watson-Crick pair. A position outside the sequence is a space (D-046).

`match_distances` gives the distance between each pair of matches in one sequence. This example uses the matches of the `MER6` seed.

```python
from seedmatch import match_distances

# One row for each pair of matches in a sequence. A sequence with one match has no row (D-011).
distances = match_distances(find_matches(guide.target(Seed.MER6), sequences))
print(distances.head())
```

The columns are `seq_id`, `gene_id`, and `distance`. The default `anchor=MatchAnchor.FIVE_PRIME_END` is the R default `fivePrimeEnd`, and it uses the `end` of each match. `MatchAnchor.THREE_PRIME_END` uses the `start`, and `MatchAnchor.CENTER` uses the center. A sequence with n matches gives n(n − 1)/2 rows, thus a sequence with many matches gives a large result.

`e_values` gives the e-value of the seed count of each gene. The e-value is the fraction of random guides whose count is at least the observed count. `e_value_null` makes the random counts of `annotate_seed_counts`. It shuffles the letters of the guide, and it counts the seed matches of each shuffle in each gene. This example uses the `de` table of the example above.

```python
from seedmatch import e_value_null, e_values

# The null: the seed count of each gene for each of 100 shuffles of the guide.
null = e_value_null(guide, sequences, rng=1, n_random=100)

# The e-value of the column "mer7m8". A small e-value shows a count that shuffles rarely reach.
scored = e_values(de, null, column="mer7m8")
print(scored[["gene_id", "mer7m8", "e_value"]].head())
```

`rng` is a `numpy.random.Generator` or an int seed, and the default `n_random` is 1000. The null is a required argument of `e_values`, because the R default does not exist. Genes match by `gene_id` (D-008). The result has the rows of `de` whose gene is in the null. `e_value_null` runs its searches on a thread pool, and `workers` sets the width of the pool. The null holds genes x `n_random` x 4 bytes. A count that is NaN, infinite, or negative raises `InputError`, because R gives it the e-value 0. This holds also for a row whose gene is not in the null, where R drops the row (D-047).

`e_values` cannot tell which count a null holds. Each null fits one count. `e_value_null` fits the count of `annotate_seed_counts`, and `wobble_bulge_null` fits the column `overall_matches` of `annotate_wobble_bulge_counts`. A null of the other count gives wrong e-values and no message. For example, the null of `e_value_null` gives e-values that are too small for `overall_matches`.

`benchmark_predictions` scores the seed counts against a gold standard of genes that go down. A gene is in the gold standard if one of its rows has `padj` at most `padj_cutoff` and `log2FoldChange` below `l2fc_cutoff`. This example uses the `de` table of the example above.

```python
from seedmatch import benchmark_predictions

# The scores are the seed counts of the column "mer7m8". The cutoffs are the R defaults: 0.15 and 0.
benchmark = benchmark_predictions(de)
print(benchmark.pr.auc_integral, benchmark.roc.auc)
```

`benchmark.pr` is a `PrCurve`, and `benchmark.roc` is a `RocCurve`. The PR curve and the ROC curve use the same labels, but the ROC curve of R does not (D-006). A row with a NaN `padj` or a NaN `log2FoldChange` does not pass the cutoffs. It is a negative unless another row of its gene passes, as in R. The negatives are the other rows of `de`, thus a filter of `de` before the call also changes the negatives. A NaN score raises `InputError`. So do a cutoff that is not a number, and a gold standard without a positive or without a negative row (D-048).

`read_mirdb` reads the miRDB table of predicted miRNA targets. This is the file that `seedmatch-fetch mirdb` downloads, and it has no header. This example keeps the rows of the rat, which has the species code `rno`.

```python
from seedmatch import read_mirdb

# The columns: mirdb_id, target_refseq_id, target_score, species, and mirna_name.
rat_targets = read_mirdb("miRDB_v6.0_prediction_result.txt.gz", species="rno")
print(rat_targets.head())
```

`species` is the text of the miRDB ID before the first `-`, and `mirna_name` is the text after it. A code without a row raises `InputError` that lists the codes of the file. The reader reads 2^18 lines at a time. The result holds about 40 bytes for each row and the distinct strings. A score that is not a finite decimal number raises `InputError`. The ID fields keep `#`, quote characters, and the text `NA` as text (D-049).

`wobble_bulge_matches` finds the sites that have G:U pairs or a bulge in the seed region. A guide U can pair with a target A or G. Thus one pattern with the letter `R` covers each set of G:U pairs. A bulge is a run of extra target bases between two seed positions. This example uses the `guide`, the `sequences`, and the `de` table of the examples above.

```python
from seedmatch import annotate_wobble_bulge_counts, wobble_bulge_matches

# One row for each distinct start of a site, with the place and the size of its bulge.
wobble_matches = wobble_bulge_matches(guide, sequences)
print(wobble_matches.head())

# The column "overall_matches": for each gene, the largest number of rows of one of its sequences (D-003).
wobble_de = annotate_wobble_bulge_counts(de, sequences, guide)
print(wobble_de[["gene_id", "overall_matches"]].head())
```

`max_bulge` is the longest run of extra bases, and its default is 4. `bulge_free` is the number of seed positions from g2 that no bulge is next to. Its default is 3, thus a bulge lies only between g5 and g6, g6 and g7, or g7 and g8. A sequence has one row for each distinct start of a match, and the row names the first pattern that matches at this start. `wobble_bulge_matches` searches the set one time for each pattern, and `workers` sets the width of the thread pool. R stops for each call of its wobble and bulge search, and its count rule counts one site more than one time (D-005). `wobble_bulge_null(guide, sequences, rng=1, n_random=100)` makes the null of `overall_matches` for `e_values`. Its shuffles are those of `e_value_null`, and its thread pool runs over the shuffles.

`compensatory_sites` finds the sites with a seed match and a second region that pairs with the sequence. A region is a window of the guide that starts after the seed window, at g9 for the `MER7M8` seed. The target of a region is the reverse complement of the window. It lies before the seed match in the same sequence. This example uses the `guide`, the `sequences`, and the `de` table of the examples above.

```python
from seedmatch import compensatory_sites, ecdf_test

# One row for each pair of a seed match and a region of one sequence.
sites = compensatory_sites(guide, sequences)
print(sites.head())

# As in the R vignette: the genes with a site against the other genes of the DE table.
site_genes = sites["gene_id"].unique()
other_genes = de.loc[~de["gene_id"].isin(site_genes), "gene_id"]
site_result = ecdf_test(de, target=site_genes, background=other_genes)
print(site_result.p_value)
```

`min_length` is the shortest region, and its default is 6. `max_length` is the longest region, and its default is each position up to the guide end. `first_position` is the first guide position of a region. R always uses g8, and `first_position=8` gives the R regions (D-056). `max_distance` is the largest number of bases between a region and its seed match. R `limit` L gives `max_distance` L - 2, thus the default 98 is the R default 100.

A region is dropped if a longer region of the same seed match contains it (D-053). The seed match can have `max_seed_mismatch` mismatches, and the region can have `max_region_mismatch`. R returns no row for each input (D-004). `workers` sets the width of the thread pool. The function takes the fields above as a `CompensatoryOptions` object, in the third position.

## Resources

The core of seedmatch reads files from a path, and the sandbox of an agent has no network access. Thus a machine with network access runs `seedmatch-fetch`, and the sandbox reads the files from a path. If you install seedmatch, you get the command `seedmatch-fetch`.

The registry in `seedmatch.resources` pins the URL, the size, and the SHA-256 digest of each file. The command downloads a file only through `https`. A download that does not agree with the registry leaves no file. If a file is present and its digest agrees, the command does not download it again. If the digest differs, the command stops with an error and keeps the file.

```bash
# The 2bit file and the EnsDb file of one genome build.
seedmatch-fetch reference rnor7 refs/rnor7

# The DESeq2 tables of the siRNA example (GEO series GSE184929).
seedmatch-fetch example sirna data/sirna

# The miRDB v6.0 predictions.
seedmatch-fetch mirdb data/mirdb

# The TRRUST table of human or mouse.
seedmatch-fetch trrust human data/trrust

# The STRING protein links of one NCBI taxon.
seedmatch-fetch string 10116 data/string

# One line for each file: the group, the key, the file name, the size in bytes, and the description.
seedmatch-fetch list
```

`python -m seedmatch.fetch` runs the same command. The files of a reference are large, and `seedmatch-fetch list` gives the size of each file.

If `DEST` is absent, the command makes it. For each file, it prints the path, a tab, and `downloaded` or `present` on standard output. Before a download, it prints the file name and the size on standard error. The exit status is 0 for success and 1 for an error. An error prints `seedmatch-fetch: <message>` on standard error. Ctrl-C ends the command with status 130. A usage error exits with status 2. A key that the registry does not have gives a message that lists the valid keys.

Use one `DEST` for each reference build. `Reference.from_directory(DEST)` opens the one `.sqlite` file and the one `.2bit` file of the directory. A directory with the EnsDb files of different builds has more than one `.sqlite` file, and `from_directory` raises `InputError`. `read_deseq2_table` reads each table of `example sirna` from `DEST`, because it accepts the gzip format. `read_mirdb` reads the miRDB file from `DEST`. No function of this version reads the TRRUST and STRING files.

## Differences from SeedMatchR

seedmatch corrects the defects of SeedMatchR that it finds. These changes are the most important for a user:

- The ranges of a search obey the same rules as the counts. In R, `vmatchPattern2` uses the indel algorithm also when `with.indels` is `FALSE`, if the pattern is longer than `max.mismatch`. Thus the ranges and the counts do not agree (D-001).
- A result for a gene uses all the sequences of the gene. R uses the first sequence for the seed count, and it counts each sequence in the report (D-003).
- `full_search_ranges` gives a result for each input. R stops if each search finds a match, or if no search finds a match (D-022).
- `wobble=True` with `fixed_pattern=True` raises `InputError`, because a fixed pattern matches the wobble letter R only to the letter R (D-026).
- A sequence from an annotation keeps its transcript ID and its gene ID. In R, `load_annotations` gives the name `NA` to each sequence of a gene union (D-002), and `tx.names` is `NULL` (D-009).
- `ecdf_test` drops the rows without a fold change for each test. R drops them for the KS and Wilcoxon tests only, and it stops for DTS and Wasserstein (D-036). `ecdf_test` raises `InputError` for an `alternative` with DTS or Wasserstein, but R ignores it without a message (D-035).
- `kmer_gene_counts` counts the genes with a k-mer. The histogram of R counts the sequences, but its label says genes (D-042).
- `benchmark_predictions` applies the `padj` cutoff. In R, `padj <= padj` compares the column with itself, thus each row with a `padj` passes. R also gives the labels to the ROC curve as the scores of the negative class (D-006).
- `alignment_text` shows the full guide, from g1 to the last base. In R, the guide line of `plotAlignment` has no g1, and the function stops for the default positions (D-045).
- The wobble and bulge search runs. In R, each call of `run_wobble_bulges_pipeline` stops with an error. The port counts each distinct start of a site one time, and a bulge is a run of extra target bases (D-005).
- The 3' compensatory search runs. R returns 0 rows for each input. Without the early return, its code stops at the first site. The port gives each pair of a seed match and a region (D-004).
- An input that R accepts in silence, but that gives a wrong result, raises `InputError`. The message gives the bad value. Examples are a gap letter in a sequence (D-017), a DESeq2 table with the wrong columns (D-021), and a NaN cutoff (D-024).

## Scope

This version has the core path:

- the FASTA reader and the sequence set
- the seeds and the targets of a guide
- the search: the counts and the ranges, with mismatches and indels
- the annotation: the EnsDb file or an Ensembl or GENCODE GTF file, the 2bit and FASTA genome files, and the feature sequences of a `Reference`
- the DE table: the reader and the filter
- the seed counts, `full_search`, `full_search_ranges`, and `seed_match_report`
- the statistics: the KS and Wilcoxon tests, DTS, the Wasserstein areas, and `ecdf_test` for two gene sets
- the PR and ROC curves of scores and 0/1 labels
- the k-mer counts of a sequence set, and the number of genes with each k-mer
- the alignment text of a seed match
- the distances between the matches of a sequence
- the e-value of a seed count against a null of shuffled guides
- the benchmark of the seed counts against a gold standard of genes that go down
- the reader of the miRDB table of predicted miRNA targets
- the wobble and bulge search: the sites with G:U pairs or a bulge, their count for each gene, and the null of this count
- the 3' compensatory search: the sites with a seed match and, before it in the sequence, a match of a 3' region of the guide
- the registry of the upstream resources, and the command `seedmatch-fetch` that downloads them

These parts of SeedMatchR are not in this version:

- the network scores
- the plots and the HTML report

## Installation

Python 3.12 or later is necessary. Install seedmatch from GitHub:

```bash
uv pip install git+https://github.com/inflexa-ai/seedmatch.git
```

To install it from a clone of this repository, use `uv pip install .`.

If numba has no writable cache directory, `import seedmatch` stops with a `RuntimeError`. Thus set `NUMBA_CACHE_DIR`, or use a writable `HOME`.

## Development

Install the dependencies:

```bash
uv sync
```

Do the checks:

```bash
uv run ruff check
uv run ruff format --check
uv run mypy
uv run pytest
```

Measure the coverage:

```bash
PYTHONPATH=tools/coverage NUMBA_DISABLE_JIT=1 uv run pytest -p serial_pool --cov
```

Without the JIT, the coverage tracer stops the pool threads on its data lock, thus the plugin `tools/coverage/serial_pool.py` runs the pool in the caller thread.

`tools/golden/README.md` tells how to make the R golden files again.

[`CONTRIBUTING.md`](CONTRIBUTING.md) gives the rules for a contribution. To report a security problem, obey [`SECURITY.md`](SECURITY.md).

## License

seedmatch uses the MIT License, the same license as SeedMatchR. See [`LICENSE`](LICENSE).
