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

Use one `DEST` for each reference build. `Reference.from_directory(DEST)` opens the one `.sqlite` file and the one `.2bit` file of the directory. A directory with the EnsDb files of different builds has more than one `.sqlite` file, and `from_directory` raises `InputError`. `read_deseq2_table` reads each table of `example sirna` from `DEST`, because it accepts the gzip format. No function of this version reads the miRDB, TRRUST, and STRING files.

## Differences from SeedMatchR

seedmatch corrects the defects of SeedMatchR that it finds. These changes are the most important for a user:

- The ranges of a search obey the same rules as the counts. In R, `vmatchPattern2` uses the indel algorithm also when `with.indels` is `FALSE`, if the pattern is longer than `max.mismatch`. Thus the ranges and the counts do not agree (D-001).
- A result for a gene uses all the sequences of the gene. R uses the first sequence for the seed count, and it counts each sequence in the report (D-003).
- `full_search_ranges` gives a result for each input. R stops if each search finds a match, or if no search finds a match (D-022).
- `wobble=True` with `fixed_pattern=True` raises `InputError`, because a fixed pattern matches the wobble letter R only to the letter R (D-026).
- A sequence from an annotation keeps its transcript ID and its gene ID. In R, `load_annotations` gives the name `NA` to each sequence of a gene union (D-002), and `tx.names` is `NULL` (D-009).
- `ecdf_test` drops the rows without a fold change for each test. R drops them for the KS and Wilcoxon tests only, and it stops for DTS and Wasserstein (D-036). `ecdf_test` raises `InputError` for an `alternative` with DTS or Wasserstein, but R ignores it without a message (D-035).
- `kmer_gene_counts` counts the genes with a k-mer. The histogram of R counts the sequences, but its label says genes (D-042).
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
- the k-mer counts of a sequence set, and the number of genes with each k-mer
- the registry of the upstream resources, and the command `seedmatch-fetch` that downloads them

These parts of SeedMatchR are not in this version:

- the PR and ROC curves
- the wobble and bulge search, the 3' compensatory search, the network scores, and the other small analyses
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
NUMBA_DISABLE_JIT=1 uv run pytest --cov
```

`tools/golden/README.md` tells how to make the R golden files again.

[`CONTRIBUTING.md`](CONTRIBUTING.md) gives the rules for a contribution. To report a security problem, obey [`SECURITY.md`](SECURITY.md).

## License

seedmatch uses the MIT License, the same license as SeedMatchR. See [`LICENSE`](LICENSE).
