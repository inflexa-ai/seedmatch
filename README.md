# seedmatch

[![Test](https://github.com/inflexa-ai/seedmatch/actions/workflows/test.yml/badge.svg)](https://github.com/inflexa-ai/seedmatch/actions/workflows/test.yml)
[![OpenSSF Scorecard](https://api.scorecard.dev/projects/github.com/inflexa-ai/seedmatch/badge)](https://scorecard.dev/viewer/?uri=github.com/inflexa-ai/seedmatch)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

seedmatch is a Python port of the R package [SeedMatchR](https://github.com/tacazares/SeedMatchR). It finds the seed matches of an siRNA guide strand in transcript sequences. It then joins the matches with the results of a differential expression (DE) analysis. Thus you can see which genes with a seed match go down.

The output is the same as the output of SeedMatchR on the core paths. R golden files in `tests/golden/` keep this true. Each change of behavior from R has an entry in [`docs/divergences.md`](docs/divergences.md).

## Offline use

The core modules never use the network. Each function reads its input from a path or an object that you give. Thus seedmatch runs in a sandbox without network access. A later version adds a separate module that downloads the resources.

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

[`docs/r-map.md`](docs/r-map.md) gives the Python name of each R function, and the changes of the arguments.

## Differences from SeedMatchR

seedmatch corrects the defects of SeedMatchR that it finds. These changes are the most important for a user:

- The ranges of a search obey the same rules as the counts. In R, `vmatchPattern2` uses the indel algorithm also when `with.indels` is `FALSE`, if the pattern is longer than `max.mismatch`. Thus the ranges and the counts do not agree (D-001).
- A result for a gene uses all the sequences of the gene. R uses the first sequence for the seed count, and it counts each sequence in the report (D-003).
- `full_search_ranges` gives a result for each input. R stops if each search finds a match, or if no search finds a match (D-022).
- An input that R accepts in silence, but that gives a wrong result, raises `InputError`. The message gives the bad value. Examples are a gap letter in a sequence (D-017), a DESeq2 table with the wrong columns (D-021), and a NaN cutoff (D-024).

## Scope

This version has the core path:

- the FASTA reader and the sequence set
- the seeds and the targets of a guide
- the search: the counts and the ranges, with mismatches and indels
- the DE table: the reader and the filter
- the seed counts, `full_search`, `full_search_ranges`, and `seed_match_report`

These parts of SeedMatchR are not in this version:

- the annotation: EnsDb, 2bit, and GTF files, and the feature sequences
- the statistics: the KS and Wilcoxon tests, DTS, the Wasserstein areas, and the PR and ROC curves
- the downloads
- the wobble and bulge search, the 3' compensatory search, the network scores, and the other small analyses
- the plots and the HTML report

## Installation

seedmatch needs Python 3.12 or later. Install it from GitHub:

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
