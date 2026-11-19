# Test data

`DESEQ2_results.tsv` and `Test_seqs.fa` are copies of the test data of SeedMatchR. The bytes are the same as in the source.

- Source repository: SeedMatchR
- Source path: `tests/testthat/testdata/`
- Source commit: `2fc0c4f`

`DESEQ2_results.tsv` holds differential expression results. `Test_seqs.fa` holds the sequences that the seed search reads.

`de_filter_cases.tsv` is not a copy. The seedmatch tests and the R script `tools/golden/R/de_and_counts.R` read it. Each row is a boundary case of `filter_res`:

- a missing value
- a value equal to a cutoff
- a value below a cutoff
- a negative fold change

## License

SeedMatchR is distributed under the MIT License. Copyright (c) 2025 Eli Lilly and Company. The `LICENSE` file in the repository root holds the full text of the license.
