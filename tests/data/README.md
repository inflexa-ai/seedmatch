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

## Annotation fixture

The files in `annotation/` are synthetic. They are not copies of SeedMatchR data. The script `tools/fixtures/annotation.py` writes each file from one specification:

- `synthetic.sqlite`: an EnsDb with schema 2.2
- `synthetic.2bit` and `synthetic.v1.2bit`: the genome in the 2bit format of version 0 and of version 1
- `synthetic.fa`: the genome in the FASTA format, with the soft-masked bases in lowercase
- `synthetic.fa.gz`: the BGZF file of `synthetic.fa`
- the `.fai` and `.gzi` files: the indexes that `samtools faidx` and `bgzip` write
- `synthetic.gtf` and `synthetic.gencode.gtf.gz`: the annotation in the Ensembl style and in the GENCODE style

To make the files again, run the script from the repository root. The script runs `bgzip` of htslib and `samtools`, thus they must be on `PATH`:

```bash
uv run python tools/fixtures/annotation.py
```

The R golden files in `tests/golden/annotation/` come from these files. `tools/golden/README.md` gives the command.

## License

SeedMatchR is distributed under the MIT License. Copyright (c) 2025 Eli Lilly and Company. The `LICENSE` file in the repository root holds the full text of the license.
