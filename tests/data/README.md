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

## Fixtures of the small analyses

`kmer_cases.tsv` and `mirdb_sample.txt` are synthetic. They are not copies of SeedMatchR data or of miRDB data. A tab separates the fields. Each line ends with a newline.

`kmer_cases.tsv` has a header line with the fields `seq_id`, `gene_id`, and `sequence`. Each row is an edge case of the k-mer count:

- `s1` has the base `N`, and `s7` has the IUPAC bases `M` and `K`.
- `s2` has lowercase bases and the base `R`.
- `s3` is shorter than most widths.
- `s4` has an empty sequence. Its line ends with a tab.
- `s5` has one repeated base, and `s8` has one repeated pair of bases.
- `s6` has a k-mer that occurs two times.
- `s1` and `s2` share the gene `g1`, and `s3` and `s6` share the gene `g2`.

`mirdb_sample.txt` has the format of the miRDB file. It has no header, and its fields are the ID, the RefSeq ID, and the score. The IDs and the RefSeq IDs are not copies of miRDB data. The rows cover these cases:

- The species `hsa`, `mmu`, `rno`, `cfa`, and `gga`.
- An ID with more than one hyphen after the species (`rno-let-7a-2-3p`).
- An ID without a hyphen (`lin4`).
- Integer scores and scores with decimals.

The R golden files in `tests/golden/analyses/` come from these files. `tools/golden/README.md` gives the command.

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
