# Test data

`DESEQ2_results.tsv` and `Test_seqs.fa` are copies of the test data of SeedMatchR. The bytes are the same as in the source.

- Source repository: SeedMatchR
- Source path: `tests/testthat/testdata/`
- Source commit: `2fc0c4f`

`SeedMatchR.NAMESPACE` is a copy of the file `NAMESPACE` of SeedMatchR at commit `2fc0c4f`. The bytes are the same as in the source. `tests/test_r_map.py` reads its `export(...)` lines.

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

## Fixtures of the non-canonical searches

`compensatory_cases.tsv` and `wobble_bulge_cases.tsv` are synthetic. They are not copies of SeedMatchR data. A tab separates the fields. Each line ends with a newline.

`compensatory_cases.tsv` has a header line with the fields `case`, `guide`, `min_length`, `max_length`, `limit`, `seed_mismatch`, `seq_id`, and `sequence`. The rows of a case are one set of sequences. The fields from `guide` to `seed_mismatch` have the same value in each row of a case. The text `NA` in `max_length` is the R default, and it stays text. The seed is `mer7m8` of the guide.

A case that ends with `_r` shows the R behavior that a divergence changes. The cases:

- `seed_past_end`: the sequence ends with the seed target without its last base. R finds a match that ends one position after the sequence. It counts the outside position as the mismatch.
- `no_site`: the sequence has no seed match, thus no region has a site.
- `repeated_name_r`: two sequences have the name `g`, and each sequence has a site. R selects the sequences by name, thus it searches the first sequence two times and gives no row for the second sequence.
- `two_seed_matches_r`: one sequence has two seed matches. R keeps the longest region of each sequence and region position, and it ignores the seed match. Thus the longer region of one seed match drops the region of the other seed match.
- `repeated_region_r`: the guide holds the region `GCAAGC` two times, at g8 to g13 and at g12 to g17. R searches the region two times and gives two equal rows. Each row lists the position of each copy.
- `region_at_g2_r`: the guide holds the region `TATAGA` at g2 to g7 and at g8 to g13. R lists the position of each copy, and the first copy is in the seed.
- `seed_mm2_r`: the seed match has two mismatches. R gives the value 1 in `mismatch_in_seed`.
- `min_above_max_r`: `min_length` 6 is above `max_length` 3. R loops over the lengths 6, 5, 4, and 3.

`wobble_bulge_cases.tsv` has a header line with the fields `case`, `guide`, and `sequence`. Each row is a guide and a sequence. Except for `five_wobbles`, the guide is `UUAUAGAGCAAGAACACUGUUUU`, the seed target is `CTCTATA`, and `G` bases are before and after the site. The cases:

- `perfect`: the site is the seed target.
- `wobble_g4`: the site has `G` in place of the `A` that pairs with g4. This is a G:U wobble.
- `mismatch_t4`: the site has `C` in place of the `T` at t4, which pairs with g5. This is a mismatch and not a wobble.
- `bulge_after_g5`: the site has one extra base between the target bases of g5 and g6.
- `bulge_after_g4`: the site has one extra base between the target bases of g4 and g5, next to g4.
- `bulge_4_after_g7`: the site has four extra bases between the target bases of g7 and g8.
- `bulge_5_after_g7`: the site has five extra bases at the same place. This is one more than the R default of `bulges_lim`.
- `one_start_two_patterns`: more than one pattern of the wobble set matches at the start of the site.
- `five_wobbles`: the guide is `AUUUUUAGCAAGAACACAGAUAA`, with `U` at g2 to g6. The site has `G` at the five target positions of these `U`.

The R golden files in `tests/golden/noncanonical/` come from these files. `tools/golden/README.md` gives the command.

## Annotation fixture

The files in `annotation/` are synthetic. They are not copies of SeedMatchR data. The script `tools/fixtures/annotation.py` writes each file from one specification:

- `synthetic.sqlite`: an EnsDb with schema 2.2, with the table `protein` (see below)
- `synthetic.2bit` and `synthetic.v1.2bit`: the genome in the 2bit format of version 0 and of version 1
- `synthetic.fa`: the genome in the FASTA format, with the soft-masked bases in lowercase
- `synthetic.fa.gz`: the BGZF file of `synthetic.fa`
- the `.fai` and `.gzi` files: the indexes that `samtools faidx` and `bgzip` write
- `synthetic.gtf` and `synthetic.gencode.gtf.gz`: the annotation in the Ensembl style and in the GENCODE style. Each `CDS` line has the attribute `protein_id` of its transcript, the protein of the table `protein`.

To make the files again, run the script from the repository root. The script runs `bgzip` of htslib and `samtools`, thus they must be on `PATH`:

```bash
uv run python tools/fixtures/annotation.py
```

The R golden files in `tests/golden/annotation/` come from these files. `tools/golden/README.md` gives the command.

The table `protein` has the columns of the table in `EnsDb.Rnorvegicus.v109.sqlite`: `tx_id`, `protein_id`, and `protein_sequence`. The column `protein_sequence` is `NULL`. The EnsDb has no `uniprot` table and no `protein_domain` table. Thus `hasProteinData` of ensembldb is `FALSE`, and the queries of ensembldb do not change.

Each transcript with a CDS has one protein. The protein ID is `ENSSYNP<n>` for the transcript number n, and n has 11 digits. The lncRNA transcript 5 has no protein. The genes 1, 2, 5, and 7 have more than one protein.

## Fixtures of the network functions

`trrust_sample.tsv` and `string_links_sample.txt` are synthetic. They are not copies of TRRUST data or of STRING data. Each line ends with a newline.

`trrust_sample.tsv` has the format of the TRRUST file. It has no header, and a tab separates the fields. The fields are the regulator, the target, the mode, and the PubMed IDs. The rows cover these cases:

- Line 1 is an edge that R loses with `header = TRUE`.
- `TFA` represses and activates `GA` (lines 1 and 6).
- Lines 2 and 9 give the same edge with other PubMed IDs.
- `TFB` regulates `GC` with the mode `Unknown`.
- `GA` is a target and a regulator (line 5), and `TFC` regulates itself (line 8).

`string_links_sample.txt` has the format of the STRING file. It has a header line, and a space separates the fields. The proteins are the proteins of the synthetic EnsDb, with its taxon 32630. The proteins 98 and 99 are not in the EnsDb. Each pair has a row for each order. The rows are in the order of `protein1`, as in the rat file.

Let Pn be `32630.ENSSYNP<n>` and Gn be `ENSSYNG<n>`, each with 11 digits. The cases:

- P1 and P2 are proteins of G1, thus the link P1-P2 joins G1 with itself.
- The links P2-P7 of G1, P4-P12 of G2, and P11-P12 of G7 come from a second protein of a gene.
- P6-P99 and P98-P99 have a protein without a gene.
- The scores 400 and 399 are at the default threshold of STRINGdb 2.22.0 (`R/rstring.R:109-110`) and below it. The score 150 is the lowest score of the rat file.
- The scores of 400 or more give the cutoff 400 at the percentile 0, 650 at 0.5, and 897.4999999999999 at 0.95.

The R golden files in `tests/golden/network/` come from these files and from the table `protein`. `tools/golden/README.md` gives the command.

## License

SeedMatchR is distributed under the MIT License. Copyright (c) 2025 Eli Lilly and Company. The `LICENSE` file in the repository root holds the full text of the license.
