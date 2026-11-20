# Divergences from SeedMatchR

This document records each change of behavior from SeedMatchR. A change is a defect fix or an improvement. A change of behavior that has no entry is a defect of the port.

Each entry has a heading in the form `## D-NNN: <title>`. The entry has these fields:

- ID: a stable ID, for example `D-001`.
- R location: the R function and the R source location.
- R behavior: the behavior of the R code, with the evidence (an R golden or an R probe).
- New behavior: the behavior of the port.
- Reason: the reason for the change.
- Test: the test that keeps the new behavior. The name or the docstring of the test contains the ID.

The IDs of the known divergences have the same order as the rows of the table in `docs/design.md`, section "Divergences (defect fixes)". A new divergence gets the next free ID.

## Entries

## D-001: Ranges obey with.indels and min.mismatch

- ID: `D-001`
- R location: `vmatchPattern2` (`R/vmatchPattern2.R:43`). Its callers are `.get_match_df` (`R/SeedMatchR.R:260`), `full_search` (`R/full_search.R:71`), and `.run_scanning` of `full_search_ranges` (`R/full_search.R:205`).
- R behavior: `vmatchPattern2` passes the algorithm `"indels"` to the C matcher in each call. Thus its ranges disagree with the counts of `vcountPattern`. These R probes with Biostrings 2.78.0 give the evidence:
  - If `with.indels` is `TRUE` and `max.mismatch` is greater than 0, `vmatchPattern` stops with "vmatchPattern() does not support indels yet" for a `DNAStringSet`. `vmatchPattern2` calls the C matcher with "indels" for each call (`R/vmatchPattern2.R:43`).
  - `vmatchPattern2("ACGTACG", DNAStringSet(c(a = "TTTTACGACGTTTT")), max.mismatch = 1, with.indels = FALSE)` returns `5-10`, a deletion match. `vmatchPattern` with the same arguments returns no range.
  - `vmatchPattern2("ACGAACG", DNAStringSet(c(a = "TTTTACGAACGTTTT")), max.mismatch = 1, min.mismatch = 1)` returns `5-11`, an exact match. `vcountPattern` with the same arguments returns 0.
  - If the pattern is not longer than `max.mismatch`, the C code sends "indels" to the mismatch search, and this search applies `min.mismatch` (`_match_pattern_XString`, `src/match_pattern.c` of Biostrings 2.78.0). For example, `vmatchPattern2("AC", DNAStringSet(c(a = "AC")), max.mismatch = 2, min.mismatch = 1)` returns `0-1` and `2-3`, but not `1-2`. Thus `vmatchPattern2` ignores `with.indels` and `min.mismatch` only for a pattern longer than `max.mismatch`.
  - `SeedMatchR` has the default `with.indels = TRUE` (`R/SeedMatchR.R:78`). Its `iranges` mode passes this value to `vmatchPattern` (`R/SeedMatchR.R:128-133`). If `max.mismatch` is greater than 0, the `iranges` mode stops with the same error.
- New behavior: `find_matches` and `count_matches` obey search rules 2, 4, 5, and 6 of `docs/design.md` for the same `SearchOptions`. If `indels` is off, the ranges are the mismatch matches, and `min_mismatch` applies. If `indels` is on and `max_mismatch` is greater than 0, the ranges are the indel matches. `SearchOptions` raises `InputError` for `min_mismatch` above 0 with `indels` on.
- Reason: A caller who sets `with.indels` or `min.mismatch` expects the ranges to obey them, as the counts do. In R, the ranges and the counts of one call can disagree.
- Test: `test_ranges_obey_the_indels_flag_D001` and `test_ranges_obey_min_mismatch_D001` in `tests/search/test_api.py`.

## D-002: Each sequence of a gene union keeps its gene ID

- ID: `D-002`
- R location: `load_annotations` (`R/load_annotations.R:103-107`). `.reduce_features` splits the features by gene (`R/load_annotations.R:207-210`). `get_3UTR_granges` calls it when `reduce.features` is `TRUE` (`R/load_annotations.R:376-378`).
- R behavior: With `reduce.features = TRUE`, the names of the features are gene IDs. `extractTranscriptSeqs` gives these names to the sequences. `tx.gtf` keeps the transcripts whose `tx_id` is one of these names (`R/load_annotations.R:103`). No transcript ID is a gene ID, thus `tx.gtf` has no row. Then `match` finds no gene ID for the names (`R/load_annotations.R:105-107`), and each sequence gets the name `NA`. This evidence shows the behavior:
  - `tests/golden/annotation/utr3_union_r_sequences.tsv` comes from the unpatched R. `name_gene` is `NA` in each row, and `name_tx` holds the gene IDs.
- New behavior: A sequence of `PerGene.UNION` has its gene ID as its ID and as its gene ID. No sequence has the name `NA`.
- Reason: A result for a gene joins the sequences with the DE table through the gene ID. A sequence with the name `NA` matches no gene.
- Test: `test_union_sequences_keep_the_gene_id_D002` in `tests/annotation/test_reference.py`.

## D-003: Gene-level results use the maximum over the sequences of a gene

- ID: `D-003`
- R location: `.get_counts_column` (`R/SeedMatchR.R:207`) and `SeedMatchReport` (`R/SeedMatchReport.R:85-99`). `SeedMatchR` calls `.get_counts_column` in `DESEQ2` mode (`R/SeedMatchR.R:94-104`).
- R behavior: The counts of `SeedMatchR` use the first sequence of a gene, but `SeedMatchReport` counts each sequence. R probes with Biostrings 2.78.0 give the evidence:
  - `match(res$gene_id, match.df$gene_id)` finds the first sequence of a gene. Thus R ignores the other sequences of the gene.
  - `seqs` has two sequences with the name `g`, `CTCTATA` and `CTCTATACTCTATACTCTATA`. `SeedMatchR(seqs = seqs, sequence = guide, res = data.frame(gene_id = "g"), seed.name = "mer7m8")` gives the count 1, but the second sequence has 3 matches.
  - The test data has one sequence for each gene. Thus the R golden (`tests/golden/seedmatchr/deseq2_mode.tsv`) and the maximum agree there.
  - `SeedMatchReport` counts each sequence, with `n()` after `group_by(best.match)` (`R/SeedMatchReport.R:85-99`).
  - `seqs` has two sequences with the name `g`, both `CTCTATAA`. `res` has one row for `g`, with `log2FoldChange` -2 and `padj` 0.01. `SeedMatchReport(res, seqs, guide)` gives 2 for `8mer` in the rows `In silico predictions`, `Expressed predictions`, and `Off-target predictions`. The percentage is 100.
- New behavior: `annotate_seed_counts` gives each gene the maximum count over the sequences of the gene. A gene without a sequence gets 0. `seed_match_report` counts each gene once, in the class of its best match over the sequences of the gene. This best match is the best match of the maximum of each count column. The `matches` table keeps one row for each sequence.
- Reason: The count of a gene must not depend on the order of its sequences. In R, the first sequence in the input decides the count. The DE table has one row for each gene. Thus a gene with more than one sequence must count once in the report.
- Test: `test_gene_count_is_the_maximum_over_its_sequences_D003` and `test_the_maximum_does_not_depend_on_the_order_of_the_sequences_D003` in `tests/report/test_annotate.py`. `test_report_counts_genes_not_sequences_D003` and `test_the_best_match_of_a_gene_is_the_best_over_its_sequences_D003` in `tests/report/test_report.py`.

## D-009: A sequence has its transcript ID, and hg38-old is not a reference

- ID: `D-009`
- R location: `load_annotations` (`R/load_annotations.R:42`, `:116`) and `load_species_anno_db` (`R/load_annotations.R:144`, `:152-154`).
- R behavior: The R source and an R probe show these defects:
  - `hg38-old` is a choice of `reference.name` (`R/load_annotations.R:42`, `:144`). The table `ah.df` has no row for it, and thus no AnnotationHub ID (`R/load_annotations.R:152-154`).
  - `tx.names` is `features.seqs$names` (`R/load_annotations.R:116`). The probe calls `load_annotations("rnor7", feature.type = "3UTR", canonical = FALSE)` on the fixture with its local EnsDb and 2bit file. It uses SeedMatchR 2.0.0, Biostrings 2.78.0, and S4Vectors 0.48.1. `$seqs` is a `DNAStringSet`. `$tx.names` is `NULL`, and so is `$seqs$names`. `mcols($seqs)` has no column. Thus the list does not give the transcript of each sequence. The names of the sequences are gene IDs by default (`return_gene_name = TRUE`).
- New behavior: There is no `hg38-old`. `Reference.from_paths` and `Reference.from_directory` take the files of a reference. `Reference.sequences` gives the transcript ID of each sequence in `ids`, and the gene ID in `gene_ids`.
- Reason: A caller must know the transcript of each sequence, for example to join the sequences with `Reference.transcripts`. A choice without an ID cannot load a reference.
- Test: `test_sequence_ids_are_the_transcript_ids_D009` in `tests/annotation/test_reference.py`.

## D-010: kmer_gene_counts returns the counts of kmer_hist

- ID: `D-010`
- R location: `kmer_hist` (`R/kmer_hist.R:32-36`).
- R behavior: `plot(hist(...))` returns `NULL`, and `NULL + abline(...)` returns `integer(0)`. Thus `kmer_hist` returns `integer(0)` and no data. An R probe with R 4.5.3 gives the evidence: `kmer_hist(<the sequences of kmer_cases.tsv>, "ACG", width = 3)`, with `pdf(NULL)` as the device, returns `integer(0)`. The note of `tests/golden/analyses/kmer_gene_counts.tsv` records this value.
- New behavior: `kmer_gene_counts` returns the gene count of each k-mer. The plot module of a later plan draws the histogram.
- Reason: R returns no data, thus a caller cannot read the counts or the target count.
- Test: `test_kmer_gene_counts_give_the_counts_of_kmer_hist_D010` in `tests/test_kmers.py`.

## D-011: A sequence with one match gives no distance

- ID: `D-011`
- R location: `match_distances` (`R/match_distances.R:33-44`).
- R behavior: For one match, `1:(length(feature_coordinates) - 1)` is `c(1, 0)`.
  - For j = 1, `(j + 1):length(feature_coordinates)` is `c(2, 1)`. Thus the loop calculates the pairs (1, 2) and (1, 1), which give `NA` and 0.
  - For j = 0, the pair has no value, and it adds nothing.
  - The distances are `c(NA, 0)`. An element without a range gives `c(NA, NA)`.
  - For a name that occurs again, the list keeps the distances of the last element (`R/match_distances.R:44`).
  - Evidence: `tests/golden/analyses/distances.tsv` with its notes, and a probe of the planner.
- New behavior: A sequence with n matches gives each of the n(n − 1)/2 pair distances. A sequence with one match gives no row. The rows of one `seq_id` are one sequence.
- Reason: A sequence with one match has no pair. The value `NA` and the self-distance 0 change each summary of the distances.
- Test: `test_a_sequence_with_one_match_gives_no_distance_D011` and `test_a_table_of_one_match_gives_no_row_D011` in `tests/test_distances.py`.

## D-012: The best match uses only the columns that the search made

- ID: `D-012`
- R location: `.get_best_match` (`R/full_search.R:96-125`). `full_search` calls it for `res_format = "counts"` (`R/full_search.R:89-91`).
- R behavior: `case_when` names each column of the default search (`R/full_search.R:98-121`). It evaluates the condition of each column, thus a column that the search did not make stops the call. R probes with Biostrings 2.78.0 and dplyr 1.2.1 on the test data (`tests/data/`) give the evidence:
  - `full_search(guide, seqs, nonseed.edit.distances = 0:1)` stops with `object 'Full.mm2' not found`.
  - `full_search(guide, seqs, seqs.to.search = c("Full", "mer7m8"))` stops with `object '18mer.mm0' not found`.
  - `full_search(guide, seqs, seqs.to.search = c("mer7m8", "mer8"))` stops with `object 'Full.mm0' not found`.
  - `full_search(guide, seqs, nonseed.edit.distances = 0:5)` makes the column `Full.mm5`, but the ranking does not name it. Thus the column does not change the best match.
- New behavior: `best_match` skips each column of the ranking that the table lacks. A table that has none of these columns gives `No.matches` for each row. The ranking and the labels are the ones of R. As in R, a column that the ranking does not name does not change the best match. The default search gives the best match of R (`tests/golden/seedmatchr/full_search_counts_default.tsv`).
- Reason: A caller who selects seeds or edit distances for `full_search` expects the best match of the columns that the search made. A missing column is not an error.
- Test: `test_best_match_skips_absent_columns_D012`, `test_best_match_of_a_frame_without_count_columns_is_no_matches_D012`, and `test_best_match_of_a_frame_without_rows_is_empty_D012` in `tests/report/test_full_search.py`.

## D-015: A seed-count column replaces a column with the same name

- ID: `D-015`
- R location: `.get_counts_column` (`R/SeedMatchR.R:214-217`). `SeedMatchR` calls it in `DESEQ2` mode.
- R behavior: `cbind(res, seed.name = match.counts)` appends the counts as the last column, and the next statement gives this column the seed name. If `res` has a column with this name, the result has two columns with the same name. `res$mer7m8` then reads the old column. The test data (`tests/data/DESEQ2_results.tsv`) has a `mer7m8` column. R probe with Biostrings 2.78.0: `names(SeedMatchR(seqs = seqs, sequence = guide, res = res, seed.name = "mer7m8"))` ends with `mer7m8, mer7m8` (`tests/golden/seedmatchr/deseq2_mode.tsv` records this in its header). After `res$mer7m8 <- 99L`, `$mer7m8` of the result is 99 for each gene, and the last column holds the counts.
- New behavior: If the table has a column with this name, `annotate_seed_counts` replaces the values of this column. The column keeps its position. A new column goes at the end of the table.
- Reason: A table with two columns of the same name is ambiguous. In pandas, `de["mer7m8"]` then returns a table, not a column. A caller who reads the count column must get the new count.
- Test: `test_existing_column_is_replaced_not_duplicated_D015` and `test_replaced_column_keeps_its_position_and_gets_the_new_values_D015` in `tests/report/test_annotate.py`.

## D-016: A FASTA record ID ends at the first white space

- ID: `D-016`
- R location: `Biostrings::readDNAStringSet`. The tests of SeedMatchR call it (`tests/testthat/helper.R:10`). The package takes the result as the argument `seqs`. `.get_counts_column` makes `gene_id = names(seqs)` (`R/SeedMatchR.R:203`) and matches it with the `gene_id` of the DE table (`R/SeedMatchR.R:207`).
- R behavior: The name of a sequence is the complete header line after `>`. R probe with Biostrings 2.78.0 on a FASTA file with the header `>ENSRNOG00000016275 chromosome:mRatBN7.2:1:1:100:1 gene_biotype:protein_coding`: `names(readDNAStringSet(file))` is `"ENSRNOG00000016275 chromosome:mRatBN7.2:1:1:100:1 gene_biotype:protein_coding"`.
- New behavior: `read_fasta` takes the ID of a record from the header text up to the first white space. The description after the white space is not kept. A header that starts with white space gives an empty ID.
- Reason: An Ensembl FASTA header has a description after the ID, and the DE table holds the bare ID. A name with the description does not match the `gene_id` of the DE table.
- Test: `test_d016_read_fasta_takes_the_id_up_to_the_first_white_space` in `tests/sequence/test_sequence_set.py`.

## D-017: Gap and mask letters in a DNA sequence cause an error

- ID: `D-017`
- R location: `Biostrings::DNAString` and `Biostrings::readDNAStringSet`. SeedMatchR makes its patterns with `DNAString` (`R/vmatchPattern2.R:24`, `R/get_seed.R:172`), and its tests read the subjects with `readDNAStringSet` (`tests/testthat/helper.R:10`).
- R behavior: Both functions accept the gap letter `-`, the mask letter `+`, and the unknown letter `.`, and they keep these letters. R probe with Biostrings 2.78.0: `as.character(DNAString("AC-GT+A."))` is `"AC-GT+A."`, and `as.character(readDNAStringSet(file))` of a record with this sequence is `"AC-GT+A."`. `DNA_ALPHABET` lists these letters after the IUPAC letters.
- New behavior: `encode_dna`, `reverse_complement`, `SequenceSet.from_records`, and `read_fasta` raise `InputError` for each of these letters. The message gives the letter and its 1-based position. For a record, the message also gives the ID of the record.
- Reason: These letters are not bases. A search over them has no defined meaning in the port, and an error at the input boundary is better than a silent result.
- Test: `test_d017_gap_mask_and_unknown_letters_are_an_error` in `tests/sequence/test_alphabet.py` and `test_d017_read_fasta_rejects_gap_mask_and_unknown_letters` in `tests/sequence/test_sequence_set.py`.

## D-018: A lowercase guide gets the wobble like an uppercase guide

- ID: `D-018`
- R location: `.get_canonical_seed` and `.get_custom_sequence`. Each function calls `stringr::str_replace_all(guide.seq, "U", "Y")` when `allow_wobbles` is true (`R/get_seed.R:106-108`, `R/get_seed.R:160-162`).
- R behavior: The replacement matches only the uppercase `U`, but `RNAString` accepts lowercase letters. R probe with Biostrings 2.78.0: `get_seed("uuauagagcaagaacacuguuuu", "mer7m8", allow_wobbles = TRUE)` gives the target `CTCTATA`. The uppercase guide gives `CTCTRTR`. R gives no message.
- New behavior: The constructor of `Guide` changes the letters to uppercase. With `wobble=True`, a lowercase guide gives the same target as the uppercase guide.
- Reason: The case of a guide letter has no meaning. A wobble option that R ignores without a message gives a wrong result.
- Test: `test_d018_lowercase_guide_gets_the_wobble_target_of_the_uppercase_guide` in `tests/test_seed.py`.

## D-019: A guide holds only A, C, G, and U

- ID: `D-019`
- R location: `get_seed`, `.get_canonical_seed`, and `.get_custom_sequence`. Each function makes an `RNAString` of the guide (`R/get_seed.R:68`, `R/get_seed.R:111`, `R/get_seed.R:165`). No function checks the letters of the guide.
- R behavior: `RNAString` accepts each IUPAC letter, thus a guide can hold an ambiguity letter. R probe with Biostrings 2.78.0: `get_seed("UUNUAGAGCAAGAACACUGUUUU", "mer7m8")` gives the target `CTCTANA`. The letter `R` at g3 gives `CTCTAYA`. R rejects `T` (`key 84 (char 'T') not in lookup table`), and the port rejects it also.
- New behavior: `Guide` raises `InputError` for each letter that is not `A`, `C`, `G`, or `U` in either case. The message gives the letter and its 1-based position. The option `wobble` is the way to search with an ambiguous base.
- Reason: A guide is a synthesized RNA with defined bases. An ambiguity letter makes a pattern that matches more than the guide, and R gives no message. An error at the input boundary prevents such a pattern.
- Test: `test_d019_guide_rejects_each_ambiguity_letter` in `tests/test_seed.py`.

## D-020: A custom window with start after stop causes an error

- ID: `D-020`
- R location: `.get_custom_sequence`. The function checks only `stop.pos <= nchar(guide.seq)` and `start.pos > 0` (`R/get_seed.R:157-158`). Then it cuts the window with `substring` (`R/get_seed.R:168`).
- R behavior: No check compares `start.pos` with `stop.pos`, and `substring` gives an empty string. R probe with Biostrings 2.78.0: `get_seed("UUAUAGAGCAAGAACACUGUUUU", "x", start.pos = 8, stop.pos = 2)` gives `Target.Seq` `''`. The window `start.pos = 24, stop.pos = 23` gives `''` also.
- New behavior: `Guide.target` raises `InputError` for a `CustomSeed` with `start` after `stop`. The message gives `start` and `stop`.
- Reason: An empty pattern has no search meaning. R gives no message.
- Test: `test_d020_custom_window_start_must_not_pass_stop` in `tests/test_seed.py`.

## D-021: read_deseq2_table rejects a table that R reads with the wrong columns

- ID: `D-021`
- R location: `download_parse_file` (`R/get_data.R:74-76`). Line 74 reads the file with `read.csv`, and line 76 names the columns.
- R behavior: `colnames(res) <- c(...)` gives 8 names. R probes with R 4.5.3 show these results:
  - If the table has more columns, R does not stop. It gives the name `NA` to each extra column. A table with 9 columns gets the names `gene_id`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`, `symbol`, and `NA`. The test data (`tests/data/DESEQ2_results.tsv`) has a ninth column, `mer7m8`.
  - If the table has fewer columns, R stops with `'names' attribute [8] must be the same length as the vector [7]`.
  - If each data row has one field more than the header, `read.csv` uses the first field as the row name. Then each name labels the next column: `gene_id` holds the `baseMean` values. R gives no message.
  - A value that is not a number in a number column makes the column text. R gives no message, and `filter_res` stops later.
- New behavior: `read_deseq2_table` raises `InputError` for each of these tables. The message gives the path. For a table that does not have exactly 8 columns, the message also gives the count found and the 8 names. For a value that is not a number, it gives the column.
- Reason: In each case, the 8 names label the wrong values, or a later step stops with a message that does not name the file.
- Test: `test_read_deseq2_table_rejects_a_table_with_more_than_8_columns_D021`, `test_read_deseq2_table_rejects_a_row_with_more_fields_than_the_header_D021`, and `test_read_deseq2_table_names_the_column_with_a_value_that_is_not_a_number_D021` in `tests/test_de.py`.

## D-022: full_search_ranges does not depend on which searches find a match

- ID: `D-022`
- R location: `full_search_ranges` (`R/full_search.R:170-183`). `mapply` collects the result of each search, and `dplyr::bind_rows` joins them.
- R behavior: `mapply` simplifies a list of results that have the same length to a matrix. A search with a match gives a data frame of 9 columns, and a search without a match gives a data frame of 4 columns. If every result has the same length, `mapply` returns a matrix of lists, and `bind_rows` stops with `Argument 1 must be a data frame or a named atomic vector`. R probes with Biostrings 2.78.0 and the guide of the test data show these results, with the original `vmatchPattern2` and with the patched one (`tools/golden/patches/vmatchPattern2.R`):
  - A subject with a match for each search stops.
  - A subject without a match for any search stops.
  - A subject with a match for each search in one sequence and no match in another sequence stops.
  - The test data gives a result, because some searches have a match and some do not.
- New behavior: `full_search_ranges` joins the results of the searches itself. It gives the matches of each search. If no search finds a match, it gives a table without rows and with the same columns.
- Reason: The success of the call must not depend on the number of searches that find a match.
- Test: `test_full_search_ranges_with_a_match_for_every_search_D022` and `test_full_search_ranges_without_any_match_is_an_empty_frame_D022` in `tests/report/test_full_search.py`.

## D-023: The options reject duplicate values and an invalid max_seed_edit_distance

- ID: `D-023`
- R location: `full_search` (`R/full_search.R:42-66`) and `full_search_ranges` (`R/full_search.R:156-162`).
- R behavior: R probes with Biostrings 2.78.0 on the test data (`tests/data/`) give these results:
  - `full_search` with a duplicate seed (`seqs.to.search = c("Full", "Full")`) stops. It stops also with a duplicate value in `nonseed.edit.distances` or in `seed.edit.distances`. The message is `Can't transform a data frame with duplicate names`, and it does not name the argument.
  - `full_search_ranges` accepts a duplicate seed, a duplicate edit distance, and a duplicate indel option. It gives each match of the duplicate search again. For the first 60 sequences, `seqs.to.search = c("Full", "mer7m8", "mer7m8")` gives 18 rows for `mer7m8.mm0.indelFALSE`. One `mer7m8` gives 9 rows.
  - `max.seed.edit.dist = -1` removes each seed search, because `filter` compares `mm` with the number. `max.seed.edit.dist = 1.5` gives the same result as `1`.
- New behavior: `FullSearchOptions` and `FullSearchRangesOptions` raise `InputError` for a duplicate seed, a duplicate edit distance, and a duplicate indel option. `FullSearchRangesOptions` raises `InputError` for a `max_seed_edit_distance` that is not an integer of 0 or more. The message gives the field and the bad value.
- Reason: A duplicate gives columns with the same name or rows with the same content, and it is likely a typo. An edit distance is a whole number. To search none of the short seeds (`mer8`, `mer7m8`, `mer7A1`, `mer6`), do not put them in `seeds`.
- Test: `test_a_duplicate_value_is_rejected_D023` and `test_max_seed_edit_distance_must_be_an_integer_of_0_or_more_D023` in `tests/report/test_full_search.py`.

## D-024: The report rejects a NaN cutoff, a cutoff that is not a number, and a DE column that is not numeric

- ID: `D-024`
- R location: `SeedMatchReport` (`R/SeedMatchReport.R:94-99`). `dplyr::filter` compares `log2FoldChange` and `padj` with the cutoffs there.
- R behavior: R does not check the cutoffs or the type of the DE columns. R probes with Biostrings 2.78.0 on the test data (`tests/data/`) give these results:
  - The default cutoffs give 42 off-target predictions in `Total`, and the row `Expressed predictions` has 125.
  - `log2fc.cut = NaN` and `padj.cut = NaN` each give 0 off-target predictions. The row `Expressed predictions` does not change.
  - `log2fc.cut = TRUE` works as the number 1, and gives 51 off-target predictions.
  - `log2fc.cut = "-1"` compares the numbers with the text `-1` as text, and gives 19. The number -1 gives 28.
  - A `log2FoldChange` column of text gives 5 off-target predictions, and a `padj` column of text gives 0. R compares the text with the cutoff as text. For example, `"-2" <= -0.7655347` is `FALSE`, and `"-0.5" <= -0.7655347` is `TRUE`.
  - A `baseMean` column of text gives 42, because R only copies this column to `matches`.
  - R has no number that a double cannot hold: R reads `1e400` as `Inf`.
- New behavior: `ReportOptions` raises `InputError` for a NaN cutoff, for a cutoff that is a bool or text, and for an int that a float cannot hold. `seed_match_report` raises `InputError` that names the column and the dtype, if `log2FoldChange`, `padj`, or `baseMean` is not numeric. A bool column and a complex column are not numeric, as for `is.numeric` of R. `filter_de` and `ecdf_test` use the same check. The port does not compare `baseMean`, but `matches` gives it with the type `float64`. A Python int can be too large for a float, thus the port rejects it.
- Reason: A NaN cutoff gives no off-target, and a bool or text cutoff gives a wrong count. A text `log2FoldChange` or `padj` column also gives a wrong count. R shows none of these in an error, and a table that a caller read as text is a likely cause.
- Test: these tests:
  - `test_a_nan_cutoff_is_an_error_D024`, `test_a_cutoff_that_is_not_an_int_or_a_float_is_an_error_D024`, `test_a_cutoff_that_a_float_cannot_hold_is_an_error_D024`, and `test_a_column_that_is_not_numeric_is_an_error_that_names_the_column_and_the_dtype_D024` in `tests/report/test_report.py`
  - `test_a_bool_and_a_complex_column_are_not_numeric_D024` in `tests/report/test_report.py`, for `seed_match_report`
  - `test_filter_de_rejects_a_bool_and_a_complex_column_D024` in `tests/test_de.py`, for `filter_de`
  - `test_ecdf_test_rejects_a_bool_and_a_complex_fold_change_column` in `tests/stats/test_ecdf.py`, for `ecdf_test`
  - `test_check_numeric_columns_rejects_a_bool_column_a_complex_column_and_text` in `tests/test_checks.py`, for the helper that each function calls

## D-025: filter_de rejects a NaN cutoff and a cutoff that is not a number

- ID: `D-025`
- R location: `filter_res` (`R/filter_res.R:43-74`). `dplyr::filter` compares `padj`, `abs(log2FoldChange)`, and `baseMean` with the cutoffs at `R/filter_res.R:64`, `:67`, and `:70`.
- R behavior: R does not do a check of the cutoffs. An R probe with SeedMatchR 2.0.0 on a table of 2 genes gives these results for each of the 3 cutoffs:
  - A `NULL` cutoff stops with the dplyr message `..1 must be of size 2 or 1, not size 0`. The message gives the comparison, for example `padj <= fdr_cutoff`.
  - A `NaN` cutoff keeps 0 rows, and R gives no message.
  - A text cutoff compares the numbers as text. For example, `basemean_cutoff = "60"` keeps 0 rows, because the text `100` is less than the text `60`. The number 60 keeps 1 row.
  - `TRUE` works as the number 1.
- New behavior: `filter_de` raises `InputError` for a NaN cutoff, and for a cutoff that is a bool, text, or `None`. It also raises `InputError` for an int that a float cannot hold. The message gives the name of the cutoff and the bad value. The checks are the same as the checks of the report cutoffs (`D-024`).
- Reason: A NaN cutoff removes each row without a message, and a bool or text cutoff gives a wrong result. R shows none of these in an error. The report rejects the same values (`D-024`), thus the two functions agree.
- Test: `test_filter_de_rejects_a_cutoff_that_is_not_a_number_D025` and `test_filter_de_checks_each_cutoff_D025` in `tests/test_de.py`.

## D-026: The wobble with a fixed pattern causes an error

- ID: `D-026`
- R location: `SeedMatchReport` (`R/SeedMatchReport.R:51,55`), with the arguments `allow_wobbles` and `fixed = T`. `SeedMatchReport` gives `fixed` to `full_search` without a change (`R/SeedMatchReport.R:64`). The R documentation of `fixed` tells the user to set `fixed` to `FALSE` with wobbles (`R/SeedMatchR.R:37`).
- R behavior: The wobble target has `R` at each position of a guide `U`. With `fixed = TRUE`, that `R` matches only the letter `R` of a sequence, not `A` or `G`. A transcript has no `R`, thus the wobble finds fewer matches than the search without the wobble. R gives no message. `SeedMatchReport` on the test data (`tests/data/`) gives these counts of in-silico predictions in `Total`:
  - 125 without the wobble.
  - 42 with `allow_wobbles = TRUE`.
  - 161 with `allow_wobbles = TRUE` and `fixed = "subject"` (`tests/golden/seedmatchr/report_wobble.tsv`).
- New behavior: `ReportOptions` and `FullSearchOptions` raise `InputError` for `wobble=True` with `fixed_pattern=True`. The default of `ReportOptions` keeps `fixed_pattern=True`, as in R, thus `fixed_pattern=False` is necessary for a wobble report. With `options.fixed_pattern`, `annotate_seed_counts` raises `InputError` for a target with an ambiguity letter. The message gives the first such letter. `count_matches` and `find_matches` keep the Biostrings behavior.
- Reason: A wobble search with a fixed pattern loses the matches at each guide `U` without a message. The R documentation gives the rule, but R does not apply it. In `count_matches` and `find_matches`, a fixed ambiguity letter is a valid Biostrings search.
- Test: `test_the_wobble_with_a_fixed_pattern_causes_an_error_D026` and `test_the_wobble_adds_in_silico_predictions_D026` in `tests/report/test_report.py`, and `test_annotate_rejects_a_wobble_target_with_a_fixed_pattern_D026` in `tests/report/test_annotate.py`.

## D-027: The ranges of a gene union come in transcript order

- ID: `D-027`
- R location: `.reduce_features` (`R/load_annotations.R:183-213`). `get_3UTR_granges` calls it with `reduce.features = TRUE` (`R/load_annotations.R:376-378`), and `load_annotations` calls `get_3UTR_granges` (`R/load_annotations.R:84-87`).
- R behavior: `GenomicRanges::reduce` gives the ranges of each gene in ascending order of the position (`R/load_annotations.R:210`). Then `GenomicFeatures::extractTranscriptSeqs` joins the ranges in the order of the list. On the minus strand, the first base of the transcript is in the range with the highest position. Thus the R sequence of a minus-strand gene with more than one range is wrong. This evidence shows the behavior:
  - The rnor7 probe gives the ranges of a minus-strand gene in ascending order. The R sequence is equal to the join in list order, and it is not equal to the join in transcript order.
  - `tests/golden/annotation/utr3_union_r_features.tsv` comes from the unpatched R. It gives the ranges of the minus-strand gene `ENSSYNG00000000002` in ascending order.
  - `tests/golden/annotation/utr3_union_features.tsv` comes from R with the patch `tools/golden/patches/reduce_features.R`. The patch reverses the ranges of each minus-strand gene. The sequence of `ENSSYNG00000000002` in `utr3_union_r_sequences.tsv` is different from its sequence in `utr3_union_sequences.tsv`.
- New behavior: `PerGene.UNION` merges the ranges of the transcripts of a gene that overlap or touch. The merged ranges come in transcript order: in ascending order on the plus strand, and in descending order on the minus strand.
- Reason: A sequence must be the transcript sequence from 5' to 3'. The order of the ranges sets the order of the join, thus the R order gives a wrong sequence for a minus-strand gene.
- Test: `test_union_gives_the_minus_strand_ranges_in_transcript_order_D027` and `test_features_match_r` in `tests/annotation/test_features.py`.

## D-028: The longest transcript of a gene comes from the transcripts with the feature

- ID: `D-028`
- R location: `.get_longest_utr` (`R/load_annotations.R:406-418`). `get_3UTR_granges` calls it with `longest.utr = TRUE` (`R/load_annotations.R:370-372`).
- R behavior: `.get_longest_utr` selects the longest transcript of each gene from `GenomicFeatures::transcriptLengths` (`R/load_annotations.R:408-411`). This table holds each transcript of the filter, also a transcript without a CDS. If a gene has no coding transcript, the selected transcript is not in the 3' UTR list. Then the subscript at `R/load_annotations.R:415` stops with "subscript contains invalid names". This evidence shows the behavior:
  - In the rnor7 probe, `longest.utr = TRUE` with `protein.coding = FALSE` stops with "subscript contains invalid names".
  - The notes of `tests/golden/annotation/utr3_longest_noncoding_features.tsv` record the same stop of the unpatched R on the fixture. The fixture gene `ENSSYNG00000000003` has only a lncRNA transcript.
  - The golden comes from R with the patch `tools/golden/patches/get_longest_utr.R`. The patch keeps only the transcripts in the 3' UTR list.
- New behavior: `PerGene.LONGEST` selects only from the transcripts in the feature result. For the 3' UTR, the 5' UTR, and the CDS, these are the transcripts with a CDS. A gene without such a transcript has no feature.
- Reason: With `protein_coding=False`, a filter keeps the genes without a coding transcript, and R stops for each such filter. The patched R gives the result of the port.
- Test: `test_longest_skips_a_gene_without_a_coding_transcript_D028` and `test_features_match_r` in `tests/annotation/test_features.py`.

## D-029: PerGene applies to each feature type

- ID: `D-029`
- R location: `load_annotations` (`R/load_annotations.R:77-91`).
- R behavior: `load_annotations` gives `longest.utr` and `reduce.features` only to `get_3UTR_granges`. For `5UTR`, `exons`, and `cds`, R ignores the two options and gives no warning. This evidence shows the behavior:
  - The R source: the `switch` at `R/load_annotations.R:82-91` calls `ensembldb::fiveUTRsByTranscript`, `ensembldb::exonsBy`, and `ensembldb::cdsBy` without the two options. The message at `R/load_annotations.R:77-79` still says that R keeps the longest UTR of each gene.
  - The fixture golden `tests/golden/annotation/exons_longest_r_features.tsv` gives the exons with `longest.utr = TRUE`. It is equal to `tests/golden/annotation/exons_all_features.tsv`, the exons without the option.
- New behavior: `PerGene.LONGEST` and `PerGene.UNION` apply to each `FeatureType`. For the exons, `LONGEST` gives the transcript of each gene with the maximum sum of the exon widths.
- Reason: An option without an effect and without a warning hides a caller error. The rule for the longest transcript and for the union does not depend on the feature type.
- Test: these tests in `tests/annotation/test_features.py`:
  - `test_per_gene_applies_to_the_exons_D029`
  - `test_longest_and_union_agree_with_the_transcripts_for_each_feature_type_D029`
  - `test_union_merges_a_range_inside_an_earlier_range_and_keeps_the_genes_apart`

## D-030: The order of a tie of two genes does not depend on the filter

- ID: `D-030`
- R location: `ensembldb::transcripts` runs the SQL `order by gene.seq_name, tx.tx_seq_start`. `load_species_anno_db` calls it for `$gtf` (`R/load_annotations.R:160`), and `load_annotations` calls `load_species_anno_db` (`R/load_annotations.R:75`). `.get_longest_utr` uses it through `transcriptLengths` (`R/load_annotations.R:408`).
- R behavior: The SQL sorts on `seq_name` and `tx_seq_start` only. A tie keeps the scan order of the SQLite query plan, and the plan depends on the filter. The R goldens give the evidence. The fixture has the transcripts `ENSSYNT00000000010`, `ENSSYNT00000000002`, and `ENSSYNT00000000001`. They start at the same base of contig 1. `ENSSYNT00000000010` belongs to a gene whose row comes first in the `gene` table. The other two belong to a gene whose row comes later, but their `tx` rows come before the `tx` row of `ENSSYNT00000000010`:
  - `tests/golden/annotation/transcripts_exons_all.tsv` uses the filter on the standard chromosomes. The tie is `ENSSYNT00000000010`, `ENSSYNT00000000002`, `ENSSYNT00000000001`.
  - `tests/golden/annotation/transcripts_exons_everything.tsv` uses no filter on the chromosomes. The tie is `ENSSYNT00000000002`, `ENSSYNT00000000001`, `ENSSYNT00000000010`.
- New behavior: The tie is `ENSSYNT00000000010`, `ENSSYNT00000000002`, `ENSSYNT00000000001` for each filter. `AnnotationDb.transcripts` and `Reference.transcripts` sort on these keys, in this order:
  - `seq_name`
  - `tx_seq_start`
  - the row order of the `gene` table
  - the row order of the `tx` table
- Reason: The order of the transcripts sets the row order of the searches, thus it must not depend on the filter. An EnsDb from ensembldb has the `tx` rows of each gene in the order of the `gene` rows. Thus the two scan orders agree in such a file, and the port gives the order of R. The fixture has a tie where the two row orders disagree. This tie lets the test show the rule of the port.
- Test: `test_tie_order_does_not_depend_on_the_filter_D030` and `test_transcripts_match_r` in `tests/annotation/test_ensdb.py`.

## D-031: TranscriptFilter rejects a support level and an Entrez ID below 1

- ID: `D-031`
- R location: `build_annotation_filter` (`R/load_annotations.R:298-330`). `load_annotations` calls it (`R/load_annotations.R:63`).
- R behavior: R limits `transcript.support` only above 5, with `stopifnot(transcript.support <= 5)`. It does not limit `entrez_id`. An R probe gives these results:
  - `transcript.support = 0` and `transcript.support = -1` select no transcript. R gives no message.
  - `transcript.support = 6` stops with `transcript.support <= 5 is not TRUE`.
  - `entrez_id = -5` selects no transcript. R gives no message.
- New behavior: `TranscriptFilter` raises `InputError` for a `support_level` below 1 and for an Entrez ID below 1. The message gives the value and the permitted range. A `support_level` above 5 raises `InputError`, as in R.
- Reason: A value below 1 cannot match a transcript. An empty result without a message hides the mistake of the caller.
- Test: `test_support_level_below_1_raises_D031` and `test_entrez_id_below_1_raises_D031` in `tests/annotation/test_filters.py`.

## D-032: The asymptotic KS p-values are accurate

- ID: `D-032`
- R location: `ecdf_stat_test` calls `stats::ks.test` (`R/deseq_fc_ecdf.R:255-257`). In R 4.5.3, `ks.test` calls the R functions `psmirnov` and `psmirnov_asymp`. The two-sided limit distribution is the C function `K2l` in `src/library/stats/src/ks.c`.
- R behavior: R uses the asymptotic p-value if `n_x * n_y >= 10000`. Let en = `n_x * n_y / (n_x + n_y)`. The asymptotic p-value of R has these defects. The R goldens in `tests/golden/stats/ks.tsv` give the evidence:
  - For `sqrt(en) * D` below 1, `K2l` keeps one term of its series, because R gives `tol = 1e-06` (`ks.c:100-111`). The p-value has an error of about 1e-5. `ks_asym_cut_series` `two.sided` (120 x 150) gives 0.3243418012381829. The limit distribution gives 0.32432903587229633.
  - `psmirnov_asymp` gives the one-sided upper tail as `1 - (-expm1(-2 * en * D^2))`. The absolute error of this value is about 1e-16, thus its relative error is above 1e-9 for a p-value below about 1e-7. An R probe with `greater` on a sample of 150 values and a sample of 120 values gives 2.6115202222953826e-09, where `exp(-2 * en * D^2)` gives 2.6115201992964189e-09 (a relative difference of 8.8e-9). Below about 1e-16, the value is 0. `ks_asym_small_p` `greater` (400 x 300), `ks_120x150_shift` `greater`, and `ks_400x300_shift` `greater` give 0. The values of `exp(-2 * en * D^2)` are 2.5313259744298458e-30, 1.2418498224782002e-58, and 6.959995697096148e-149.
  - R finds D with a cumulative sum of 1/m and -1/n. Thus a rounding error can put D above 1. `psmirnov` gives the p-value 0 for D above 1. `ks_100x100_shift` has D = 1.0000000000000007, and `two.sided` and `greater` give 0. The limit values are 7.440151952041566e-44 and 3.720075976020836e-44.
- New behavior: The statistic comes from integer counts. A checkpoint is the end of a run of equal pooled values. Let i and j be the counts of x values and y values up to a checkpoint. The gap of the two ECDFs, scaled by `m * n`, is the integer `i * n - j * m`. Let d be the maximum of the gap for `GREATER`, of the negative gap for `LESS`, and of the absolute gap for `TWO_SIDED`. D is d divided by `m * n`, and d is at most `m * n`. Thus D is never above 1. The two-sided p-value is `scipy.special.kolmogorov(sqrt(en) * D)`. The one-sided p-value is `exp(-2 * en * D^2)`.
- Reason: R gives 0 for a p-value that a float64 can hold, for example 1e-43. A p-value of 0 hides the size of an effect. Below about 1e-7, the one-sided p-value of R also loses relative precision. The cut series changes the fifth digit of a p-value. Each defect is a numeric error of R, not a property of the test.
- Test: `test_ks_asymptotic_two_sided_uses_the_limit_distribution_D032`, `test_ks_asymptotic_one_sided_keeps_a_small_p_value_D032`, `test_ks_keeps_a_p_value_that_r_rounds_to_zero_D032`, and `test_ks_statistic_cannot_exceed_1_where_the_sum_of_r_does_D032` in `tests/stats/test_ks.py`. `test_ks_matches_r` in the same file covers each other row of `ks.tsv`.

## D-033: The Wilcoxon test gives the p-value 1 when each value is equal

- ID: `D-033`
- R location: `ecdf_stat_test` calls `stats::wilcox.test` (`R/deseq_fc_ecdf.R:274-276`). In R 4.5.3, `wilcox.test.default` uses the normal approximation if two values are equal.
- R behavior: If each pooled value is equal, each rank is the average rank, and z is 0 before the continuity correction. The variance SIGMA^2 is 0, thus R divides by 0. The golden `wilcoxon_all_equal` in `tests/golden/stats/wilcoxon.tsv` (5 values and 7 values, each equal to 1) gives the evidence:
  - `two.sided`: the correction `sign(z) * 0.5` is 0, thus z is 0/0, and the p-value is NaN.
  - `less`: z is 0.5/0 = Inf, and the p-value is 1.
  - `greater`: z is -0.5/0 = -Inf, and the p-value is 1.
- New behavior: `wilcoxon_test` gives the p-value 1 for each alternative if each pooled value is equal. The statistic is W, as in R, and the method is `ASYMPTOTIC`. The results for `LESS` and `GREATER` do not change. Only the two-sided result differs from R.
- Reason: Two samples of equal values show no difference, thus 1 is the correct p-value. A NaN p-value fails each comparison with a threshold. R already gives 1 for the one-sided tests.
- Test: `test_wilcoxon_with_all_values_equal_gives_1_D033` in `tests/stats/test_wilcoxon.py`. `test_wilcoxon_matches_r` in the same file covers each other row of `wilcoxon.tsv`.

## D-034: The DTS test rejects an infinite value, an empty sample, and a range that overflows

- ID: `D-034`
- R location: `ecdf_stat_test` calls `twosamples::dts_test` (`R/deseq_fc_ecdf.R:259-264`). The version of twosamples is 2.0.1.
- R behavior: `dts_test` does not check its samples. Let `a` be `c(0.5, -0.2, 1.1, 0.3, -0.7)` and `b` be `c(0.1, -0.4, 0.9, 0.2, -1, 0.6, -0.3)`. R probes with R 4.5.3 and `set.seed(1)` give the evidence:
  - `dts_test(c(a, Inf), b)`, `dts_test(c(a, -Inf), b)`, and `dts_test(a, c(b, Inf))` give the statistic `Inf` and the p-value 1.
  - `dts_test(numeric(0), b)` gives the statistic 3.9986665900519656 and the p-value 1. `dts_test(a, numeric(0))` gives 3.3254468626107201 and the p-value 1.
  - A NaN or an NA makes the function stop with the message `missing value where TRUE/FALSE needed`. The port also raises an error, and its message gives the sample and the index.
  - Finite values near the limit of float64 overflow. With `set.seed(1)`, `dts_test(c(-1e308, rnorm(5)), c(1e308, rnorm(5)))` gives the statistic `Inf` and the p-value 1. The same holds for `c(-8e307, rnorm(5))` and `c(8e307, rnorm(5))`, with the range 1.6e308.
- New behavior: `dts_test` raises `InputError` for an empty sample and for an infinite value, `Inf` or `-Inf`. The message gives the sample and the index of the value. It also raises `InputError` for a pooled range above the largest float divided by 2 n, where n is the pooled size. A range that overflows to `Inf` is above this limit. The message is `x and y span from -1e+308 to 1e+308; expected a range of at most` and the limit. `ecdf_test` gives the names `target` and `background`.
- Reason: An infinite value makes a statistic that does not measure a difference, and the p-value 1 states that the samples agree. An empty sample has no distribution to compare. R gives a result and no message, thus the mistake of the caller stays hidden. `ks_test` and `wilcoxon_test` accept an infinite value, because they use only the order of the values. Each term of the statistic is at most the width of its gap times n, thus the statistic is at most the range times n. The limit keeps the statistic at half of the largest float. The port multiplies the width by a factor of at most n / (n - 1). It does not multiply by n^1.5 first. Thus no product overflows below the limit.
- Test: `test_dts_rejects_an_infinite_value_D034`, `test_dts_rejects_an_empty_sample_D034`, `test_dts_rejects_a_range_that_overflows_D034`, `test_dts_names_the_limit_of_the_range_as_the_largest_float_over_2_n_D034`, and `test_dts_is_finite_and_scales_with_the_values_up_to_the_limit_D034` in `tests/stats/test_dts.py`. `test_ecdf_test_rejects_a_range_that_a_float_cannot_hold_with_the_names_of_the_sets` in `tests/stats/test_ecdf.py` covers the names `target` and `background`.

## D-035: ecdf_test rejects an alternative for DTS and Wasserstein

- ID: `D-035`
- R location: `ecdf_stat_test` passes `alternative` only to `stats::ks.test` (`R/deseq_fc_ecdf.R:255-257`) and to `stats::wilcox.test` (`R/deseq_fc_ecdf.R:274-276`). The branches of `twosamples::dts_test` (`R/deseq_fc_ecdf.R:259-264`) and `wass_dist` (`R/deseq_fc_ecdf.R:266-271`) do not use it.
- R behavior: R ignores `alternative` for `DTS` and `Wass`, and it gives no message. An R probe with R 4.5.3 and `set.seed(1)` on `tests/data/DESEQ2_results.tsv` gives the evidence. The target is the genes with `mer7m8 >= 1`, and the background is the genes with `mer7m8 == 0`:
  - `stats_test = "DTS"` gives the statistic 74.911568929640467 and the p-value 0.00025 for `greater`, `less`, and `two.sided`.
  - `stats_test = "Wass"` gives the statistic 2.1960000000000002 for `greater`, `less`, and `two.sided`.
  - `stats_test = "KS"` gives a different result for each `alternative`. The statistic is 0.51035196687370576 for `greater` and 1.1015494072452725e-15 for `less`.
- New behavior: `ecdf_test` raises `InputError` for an `alternative` that is not `None`, if the test is `DTS` or `WASSERSTEIN`. The message names the value and the test. The result has `alternative=None` for `DTS` and `WASSERSTEIN`.
- Reason: An argument that has no effect, and no message, hides a mistake of the caller. A caller who asks for `LESS` with DTS gets a result that does not depend on `LESS`.
- Test: `test_ecdf_test_rejects_an_alternative_for_dts_and_wasserstein_D035` in `tests/stats/test_ecdf.py`. `test_ecdf_test_matches_r` in the same file covers the rows of `ecdf.tsv` for `DTS` and `WASSERSTEIN`, which have no alternative.

## D-036: ecdf_test drops the rows without a fold change for each test

- ID: `D-036`
- R location: `ecdf_stat_test` selects the samples with `res$log2FoldChange[res$gene_id %in% unlist(target_list)]` (`R/deseq_fc_ecdf.R:247-252`). It then calls the test (`R/deseq_fc_ecdf.R:254-277`).
- R behavior: The samples keep the NA values, and each test reacts in its own way. An R probe with R 4.5.3 gives the evidence. The table is `tests/data/DESEQ2_results.tsv` with 5 more rows whose `log2FoldChange` is NA: 2 rows for the target and 3 rows for the background:
  - `stats::ks.test` and `stats::wilcox.test` drop the NA values. The result is equal to the result for the table without the 5 rows. The KS statistic is 0.51035196687370576, and the p-value is 0.00010377312440832125. The Wilcoxon statistic is 748, and the p-value is 0.99961672019886205.
  - `twosamples::dts_test` stops with the message `missing value where TRUE/FALSE needed`.
  - `wass_dist` stops with the message `'from' must be a finite number`.
  - Direct calls with NA values give the same reactions. For example, `ks.test(c(x, NA), c(y, NA, NA))` is equal to `ks.test(x, y)`.
  - A sample without a value after the drop makes `ks.test` stop with `not enough 'x' data`, and `wilcox.test` stop with `not enough (non-missing) 'x' observations`.
- New behavior: `ecdf_test` drops each row whose `log2FoldChange` is NaN, for each test. The fields `n_target` and `n_background` count the rows that remain. A set with no row after the drop raises `InputError` that names the set and gives the count 0. The functions `ks_test`, `wilcoxon_test`, `dts_test`, and `wasserstein_areas` drop no value. Each of them raises `InputError` for NaN, and the message names the sample and the index.
- Reason: A row without a fold change has no value to compare, thus the result must not depend on the test. A function that takes values must not drop a value in silence, because the number of values changes the p-value. The drop belongs to the selection of the rows of a table, in `ecdf_test`.
- Test: `test_ecdf_test_drops_rows_without_a_fold_change_D036` and `test_ecdf_test_rejects_a_set_whose_rows_have_no_fold_change` in `tests/stats/test_ecdf.py`. `test_the_direct_tests_raise_for_nan_where_ecdf_test_drops_the_row_D036` in the same file covers `ks_test`, `wilcoxon_test`, `dts_test`, and `wasserstein_areas`. Each module of a direct function also has its own test: `test_ks_rejects_nan`, `test_wilcoxon_rejects_nan`, `test_dts_rejects_nan`, and `test_wasserstein_rejects_nan_and_an_infinite_value`.

## D-037: The statistical tests give no warning for ties

- ID: `D-037`
- R location: `ecdf_stat_test` calls `stats::ks.test` (`R/deseq_fc_ecdf.R:255-257`) and `stats::wilcox.test` (`R/deseq_fc_ecdf.R:274-276`). `wass_dist` calls `stats::approx` for each sample (`R/wass_dist.R:70-71`).
- R behavior: R warns if ties make a p-value approximate. `approx` warns if it merges equal x values. R probes with R 4.5.3 give the evidence:
  - `ks.test(rep(c(1, 2, 3), 40), rep(c(2, 3, 4), 50))` uses the asymptotic method and gives the warning `p-value will be approximate in the presence of ties`. The exact method gives no warning, for example `ks.test(c(1, 2, 2, 3), c(2, 3, 4))`.
  - `wilcox.test(c(1, 2, 3, 3), c(3, 4, 5))` gives the warning `cannot compute exact p-value with ties`.
  - Let `x` be `c(1, 2, 2, 3, 4, 5, 6, 7, 8, 9)`. Then `approx(x, ecdf(x)(x), xout = c(1.5, 2.5), yleft = 0, yright = 1)` gives the warning `collapsing to unique 'x' values`. `wass_dist` gives this warning for each sample that has equal values.
- New behavior: `ks_test`, `wilcoxon_test`, and `wasserstein_areas` give no warning. The field `method` of a test result is `EXACT` or `ASYMPTOTIC`, thus the caller can read how the p-value was made. `wasserstein_areas` merges the equal values, as `approx` does.
- Reason: A caller that does many tests gets one warning for each call with ties. The result already holds this information, in a field that a program can read. The areas use the merged values, as in R. Thus the warning does not say that the areas are approximate.
- Test: `test_ks_gives_no_warning_for_ties_D037` in `tests/stats/test_ks.py`, `test_wilcoxon_gives_no_warning_for_ties_D037` in `tests/stats/test_wilcoxon.py`, and `test_wasserstein_gives_no_warning_for_ties_D037` in `tests/stats/test_wasserstein.py`.

## D-038: The Wasserstein areas reject a grid of fewer than 2 points and a range that overflows

- ID: `D-038`
- R location: `wass_dist` builds the grid with `seq` (`R/wass_dist.R:69`). It integrates with `caTools::trapz` (`R/wass_dist.R:79-81`).
- R behavior: R probes with R 4.5.3 give the evidence. Let `a` be `c(0.5, -0.2, 1.1, 0.3, -0.7, 0.9, -0.4, 0.1, 0.6, -0.9)`. Let `b` be `c(0.1, -0.4, 0.9, 0.2, -1, 0.6, -0.3, 0.05, 0.45, -0.15, 0.8, -0.6)`.
  - `wass_dist(a, b, n.interp = 1)` gives `NA` for each area, and no message.
  - `wass_dist(a, b, n.interp = 0)` stops with the message `only 0's may be mixed with negative subscripts`.
  - `wass_dist(a, b, n.interp = -5)` stops with the message `'length.out' must be a non-negative number`.
  - `wass_dist(c(-1e308, a), c(1e308, b))` gives the areas 1.74825173985846e+307, 1.74825173985846e+307, and 0, and no message. The width of the range is `Inf` in float64.
  - `wass_dist(seq(-6e307, -6.6e307, length.out = 10), seq(6e307, 6.6e307, length.out = 10), n.interp = 100)` gives the areas `Inf`, `Inf`, and 0, and no message. The range is 1.32e308, which float64 can hold, but the sum of `trapz` before the division by 2 is about 2.4e308.
- New behavior: `wasserstein_areas` raises `InputError` for an `n_interp` below 2, with the message `n_interp 1 is less than 2; expected 2 or more`. It raises `InputError` if the range of the pooled values is above half of the largest float, 8.988465674311579e+307. The message is `seed and non_seed span from -1e+308 to 1e+308; expected a range of at most 8.988465674311579e+307`. `ecdf_test` gives the names `target` and `background`.
- Reason: An `NA` area or a number from an overflow hides the mistake of the caller. A grid of 1 point has no trapezoid, and the area of an overflowed range has no meaning. Each product of the trapezoid rule is a width times a sum of two heights from 0 to 1. Thus the sum is at most twice the range, and a range above half of the largest float can overflow it. The messages of R for 0 and for a negative value do not name `n.interp`.
- Test: `test_wasserstein_rejects_n_interp_below_2_D038`, `test_wasserstein_rejects_a_range_that_a_float_cannot_hold_D038`, and `test_wasserstein_rejects_a_finite_range_whose_sum_overflows_D038` in `tests/stats/test_wasserstein.py`. `test_wasserstein_accepts_a_range_just_below_the_limit_D038` in the same file shows that the limit is not too low. `test_ecdf_test_rejects_a_range_that_a_float_cannot_hold_with_the_names_of_the_sets` in `tests/stats/test_ecdf.py` covers the names `target` and `background`.

## D-039: assign_groups warns only if two groups share a gene

- ID: `D-039`
- R location: `.check_gene_list_overlap` (`R/deseq_fc_ecdf.R:165-169`). `deseq_fc_ecdf` calls it (`R/deseq_fc_ecdf.R:97`) before `.assign_groups` (`R/deseq_fc_ecdf.R:171-182`).
- R behavior: R warns `Gene lists overlap!` if `any(duplicated(unlist(gene_lists)))`. Thus R also warns for a gene that occurs two times in one list. An R probe with R 4.5.3 gives the evidence:
  - `list(a = c("g1", "g1"), b = c("g2"))` gives the warning `Gene lists overlap!`, but no two lists share a gene.
  - `list(a = "g1", b = "g1")` gives the same warning, and here two lists share a gene.
  - `.assign_groups` labels the genes of `list(a = c("g1", "g1", "g2"), b = c("g2", "g3"))` as `a` for `g1` and `b` for `g2` and `g3`. The last list wins. The second `g1` in `a` does not change a label.
- New behavior: `assign_groups` warns with `SeedMatchWarning` only if two groups share a gene. The message names the first two groups that share a gene, and one shared gene. A gene that occurs two times in one group gives no warning. The last group wins, as in R.
- Reason: The warning tells the caller that the label of a gene depends on the order of the groups. A gene that occurs two times in one group keeps its label, thus the R warning is false for it. The message names the groups and a gene, thus the caller can find the overlap.
- Test: `test_assign_groups_gives_no_warning_for_a_repeat_inside_one_group_D039` and `test_assign_groups_gives_the_last_group_and_warns_for_an_overlap` in `tests/stats/test_ecdf.py`.

## D-040: assign_groups keeps a group named NA apart from the genes in no group

- ID: `D-040`
- R location: `.assign_groups` (`R/deseq_fc_ecdf.R:171-182`) and the `factor` call of `deseq_fc_ecdf` (`R/deseq_fc_ecdf.R:109`).
- R behavior: `.assign_groups` writes the text `"NA"` as the group of each gene in no group. The `factor` call turns this text into a level only if a group has the name `NA`. Then this level holds the genes of that group and each gene in no group. An R probe with R 4.5.3 gives the evidence. The table has the genes `g1`, `g2`, `g3`, `g4`, `g2`, `g5`, and `g6` in this order. The groups are `list(a = c("g1"), "NA" = c("g2"))`:
  - The levels are `a` and `NA`.
  - The label of `g1` is `a`. The label of each other row is `NA`, also for `g3`, `g4`, `g5`, and `g6`, which are in no group.
  - Without a group named `NA`, the same rows have no value, as in the probe `list(a = c("g1"), b = c("g2"))`.
- New behavior: `assign_groups` gives no value to a gene in no group, also if a group has the name `NA`. A group named `NA` holds only its own genes.
- Reason: A gene in no group is not a member of a group named `NA`. The mix gives wrong group sizes and a wrong comparison of the groups. The name of a group is data of the caller, thus it must not change the meaning of a gene in no group.
- Test: `test_assign_groups_keeps_a_group_named_na_apart_D040` in `tests/stats/test_ecdf.py`.

## D-041: the registry pins the STRING version

- ID: `D-041`
- R location: `get_stringPPI_to_stringGGI` (`R/get_stringPPI_to_stringGGI.R:15`) and `STRINGdb$new` (STRINGdb 2.22.0, `R/rstring.R:145-151`).
- R behavior: `get_stringPPI_to_stringGGI` calls `STRINGdb::STRINGdb$new(species = sp_number)` without a `version`. In STRINGdb 2.22.0 (Bioconductor 3.22), `STRINGdb$new` then asks `https://string-db.org/api/tsv-no-header/version` and uses the answer. On 2026-11-20, the answer is 12.5. Thus the network of R depends on the date of the run. Evidence: the source of STRINGdb 2.22.0 and the answer of the API on 2026-11-20. No R run.
- New behavior: The registry pins STRING 12.5. `seedmatch-fetch string TAXON DEST` downloads `<TAXON>.protein.links.v12.5.txt.gz` and makes sure that the size and the SHA-256 digest agree with the registry.
- Reason: A pinned version and its digest give the same file at each run. A new STRING release does not change a result, because the registry must change first.
- Test: `test_the_string_urls_name_the_pinned_version_D041` in `tests/resources/test_registry.py`.

## D-042: kmer_gene_counts counts genes, not sequences

- ID: `D-042`
- R location: `kmer_hist` (`R/kmer_hist.R:24-28`).
- R behavior: `colSums` of the presence matrix counts the sequences with the k-mer. The axis label is `Number of genes with k-mer` (`R/kmer_hist.R:32`), and the title is `K-mer gene frequencies`. A gene with more than one sequence counts more than once. Evidence: the case `cases_k2_s1` of `tests/golden/analyses/kmer_gene_counts.tsv`. The 2-mer `AC` is in the sequences `s1`, `s2`, `s3`, and `s6`, which belong to the genes `g1` and `g2`. R gives 4, and the port gives 2.
- New behavior: `kmer_gene_counts` counts the distinct gene IDs. A set with one sequence for each gene gives the R values.
- Reason: The label of R names genes. A result for a gene uses all the sequences of the gene (`D-003`).
- Test: `test_kmer_gene_counts_count_genes_not_sequences_D042` in `tests/test_kmers.py`. `test_kmer_gene_counts_match_r_for_one_sequence_per_gene` in the same file shows that the R values agree if each gene has one sequence.

## D-043: kmer_counts rejects a k above 12

- ID: `D-043`
- R location: `sequence_kmer_counts` calls `oligonucleotideFrequency` (`R/kmer_hist.R:63-66`).
- R behavior: R accepts `width = 13`. It makes a dense integer matrix with 4^13 = 67,108,864 columns for each sequence. Evidence: the note of `tests/golden/analyses/kmer_counts.tsv`. `dim` of `oligonucleotideFrequency(DNAStringSet("ACGTACGTACGTAC"), width = 13)` is 1 x 67108864.
- New behavior: `kmer_counts` raises `InputError` for a `k` above 12, with the message `k 13 is out of range; expected an integer from 1 to 12`. A `k` below 1 raises the same message. The matrix is sparse.
- Reason: The labels of `kmers()`, `to_frame()`, and `kmer_gene_counts` grow as 4^k. 4^12 is 16,777,216 labels. `kmers()` holds about 1 GB, and `kmer_gene_counts` takes about 3 s and peaks near 2 GB.
- Test: `test_kmer_counts_reject_k_above_12_D043` in `tests/test_kmers.py`. `test_kmer_counts_accept_k_12` in the same file covers the limit.

## D-044: kmer_counts gives a row for each sequence, also for a repeated ID

- ID: `D-044`
- R location: `sequence_kmer_counts` (`R/kmer_hist.R:71`).
- R behavior: `rownames(olig.freq) = names(seqs)` stops with `duplicate 'row.names' are not allowed`. Evidence: the note of `tests/golden/analyses/kmer_counts.tsv`.
- New behavior: `kmer_counts` gives one row for each sequence, in the order of the set.
- Reason: A row of the matrix is a sequence, thus a repeated ID does not make a count ambiguous.
- Test: `test_kmer_counts_keep_each_sequence_with_a_repeated_id_D044` in `tests/test_kmers.py`.

## D-045: alignment_text shows the full guide

- ID: `D-045`
- R location: `plotAlignment` (`R/plotAlignment.R:28-57`), `get_seed` (`R/get_seed.R:121`), and `.matchStrings` (`R/plotAlignment.R:84`).
- R behavior: `Full` is a canonical seed, thus `get_seed(guide.seq, "Full", start.pos = 1, stop.pos = guide.len)` ignores the positions. Its target is g2 to gN, with N - 1 letters, but the window has N letters.
  - Each call with the default positions stops with `nchar(s1) == nchar(s2) is not TRUE`. Evidence: the notes of `tests/golden/analyses/alignment.tsv`.
  - With `seq.start.pos = 1`, the call works, but the guide line omits g1.
  - A row of `SeedMatchR(res.format = "data.frame")` has no column `names`, thus R stops with `subscript out of bounds`. Evidence: the notes of `tests/golden/analyses/alignment.tsv`.
- New behavior: The guide line holds g1 to gN, and `match` gives the sequence by `seq_id`.
- Reason: The R documentation calls `guide.seq` "The full siRNA guide sequence". `tests/testthat/helper.R:27-45` of SeedMatchR uses strings of 23 letters for the guide of 23 letters. The function must work for its documented input.
- Test: `test_alignment_text_shows_the_full_guide_D045` in `tests/test_alignment.py`. `test_alignment_text_matches_r` in the same file compares the text with the R text of the patched function (`tools/golden/patches/plotAlignment.R`).

## D-046: alignment_text shows a window past a sequence end

- ID: `D-046`
- R location: `plotAlignment` (`R/plotAlignment.R:43-45`, `XVector::subseq`).
- R behavior: A window that starts before 1 stops with `solving row 1: the supplied start/end lead to a negative width`. A window that starts at 0 stops with `solving row 1: 'allow.nonnarrowing' is FALSE and the supplied start (0) is < 1`. A window that ends after the sequence stops with `solving row 1: 'allow.nonnarrowing' is FALSE and the supplied end (147) is > refwidth`. Evidence: the notes of `tests/golden/analyses/alignment.tsv`.
- New behavior: The target line and the pair line have a space for each position outside the sequence. Each pad must be at most the length of the sequence of the match, and a longer pad raises `InputError`. A pad above the length shows no base, and R stops for such a window too. A short sequence, also an empty sequence, accepts only a pad up to its length.
- Reason: A seed match in the first N - b bases of a 3' UTR has no alignment in R. Here N is the length of the guide, and b is the last guide position of the seed window. `find_matches` also gives mismatch matches that start before 1 (search rule 2). The bound limits the window to N + 2L bases for a sequence of length L. It never rejects a call that R accepts.
- Test: `test_alignment_text_puts_spaces_outside_the_sequence_D046`, `test_alignment_text_accepts_a_mismatch_match_before_the_start_D046`, `test_alignment_text_accepts_a_mismatch_match_after_the_end_D046`, and `test_alignment_text_rejects_a_pad_longer_than_the_sequence_D046` in `tests/test_alignment.py`.
