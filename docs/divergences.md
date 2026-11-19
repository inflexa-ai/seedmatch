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
- New behavior: `ReportOptions` raises `InputError` for a NaN cutoff, for a cutoff that is a bool or text, and for an int that a float cannot hold. `seed_match_report` raises `InputError` that names the column and the dtype, if `log2FoldChange`, `padj`, or `baseMean` is not numeric. The port does not compare `baseMean`, but `matches` gives it with the type `float64`. A Python int can be too large for a float, thus the port rejects it.
- Reason: A NaN cutoff gives no off-target, and a bool or text cutoff gives a wrong count. A text `log2FoldChange` or `padj` column also gives a wrong count. R shows none of these in an error, and a table that a caller read as text is a likely cause.
- Test: `test_a_nan_cutoff_is_an_error_D024`, `test_a_cutoff_that_is_not_an_int_or_a_float_is_an_error_D024`, `test_a_cutoff_that_a_float_cannot_hold_is_an_error_D024`, and `test_a_column_that_is_not_numeric_is_an_error_that_names_the_column_and_the_dtype_D024` in `tests/report/test_report.py`.

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
