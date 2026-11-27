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
- New behavior: `annotate_seed_counts` and `annotate_wobble_bulge_counts` give each gene the maximum count over the sequences of the gene. A gene without a sequence gets 0. `seed_match_report` counts each gene once, in the class of its best match over the sequences of the gene. This best match is the best match of the maximum of each count column. The `matches` table keeps one row for each sequence.
- Reason: The count of a gene must not depend on the order of its sequences. In R, the first sequence in the input decides the count. The DE table has one row for each gene. Thus a gene with more than one sequence must count once in the report.
- Test: `test_gene_count_is_the_maximum_over_its_sequences_D003` and `test_the_maximum_does_not_depend_on_the_order_of_the_sequences_D003` in `tests/report/test_annotate.py`. `test_report_counts_genes_not_sequences_D003` and `test_the_best_match_of_a_gene_is_the_best_over_its_sequences_D003` in `tests/report/test_report.py`. `test_annotate_wobble_bulge_counts_give_the_gene_maximum_D003` and `test_wobble_bulge_null_takes_the_gene_maximum_of_the_sequences_D003` in `tests/noncanonical/test_wobble.py`.

## D-004: The 3' compensatory search runs and returns the sites

- ID: `D-004`
- R location: `get_compensatory_sequences` (`R/get_compensatory_sequences.R:57-77`, `:34`, `:99`).
- R behavior: The frame of `:57-66` has no row, thus the test of `:67` is true. The function returns 0 rows with the message `No compensatory regions found.` for each input (note `original_vignette` of `tests/golden/noncanonical/compensatory.tsv`).
  - Without the early return, `:99` stops at the first site. `gsub` at `:34` makes the seed a `character` value, and `reverseComplement` has no method for it (note `line_99`).
- New behavior: `compensatory_sites` runs the search. Each row is a pair of a seed match and a region. The golden of the patched function (`tools/golden/patches/get_compensatory_sequences.R`) gives the R rows of each case. With `first_position=8`, the port gives the R rows of the cases on `Test_seqs.fa`.
- Reason: R gives no result for each input.
- Test: `test_compensatory_sites_match_r_D004` in `tests/noncanonical/test_compensatory.py`.

## D-005: The wobble and bulge search runs, with one pattern for the wobbles and extra bases for a bulge

- ID: `D-005`
- R location: `run_wobble_bulges_pipeline` (`R/run_wobble_bulges_pipeline.R:18-62`), `.get_all_combination_w_b` (`:75-129`), `.get_target_bulged_guides` (`:160-173`), `.sequence_with_tb` (`:185-198`), and `check_bulges_nt` (`:208-225`).
- R behavior: Each call stops, and the rules of the pipeline are wrong in other places too. The golden file `tests/golden/noncanonical/wobble_bulge_r.tsv` and its notes give the evidence. R 4.5.3, Biostrings 2.78.0, and gtools 3.9.5 made the file.
  - `.get_target_bulged_guides` sets `bulge_allowance` in its own frame (`:167`). `.sequence_with_tb` reads it at `:192` and `:197`, and its environment is the namespace. Thus each call stops with `object 'bulge_allowance' not found` (note `pipeline_error`). The call with `bulges_lim = 1` stops with the same message (note `pipeline_error_bulges_lim_1`).
  - The call at `:27` gives `bulges_lim` to the parameter `guide_seq_length` and the guide length to `bulges_lim` (`:160`, note `formals_get_target_bulged_guides`).
  - `.get_all_combination_w_b` sets exactly 4 U positions to `Y` (`:82-93`). A guide with fewer than 4 U stops (`:101`, notes `wobble_seeds_2_u` and `wobble_seeds_0_u`). No seed has `R` at 5 positions (note `wobble_seeds_5_u_in_seed`).
  - `.sequence_with_tb` replaces seed letters with `N` (`:197`). Thus a "bulge" of R matches a mismatch, and no R pattern matches an extra base (the cases `mismatch_t4` and `bulge_after_g5`).
  - `check_bulges_nt` looks for an `N` in the first `nt` letters of the target (`:213`). They pair with g8, g7, and g6, not with g2 to g4 (note `check_bulges_nt`). The R documentation of `nt` names the start of the siRNA and the wobbles (`:11`), and the bulges (`:203`). The R code applies `nt` only to the bulges (`:28`, `:213`).
  - The counts add each bulged string of each wobble seed (`:35-55`). The columns `r_rule_bulges_lim_1` and `r_rule_bulges_lim_4` of the golden file give this count rule, as if the call did not stop. One site gives 12 for `bulges_lim = 1` and 24 for 4 (case `perfect`). A mismatch gives 4 and 12 (case `mismatch_t4`). Each case with an extra target base gives 0.
- New behavior: `wobble_bulge_matches` and `annotate_wobble_bulge_counts` run. The rules follow.
  - The seed is `Seed.MER7M8`, as in R (`:113`). `Guide.target(Seed.MER7M8, wobble=True)` has `R` at each target position that pairs with a guide `U`. Thus one pattern covers each set of G:U pairs.
  - A bulge is 1 to `max_bulge` extra target bases between two adjacent seed positions. The pattern has the letter `N` for each extra base. A site has one bulge at most. `max_bulge` is at most 19993, because a search accepts a pattern of 20000 letters at most.
  - No bulge is next to the `bulge_free` seed positions from g2. With the default 3, a bulge can be only between g5 and g6, g6 and g7, or g7 and g8. As in the R code, `bulge_free` applies only to the bulges, thus a G:U pair can be at each seed position.
  - The search uses `fixed = FALSE` of R for the pattern and the subject (`:39`, `:42`). Thus a subject `N` matches each pattern letter.
  - A sequence has one row for each distinct start of a match. The row takes the first pattern in pattern order that matches at this start. Two patterns with different starts and one end thus count two times, for example `ENSRNOG00000016275` at 69 and 73.
  - `annotate_wobble_bulge_counts` gives each sequence its number of rows, and each gene the maximum over its sequences (`D-003`). It gives only the column `overall_matches` (`:51-55`), not one column for each wobble seed.
- Reason: R gives no result. The R count rule counts one site more than one time, counts mismatches as bulges, and protects the wrong end of the seed.
- Test: `test_wobble_bulge_matches_of_the_cases_D005`, `test_each_r_call_stops_D005`, `test_one_site_counts_once_D005`, `test_a_mismatch_is_not_a_bulge_D005`, `test_an_extra_target_base_is_a_bulge_D005`, `test_no_bulge_is_next_to_g2_to_g4_D005`, `test_each_set_of_wobbles_is_in_one_pattern_D005`, and `test_each_distinct_start_counts_once_D005` in `tests/noncanonical/test_wobble.py`.

## D-006: The benchmark applies the padj cutoff, and the ROC curve uses the labels

- ID: `D-006`
- R location: `benchmark_preds` (`R/benchmark.R:22-38`).
- R behavior: The notes of `tests/golden/analyses/benchmark.tsv` give the evidence for the defects below.
  - In `dplyr::filter(padj <= padj, ...)` (`R/benchmark.R:23`), both names are the column `padj` of the data mask. Thus each row with a `padj` passes, and the argument `padj` has no effect. On the test data, `benchmark_preds(res, padj = 1e-10)` keeps 63 genes. The cutoff keeps 24 genes.
  - `roc.curve(scores.class0 = preds, scores.class1 = res$GS.label)` (`R/benchmark.R:36-38`) takes the labels as the scores of the negative class. Thus R compares the seed counts with the numbers 0 and 1. These 63 genes give the case `benchmark_1`, which has the same PR area. For them, the ROC area is 0.37626790674448363 in R, but 0.59195402298850586 with the labels as weights, as `pr.curve` gets them.
- New behavior: `benchmark_predictions` applies `padj_cutoff`, and `roc_curve` takes the labels as `pr_curve` does.
- Reason: A cutoff without an effect gives a wrong gold standard. The ROC area of R does not measure the predictions.
- Test: `test_the_padj_cutoff_applies_D006` and `test_the_roc_curve_uses_the_labels_D006` in `tests/stats/test_benchmark.py`. `test_benchmark_matches_r` in the same file compares the areas and the curves with the golden files of the patched R function (`tools/golden/patches/benchmark_preds.R`).

## D-007: trrust_scores takes the messages from the regulators

- ID: `D-007`
- R location: `message_passing` (`R/get_network_matches.R:51-122`), called by `get_network_matches` (`:35-39`).
- R behavior: The case `mean` of `tests/golden/network/trrust_matches.tsv` and R probes give the evidence for the defects below.
  - `igraph::neighbors(graph, v)` has the default mode `out`, thus the messages of a gene come from its targets (`:59`).
  - `mean(aggregated_message_rep, as.numeric(features_rep[v]))` gives the own flag of the gene as the argument `trim`, thus R ignores the flag (`:84-85`). The probes `mean(3, 1)` = 3 and `mean(0, 1)` = 0 show this. `length(neighbors != 0)` is the number of neighbors (`:61`).
  - The mode of a message comes from each TRRUST row whose target is the neighbor. `Repression` wins over `Activation` (`:71-79`). `features_act` reads `features_rep[v]` after its change (`:85`, `:110`).
  - `degree_importance` divides by `igraph::degree`, the total degree with each row (`:99`, `:103`, `:109-110`).
  - R gives `GA` the repression value 1, because its target `TFA` is a seed.
- New behavior: `trrust_scores` takes the messages of a gene from its regulators.
  - The score of a gene in a channel is the sum of two terms, divided by 2. The terms are the mean of the weighted seed flags of its regulators in the channel, and its own seed flag.
  - The mode of an edge selects its channel. An `Unknown` edge is in no channel. A gene without a regulator in a channel has the score NaN.
  - With `Aggregator.DEGREE_IMPORTANCE`, the weight of a regulator is 1 divided by the number of its distinct targets, over each mode.
- Reason: A TRRUST edge goes from a regulator to its target. A seed match in a regulator can change the expression of its targets, not of its regulators.
- Test: `test_trrust_scores_use_the_regulators_D007`, `test_trrust_scores_give_the_degree_importance_D007`, and `test_trrust_scores_count_an_unknown_edge_in_the_out_degree_D007` in `tests/network/test_trrust.py`.

## D-008: e_values has no default null, and genes match by gene_id

- ID: `D-008`
- R location: `generate_e_value` (`R/generate_e_value.R:12-24`).
- R behavior: The default `random_siRNA_results_overall` does not exist. Thus a call without a null stops with `object 'random_siRNA_results_overall' not found`. SeedMatchR has no function that makes the null.
  - Genes match by `symbol` (`:14`, `:18-20`).
  - The function reads the column `overall_matches` (`:19`). This is the column of `run_wobble_bulges_pipeline` (`R/run_wobble_bulges_pipeline.R:51-55`).
  - The e-value is the fraction of the null columns with `log(null) >= log(observed)`, and `log(0)` is `-Inf`.
  - Evidence: the notes of `tests/golden/analyses/evalue.tsv`.
- New behavior: The null is a required argument. `e_value_null` makes it from `n_random` shuffles of the guide (default 1000), with the random generator of the caller. Genes match by `gene_id`.
  - `e_value_null` counts the matches of `count_matches`, thus it is the null of `annotate_seed_counts`.
  - `wobble_bulge_null` makes the null of the column `overall_matches` of `annotate_wobble_bulge_counts`. For the same `rng`, shuffle j is the same guide as in `e_value_null`.
- Reason: Without a null, R stops for each call. A symbol can name two genes, and a DE table can have no symbol.
- Test: `test_e_values_match_genes_by_gene_id_D008`, `test_e_value_null_counts_each_shuffle_like_annotate_seed_counts_D008`, `test_e_value_null_uses_the_default_options_of_annotate_seed_counts_D008`, `test_each_shuffle_keeps_the_letters_of_the_guide_D008`, and `test_the_shuffles_are_the_permutations_of_the_generator_in_order_D008` in `tests/test_evalue.py`. `test_e_values_match_r` in the same file compares the e-values with the R values of the golden file. `test_wobble_bulge_null_counts_each_shuffle_like_annotate_wobble_bulge_counts_D008` and `test_wobble_bulge_null_draws_the_shuffles_of_e_value_null_D008` in `tests/noncanonical/test_wobble.py`.

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
- New behavior: `kmer_gene_counts` returns the gene count of each k-mer. `plot_kmer_histogram` draws the histogram (`D-066`).
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
  - Evidence: `tests/golden/analyses/distances.tsv` with its notes, and an R probe with R 4.5.3 and SeedMatchR 2.0.0.
  - The probe calls `match_distances` with a `GRangesList` whose first element and last element have the name `a`. The first element has the starts 10 and 30. The last element has the starts 1, 5, and 9. Each range has the width 7.
  - The result is only `$a` = 4 8 4, the distances of the last element.
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

## D-013: string_gene_network maps the proteins with the annotation of the Reference

- ID: `D-013`
- R location: `get_stringPPI_to_stringGGI` (`R/get_stringPPI_to_stringGGI.R:32-60`).
- R behavior: `biomaRt::getBM` gives the gene of each protein. The rat query uses the current Ensembl release of BioMart. The human query uses the archive `grch37.ensembl.org` (`:32-39`).
  - For a value of `sp` other than `rat` and `human`, `mart` has no value. Evidence: the source. There is no R run, because biomaRt is not in `tools/golden/.rlib`, and the tests do not use the internet.
- New behavior: `string_gene_network` maps the protein IDs with the tables `protein` and `tx` of the annotation of the `Reference`. For a GTF reference, `load_gtf` makes the `protein` table from the attribute `protein_id` of the `CDS` lines.
- Reason: The map then uses the Ensembl release of the sequences, and it uses no internet. A GRCh37 map does not agree with a GRCh38 reference.
- Test: `test_string_gene_network_maps_with_the_protein_table_D013` and `test_string_gene_network_maps_with_the_protein_ids_of_a_gtf_D013` in `tests/network/test_string_network.py`. The `test_protein_genes_*` tests in `tests/annotation/test_ensdb.py` check the method of the EnsDb.

## D-014: The documentation uses only functions that exist

- ID: `D-014`
- R location: the vignettes `vignettes/articles/SeedMatchR_miRNA_analysis.Rmd` (`:64`, `:76`, `:178`) and `vignettes/articles/SeedMatchR_build_PPI_network.Rmd` (`:30`).
- R behavior: The miRNA vignette calls `filter_deseq` (`:64`), `get_feature_seqs` (`:76`), and `de_fc_ecdf` (`:178`). The network vignette calls `get_network_matches` without `organism` (`:30`). An R probe with R 4.5.3 and SeedMatchR 2.0.0 gives the evidence:
  - `exists(name, envir = asNamespace("SeedMatchR"), inherits = FALSE)` is `FALSE` for `filter_deseq`, `de_fc_ecdf`, and `get_feature_seqs`. It is `TRUE` for `filter_res` and `deseq_fc_ecdf`.
  - `get_network_matches(c("A"), aggregator = "mean")` stops with `argument "organism" is missing, with no default`.
  - The export list of SeedMatchR 2.0.0 (`tests/data/SeedMatchR.NAMESPACE`) has `deseq_fc_ecdf` (`:25`), `filter_res` (`:28`), and `get_network_matches` (`:35`). It has no line for `filter_deseq`, `de_fc_ecdf`, or `get_feature_seqs`.
- New behavior: The examples of `README.md` use only the public names of seedmatch, and a test runs each Python example. Each `seedmatch.` name and each `seedmatch-fetch` command of `src/seedmatch/r-map.md` exists.
- Reason: An agent copies the examples of the documentation. An example that calls an absent function stops with an error, and the agent cannot use it.
- Test: `test_each_python_name_of_the_r_map_exists_D014` in `tests/test_r_map.py`. The tests of `tests/test_readme.py` run each Python example of `README.md`.

## D-015: A seed-count column replaces a column with the same name

- ID: `D-015`
- R location: `.get_counts_column` (`R/SeedMatchR.R:214-217`). `SeedMatchR` calls it in `DESEQ2` mode.
- R behavior: `cbind(res, seed.name = match.counts)` appends the counts as the last column, and the next statement gives this column the seed name. If `res` has a column with this name, the result has two columns with the same name. `res$mer7m8` then reads the old column. The test data (`tests/data/DESEQ2_results.tsv`) has a `mer7m8` column. R probe with Biostrings 2.78.0: `names(SeedMatchR(seqs = seqs, sequence = guide, res = res, seed.name = "mer7m8"))` ends with `mer7m8, mer7m8` (`tests/golden/seedmatchr/deseq2_mode.tsv` records this in its header). After `res$mer7m8 <- 99L`, `$mer7m8` of the result is 99 for each gene, and the last column holds the counts.
- New behavior: If the table has a column with this name, `annotate_seed_counts` and `annotate_wobble_bulge_counts` replace the values of this column. The column keeps its position. A new column goes at the end of the table.
- Reason: A table with two columns of the same name is ambiguous. In pandas, `de["mer7m8"]` then returns a table, not a column. A caller who reads the count column must get the new count. A `de` that has this column twice raises `InputError` (`D-051`).
- Test: `test_existing_column_is_replaced_not_duplicated_D015` and `test_replaced_column_keeps_its_position_and_gets_the_new_values_D015` in `tests/report/test_annotate.py`. `test_annotate_wobble_bulge_counts_replace_a_column_with_the_name_in_place_D015` in `tests/noncanonical/test_wobble.py`.

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
- New behavior: `encode_dna`, `reverse_complement`, `SequenceSet.from_records`, and `read_fasta` raise `InputError` for each of these letters. The message gives the letter and its 1-based position. For a record, the message also gives the ID of the record. For `read_fasta`, the message also gives the path and the line, and the position is in the line.
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

## D-047: e_values rejects a count that is not valid, and gives each row its own e-value

- ID: `D-047`
- R location: `generate_e_value` (`R/generate_e_value.R:16-21`).
- R behavior: A count of `NA` gives the e-value 0. A negative count gives 0 and the warning `NaNs produced`.
  - The R probes use R 4.5.3 and SeedMatchR 2.0.0. In the first probe, the null row is `a = 1, 2, 3, 0`. `generate_e_value(data.frame(symbol = "a", overall_matches = Inf), null)` gives 0, because no finite value of the null is at least `Inf`. The count `-Inf` gives 0 and the warning `NaNs produced`.
  - R probe, with the observed count `a = 1`: the null row `a = 1, 2, 3, 0` gives 0.75. If the first cell is `NA`, a negative number, or `-Inf`, R gives 0.5. R does not count the cell, but it keeps the cell in the denominator. A negative number and `-Inf` also give the warning `NaNs produced`. If the first cell is `Inf`, R gives 0.75, because R counts the cell.
  - A gene with two rows in the null gives an e-value above 1. For the null rows `a = 1, 2, 3, 0` and `a = 5, 5, 5, 5`, the count `a = 2` gives 1.5.
  - Two rows of one gene in `result_output` both get the value of a recycled comparison. For the null row `a = 1, 2, 3, 0`, the counts `a = 1` and `a = 3` both get 0.5, but one by one they get 0.75 and 0.25.
  - R drops a row whose gene is not in the null before it reads the count (`R/generate_e_value.R:14`). Thus such a row gives no e-value and no message.
  - R probe, with the null row `a = 1, 2, 3, 0`: `generate_e_value(data.frame(symbol = c("a", "zz"), overall_matches = c(1, value)), null)` gives only the row `a`, with the e-value 0.75. This holds for `value` equal to `NA`, -1, `Inf`, and `-Inf`, and R gives no warning.
  - Evidence: the notes of `tests/golden/analyses/evalue.tsv`, and the R probes above.
- New behavior: `e_values` raises `InputError` for a count that is NaN, infinite, or negative. This rule covers `de` and the null. It covers also a row of `de` whose gene is not in the null. It raises `InputError` also for a gene with two rows in the null. The message names the column or the null, the value (`NaN`, `inf`, `-inf`, or the negative number), and the gene. A count that is not an integer, for example 1.5, is valid. Each row gets its own e-value.
- Reason: The e-value 0 is the strongest result, thus a count that is not known must not give it. A count that is not finite is not a count. A fraction above 1 has no meaning. A count that is not valid shows a broken input, also in a row that R drops.
- Test: `test_e_values_reject_a_nan_count_D047`, `test_e_values_reject_a_negative_count_D047`, `test_e_values_reject_an_infinite_count_D047`, `test_e_values_reject_a_repeated_gene_of_the_null_D047`, and `test_each_row_of_a_gene_gets_its_own_e_value_D047` in `tests/test_evalue.py`. `test_e_values_reject_an_invalid_count_for_a_gene_that_is_not_in_the_null_D047` shows the rule for a row outside the null. The tests `test_e_values_reject_a_nan_count_in_the_null_D047`, `test_e_values_reject_a_negative_count_in_the_null_D047`, and `test_e_values_reject_an_infinite_count_in_the_null_D047` in the same file cover the counts of the null.

## D-048: The benchmark rejects a NaN score, a cutoff that is not a number, and a gold standard with one class

- ID: `D-048`
- R location: `benchmark_preds` (`R/benchmark.R:22-38`) and `pr.curve` and `roc.curve` of PRROC 1.4.
- R behavior: R gives no message, or a message that does not name the cause. The notes of `tests/golden/analyses/benchmark.tsv` give some of the evidence. R probes give the rest. A probe uses R 4.5.3, SeedMatchR 2.0.0, dplyr 1.2.1, and PRROC 1.4 on the test data. No golden file keeps the probe results.
  - To do a probe again, run `source("tools/golden/R/common.R")` and `res <- load_test_data()$res` from the repository root.
  - The area probe is `with_seedmatchr_patch("tools/golden/patches/benchmark_preds.R", SeedMatchR::benchmark_preds(res, padj = NaN))`. Read the integral area, the Davis-Goadrich area, and the ROC area as `$PRC$auc.integral`, `$PRC$auc.davis.goadrich`, and `$ROC$auc`.
  - The gene probe is `sum(res$gene_id %in% dplyr::pull(dplyr::filter(res, padj <= !!padj, log2FoldChange < !!l2fc), gene_id))`. Set the variables `padj` and `l2fc` to the cutoffs first.
  - A NaN score in one row stops the function with `NAs are not allowed in subscripted assignments` (the notes).
  - A NaN cutoff keeps no gene, because each comparison with NaN fails. The area probe with `padj = NaN` and with `l2fc = NaN` gives NaN for the integral area, the Davis-Goadrich area, and the ROC area.
  - A gold standard without a positive (`padj = 0`) gives NaN for the integral area, the Davis-Goadrich area, and the ROC area (the notes).
  - A gold standard without a negative (`padj = 1` and `l2fc = Inf`) gives 1 for the integral area, 1 for the Davis-Goadrich area, and NaN for the ROC area (the notes).
  - A text cutoff compares the values as text. The gene probe with `padj = "0.15"` keeps no gene. With `padj = 1e-10`, `l2fc = "-1"` keeps 5 genes, but `l2fc = -1` keeps 19 genes.
  - A bool cutoff is the number 1. The gene probe with `padj = 0.15` and `l2fc = TRUE` keeps 68 genes, and `l2fc = 0` keeps 63 genes.
- New behavior: `benchmark_predictions` raises `InputError` for a NaN score, and for a cutoff that is a bool, text, or NaN. It also raises `InputError` for a gold standard without a positive row and for a gold standard without a negative row. A score of infinity is valid. A row with a NaN `padj` or a NaN `log2FoldChange` does not pass the cutoffs. It is a negative unless another row of its gene passes, as in R (the note on `benchmark_na`).
- Reason: An area from a broken input looks like a result, and no message tells the caller.
- Test: `test_benchmark_rejects_a_nan_score_D048`, `test_benchmark_rejects_a_missing_value_of_a_nullable_integer_score_D048`, `test_benchmark_rejects_a_nan_cutoff_D048`, `test_benchmark_rejects_a_cutoff_that_is_not_a_number_D048`, `test_benchmark_rejects_a_gold_standard_without_a_positive_D048`, and `test_benchmark_rejects_a_gold_standard_without_a_negative_D048` in `tests/stats/test_benchmark.py`.

## D-049: read_mirdb rejects a score that is not a number, and keeps the ID text as it is

- ID: `D-049`
- R location: `get_example_data` with `example.type = "mirna"` (`R/get_data.R:43-46`, `read.table` with the defaults).
- R behavior: `read.table` converts each column and gives no message. The evidence is in the notes of `tests/golden/analyses/mirdb.tsv`.
  - A score `abc` makes the column `target.score` text, with the values `50` and `abc`. Then the comparison `mirdb.rat[,3] >= 90` of the vignette (`vignettes/articles/SeedMatchR_miRNA_analysis.Rmd:149`) compares text.
  - An empty score gives `NA`.
  - The ID text `NA` becomes `NA`, the R value for no data. The species is then `NA`, and the miRNA name is empty.
  - A `#` in a field starts a comment, and R stops with `line 2 did not have 3 elements`.
  - A quote character in a field loses rows, and R gives no message. For a file of 4 lines with a quote character in line 2, R gives 2 rows, with the IDs `hsa-c` and `hsa-d`. The default quote of `read.table` includes the quote character.
  - A line with 2 fields stops R with `line 2 did not have 3 elements`. The port also raises `InputError` for it, thus this case is no change.
  - R probes use R 4.5.3 and the call of `R/get_data.R:43-46`. The lines are in a file `path`, and the call is `m <- utils::read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)`. No golden file keeps the probe results.
  - R reads a line of any length. After `writeLines(c("hsa-a\tNM_1\t50", paste0("hsa-", strrep("x", 5000), "\tNM_2\t60")), path)`, `nrow(m)` is 2.
  - R reads a score with spaces before it, after it, or on each side as a number, and the port does the same. After `writeLines(c("hsa-a\tNM_1\t50", "hsa-b\tNM_2\t 60 "), path)`, `class(m$V3)` is `numeric`. A score with a space inside (`6 0`) makes the column text.
  - R reads the scores `0x10`, `Inf`, `NaN`, and `1e999` as numbers. After `writeLines(c("hsa-a\tNM_1\t50", paste0("hsa-b\tNM_2\t", score)), path)`, `m$V3[2]` is 16, `Inf`, `NaN`, and `Inf`.
  - A NUL character in a line gives the warning `line 2 appears to contain embedded nulls`. R ends a score at the NUL: the score `6<NUL>0` gives the number 6. In an ID or a RefSeq ID, R stops with `line 2 did not have 3 elements`. The file is `writeBin(c(charToRaw("hsa-a\tNM_1\t50\n"), charToRaw("rno-b\tNM_2\t6"), as.raw(0), charToRaw("0\n")), path)`.
- New behavior: `read_mirdb` raises `InputError` for a score that is not a finite decimal number. This includes an empty score, `NA`, `NaN`, and `Inf`. The ID fields keep `#`, quote characters, and `NA` as text, and each line is one row. A line of more than 1024 characters with its line end raises `InputError`, and R reads such a line. A line with a NUL character raises `InputError` that names the line.
- Reason: A text score gives wrong comparisons and no message. A row that R loses is a loss of data without a message. The miRDB format has no quotes and no comments. R cuts a score at a NUL without an error, and a cut score is wrong data. The limit on the line bounds the memory of a chunk of 2^18 lines in bytes. The longest line of the miRDB 6.0 file of `seedmatch-fetch mirdb` has 51 characters with its line end.
- Test: `test_read_mirdb_rejects_a_score_that_is_not_a_number_D049`, `test_read_mirdb_keeps_a_hash_a_quote_and_the_text_na_D049`, `test_read_mirdb_rejects_a_line_above_the_limit_D049`, and `test_read_mirdb_rejects_a_nul_character_D049` in `tests/test_mirdb.py`. In the same file, `test_read_mirdb_keeps_each_line_that_a_quote_character_would_join_D049` shows that the port keeps the 4 lines of the last note.

## D-050: auc_integral loses fewer digits than the difference of two logarithms

- ID: `D-050`
- R location: `compute.pr` of PRROC 1.4, which `pr.curve` calls. The expression for `v` is at lines 71 and 72 of the printed function `PRROC:::compute.pr`. Line 1 is the line that starts with `function`.
- R behavior: `v <- (tp/r.fg - tp.prev/r.fg - b/a * (log(a * tp/r.fg + b) - log(a * tp.prev/r.fg + b)))/a` gives the area of a step as the difference of two logarithms. If the negatives far outnumber the positives, `b / a` is large and the two logarithms are almost equal. Then the difference loses digits. R probe with R 4.5.3 and PRROC 1.4. The exact value is an evaluation of the same formula in decimal arithmetic with 50 digits.
  - Samples 1 and 2 have 200,000 scores. The generator `numpy.random.default_rng(20261125)` chooses the k positions of the positives with `choice(200_000, k, replace=False)`. The score is `normal(size=200_000) + 0.8 * label`.
  - Sample 1 has k = 11. `pr.curve(scores.class0 = score, weights.class0 = label)$auc.integral` is 0.00021605166950061649, and the exact value is 0.00021605166048556729. The relative error of R is 4.2e-8.
  - Sample 2 has k = 2044. R gives 0.032534101333943015, and the exact value is 0.032534101333588007. The relative error of R is 1.1e-11.
  - Sample 3 has 1,000,000 scores `numpy.random.default_rng(4711).normal(size=1_000_000)`, and its one positive is the lowest score. R gives 4.9981627192607903e-07, and the exact value is 5.0000016666674998e-07. The relative error of R is 3.7e-4.
- New behavior: `pr_curve` computes `ln(a * r_k + b) - ln(a * r_{k-1} + b)` as `log1p(a * (r_k - r_{k-1}) / c)`. Here c is the value of `a * r + b` at the step before, that is (TP + FP) / P of that step. The formula is the same in exact arithmetic. The port loses fewer digits, but it still subtracts a term from the width of the step, thus a smaller error remains. On the same samples, the relative error is 2.4e-13 for sample 1, 0 for sample 2, and 3.2e-10 for sample 3.
- Reason: The area loses digits that the double format can partly keep. The R value thus differs from the port by up to the error of R. The goldens of small samples agree within 1e-16.
- Test: `test_auc_integral_keeps_its_precision_for_few_positives_D050` compares samples 1 and 2 with the exact value, within a relative 1e-12. `test_auc_integral_has_a_small_error_for_one_positive_in_a_million_scores_D050` compares sample 3 within a relative 1e-9. Each test is in `tests/stats/test_curves.py`.

## D-051: A repeated column that a function uses raises InputError

- ID: `D-051`
- R location: `ecdf_stat_test` (`R/deseq_fc_ecdf.R:248`, `:252`) and `.assign_groups` (`R/deseq_fc_ecdf.R:176`), which `deseq_fc_ecdf` calls (`:100`, `:124`). Also `.get_counts_column` (`R/SeedMatchR.R:193`, `:207`, `:214-217`). The calls of `dplyr::filter` on the table are in `filter_res` (`R/filter_res.R:55-70`), `SeedMatchReport` (`R/SeedMatchReport.R:70`), `generate_e_value` (`R/generate_e_value.R:14`), and `benchmark_preds` (`R/benchmark.R:23`).
- R behavior: `res$log2FoldChange` and `res$gene_id` take the first column of that name. R gives no message. `dplyr::filter` stops for each repeated name. The R probes use R 4.5.3, SeedMatchR 2.0.0, and dplyr 1.2.1. No golden file keeps the probe results.
  - The table is `res <- data.frame(gene_id = paste0("g", 0:19), log2FoldChange = rnorm(20))` after `set.seed(1)`. The sets are `t <- res$gene_id[1:10]` and `b <- res$gene_id[11:20]`.
  - `dup <- cbind(res, data.frame(log2FoldChange = rep(5, 20)))` adds a second `log2FoldChange` column. `ecdf_stat_test(dup, t, b)$p.value` is 0.4195804, the value for `res`. `deseq_fc_ecdf(dup, list(a = t, b = b), stats_test = "KS")` runs.
  - `cbind(data.frame(log2FoldChange = rep(5, 20)), res)` puts the repeated column first, and the p-value is 1. Thus the result of R depends on the column order.
  - The table `res4` is `data.frame(gene_id = c("g1", "g2", "g3", "g4"), log2FoldChange = c(0.5, -1, 0.2, 0.1))`. `dup4 <- cbind(res4, data.frame(gene_id = c("x1", "x2", "x3", "x4")))` adds a second `gene_id` column. `.assign_groups(dup4, list(a = c("g1", "g2"), b = c("g3", "g4")))$Group` is `a a b b`. R uses the first `gene_id`.
  - `SeedMatchR(seqs = seqs, sequence = guide, res = dup4, res.format = "DESEQ2", shared_genes = FALSE)` gives the count column and keeps each `gene_id` column. The guide is `UUAUAGAGCAAGAACACUGUUUU`, and the seed is `mer7m8` with the target `CTCTATA`.
  - The sequences `g1`, `g2`, and `g3` of `seqs` are `GGGGGGGGGGGGGGGGGGGGAAACAGTGTTCTTGCTCTATAAGGGGG`, `CCCCCCCCCCCCCCCCCCCCCCCCC`, and `TTTTTTTTTTTTTTTTTTTTTT`. Only `g1` has a match. The count column is 1 0 0 0, thus R counts with the first `gene_id`. If `x1` to `x4` is the first `gene_id` column, the count column is 0 0 0 0.
  - With the default `shared_genes = TRUE`, the same call stops in `dplyr::filter` (`R/SeedMatchR.R:193`) with `Can't transform a data frame with duplicate names.`
  - For `res <- data.frame(gene_id = c("g1", "g2", "g3", "g4"), mer7m8 = c(5, 6, 7, 8), mer7m8 = c(9, 10, 11, 12), check.names = FALSE)` and `shared_genes = FALSE`, the names are `gene_id`, `mer7m8`, `mer7m8`, and `mer7m8`. The old columns keep their values, and the new column is 1 0 0 0. With `shared_genes = TRUE`, R stops with the message above.
  - `filter_res`, `SeedMatchReport`, `generate_e_value`, and `benchmark_preds` stop with the same message for a repeated column that they do not read. Each probe table has a column `extra` that occurs twice. `SeedMatchR` with `shared_genes = TRUE` stops for such a table too. With `shared_genes = FALSE`, it accepts the table.
- New behavior: A function raises `InputError` for a repeated column that it reads or writes. For a column `log2FoldChange` that occurs twice, the message is `DE table has 2 columns 'log2FoldChange'; expected one column for each name`.
  - `ecdf_test` and `assign_groups` raise it for a repeated `gene_id` or `log2FoldChange` column. R uses the first column.
  - `annotate_seed_counts` and `annotate_wobble_bulge_counts` raise it for a repeated `gene_id` column and for a repeated column with the name of the count column. R stops for `shared_genes = TRUE`. For `FALSE`, R uses the first `gene_id`, and it appends the counts after the old columns.
  - `filter_de`, `e_values`, `benchmark_predictions`, and `seed_match_report` raise it for a repeated column that they read. R stops in `filter_res`, `generate_e_value`, `benchmark_preds`, and `SeedMatchReport`.
  - A repeated column that the function does not read is valid. A function that returns the rows of `de` keeps the column. R stops for such a table in `filter_res`, `generate_e_value`, `benchmark_preds`, `SeedMatchReport`, and `SeedMatchR` with `shared_genes = TRUE`.
- Reason: A repeated column is ambiguous. The result of R depends on the column order, or R stops. A repeated count column gives a wrong count in pandas, and no message. The message names the column, thus the caller can correct the table.
- Test: `test_ecdf_test_rejects_a_repeated_fold_change_column_D051` and `test_assign_groups_rejects_a_repeated_gene_id_column_D051` in `tests/stats/test_ecdf.py`. `test_a_repeated_gene_id_column_is_an_error_D051`, `test_a_repeated_count_column_is_an_error_D051`, and `test_a_repeated_column_that_is_not_read_is_kept_D051` in `tests/report/test_annotate.py`. `test_annotate_wobble_bulge_counts_reject_a_repeated_gene_id_or_count_column_D051` in `tests/noncanonical/test_wobble.py`. `test_filter_de_rejects_a_repeated_padj_column_D051` in `tests/test_de.py`. The tests of `check_columns` and `check_unique_columns` with `_D051` in their names are in `tests/test_checks.py`. One of them shows that a repeated column that is not required is valid.

## D-052: Each sequence is searched, also for a repeated ID

- ID: `D-052`
- R location: `get_compensatory_sequences` (`R/get_compensatory_sequences.R:36-41`, `:90`, and `:106-108`).
- R behavior: R selects the sequences with a seed match by name (`sequences[names]`). For a name that occurs two times, this gives the first sequence two times (note `repeated_name_r`).
  - The merge by `gene_name` (`:90`) then pairs each seed row with each region row of the name.
  - The case `repeated_name_r` has two sequences with the name `g`, and each sequence has a site. R gives 4 rows, and each row belongs to the first sequence. The second sequence has no row.
  - Evidence: the cases and the notes of `tests/golden/noncanonical/compensatory.tsv`, which come from the R function with the patch of `D-004`.
- New behavior: Each sequence of the set is searched. A pair holds a seed match and a region of the same sequence, by the index of the sequence in the set.
- Reason: A `SequenceSet` can hold two sequences with one ID, for example with `PerGene.EACH_TRANSCRIPT` and gene IDs as names.
- Test: `test_each_sequence_is_searched_for_a_repeated_id_D052` in `tests/noncanonical/test_compensatory.py`.

## D-053: A region is dropped only for a longer region of the same seed match

- ID: `D-053`
- R location: `get_compensatory_sequences` (`R/get_compensatory_sequences.R:123-125`).
- R behavior: Two `dplyr` filters keep the longest region of each transcript and region end, and then of each transcript and region start. The groups ignore the seed match.
  - Thus a longer region of one seed match drops the region of another seed match.
  - The case `two_seed_matches_r` has two seed matches in one sequence. R gives one row, and the seed match at 13 has no row.
  - Evidence: the cases and the notes of `tests/golden/noncanonical/compensatory.tsv`, which come from the R function with the patch of `D-004`.
- New behavior: A pair is dropped if a pair of the same seed match has a region that contains its region and is longer. Pairs with the same region and different guide windows stay.
- Reason: Each seed match is a site of its own. On the R goldens without this case, the R rule and the port rule give the same rows.
- Test: `test_a_region_of_one_seed_match_does_not_drop_a_region_of_another_D053` and `test_a_region_inside_a_longer_region_of_its_seed_match_is_dropped` in `tests/noncanonical/test_compensatory.py`.

## D-054: Each row gives the guide positions of its region and the mismatches of its seed match

- ID: `D-054`
- R location: `relative_position_seq` (`R/get_compensatory_sequences.R:3-8`), `:127-128`, and `:106-108`.
- R behavior: `relative_position_seq` searches the full guide and joins each start with a comma.
  - A region that occurs two times in the guide gives two equal rows with the text `8,12` (case `repeated_region_r`).
  - A region at g8 that also occurs at g2 gives the text `2,8` (case `region_at_g2_r`, note `relative_position_seq_g6`).
  - `mismatch_in_seed` is 1 if the text of the seed match differs from the target (`:106-108`). Thus a match with 2 mismatches gives 1 (case `seed_mm2_r`).
  - Evidence: the cases and the notes of `tests/golden/noncanonical/compensatory.tsv`, which come from the R function with the patch of `D-004`.
- New behavior: `guide_start` and `guide_end` are the integers of the window of the row. `seed_mismatches` is the number of mismatches, and a position outside the sequence is a mismatch. For `max_seed_mismatch` of 0 or 1, it equals the R value.
- Reason: A text of positions is not a number, and the position of g2 is not in the region. The flag hides the number of mismatches.
- Test: `test_each_row_gives_the_guide_positions_of_its_region_D054` and `test_seed_mismatches_is_the_number_of_mismatches_D054` in `tests/noncanonical/test_compensatory.py`.

## D-055: CompensatoryOptions rejects lengths and distances that select no region

- ID: `D-055`
- R location: `get_compensatory_sequences` (`R/get_compensatory_sequences.R:46-54` and `:92`).
- R behavior: `compensatory_region_length_min = 6` with `compensatory_region_length_max = 3` gives the lengths 6, 5, 4, and 3, and R returns 3 rows of length 6 (case `min_above_max_r`).
  - A minimum above the 16 positions of g8 to g23 searches the region g8 to g23 three times (note `min_17`).
  - `limit = 2` gives no row. R keeps `seed_start > end + 1` and `seed_start < end + 2`, and no integer obeys these two tests.
  - Evidence: the cases and the notes of `tests/golden/noncanonical/compensatory.tsv`, which come from the R function with the patch of `D-004`. No case has `limit = 2`, and the filter at `:92` gives the `limit = 2` result.
- New behavior: `InputError` for `max_length` below `min_length`, for a length above the number of positions of the region, and for a `max_distance` below 1. A longest region of more letters than the pattern limit of the search also raises `InputError`, before the first search.
- Reason: R gives a result for other lengths than the arguments, and it gives no message.
- Test: `test_compensatory_options_reject_lengths_that_select_no_region_D055`, `test_a_short_guide_has_no_room_for_the_region_D055`, `test_the_message_names_one_position_in_the_singular_D055`, and `test_a_window_above_the_pattern_limit_raises_before_a_search_D055` in `tests/noncanonical/test_compensatory.py`.

## D-056: The compensatory region starts after the seed

- ID: `D-056`
- R location: `get_compensatory_sequences` (`R/get_compensatory_sequences.R:44`).
- R behavior: The region is always g8 to the end of the guide. With the `mer7m8` seed, g8 pairs with the seed match and with the region of one row. In the golden `vignette`, 10 of the 104 rows have a region from g8.
- New behavior: The default region starts after the seed window, at g9 for `Seed.MER7M8` and at g8 for `Seed.MER6`. `first_position=8` gives the R region.
- Reason: One guide base pairs with one target base.
- Test: `test_the_region_starts_after_the_seed_D056`, `test_the_default_region_start_follows_the_seed_window_D056`, and `test_the_default_region_starts_after_a_custom_seed_D056` in `tests/noncanonical/test_compensatory.py`.

## D-057: read_trrust keeps the first line of the TRRUST file

- ID: `D-057`
- R location: `get_network_matches` (`R/get_network_matches.R:25`, `:27`).
- R behavior: `read.table(..., header = TRUE, sep = "\t")` reads line 1 as the column names, but the TRRUST file has no header.
  - A planner probe with R 4.5.3 on the human file gives 9395 rows, not 9396. R loses the edge `AATF BAX Repression 22909821`.
  - A `'` or a `#` stops R with `line 1 did not have 4 elements`. The text `NA` becomes the value `NA`. The notes of `tests/golden/network/trrust_rows.tsv` give the evidence.
- New behavior: Each line is an edge. The fields keep quote characters, `#`, and `NA` as text. A bad mode, an empty gene, another number of fields, a NUL character, and a line of more than 1024 characters raise `InputError`.
- Reason: Line 1 is an edge of the data. A gene symbol is text, and the TRRUST format has no quotes and no comments. The limit bounds the memory of a chunk.
- Test: `test_read_trrust_keeps_the_first_line_D057`, `test_read_trrust_keeps_quote_hash_and_na_text_D057`, `test_read_trrust_keeps_a_double_quote_character_D057`, and `test_read_trrust_rejects_a_bad_line_D057` in `tests/network/test_trrust.py`.

## D-058: trrust_scores warns if no seed gene is in the network

- ID: `D-058`
- R location: `get_network_matches` (`R/get_network_matches.R:36-37`), and the call of the vignette (`vignettes/articles/SeedMatchR_build_PPI_network.Rmd:28-30`).
- R behavior: If no gene of `ENSG_Genelist` is a vertex, each value is 0, and R gives no message. The vignette gives `res$gene_id`, the Ensembl gene IDs, but TRRUST names each gene by its symbol. The case `mean_no_seed` of `tests/golden/network/trrust_matches.tsv` gives the evidence.
- New behavior: `trrust_scores` warns with `SeedMatchWarning`. The message names the first gene of the network as an example.
- Reason: IDs of another type give scores that look valid. The function warns and does not raise an error, because a valid list of genes can have no gene of the network.
- Test: `test_trrust_scores_warn_without_a_seed_gene_in_the_network_D058` in `tests/network/test_trrust.py`.

## D-059: read_string_links rejects a bad line of the STRING file

- ID: `D-059`
- R location: `get_stringPPI_to_stringGGI` (`R/get_stringPPI_to_stringGGI.R:15-16`), which calls `load` of STRINGdb 2.22.0 through `get_graph`.
- R behavior: `load` reads the file with `read.table(temp, sep = " ", header=TRUE, stringsAsFactors=FALSE, fill = TRUE)` (`R/rstring.R:1064`). It keeps the rows with `PPI$combined_score >= score_threshold` (`:1067`).
  - A text score makes the column text, and `>=` then compares text. A planner probe with R 4.5.3 keeps `abc` and `99`, but not `1000`.
  - The result of R depends on the position of the bad line, because `read.table` finds the number of columns from the first five lines. Each case below is a line in the first five lines. The notes of `tests/golden/network/string_links.tsv` also give a late `'` and a late line with 4 fields.
  - A line with 2 fields or a `#` gives a row of `NA` values in the subset `>= 400`. A `'` gives a row of `NA` values if a later line closes the quote. If no line closes it, R reads no row and warns `incomplete final line found by readTableHeader`.
  - A line with 4 fields moves the columns, and R keeps the row.
  - A score of 450.5 or 5000 gives a numeric column, and R keeps the row. A score of -5 is below the threshold, and R drops it.
  - A late line with 4 fields gives two rows, and the second row has `NA` values. A late `'` joins the next lines into one row, and R warns `EOF within quoted string`.
- New behavior: Line 1 must be the header `protein1 protein2 combined_score`. Each other line has 3 fields with one space between them. The first two fields hold a protein ID, and the third field holds a score with an integer value from 0 to 1000.
  - Another line, an empty protein ID, a NUL character, and a line of more than 1024 characters raise `InputError` with the line number. `#` and quote characters are text.
- Reason: STRING gives each combined score as an integer from 0 to 1000. A text score gives a wrong selection and no message. The limit bounds the memory of a chunk.
- Test: `test_read_string_links_reject_a_bad_line_D059`, `test_read_string_links_keep_hash_and_quote_text_D059`, `test_read_string_links_reject_a_long_line_D059`, and `test_read_string_links_reject_a_nul_character_D059` in `tests/network/test_string_links.py`.

## D-060: string_gene_network joins the links of each protein of a gene

- ID: `D-060`
- R location: `get_stringPPI_to_stringGGI` (`R/get_stringPPI_to_stringGGI.R:62-64`).
- R behavior: `adj_matrix[!duplicated(newnames), !duplicated(newnames)]` keeps only the row and the column of the first protein of each gene. R drops the links of the other proteins, and then the genes without a link. R keeps a link of a protein with itself as a loop of the gene.
  - Evidence: `tests/golden/network/string_ppi.tsv`, case `p50`, has the pairs G1-G4 and G2-G6. Gn is the fixture gene `ENSSYNG<n>`, with n as 11 digits. The case does not have G1-G5 of the second protein of G1, or G2-G8 of the second protein of G2.
  - Evidence: on the rat file of STRING 12.5 and the EnsDb 109, 96 genes of the default network have more than one protein. Each gene pair of R is a pair of the port, and the port has 1,112 rows more.
  - Evidence: the steps of `R/get_stringPPI_to_stringGGI.R:22-64` on the links P1-P1 and P1-P6 give the entry `G1 G1 = 1` in the matrix of the gene pairs (R 4.5.3).
- New behavior: The links of a gene are the links of each of its proteins. A gene pair gets the highest score of its links. A link of two proteins of one gene gives no pair. A link of a protein with itself gives no pair.
- Reason: A gene-level result uses each sequence of the gene (`D-003`). In the same way, the network uses each protein of the gene.
- Test: `test_string_gene_network_joins_the_proteins_of_a_gene_D060` and `test_string_gene_network_gives_a_pair_the_highest_score_of_its_links_D060` in `tests/network/test_string_network.py`.

## D-061: string_gene_network drops a protein without a gene

- ID: `D-061`
- R location: `get_stringPPI_to_stringGGI` (`R/get_stringPPI_to_stringGGI.R:43-60`).
- R behavior: The name of a protein is the text of its ID between the first and the second `.` (`:43-44`). A protein that BioMart does not map keeps this name (`:56-58`). Thus the result has gene IDs and protein IDs. If no protein maps, the result has only protein IDs, and R gives no message. R drops the text after a second `.`.
  - Evidence: `tests/golden/network/string_nodes.tsv`, case `p0`, has the nodes `ENSSYNP00000000098` and `ENSSYNP00000000099`.
  - Evidence: on the rat file of STRING 12.5 and the EnsDb 109, 628 of the 8,809 proteins of the default network have no gene.
  - Evidence: the expression of `:43-44` gives `AT1G01010` for the ID `3702.AT1G01010.1` (R 4.5.3).
- New behavior: The network has only gene pairs, and `GeneNetwork.unmapped_proteins` lists each protein without a gene. The protein ID is the text after the first `.`, also if the ID has a second `.`. An ID without `.`, with no text before it, or with no text after it raises `InputError`. Links that join no two genes also raise `InputError`.
- Reason: A protein ID is not a gene ID, thus no row of `de` can match it.
- Test: `test_string_gene_network_drops_the_unmapped_proteins_D061`, `test_string_gene_network_gives_no_pair_for_a_link_of_an_unmapped_protein_D061`, `test_string_gene_network_keeps_the_text_after_the_first_dot_of_a_protein_ID_D061`, `test_string_gene_network_rejects_links_without_a_gene_pair_D061`, and `test_string_gene_network_rejects_a_protein_id_without_a_taxon_prefix_D061` in `tests/network/test_string_network.py`.

## D-062: string_scores takes the highest count of a gene, and rejects a bad count

- ID: `D-062`
- R location: `create_graph` (`R/create_graph.R:13-23`).
- R behavior: `f[gene] = res$overall_matches[which(res$gene_id == gene)]` assigns the counts of all rows of a gene to one value. `tests/golden/network/string_features.tsv` gives the evidence.
  - For a gene with more than one row in `res`, R keeps the first count and warns `number of items to replace is not a multiple of replacement length`. The case `repeated` gives the gene with the counts 1 and 3 the feature 1.
  - A count of `NA` stays `NA`, and R accepts a negative count and `Inf`. The cases `na_count` and `bad_counts` give the evidence.
  - If no gene of `res` is in the network, each feature is 0, and R gives no message. The case `none_in_network` gives the evidence.
- New behavior: The seed count of a gene is the highest count of its rows in `de`. A count that is NaN, infinite, or negative raises `InputError`, also in a row of a gene outside the network.
  - If no gene of `de` is in the network, `string_scores` warns with `SeedMatchWarning`. The message names the first gene of the network as an example.
- Reason: A gene-level result uses the maximum count (`D-003`). A count of NaN, below 0, or infinite is not a count of matches. IDs of another type give scores that look valid, thus the function warns (`D-058`).
- Test: `test_string_scores_take_the_highest_count_of_a_gene_D062`, `test_string_scores_reject_a_bad_count_D062`, and `test_string_scores_warn_without_a_de_gene_in_the_network_D062` in `tests/network/test_string_scores.py`.

## D-063: plot_ecdf labels each p-value

- ID: `D-063`
- R location: `.plot_ecdf` (`R/deseq_fc_ecdf.R:312-324`), called by `deseq_fc_ecdf` (`:133-145`).
- R behavior: The branches compare `p < .01 && p >= .001`, `p > .01`, and `p < .001`. No branch takes p = 0.01, thus R stops with `object 'pvalue_label' not found`.
  - A NaN or NA p-value stops R with `missing value where TRUE/FALSE needed`.
  - A DTS p-value of 20 of 2000 permutations is exactly 0.01.
  - Evidence: `tests/golden/plots/ecdf_labels.tsv`.
- New behavior: The label of p = 0.01 is `0.01`, and NaN gives `NaN`. The other labels are the labels of R.
- Reason: Each result of `ecdf_test` must give a figure.
- Test: `test_plot_ecdf_labels_each_p_value_D063` and `test_plot_ecdf_labels_the_p_value_of_exactly_one_hundredth_D063` in `tests/plot/test_plot_ecdf.py`.

## D-064: plot_ecdf draws only the groups, and its legend counts the values of each line

- ID: `D-064`
- R location: `deseq_fc_ecdf` (`R/deseq_fc_ecdf.R:100-121`), and `.plot_ecdf` (`:326-336`).
- R behavior: A row of `res` in no gene list has the group `NA`. `stat_ecdf` draws a gray line for it, with the legend key `NA`.
  - Evidence: the case `no_group_r` of `tests/golden/plots/ecdf_lines.tsv`.
- R behavior: The count of a legend label is the length of the gene list (`:112`). Thus a gene that `res` does not have, a repeated gene, and a row without a fold change add to the count.
  - Evidence: the notes `legend absent_and_repeated` (`Background: 164`) and `legend na_fold_change` of `tests/golden/plots/ecdf_lines.tsv`.
- R behavior: The labels come in the order of `factor_order`, but the counts come in the order of `gene_lists` (`:112`, `:121`). Thus another order gives each group the count of another group.
  - Evidence: the note `legend factor_order_reversed` (`mer7m8: 161`) of `tests/golden/plots/ecdf_lines.tsv`.
- New behavior: The figure has one line for each key of `groups`, in the order of the keys. The count of a label is the number of finite `log2FoldChange` values of its line.
- Reason: A legend label must describe its line.
- Test: `test_plot_ecdf_draws_only_the_groups_D064`, `test_plot_ecdf_counts_the_values_of_each_line_D064`, and `test_plot_ecdf_follows_the_order_of_the_keys_D064` in `tests/plot/test_plot_ecdf.py`.

## D-065: plot_wasserstein draws the ECDF that the areas use

- ID: `D-065`
- R location: `wass_dist` (`R/wass_dist.R:89-93`).
- R behavior: `stat_ecdf(geom = "line")` adds a point at `-Inf` with 0 and a point at `Inf` with 1. `geom_line` joins them to the values.
  - Evidence: `tests/golden/plots/wasserstein_ecdf.tsv`.
- R behavior: Thus the line goes from 0 at the left edge of the panel to the height of the first value, with a slope. The slope changes with the width of the panel. The areas use 0 below the first value (`yleft = 0`, `:70-71`).
- New behavior: Each line is the interpolated ECDF of `wasserstein_areas`. It is 0 up to the first value, and then it goes up to the height of the first value. Straight lines join the values, and the line is 1 after the last value.
- Reason: The figure must show the functions that the areas integrate.
- Test: `test_plot_wasserstein_draws_the_ecdf_of_the_areas_D065` in `tests/plot/test_plot_wasserstein.py`.

## D-066: plot_kmer_histogram draws one figure with bins of integer counts and the target line

- ID: `D-066`
- R location: `kmer_hist` (`R/kmer_hist.R:30-34`).
- R behavior: `plot(graphics::hist(...))` draws the histogram two times. Page 1 has the title and the label of `kmer_hist`. Page 2 has `Histogram of gene.counts` and `gene.counts`, and `abline` draws the target line on page 2 only.
  - Evidence: the notes `pages` and `title calls` of `tests/golden/plots/kmer_hist.tsv`.
- R behavior: `breaks = 50` gives the breaks of `pretty`, and each bin is closed on the right. The first bin also holds its lower break. Thus for the breaks 0, 2, ..., 88, the first bin holds the counts 0, 1, and 2, and each other bin holds two counts.
- R behavior: on `Test_seqs.fa`, the first bar of R is 3394, but only 2142 k-mers have 0 or 1 genes. For other counts, `pretty` gives more than 50 bins, or bins of a fractional width. For example, the range 1 to 7 gives 60 bins of width 0.1.
- R behavior: a target that is not a k-mer of the counts gives `NA`, for example a lowercase target or a target of another length. R then draws no line and gives no message.
  - Evidence: the target notes of `tests/golden/plots/kmer_hist.tsv`.
- New behavior: `kmer_histogram` gives bins of w integer counts each. w is the first of 1, 2, 5, 10, 20, 50, and so on, that gives at most `max_bins` bins.
- New behavior: `plot_kmer_histogram` draws one Axes with the bars, the title, and the target line. A target that is not a k-mer raises `InputError`.
- Reason: Each bar must hold the same number of integer counts, and the target line must be on the histogram.
- Test: `test_kmer_histogram_bins_the_integer_counts_D066` in `tests/test_kmers.py`. `test_plot_kmer_histogram_draws_the_bins_and_the_target_D066`, `test_plot_kmer_histogram_titles_the_first_page_of_r_D066`, and `test_plot_kmer_histogram_rejects_a_target_that_is_not_a_kmer_D066` in `tests/plot/test_plot_kmers.py`.

## D-067: plot_alignment names the guide line siRNA

- ID: `D-067`
- R location: `plotAlignment` (`R/plotAlignment.R:63`), and `.create_alignment_df` (`:135-137`).
- R behavior: the figure names the guide line `miRNA`, but the data.frame names it `siRNA `.
  - Evidence: `tests/golden/plots/alignment_plot.tsv` and its note `data.frame row names`.
- New behavior: the names are `siRNA `, six spaces, and `target`, as in the data.frame.
- Reason: the guide of SeedMatchR is the guide strand of an siRNA. One thing must have one name.
- Test: `test_plot_alignment_shows_the_r_lines_D067` in `tests/plot/test_plot_alignment.py`.

## D-068: plot_ranges draws each range at the position of its row, with the options

- ID: `D-068`
- R location: `plotIRanges` (`R/plotIRanges.R:3-14`) and `plotRangesDF` (`:16-21`).
- R behavior: `plotRangesDF` draws only the first 6 rows, because of `head()` (`:20`). It ignores `xlim`, `main`, and `col`.
  - Evidence: the case `eight_rows` and the note `options` of `tests/golden/plots/ranges_plot.tsv`.
- R behavior: the y value of a row is `as.integer(rownames(ir))` (`:17-18`). Thus a subset keeps the old positions of its rows, and text row names give NA.
  - Evidence: the cases `subset` and `text_row_names` of `tests/golden/plots/ranges_plot.tsv`.
- R behavior: the bar of a range of one base has the width 0, because `geom_rect` goes from `start` to `end`. `plotIRanges` skips such a range with a warning.
  - Evidence: the case `one_base` and the note `plotIRanges one_base warning` of `tests/golden/plots/ranges_plot.tsv`.
- New behavior: `plot_ranges` draws one bar for each row, at the position of the row in the table. A bar goes from `start - 0.5` to `end + 0.5`, and each option has an effect.
- Reason: the figure must show each range, and a range of one base covers one position.
- Test: `test_plot_ranges_draws_each_row_D068`, `test_plot_ranges_uses_the_options_D068`, `test_plot_ranges_puts_each_row_at_its_position_D068`, and `test_plot_ranges_shows_a_range_of_one_base_D068` in `tests/plot/test_plot_ranges.py`.

## D-069: plot_benchmark colors the curve of one distinct score and of an infinite score

- ID: `D-069`
- R location: `plot.PRROC` and `getColor` of PRROC 1.4, for the curves of `benchmark_preds` (`R/benchmark.R:32-40`).
- R behavior: for one distinct score s, the color range goes from 0.9 s to s + 0.09 s. For s = 0, each color is NA.
  - For s = 0 and for s < 0, `image` stops R with `increasing 'x' and 'y' values expected`. Evidence: the notes `one score 0` and `one score -3` of `tests/golden/plots/benchmark_plot.tsv`.
- R behavior: an infinite score gives the color NA to its threshold, and R stops with `'to' must be a finite number`.
  - Evidence: the note `an Inf score` of `tests/golden/plots/benchmark_plot.tsv`.
- New behavior: the color range comes from the finite thresholds, and an infinite threshold gets the color of the nearest end of the scale.
- New behavior: without two distinct finite thresholds, each segment gets the middle color `#00ff60`, and the color bar goes from c - 0.5 to c + 0.5. If c - 0.5 and c + 0.5 round to c, matplotlib widens the color bar by 10% of |c| on each side (0.9 c to 1.1 c for c above 0, and 1.1 c to 0.9 c for c below 0). For s = 2, R gives `#00FF85` (the case `one_score_2`).
- New behavior: a finite threshold outside -9e305 to 9e305 raises `InputError`. Inside this bound, each position on the scale and each boundary of the color bar is finite.
- Reason: each result of `benchmark_predictions` must give a figure. A seed count of 0 for each gene is a valid input.
- Test: `test_plot_benchmark_colors_one_distinct_score_D069`, `test_plot_benchmark_colors_an_infinite_score_D069`, `test_plot_benchmark_colors_a_negative_infinite_score_D069`, `test_plot_benchmark_rejects_a_threshold_that_is_too_large_D069`, and `test_plot_benchmark_takes_a_threshold_just_inside_the_bound_D069` in `tests/plot/test_plot_benchmark.py`.

## D-070: report_html labels the windows from the guide, and gives no fixed element ID

- ID: `D-070`
- R location: `.make_table` (`R/SeedMatchReport.R:142-196`), called by `SeedMatchReport` (`:135`).
- R behavior: the spanner labels are fixed texts: `Full Guide Strand (g2:g23)` (`:158`), `18-mer (g2:g19)` (`:162`), and `15-mer (g2:g19)` (`:166`).
  - The `15mer` window is g2 to g16 (`R/get_seed.R:121-123`), thus the 15-mer label is wrong. For a guide that does not have 23 bases, the `Full` label is wrong.
  - Evidence: the notes `spanner 15mer` and `spanner Full of a 21-base guide` of `tests/golden/plots/report_table.tsv`.
- R behavior: `gt::gt(..., id = "mygt")` (`:145`) gives each table the element ID `mygt`, and each CSS rule selects `#mygt`. Thus two tables in one HTML page have the same ID and the same rules.
  - Evidence: the notes `id` and `css selector` of `tests/golden/plots/report_table.tsv`.
- New behavior: `seed_window` of `Seed.FULL`, `Seed.MER18`, and `Seed.MER15` for the guide gives the labels. The table has the class `seedmatch-report` and no ID, and each CSS rule selects this class.
- Reason: a label must give the window of the search, and an ID must be unique in a page.
- Test: `test_report_html_labels_the_windows_from_the_guide_D070` and `test_report_html_has_no_fixed_id_D070` in `tests/report/test_report_html.py`.

## D-071: read_fasta skips a byte order mark at the start of the file

- ID: `D-071`
- R location: `Biostrings::readDNAStringSet`. The tests of SeedMatchR call it (`tests/testthat/helper.R:10`).
- R behavior: a FASTA file that starts with a UTF-8 byte order mark stops the read. R probe with R 4.5.3 and Biostrings 2.78.0: `readDNAStringSet` of a file with the bytes `EF BB BF` before `>a` stops with `reading FASTA file ...: ">" expected at beginning of line 1`.
- New behavior: `read_fasta` skips a byte order mark at the start of the file. A U+FEFF after the start is not a DNA letter, thus it raises `InputError`. An incomplete mark is not valid UTF-8, thus it also raises `InputError`.
- Reason: some text editors write the mark at the start of a file. The table readers of the port also skip it at the start of a file, as R does for a table file.
- Test: `test_d071_read_fasta_skips_a_byte_order_mark_at_the_start_of_the_file` in `tests/sequence/test_sequence_set.py`.

## D-072: A custom seed with the name of a canonical seed causes an error

- ID: `D-072`
- R location: `get_seed`. The function compares `seed.name` with the canonical names before it reads `start.pos` and `stop.pos` (`R/get_seed.R:55-58`).
- R behavior: a canonical name gives the canonical seed, and R ignores the window. R probe with R 4.5.3 and Biostrings 2.78.0: `get_seed("UUAUAGAGCAAGAACACUGUUUU", "mer8", start.pos = 2, stop.pos = 5)` gives the `Target.Seq` `CTCTATAA` of `mer8`.
- New behavior: `CustomSeed` raises `InputError` for the name of a canonical seed. The message names the `Seed` member. The comparison uses the case, as in R. Thus `CustomSeed("MER8", 2, 8)` is valid.
- Reason: a result names its rows with the seed name. Thus a reader cannot tell the rows of a custom window with a canonical name from the rows of the canonical seed. R gives no message.
- Test: `test_d072_custom_seed_rejects_the_name_of_a_canonical_seed` in `tests/test_seed.py`.

## D-073: read_deseq2_table rejects a quoted field without its closing quote, and a long field

- ID: `D-073`
- R location: `download_parse_file` (`R/get_data.R:74`), which reads the file with `read.csv`.
- R behavior: R probes with R 4.5.3 show these results:
  - If the last row starts a quoted field and no line closes it, `read.csv` gives no row. This includes each complete row before it. R only warns `incomplete final line found by readTableHeader`.
  - `read.csv` reads a field of 200,000 characters.
- New behavior: `read_deseq2_table` raises `InputError` for a quoted field without its closing quote. The message gives the first line of the row. A field of more than `csv.field_size_limit()` characters raises `InputError` with the line. The default limit is 131,072 characters.
- Reason: a quote without its end loses each row, and R gives only a warning. The `csv` module of Python limits the length of a field. The reader keeps this limit, because a change of the limit applies to the full process. A gene ID or a gene symbol has fewer characters.
- Test: `test_read_deseq2_table_rejects_a_quoted_field_without_its_closing_quote_D073` and `test_read_deseq2_table_names_the_line_of_a_field_above_the_limit_of_the_csv_module_D073` in `tests/test_de.py`.

## D-074: read_fasta reads bzip2 and xz, and it rejects damaged compressed data

- ID: `D-074`
- R location: `Biostrings::readDNAStringSet`. The tests of SeedMatchR call it (`tests/testthat/helper.R:10`).
- R behavior: R probes with R 4.5.3 and Biostrings 2.78.0 show these results:
  - `readDNAStringSet` reads a gzip FASTA file. The first bytes of the file, not the suffix, show the gzip format.
  - For a bzip2 file, R stops with `has unsupported type: bzfile`. For an xz file, R stops with `has unsupported type: xzfile`.
  - For a gzip file that stops after 30 bytes, R gives the bases before the damage, with no error. Thus the record `a` gets the bases `AC`, not `ACGTGG`.
- New behavior: `read_fasta` reads a gzip, bzip2, or xz file. The first bytes of the file show the format. The file can hold more than one compressed stream. Damaged or truncated compressed data raises `InputError` with the path, also in a later stream.
- Reason: the table readers of the port read the same formats, as `read.csv` does. A truncated file gives sequences that are too short, and R gives no message.
- Test: `test_d074_read_fasta_reads_bzip2_and_xz` and `test_d074_read_fasta_names_the_file_of_damaged_compressed_data` in `tests/sequence/test_sequence_set.py`.
