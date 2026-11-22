# Map from R names to Python names

This table maps each SeedMatchR function to its Python replacement. The third cell of a row links to the notes of the row, after the table. A test makes sure that each exported R function has a row. If a column gets a new name or a value gets a new type, the notes record the change. Such a change is not a divergence. A change of behavior goes in `docs/divergences.md`.

| R function | Python | Notes |
| --- | --- | --- |
| `Biostrings::readDNAStringSet` | `seedmatch.sequence.read_fasta` | [Notes](#biostringsreaddnastringset) |
| `Biostrings::reverseComplement` | `seedmatch.sequence.reverse_complement` | [Notes](#biostringsreversecomplement) |
| `get_seed`, `.get_canonical_seed`, `.get_custom_sequence` | `seedmatch.Guide.target` | [Notes](#get_seed-get_canonical_seed-get_custom_sequence) |
| `Biostrings::vcountPattern` | `seedmatch.count_matches` | [Notes](#biostringsvcountpattern) |
| `Biostrings::vmatchPattern` | `seedmatch.find_matches` | [Notes](#biostringsvmatchpattern) |
| `vmatchPattern2` | `seedmatch.find_matches` | [Notes](#vmatchpattern2) |
| `.normalize_fixed_arg`, `.get_algo`, `.valid_algos` | `seedmatch.SearchOptions`, `seedmatch.search.MAX_PATTERN_LENGTH` | [Notes](#normalize_fixed_arg-get_algo-valid_algos) |
| `filter_res` | `seedmatch.filter_de` | [Notes](#filter_res) |
| `download_parse_file` | `seedmatch.read_deseq2_table` | [Notes](#download_parse_file) |
| `get_example_data` | `seedmatch-fetch example sirna`, `seedmatch-fetch mirdb`, `seedmatch.read_deseq2_table`, `seedmatch.read_mirdb` | [Notes](#get_example_data) |
| `get_example_data` with `example.type = "mirna"` (the table), `load_example_data("mirna")` | `seedmatch.read_mirdb` | [Notes](#get_example_data-with-exampletype--mirna-the-table-load_example_datamirna) |
| `load_example_data` with `example.type = "sirna"` | `seedmatch.read_deseq2_table` | [Notes](#load_example_data-with-exampletype--sirna) |
| `SeedMatchR` with `res.format = "DESEQ2"` | `seedmatch.annotate_seed_counts` | [Notes](#seedmatchr-with-resformat--deseq2) |
| `.get_counts_column` | `seedmatch.annotate_seed_counts` | [Notes](#get_counts_column) |
| `full_search` with `res_format = "counts"` | `seedmatch.full_search` | [Notes](#full_search-with-res_format--counts) |
| `full_search_ranges` | `seedmatch.full_search_ranges` | [Notes](#full_search_ranges) |
| `.run_scanning` | `seedmatch.full_search_ranges` | [Notes](#run_scanning) |
| `.get_best_match` | `seedmatch.best_match` | [Notes](#get_best_match) |
| `.get_match_seq`, `.get_alignment` | no counterpart | [Notes](#get_match_seq-get_alignment) |
| `SeedMatchR` with `res.format = "data.frame"` | `seedmatch.find_matches` | [Notes](#seedmatchr-with-resformat--dataframe) |
| `.get_match_df` | `seedmatch.find_matches` | [Notes](#get_match_df) |
| `SeedMatchReport` | `seedmatch.seed_match_report` | [Notes](#seedmatchreport) |
| `.make_table` | `seedmatch.report_html` | [Notes](#make_table) |
| `load_annotations` | `seedmatch.Reference.sequences` | [Notes](#load_annotations) |
| `load_species_anno_db` | `seedmatch.Reference.from_paths`, `seedmatch.Reference.from_directory`, `seedmatch-fetch reference` | [Notes](#load_species_anno_db) |
| `build_annotation_filter` | `seedmatch.TranscriptFilter` | [Notes](#build_annotation_filter) |
| `get_3UTR_granges` | `seedmatch.Reference.features` | [Notes](#get_3utr_granges) |
| `.get_longest_utr` | `seedmatch.PerGene.LONGEST` | [Notes](#get_longest_utr) |
| `.reduce_features` | `seedmatch.PerGene.UNION` | [Notes](#reduce_features) |
| `ensembldb::transcripts` | `seedmatch.Reference.transcripts` | [Notes](#ensembldbtranscripts) |
| `GenomicFeatures::extractTranscriptSeqs` | `seedmatch.annotation.extract.extract_sequences` | [Notes](#genomicfeaturesextracttranscriptseqs) |
| `rtracklayer::TwoBitFile` | `seedmatch.annotation.TwoBitGenome`, `seedmatch.annotation.open_genome` | [Notes](#rtracklayertwobitfile) |
| `ecdf_stat_test` | `seedmatch.ecdf_test` | [Notes](#ecdf_stat_test) |
| `deseq_fc_ecdf` | `seedmatch.ecdf_test`, `seedmatch.assign_groups`, `seedmatch.plot.plot_ecdf` | [Notes](#deseq_fc_ecdf) |
| `.check_gene_list_overlap` | `seedmatch.assign_groups` | [Notes](#check_gene_list_overlap) |
| `.assign_groups` | `seedmatch.assign_groups` | [Notes](#assign_groups) |
| `.plot_ecdf` | `seedmatch.plot.plot_ecdf` | [Notes](#plot_ecdf) |
| `SeedMatchR.theme`, `SeedMatchR.palette` | `seedmatch.plot.PALETTE` | [Notes](#seedmatchrtheme-seedmatchrpalette) |
| `wass_dist` | `seedmatch.wasserstein_areas`, `seedmatch.wasserstein_difference`, `seedmatch.plot.plot_wasserstein` | [Notes](#wass_dist) |
| `stats::ks.test` | `seedmatch.stats.ks_test` | [Notes](#statskstest) |
| `stats::wilcox.test` | `seedmatch.stats.wilcoxon_test` | [Notes](#statswilcoxtest) |
| `twosamples::dts_test` | `seedmatch.stats.dts_test` | [Notes](#twosamplesdts_test) |
| `PRROC::pr.curve` (with `weights.class0` and `curve = TRUE`) | `seedmatch.stats.pr_curve` | [Notes](#prrocprcurve-with-weightsclass0-and-curve--true) |
| `PRROC::roc.curve` (with the same arguments) | `seedmatch.stats.roc_curve` | [Notes](#prrocroccurve-with-the-same-arguments) |
| `benchmark_preds` | `seedmatch.benchmark_predictions`, `seedmatch.plot.plot_benchmark` | [Notes](#benchmark_preds) |
| `plot.PRROC` (PRROC 1.4, for the curves of `benchmark_preds`) | `seedmatch.plot.plot_benchmark` | [Notes](#plotprroc-prroc-14-for-the-curves-of-benchmark_preds) |
| `sequence_kmer_counts` | `seedmatch.kmer_counts`, `seedmatch.KmerCounts.to_frame` | [Notes](#sequence_kmer_counts) |
| `kmer_hist` | `seedmatch.kmer_gene_counts`, `seedmatch.kmer_histogram`, `seedmatch.plot.plot_kmer_histogram` | [Notes](#kmer_hist) |
| `plotAlignment` with `output.type = "data.frame"` | `seedmatch.alignment_text` | [Notes](#plotalignment-with-outputtype--dataframe) |
| `plotAlignment` with `output.type = "ggplot"` | `seedmatch.plot.plot_alignment` | [Notes](#plotalignment-with-outputtype--ggplot) |
| `plotIRanges`, `plotRangesDF` | `seedmatch.plot.plot_ranges` | [Notes](#plotiranges-plotrangesdf) |
| `.matchStrings`, `.make_space`, `.create_alignment_df` | `seedmatch.alignment_text` | [Notes](#matchstrings-make_space-create_alignment_df) |
| `match_distances`, `.get_centers`, `.get_fivePrimeEnd`, `.get_threePrimeEnd` | `seedmatch.match_distances` | [Notes](#match_distances-get_centers-get_fiveprimeend-get_threeprimeend) |
| `generate_e_value` | `seedmatch.e_values`, `seedmatch.e_value_null`, `seedmatch.wobble_bulge_null` | [Notes](#generate_e_value) |
| `get_network_matches` (the TRRUST URLs) | `seedmatch-fetch trrust` | [Notes](#get_network_matches-the-trrust-urls) |
| `get_network_matches` | `seedmatch.read_trrust`, `seedmatch.trrust_scores` | [Notes](#get_network_matches) |
| `message_passing`, `as_dataframe_mp_matches` | `seedmatch.trrust_scores` | [Notes](#message_passing-as_dataframe_mp_matches) |
| `get_stringPPI_to_stringGGI` (the STRING download) | `seedmatch-fetch string` | [Notes](#get_stringppi_to_stringggi-the-string-download) |
| `get_stringPPI_to_stringGGI` | `seedmatch.read_string_links`, `seedmatch.string_gene_network` | [Notes](#get_stringppi_to_stringggi) |
| `create_graph` | `seedmatch.string_scores` | [Notes](#create_graph) |
| `run_wobble_bulges_pipeline` | `seedmatch.wobble_bulge_matches`, `seedmatch.annotate_wobble_bulge_counts` | [Notes](#run_wobble_bulges_pipeline) |
| `.get_all_combination_w_b`, `.find_factorial_internal`, `.get_target_bulged_guides`, `.sequence_with_tb`, `check_bulges_nt` | `seedmatch.wobble_bulge_patterns` | [Notes](#get_all_combination_w_b-find_factorial_internal-get_target_bulged_guides-sequence_with_tb-check_bulges_nt) |
| `get_compensatory_sequences` | `seedmatch.compensatory_sites` | [Notes](#get_compensatory_sequences) |
| `relative_position_seq`, `.gc_content` | the columns `guide_start`, `guide_end`, and `gc_content` of `seedmatch.compensatory_sites` | [Notes](#relative_position_seq-gc_content) |
| `.update_sequence_at` | no counterpart | [Notes](#update_sequence_at) |

## Notes

### `Biostrings::readDNAStringSet`

- Returns a `SequenceSet`.
- The ID of a record is the header text up to the first white space (`D-016`).
- The reader changes lowercase letters to uppercase.
- A letter that is not an IUPAC DNA letter raises `InputError`.
- R accepts `-`, `+`, and `.`, but the port does not (`D-017`).
- The reader skips a byte order mark at the start of the file, but R stops (`D-071`).
- As in R, the reader reads a gzip file, and the first bytes of the file show the format.
- The reader also reads bzip2 and xz, but R stops. Damaged or truncated compressed data raises `InputError`, but R reads a truncated gzip file up to the damage (`D-074`).

### `Biostrings::reverseComplement`

- Takes and returns a `str`.
- The result has uppercase letters.
- Each IUPAC DNA letter has its complement, so `YYAYAGAG` gives `CTCTRTRR`, as in R.
- A letter that is not an IUPAC DNA letter raises `InputError` (`D-017`).

### `get_seed`, `.get_canonical_seed`, `.get_custom_sequence`

- `Guide(sequence)` holds the guide as uppercase RNA.
- `Guide.target(seed, wobble=False)` returns a `Target` with the fields `name`, `sequence`, and `seed_rna`.
- `sequence` is `Target.Seq` of R, and `seed_rna` is `Seed.Seq.RNA`.
- A `Seed` member replaces a canonical `seed.name`.
- A `CustomSeed(name, start, stop)` replaces `seed.name` with `start.pos` and `stop.pos`.
- `wobble` replaces `allow_wobbles`.
- The wobble target has `R` at each `U` of the guide, and only a search with `fixed_pattern=False` matches an `A` or a `G` there (`D-026`).
- A `seed` that is not a `Seed` or a `CustomSeed`, for example the text `"mer7m8"`, raises `InputError`.
- A guide that is too short for the seed raises `InputError` that names the seed and the length.
- R stops with `nchar(guide.seq) >= 8 is not TRUE` or `stop.pos <= nchar(guide.seq) is not TRUE`.
- A lowercase guide gets the wobble (`D-018`).
- The guide takes only `A`, `C`, `G`, and `U` (`D-019`).
- A `CustomSeed` with `start` after `stop` raises `InputError` (`D-020`).
- `Seed(name)` raises `InputError` for a name that is not the name of a canonical seed.
- A `CustomSeed` with the name of a canonical seed raises `InputError`, but R uses the canonical seed and ignores the window (`D-072`).

### `Biostrings::vcountPattern`

- Takes a pattern (a `str` or a `Target`), a `SequenceSet`, and a `SearchOptions`.
- Returns an `int64` array with one count for each sequence.
- `max.mismatch`, `min.mismatch`, and `with.indels` become `max_mismatch`, `min_mismatch`, and `indels` of `SearchOptions`.
- The two entries of `fixed` become `fixed_pattern` and `fixed_subject`.
- The default of `indels` is `False`, as in Biostrings, but `SeedMatchR` has the default `with.indels = TRUE`.
- An empty pattern raises `InputError`, and so does a letter that is not an IUPAC DNA letter (`D-017`).

### `Biostrings::vmatchPattern`

- Returns a `pandas.DataFrame`, not an `MIndex`.
- One row is one match.
- The columns are `seq_id`, `gene_id`, `start`, `end`, and `width`.
- The rows come in the order of the sequences.
- A mismatch match can start before 1 or end after the sequence length, as in R.
- The options are the options of `count_matches`.
- `vmatchPattern` stops for a `DNAStringSet` with `with.indels = TRUE` and `max.mismatch > 0`, but `find_matches` does not (`D-001`).
- `seedmatch.search.match_ranges` gives the same matches as three `int64` arrays: the index of the sequence in the set, the start, and the end.

### `vmatchPattern2`

- The R data frame has the columns `start`, `end`, `width`, and `names`.
- `seq_id` replaces `names`, and `gene_id` is new.
- A sequence without a match has no row, as in R.
- The ranges obey `indels` and `min_mismatch` (`D-001`).

### `.normalize_fixed_arg`, `.get_algo`, `.valid_algos`

- These helpers make sure that the arguments of `vmatchPattern2` are valid (`R/vmatchPattern2.R:29-33`).
- `.normalize_fixed_arg` takes `fixed` as one or two flags or as names, and `SearchOptions` has the two flags `fixed_pattern` and `fixed_subject`.
- For an empty pattern or a pattern of more than 20000 letters, `.valid_algos` stops (`:131-134`), and `count_matches` and `find_matches` raise `InputError`.
- For `min.mismatch` above 0 with indels and `max.mismatch` above 0, `.valid_algos` stops (`:135-137`).
- `SearchOptions` raises `InputError` for `min_mismatch` above 0 with `indels` on.
- `vmatchPattern2` gets an algorithm from `.get_algo`, but it gives `"indels"` to the C code (`:43`, `D-001`).

### `filter_res`

- Takes a `DataFrame` and returns a new `DataFrame` with a new index.
- `res` becomes `de`.
- The cutoffs `fdr_cutoff`, `fc_cutoff`, and `basemean_cutoff` have the R names and the R defaults, but they are keyword-only.
- The result keeps the order of the rows.
- A row with a missing `padj`, `pvalue`, `log2FoldChange`, or `baseMean` is removed.
- A table without a required column raises `InputError` that lists the missing columns in sorted order.
- R stops with `"padj" %in% colnames(res) is not TRUE`.
- A filter column with a dtype that is not numeric raises `InputError` that names the column and the dtype.
- A cutoff is an int or a float.
- A bool, text, `None`, and NaN raise `InputError` (`D-025`).

### `download_parse_file`

- Reads a table from a path, because the core modules do not use the network.
- A URL raises `InputError`.
- The file can be gzip, bzip2, or xz.
- As in `read.csv`, the first bytes of the file show the compression.
- The function replaces the header with the DESeq2 column names `gene_id`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`, and `symbol` (`R/get_data.R:74-76`).
- `gene_id` and `symbol` are text, and the other columns have the type `float64`.
- R reads a column of digits as numbers, but the port reads `gene_id` as text.
- As in `read.csv`, only `NA` is missing in a text column, and `NA` or an empty field is missing in a number column.
- A table that does not have exactly 8 columns raises `InputError` that gives the count found (`D-021`). These tables also raise `InputError`: an empty file, a row with more fields than the header, and a number column with a value that is not a number.
- As in `read.csv`, a blank line is no row, and after the header, a line `""` is also no row. A line of spaces is a row.
- A quoted field without its closing quote raises `InputError`, but R reads no row and only warns. A field above the limit of the `csv` module raises `InputError`, but R reads it (`D-073`).

### `get_example_data`

- `example.type = "sirna"` becomes `seedmatch-fetch example sirna DEST`, and `"mirna"` becomes `seedmatch-fetch mirdb DEST`.
- The command downloads the files, because the core modules do not use the network.
- R writes to `tools::R_user_dir("SeedMatchR", "data")` (`R/get_data.R:16`), but the command writes only to `DEST`.
- The `sirna` files are the GEO URLs of `R/get_data.R:23`, `:27`, and `:31`.
- R parses each file with `download_parse_file` and saves the tables in `Schlegel_2022.Rdata`.
- The command keeps the files as they are, and `read_deseq2_table` reads each file.
- The `mirdb` file is the URL of `R/get_data.R:39`.
- R saves it as `miRDB_v6_prediction_results.txt.gz` (`R/get_data.R:37`), but the command keeps the name of the URL, `miRDB_v6.0_prediction_result.txt.gz`.
- `seedmatch.read_mirdb` reads the miRDB table (see the next row).
- R does no check of the size or the digest, but the command makes sure that the size and the digest agree with the registry.

### `get_example_data` with `example.type = "mirna"` (the table), `load_example_data("mirna")`

- R downloads the file, reads it with `read.table`, and saves `mirdb.Rdata` (`R/get_data.R:39-54`).
- `load_example_data("mirna")` loads this file (`R/get_data.R:101-104`).
- The port reads the file of `seedmatch-fetch mirdb` and saves nothing.
- `miRDB.ID`, `target.REFSEQ.ID`, `target.score`, `species`, and `miRNA.name` become `mirdb_id`, `target_refseq_id`, `target_score`, `species`, and `mirna_name`.
- The keyword `species` is new. It keeps the rows of one species code, and a code without a row raises `InputError`.
- A score that is not a finite number raises `InputError`, and the ID fields keep `#`, quote characters, and `NA` as text (`D-049`).

### `load_example_data` with `example.type = "sirna"`

- `"sirna"` is the default.
- R loads the DESeq2 tables that `get_example_data("sirna")` saves in `Schlegel_2022.Rdata`, and it returns them as a list (`R/get_data.R:96-99`).
- `seedmatch-fetch example sirna DEST` downloads the files, and `read_deseq2_table` reads each file.
- The port saves no `.Rdata` file.

### `SeedMatchR` with `res.format = "DESEQ2"`

- `res` becomes `de`, and `seqs` becomes `sequences`.
- A `Target` replaces `sequence` and `seed.name` (`Guide.target`).
- `col.name` becomes `column`, and the default is the name of the target.
- `max.mismatch`, `with.indels`, and `fixed` become `options`.
- The default is `SearchOptions(indels=True)`, the search with the R defaults.
- `shared_genes` keeps its meaning.
- With `options.fixed_pattern`, a target with an ambiguity letter, for example a wobble target, raises `InputError` (`D-026`).
- The result keeps the rows of the genes with a sequence, in the order of the input, and has a new index.
- The function returns a new `DataFrame` and does not change `de`.
- The count column has the type `int64`, but R gives a `double`.
- The count of a gene is the maximum over its sequences (`D-003`).
- A column with the same name is replaced (`D-015`).
- A `de` without a `gene_id` column raises `InputError`.
- A repeated `gene_id` column also raises `InputError`, but R uses the first column if `shared_genes` is `FALSE`.
- A repeated column with the name of the count column raises `InputError` too (`D-051`).
- A `gene_id` column that is not text raises `InputError`, because no gene would match.
- `column` cannot be `gene_id`.

### `.get_counts_column`

- `SeedMatchR` with `res.format = "DESEQ2"` calls this helper (`R/SeedMatchR.R:95-104`), thus the row of that mode gives the arguments.
- A direct call has the defaults `with.indels = F` and `fixed = F` (`R/SeedMatchR.R:190-191`).
- `options=SearchOptions(fixed_pattern=False, fixed_subject=False)` gives these defaults.
- The default of `annotate_seed_counts` is the default of `SeedMatchR`.
- `seed.name` and `col.name` become `column`.
- The count of a gene is the maximum over its sequences (`D-003`), and the counts replace a column with the same name (`D-015`).

### `full_search` with `res_format = "counts"`

- `guide.seq` becomes a `Guide`, `seqs` becomes `sequences`, and `group.name` becomes `group`.
- The other arguments become the fields of `FullSearchOptions`, and the defaults are the R defaults.
- `seqs.to.search` becomes `seeds` (a tuple of `Seed`).
- The default is `DEFAULT_SEEDS`, the seeds of R in the order of R.
- `nonseed_edit_distances` applies to the seeds of `NONSEED_SEEDS`: `Seed.FULL`, `Seed.MER18`, and `Seed.MER15`, as in `R/full_search.R:44`.
- `nonseed.edit.distances` and `seed.edit.distances` become `nonseed_edit_distances` and `seed_edit_distances`.
- `indel.bool` becomes `indels`, and `allow_wobbles` becomes `wobble`.
- `fixed` becomes `fixed_pattern` and `fixed_subject`. These are off by default, unlike in `SearchOptions`.
- `wobble=True` with `fixed_pattern=True` raises `InputError` (`D-026`).
- R passes `with.indel` (a partial name) to `vcountPattern` (`R/full_search.R:61`), and Biostrings reads it as `with.indels`.
- The keyword `workers` sets the width of the thread pool, and the result does not depend on it.
- The count columns keep the R names `<seed>.mm<distance>` and have the type `int64`.
- `seqs`, `Group`, and `best.match` become `seq_id`, `group`, and `best_match`.
- A new column, `gene_id`, follows `seq_id`.
- `FullSearchOptions` raises `InputError` for an invalid option, for a duplicate value (`D-023`), and for a set instead of a tuple or a list.
- A guide that is too short for a seed raises `InputError` that names the seed and the length, before any search.
- `res_format = "iranges"` has no counterpart.
- `full_search_ranges` does not replace it, because the functions differ.
- The `iranges` mode of R uses `allow_wobbles`. It takes the separate lists `nonseed.edit.distances` and `seed.edit.distances` and the flag `indel.bool`. It names a search `<seed>.mm<distance>`, without an indel flag. Its columns are `start`, `end`, `width`, `names`, `query`, and `Group`. Its rows go by seed first, then by distance.
- `full_search_ranges` has no wobble option. It takes a list `edit_distances` with the limit `max_seed_edit_distance`, and a list `indel_options`. Its rows go by distance first, then by seed, then by flag.
- To get the `iranges` result, call `find_matches` for each seed and distance, with a target from `Guide.target`.

### `full_search_ranges`

- The arguments `indel.bool.opt`, `edit.dist.list`, `seqs.to.search`, `fixed.opt`, and `max.seed.edit.dist` become `indel_options`, `edit_distances`, `seeds`, `fixed_pattern` and `fixed_subject`, and `max_seed_edit_distance` of `FullSearchRangesOptions`.
- The rows follow the `expand.grid` order of R.
- The ranges obey `indels` (`D-001`).
- `seqnames`, `TargetSeq`, `EditDistance`, `AllowIndels`, and `Group` become `seq_id`, `target_seq`, `edit_distance`, `indels`, and `group`.
- `indels` has the type `bool`, and `edit_distance` has the type `int64`.
- The column `names` repeats `seqnames`, thus the port does not keep it.
- A new column, `gene_id`, follows `seq_id`.
- The column `query` has the R values.
- The function gives a result for any number of searches with a match (`D-022`).
- `workers` and the checks of the options are as in `full_search`.
- `max_seed_edit_distance` must be an integer of 0 or more (`D-023`).

### `.run_scanning`

- `full_search_ranges` calls this helper for each search (`R/full_search.R:170-179`). It names a search `<seed>.mm<distance>.indel<flag>`, and it adds the columns `query`, `TargetSeq`, `seqnames`, `EditDistance`, and `AllowIndels` (`:191-223`).
- The row of `full_search_ranges` gives the Python columns.

### `.get_best_match`

- Takes the table of counts and returns a `Series` named `best_match` with the index of the table.
- The ranking and the labels are the R ones.
- A column that the table lacks is skipped (`D-012`).

### `.get_match_seq`, `.get_alignment`

- No SeedMatchR function calls them.
- The calls of `get_match_seq` are comments (`R/SeedMatchR.R:276`, `R/full_search.R:217`).
- `.get_alignment` reads the column `MatchString`, and no code makes this column (`R/full_search.R:231-238`).
- `alignment_text` gives the alignment of a match.

### `SeedMatchR` with `res.format = "data.frame"`

- `target.seq`, `seqs`, `max.mismatch`, `with.indels`, and `fixed` become the arguments of `find_matches` (`Guide.target` makes the target).
- The search with the R defaults uses `SearchOptions(indels=True)`.
- The columns `seqnames`, `start`, `end`, and `width` become `seq_id`, `gene_id`, `start`, `end`, and `width`.
- `match.df` has no argument: `pandas.concat([new, old])` gives the rows of R, with the new rows first.
- The columns of `sirna.name` (`sirna`, `query`, `TargetSeq`, `EditDistance`, and `AllowIndels`) are not made.
- `res.format = "granges"` and `"iranges"` have no counterpart.

### `.get_match_df`

- `SeedMatchR` with `res.format = "data.frame"` or `"granges"` calls this helper (`R/SeedMatchR.R:110-117`), thus the row of that mode gives the arguments.
- A direct call has the defaults `with.indels = TRUE` and `fixed = TRUE` (`R/SeedMatchR.R:257-258`), thus `SearchOptions(indels=True)`.
- `seed.name` and `sirna.name` have no counterpart.
- R uses them only for the columns `sirna`, `query`, `TargetSeq`, `EditDistance`, and `AllowIndels` (`:272-279`), and the port does not make these columns.
- The ranges obey `indels` (`D-001`).

### `SeedMatchReport`

- `res` becomes `de`, `seqs` becomes `sequences`, and `guide.seq` becomes a `Guide`.
- `group.name` becomes `group`.
- The other arguments become the fields of `ReportOptions`, and the defaults are the R defaults.
- `indel.bool` becomes `indels`, and `allow_wobbles` becomes `wobble`.
- `log2fc.cut` and `padj.cut` become `log2fc_cut` and `padj_cut`, and `seed.edit.distances` becomes `seed_edit_distances` (a tuple).
- `fixed` becomes `fixed_pattern` and `fixed_subject`.
- The default of `fixed_pattern` is `True`, as in R, but `wobble=True` is permitted only with `fixed_pattern=False` (`D-026`).
- As in R, the edit distances of `Full`, `18mer`, and `15mer` are not options, and they are 0 to 4.
- `ReportOptions` raises `InputError` for an invalid field, as `FullSearchOptions` does.
- A cutoff is an int or a float.
- A bool, text, and NaN are not valid (`D-024`).
- The function returns a `SeedMatchReport` (a frozen dataclass) with the fields `report` and `matches`, not a list.
- `report_html(report, guide)` gives the HTML of the `gt` table (`table`) (`D-070`).
- In `report`, `Category` and `Group` become `category` and `group`.
- The category columns and `Total` have the type `float64`, because one column holds counts and a percentage, as in R.
- `REPORT_CATEGORIES` holds the category names in the R order (`expected.names`).
- In `matches`, `seqs`, `Group`, `best.match`, and `IN.DESEQ2` become `seq_id`, `group`, `best_match`, and `detected`.
- `detected` has the type `bool`, not `Detected` or `Undetected`.
- A new column, `gene_id`, follows `seq_id`.
- `matches` keeps one row for each sequence, in the order of the set.
- The counts of `report` are gene counts, and the best match of a gene is the best over its sequences (`D-003`).
- For a gene with more than one row in `de`, the first row counts, as in R.
- `de` must have the columns `gene_id`, `log2FoldChange`, `padj`, and `baseMean`.
- A missing column raises `InputError` that lists the missing columns in sorted order.
- A `gene_id` column that is not text raises `InputError`.
- A value column that is not numeric also raises `InputError` (`D-024`).

### `.make_table`

- R returns a `gt_tbl` object, and the port returns HTML text.
- The columns, the formats, the title, and the subtitle are the same.
- The labels of the windows come from the guide (`D-070`).
- The color of a percent cell is a linear blend from white to red, with black text.
- `gt` blends in another color space and gives white text on a dark cell. 0 and 100 give the colors of `gt` (`D-070`).
- The port writes the minus sign as `&#8722;` and omits the space that `gt` writes before `Inf` in a count column that also has `-Inf`.
- `gt` writes `NA` for NA and `NaN` for NaN, and the port writes `NaN` for each missing value.

### `load_annotations`

- `load_annotations` returns a list.
- `$seqs` becomes `Reference.sequences`, a `SequenceSet`.
- `$features` becomes `Reference.features`, a `Features` object.
- `$gtf` becomes `Reference.transcripts`.
- `$txdb` and `$dna` become `Reference.annotation` and `Reference.genome`.
- `$filter` becomes a `TranscriptFilter`.
- `$tx.names` has no counterpart, because it is `NULL` in R (`D-009`).
- `reference.name` has no counterpart, because the port reads local files (see `load_species_anno_db`).
- `feature.type` becomes `feature`, a `FeatureType`: `UTR3`, `UTR5`, `EXONS`, or `CDS`.
- The filter arguments become the fields of a `TranscriptFilter` (see `build_annotation_filter`), and `Reference.sequences` gives the default filter of `load_annotations`.
- `longest.utr` and `reduce.features` become `per_gene`.
- `PerGene` replaces the flags `longest.utr` and `reduce.features`, and it applies to each feature type (`D-029`).
- `min.feature.width` becomes `min_width`, with the default 8.
- A negative `min_width` raises `InputError`.
- `return_gene_name` has no counterpart.
- Each sequence keeps its transcript ID in `ids` and its gene ID in `gene_ids` (`D-009`).
- With `PerGene.UNION`, `ids` holds the gene IDs, and R gives the name `NA` to each sequence (`D-002`).
- A transcript without the feature has no range, but R gives it a range of width 0.
- The sequence is empty in R and in the port, and `min_width` removes it.

### `load_species_anno_db`

- R reads the EnsDb and the 2bit file from AnnotationHub, thus it uses the network.
- The port reads the files that you name.
- `seedmatch-fetch reference BUILD DEST` downloads the 2bit file and the EnsDb file of a build, and `from_directory(DEST)` opens them.
- `BUILD` is a value of `reference` in `R/load_annotations.R:152`, for example `rnor7`.
- The registry pins the size and the SHA-256 digest of each file.
- Use one `DEST` for each build, because `from_directory` accepts only one EnsDb file.
- `reference.name` becomes the arguments `ensdb` and `genome` of `from_paths`.
- R has no ID for the choice `hg38-old` (`D-009`), and the port has no such choice.
- `from_directory(path)` finds the one `.sqlite` file and the one genome file of a directory.
- The genome file ends with `.2bit` or `.fa.gz`.
- An absent file, or more than one candidate, raises `InputError` that names the directory and the candidates.
- The genome is a 2bit file or a FASTA file, and the first bytes of the file show the format.
- A `Reference` raises `InputError` if no contig of the annotation is in the genome, or if a contig has a different length in the genome.
- `filter` has no argument.
- `Reference.transcripts(transcript_filter)` replaces `$gtf` for all transcripts of the filter.
- A `Reference` holds open files, thus `close()` or the `with` statement releases them.

### `build_annotation_filter`

- `TranscriptFilter` is a frozen dataclass.
- `standard.chroms`, `remove.na.rows`, `protein.coding`, and `canonical` become `standard_chromosomes`, `require_tx_id`, `protein_coding`, and `canonical`.
- Each is a `bool`.
- `transcript.support` becomes `support_level`, an integer from 1 to 5.
- `gene_id`, `tx_id`, `symbol`, and `entrez_id` become `gene_ids`, `tx_ids`, `symbols`, and `entrez_ids`.
- Each takes a sequence.
- A `str` or an empty sequence raises `InputError`.
- `tx_biotypes` and `gene_biotypes` replace the free `add.filter`.
- The defaults are the defaults of `load_annotations`. Thus `canonical` is `True`, but the default of `build_annotation_filter` is `FALSE`.
- A `support_level` below 1 and an Entrez ID below 1 raise `InputError` (`D-031`).
- The filter is a condition of the SQL query on the EnsDb.
- `message` has no counterpart.

### `get_3UTR_granges`

- `Reference.features(FeatureType.UTR3, transcript_filter, per_gene=...)` returns a `Features` object, not a `GRangesList`.
- Feature i has the ranges `starts[k]` to `ends[k]` for k in `offsets[i]:offsets[i + 1]`.
- The coordinates are 1-based and inclusive.
- `Features.to_frame()` gives one row for each range.
- `txdb` and `filter` become the `Reference` and the `transcript_filter`.
- `longest.utr` and `reduce.features` become `per_gene`.
- The ranges of a feature come in transcript order, thus they descend on the minus strand.
- A transcript without a 3' UTR has no range, but R gives it a range of width 0.
- `Reference.features` also takes `FeatureType.UTR5`, `EXONS`, and `CDS`, which replace `fiveUTRsByTranscript`, `exonsBy`, and `cdsBy`.

### `.get_longest_utr`

- `Reference.features` and `Reference.sequences` take `per_gene=PerGene.LONGEST`.
- The result has one transcript for each gene, in the order of `gene_id`. This transcript has the feature of the maximum width, for example the longest 3' UTR for `UTR3`.
- A tie goes to the first transcript in the order of `Reference.transcripts`.
- `transcriptLengths` has no counterpart.
- Only a transcript with the feature is a candidate, and a gene without such a transcript has no feature (`D-028`).
- `LONGEST` applies to each feature type (`D-029`).

### `.reduce_features`

- `per_gene=PerGene.UNION` merges the ranges of the transcripts of a gene that overlap or touch.
- The result has one feature for each gene, in the order of `gene_id`.
- The `ids` and the `gene_ids` of the features are the gene IDs (`D-002`).
- The ranges come in transcript order: from the lowest position to the highest on the plus strand, and from the highest to the lowest on the minus strand (`D-027`).
- `GenomicRanges::reduce` has no counterpart.
- `UNION` applies to each feature type (`D-029`).

### `ensembldb::transcripts`

- Returns a `pandas.DataFrame` with one row for each transcript that passes the filter, not a `GRanges`.
- The columns are `tx_id`, `gene_id`, `gene_name`, `contig`, `start`, `end`, `strand`, `tx_biotype`, `gene_biotype`, `cds_start`, `cds_end`, `support_level`, and `canonical`.
- `seqnames`, `tx_cds_seq_start`, `tx_cds_seq_end`, `tx_support_level`, and `tx_is_canonical` become `contig`, `cds_start`, `cds_end`, `support_level`, and `canonical`.
- `strand` is `+` or `-`.
- A column that the file lacks has no value, and neither has the CDS of a transcript without a CDS.
- The rows come in the order of `contig` and `start`.
- A tie has a fixed order (`D-030`).
- `Reference.transcripts` replaces `$gtf` of `load_annotations`.
- R keeps in `$gtf` only the transcripts that have a sequence, but `Reference.transcripts` gives all transcripts of the filter.

### `GenomicFeatures::extractTranscriptSeqs`

- `Reference.sequences` calls `extract_sequences(features, genome)`. It returns a `SequenceSet`, not a `DNAStringSet`. It joins the ranges of a feature in list order.
- On the minus strand, it uses the reverse complement of each range.
- `Genome.resolve` maps each contig of the annotation to a contig of the genome. It tries the same name, the name without `chr`, the name with `chr`, and the other name of the mitochondrial contig.
- A contig that the genome lacks raises `InputError`.
- The sequence has uppercase letters, also for the soft-masked bases.
- A feature without a range has an empty sequence.
- R gives a transcript without the feature a range of width 0, and the sequence is empty in R and in the port.

### `rtracklayer::TwoBitFile`

- `TwoBitGenome(path)` reads a 2bit file of version 0 or 1, in each byte order.
- `rtracklayer` reads only version 0.
- The reader applies the N blocks, and it decodes only the bytes of the requested range.
- `seedmatch.annotation.Genome` is the type of each genome reader.
- `Genome.codes(contig, start, end)` gives the masks of a range with 1-based inclusive coordinates.
- `Genome.sequence` gives uppercase text.
- `Genome.contigs` maps each contig to its length, as `seqlengths` does.
- `open_genome(path)` reads the first bytes of the file, and it gives a `TwoBitGenome` or a `FastaGenome`.
- `seedmatch.annotation.FastaGenome` reads a plain FASTA file with its `.fai`, or a BGZF file with its `.fai` and `.gzi`.
- A URL raises `InputError`.
- The readers write no file.

### `ecdf_stat_test`

- `res` becomes `de`.
- `target_list` and `background_list` become `target` and `background`.
- Each is an iterable of gene IDs, for example a `Series`, a `list`, or a `set`.
- A `str` raises `InputError`, because it is one ID, and so does a gene ID that is not text.
- `stats_test` becomes `test`, an `EcdfTest` member: `KS`, `WILCOXON`, `WASSERSTEIN`, or `DTS`.
- `Wilcoxen` of R becomes `EcdfTest.WILCOXON`.
- `alternative` takes an `Alternative` member, not text.
- `None` gives `GREATER`, the default of R.
- An `alternative` with `DTS` or `WASSERSTEIN` raises `InputError` (`D-035`).
- `ecdf_test` drops each row with a NaN `log2FoldChange`, for each test (`D-036`).
- A set with no row after the drop raises `InputError` that names the set.
- R stops with a message such as `not enough 'x' data`.
- `DTS` must have `rng`, and `n_boot` replaces `nboots`.
- `n_interp` replaces `n.interp`.
- `n_boot` and `n_interp` default to `None`, which gives the default of the test: 2000 and 10000.
- A value that is not `None` for an option of another test raises `InputError`.
- `DTS` and `WASSERSTEIN` raise `InputError` for an infinite `log2FoldChange`, and the message names the set and the gene ID.
- R gives the statistic `Inf` for `DTS` (`D-034`), and it stops for `WASSERSTEIN`.
- A pooled range that a float cannot hold also raises `InputError` that names `target` and `background` (`D-034`, `D-038`).
- `WASSERSTEIN` also raises `InputError` for a set with fewer than 10 values or fewer than 2 distinct values, as `wass_dist` stops.
- The function returns an `EcdfTestResult`, not the object of the test.
- `$statistic` and `$p.value` become `statistic` and `p_value`.
- The other fields are `test`, `method`, `alternative`, `n_target`, `n_background`, and `areas`.
- For `WASSERSTEIN`, `statistic` is `neg.cens.auc`, as in R.
- `p_value` and `method` are `None`, and `areas` holds the `WassersteinAreas`.
- A table without `gene_id` or `log2FoldChange` raises `InputError` that lists the missing columns in sorted order.
- A table in which one of these columns occurs more than once also raises `InputError`, but R uses the first column (`D-051`).
- A `log2FoldChange` column that is not numeric raises `InputError`.
- A column of bool or complex values raises `InputError` too, because `is.numeric` of R is false for it.

### `deseq_fc_ecdf`

- `stats_test` and `alternative` become `test` and `alternative` of `ecdf_test`.
- `$stats` becomes the `EcdfTestResult` of `ecdf_test`.
- `gene_lists` becomes `groups` of `assign_groups`, a mapping from a group name to gene IDs.
- `null_name` and `target_name` have no counterpart.
- The caller gives `target` and `background` to `ecdf_test`.
- `factor_order` has no counterpart.
- The order of the keys of `groups` is the order of the categories.
- `$plot` becomes the `Figure` of `plot_ecdf(de, groups, result)`.
- `title`, `palette`, and `x_lims` become `title`, `palette`, and `x_limits`.
- `height` and `width` become `size`.
- `output_filename` and `dpi` have no counterpart.
- Call `Figure.savefig(path, dpi=...)`.
- R writes a PDF for each file name, also for `.png` (`R/deseq_fc_ecdf.R:357`).
- `theme` has no counterpart, because R ignores it (`:349`).
- The p-value labels obey `D-063`, and the lines and the legend obey `D-064`.

### `.check_gene_list_overlap`

- `assign_groups` does the check. It warns with `SeedMatchWarning` if two groups share a gene.
- The message names the two groups and one shared gene.
- R warns `Gene lists overlap!` also for a gene that occurs two times in one list (`D-039`).

### `.assign_groups`

- `assign_groups(de, groups)` returns a `Series` named `group` with the index of `de`. It does not add a column to a copy of `res`.
- `Group` becomes `group`.
- The dtype is categorical, and the categories are the keys of `groups` in their order, as the `factor` call of `deseq_fc_ecdf` makes them.
- A row whose gene is in no group has no value, not the text `NA`.
- If a gene is in more than one group, the last group wins, as in R.
- A group must be an iterable of gene IDs.
- A `str` raises `InputError`.
- A table without `gene_id` raises `InputError` that lists the missing columns, and so does a `gene_id` column that is not text.
- A repeated `gene_id` column also raises `InputError`, but R uses the first column (`D-051`).

### `.plot_ecdf`

- `res` with the column `Group` becomes `de` and `groups`.
- `test_results` becomes `result`, an `EcdfTestResult`.
- `legend_labels` has no counterpart, because the legend counts the values of each line (`D-064`).
- The other arguments are the arguments of the row of `deseq_fc_ecdf`.
- R warns `Removed N rows containing non-finite outside the scale range`, and the port warns for each group.
- R gives `1e+05` for a statistic of 100000, and the port gives `100000`.
- A value of `x_limits` or of the fold changes outside -9e306 to 9e306, or a range wider than 9e306, raises `InputError`.
- Name `D-063`.

### `SeedMatchR.theme`, `SeedMatchR.palette`

- `PALETTE` has the colors of R.
- The figure of `plot_ecdf` has the colors of the theme, and sizes near its sizes.
- There is no theme object.

### `wass_dist`

- `wasserstein_areas(seed, non_seed, *, n_interp=10000)`.
- `non.seed` becomes `non_seed`, and `n.interp` becomes `n_interp`.
- The function returns a `WassersteinAreas`, not a list.
- `abs.auc`, `neg.cens.auc`, and `pos.cens.auc` become the fields `absolute`, `negative_censored`, and `positive_censored`.
- `plot_wasserstein(seed, non_seed, n_interp=...)` draws `dens.plot`, `cdf.plot`, and `diff.plot` as the Axes of one figure, in this order.
- The title of `diff.plot` is the title of the figure.
- `wasserstein_difference` gives the data of `diff.plot`.
- The density is the exact Gaussian sum with the bandwidth of `bw.nrd0`.
- `density()` of R uses an FFT approximation, and for the test data the difference is at most 0.0014 of the maximum.
- The ECDF lines are the ECDFs of the areas (`D-065`).
- Each area has 3 decimals, rounded with the rule of R.
- A sample must have at least 10 values and at least 2 distinct values, as in R.
- NaN and an infinite value raise `InputError`, and R stops for them.
- An `n_interp` below 2 and a range above half of the largest float raise `InputError` (`D-038`).
- In `plot_wasserstein`, a pooled value outside -9e306 to 9e306, a pooled range wider than 9e306, and a density bandwidth below about 4.43e-308 raise `InputError`.
- The function gives no warning for equal values (`D-037`).

### `stats::ks.test`

- `ks_test(x, y, *, alternative=Alternative.TWO_SIDED)` replaces the two-sample call.
- The one-sample test and the arguments `exact`, `simulate.p.value`, and `B` have no counterpart.
- The function uses the rule of R for the p-value: exact if the product of the two sizes is below 10000, and else asymptotic.
- `alternative` takes an `Alternative` member, not text.
- `two.sided`, `less`, and `greater` become `TWO_SIDED`, `LESS`, and `GREATER`.
- The function returns a `seedmatch.stats.TwoSampleResult` with the fields `statistic`, `p_value`, `method`, `alternative`, `n_x`, and `n_y`.
- `method` is `PValueMethod.EXACT` or `PValueMethod.ASYMPTOTIC`.
- The asymptotic p-values are accurate (`D-032`).
- NaN raises `InputError`, but R drops an NA value (`D-036`).
- The function gives no warning for ties (`D-037`).

### `stats::wilcox.test`

- `wilcoxon_test(x, y, *, alternative=Alternative.TWO_SIDED)` replaces the two-sample call with the default arguments of R.
- The arguments `paired`, `mu`, `exact`, `correct`, `conf.int`, and `conf.level` have no counterpart.
- The function uses the rule of R: an exact p-value if both sizes are below 50 and no two values are equal, and else the normal approximation with the continuity correction.
- `statistic` is W.
- `alternative` and the result are as for `ks_test`.
- If each value is equal, the p-value is 1 (`D-033`).
- NaN raises `InputError`, but R drops an NA value (`D-036`).
- The function gives no warning for ties (`D-037`).

### `twosamples::dts_test`

- `dts_test(x, y, *, rng, n_boot=2000)`.
- `a` and `b` become `x` and `y`, and `nboots` becomes `n_boot`.
- The caller must give the new argument `rng`. It is a `numpy.random.Generator` or a non-negative int seed, and it replaces the random state of R (`set.seed`).
- `p`, `keep.boots`, and `keep.samples` have no counterpart.
- The function returns a `TwoSampleResult`, not a named vector.
- `method` is `PValueMethod.PERMUTATION`, and `alternative` is `None`.
- The statistic is the same as in R.
- The p-value agrees with R only within the random error of the permutations.
- An infinite value, an empty sample, and a pooled range that overflows raise `InputError` (`D-034`).
- NaN raises `InputError`, as the stop of R does (`D-036`).

### `PRROC::pr.curve` (with `weights.class0` and `curve = TRUE`)

- `pr_curve(scores, labels)`.
- `scores.class0` becomes `scores`, and `weights.class0` becomes `labels`, with the values 0 and 1.
- A bool is valid, and a label of 1 or `True` is a positive.
- Soft labels, `scores.class1`, `weights.class1`, `sorted`, `minStepSize`, `max.compute`, `min.compute`, `rand.compute`, and `dg.compute` have no counterpart.
- The step of the recall keeps the default of PRROC.
- `auc.integral` and `auc.davis.goadrich` become `auc_integral` and `auc_davis_goadrich`.
- `auc_integral` loses fewer digits than PRROC if the negatives far outnumber the positives (`D-050`).
- The labels are hard, thus `auc_davis_goadrich` is always a number.
- The matrix `curve` becomes a `DataFrame` with the columns `recall`, `precision`, and `threshold`.
- The function returns a `PrCurve`, which also has `n_positive` and `n_negative`.
- NaN in `scores`, a label that is not 0, 1, `False`, or `True`, labels of another length, and labels without a positive raise `InputError`.
- An infinite score is valid.

### `PRROC::roc.curve` (with the same arguments)

- `roc_curve(scores, labels)`.
- `scores.class0` becomes `scores`, and `weights.class0` becomes `labels`, as for `pr_curve`.
- Soft labels, `scores.class1`, `weights.class1`, `sorted`, `max.compute`, `min.compute`, and `rand.compute` have no counterpart.
- `auc` keeps its name.
- The matrix `curve` becomes a `DataFrame` with the columns `fpr`, `tpr`, and `threshold`.
- The threshold column keeps the rule of PRROC: the first row has the highest score, and the last row repeats the lowest score.
- The function returns a `RocCurve`, which also has `n_positive` and `n_negative`.
- NaN in `scores` raises `InputError`, and so does a label that is not 0, 1, `False`, or `True`.
- Labels of another length, and labels without a positive or without a negative, raise it too.
- An infinite score is valid.

### `benchmark_preds`

- `benchmark_predictions(de, *, column="mer7m8", padj_cutoff=0.15, l2fc_cutoff=0.0)`.
- `res` becomes `de`, and `preds_col` becomes `column`.
- `padj` and `l2fc` become `padj_cutoff` and `l2fc_cutoff`.
- `column`, `padj_cutoff`, and `l2fc_cutoff` are keyword-only, and their defaults are the R defaults.
- The function returns a `Benchmark` with the fields `pr` (a `PrCurve`) and `roc` (a `RocCurve`), not a list with `PRC` and `ROC`.
- `padj_cutoff` has an effect, and the ROC curve uses the labels as the PR curve does (`D-006`).
- A row is a positive if its `gene_id` is in the gold standard, as `%in%` does in R.
- A row with a NaN `padj` or a NaN `log2FoldChange` does not pass the cutoffs. It is a negative unless another row of its gene passes, as in R.
- A `gene_id` that is `None`, NaN, or `pd.NA` is one gene, as `NA` is in R.
- A NaN score raises `InputError`.
- So do a cutoff that is a bool, text, or NaN, and a gold standard without a positive or without a negative row (`D-048`).
- R documents the result as figures (`R/benchmark.R:11-13`). But it returns the curve objects of PRROC (`:32-40`), and it draws nothing.
- In R, `plot(res$PRC)` and `plot(res$ROC)` draw the figures.
- `plot_benchmark(benchmark)` draws the PR curve and the ROC curve in one figure.

### `plot.PRROC` (PRROC 1.4, for the curves of `benchmark_preds`)

- `plot_benchmark(benchmark, size=...)` takes a `Benchmark`, not one curve object.
- The figure has the defaults of `plot.PRROC`, and the other arguments of `plot.PRROC` have no counterpart.
- `size` is the width and the height of the figure in inches.
- Each curve has its own Axes and its own color bar.
- The PR area is `auc_integral`, the default `auc.type` of `plot.PRROC`.
- For many negatives, it differs from PRROC (`D-050`).
- For some AUC values, R gives a short scientific form, for example `1e-04` for 0.0001.
- The port gives `0.0001`.
- A curve of one distinct score, or with an infinite threshold, gives a figure (`D-069`).
- A finite threshold outside -9e305 to 9e305 raises `InputError` (`D-069`).

### `sequence_kmer_counts`

- `seqs` becomes `sequences`, a `SequenceSet`.
- `width` becomes `k`, and `step` keeps its name.
- `k` and `step` are keyword-only, and their defaults are the R defaults 7 and 1.
- `kmer_counts` returns a `KmerCounts`. Its `matrix` is a sparse `csr_array` with one row for each sequence and `4**k` columns in the order of Biostrings.
- `tx2gene` has no counterpart, because each sequence of a `SequenceSet` has its gene ID.
- R gives `NA` to a name that `tx2gene` lacks, and the IDs `1`, `2`, and so on to a set without names.
- `out.path` and `path.basename` have no counterpart, because the core writes no file.
- Use `to_frame().to_csv(path, sep="\t")`.
- `to_frame()` gives the data frame of R with `seq_id` for `tx_id`, and with `int64` for the integers of R.
- `simplify.as = "list"` and `"collapsed"` have no counterpart.
- Use a row of `matrix`, or `matrix.sum(axis=0)`.
- R stops with `'width' must be an integer >= 1` and `'step' must be a positive integer`.
- The port raises `InputError` for a `k` outside 1 to 12 (`D-043`) and for a `step` below 1.
- R stops for a repeated name, but the port gives each sequence its own row (`D-044`).

### `kmer_hist`

- `seqs`, `width`, and `step` become the arguments of `kmer_counts`, and `kmer_gene_counts` takes its result.
- `target.seq` has no argument, because `result[target.sequence]` gives the target count.
- `plot_kmer_histogram(gene_counts, target.sequence, title=..., max_bins=..., size=...)` draws the histogram, and `kmer_histogram(gene_counts, max_bins=...)` gives its bins (`D-066`).
- `max_bins` is the most bins, and it replaces `breaks = 50`.
- `size` is the width and the height of the figure in inches.
- R gives `lw = 3` to `abline`, and R reads it as `lwd`.
- The line of the port has the width 2.25 points, the width of the R line.
- The function returns the gene count of each k-mer, as a `Series`.
- R returns `integer(0)` (`D-010`).
- The function counts the distinct genes, but R counts the sequences (`D-042`).

### `plotAlignment` with `output.type = "data.frame"`

- `match.info` becomes `match`, a row of `find_matches`.
- R reads the column `names`, and the port reads `seq_id`.
- `guide.seq` becomes a `Guide`, and `seqs` becomes `sequences`.
- `seq.start.pos` and `seq.stop.pos` become `seed`, a `Seed` or a `CustomSeed(name, start, stop)`.
- `seedmatch.seed.seed_window` gives the positions.
- The `A` of `MER8` and `MER7A1` takes g1, thus an R call for `mer8` must use `seq.start.pos = 1`.
- `left.pad` and `right.pad` become `left_pad` and `right_pad`.
- `aln.name` has no counterpart.
- The function returns the three values of the column `alignment` as one `str` with line ends between the lines.
- The guide line shows the full guide (`D-045`), and a position outside the sequence is a space (`D-046`).
- A match whose width differs from the width of the window raises `InputError`.
- A search with `indels` on can give such a match.
- R stops with `nchar(s1) == nchar(s2) is not TRUE`.
- A pad that is longer than the sequence of the match also raises `InputError` (`D-046`).
- `seedmatch.plot.plot_alignment` draws the figure of `output.type = "ggplot"`.

### `plotAlignment` with `output.type = "ggplot"`

- `plot_alignment(text)` takes the text of `alignment_text`, thus the arguments of R go to `alignment_text`.
- `aln.name` has no counterpart.
- The figure shows each line of the text after its name.
- R names the guide line `miRNA` in the figure, and `siRNA ` in the data.frame.
- The port uses `siRNA ` (`D-067`).
- The guide line shows the full guide (`D-045`), and a position outside the sequence is a space (`D-046`).
- A text with another number of lines raises `InputError`.

### `plotIRanges`, `plotRangesDF`

- `ir` becomes `ranges`, a table with the integer columns `start` and `end`, for example a result of `find_matches`.
- A data.frame stops `plotIRanges` with `invalid time series parameters specified (1)`, because it takes only an `IRanges`.
- `xlim`, `main`, and `col` become `x_limits`, `title`, and `color`.
- The default title is `Ranges`, not the default `IRanges Plot` of `plotIRanges` or `Ranges DF Plot` of `plotRangesDF`.
- `bar_height` keeps its name.
- R accepts each value of `bar_height`, but the port accepts only a value above 0 and at most 1.
- A higher bar covers the bar of the next row.
- A `start` or an `end` outside -2^52 to 2^52 raises `InputError`.
- A limit outside -9e306 to 9e306, or an `x_limits` range wider than 9e306, raises `InputError`.
- `size` is the width and the height of the figure in inches.
- `...` has no counterpart.
- `NAMESPACE` does not export the R functions.
- Name `D-068`.

### `.matchStrings`, `.make_space`, `.create_alignment_df`

- These helpers make the lines.
- `alignment_text` does their work.
- The `gsub` calls of `.matchStrings` (`R/plotAlignment.R:89-90`) change no line, because the pair line holds only the bar character and spaces.

### `match_distances`, `.get_centers`, `.get_fivePrimeEnd`, `.get_threePrimeEnd`

- `gr_list` becomes `matches`, a table of `find_matches`.
- The rows of one `seq_id` replace an element of the list.
- `feature.pos` becomes `anchor`, a `MatchAnchor`. It is keyword-only, and its default is `FIVE_PRIME_END`, the R default `fivePrimeEnd`.
- `THREE_PRIME_END` and `CENTER` are `threePrimeEnd` and `center`.
- The function returns a `DataFrame` with `seq_id`, `gene_id`, and `distance`, not a list of vectors.
- `distance` has the type `float64`, because `CENTER` gives half values.
- A list without names stops in R with `attempt to select less than one element in OneIndex`.
- The port uses `seq_id` as the name.
- A sequence with one match gives no distance, but R gives `NA` and 0 (`D-011`).

### `generate_e_value`

- `result_output` becomes `de`, and `random_seq_results` becomes `null`, a required argument.
- The column `overall_matches` becomes the argument `column`. This column is the count of `run_wobble_bulges_pipeline` (`R/generate_e_value.R:19`, `R/run_wobble_bulges_pipeline.R:51-55`).
- The null is a `DataFrame` with the gene IDs as its index and one column for each random guide.
- The R null has a `symbol` column.
- Genes match by `gene_id`, not by `symbol` (`D-008`).
- R has no function that makes the null, and its default null does not exist.
- `e_value_null(guide, sequences, *, rng, n_random=1000, seed=Seed.MER7M8, options=SearchOptions(indels=True), workers=None)` makes the null of the column of `annotate_seed_counts` from `n_random` shuffles of the guide.
- `e_value_null` counts the matches of `count_matches` for each shuffle.
- `wobble_bulge_null(guide, sequences, options=WobbleBulgeOptions(), *, rng, n_random=1000, workers=None)` makes the null of the column `overall_matches` of `annotate_wobble_bulge_counts`.
- `wobble_bulge_null` counts the rows of `wobble_bulge_matches` for each shuffle.
- Shuffle j is the same guide in both functions for the same `rng`.
- A null of another count gives wrong e-values and no message.
- `rng` is a `numpy.random.Generator` or a non-negative int.
- R compares the logarithms of the counts. This gives the same result as the counts.
- A count that is NaN, infinite, or negative raises `InputError`, and each row of a gene gets its own e-value (`D-047`).

### `get_network_matches` (the TRRUST URLs)

- If `network_data` is `NULL`, R reads `trrust_rawdata.human.tsv` or `trrust_rawdata.mouse.tsv` from `https://www.grnpedia.org/trrust/data/` (`R/get_network_matches.R:21-22`).
- `seedmatch-fetch trrust ORGANISM DEST` downloads the file of one organism.
- `ORGANISM` is `human` or `mouse`, as in the `organism` argument of R.
- The file keeps the name of the URL.
- `read_trrust` reads the file.

### `get_network_matches`

- `read_trrust(path)` reads the TRRUST file and gives the columns `source`, `target`, `edge_type`, and `evidence`.
- `read.table` of R reads line 1 as the column names, but the file has no header (`D-057`).
- `trrust_scores(network, seed_genes, *, aggregator=Aggregator.MEAN)` gives the scores.
- `ENSG_Genelist` becomes `seed_genes`, and `aggregator` becomes an `Aggregator`.
- `network_data` becomes `network`, an argument without a default, because the core does not download.
- `organism` has no counterpart.
- R stops with `Organism not supported` for an `organism` that is not `human` or `mouse`.
- The messages come from the regulators, not from the targets (`D-007`).
- A seed gene list with no gene of the network gives a `SeedMatchWarning` (`D-058`).

### `message_passing`, `as_dataframe_mp_matches`

- `ENSG` becomes `gene`.
- `repression_matches` and `activation_matches` become `repression_score` and `activation_score`.
- The new columns are `seed`, `repressors`, and `activators`.
- An `agg` that is not `mean` or `degree_importance` makes R return the seed flags without a message.
- `aggregator` raises `InputError` if it is not an `Aggregator`.
- The score of a gene comes from its regulators (`D-007`).

### `get_stringPPI_to_stringGGI` (the STRING download)

- `STRINGdb::STRINGdb$new(species = sp_number)` sets the STRING version (`R/get_stringPPI_to_stringGGI.R:15`), and `string_db$get_graph()` downloads the network (`:16`).
- The call gets no `version`, thus STRINGdb 2.22.0 uses the current STRING version at the time of the run (`R/rstring.R:145-151` of STRINGdb 2.22.0).
- `seedmatch-fetch string TAXON DEST` downloads `<TAXON>.protein.links.v12.5.txt.gz` of the pinned STRING version 12.5 (`D-041`).
- `TAXON` is the NCBI taxon, the `sp_number` of R.
- The registry pins the size and the SHA-256 digest.
- `read_string_links` reads the file, and `string_gene_network` maps the protein IDs.

### `get_stringPPI_to_stringGGI`

- `read_string_links(path, *, score_threshold=400, percentile=0.95)` reads the links file and keeps the links above the percentile.
- `p` becomes `percentile`, and `score_threshold` is new, with the default 400 of STRINGdb.
- `sp_number` is the `TAXON` of `seedmatch-fetch string`, which selects the file.
- `read_string_links` returns a `StringLinks` with the table `links` and the `cutoff`, the value of `quantile` in R.
- `string_gene_network(links, reference)` takes `StringLinks.links`, maps the protein IDs, and `sp` becomes the `Reference` of the same species and build.
- The result is a `GeneNetwork`, not a sparse matrix.
- `edges` has a row for each order of each pair, with the score of the pair.
- The R matrix holds 2 for each pair, the count of its parallel edges.
- The divergences are `D-013`, `D-059`, `D-060`, and `D-061`.

### `create_graph`

- `string_scores(network, de, *, column, aggregator=Aggregator.MEAN)`.
- `string_ppi` becomes `network`, the `edges` of a `GeneNetwork`.
- `res` becomes `de`, and the column `overall_matches` becomes the argument `column`.
- R gives an igraph object with the seed count of each gene as the vertex `feature`.
- The port gives a table with the columns `gene_id`, `seed_count`, `neighbors`, and `score`.
- `seed_count` is the `feature`, and `neighbors` and `score` are new.
- A gene with more than one row in `de` gets the highest count.
- A count that is NaN, infinite, or negative raises `InputError`.
- If no gene of `de` is a gene of the network, the function warns with `SeedMatchWarning` (`D-062`).

### `run_wobble_bulges_pipeline`

- `guide.seq` becomes a `Guide`, `features` becomes `sequences`, and `res` becomes `de` of `annotate_wobble_bulge_counts`.
- R ignores `res`.
- `bulges_lim` and `nt` become `max_bulge` and `bulge_free` of `WobbleBulgeOptions`, with the R defaults 4 and 3.
- `silent` has no counterpart, because the port prints nothing.
- R gives one column `SeedMatches_<j>` for each wobble seed, and `overall_matches`.
- The port has one pattern for each set of wobbles, thus it gives only `overall_matches`.
- `transcript` becomes `seq_id` and `gene_id`.
- `wobble_bulge_matches` gives one row for each distinct start of a site, and `annotate_wobble_bulge_counts` gives the count of each gene (`D-005`).

### `.get_all_combination_w_b`, `.find_factorial_internal`, `.get_target_bulged_guides`, `.sequence_with_tb`, `check_bulges_nt`

- One IUPAC pattern replaces the wobble seeds.
- The bulged patterns insert bases, and `bulge_free` replaces `check_bulges_nt` (`D-005`).
- `wobble_bulge_patterns(guide, options)` returns a table with the columns `pattern`, `bulge_after`, and `bulge_size`.

### `get_compensatory_sequences`

- The signature is `compensatory_sites(guide, sequences, options=CompensatoryOptions(), *, seed=Seed.MER7M8, workers=None)`, and `options` is in position, as in `wobble_bulge_matches`.
- `guide.seq` becomes a `Guide`, and `sequences` keeps its name.
- `seed` becomes a `SeedSpec`, and the default is `Seed.MER7M8`, the seed of the R vignette.
- R takes any target text.
- The other arguments become the fields of `CompensatoryOptions`.
- `compensatory_region_length_min` and `compensatory_region_length_max` become `min_length` and `max_length`, and `NA` becomes `None`.
- `mismatch_` and `mismatch_compensatory_region_` become `max_seed_mismatch` and `max_region_mismatch`.
- `limit = L` becomes `max_distance = L - 2`.
- `first_position` is new (`D-056`).
- The R columns map to the port columns.
- `TX_ID` is `seq_id`, `mismatch_in_seed` is `seed_mismatches`, and `CompensatoryLength` is `region_length`.
- `Compensatory_Region_in_Transcript` is `region_target`, `Distance` is `distance`, and `GC_Content` is `gc_content`.
- A new column, `gene_id`, follows `seq_id`.
- R gives the transcript coordinates of the seed and of the region in swapped order.
- `Seed_Start` is the end of the seed match, and `Seed_End` is its start. Thus `Seed_End` and `Seed_Start` become `seed_start` and `seed_end`.
- `Compensatory_End` and `Compensatory_Start` become `region_start` and `region_end`.
- `relative_position_guide_seq_start` and `relative_position_guide_seq_end` become `guide_start` and `guide_end`.
- `Seed_in_Transcript` and `Seed_Region_in_Guide` have no column, because each row holds `guide.target(seed).sequence` and its reverse complement.
- `Compensatory_Region_in_Guide` has no column. It is the DNA text of `guide.sequence[guide_start - 1 : guide_end]`.
- The gene list of the R vignette, `unique(comp_matches$TX_ID)`, becomes `sites["gene_id"].unique()`.
- `workers` is new. It sets the width of the thread pool.
- The function runs the search (`D-004`), and it searches each sequence of the set (`D-052`). It drops a region only for a longer region of the same seed match (`D-053`). It gives the guide positions as integers and the number of seed mismatches (`D-054`). It raises `InputError` for lengths and distances that select no region (`D-055`).

### `relative_position_seq`, `.gc_content`

- `compensatory_sites` gives the window of the guide for each row, and the percentage of `G` and `C` of the region (`D-054`).

### `.update_sequence_at`

- No SeedMatchR function calls it.
