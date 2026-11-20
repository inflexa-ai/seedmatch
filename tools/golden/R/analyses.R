source("tools/golden/R/common.R")

PRROC_NOTE <- paste0("PRROC: ", packageVersion("PRROC"))
PATCH_DIR <- "tools/golden/patches"

write_analysis <- function(df, path, notes = character(0)) {
  stop_on_nan(df)
  write_golden(df, file.path("analyses", path), "analyses.R", notes)
}

data <- load_test_data()
seqs <- data$seqs
res <- data$res
seq_widths <- setNames(Biostrings::width(seqs), names(seqs))
target_of <- function(seed, ...) SeedMatchR::get_seed(GUIDE, seed, ...)$Target.Seq

# ---- k-mers ----
cases_fixture <- utils::read.delim(
  file.path(DATA_DIR, "kmer_cases.tsv"), colClasses = "character", na.strings = character(0)
)
stopifnot(identical(cases_fixture$sequence[4], ""))
case_seqs <- Biostrings::DNAStringSet(setNames(cases_fixture$sequence, cases_fixture$seq_id))
case_tx2gene <- data.frame(tx_id = cases_fixture$seq_id, gene_id = cases_fixture$gene_id)

test_names <- names(seqs)
test_tx2gene <- data.frame(tx_id = test_names, gene_id = test_names)

kmer_cases <- list(
  test_data_first10_k7_s1 = list(seqs = seqs[1:10], tx2gene = test_tx2gene, width = 7, step = 1),
  test_data_k3_s2 = list(seqs = seqs, tx2gene = test_tx2gene, width = 3, step = 2),
  cases_k1_s1 = list(seqs = case_seqs, tx2gene = case_tx2gene, width = 1, step = 1),
  cases_k2_s1 = list(seqs = case_seqs, tx2gene = case_tx2gene, width = 2, step = 1),
  cases_k3_s2 = list(seqs = case_seqs, tx2gene = case_tx2gene, width = 3, step = 2),
  cases_k4_s3 = list(seqs = case_seqs, tx2gene = case_tx2gene, width = 4, step = 3)
)

kmer_frames <- lapply(kmer_cases, function(case) {
  SeedMatchR::sequence_kmer_counts(case$seqs, case$tx2gene, width = case$width, step = case$step)
})

# One row for each count above 0, in the order of the rows of the frame and then of its columns.
nonzero_counts <- function(id, frame) {
  counts <- as.matrix(frame[, -(1:2)])
  cells <- which(counts > 0, arr.ind = TRUE)
  cells <- cells[order(cells[, "row"], cells[, "col"]), , drop = FALSE]
  data.frame(
    case = id, seq_id = frame$tx_id[cells[, "row"]], gene_id = frame$gene_id[cells[, "row"]],
    kmer = colnames(counts)[cells[, "col"]], count = counts[cells], stringsAsFactors = FALSE
  )
}
kmer_counts <- do.call(rbind, Map(nonzero_counts, names(kmer_frames), kmer_frames))

duplicate_names_error <- error_of(suppressWarnings(SeedMatchR::sequence_kmer_counts(
  Biostrings::DNAStringSet(c(a = "ACGT", a = "ACGA")), data.frame(tx_id = "a", gene_id = "g"), width = 2
)))
wide_dim <- dim(Biostrings::oligonucleotideFrequency(Biostrings::DNAStringSet("ACGTACGTACGTAC"), width = 13))
stopifnot(wide_dim == c(1, 4^13))
write_analysis(
  kmer_counts, "kmer_counts.tsv",
  c(
    "the rows are the counts above 0 of sequence_kmer_counts, in the order of its rows and then of its columns",
    "the cases test_data_* read tests/data/Test_seqs.fa, and tx2gene maps each name to itself",
    "the cases cases_* read tests/data/kmer_cases.tsv, and tx2gene maps seq_id to gene_id",
    paste0("sequence_kmer_counts for a set with a repeated name stops: ", duplicate_names_error),
    paste0(
      "dim of oligonucleotideFrequency(DNAStringSet(\"ACGTACGTACGTAC\"), width = 13) is ",
      wide_dim[1], " x ", wide_dim[2]
    )
  )
)

write_analysis(
  kmer_frames$cases_k2_s1, "kmer_frame.tsv",
  "the data frame of sequence_kmer_counts for cases_k2_s1, as R gives it: tx_id, gene_id, and one column for each 2-mer"
)

# The presence matrix of R/kmer_hist.R:19-28.
gene_counts <- function(id, set, width) {
  presence <- Biostrings::oligonucleotideFrequency(x = set, width = width, step = 1, simplify.as = "matrix")
  presence[presence > 0] <- 1
  n <- colSums(presence)
  data.frame(case = id, kmer = names(n), n = unname(n), stringsAsFactors = FALSE)
}
kmer_gene_counts <- rbind(
  gene_counts("test_data_k7_s1", seqs, 7),
  gene_counts("cases_k2_s1", case_seqs, 2)
)

mer7m8_target <- as.character(target_of("mer7m8"))
stopifnot(mer7m8_target == "CTCTATA")
mer7m8_row <- kmer_gene_counts$case == "test_data_k7_s1" & kmer_gene_counts$kmer == mer7m8_target
stopifnot(sum(mer7m8_row) == 1)

grDevices::pdf(NULL)
hist_value <- SeedMatchR::kmer_hist(case_seqs, "ACG", width = 3)
invisible(grDevices::dev.off())
stopifnot(identical(hist_value, integer(0)))
write_analysis(
  kmer_gene_counts, "kmer_gene_counts.tsv",
  c(
    "n is the number of sequences with the k-mer: colSums of the presence matrix of R/kmer_hist.R:19-28, with step 1",
    paste0(
      "n of the mer7m8 target ", mer7m8_target, " for test_data_k7_s1 is ", kmer_gene_counts$n[mer7m8_row]
    ),
    "the value of kmer_hist(cases_*, \"ACG\", width = 3), with pdf(NULL) as the device, is integer(0)"
  )
)

# ---- Alignment ----
plot_alignment_patch <- file.path(PATCH_DIR, "plotAlignment.R")
guide_length <- nchar(GUIDE)

# The window of plotAlignment, in positions of the sequence: the match, the part of the guide outside
# the seed, and the pads.
window_of <- function(matches, start_pos, stop_pos, left_pad = 0, right_pad = 0) {
  list(
    first = matches$start - (guide_length - stop_pos) - left_pad,
    last = matches$end + (start_pos - 1) + right_pad
  )
}
window_fits <- function(matches, ...) {
  window <- window_of(matches, ...)
  window$first >= 1 & window$last <= unname(seq_widths[matches$names])
}

alignment_rows <- list()
add_alignments <- function(prefix, matches, seed, start_pos, stop_pos, left_pad = 0, right_pad = 0) {
  stopifnot(nrow(matches) > 0, all(window_fits(matches, start_pos, stop_pos, left_pad, right_pad)))
  for (i in seq_len(nrow(matches))) {
    aligned <- SeedMatchR::plotAlignment(
      matches[i, ], GUIDE, seqs,
      left.pad = left_pad, right.pad = right_pad, seq.start.pos = start_pos, seq.stop.pos = stop_pos
    )
    stopifnot(nrow(aligned) == 3, length(unique(nchar(aligned$alignment))) == 1)
    alignment_rows[[length(alignment_rows) + 1]] <<- data.frame(
      case = paste0(prefix, "_", i), seq_id = matches$names[i], start = matches$start[i],
      end = matches$end[i], seed = seed, start_pos = start_pos, stop_pos = stop_pos,
      left_pad = left_pad, right_pad = right_pad, line = 1:3, text = aligned$alignment,
      stringsAsFactors = FALSE
    )
  }
}

mer7m8_matches <- SeedMatchR::vmatchPattern2(target_of("mer7m8"), seqs)
mer8_matches <- SeedMatchR::vmatchPattern2(target_of("mer8"), seqs)
mer6_matches <- SeedMatchR::vmatchPattern2(target_of("mer6"), seqs)
custom_matches <- SeedMatchR::vmatchPattern2(target_of("custom", start.pos = 3, stop.pos = 12), seqs)
# Biostrings::vmatchPattern counts substitutions only. A window must end inside its sequence.
full_mismatch_matches <- function(max_mismatch) {
  found <- as.data.frame(unlist(Biostrings::vmatchPattern(target_of("Full"), seqs, max.mismatch = max_mismatch)))
  found[window_fits(found, 2, guide_length), ]
}
full_max_mismatch <- 2
full_matches <- full_mismatch_matches(full_max_mismatch)
while (nrow(full_matches) == 0) {
  full_max_mismatch <- full_max_mismatch + 1
  full_matches <- full_mismatch_matches(full_max_mismatch)
}
rownames(full_matches) <- NULL

first_name <- names(seqs)[1]
first_width <- seq_widths[[first_name]]
match_in_first <- function(start, end) data.frame(start = start, end = end, width = end - start + 1, names = first_name)

with_seedmatchr_patch(plot_alignment_patch, {
  add_alignments("mer7m8", mer7m8_matches, "mer7m8", 2, 8)
  add_alignments("mer7m8_pad", mer7m8_matches[1:3, ], "mer7m8", 2, 8, left_pad = 3, right_pad = 2)
  add_alignments("mer8", mer8_matches, "mer8", 1, 8)
  add_alignments("mer6", mer6_matches[1:5, ], "mer6", 2, 7)
  add_alignments("full_mm2", full_matches, "Full", 2, guide_length)
  add_alignments("custom", custom_matches, "custom", 3, 12)
  # These windows lie at the edges of the first sequence. The positions are not results of a search.
  add_alignments("edge_first", match_in_first(guide_length - 7, guide_length - 1), "mer7m8", 2, 8)
  add_alignments("edge_last", match_in_first(first_width - 7, first_width - 1), "mer7m8", 2, 8)
})
alignment <- do.call(rbind, alignment_rows)

indel_matches <- SeedMatchR::vmatchPattern2(target_of("mer7m8"), seqs, max.mismatch = 1, with.indels = TRUE)
indel_row <- indel_matches[which(indel_matches$width != 7)[1], ]
data_frame_row <- SeedMatchR::SeedMatchR(
  seqs = seqs, sequence = GUIDE, res.format = "data.frame", seed.name = "mer7m8"
)[1, ]
patched_error <- function(...) {
  with_seedmatchr_patch(plot_alignment_patch, error_of(SeedMatchR::plotAlignment(..., GUIDE, seqs)))
}

alignment_notes <- c(
  paste0("the mer7m8 target is ", mer7m8_target, ", and the guide is ", GUIDE),
  "text is the column alignment of plotAlignment of tools/golden/patches/plotAlignment.R, with each space",
  paste0(
    "the unpatched plotAlignment stops for mer7m8_1: ",
    error_of(SeedMatchR::plotAlignment(mer7m8_matches[1, ], GUIDE, seqs))
  ),
  paste0(
    "the unpatched plotAlignment stops for a row of SeedMatchR(res.format = \"data.frame\"): ",
    error_of(SeedMatchR::plotAlignment(data_frame_row, GUIDE, seqs))
  ),
  paste0(
    "the patched plotAlignment stops for a window that starts before 1 (match 5 to 11, window from -10): ",
    patched_error(match_in_first(5, 11))
  ),
  paste0(
    "the patched plotAlignment stops for a window that starts at 0 (match 15 to 21): ",
    patched_error(match_in_first(15, 21))
  ),
  paste0(
    "the patched plotAlignment stops for a window that ends after its sequence (the sequence has ",
    first_width, " bases, and the match ends at the last base): ",
    patched_error(match_in_first(first_width - 6, first_width))
  ),
  paste0(
    "the patched plotAlignment stops for a match with an indel (width ", indel_row$width, "): ",
    patched_error(indel_row)
  ),
  paste0(
    "edge_first_1 has the window 1 to ", guide_length, " (match ", guide_length - 7, " to ", guide_length - 1,
    " of the first sequence). edge_last_1 has a window that ends at the last base (match ",
    first_width - 7, " to ", first_width - 1, " of ", first_width, " bases). These matches are not search results"
  )
)
if (full_max_mismatch != 2) {
  alignment_notes <- c(
    alignment_notes,
    paste0("the Full search has no match with max.mismatch 2, thus full_mm2_* uses max.mismatch ", full_max_mismatch)
  )
}
write_analysis(alignment, "alignment.tsv", alignment_notes)

# ---- Distances ----
ranges_by_name <- function(names, start, width) {
  split(IRanges::IRanges(start = start, width = width), factor(names, levels = unique(names)))
}

# The vmatchPattern2 of SeedMatchR gives the same indel ranges as the patched function.
patched_vmatch <- new.env()
sys.source(file.path(PATCH_DIR, "vmatchPattern2.R"), envir = patched_vmatch)
stopifnot(identical(
  indel_matches,
  patched_vmatch$vmatchPattern2(target_of("mer7m8"), seqs, max.mismatch = 1, with.indels = TRUE)
))

distance_cases <- list(
  mer6 = ranges_by_name(mer6_matches$names, mer6_matches$start, mer6_matches$width),
  mer7m8_indel = ranges_by_name(indel_matches$names, indel_matches$start, indel_matches$width),
  synthetic = ranges_by_name(
    c("a", "a", "a", "b", "c", "c"), c(10, 30, 35, 5, 100, 50), c(7, 7, 8, 7, 7, 7)
  )
)
stopifnot(length(unique(indel_matches$width)) > 1)

distance_ranges <- do.call(rbind, lapply(names(distance_cases), function(id) {
  do.call(rbind, lapply(names(distance_cases[[id]]), function(name) {
    ranges <- distance_cases[[id]][[name]]
    data.frame(
      case = id, seq_id = name, start = IRanges::start(ranges), end = IRanges::end(ranges),
      width = IRanges::width(ranges), stringsAsFactors = FALSE
    )
  }))
}))

ANCHORS <- c("fivePrimeEnd", "threePrimeEnd", "center")
distances <- do.call(rbind, lapply(names(distance_cases), function(id) {
  do.call(rbind, lapply(ANCHORS, function(anchor) {
    by_name <- SeedMatchR::match_distances(distance_cases[[id]], feature.pos = anchor)
    stopifnot(identical(names(by_name), names(distance_cases[[id]])))
    data.frame(
      case = id, anchor = anchor, seq_id = rep(names(by_name), lengths(by_name)),
      distance = as.double(unlist(by_name, use.names = FALSE)), stringsAsFactors = FALSE
    )
  }))
}))

empty_element <- SeedMatchR::match_distances(
  split(IRanges::IRanges(start = 5, width = 7), factor("a", levels = c("a", "b")))["b"],
  feature.pos = "center"
)$b
stopifnot(identical(empty_element, c(NA_real_, NA_real_)))
single_element <- SeedMatchR::match_distances(
  split(IRanges::IRanges(start = 5, width = 7), factor("a", levels = c("a", "b")))["a"],
  feature.pos = "center"
)$a
stopifnot(identical(single_element, c(NA_real_, 0)))

mer6_names <- table(mer6_matches$names)
write_analysis(
  distance_ranges, "distances_ranges.tsv",
  c(
    "the ranges that the script gives to match_distances, in list order",
    "the script splits the ranges by name with factor(names, levels = unique(names))",
    "mer6: the matches of vmatchPattern2 for the mer6 target. Some names have one match and some names have more",
    "mer7m8_indel: vmatchPattern2 for the mer7m8 target with max.mismatch = 1 and with.indels = TRUE. The widths differ",
    "synthetic: a (starts 10, 30, 35 and widths 7, 7, 8), b (start 5 and width 7), and c (starts 100, 50 and width 7)"
  )
)
write_analysis(
  distances, "distances.tsv",
  c(
    "distance holds the values of the R list of match_distances, in list order and in the order of each vector",
    "an element with one range gives NA and 0",
    paste0(
      "match_distances of a list element without a range gives ",
      paste(empty_element, collapse = " ")
    ),
    "the cases mer6 and mer7m8_indel have the same vmatchPattern2 ranges as tools/golden/patches/vmatchPattern2.R"
  )
)

# ---- E-value ----
set.seed(20261122)
genes <- sprintf("g%02d", 1:40)
observed <- sample(0:4, 40, replace = TRUE)
observed[1] <- 0L
null_counts <- matrix(sample(0:5, 1000, replace = TRUE), nrow = 40, ncol = 25)
# The rows of the null have the reverse order of the genes, thus a match by row position fails.
null_table <- data.frame(symbol = rev(genes), null_counts)
observed_table <- data.frame(symbol = genes, overall_matches = observed)
evalue_result <- SeedMatchR::generate_e_value(observed_table, null_table)
stopifnot(identical(evalue_result$symbol, genes), evalue_result$e_value[1] == 1)

evalue_observed <- data.frame(gene_id = genes, overall_matches = observed)
evalue_null <- data.frame(
  gene_id = rep(null_table$symbol, times = 25),
  shuffle = rep(1:25, each = 40),
  count = as.integer(null_counts)
)
evalue_null <- evalue_null[order(match(evalue_null$gene_id, null_table$symbol), evalue_null$shuffle), ]
stopifnot(nrow(evalue_null) == 1000)

small_null <- data.frame(symbol = c("a", "b"), r1 = c(1, 0), r2 = c(2, 1), r3 = c(3, 0), r4 = c(0, 0))
evalue_of <- function(symbol, overall_matches, null = small_null) {
  SeedMatchR::generate_e_value(data.frame(symbol = symbol, overall_matches = overall_matches), null)
}
negative_warning <- NULL
negative_value <- withCallingHandlers(
  evalue_of("a", -1)$e_value,
  warning = function(w) {
    negative_warning <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  }
)
duplicate_observed <- evalue_of(c("a", "a"), c(1, 3))
doubled_null <- rbind(small_null, data.frame(symbol = "a", r1 = 5, r2 = 5, r3 = 5, r4 = 5))
stopifnot(
  evalue_of("a", NA_real_)$e_value == 0,
  negative_value == 0, identical(negative_warning, "NaNs produced"),
  duplicate_observed$e_value == c(0.5, 0.5),
  evalue_of("a", 1)$e_value == 0.75, evalue_of("a", 3)$e_value == 0.25,
  evalue_of("a", 2, doubled_null)$e_value == 1.5
)
missing_error <- error_of(evalue_of("z", 1))
partial_observed <- evalue_of(c("a", "z"), c(1, 1))
stopifnot(identical(partial_observed$symbol, "a"), partial_observed$e_value == 0.75)
default_null_error <- error_of(SeedMatchR::generate_e_value(data.frame(symbol = "a", overall_matches = 1)))

write_analysis(evalue_observed, "evalue_observed.tsv", "the observed counts: sample(0:4, 40, replace = TRUE), then g01 is 0")
write_analysis(
  evalue_null, "evalue_null.tsv",
  c(
    "the null: a 40 x 25 matrix of sample(0:5, 1000, replace = TRUE), filled by column",
    "the rows have the reverse order of the genes g01 to g40, thus a match by row position fails",
    "the rows of this file follow the rows of the null, and then shuffle 1 to 25"
  )
)
write_analysis(
  data.frame(
    gene_id = evalue_result$symbol, overall_matches = evalue_result$overall_matches,
    e_value = evalue_result$e_value
  ),
  "evalue.tsv",
  c(
    "e_value is the result of generate_e_value for the observed counts and the null, in the row order of the result",
    paste0("generate_e_value without a null stops: ", default_null_error),
    "the next notes use the null with the genes a and b: a has the counts 1, 2, 3, 0 and b has the counts 0, 1, 0, 0",
    "the observed count NA gives the e-value 0",
    paste0("the observed count -1 gives the e-value 0, and the warning ", negative_warning),
    paste0(
      "the observed rows a = 1 and a = 3 both get ", duplicate_observed$e_value[1],
      ", because the comparison recycles the two counts over the null row. One by one they get ",
      evalue_of("a", 1)$e_value, " and ", evalue_of("a", 3)$e_value
    ),
    paste0(
      "a null with two rows for a (the second row 5, 5, 5, 5) gives the e-value ",
      evalue_of("a", 2, doubled_null)$e_value, " for a = 2"
    ),
    paste0("an observed table without a gene of the null stops: ", missing_error),
    paste0(
      "an observed table with the genes a and z, where z is not in the null, gives ", nrow(partial_observed),
      " row, for ", partial_observed$symbol, ". R drops the row of z"
    )
  )
)

# ---- Curves ----
set.seed(20261123)
curve_cases <- list()
add_curve_case <- function(id, score, label) {
  stopifnot(length(score) == length(label), any(label == 1), any(label == 0))
  curve_cases[[id]] <<- list(score = score, label = label)
}
score <- sample(0:3, 40, replace = TRUE)
add_curve_case("ties_small", score, rbinom(40, 1, 0.4))
score <- rnorm(60)
add_curve_case("continuous", score, rbinom(60, 1, 0.5))
score <- sample(0:5, 400, replace = TRUE)
label <- rbinom(400, 1, 0.3 + 0.1 * score)
stopifnot(sum(label) >= 100)
add_curve_case("large", score, label)
add_curve_case("one_value", rep(1, 10), rep(c(0, 1), 5))
add_curve_case("perfect", 5:1, c(1, 1, 0, 0, 0))
add_curve_case("first_negative", c(9, 9, 5, 5, 1, 1), c(0, 0, 1, 0, 1, 1))
add_curve_case("inverse", 1:6, c(1, 1, 1, 0, 0, 0))

pr_of <- function(case) PRROC::pr.curve(scores.class0 = case$score, weights.class0 = case$label, curve = TRUE)
roc_of <- function(case) PRROC::roc.curve(scores.class0 = case$score, weights.class0 = case$label, curve = TRUE)

curve_samples <- list()
pr_summary <- list()
pr_points <- list()
roc_summary <- list()
roc_points <- list()
collect_curves <- function(id, case, summary = TRUE) {
  pr <- pr_of(case)
  roc <- roc_of(case)
  curve_samples[[id]] <<- data.frame(case = id, score = case$score, label = case$label, stringsAsFactors = FALSE)
  pr_points[[id]] <<- data.frame(
    case = id, row = seq_len(nrow(pr$curve)), recall = pr$curve[, 1], precision = pr$curve[, 2],
    threshold = pr$curve[, 3], stringsAsFactors = FALSE
  )
  roc_points[[id]] <<- data.frame(
    case = id, row = seq_len(nrow(roc$curve)), fpr = roc$curve[, 1], tpr = roc$curve[, 2],
    threshold = roc$curve[, 3], stringsAsFactors = FALSE
  )
  if (summary) {
    pr_summary[[id]] <<- data.frame(
      case = id, auc_integral = pr$auc.integral, auc_davis_goadrich = pr$auc.davis.goadrich,
      stringsAsFactors = FALSE
    )
    roc_summary[[id]] <<- data.frame(case = id, auc = roc$auc, stringsAsFactors = FALSE)
  }
  invisible(list(pr = pr, roc = roc))
}
for (id in names(curve_cases)) collect_curves(id, curve_cases[[id]])

# ---- Benchmark ----
benchmark_patch <- file.path(PATCH_DIR, "benchmark_preds.R")
# The labels of the filter of tools/golden/patches/benchmark_preds.R.
benchmark_labels <- function(res, padj, l2fc) {
  de_genes <- dplyr::pull(dplyr::filter(res, padj <= !!padj, log2FoldChange < l2fc), gene_id)
  as.integer(res$gene_id %in% de_genes)
}

# benchmark_na has the filter of benchmark_2 and a padj of NA for two genes that pass this filter.
na_genes <- head(dplyr::pull(dplyr::filter(res, padj <= 1e-10, log2FoldChange < 0), gene_id), 2)
stopifnot(length(na_genes) == 2)
na_res <- res
na_res$padj[na_res$gene_id %in% na_genes] <- NA

benchmark_cases <- list(
  benchmark_1 = list(padj = 0.15, l2fc = 0, res = res),
  benchmark_2 = list(padj = 1e-10, l2fc = 0, res = res),
  benchmark_3 = list(padj = 1e-10, l2fc = -1, res = res),
  benchmark_na = list(padj = 1e-10, l2fc = 0, res = na_res)
)
benchmark_rows <- list()
for (id in names(benchmark_cases)) {
  args <- benchmark_cases[[id]]
  patched <- with_seedmatchr_patch(
    benchmark_patch,
    SeedMatchR::benchmark_preds(args$res, preds_col = "mer7m8", padj = args$padj, l2fc = args$l2fc)
  )
  case <- list(score = args$res$mer7m8, label = benchmark_labels(args$res, args$padj, args$l2fc))
  curves <- collect_curves(id, case, summary = FALSE)
  # The labels of the script are the labels of the patched function.
  stopifnot(
    identical(curves$pr$curve, patched$PRC$curve), identical(curves$roc$curve, patched$ROC$curve),
    curves$pr$auc.integral == patched$PRC$auc.integral
  )
  benchmark_rows[[id]] <- data.frame(
    case = id, column = "mer7m8", padj = args$padj, l2fc = args$l2fc,
    n_positive = sum(case$label), n_negative = sum(1 - case$label),
    auc_integral = patched$PRC$auc.integral, auc_davis_goadrich = patched$PRC$auc.davis.goadrich,
    roc_auc = patched$ROC$auc, stringsAsFactors = FALSE
  )
}
benchmark <- do.call(rbind, benchmark_rows)
stopifnot(benchmark$n_positive[1:2] == c(63, 24))
n_positive_of <- setNames(benchmark$n_positive, benchmark$case)
# R makes a row with a padj of NA a negative.
stopifnot(
  all(benchmark_labels(res, 1e-10, 0)[res$gene_id %in% na_genes] == 1),
  n_positive_of[["benchmark_na"]] == n_positive_of[["benchmark_2"]] - 2
)

unpatched <- SeedMatchR::benchmark_preds(res, padj = 1e-10)
unpatched_labels <- as.integer(res$gene_id %in% dplyr::pull(dplyr::filter(res, padj <= padj, log2FoldChange < 0), gene_id))
stopifnot(
  identical(
    PRROC::pr.curve(scores.class0 = res$mer7m8, weights.class0 = unpatched_labels, curve = TRUE)$curve,
    unpatched$PRC$curve
  )
)

nan_res <- res
nan_res$mer7m8 <- as.numeric(nan_res$mer7m8)
nan_res$mer7m8[1] <- NaN
nan_error <- with_seedmatchr_patch(
  benchmark_patch, error_of(SeedMatchR::benchmark_preds(nan_res, preds_col = "mer7m8", padj = 1e-10))
)
no_positive <- with_seedmatchr_patch(
  benchmark_patch, SeedMatchR::benchmark_preds(res, preds_col = "mer7m8", padj = 0)
)
stopifnot(sum(benchmark_labels(res, 0, 0)) == 0)
no_negative <- with_seedmatchr_patch(
  benchmark_patch, SeedMatchR::benchmark_preds(res, preds_col = "mer7m8", padj = 1, l2fc = Inf)
)
stopifnot(all(benchmark_labels(res, 1, Inf) == 1))

pr_file_notes <- c(
  PRROC_NOTE,
  "the call is PRROC::pr.curve(scores.class0 = score, weights.class0 = label, curve = TRUE)",
  "the cases ties_small, continuous, and large use set.seed(20261123), and the other cases are fixed values"
)
write_analysis(
  do.call(rbind, curve_samples), "curve_samples.tsv",
  c(
    "the scores and the labels (1 is a positive) of the cases of pr.tsv and roc.tsv, and of the cases of benchmark.tsv",
    "the benchmark cases use the column mer7m8 of tests/data/DESEQ2_results.tsv as the score"
  )
)
write_analysis(do.call(rbind, pr_summary), "pr.tsv", pr_file_notes)
write_analysis(
  do.call(rbind, pr_points), "pr_curve.tsv",
  c(pr_file_notes, "the rows of the curve of the cases of pr.tsv and of the benchmark cases, in the order of the matrix of PRROC")
)
write_analysis(
  do.call(rbind, roc_summary), "roc.tsv",
  c(PRROC_NOTE, "the call is PRROC::roc.curve(scores.class0 = score, weights.class0 = label, curve = TRUE)")
)
write_analysis(
  do.call(rbind, roc_points), "roc_curve.tsv",
  c(
    PRROC_NOTE,
    "the rows of the curve of the cases of roc.tsv and of the benchmark cases, in the order of the matrix of PRROC"
  )
)
write_analysis(
  benchmark, "benchmark.tsv",
  c(
    PRROC_NOTE,
    "the result of benchmark_preds of tools/golden/patches/benchmark_preds.R for load_test_data()$res and preds_col = mer7m8",
    "n_positive and n_negative count the labels of the filter of the patch. The curves are in pr_curve.tsv and roc_curve.tsv",
    paste0(
      "the unpatched benchmark_preds(res, padj = 1e-10) keeps ", sum(unpatched_labels),
      " genes, and the patched function keeps ", benchmark$n_positive[2]
    ),
    paste0(
      "the unpatched benchmark_preds(res, padj = 1e-10) gives the PR auc.integral ",
      double_to_string(unpatched$PRC$auc.integral), " and the ROC auc ", double_to_string(unpatched$ROC$auc)
    ),
    paste0("the patched function stops for a NaN score in one row: ", nan_error),
    paste0(
      "the patched function for padj = 0 (no positive) gives the PR auc.integral ",
      no_positive$PRC$auc.integral, ", the PR auc.davis.goadrich ", no_positive$PRC$auc.davis.goadrich,
      ", and the ROC auc ", no_positive$ROC$auc
    ),
    paste0(
      "the patched function for padj = 1 and l2fc = Inf (no negative) gives the PR auc.integral ",
      no_negative$PRC$auc.integral, ", the PR auc.davis.goadrich ", no_negative$PRC$auc.davis.goadrich,
      ", and the ROC auc ", no_negative$ROC$auc
    ),
    paste0(
      "benchmark_na has the filter of benchmark_2 and a padj of NA for ", paste(na_genes, collapse = " and "),
      ". R makes these rows negatives: n_positive is ", n_positive_of[["benchmark_na"]], " instead of ",
      n_positive_of[["benchmark_2"]]
    )
  )
)

# ---- miRDB ----
# The lines of R/get_data.R:43-52, for a file.
mirdb_table <- function(path) {
  mirdb <- utils::read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  mirdb[, c("species", "miRNA.name")] <- stringr::str_split_fixed(mirdb[, 1], "-", 2)
  colnames(mirdb) <- c("miRDB.ID", "target.REFSEQ.ID", "target.score", "species", "miRNA.name")
  mirdb
}
mirdb_of_lines <- function(lines) {
  path <- tempfile(fileext = ".txt")
  on.exit(unlink(path))
  writeLines(lines, path)
  mirdb_table(path)
}

mirdb <- mirdb_table(file.path(DATA_DIR, "mirdb_sample.txt"))
stopifnot(nrow(mirdb) == 12, is.double(mirdb$target.score))
lin4 <- mirdb[mirdb$miRDB.ID == "lin4", ]
stopifnot(lin4$species == "lin4", lin4$miRNA.name == "")

text_score <- mirdb_of_lines(c("hsa-a\tNM_1\t50", "hsa-b\tNM_2\tabc"))
empty_score <- mirdb_of_lines(c("hsa-a\tNM_1\t50", "hsa-b\tNM_2\t"))
na_id <- mirdb_of_lines(c("NA\tNM_1\t50", "hsa-b\tNM_2\t10"))
short_line_error <- error_of(mirdb_of_lines(c("hsa-a\tNM_1\t50", "hsa-b\tNM_2")))
hash_error <- error_of(mirdb_of_lines(c("hsa-a\tNM_1\t50", "hsa-b\tNM#_2\t10")))
quote_lines <- c("hsa-a\tNM_1\t50", "hsa-b\tNM'_2\t10", "hsa-c\tNM_3\t11", "hsa-d\tNM_4\t12")
quote_table <- suppressWarnings(mirdb_of_lines(quote_lines))
stopifnot(
  is.character(text_score$target.score), is.na(empty_score$target.score[2]),
  is.na(na_id$miRDB.ID[1]), is.na(na_id$species[1])
)

write_analysis(
  mirdb, "mirdb.tsv",
  c(
    "the table of R/get_data.R:43-52 for tests/data/mirdb_sample.txt, with the R column names",
    paste0(
      "a score abc gives the column target.score the type ", class(text_score$target.score),
      ", with the values ", paste(text_score$target.score, collapse = " and ")
    ),
    paste0("an empty score gives ", empty_score$target.score[2]),
    paste0(
      "the ID NA gives the value ", na_id$miRDB.ID[1], " for the ID, the species ", na_id$species[1],
      ", and the miRNA name \"", na_id$miRNA.name[1], "\""
    ),
    paste0("a line with 2 fields stops: ", short_line_error),
    paste0("a line with a field with a hash character stops: ", hash_error),
    paste0(
      "a quote character in a field of line 2 of ", length(quote_lines), " lines gives ", nrow(quote_table),
      " rows with the IDs ", paste(quote_table$miRDB.ID, collapse = " and "),
      ", because the default quote of read.table includes the quote character"
    )
  )
)
