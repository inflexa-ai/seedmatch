source("tools/golden/R/common.R")

suppressPackageStartupMessages(library(Matrix))

# R quotes a name in an error message with curly quotes in a UTF-8 locale. The golden files hold ASCII only.
options(useFancyQuotes = FALSE)

TRRUST_SAMPLE <- file.path(DATA_DIR, "trrust_sample.tsv")
STRING_SAMPLE <- file.path(DATA_DIR, "string_links_sample.txt")
ENSDB <- file.path(DATA_DIR, "annotation", "synthetic.sqlite")
TRRUST_COLUMNS <- c("source", "target", "edge_type", "evidence")
STRING_THRESHOLD <- 400
STRING_HEADER <- "protein1 protein2 combined_score"
VERSIONS <- paste0(
  "igraph: ", packageVersion("igraph"), ", Matrix: ", packageVersion("Matrix"),
  ", DBI: ", packageVersion("DBI"), ", RSQLite: ", packageVersion("RSQLite")
)

write_network <- function(df, path, notes = character(0)) {
  stop_on_nan(df)
  write_golden(df, file.path("network", path), "network.R", c(VERSIONS, notes))
}

same_double <- function(x, expected) identical(double_to_string(x), expected)
gene <- function(n) sprintf("ENSSYNG%011d", n)
protein <- function(n) sprintf("ENSSYNP%011d", n)

# ---- TRRUST rows ----
trrust_rows <- utils::read.table(
  TRRUST_SAMPLE, header = FALSE, sep = "\t", colClasses = "character", quote = "", comment.char = ""
)
colnames(trrust_rows) <- TRRUST_COLUMNS
stopifnot(
  nrow(trrust_rows) == 10,
  identical(unlist(trrust_rows[1, ], use.names = FALSE), c("TFA", "GA", "Repression", "1001")),
  identical(trrust_rows$evidence[9], "1009;1010")
)

# The call of R/get_network_matches.R:25.
read_trrust_like_r <- function(path) {
  utils::read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

header_true <- read_trrust_like_r(TRRUST_SAMPLE)
header_true_names <- names(header_true)
# The names of R/get_network_matches.R:30.
colnames(header_true) <- TRRUST_COLUMNS
stopifnot(nrow(header_true) == 9, identical(header_true_names, c("TFA", "GA", "Repression", "X1001")))

# The result of the call for a file with the line of a regulator and then the lines of a case.
trrust_case <- function(lines) {
  path <- file.path(tempdir(), "trrust_case.tsv")
  writeLines(c("TF1\tG1\tRepression\t1", lines), path)
  tryCatch(read_trrust_like_r(path), error = root_message)
}

quote_lines <- c("TF'2\tG2\tActivation\t2", "TF3\tG3\tUnknown\t3", "TF4\tG'4\tActivation\t4")
hash_lines <- c("TF#2\tG2\tActivation\t2", "TF3\tG3\tUnknown\t3")
fields_lines <- c("TF2\tG2\tActivation", "TF3\tG3\tUnknown\t3")
na_lines <- c("NA\tG2\tActivation\t2", "TF3\tNA\tUnknown\t3")
quote_result <- trrust_case(quote_lines)
hash_result <- trrust_case(hash_lines)
fields_result <- trrust_case(fields_lines)
na_result <- trrust_case(na_lines)
colnames(na_result) <- TRRUST_COLUMNS
four_elements <- "line 1 did not have 4 elements"
stopifnot(
  identical(quote_result, four_elements),
  identical(hash_result, four_elements),
  identical(fields_result, four_elements),
  identical(which(is.na(na_result$source)), 1L),
  identical(which(is.na(na_result$target)), 2L),
  identical(na_result$source[2], "TF3")
)

# A note shows each tab of a line as \t.
shown <- function(lines) paste(gsub("\t", "\\t", lines, fixed = TRUE), collapse = " | ")
case_note <- function(id, lines, result) {
  paste0(
    id, ": the same call on a file with the lines ", shown(c("TF1\tG1\tRepression\t1", lines)),
    " stops with: ", result
  )
}

write_network(
  trrust_rows, "trrust_rows.tsv",
  c(
    paste0(
      "header_true: the call of R/get_network_matches.R:25 on trrust_sample.tsv gives ", nrow(header_true),
      " rows, and the names before :30 are ", paste(header_true_names, collapse = " | ")
    ),
    case_note("quote", quote_lines, quote_result),
    case_note("hash", hash_lines, hash_result),
    case_note("three_fields", fields_lines, fields_result),
    paste0(
      "na_text: the same call on a file with the lines ", shown(c("TF1\tG1\tRepression\t1", na_lines)),
      " gives the value NA, not the text: source NA in data row ", which(is.na(na_result$source)),
      " and target NA in data row ", which(is.na(na_result$target))
    )
  )
)

# ---- TRRUST matches ----
network_matches <- function(seeds, aggregator, network) {
  SeedMatchR::get_network_matches(seeds, organism = "human", aggregator = aggregator, network_data = network)
}

SEEDS <- c("TFA", "GC", "TFC")
NO_SEEDS <- c("ENSG00000000001", "ENSG00000000002")
match_cases <- list(
  mean = list(seeds = SEEDS, aggregator = "mean", network = trrust_rows),
  degree_importance = list(seeds = SEEDS, aggregator = "degree_importance", network = trrust_rows),
  mean_header_true = list(seeds = SEEDS, aggregator = "mean", network = header_true),
  degree_importance_header_true = list(
    seeds = SEEDS, aggregator = "degree_importance", network = header_true
  ),
  mean_no_seed = list(seeds = NO_SEEDS, aggregator = "mean", network = trrust_rows)
)
matches <- do.call(rbind, lapply(names(match_cases), function(id) {
  case <- match_cases[[id]]
  cbind(case = id, network_matches(case$seeds, case$aggregator, case$network))
}))

match_value <- function(case, gene, column) {
  matches[[column]][matches$case == case & matches$ENSG == gene]
}
no_seed <- matches[matches$case == "mean_no_seed", ]
stopifnot(
  match_value("mean", "GA", "repression_matches") == 1,
  match_value("mean", "GA", "activation_matches") == 0,
  match_value("mean", "GC", "repression_matches") == 1,
  match_value("mean", "GC", "activation_matches") == 1,
  same_double(match_value("degree_importance", "TFA", "repression_matches"), "0.20000000000000001"),
  same_double(match_value("degree_importance", "TFA", "activation_matches"), "0.040000000000000001"),
  match_value("degree_importance_header_true", "TFA", "repression_matches") == 0.25,
  match_value("degree_importance_header_true", "TFA", "activation_matches") == 0.0625,
  nrow(no_seed) == 7,
  all(no_seed$repression_matches == 0),
  all(no_seed$activation_matches == 0)
)

# The graph of R/get_network_matches.R:35.
trrust_graph <- igraph::graph_from_data_frame(trrust_rows[, 1:2], directed = TRUE)
tfa_targets <- names(igraph::neighbors(trrust_graph, "TFA"))
median_matches <- network_matches(SEEDS, "median", trrust_rows)
median_flags <- as.numeric(median_matches$ENSG %in% SEEDS)
organism_error <- error_of(SeedMatchR::get_network_matches("TFA", network_data = trrust_rows))
stopifnot(
  identical(tfa_targets, c("GA", "GA", "GB", "GB")),
  identical(median_matches$ENSG, names(igraph::V(trrust_graph))),
  identical(unname(median_matches$repression_matches), median_flags),
  identical(unname(median_matches$activation_matches), median_flags),
  identical(organism_error, "argument \"organism\" is missing, with no default")
)

write_network(
  matches, "trrust_matches.tsv",
  c(
    paste0(
      "out_neighbors_tfa: names(igraph::neighbors(g, \"TFA\")) for the graph of R/get_network_matches.R:35 ",
      "gives ", paste(tfa_targets, collapse = " ")
    ),
    paste0(
      "median: get_network_matches(seeds, organism = \"human\", aggregator = \"median\", network_data = ...) ",
      "runs no branch, thus each value is the seed flag: ",
      paste(median_matches$ENSG, median_flags, collapse = ", ")
    ),
    paste0("organism_missing: get_network_matches without organism stops with: ", organism_error)
  )
)

# ---- Quantiles ----
string_frame <- utils::read.table(
  STRING_SAMPLE, sep = " ", header = TRUE, stringsAsFactors = FALSE, fill = TRUE
)
link_scores <- string_frame$combined_score[string_frame$combined_score >= STRING_THRESHOLD]
stopifnot(length(link_scores) == 22)

quantile_cases <- list(
  seq21_p95 = list(values = 1:21, probability = 0.95, expected = "20"),
  ties_p30 = list(values = c(5, 5, 7), probability = 0.3, expected = "5"),
  one_p95 = list(values = 1, probability = 0.95, expected = "1"),
  ties_high_p95 = list(
    values = c(400, 940, 940, 999), probability = 0.95, expected = "990.14999999999998"
  ),
  two_third = list(values = c(1, 2), probability = 1 / 3, expected = "1.3333333333333333"),
  seq10_p10 = list(values = 1:10, probability = 0.1, expected = "1.8999999999999999"),
  links_p0 = list(values = link_scores, probability = 0, expected = "400"),
  links_p50 = list(values = link_scores, probability = 0.5, expected = "650"),
  links_p95 = list(values = link_scores, probability = 0.95, expected = "897.49999999999989"),
  max_p100 = list(values = c(3, 9, 4), probability = 1, expected = "9")
)
quantile_rows <- do.call(rbind, lapply(names(quantile_cases), function(id) {
  case <- quantile_cases[[id]]
  result <- stats::quantile(case$values, case$probability, names = FALSE)
  stopifnot(same_double(result, case$expected))
  data.frame(
    case = id, values = paste(case$values, collapse = ","), probability = case$probability,
    quantile = result
  )
}))

write_network(quantile_rows, "r_quantiles.tsv")

# ---- STRING links ----
# get_stringPPI_to_stringGGI runs on the local files. The functions local_stringdb_new, local_use_mart,
# and local_get_bm replace STRINGdb$new, useMart, and getBM. They give what STRINGdb 2.22.0 and biomaRt
# give for the defaults of the function. They do not copy the source of STRINGdb.
# The read.table call of load (R/rstring.R:1032 of STRINGdb 2.22.0) is at :1064.
read_string_like_r <- function(path) {
  utils::read.table(path, sep = " ", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
}
# STRINGdb keeps the rows with a score at or above the threshold (:1067), and the default is 400 (:109-110).
keep_default_threshold <- function(links) links[links$combined_score >= STRING_THRESHOLD, ]

# get_graph (R/rstring.R:692) calls load, and load returns an undirected graph from the kept rows (:1073).
local_stringdb_new <- function(species) {
  list(get_graph = function() {
    kept <- keep_default_threshold(read_string_like_r(STRING_SAMPLE))
    igraph::graph_from_data_frame(kept, FALSE)
  })
}

local_use_mart <- function(...) NULL

local_get_bm <- function(attributes, filters, values, mart) {
  con <- DBI::dbConnect(RSQLite::SQLite(), ENSDB, flags = RSQLite::SQLITE_RO)
  on.exit(DBI::dbDisconnect(con))
  mapping <- DBI::dbGetQuery(con, paste(
    "SELECT tx.gene_id AS ensembl_gene_id, protein.protein_id AS ensembl_peptide_id",
    "FROM protein JOIN tx ON tx.tx_id = protein.tx_id ORDER BY protein.protein_id"
  ))
  mapping[mapping$ensembl_peptide_id %in% values, ]
}

replace_text <- function(source, text, replacement, count) {
  stopifnot(lengths(regmatches(source, gregexpr(text, source, fixed = TRUE))) == count)
  gsub(text, replacement, source, fixed = TRUE)
}

patched_source <- paste(deparse(SeedMatchR::get_stringPPI_to_stringGGI), collapse = "\n")
patched_source <- replace_text(patched_source, "STRINGdb::STRINGdb$new", "local_stringdb_new", 1)
patched_source <- replace_text(patched_source, "biomaRt::useMart", "local_use_mart", 2)
patched_source <- replace_text(patched_source, "biomaRt::getBM", "local_get_bm", 1)
string_to_genes <- eval(parse(text = patched_source))
environment(string_to_genes) <- globalenv()

# igraph 2.3.3 warns that subgraph.edges is deprecated. The script keeps the message as a note.
subgraph_warnings <- character(0)
without_subgraph_warning <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("subgraph.edges", conditionMessage(w), fixed = TRUE)) {
      subgraph_warnings <<- c(subgraph_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  })
}

string_selected <- keep_default_threshold(read_string_like_r(STRING_SAMPLE))
string_graph <- local_stringdb_new(32630)$get_graph()
edge_scores <- igraph::E(string_graph)$combined_score
stopifnot(
  nrow(string_selected) == 22,
  igraph::ecount(string_graph) == 22,
  !igraph::is_directed(string_graph)
)

# The entries above 0 of a matrix, in the order of the row and then of the column.
matrix_entries <- function(ppi) {
  entries <- Matrix::summary(methods::as(ppi, "generalMatrix"))
  entries <- entries[entries$x > 0, ]
  entries <- entries[order(entries$i, entries$j), ]
  data.frame(
    gene1 = rownames(ppi)[entries$i], gene2 = colnames(ppi)[entries$j], value = entries$x,
    stringsAsFactors = FALSE
  )
}

string_cases <- list(p0 = 0, p50 = 0.5, p95 = 0.95)
expected_cutoffs <- c(p0 = "400", p50 = "650", p95 = "897.49999999999989")
expected_counts <- list(p0 = c(20, 2), p50 = c(10, 2), p95 = c(2, 0))
expected_nodes <- list(
  p0 = c(gene(c(1, 2, 4, 5, 6, 7, 8)), protein(c(98, 99))),
  p50 = gene(c(1, 2, 4, 6)),
  p95 = gene(c(1, 4))
)
expected_entries <- c(p0 = 12, p50 = 4, p95 = 2)

cutoff_rows <- list()
link_rows <- list()
node_rows <- list()
entry_rows <- list()
ppi_of <- list()
for (id in names(string_cases)) {
  p <- string_cases[[id]]
  cutoff <- unname(stats::quantile(edge_scores, p))
  above <- string_selected[string_selected$combined_score > cutoff, ]
  subgraph <- without_subgraph_warning(
    igraph::subgraph.edges(string_graph, igraph::E(string_graph)[combined_score > cutoff])
  )
  ppi <- without_subgraph_warning(string_to_genes(32630, p = p, sp = "rat"))
  entries <- matrix_entries(ppi)
  stopifnot(
    same_double(cutoff, expected_cutoffs[[id]]),
    identical(c(sum(edge_scores > cutoff), sum(edge_scores == cutoff)), as.integer(expected_counts[[id]])),
    nrow(above) == igraph::ecount(subgraph),
    identical(rownames(ppi), expected_nodes[[id]]),
    nrow(entries) == expected_entries[[id]],
    all(entries$value == 2)
  )
  cutoff_rows[[id]] <- data.frame(
    case = id, percentile = p, cutoff = cutoff, n_above = sum(edge_scores > cutoff),
    n_equal = sum(edge_scores == cutoff)
  )
  link_rows[[id]] <- cbind(case = id, above[c("protein1", "protein2", "combined_score")])
  node_rows[[id]] <- data.frame(case = id, position = seq_along(rownames(ppi)), name = rownames(ppi))
  entry_rows[[id]] <- cbind(case = id, entries)
  ppi_of[[id]] <- ppi
}
stopifnot(length(subgraph_warnings) >= 1)
subgraph_message <- one_line(strsplit(subgraph_warnings[1], "\n", fixed = TRUE)[[1]][1])
stopifnot(grepl("^`subgraph.edges\\(\\)` was deprecated in igraph", subgraph_message))

# The errors of the percentile.
percentile_1_error <- error_of(without_subgraph_warning(string_to_genes(32630, p = 1, sp = "rat")))
percentile_high_error <- error_of(stats::quantile(edge_scores, 1.2))
percentile_na_quantile <- stats::quantile(edge_scores, NA)
percentile_na_error <- error_of(without_subgraph_warning(string_to_genes(32630, p = NA, sp = "rat")))
stopifnot(
  identical(percentile_1_error, "non-character argument"),
  identical(percentile_high_error, "'probs' outside [0,1]"),
  is.na(percentile_na_quantile),
  identical(percentile_na_error, "Unknown edge selected")
)

# The read of a file with the header, the good rows, and then the lines of a case. The warning of
# read.table names the file, thus the script replaces the path with <file>.
read_string_case <- function(before, lines) {
  path <- file.path(tempdir(), "string_case.txt")
  writeLines(c(STRING_HEADER, before, lines), path)
  warnings <- character(0)
  links <- withCallingHandlers(
    read_string_like_r(path),
    warning = function(w) {
      warnings <<- c(warnings, gsub(path, "<file>", conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning")
    }
  )
  has_row_names <- is.character(attr(links, "row.names"))
  with_row_names <- function(rows, frame) {
    if (has_row_names) paste0("row name ", rownames(frame), ", ", rows) else rows
  }
  text <- function(x) ifelse(is.na(x), "NA", as.character(x))
  # A text value shows its quotes and its escapes, because a garbled value can hold a line end.
  cell <- function(x) if (is.character(x)) ifelse(is.na(x), "NA", encodeString(x, quote = "\"")) else text(x)
  kept <- keep_default_threshold(links)
  after <- links[-seq_along(before), , drop = FALSE]
  list(
    class = class(links$combined_score), n_rows = nrow(links),
    after = with_row_names(paste(cell(after$protein1), cell(after$protein2), cell(after$combined_score)), after),
    kept = with_row_names(paste(text(kept$protein1), text(kept$protein2), text(kept$combined_score)), kept),
    warnings = warnings
  )
}

good_rows <- c("9.A 9.B 500", "9.B 9.A 500")
late_good_rows <- c(
  good_rows, "9.G 9.H 500", "9.H 9.G 500", "9.I 9.J 500", "9.J 9.I 500", "9.K 9.L 500", "9.L 9.K 500"
)
na_row <- "NA NA NA"
read_cases <- list(
  text_score = list(
    lines = c("9.C 9.D abc", "9.D 9.C 1000", "9.E 9.F 99"), class = "character", n_rows = 5,
    after = c("\"9.C\" \"9.D\" \"abc\"", "\"9.D\" \"9.C\" \"1000\"", "\"9.E\" \"9.F\" \"99\""),
    kept = c(good_rows, "9.C 9.D abc", "9.E 9.F 99")
  ),
  two_fields = list(
    lines = "9.C 9.D", class = "integer", n_rows = 3, after = "\"9.C\" \"9.D\" NA",
    kept = c(good_rows, na_row)
  ),
  hash = list(
    lines = "9.C#x 9.D 700", class = "integer", n_rows = 3, after = "\"9.C\" \"\" NA",
    kept = c(good_rows, na_row)
  ),
  quote_pair = list(
    lines = c("9.C'x 9.D 700", "9.D 9.C'x 700"), class = "integer", n_rows = 3,
    after = "\"9.Cx 9.D 700\\n9.D 9.Cx\" \"700\" NA", kept = c(good_rows, na_row)
  ),
  quote_single = list(
    lines = "9.C'x 9.D 700", class = "logical", n_rows = 0, after = character(0), kept = character(0),
    warnings = "incomplete final line found by readTableHeader on '<file>'"
  ),
  four_fields = list(
    lines = "9.C 9.D 700 9.E", class = "character", n_rows = 3,
    after = "row name 9.C, \"9.D\" 700 \"9.E\"", kept = "row name 9.C, 9.D 700 9.E"
  ),
  negative_and_large = list(
    lines = c("9.C 9.D -5", "9.D 9.C 5000"), class = "integer", n_rows = 4,
    after = c("\"9.C\" \"9.D\" -5", "\"9.D\" \"9.C\" 5000"), kept = c(good_rows, "9.D 9.C 5000")
  ),
  decimal = list(
    lines = "9.C 9.D 450.5", class = "numeric", n_rows = 3, after = "\"9.C\" \"9.D\" 450.5",
    kept = c(good_rows, "9.C 9.D 450.5")
  ),
  quote_unclosed_late = list(
    before = late_good_rows, lines = c("9.C'x 9.D 700", "9.E 9.F 800"), class = "integer", n_rows = 9,
    after = "\"9.Cx 9.D 700\\n9.E 9.F 800\\n\" \"\" NA", kept = c(late_good_rows, na_row),
    warnings = "EOF within quoted string"
  ),
  four_fields_late = list(
    before = late_good_rows, lines = c("9.C 9.D 700 9.E", "9.F 9.G 800"), class = "integer", n_rows = 11,
    after = c("\"9.C\" \"9.D\" 700", "\"9.E\" \"\" NA", "\"9.F\" \"9.G\" 800"),
    kept = c(late_good_rows, "9.C 9.D 700", na_row, "9.F 9.G 800")
  )
)
read_notes <- character(0)
for (id in names(read_cases)) {
  case <- read_cases[[id]]
  before <- if (is.null(case$before)) good_rows else case$before
  result <- read_string_case(before, case$lines)
  stopifnot(
    identical(result$class, case$class),
    result$n_rows == case$n_rows,
    identical(result$after, case$after),
    identical(result$kept, case$kept),
    identical(result$warnings, if (is.null(case$warnings)) character(0) else case$warnings)
  )
  shown_rows <- function(rows) if (length(rows) == 0) "none" else paste(rows, collapse = " | ")
  read_notes <- c(read_notes, paste0(
    "read_", id, ": the read.table call of STRINGdb 2.22.0 (R/rstring.R:1064) with the lines ",
    paste(before, collapse = " | "), " and then ", paste(case$lines, collapse = " | "), " gives ",
    result$n_rows, " rows and the class ", result$class, " for combined_score. The rows after the first ",
    length(before), " lines are: ", shown_rows(result$after), ". The subset combined_score >= ",
    STRING_THRESHOLD, " (R/rstring.R:1067) keeps the rows: ", shown_rows(result$kept),
    if (length(result$warnings) > 0) paste0(". The warning is: ", result$warnings) else ""
  ))
}

cutoff_notes <- c(
  paste0("subgraph_edges: igraph ", packageVersion("igraph"), " warns with: ", subgraph_message),
  paste0(
    "percentile_1: get_stringPPI_to_stringGGI with p = 1 stops with: ", percentile_1_error,
    ", because no edge is above the cutoff ", double_to_string(max(edge_scores))
  ),
  paste0("percentile_1_2: stats::quantile(edge.scores, 1.2) stops with: ", percentile_high_error),
  paste0(
    "percentile_na: stats::quantile(edge.scores, NA) gives NA, and get_stringPPI_to_stringGGI with p = NA ",
    "stops with: ", percentile_na_error
  )
)
write_network(do.call(rbind, cutoff_rows), "string_cutoffs.tsv", cutoff_notes)
write_network(do.call(rbind, link_rows), "string_links.tsv", read_notes)
write_network(do.call(rbind, node_rows), "string_nodes.tsv")
write_network(do.call(rbind, entry_rows), "string_ppi.tsv")

# ---- Seed counts ----
node_names <- rownames(ppi_of$p0)
feature_cases <- list(
  repeated = list(
    res = data.frame(gene_id = gene(c(1, 2, 2, 99)), overall_matches = c(2, 1, 3, 5)),
    expected = setNames(c(2, 1), gene(c(1, 2)))
  ),
  na_count = list(
    res = data.frame(gene_id = gene(6), overall_matches = NA_real_),
    expected = setNames(NA_real_, gene(6))
  ),
  bad_counts = list(
    res = data.frame(gene_id = gene(c(1, 4)), overall_matches = c(-1, Inf)),
    expected = setNames(c(-1, Inf), gene(c(1, 4)))
  ),
  none_in_network = list(
    res = data.frame(gene_id = gene(99), overall_matches = 4), expected = numeric(0)
  )
)

feature_rows <- list()
feature_warnings <- list()
for (id in names(feature_cases)) {
  case <- feature_cases[[id]]
  warnings <- character(0)
  graph <- withCallingHandlers(
    SeedMatchR::create_graph(ppi_of$p0, case$res),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expected <- setNames(rep(0, length(node_names)), node_names)
  expected[names(case$expected)] <- case$expected
  stopifnot(
    identical(names(igraph::V(graph)), node_names),
    identical(unname(igraph::V(graph)$feature), unname(expected))
  )
  feature_rows[[id]] <- data.frame(case = id, name = node_names, feature = igraph::V(graph)$feature)
  feature_warnings[[id]] <- warnings
}

# The graph of the case repeated.
repeated_graph <- suppressWarnings(SeedMatchR::create_graph(ppi_of$p0, feature_cases$repeated$res))
repeated_neighbors <- names(igraph::neighbors(repeated_graph, gene(1)))
no_counts_error <- error_of(SeedMatchR::create_graph(ppi_of$p0, data.frame(gene_id = gene(1))))
repeated_warning <- "number of items to replace is not a multiple of replacement length"
stopifnot(
  identical(feature_warnings$repeated, repeated_warning),
  all(lengths(feature_warnings[names(feature_warnings) != "repeated"]) == 0),
  igraph::is_directed(repeated_graph),
  igraph::ecount(repeated_graph) == 24,
  identical(repeated_neighbors, c(gene(4), gene(4))),
  identical(no_counts_error, "replacement has length zero")
)

write_network(
  do.call(rbind, feature_rows), "string_features.tsv",
  c(
    paste0("repeated: create_graph warns with: ", repeated_warning, ". It keeps the first count of a gene"),
    paste0(
      "repeated_graph: igraph::is_directed ", igraph::is_directed(repeated_graph), ", igraph::ecount ",
      igraph::ecount(repeated_graph), ", and the neighbors of ", gene(1), " are ",
      paste(repeated_neighbors, collapse = " ")
    ),
    paste0(
      "no_counts: create_graph(ppi, data.frame(gene_id = \"", gene(1), "\")), a res without the column ",
      "overall_matches for a gene of the network, stops with: ", no_counts_error
    )
  )
)

# The golden files hold ASCII only.
for (file in list.files(file.path(GOLDEN_DIR, "network"), full.names = TRUE)) {
  bytes <- readBin(file, "raw", file.size(file))
  stopifnot(all(as.integer(bytes) < 128))
}
