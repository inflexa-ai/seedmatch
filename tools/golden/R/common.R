SEEDMATCHR_COMMIT <- "2fc0c4f"
GOLDEN_DIR <- "tests/golden"
DATA_DIR <- "tests/data"
GUIDE <- "UUAUAGAGCAAGAACACUGUUUU"

load_test_data <- function() {
  list(
    res = utils::read.table(file.path(DATA_DIR, "DESEQ2_results.tsv"), sep = "\t", header = 1),
    seqs = Biostrings::readDNAStringSet(file.path(DATA_DIR, "Test_seqs.fa")),
    guide = GUIDE
  )
}

# write.table cuts a double to 15 significant digits. 17 significant digits are exact for a correctly rounded reader such as Python float().
double_to_string <- function(x) {
  sprintf("%.17g", x)
}

# write_golden writes a NaN as text. A script that has no NaN in a golden file stops at one.
stop_on_nan <- function(df) {
  stopifnot(!any(vapply(df, function(col) is.double(col) && any(is.nan(col)), logical(1))))
}

# golden in tests/helpers.py splits a line at each tab, reads a line that starts with # as a note, and
# reads the text NA as a missing value. write_golden stops on a value that this reader changes.
stop_on_unreadable_text <- function(df, notes) {
  text <- lapply(df, function(col) if (is.character(col) || is.factor(col)) as.character(col)[!is.na(col)])
  columns_with <- function(test) names(Filter(function(col) any(test(col)), text))
  if (any(grepl("[\t\r\n]", names(df)))) stop("a column name holds a tab or a line break")
  broken <- columns_with(function(col) grepl("[\t\r\n]", col))
  if (length(broken) > 0) stop("a text cell of the columns ", toString(broken), " holds a tab or a line break")
  if (any(grepl("[\r\n]", notes))) stop("a note holds a line break")
  na_text <- columns_with(function(col) col == "NA")
  if (length(na_text) > 0) {
    stop("a text cell of the columns ", toString(na_text), " is the text NA, which the reader reads as a missing value; use NA")
  }
  starts <- c(names(df)[1], if (ncol(df) > 0) as.character(df[[1]]))
  if (any(startsWith(starts[!is.na(starts)], "#"))) stop("a line starts with #, which the reader reads as a note")
}

# install.R writes the commit of ../SeedMatchR into the installed package. Each golden file names
# SEEDMATCHR_COMMIT, thus write_golden stops if the installed package comes from another commit.
stop_on_other_commit <- function() {
  file <- system.file("COMMIT", package = "SeedMatchR")
  commit <- if (nzchar(file)) readLines(file, n = 1, warn = FALSE) else character(0)
  if (length(commit) == 0 || !startsWith(commit, SEEDMATCHR_COMMIT)) {
    stop(
      "the installed SeedMatchR has the commit ", if (length(commit) > 0) commit else "unknown",
      ", expected a commit that starts with ", SEEDMATCHR_COMMIT,
      "; run tools/golden/run.sh tools/golden/install.R"
    )
  }
}

write_golden <- function(df, path, script, notes = character(0)) {
  stop_on_other_commit()
  stop_on_unreadable_text(df, notes)
  df[] <- lapply(df, function(col) if (is.double(col) && !is.object(col)) double_to_string(col) else col)
  out <- file.path(GOLDEN_DIR, path)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  # The gzip header of gzfile has no time stamp and no file name, thus each run writes the same bytes.
  con <- if (endsWith(out, ".gz")) gzfile(out, "w") else file(out, "w")
  on.exit(close(con))
  writeLines(c(
    paste0("# generated-by: tools/golden/R/", script),
    paste0("# R: ", getRversion()),
    paste0("# Biostrings: ", packageVersion("Biostrings")),
    paste0("# SeedMatchR: ", packageVersion("SeedMatchR"), " (", SEEDMATCHR_COMMIT, ")"),
    if (length(notes) > 0) paste0("# ", notes)
  ), con)
  write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(out)
}

# Put each function of the patch file into the SeedMatchR namespace, evaluate code, and then put
# back the original functions. The patched functions get the namespace as the parent environment,
# thus they find the imports of SeedMatchR, for example %>%.
with_seedmatchr_patch <- function(patch, code) {
  ns <- asNamespace("SeedMatchR")
  patched <- new.env(parent = ns)
  sys.source(patch, envir = patched)
  names <- ls(patched, all.names = TRUE)
  originals <- mget(names, envir = ns)
  for (name in names) utils::assignInNamespace(name, patched[[name]], ns = "SeedMatchR")
  on.exit(for (name in names) utils::assignInNamespace(name, originals[[name]], ns = "SeedMatchR"))
  code
}

# A note is one comment line, thus a message of R with line breaks becomes one line.
one_line <- function(text) gsub("[[:space:]]+", " ", text)

# dplyr chains its errors. The message of the first error in the chain is the message of the cause.
root_message <- function(e) {
  while (!is.null(e$parent)) e <- e$parent
  one_line(conditionMessage(e))
}

# The message of the error of expr. The script stops if expr gives no error.
error_of <- function(expr) {
  message <- tryCatch({
    expr
    NULL
  }, error = root_message)
  if (is.null(message)) stop("an expected error did not occur")
  message
}
