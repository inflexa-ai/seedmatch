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

write_golden <- function(df, path, script, notes = character(0)) {
  df[] <- lapply(df, function(col) if (is.double(col) && !is.object(col)) double_to_string(col) else col)
  out <- file.path(GOLDEN_DIR, path)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  con <- file(out, "w")
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
