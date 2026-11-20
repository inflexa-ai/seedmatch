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
