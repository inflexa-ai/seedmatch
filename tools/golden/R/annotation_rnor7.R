if (!requireNamespace("openssl", quietly = TRUE)) {
  stop("the R package openssl is not installed; annotation_rnor7.R needs it for the digests")
}

source("tools/golden/R/common.R")
source("tools/golden/R/annotation_common.R")

SCRIPT <- "annotation_rnor7.R"
ENSDB <- "EnsDb.Rnorvegicus.v109.sqlite"
TWOBIT <- "Rattus_norvegicus.mRatBN7.2.dna_sm.toplevel.2bit"

refs <- Sys.getenv("SEEDMATCH_TEST_REFS")
if (!nzchar(refs)) stop("SEEDMATCH_TEST_REFS is not set; set it to the directory that holds rnor7/")
rnor7 <- file.path(refs, "rnor7")

txdb <- ensembldb::EnsDb(file.path(rnor7, ENSDB))
dna <- rtracklayer::TwoBitFile(file.path(rnor7, TWOBIT))
utils::assignInNamespace("load_species_anno_db", local_anno_db(txdb, dna), ns = "SeedMatchR")

# The feature sequences of rnor7 are too large for a golden file. tests/annotation/test_refs.py
# computes the same digest from the uppercase text of each sequence.
digests_table <- function(result) {
  digests <- sequences_table(result)
  digests$digest <- substr(as.character(openssl::sha256(digests$sequence)), 1, 16)
  digests$sequence <- NULL
  digests
}

for (name in c("utr3", "utr5", "exons", "cds", "utr3_longest", "utr3_union")) {
  started <- Sys.time()
  run <- run_config(name, ANNOTATION_CONFIGS[[name]])
  digests <- digests_table(run$result)
  notes <- c(run$notes, paste0("reference: rnor7/", ENSDB, ", rnor7/", TWOBIT))
  write_golden(digests, paste0("annotation/rnor7/", name, ".tsv.gz"), SCRIPT, notes)
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  cat(name, ": ", nrow(digests), " rows, ", round(seconds), " s\n", sep = "")
}
