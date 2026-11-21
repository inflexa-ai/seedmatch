# Measure the wall time and the memory of SeedMatchR::SeedMatchReport.
# The input is the FASTA file of tools/perf/utr3_fasta.py and the DE table of the R goldens.
# Run it from the repository root, with the file in PERF_FASTA and the number of calls in PERF_CALLS:
# PERF_FASTA=tmp/perf/rnor7_utr3.fa /usr/bin/time -l tools/golden/run.sh tools/perf/report.R
source("tools/golden/R/common.R")

# The resident memory of this process now. /usr/bin/time gives the peak of the process.
rss_mib <- function() {
  as.numeric(system(sprintf("ps -o rss= -p %d", Sys.getpid()), intern = TRUE)) / 1024
}

fasta <- Sys.getenv("PERF_FASTA")
calls <- as.integer(Sys.getenv("PERF_CALLS", "3"))
stopifnot(nzchar(fasta), calls >= 1)

seqs <- Biostrings::readDNAStringSet(fasta)
data <- load_test_data()
cat(sprintf("%d sequences; resident memory before the first call: %.0f MiB\n", length(seqs), rss_mib()))
for (call in seq_len(calls)) {
  started <- proc.time()[["elapsed"]]
  result <- SeedMatchR::SeedMatchReport(data$res, seqs, data$guide)
  seconds <- proc.time()[["elapsed"]] - started
  cat(sprintf("call %d: %.2f s, resident memory after the call: %.0f MiB\n", call, seconds, rss_mib()))
}
report <- result$report
cat("In silico predictions, Total:", report$Total[report$Category == "In silico predictions"], "\n")
versions <- c("SeedMatchR", "Biostrings", "dplyr", "gt")
cat(paste(versions, vapply(versions, function(p) as.character(packageVersion(p)), "")), sep = ", ")
cat(", R", as.character(getRversion()), "\n")
