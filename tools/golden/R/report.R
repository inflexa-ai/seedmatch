source("tools/golden/R/common.R")

data <- load_test_data()

# SeedMatchReport calls full_search, which counts with vcountPattern, so the original package
# runs and no patch of vmatchPattern2 is necessary.
# The wobble target has R at each U of the guide, and the default fixed = TRUE matches R only to R
# (D-026). fixed = "subject" is fixed_pattern=False, fixed_subject=True (Biostrings normargFixed).
variants <- list(
  default = list(),
  wobble = list(allow_wobbles = TRUE, fixed = "subject")
)
for (variant in names(variants)) {
  args <- c(list(res = data$res, seqs = data$seqs, guide.seq = data$guide), variants[[variant]])
  result <- do.call(SeedMatchR::SeedMatchReport, args)
  write_golden(result$report, paste0("seedmatchr/report_", variant, ".tsv"), "report.R")
  write_golden(result$matches, paste0("seedmatchr/report_matches_", variant, ".tsv"), "report.R")
}
