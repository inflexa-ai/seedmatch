source("tools/golden/R/common.R")
source("tools/golden/patches/vmatchPattern2.R")

data <- load_test_data()

patch_note <- "vmatchPattern2 is the patched function of tools/golden/patches/vmatchPattern2.R"

variants <- list(
  default = list(),
  wobble = list(allow_wobbles = TRUE),
  indels = list(indel.bool = TRUE)
)
for (variant in names(variants)) {
  args <- c(list(guide.seq = data$guide, seqs = data$seqs), variants[[variant]])
  counts <- do.call(SeedMatchR::full_search, args)
  write_golden(counts, paste0("seedmatchr/full_search_counts_", variant, ".tsv"), "full_search.R")
}

# full_search_ranges and the data.frame mode of SeedMatchR look up vmatchPattern2 in the namespace.
# The counts above do not call it, so they come from the original package.
assignInNamespace("vmatchPattern2", vmatchPattern2, ns = "SeedMatchR")

ranges <- SeedMatchR::full_search_ranges(data$guide, data$seqs)
write_golden(ranges, "seedmatchr/full_search_ranges.tsv", "full_search.R", patch_note)

# The calls of tests/testthat/helper.R:49 and :51 in SeedMatchR. The second call puts its rows first.
sample_seqs <- data$seqs[names(data$seqs) %in% "ENSRNOG00000016275"]
partial <- SeedMatchR::SeedMatchR(
  sequence = data$guide, seqs = sample_seqs, res.format = "data.frame",
  seed.name = "Partial", start.pos = 2, stop.pos = 19, max.mismatch = 3
)
both <- SeedMatchR::SeedMatchR(
  sequence = data$guide, seqs = sample_seqs, res.format = "data.frame",
  seed.name = "mer7m8", max.mismatch = 1, match.df = partial
)
write_golden(both, "seedmatchr/dataframe_mode.tsv", "full_search.R", patch_note)
