source("tools/golden/R/common.R")
suppressPackageStartupMessages(library(Biostrings))

iupac <- names(IUPAC_CODE_MAP)
grid <- expand.grid(
  subject = iupac,
  pattern = iupac,
  fixed_subject = c(TRUE, FALSE),
  fixed_pattern = c(TRUE, FALSE),
  stringsAsFactors = FALSE
)
grid <- grid[, c("pattern", "subject", "fixed_pattern", "fixed_subject")]
grid$match <- mapply(
  function(p, s, fp, fs) countPattern(p, DNAString(s), fixed = c(pattern = fp, subject = fs)) > 0,
  grid$pattern, grid$subject, grid$fixed_pattern, grid$fixed_subject,
  USE.NAMES = FALSE
)

write_golden(grid, "search/iupac_compat.tsv", "iupac_table.R")
