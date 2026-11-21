source("tools/golden/R/common.R")
source("tools/golden/R/annotation_common.R")

write_reference_goldens(
  build = "rnor7",
  ensdb = "EnsDb.Rnorvegicus.v109.sqlite",
  twobit = "Rattus_norvegicus.mRatBN7.2.dna_sm.toplevel.2bit",
  configs = c("utr3", "utr5", "exons", "cds", "utr3_longest", "utr3_union"),
  script = "annotation_rnor7.R"
)
