source("tools/golden/R/common.R")
source("tools/golden/R/annotation_common.R")

write_reference_goldens(
  build = "hg38",
  ensdb = "EnsDb.Hsapiens.v109.sqlite",
  twobit = "Homo_sapiens.GRCh38.dna.primary_assembly.2bit",
  configs = c("utr3", "utr5", "exons", "cds", "utr3_longest", "utr3_union", "entrez", "tsl1"),
  script = "annotation_hg38.R"
)
