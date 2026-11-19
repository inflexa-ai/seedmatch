source("tools/golden/R/common.R")

data <- load_test_data()
report <- SeedMatchR::SeedMatchReport(res = data$res, seqs = data$seqs, guide.seq = data$guide)$report

row <- report[report$Category == "In silico predictions", ]
if (nrow(row) != 1 || row$Total != 125 || row$`15mer-D4` != 73) {
  stop("In silico predictions: Total ", row$Total, " and 15mer-D4 ", row$`15mer-D4`, ", expected 125 and 73")
}

seedmatchr <- as.character(packageVersion("SeedMatchR"))
biostrings <- as.character(packageVersion("Biostrings"))
if (seedmatchr != "2.0.0" || !startsWith(biostrings, "2.78.")) {
  stop("SeedMatchR ", seedmatchr, " and Biostrings ", biostrings, ", expected 2.0.0 and 2.78.x")
}

cat("golden environment OK: SeedMatchR ", seedmatchr, ", Biostrings ", biostrings, "\n", sep = "")
