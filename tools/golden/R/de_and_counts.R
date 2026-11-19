source("tools/golden/R/common.R")

data <- load_test_data()

cases <- utils::read.table(
  file.path(DATA_DIR, "de_filter_cases.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE
)
cutoff_sets <- list(c(1, 0, 30), c(0.05, 0.5, 10))
kept <- do.call(rbind, lapply(cutoff_sets, function(cutoffs) {
  filtered <- SeedMatchR::filter_res(
    cases, fdr_cutoff = cutoffs[1], fc_cutoff = cutoffs[2], basemean_cutoff = cutoffs[3]
  )
  data.frame(
    fdr_cutoff = cutoffs[1], fc_cutoff = cutoffs[2], basemean_cutoff = cutoffs[3],
    gene_id = filtered$gene_id, stringsAsFactors = FALSE
  )
}))
write_golden(kept, "de/filter_res.tsv", "de_and_counts.R")

seeds <- c("mer7m8", "mer8", "mer6", "mer7A1", "Full")
counts <- list()
notes <- character(0)
for (seed in seeds) {
  result <- SeedMatchR::SeedMatchR(
    seqs = data$seqs, sequence = data$guide, res = data$res, seed.name = seed
  )
  counts[[seed]] <- data.frame(
    gene_id = result$gene_id, seed = seed, count = result[[ncol(result)]], stringsAsFactors = FALSE
  )
  if (seed %in% colnames(data$res)) {
    notes <- c(notes, paste0("columns of the ", seed, " result: ", paste(names(result), collapse = ",")))
  }
}
write_golden(do.call(rbind, counts), "seedmatchr/deseq2_mode.tsv", "de_and_counts.R", notes)

absent_gene <- "ENSRNOG_NO_SEQUENCE"
absent_row <- data$res[1, ]
absent_row$gene_id <- absent_gene
res_with_absent <- rbind(data$res[1:50, ], absent_row, data$res[51:nrow(data$res), ])
by_sharing <- lapply(c(TRUE, FALSE), function(shared_genes) {
  result <- SeedMatchR::SeedMatchR(
    seqs = data$seqs, sequence = data$guide, res = res_with_absent, seed.name = "mer7m8",
    shared_genes = shared_genes
  )
  data.frame(
    shared_genes = shared_genes, gene_id = result$gene_id, count = result[[ncol(result)]],
    stringsAsFactors = FALSE
  )
})
write_golden(
  do.call(rbind, by_sharing), "seedmatchr/deseq2_mode_absent_gene.tsv", "de_and_counts.R",
  paste0("the gene ", absent_gene, " has no sequence and follows row 50 of the test DE table")
)
