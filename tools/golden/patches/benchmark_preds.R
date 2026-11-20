# A patched benchmark_preds for the golden files. The version in SeedMatchR (R/benchmark.R:23)
# filters with padj <= padj. Both sides are the column padj, thus the filter ignores the argument
# padj and keeps each row that has a padj (D-006). This version compares the column with the
# argument (padj <= !!padj). It also gives roc.curve the labels as weights.class0, as pr.curve gets
# them (R/benchmark.R:37 passes them as scores.class1). The other lines are the same as in the
# original.
benchmark_preds <- function(res,
                            preds_col = "mer7m8",
                            padj = .15,
                            l2fc = 0) {
  de_genes <- res %>%
    dplyr::filter(padj <= !!padj, log2FoldChange < l2fc) %>%
    dplyr::select(gene_id) %>%
    dplyr::distinct() %>%
    dplyr::pull()

  res$GS.label <- ifelse(res$gene_id %in% de_genes, 1, 0)

  preds <- res %>% dplyr::select(preds_col) %>% dplyr::pull()

  pr <- PRROC::pr.curve(scores.class0 = preds,
                        weights.class0 = res$GS.label,
                        curve = TRUE)

  roc <- PRROC::roc.curve(scores.class0 = preds,
                          weights.class0 = res$GS.label,
                          curve = TRUE)

  return(list("PRC" = pr,
              "ROC" = roc))
}
