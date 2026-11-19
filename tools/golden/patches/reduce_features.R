# A patched .reduce_features for the golden files. The version in SeedMatchR
# (R/load_annotations.R:183-213) gives the ranges of each gene in the order of the position, as
# GenomicRanges::reduce gives them. On the minus strand, extractTranscriptSeqs then joins the ranges
# in the wrong order (D-027). This version reverses the ranges of each element on the minus strand
# after reduce, thus the ranges come in transcript order. The other lines are the same as in the
# original.
.reduce_features <- function(gr.list, txdb, filter){
  # Create a df to map tx_id, tx_name, and gene_id
  tx2gene <- S4Vectors::mcols(GenomicFeatures::transcripts(txdb,
                                                           columns=c("tx_id",
                                                                     "tx_name",
                                                                     "gene_id"),
                                                           filter = filter
  ))

  tx2gene$gene_id <- as.character(tx2gene$gene_id)

  # Unlist the GTF so that they are not longer grouped by transcript
  all_features <- BiocGenerics::unlist(gr.list)

  ## Add "tx_name" and "gene_id" metadata columns to feature specific gtf
  m <- match(names(all_features), tx2gene$tx_name)

  mcol.df <- S4Vectors::mcols(all_features)

  S4Vectors::mcols(all_features) <- cbind(mcol.df, tx2gene[m, -1L, drop=FALSE])

  names(all_features) <- NULL

  # Split the GTF into a list by gene_id
  features_by_genes <- GenomicRanges::split(all_features, all_features$gene_id)

  # Reduce the overlapping genomic ranges to only have unique intervals
  gtf = GenomicRanges::reduce(features_by_genes)

  minus <- any(BiocGenerics::strand(gtf) == "-")
  gtf <- S4Vectors::revElements(gtf, minus)

  return(gtf)
}
