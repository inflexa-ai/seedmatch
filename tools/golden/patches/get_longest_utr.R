# A patched .get_longest_utr for the golden files. The version in SeedMatchR
# (R/load_annotations.R:406-418) selects the longest transcript of each gene from transcriptLengths,
# which also holds the transcripts that are not in utr.grlist. If a gene has no coding transcript,
# the selected transcript is not in utr.grlist, and the function stops with "subscript contains
# invalid names" (D-028). This version keeps only the rows of tx_lens whose tx_id is in
# names(utr.grlist). The other lines are the same as in the original.
.get_longest_utr <- function(utr.grlist, txdb, filter = NULL){
  # Extract transcript lengths for each gene
  tx_lens <- GenomicFeatures::transcriptLengths(txdb, with.utr3_len=TRUE, filter = filter)

  tx_lens <- tx_lens[tx_lens$tx_id %in% names(utr.grlist), ]

  # Find the longest transcript per gene
  longest_utr_per_gene = tx_lens %>% dplyr::group_by(gene_id) %>% dplyr::summarize(longest_UTR = max(utr3_len),
                                                                            tx_id = tx_id[which.max(utr3_len)])

  # Subset annotations for the longest transcript
  utr.grlist = utr.grlist[longest_utr_per_gene$tx_id]

  return(utr.grlist)
}
