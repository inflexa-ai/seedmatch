# A patched get_compensatory_sequences for the golden files. The version in SeedMatchR
# (R/get_compensatory_sequences.R:26-133) always returns an empty data frame (D-004). Its test of
# nrow(compensatory_regions_captured) == 0 at :67 reads the frame of :57-66, and this frame is
# always empty. This version removes the test, its early return, the else, and the brace that
# closes the else (:67-77, :130). It keeps the body of the else (:78-129). After this change, :99
# stops for the first match, because gsub at :34 makes the seed a character string and
# reverseComplement has no method for it. Thus this version gives :99 the seed as a DNAString. The
# other lines are the same as in the original.
get_compensatory_sequences = function(guide.seq,
                                      seed,
                                      sequences,
                                      compensatory_region_length_min = 6,
                                      compensatory_region_length_max = NA,
                                      limit=100,
                                      mismatch_ = 1,
                                      mismatch_compensatory_region_ = 0) {
  seed = gsub("U", "T", seed)

  non_can_seed_match_transcripts <- names(sequences)[which(Biostrings::vcountPattern(seed,
                                                                                     sequences,
                                                                                     max.mismatch = mismatch_,
                                                                                     with.indel = FALSE) > 0)]

  non_can_seed_match_transcripts <- sequences[non_can_seed_match_transcripts]

  stopifnot(is.character(guide.seq), nchar(guide.seq) >= 8)
  comp_seq <- substr(guide.seq, 8, nchar(guide.seq))

  if(is.na(compensatory_region_length_max)){
    compensatory_region_length_max <- nchar(comp_seq)
  }

  # All possible substrings of indexes of the compensatory region
  comp_seqs <- c()
  for (comp_len in compensatory_region_length_min:compensatory_region_length_max){
    comp_seqs <- c(comp_seqs, sapply(1:(nchar(comp_seq) - comp_len + 1), function(i){ return(gsub(" ","",substr(comp_seq,i,i+comp_len-1)))}))
  }

  # Data Frame Initialization
  compensatory_regions_captured <- data.frame("Seed" = character(),
                                             "Compensatory_Region" = character(),
                                             "TX_ID" = character(),
                                             "Distance" = numeric(),
                                             "Seed_Region_in_Guide" = character(),
                                             "Seed_Start" = numeric(),
                                             "Seed_End" = numeric(),
                                             "Compensatory_Start" = numeric(),
                                             "Compensatory_End" = numeric(),
                                             "Compensatory_Region_in_Guide" = character())
  seed_matches = as.data.frame(unlist(Biostrings::vmatchPattern(as.character(seed), non_can_seed_match_transcripts, max.mismatch = mismatch_)))
  seed_matches <- stats::setNames(seed_matches, c("seed_start","seed_end","seed_width", "gene_name"))

  # LOOP OVER COMPENSATORY SEQS
  k <- 1
  for (seq in comp_seqs){

    rev_seq <- Biostrings::reverseComplement(Biostrings::DNAString(gsub("U", "T", seq)))
    tmp = Biostrings::vmatchPattern(as.character(rev_seq), non_can_seed_match_transcripts, max.mismatch = mismatch_compensatory_region_)  # subset of non_can_seed_matches matching kmer
    rev_seq_matches <- as.data.frame(unlist(tmp))
    rev_seq_matches <- stats::setNames(rev_seq_matches, c("start","end","width", "gene_name"))

    merged_matches <- merge(seed_matches, rev_seq_matches, by = "gene_name")
    # Reduced merges
    merged_matches <- merged_matches[ merged_matches$seed_start > merged_matches$end + 1 & merged_matches$seed_start < merged_matches$end + limit,]

    if (nrow(merged_matches) > 0){
      merged_matches$Seed_in_Transcript = as.character(seed)
      merged_matches$Compensatory_Region_in_Transcript = as.character(rev_seq)
      merged_matches$TX_ID = merged_matches$gene_name
      merged_matches$Distance = merged_matches$seed_start - merged_matches$end - 1
      merged_matches$Seed_Region_in_Guide = as.character(Biostrings::reverseComplement(Biostrings::DNAString(seed)))
      merged_matches$Seed_Start = merged_matches$seed_end
      merged_matches$Seed_End = merged_matches$seed_start
      merged_matches$Compensatory_Start = merged_matches$end
      merged_matches$Compensatory_End = merged_matches$start
      merged_matches$Compensatory_Region_in_Guide = as.character(Biostrings::reverseComplement(rev_seq))

      merged_matches$mismatch_in_seed = as.numeric(as.character(seed) != as.character(Biostrings::subseq(non_can_seed_match_transcripts[merged_matches$gene_name],
                                                                                                    start=merged_matches$seed_start,
                                                                                                    end=pmin(Biostrings::width(non_can_seed_match_transcripts[merged_matches$gene_name]), merged_matches$seed_end))))
      compensatory_regions_captured = rbind(compensatory_regions_captured, merged_matches)
    }

    k <- k + 1
  }

  # Clean up extra columns
  compensatory_regions_captured <- compensatory_regions_captured[c("Seed_in_Transcript", "Compensatory_Region_in_Transcript", "TX_ID", "Distance", "Seed_Region_in_Guide", "Seed_Start",
                                    "Seed_End", "Compensatory_Start", "Compensatory_End", "Compensatory_Region_in_Guide","mismatch_in_seed")]


    # DATA SUMMARY
    compensatory_regions_captured$GC_Content = sapply(compensatory_regions_captured$Compensatory_Region_in_Transcript, .gc_content)
    compensatory_regions_captured$CompensatoryLength = sapply(compensatory_regions_captured$Compensatory_Region_in_Transcript, nchar)
    compensatory_regions_captured = compensatory_regions_captured %>% dplyr::group_by(TX_ID, Compensatory_Start) %>%  dplyr::filter(CompensatoryLength == max(CompensatoryLength)) %>% dplyr::ungroup()
    compensatory_regions_captured$check_d_tl = compensatory_regions_captured$Compensatory_Start - compensatory_regions_captured$CompensatoryLength
    compensatory_regions_captured = compensatory_regions_captured %>% dplyr::group_by(TX_ID, check_d_tl) %>%  dplyr::filter(CompensatoryLength == max(CompensatoryLength)) %>% dplyr::ungroup()
    compensatory_regions_captured$check_d_tl = NULL
    compensatory_regions_captured$relative_position_guide_seq_start = sapply(compensatory_regions_captured$Compensatory_Region_in_Guide,function(x){return(relative_position_seq(x,guide.seq)$start)})
    compensatory_regions_captured$relative_position_guide_seq_end = sapply(compensatory_regions_captured$Compensatory_Region_in_Guide,function(x){return(relative_position_seq(x,guide.seq)$end)})
    return(compensatory_regions_captured)


}
