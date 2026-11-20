# A patched plotAlignment for the golden files. The version in SeedMatchR (R/plotAlignment.R:32)
# takes the guide target from get_seed(guide.seq, "Full", ...). That target starts at g2, thus it
# has one base less than the guide, and .matchStrings stops because the two strings differ in
# length (D-045). This version takes the guide target from .get_custom_sequence(guide.seq, 1,
# guide.len), thus the guide string has the full guide, g1 to gN. The other lines are the same as
# in the original.
plotAlignment <- function(match.info,
                          guide.seq,
                          seqs,
                          aln.name = "alignment",
                          left.pad = 0,
                          right.pad = 0,
                          seq.start.pos = 2,
                          seq.stop.pos = 8,
                          output.type = c("data.frame", "ggplot")){
  x <- y <- label <- NULL

  output.type = match.arg(output.type)

  feature.seq = seqs[names(seqs) %in% match.info$names]

  guide.len = nchar(guide.seq)

  guide.target = .get_custom_sequence(guide.seq, 1, guide.len)

  sirna.target.str = paste0(.make_space(left.pad),
                            as.character(guide.target),
                            .make_space(right.pad))

  sirna.str.base = as.character(Biostrings::complement(guide.target))

  d1 = (guide.len - seq.stop.pos)
  d2 = (seq.start.pos - 1)

  matching.seq = XVector::subseq(feature.seq,
                        match.info[1,"start"] - d1 - left.pad,
                        match.info[1,"end"] + d2 + right.pad)

  matching.str = as.character(matching.seq)

  match_char = .matchStrings(sirna.target.str,
                             matching.str)

  aln.df = .create_alignment_df(sirna.str.base,
                                match_char,
                                matching.str,
                                aln.name,
                                left.pad,
                                right.pad)

  if (output.type == "data.frame"){
    return(aln.df)

  } else {
    d2 <- paste0(paste(c("miRNA ","      ","target"), aln.df$alignment), collapse="\n")

    aln.plot <- ggplot2::ggplot(data.frame("x"=1, "y"=1, "label"=d2), ggplot2::aes(x,y,label=label)) +
      ggplot2::theme_void() + ggplot2::geom_text(family="mono", fontface="bold")

    return(aln.plot)

  }
}
