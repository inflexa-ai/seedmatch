# A patched vmatchPattern2 for the golden files. The version in SeedMatchR (R/vmatchPattern2.R:43)
# sends the algorithm "indels" to the C matcher for each call, so its ranges ignore with.indels and
# min.mismatch (docs/divergences.md, D-001). This version gives the ranges that agree with the
# counts of vcountPattern. The data frame has the shape of the original: the columns start, end,
# width, and names, with one row for each match and no row for a sequence without a match.
vmatchPattern2 <- function(pattern,
                           subject,
                           max.mismatch = 0,
                           min.mismatch = 0,
                           with.indels = FALSE,
                           fixed = TRUE,
                           algorithm = "auto") {
  pattern <- Biostrings::DNAString(pattern)

  if (with.indels) {
    # vmatchPattern does not support indels for a DNAStringSet, so each sequence has its own search.
    views <- lapply(seq_along(subject), function(i) {
      Biostrings::matchPattern(pattern, subject[[i]],
        max.mismatch = max.mismatch, min.mismatch = min.mismatch,
        with.indels = TRUE, fixed = fixed
      )
    })
    starts <- lapply(views, function(v) as.integer(BiocGenerics::start(v)))
    widths <- lapply(views, function(v) as.integer(BiocGenerics::width(v)))
  } else {
    index <- Biostrings::vmatchPattern(pattern, subject,
      max.mismatch = max.mismatch, min.mismatch = min.mismatch,
      with.indels = FALSE, fixed = fixed
    )
    starts <- lapply(Biostrings::startIndex(index), as.integer)
    widths <- lapply(starts, function(s) rep(length(pattern), length(s)))
  }

  unlisted_ans <- IRanges::IRanges(
    start = unlist(starts, use.names = FALSE),
    width = unlist(widths, use.names = FALSE)
  )

  ir.list <- IRanges::relist(unlisted_ans, starts)

  names(ir.list) <- names(subject)

  ir.list <- ir.list[!sapply(ir.list, function(x) length(x) == 0)]

  as.data.frame(unlist(ir.list))
}
