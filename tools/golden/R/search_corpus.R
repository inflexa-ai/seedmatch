source("tools/golden/R/common.R")
suppressPackageStartupMessages(library(Biostrings))

option_grid <- function(max_mms, indels_values) {
  option_sets <- list()
  for (max_mm in max_mms) {
    for (indels in indels_values) {
      for (min_mm in if (indels) 0L else 0:1) {
        for (fixed_pattern in c(TRUE, FALSE)) {
          for (fixed_subject in c(TRUE, FALSE)) {
            option_sets[[length(option_sets) + 1]] <- data.frame(
              max_mm = max_mm, min_mm = min_mm, indels = indels,
              fixed_pattern = fixed_pattern, fixed_subject = fixed_subject
            )
          }
        }
      }
    }
  }
  option_sets
}

search <- function(fun, case, opts) {
  fun(
    case$pattern, DNAString(case$subject),
    max.mismatch = opts$max_mm, min.mismatch = opts$min_mm, with.indels = opts$indels,
    fixed = c(pattern = opts$fixed_pattern, subject = opts$fixed_subject)
  )
}

describe <- function(opts) {
  paste(names(opts), vapply(opts, as.character, ""), sep = "=", collapse = " ")
}

write_corpus <- function(name, option_sets) {
  cases <- utils::read.delim(
    file.path(GOLDEN_DIR, "search", paste0(name, "_cases.tsv")),
    comment.char = "#", colClasses = "character", na.strings = character(0)
  )
  matches <- list()
  counts <- list()
  skipped <- character(0)
  for (opts in option_sets) {
    errors <- character(0)
    for (i in seq_len(nrow(cases))) {
      case <- cases[i, ]
      found <- tryCatch(
        list(views = search(matchPattern, case, opts), count = search(countPattern, case, opts)),
        error = function(e) conditionMessage(e)
      )
      if (is.character(found)) {
        errors[case$case_id] <- found
        next
      }
      counts[[length(counts) + 1]] <- cbind(data.frame(case_id = case$case_id), opts, count = found$count)
      if (length(found$views) > 0) {
        matches[[length(matches) + 1]] <- cbind(
          data.frame(case_id = case$case_id), opts,
          start = start(found$views), end = end(found$views), row.names = NULL
        )
      }
    }
    for (message in unique(errors)) {
      ids <- names(errors)[errors == message]
      who <- if (length(ids) == nrow(cases)) "each case" else paste("cases", paste(ids, collapse = ","))
      skipped <- c(skipped, paste0("skipped: ", describe(opts), " for ", who, ": ", message))
    }
  }
  for (table in c("matches", "counts")) {
    rows <- do.call(rbind, if (table == "matches") matches else counts)
    write_golden(rows, paste0("search/", name, "_", table, ".tsv"), "search_corpus.R", skipped)
  }
}

write_corpus("corpus", option_grid(0:4, c(FALSE, TRUE)))
write_corpus("seed_corpus", option_grid(1:4, TRUE))
