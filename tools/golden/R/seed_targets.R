source("tools/golden/R/common.R")

guides <- c("UUAUAGAGCAAGAACACUGUUUU", "UGAGGUAGUAGGUUGUAUAGUU", "AUUUUUUUUUUUUUUUUUUUA")
seeds <- c("mer7m8", "mer8", "mer6", "mer7A1", "Full", "18mer", "15mer")
windows <- list(c(2, 19), c(1, 23))

row_of <- function(guide, seed, start, stop, wobble, result) {
  data.frame(
    guide = guide,
    seed = seed,
    start = if (is.na(start)) "" else as.character(start),
    stop = if (is.na(stop)) "" else as.character(stop),
    wobble = wobble,
    target = as.character(result[["Target.Seq"]]),
    seed_rna = as.character(result[["Seed.Seq.RNA"]]),
    stringsAsFactors = FALSE
  )
}

rows <- list()
for (guide in guides) {
  for (seed in seeds) {
    for (wobble in c(FALSE, TRUE)) {
      result <- SeedMatchR::get_seed(guide, seed.name = seed, allow_wobbles = wobble)
      rows[[length(rows) + 1]] <- row_of(guide, seed, NA, NA, wobble, result)
    }
  }
}
for (window in windows) {
  for (wobble in c(FALSE, TRUE)) {
    result <- SeedMatchR::get_seed(
      guides[1], seed.name = "custom", start.pos = window[1], stop.pos = window[2], allow_wobbles = wobble
    )
    rows[[length(rows) + 1]] <- row_of(guides[1], "custom", window[1], window[2], wobble, result)
  }
}

write_golden(do.call(rbind, rows), "seed/targets.tsv", "seed_targets.R")
