source("tools/golden/R/common.R")

PATCH <- "tools/golden/patches/get_compensatory_sequences.R"

# R quotes a name in an error message with curly quotes in a UTF-8 locale. The golden files hold ASCII only.
options(useFancyQuotes = FALSE)

write_noncanonical <- function(df, path, notes = character(0)) {
  stop_on_nan(df)
  write_golden(df, file.path("noncanonical", path), "noncanonical.R", notes)
}

R_COLUMNS <- c(
  "Seed_in_Transcript", "Compensatory_Region_in_Transcript", "TX_ID", "Distance", "Seed_Region_in_Guide",
  "Seed_Start", "Seed_End", "Compensatory_Start", "Compensatory_End", "Compensatory_Region_in_Guide",
  "mismatch_in_seed", "GC_Content", "CompensatoryLength", "relative_position_guide_seq_start",
  "relative_position_guide_seq_end"
)

data <- load_test_data()
seqs <- data$seqs
seed_target <- function(guide, seed) SeedMatchR::get_seed(guide, seed)$Target.Seq

# ---- Compensatory search ----
patched_compensatory <- function(guide, seed, sequences, ...) {
  with_seedmatchr_patch(PATCH, suppressMessages(
    SeedMatchR::get_compensatory_sequences(guide, seed_target(guide, seed), sequences, ...)
  ))
}

test_data_cases <- list(
  vignette = list(seed = "mer7m8"),
  seed_mm0 = list(seed = "mer7m8", mismatch_ = 0),
  min8_max10 = list(
    seed = "mer7m8", compensatory_region_length_min = 8, compensatory_region_length_max = 10
  ),
  limit40 = list(seed = "mer7m8", limit = 40),
  region_mm1 = list(
    seed = "mer7m8", compensatory_region_length_min = 9, mismatch_compensatory_region_ = 1
  ),
  mer6 = list(seed = "mer6")
)

rows <- list()
for (id in names(test_data_cases)) {
  arguments <- test_data_cases[[id]]
  result <- do.call(
    patched_compensatory, c(list(GUIDE, arguments$seed, seqs), arguments[names(arguments) != "seed"])
  )
  rows[[id]] <- cbind(case = id, as.data.frame(result)[R_COLUMNS])
}

cases <- utils::read.delim(
  file.path(DATA_DIR, "compensatory_cases.tsv"), colClasses = "character", na.strings = character(0)
)
notes <- character(0)
for (id in unique(cases$case)) {
  part <- cases[cases$case == id, ]
  stopifnot(all(vapply(
    part[c("guide", "min_length", "max_length", "limit", "seed_mismatch")],
    function(col) length(unique(col)) == 1, logical(1)
  )))
  one <- part[1, ]
  case_seqs <- Biostrings::DNAStringSet(setNames(part$sequence, part$seq_id))
  result <- tryCatch(
    patched_compensatory(
      one$guide, "mer7m8", case_seqs,
      compensatory_region_length_min = as.integer(one$min_length),
      compensatory_region_length_max = if (one$max_length == "NA") NA else as.integer(one$max_length),
      limit = as.integer(one$limit), mismatch_ = as.integer(one$seed_mismatch)
    ),
    error = root_message
  )
  if (is.character(result)) {
    notes <- c(notes, paste0(id, ": the patched function stops with: ", result))
  } else {
    rows[[id]] <- cbind(case = id, as.data.frame(result)[R_COLUMNS])
  }
}

# The unpatched function for the vignette case.
unpatched_messages <- character(0)
unpatched <- withCallingHandlers(
  SeedMatchR::get_compensatory_sequences(GUIDE, seed_target(GUIDE, "mer7m8"), seqs),
  message = function(m) {
    unpatched_messages <<- c(unpatched_messages, trimws(conditionMessage(m)))
    invokeRestart("muffleMessage")
  }
)
stopifnot(nrow(unpatched) == 0, identical(unpatched_messages, "No compensatory regions found."))

line_99_error <- error_of(
  Biostrings::reverseComplement(gsub("U", "T", seed_target(GUIDE, "mer7m8")))
)

repeated_sequence <- cases$sequence[cases$case == "repeated_name_r"]
repeated_set <- Biostrings::DNAStringSet(setNames(repeated_sequence, c("g", "g")))
repeated_selected <- sum(as.character(repeated_set[c("g", "g")]) == repeated_sequence[1])
stopifnot(repeated_selected == 2)

region_at_g2_guide <- cases$guide[cases$case == "region_at_g2_r"][1]
stopifnot(identical(region_at_g2_guide, "UUAUAGAUAUAGACACUGUUUUA"))
relative <- SeedMatchR:::relative_position_seq("TATAGA", region_at_g2_guide)
stopifnot(identical(relative, list(start = "2,8", end = "7,13")))

# The loop of R/get_compensatory_sequences.R:52-54 for a minimum above the length of g8 to gN. The
# maximum is that length, thus the loop runs min:max = 17:16. The length 17 gives the starts 1 and 0,
# and the length 16 gives the start 1.
region <- substr(GUIDE, 8, nchar(GUIDE))
region_min <- 17
region_substrings <- unlist(lapply(region_min:nchar(region), function(length) {
  vapply(
    1:(nchar(region) - length + 1), function(i) substr(region, i, i + length - 1), character(1)
  )
}))
stopifnot(identical(region_substrings, rep(region, 3)))

gc_content <- c(SeedMatchR:::.gc_content("GCAT"), SeedMatchR:::.gc_content(""))
stopifnot(gc_content[1] == 50, is.nan(gc_content[2]))

write_noncanonical(
  do.call(rbind, rows), "compensatory.tsv",
  c(
    paste0("dplyr: ", packageVersion("dplyr")),
    "patch: tools/golden/patches/get_compensatory_sequences.R (D-004)",
    notes,
    paste0(
      "original_vignette: the unpatched function gives ", nrow(unpatched), " rows and the message: ",
      unpatched_messages
    ),
    paste0("line_99: reverseComplement of the seed after gsub stops with: ", line_99_error),
    paste0(
      "repeated_name_r: the set of the names with a seed match, c(\"g\", \"g\"), selects the first sequence ",
      repeated_selected, " times"
    ),
    paste0(
      "relative_position_seq_g6: relative_position_seq(\"TATAGA\", \"", region_at_g2_guide,
      "\"), with the guide of region_at_g2_r, gives start ", relative$start, " and end ", relative$end
    ),
    paste0(
      "min_17: the substrings of comp_seqs for min 17 are ", paste(region_substrings, collapse = ",")
    ),
    paste0(
      "gc_content: .gc_content(\"GCAT\") = ", gc_content[1], ", .gc_content(\"\") = ", gc_content[2]
    )
  )
)

# ---- Wobble and bulge pipeline ----
seed_combinations <- function(guide) {
  seeds <- SeedMatchR:::.get_all_combination_w_b(guide, silent = TRUE)$seed_combinations
  unname(vapply(seeds, as.character, character(1)))
}

# The pipeline prints a line also with silent = TRUE (R/run_wobble_bulges_pipeline.R:22).
pipeline_error <- function(...) {
  error_of(utils::capture.output(
    SeedMatchR::run_wobble_bulges_pipeline(GUIDE, res = data$res, features = seqs, silent = TRUE, ...)
  ))
}
bulge_allowance_error <- "object 'bulge_allowance' not found"
pipeline_errors <- c(pipeline_error(), pipeline_error(bulges_lim = 1))
stopifnot(pipeline_errors == bulge_allowance_error)

seeds_5_u <- seed_combinations("AUUUUUAGCAAGAACACAGAUAA")
bulges_nt <- SeedMatchR::check_bulges_nt(Biostrings::RNAStringSet(c("CNUCUAUA", "CUCUANUA", "CUCNUAUA")), 3)
stopifnot(as.character(bulges_nt) == c("CUCUANUA", "CUCNUAUA"))

wobble_notes <- c(
  paste0("pipeline_error: run_wobble_bulges_pipeline(GUIDE, res, seqs) stops with: ", pipeline_errors[1]),
  paste0("pipeline_error_bulges_lim_1: with bulges_lim = 1, it stops with: ", pipeline_errors[2]),
  paste0("wobble_seeds_guide: ", paste(seed_combinations(GUIDE), collapse = ",")),
  paste0(
    "wobble_seeds_4_u: AUAUAGAGCAAGAACACAGAUAU (4 U) gives only ",
    paste(seed_combinations("AUAUAGAGCAAGAACACAGAUAU"), collapse = ",")
  ),
  paste0(
    "wobble_seeds_2_u: AUAGAGAGCAAGAACACAGAUAA (2 U) stops with: ",
    error_of(seed_combinations("AUAGAGAGCAAGAACACAGAUAA"))
  ),
  paste0(
    "wobble_seeds_0_u: CAGAGAGCAAGAACACAGACAAC (no U) stops with: ",
    error_of(suppressWarnings(seed_combinations("CAGAGAGCAAGAACACAGACAAC")))
  ),
  paste0("wobble_seeds_5_u_in_seed: AUUUUUAGCAAGAACACAGAUAA gives ", paste(seeds_5_u, collapse = ",")),
  paste0(
    "formals_get_target_bulged_guides: ", paste(names(formals(SeedMatchR:::.get_target_bulged_guides)), collapse = ","),
    "; the pipeline passes (x, bulges_lim, guide_seq_length)"
  ),
  paste0(
    "check_bulges_nt: of CNUCUAUA, CUCUANUA, and CUCNUAUA, nt = 3 keeps ",
    paste(as.character(bulges_nt), collapse = ",")
  )
)

# R stops before it counts, because .sequence_with_tb reads bulge_allowance from the global
# environment and the pipeline gives .get_target_bulged_guides its arguments in the wrong order.
# The script sets the variable for the time of one call and gives the arguments in the documented
# order, thus it shows the count rule of R/run_wobble_bulges_pipeline.R:35-55.
with_bulge_allowance <- function(allowance, code) {
  assign("bulge_allowance", allowance, envir = globalenv())
  on.exit(rm("bulge_allowance", envir = globalenv()))
  code
}

bulged_strings <- function(seed_rna, guide_length, bulges_lim) {
  bases <- strsplit(seed_rna, "")[[1]]
  strings <- unlist(lapply(seq_len(bulges_lim), function(allowance) {
    with_bulge_allowance(allowance, vapply(
      allowance:(length(bases) - allowance - 1),
      function(start) SeedMatchR:::.sequence_with_tb(start, bases, guide_length),
      character(1)
    ))
  }))
  Biostrings::RNAStringSet(strings)
}

r_rule_count <- function(guide, sequence, bulges_lim, nt = 3) {
  subject <- Biostrings::DNAStringSet(sequence)
  seeds <- gsub("T", "U", seed_combinations(guide))
  sum(vapply(seeds, function(seed) {
    strings <- SeedMatchR::check_bulges_nt(bulged_strings(seed, nchar(guide), bulges_lim), nt)
    sum(vapply(as.character(strings), function(string) {
      Biostrings::vcountPattern(Biostrings::DNAString(gsub("U", "T", string)), subject, fixed = FALSE)
    }, integer(1)))
  }, integer(1)))
}

wobble_cases <- utils::read.delim(
  file.path(DATA_DIR, "wobble_bulge_cases.tsv"), colClasses = "character", na.strings = character(0)
)
for (limit in c(1, 4)) {
  wobble_cases[[paste0("r_rule_bulges_lim_", limit)]] <- unname(mapply(
    r_rule_count, wobble_cases$guide, wobble_cases$sequence, MoreArgs = list(bulges_lim = limit)
  ))
}
stopifnot(!exists("bulge_allowance", envir = globalenv()))

write_noncanonical(
  wobble_cases, "wobble_bulge_r.tsv",
  c(paste0("gtools: ", packageVersion("gtools")), wobble_notes)
)
