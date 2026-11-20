source("tools/golden/R/common.R")

N_BOOT <- 2000
ALTERNATIVES <- c("two.sided", "less", "greater")
VERSIONS <- c(
  paste0("twosamples: ", packageVersion("twosamples")),
  paste0("caTools: ", packageVersion("caTools"))
)

write_stats <- function(df, path, notes = character(0)) {
  # R gives NaN where a statistic is 0/0. The golden files write each missing value as NA.
  df[] <- lapply(df, function(col) if (is.double(col)) replace(col, is.nan(col), NA_real_) else col)
  write_golden(df, file.path("stats", path), "stats.R", c(VERSIONS, notes))
}

# One random stream feeds each case in the order of the script. A new case that draws moves each later
# case, thus a run changes the goldens of the later sections. Only a case without a draw keeps the other
# cases. A case at the end of the Wass section moves no other case.
set.seed(20261120)

cases <- list()
add_case <- function(test, id, x, y) {
  stopifnot(!(id %in% names(cases)))
  cases[[id]] <<- list(test = test, x = x, y = y)
}

# A case with ties rounds the values to 1 decimal, and at least one value occurs in both samples.
with_ties <- function(x, y) {
  x <- round(x, 1)
  y <- round(y, 1)
  if (!any(x %in% y)) {
    x[1] <- y[1]
  }
  list(x = x, y = y)
}

add_pair <- function(test, m, n, shift = 0.3, sd = 1) {
  x <- rnorm(m)
  y <- rnorm(n, mean = shift, sd = sd)
  stopifnot(anyDuplicated(c(x, y)) == 0)
  id <- paste0(test, "_", m, "x", n)
  add_case(test, paste0(id, "_plain"), x, y)
  tied <- with_ties(x, y)
  add_case(test, paste0(id, "_ties"), tied$x, tied$y)
}

scaled_statistic <- function(x, y) {
  sqrt(length(x) * length(y) / (length(x) + length(y))) * unname(stats::ks.test(x, y)$statistic)
}

# ---- KS ----
for (size in list(c(1, 12), c(5, 7), c(30, 40), c(99, 101), c(100, 100), c(120, 150), c(400, 300))) {
  add_pair("ks", size[1], size[2])
}
for (size in list(c(99, 101), c(100, 100), c(120, 150), c(400, 300))) {
  x <- rnorm(size[1])
  y <- rnorm(size[2], mean = 6)
  add_case("ks", paste0("ks_", size[1], "x", size[2], "_shift"), x, y)
}
add_case("ks", "ks_all_equal", rep(2.5, 5), rep(2.5, 7))
add_case("ks", "ks_inf", c(rnorm(5), Inf), rnorm(8))

# The grid of shifts has no random draw. The two samples are the quantiles of a normal distribution.
quantiles <- function(n, shift) qnorm((seq_len(n) - 0.5) / n) + shift
shifts <- seq(0, 1, by = 0.001)
scaled <- vapply(shifts, function(s) scaled_statistic(quantiles(120, 0), quantiles(150, s)), numeric(1))
shift <- shifts[which.min(abs(scaled - 0.95))]
stopifnot(scaled_statistic(quantiles(120, 0), quantiles(150, shift)) > 0.9)
stopifnot(scaled_statistic(quantiles(120, 0), quantiles(150, shift)) < 1)
add_case("ks", "ks_asym_cut_series", quantiles(120, 0), quantiles(150, shift))

add_case("ks", "ks_asym_small_p", rnorm(400), rnorm(300, mean = 1))

# ---- Wilcoxon ----
for (size in list(c(1, 5, 0.5), c(10, 12, -0.5), c(49, 49, 0.3), c(49, 50, -0.3), c(50, 49, 0.3), c(60, 80, 0.3))) {
  add_pair("wilcoxon", size[1], size[2], shift = size[3])
}
add_case("wilcoxon", "wilcoxon_all_equal", rep(1, 5), rep(1, 7))
add_case("wilcoxon", "wilcoxon_inf", c(rnorm(5), Inf), rnorm(8))

# ---- DTS ----
for (size in list(c(5, 7), c(30, 40), c(200, 300))) {
  add_pair("dts", size[1], size[2])
}

# ---- Wasserstein ----
for (size in list(c(10, 10, 0.5), c(10, 50, -0.4), c(40, 300, 0.3))) {
  add_pair("wass", size[1], size[2], shift = size[3], sd = 1.3)
}
add_case("wass", "wass_seed_one_value", rep(0.5, 10), rnorm(20))
add_case("wass", "wass_non_seed_one_value", rnorm(10), rep(0.5, 20))

# ---- samples.tsv ----
samples <- do.call(rbind, lapply(names(cases), function(id) {
  case <- cases[[id]]
  data.frame(
    case = id,
    sample = rep(c("x", "y"), c(length(case$x), length(case$y))),
    value = c(case$x, case$y),
    stringsAsFactors = FALSE
  )
}))
write_stats(samples, "samples.tsv", "the samples x and y of each case, in the order of the values")

ids_of <- function(test) names(Filter(function(case) case$test == test, cases))

# ---- ks.tsv ----
ks_rows <- list()
for (id in ids_of("ks")) {
  for (alternative in ALTERNATIVES) {
    scored <- suppressWarnings(
      stats::ks.test(cases[[id]]$x, cases[[id]]$y, alternative = alternative)
    )
    method <- if (startsWith(scored$method, "Exact")) {
      "exact"
    } else if (startsWith(scored$method, "Asymptotic")) {
      "asymptotic"
    } else {
      stop("ks.test method not known: ", scored$method)
    }
    ks_rows[[length(ks_rows) + 1]] <- data.frame(
      case = id, alternative = alternative, statistic = unname(scored$statistic),
      p_value = scored$p.value, method = method, stringsAsFactors = FALSE
    )
  }
}
ks <- do.call(rbind, ks_rows)

small_p <- ks[ks$case == "ks_asym_small_p" & ks$alternative == "greater", ]
en <- 400 * 300 / (400 + 300)
stopifnot(small_p$p_value == 0, exp(-2 * en * small_p$statistic^2) < 1e-17)
stopifnot(all(ks$method[ks$case %in% c("ks_99x101_plain", "ks_99x101_ties")] == "exact"))
stopifnot(all(ks$method[ks$case %in% c("ks_100x100_plain", "ks_100x100_ties")] == "asymptotic"))
write_stats(
  ks, "ks.tsv",
  "statistic: D, D^+ for greater, D^- for less. method: exact or asymptotic, from the method text of ks.test"
)

# ---- wilcoxon.tsv ----
wilcoxon_rows <- list()
for (id in ids_of("wilcoxon")) {
  for (alternative in ALTERNATIVES) {
    scored <- suppressWarnings(
      stats::wilcox.test(cases[[id]]$x, cases[[id]]$y, alternative = alternative)
    )
    method <- if (grepl("continuity correction", scored$method, fixed = TRUE)) {
      "asymptotic"
    } else if (grepl("exact", scored$method, fixed = TRUE)) {
      "exact"
    } else {
      stop("wilcox.test method not known: ", scored$method)
    }
    wilcoxon_rows[[length(wilcoxon_rows) + 1]] <- data.frame(
      case = id, alternative = alternative, statistic = unname(scored$statistic),
      p_value = scored$p.value, method = method, stringsAsFactors = FALSE
    )
  }
}
wilcoxon <- do.call(rbind, wilcoxon_rows)
write_stats(
  wilcoxon, "wilcoxon.tsv",
  c(
    "statistic: W. method: exact, or asymptotic for the normal approximation with the continuity correction",
    "NA in p_value: R gives NaN, because each value is equal and z is 0/0"
  )
)

# ---- dts.tsv ----
dts_rows <- lapply(ids_of("dts"), function(id) {
  set.seed(1)
  scored <- twosamples::dts_test(cases[[id]]$x, cases[[id]]$y, nboots = N_BOOT)
  data.frame(
    case = id, statistic = unname(scored[[1]]), p_value = unname(scored[[2]]), n_boot = N_BOOT,
    stringsAsFactors = FALSE
  )
})
write_stats(do.call(rbind, dts_rows), "dts.tsv", "set.seed(1) comes before each dts_test")

# ---- wasserstein.tsv ----
wass_rows <- list()
wass_notes <- character(0)
for (id in ids_of("wass")) {
  for (n_interp in c(2, 100, 10000)) {
    areas <- tryCatch(
      suppressWarnings(SeedMatchR::wass_dist(cases[[id]]$x, cases[[id]]$y, n.interp = n_interp)),
      error = function(e) conditionMessage(e)
    )
    if (is.character(areas)) {
      if (n_interp == 2) {
        wass_notes <- c(wass_notes, paste0("R stops for ", id, ": ", areas))
      }
      next
    }
    wass_rows[[length(wass_rows) + 1]] <- data.frame(
      case = id, n_interp = n_interp, abs_auc = areas$abs.auc, neg_cens_auc = areas$neg.cens.auc,
      pos_cens_auc = areas$pos.cens.auc, stringsAsFactors = FALSE
    )
  }
}
stopifnot(length(wass_notes) == 2)
write_stats(
  do.call(rbind, wass_rows), "wasserstein.tsv",
  c(
    "case x is the seed sample and case y is the non-seed sample. R rounds each area to 3 decimals",
    wass_notes,
    "the cases that stop have no row"
  )
)

# ---- round.tsv ----
boundary <- c(0.0005, 0.0015, 0.0025, 0.0035, 0.1235, 2.675, 1.0005, -0.0025, -0.1235)
grid <- (0:4000) / 4000 + 0.0005
# The random values of round.tsv have their own seed. Thus no case moves them.
set.seed(20261120)
random <- runif(200, 0, 3)
rounded <- function(x, digits) data.frame(x = x, digits = digits, r = round(x, digits))
write_stats(
  rbind(
    rounded(boundary, 3), rounded(grid, 3), rounded(random, 3),
    rounded(boundary, 2), rounded(boundary, 0)
  ),
  "round.tsv", "r is round(x, digits) of R 4.5.3"
)

# ---- ecdf.tsv ----
data <- load_test_data()

ecdf_stats <- c(KS = "KS", Wilcoxon = "Wilcoxen", Wass = "Wass", DTS = "DTS")
ecdf_rows <- function(input, res, target, background, scope) {
  fold_changes <- function(genes) res$log2FoldChange[res$gene_id %in% genes]
  counts <- c(sum(!is.na(fold_changes(target))), sum(!is.na(fold_changes(background))))
  rows <- list()
  for (test in names(ecdf_stats)) {
    # R ignores the alternative for Wass and DTS, thus their rows have NA.
    alternatives <- if (test %in% c("KS", "Wilcoxon")) ALTERNATIVES else NA_character_
    for (alternative in alternatives) {
      key <- if (is.na(alternative)) test else paste(test, alternative)
      if (!is.null(scope) && !(key %in% scope)) {
        next
      }
      args <- list(res, target, background, stats_test = ecdf_stats[[test]])
      if (test %in% c("KS", "Wilcoxon")) {
        args$alternative <- alternative
      }
      if (test == "DTS") {
        set.seed(1)
      }
      scored <- suppressWarnings(do.call(SeedMatchR::ecdf_stat_test, args))
      rows[[length(rows) + 1]] <- data.frame(
        input = input, test = test, alternative = alternative,
        statistic = unname(scored$statistic),
        p_value = if (is.null(scored$p.value)) NA_real_ else unname(scored$p.value),
        n_target = counts[1], n_background = counts[2], stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

res <- data$res
target <- res$gene_id[res$mer7m8 >= 1]
background <- res$gene_id[res$mer7m8 == 0]
test_data <- ecdf_rows("test_data", res, target, background, NULL)
swapped <- ecdf_rows(
  "test_data_swapped", res, background, target, c("Wilcoxon less", "Wass")
)

readme_res <- SeedMatchR::SeedMatchR(
  seqs = data$seqs, sequence = data$guide,
  res = SeedMatchR::filter_res(res[, colnames(res) != "mer7m8"]), seed.name = "mer7m8"
)
readme <- ecdf_rows(
  "readme", readme_res, readme_res$gene_id[readme_res$mer7m8 >= 1],
  readme_res$gene_id[readme_res$mer7m8 == 0], NULL
)

ks_greater <- test_data[test_data$test == "KS" & test_data$alternative == "greater", ]
stopifnot(abs(ks_greater$statistic - 0.51035) < 5e-6, abs(ks_greater$p_value - 0.0001037731) < 5e-10)
wilcoxon_less <- swapped[swapped$test == "Wilcoxon", ]
stopifnot(wilcoxon_less$statistic == 2150, abs(wilcoxon_less$p_value - 0.9996) < 5e-5)
stopifnot(swapped$statistic[swapped$test == "Wass"] == 0)

write_stats(
  rbind(test_data, swapped, readme), "ecdf.tsv",
  c(
    "input: test_data is DESEQ2_results.tsv, with the target mer7m8 >= 1 and the background mer7m8 == 0",
    "input: test_data_swapped has the target and the background of test_data in the other order",
    "input: readme is filter_res with the defaults, then the mer7m8 counts of SeedMatchR, with the same sets",
    "test: Wilcoxon is the stats_test Wilcoxen of R. Wass and DTS ignore the alternative in R, thus their alternative is NA",
    "statistic of Wass is neg.cens.auc, and its p_value is NA. n_target and n_background count the fold changes that are not NA"
  )
)
