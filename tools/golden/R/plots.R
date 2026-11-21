source("tools/golden/R/common.R")

# ggplot_build and hist open a device. Without this one, R writes Rplots.pdf into the working directory.
grDevices::pdf(NULL)
suppressPackageStartupMessages(library(ggplot2))

PATCH_DIR <- "tools/golden/patches"
VERSIONS <- c(
  paste0("ggplot2: ", packageVersion("ggplot2")),
  paste0("gt: ", packageVersion("gt")),
  paste0("scales: ", packageVersion("scales")),
  paste0("PRROC: ", packageVersion("PRROC"))
)

write_plot <- function(df, path, notes = character(0)) {
  stop_on_nan(df)
  write_golden(df, file.path("plots", path), "plots.R", c(VERSIONS, notes))
}

note <- function(key, value) {
  stopifnot(length(value) == 1)
  paste0(key, ": ", gsub("\n", " | ", value, fixed = TRUE))
}
list_of <- function(x) paste(x, collapse = ",")
fields <- function(...) paste(c(...), collapse = " / ")
or_none <- function(x) if (is.null(x)) "none" else paste(x, collapse = ",")

data <- load_test_data()
res <- data$res
seqs <- data$seqs
seed_genes <- res$gene_id[res$mer7m8 >= 1]
background <- res$gene_id[res$mer7m8 == 0]
base_lists <- list(Background = background, mer7m8 = seed_genes)

# The graphics functions that the plots call. A tracer runs in the frame of the traced function and
# records its arguments in traced_calls, which each tracer finds through the global environment.
traced_calls <- new.env()
reset_traced <- function() {
  for (name in c("title", "arrows", "segments", "image", "axis")) traced_calls[[name]] <- list()
}
reset_traced()
TRACERS <- list(
  title = quote(traced_calls$title[[length(traced_calls$title) + 1]] <- list(main = main, xlab = xlab, ylab = ylab)),
  arrows = quote(traced_calls$arrows[[length(traced_calls$arrows) + 1]] <- list(x0 = x0, y0 = y0, x1 = x1, col = col)),
  segments = quote(traced_calls$segments[[length(traced_calls$segments) + 1]] <- list(
    x0 = x0, y0 = y0, x1 = x1, y1 = y1, col = col, lwd = lwd, usr = par("usr")
  )),
  image.default = quote(traced_calls$image[[length(traced_calls$image) + 1]] <- list(x = x, y = y, col = col)),
  axis = quote(traced_calls$axis[[length(traced_calls$axis) + 1]] <- list(side = side))
)

with_tracers <- function(names, code) {
  graphics_ns <- asNamespace("graphics")
  for (name in names) suppressMessages(trace(name, tracer = TRACERS[[name]], where = graphics_ns, print = FALSE))
  on.exit(for (name in names) suppressMessages(untrace(name, where = graphics_ns)))
  code
}

with_warnings <- function(expr) {
  seen <- character(0)
  value <- withCallingHandlers(expr, warning = function(w) {
    seen <<- c(seen, one_line(conditionMessage(w)))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = unique(seen))
}

title_text <- function(call) {
  show <- function(x) if (is.null(x)) "NULL" else x
  paste0("main=", show(call$main), ", xlab=", show(call$xlab), ", ylab=", show(call$ylab))
}

# ---- ECDF ----
legend_keys <- function(pl) as.character(suppressWarnings(get_guide_data(pl, "colour"))$.value)
legend_of <- function(pl) {
  guide <- suppressWarnings(get_guide_data(pl, "colour"))
  paste(paste0(guide$.label, " [", guide$colour, "]"), collapse = " / ")
}
ecdf_plot <- function(res, lists, ...) suppressWarnings(SeedMatchR::deseq_fc_ecdf(res, lists, ...))$plot

# The rows of the layer in the order of a step line: by group, then by x.
ecdf_rows <- function(case, pl) {
  layer <- ggplot_build(pl)$data[[1]]
  keys <- legend_keys(pl)
  stopifnot(length(keys) == max(layer$group))
  layer <- layer[order(layer$group, layer$x), ]
  data.frame(case = case, group = keys[layer$group], x = layer$x, y = layer$y, stringsAsFactors = FALSE)
}

test_data_plot <- ecdf_plot(res, base_lists)
test_data_build <- ggplot_build(test_data_plot)

no_group_res <- rbind(res, res[1, ])
no_group_res$gene_id[nrow(no_group_res)] <- "ENSNOGROUP"
no_group_res$log2FoldChange[nrow(no_group_res)] <- 5
no_group_plot <- ecdf_plot(no_group_res, base_lists)

reversed_plot <- ecdf_plot(res, base_lists, factor_order = c("mer7m8", "Background"))
repeated_plot <- ecdf_plot(res, list(Background = c(background, "ENSABSENT1", "ENSABSENT2", background[1]), mer7m8 = seed_genes))

na_res <- res
na_res$log2FoldChange[na_res$gene_id == seed_genes[1]] <- NA
na_plot <- ecdf_plot(na_res, base_lists)
na_built <- with_warnings(ggplot_build(na_plot))
na_layer <- na_built$value$data[[1]]

big <- data.frame(gene_id = sprintf("G%05d", 1:12345), log2FoldChange = seq(-1, 1, length.out = 12345))
big_plot <- ecdf_plot(big, list(Background = big$gene_id[1:12000], target = big$gene_id[12001:12345]))

subtitle_of <- function(stats_test) ecdf_plot(res, base_lists, stats_test = stats_test)$labels$subtitle
set.seed(1)
dts_subtitle <- subtitle_of("DTS")

bottom_plot <- ecdf_plot(res, base_lists, theme = ggplot2::theme(legend.position = "bottom"))

lists10 <- split(res$gene_id, rep(1:10, length.out = nrow(res)))
names(lists10) <- paste0("g", 1:10)
ten_groups_plot <- ecdf_plot(res, lists10)
ten_groups_error <- error_of(ggplot_build(ten_groups_plot))

png_name <- file.path(tempdir(), "ecdf.png")
invisible(SeedMatchR::deseq_fc_ecdf(res, base_lists, output_filename = png_name))
png_magic <- rawToChar(readBin(png_name, "raw", 5))
unlink(png_name)

layer_of_plot <- test_data_plot$layers[[1]]
panel <- test_data_build$layout$panel_params[[1]]
plot_theme <- SeedMatchR:::SeedMatchR.theme
theme_note <- fields(
  paste0("text colour ", plot_theme$text$colour, " size ", plot_theme$text$size),
  paste0("title size ", plot_theme$title$size),
  paste0(
    "panel.grid.major colour ", plot_theme$panel.grid.major$colour, " linetype ", plot_theme$panel.grid.major$linetype,
    " linewidth ", plot_theme$panel.grid.major$linewidth
  ),
  paste0("panel.background fill ", plot_theme$panel.background$fill, " colour ", plot_theme$panel.background$colour),
  paste0("legend.box.background fill ", plot_theme$legend.box.background$fill, " colour ", plot_theme$legend.box.background$colour),
  paste0(
    "legend.position ", plot_theme$legend.position, " at ", list_of(plot_theme$legend.position.inside),
    " justification ", list_of(plot_theme$legend.justification)
  ),
  paste0("legend.title face ", plot_theme$legend.title$face, " size ", plot_theme$legend.title$size)
)

write_plot(
  rbind(ecdf_rows("test_data", test_data_plot), ecdf_rows("no_group_r", no_group_plot)), "ecdf_lines.tsv",
  c(
    "the rows of the layer data of ggplot_build, by group and then by x, as the step line goes",
    "group is the legend key of the row. The key NA is a gene in no list",
    "test_data is deseq_fc_ecdf(res, list(Background = background, mer7m8 = seed_genes)) for the test data",
    "no_group_r has one more row of res: row 1 with the gene_id ENSNOGROUP and the log2FoldChange 5",
    note("legend test_data", legend_of(test_data_plot)),
    note("legend no_group_r", legend_of(no_group_plot)),
    note("legend factor_order_reversed", legend_of(reversed_plot)),
    note("legend absent_and_repeated", legend_of(repeated_plot)),
    note("legend na_fold_change", legend_of(na_plot)),
    note("na_fold_change rows of mer7m8", sum(legend_keys(na_plot)[na_layer$group] == "mer7m8")),
    note("na_fold_change build warning", fields(na_built$warnings)),
    note("legend big_counts", legend_of(big_plot)),
    note("subtitle KS", test_data_plot$labels$subtitle),
    note("subtitle Wilcoxen", subtitle_of("Wilcoxen")),
    note("subtitle Wass", subtitle_of("Wass")),
    note("subtitle DTS", dts_subtitle),
    note("labels", fields(
      paste0("title=", test_data_plot$labels$title), paste0("x=", test_data_plot$labels$x),
      paste0("y=", test_data_plot$labels$y)
    )),
    note("layer", fields(
      paste0("geom ", class(layer_of_plot$geom)[1]), paste0("stat ", class(layer_of_plot$stat)[1]),
      paste0("linewidth ", layer_of_plot$aes_params$linewidth),
      paste0("colours ", list_of(unique(test_data_build$data[[1]]$colour)))
    )),
    note("coord limits", list_of(test_data_plot$coordinates$limits$x)),
    note("panel range", fields(paste0("x ", list_of(panel$x.range)), paste0("y ", list_of(panel$y.range)))),
    note("theme argument", paste0("theme = theme(legend.position = \"bottom\") gives legend.position ", bottom_plot$theme$legend.position)),
    note("ten groups", ten_groups_error),
    note("ggsave to a png name", png_magic),
    note("theme", theme_note),
    note("palette", list_of(SeedMatchR:::SeedMatchR.palette))
  )
)

# ---- ECDF labels ----
p_values <- list(0.05, 0.0149, 0.0105, 0.01, 0.0099999, 0.00999, 0.005, 0.0025, 0.001, 0.000999, 1e-10, 0.996, 1, NaN, NA, NULL)
label_res <- res
label_res$Group <- factor(ifelse(label_res$mer7m8 >= 1, "mer7m8", "Background"), levels = c("Background", "mer7m8"))
p_value_text <- function(p) {
  if (is.null(p)) "NULL" else if (is.nan(p)) "NaN" else if (is.na(p)) NA_character_ else double_to_string(p)
}
p_value_label <- function(p) {
  tryCatch({
    pl <- SeedMatchR::.plot_ecdf(
      res = label_res, title = "t", legend_labels = c("a", "b"),
      test_results = list(statistic = c(D = 0.123456), p.value = p)
    )
    sub("^p-value: ", "", strsplit(pl$labels$subtitle, "\n", fixed = TRUE)[[1]][2])
  }, error = function(e) paste0("ERROR: ", root_message(e)))
}
write_plot(
  data.frame(
    p_value = vapply(p_values, p_value_text, character(1)), label = vapply(p_values, p_value_label, character(1)),
    stringsAsFactors = FALSE
  ),
  "ecdf_labels.tsv",
  c(
    "label is the text after p-value: in the subtitle of .plot_ecdf, with test_results = list(statistic = c(D = 0.123456), p.value = p)",
    "p_value is the text of double_to_string, or NaN, NA, or NULL"
  )
)

# ---- Wasserstein ----
stats_samples <- utils::read.delim(
  file.path(GOLDEN_DIR, "stats/samples.tsv"), comment.char = "#", colClasses = "character"
)
sample_of <- function(case, sample) {
  text <- stats_samples$value[stats_samples$case == case & stats_samples$sample == sample]
  stopifnot(identical(sprintf("%.17g", as.numeric(text)), text))
  as.numeric(text)
}
wass_cases <- list(
  test_data = list(seed = res$log2FoldChange[res$mer7m8 >= 1], non_seed = res$log2FoldChange[res$mer7m8 == 0]),
  wass_40x300_ties = list(seed = sample_of("wass_40x300_ties", "x"), non_seed = sample_of("wass_40x300_ties", "y"))
)

GROUPS <- c("non-seed", "seed")
wass_ecdf <- list()
wass_density <- list()
wass_difference <- list()
wass_notes <- list()
wass_descriptions <- list()
for (case in names(wass_cases)) {
  seed <- wass_cases[[case]]$seed
  non_seed <- wass_cases[[case]]$non_seed
  full <- suppressWarnings(SeedMatchR::wass_dist(seed, non_seed))
  coarse <- suppressWarnings(SeedMatchR::wass_dist(seed, non_seed, n.interp = 100))

  cdf_build <- ggplot_build(full$cdf.plot)
  dens_build <- ggplot_build(full$dens.plot)
  diff_build <- ggplot_build(coarse$diff.plot)
  cdf_layer <- cdf_build$data[[1]]
  dens_layer <- dens_build$data[[1]]
  diff_layer <- diff_build$data[[1]]
  stopifnot(
    identical(legend_keys(full$cdf.plot), GROUPS),
    identical(as.character(get_guide_data(full$dens.plot, "fill")$.value), GROUPS)
  )
  cdf_layer <- cdf_layer[order(cdf_layer$group, cdf_layer$x), ]
  dens_layer <- dens_layer[order(dens_layer$group, dens_layer$x), ]

  wass_ecdf[[case]] <- data.frame(
    case = case, group = GROUPS[cdf_layer$group], x = cdf_layer$x, y = cdf_layer$y, stringsAsFactors = FALSE
  )
  wass_density[[case]] <- data.frame(
    case = case, group = GROUPS[dens_layer$group], x = dens_layer$x, density = dens_layer$density,
    stringsAsFactors = FALSE
  )
  wass_difference[[case]] <- data.frame(
    case = case, n_interp = 100L, x = diff_layer$x, difference = diff_layer$y, stringsAsFactors = FALSE
  )

  x_ranges <- lapply(list(cdf_build, dens_build, diff_build), function(b) b$layout$panel_params[[1]]$x.range)
  stopifnot(identical(x_ranges[[1]], x_ranges[[2]]), identical(x_ranges[[1]], x_ranges[[3]]))
  wass_notes[[case]] <- c(
    note(paste("bandwidth", case, "non-seed"), format(stats::bw.nrd0(non_seed), digits = 17)),
    note(paste("bandwidth", case, "seed"), format(stats::bw.nrd0(seed), digits = 17)),
    note(paste("title", case, "n_interp 100"), coarse$diff.plot$labels$title),
    note(paste("title", case, "n_interp 10000"), full$diff.plot$labels$title),
    note(paste("panel x range", case), list_of(x_ranges[[1]]))
  )

  scale_name <- function(build, axis) build$layout[[paste0("panel_scales_", axis)]][[1]]$name
  wass_descriptions[[case]] <- c(
    note("ecdf plot", fields(
      paste0("geom ", class(full$cdf.plot$layers[[1]]$geom)[1]), paste0("stat ", class(full$cdf.plot$layers[[1]]$stat)[1]),
      paste0("linewidth ", list_of(unique(cdf_layer$linewidth))),
      paste0("colours ", paste(GROUPS, unique(cdf_layer$colour[order(cdf_layer$group)]), collapse = ", ")),
      paste0("x ", scale_name(cdf_build, "x")), paste0("y ", scale_name(cdf_build, "y")),
      paste0("title ", or_none(full$cdf.plot$labels$title))
    )),
    note("density plot", fields(
      paste0("geom ", class(full$dens.plot$layers[[1]]$geom)[1]), paste0("stat ", class(full$dens.plot$layers[[1]]$stat)[1]),
      paste0("fills ", paste(GROUPS, unique(dens_layer$fill[order(dens_layer$group)]), collapse = ", ")),
      paste0("alpha ", list_of(unique(dens_layer$alpha))), paste0("outline ", list_of(unique(dens_layer$colour))),
      paste0("linewidth ", list_of(unique(dens_layer$linewidth))),
      paste0("x ", scale_name(dens_build, "x")), paste0("y ", scale_name(dens_build, "y")),
      paste0("title ", or_none(full$dens.plot$labels$title))
    )),
    note("difference plot", fields(
      paste0("geom ", class(coarse$diff.plot$layers[[1]]$geom)[1]), paste0("colour ", list_of(unique(diff_layer$colour))),
      paste0("linewidth ", list_of(unique(diff_layer$linewidth))),
      paste0("x ", scale_name(diff_build, "x")), paste0("y ", scale_name(diff_build, "y"))
    ))
  )
}
stopifnot(identical(wass_descriptions[[1]], wass_descriptions[[2]]))
wass_pointer <- "the notes of the Wasserstein plots are in wasserstein_density.tsv"

write_plot(
  do.call(rbind, wass_ecdf), "wasserstein_ecdf.tsv",
  c(
    "the rows of the layer data of cdf.plot of wass_dist(seed, non.seed), by group and then by x",
    "test_data is the fold changes of the test data: seed is mer7m8 >= 1, non-seed is mer7m8 == 0",
    "wass_40x300_ties is the case of tests/golden/stats/samples.tsv: x is the seed sample, y is the non-seed sample",
    wass_pointer
  )
)
write_plot(
  do.call(rbind, wass_density), "wasserstein_density.tsv",
  c(
    "the rows of the layer data of dens.plot of wass_dist(seed, non.seed), by group and then by x",
    "the cases are the cases of wasserstein_ecdf.tsv",
    unlist(wass_notes, use.names = FALSE),
    wass_descriptions[[1]]
  )
)
write_plot(
  do.call(rbind, wass_difference), "wasserstein_difference.tsv",
  c(
    "the rows of the layer data of diff.plot of wass_dist(seed, non.seed, n.interp = 100)",
    "the cases are the cases of wasserstein_ecdf.tsv",
    wass_pointer
  )
)

# ---- K-mer histogram ----
# The calls of R/kmer_hist.R:19-28.
olig.freq.mat <- Biostrings::oligonucleotideFrequency(x = seqs, width = 7, step = 1, simplify.as = "matrix")
olig.freq.bin <- olig.freq.mat
olig.freq.bin[olig.freq.bin > 0] <- 1
gene.counts <- colSums(olig.freq.bin)

counts_table <- as.data.frame(table(gene.counts))
hist_value <- graphics::hist(gene.counts, breaks = 50, plot = FALSE)

kmer_pdf <- file.path(tempdir(), "kmer.pdf")
reset_traced()
grDevices::pdf(kmer_pdf, compress = FALSE)
kmer_return <- with_tracers("title", SeedMatchR::kmer_hist(seqs, "CTCTATA"))
invisible(grDevices::dev.off())
kmer_titles <- traced_calls$title
pdf_lines <- readLines(kmer_pdf, warn = FALSE)
unlink(kmer_pdf)
page_count <- sum(grepl("/Type /Page ", pdf_lines, fixed = TRUE, useBytes = TRUE))
line_widths <- grep("^[0-9.]+ w$", pdf_lines, value = TRUE, useBytes = TRUE)

write_plot(
  data.frame(
    gene_count = as.integer(as.character(counts_table$gene.counts)), kmers = counts_table$Freq,
    stringsAsFactors = FALSE
  ),
  "kmer_hist.tsv",
  c(
    "gene_count is the number of genes with a 7-mer, and kmers is the number of 7-mers with this gene count: as.data.frame(table(gene.counts))",
    "gene.counts is colSums of the presence matrix of R/kmer_hist.R:19-28, for Test_seqs.fa with width 7",
    note("hist breaks", list_of(hist_value$breaks)),
    note("hist counts", list_of(hist_value$counts)),
    vapply(
      c("CTCTATA", "CUCUAUA", "ctctata", "CTCTATAA", "NOTAKMER"),
      function(target) note(paste("target", target), gene.counts[target][[1]]), character(1), USE.NAMES = FALSE
    ),
    note("return value", paste(deparse(kmer_return), collapse = "")),
    note("pages", page_count),
    note("title calls", do.call(fields, lapply(kmer_titles, title_text))),
    # abline draws after the histogram, thus its width is the last width in the file.
    note("abline width", tail(line_widths, 1))
  )
)

# ---- Alignment ----
alignment_cases <- utils::read.delim(
  file.path(GOLDEN_DIR, "analyses/alignment.tsv"), comment.char = "#", colClasses = "character", quote = ""
)
alignment_row <- alignment_cases[alignment_cases$case == "mer7m8_1", ][1, ]
alignment_match <- data.frame(
  start = as.integer(alignment_row$start), end = as.integer(alignment_row$end), names = alignment_row$seq_id
)
with_seedmatchr_patch(file.path(PATCH_DIR, "plotAlignment.R"), {
  alignment_df <- SeedMatchR::plotAlignment(alignment_match, GUIDE, seqs, output.type = "data.frame")
  alignment_plot <- SeedMatchR::plotAlignment(alignment_match, GUIDE, seqs, output.type = "ggplot")
})
alignment_layer <- ggplot_build(alignment_plot)$data[[1]]
alignment_lines <- strsplit(alignment_layer$label, "\n", fixed = TRUE)[[1]]
stopifnot(length(alignment_lines) == 3)

write_plot(
  data.frame(case = "mer7m8_1", line = 1:3, text = alignment_lines, stringsAsFactors = FALSE),
  "alignment_plot.tsv",
  c(
    paste0(
      "text is the line of the label of the layer of plotAlignment(match, GUIDE, seqs, output.type = \"ggplot\") of ",
      "tools/golden/patches/plotAlignment.R. The match is ", alignment_row$seq_id, " from ",
      alignment_row$start, " to ", alignment_row$end, ", the case mer7m8_1 of analyses/alignment.tsv"
    ),
    note("text", fields(
      paste0("family ", alignment_layer$family), paste0("face ", alignment_layer$fontface),
      paste0("size ", alignment_layer$size), paste0("hjust ", alignment_layer$hjust),
      paste0("vjust ", alignment_layer$vjust), paste0("colour ", alignment_layer$colour)
    )),
    note("theme", fields(
      paste0("panel.background ", class(alignment_plot$theme$panel.background)[1]),
      paste0("equal to theme_void() ", isTRUE(all.equal(alignment_plot$theme, ggplot2::theme_void())))
    )),
    note("data.frame row names", paste0("[", rownames(alignment_df), "]", collapse = ",")),
    note("unpatched", error_of(
      SeedMatchR::plotAlignment(alignment_match, GUIDE, seqs, output.type = "ggplot")
    ))
  )
)

# ---- Ranges ----
eight <- data.frame(start = seq(10, 150, by = 20))
eight$end <- eight$start + 6
named_rows <- eight[1:3, ]
rownames(named_rows) <- c("a", "b", "c")

ranges_of <- function(ir, ...) {
  made <- with_warnings(SeedMatchR:::plotRangesDF(ir, ...))
  built <- with_warnings(ggplot_build(made$value))
  list(layer = built$value$data[[1]], warnings = unique(c(made$warnings, built$warnings)))
}
range_cases <- list(
  eight_rows = eight, subset = eight[c(3, 7, 8), ], text_row_names = named_rows,
  one_base = data.frame(start = c(5, 20), end = c(10, 20))
)
range_results <- lapply(range_cases, ranges_of)
range_rows <- do.call(rbind, lapply(names(range_results), function(case) {
  layer <- range_results[[case]]$layer
  data.frame(case = case, xmin = layer$xmin, xmax = layer$xmax, ymin = layer$ymin, ymax = layer$ymax, stringsAsFactors = FALSE)
}))

options_plot <- SeedMatchR:::plotRangesDF(eight, xlim = c(0, 1000), main = "My title", col = "red")
options_layer <- ggplot_build(options_plot)$data[[1]]
tall_layer <- ranges_of(eight, bar_height = 0.8)$layer

iranges <- IRanges::IRanges(start = c(5, 20, 40), end = c(10, 20, 60))
reset_traced()
iranges_run <- with_tracers(c("arrows", "title"), with_warnings(SeedMatchR:::plotIRanges(iranges)))
iranges_arrows <- vapply(
  traced_calls$arrows, function(a) paste0(a$x0, "-", a$x1, "@", a$y0, " col=", a$col), character(1)
)
iranges_titles <- traced_calls$title
empty_error <- error_of(suppressWarnings(SeedMatchR:::plotIRanges(IRanges::IRanges())))
data_frame_error <- error_of(suppressWarnings(SeedMatchR:::plotIRanges(eight)))

write_plot(
  range_rows, "ranges_plot.tsv",
  c(
    "the rows of the layer data of ggplot_build(SeedMatchR:::plotRangesDF(ir)), for the cases below",
    "eight is a table with the starts 10, 30, ..., 150 and the ends start + 6",
    "eight_rows is eight. subset is eight[c(3, 7, 8), ]. text_row_names is eight[1:3, ] with the row names a, b, and c",
    "one_base is data.frame(start = c(5, 20), end = c(10, 20))",
    note("options", fields(
      paste0("title ", or_none(options_plot$labels$title)), paste0("coord limits ", or_none(options_plot$coordinates$limits$x)),
      paste0("layer rows ", nrow(options_layer))
    )),
    note("options call", "plotRangesDF(eight, xlim = c(0, 1000), main = \"My title\", col = \"red\")"),
    note("fill", list_of(unique(range_results$eight_rows$layer$fill))),
    note("bar_height 0.8", fields(paste0("ymin ", tall_layer$ymin[1]), paste0("ymax ", tall_layer$ymax[1]))),
    note("text_row_names warning", fields(range_results$text_row_names$warnings)),
    note("plotIRanges arrows", fields(iranges_arrows)),
    note("plotIRanges titles", do.call(fields, lapply(iranges_titles, title_text))),
    note("plotIRanges one_base warning", fields(iranges_run$warnings)),
    note("plotIRanges empty", empty_error),
    note("plotIRanges data.frame", data_frame_error)
  )
)

# ---- PR and ROC ----
LABELS <- rep(c(1, 0), 10)
run_plot <- function(curve) {
  reset_traced()
  error <- tryCatch({
    suppressWarnings(with_tracers(c("segments", "title", "image.default", "axis"), plot(curve)))
    NULL
  }, error = root_message)
  list(error = error, segments = traced_calls$segments, title = traced_calls$title, image = traced_calls$image, axis = traced_calls$axis)
}
segment_rows <- function(case, curve_name, run) {
  stopifnot(length(run$segments) == 1)
  s <- run$segments[[1]]
  n <- length(s$x0)
  data.frame(
    case = case, curve = curve_name, segment = seq_len(n), x0 = s$x0, y0 = s$y0, x1 = s$x1, y1 = s$y1,
    color = rep_len(s$col, n), stringsAsFactors = FALSE
  )
}
color_bar <- function(run) {
  image <- run$image[[1]]
  fields(paste0("first y ", image$y[1]), paste0("last y ", image$y[length(image$y)]), paste0("colors ", length(image$col)))
}
curve_of <- function(scores) PRROC::pr.curve(scores.class0 = scores, weights.class0 = LABELS, curve = TRUE)

benchmark <- with_seedmatchr_patch(
  file.path(PATCH_DIR, "benchmark_preds.R"), SeedMatchR::benchmark_preds(res, padj = 0.15, l2fc = 0)
)
pr_run <- run_plot(benchmark$PRC)
roc_run <- run_plot(benchmark$ROC)
unpatched_run <- run_plot(SeedMatchR::benchmark_preds(res)$ROC)
one_score_2_run <- run_plot(curve_of(rep(2, 20)))
stopifnot(is.null(pr_run$error), is.null(roc_run$error), is.null(one_score_2_run$error))
stopifnot(
  identical(pr_run$segments[[1]]$lwd, roc_run$segments[[1]]$lwd),
  identical(pr_run$segments[[1]]$usr, roc_run$segments[[1]]$usr),
  identical(pr_run$axis, roc_run$axis)
)

score_scale <- eval(formals(PRROC:::plot.PRROC)$scale.color)
error_run <- function(scores) {
  run <- run_plot(curve_of(scores))
  fields(
    paste0("PR error ", run$error),
    paste0("segment colors ", list_of(unique(rep_len(run$segments[[1]]$col, length(run$segments[[1]]$x0)))))
  )
}

write_plot(
  rbind(
    segment_rows("benchmark_1", "PR", pr_run), segment_rows("benchmark_1", "ROC", roc_run),
    segment_rows("one_score_2", "PR", one_score_2_run)
  ),
  "benchmark_plot.tsv",
  c(
    "the arguments of graphics::segments in plot(res$PRC) and plot(res$ROC) of PRROC::plot.PRROC, with segment the 1-based position in the call",
    "benchmark_1 is benchmark_preds(res, padj = 0.15, l2fc = 0) of tools/golden/patches/benchmark_preds.R, for the test data",
    "one_score_2 is pr.curve(scores.class0 = rep(2, 20), weights.class0 = rep(c(1, 0), 10), curve = TRUE), with only the curve PR",
    note("title PR", pr_run$title[[1]]$main),
    note("title ROC", roc_run$title[[1]]$main),
    note("labels PR", fields(paste0("xlab ", pr_run$title[[1]]$xlab), paste0("ylab ", pr_run$title[[1]]$ylab))),
    note("labels ROC", fields(paste0("xlab ", roc_run$title[[1]]$xlab), paste0("ylab ", roc_run$title[[1]]$ylab))),
    note("lwd", pr_run$segments[[1]]$lwd),
    note("usr", list_of(pr_run$segments[[1]]$usr)),
    note("axis sides", list_of(vapply(pr_run$axis, function(a) a$side, numeric(1)))),
    note("scale", do.call(fields, lapply(c(1, 50, 53, 100), function(i) paste(i, score_scale[i])))),
    note("color bar PR", color_bar(pr_run)),
    note("color bar ROC", color_bar(roc_run)),
    note("title ROC unpatched", unpatched_run$title[[1]]$main),
    note("one_score_2 color bar", color_bar(one_score_2_run)),
    note("one score 0", error_run(rep(0, 20))),
    note("one score -3", error_run(rep(-3, 20))),
    note("an Inf score", error_run(c(Inf, rep(1, 9), rep(0, 10))))
  )
)

# ---- Report table ----
report_table <- function(guide) {
  SeedMatchR::SeedMatchReport(res = res, seqs = seqs, guide.seq = guide)$table
}
cells_of <- function(html, pattern, group) {
  found <- regmatches(html, gregexpr(pattern, html))[[1]]
  matrix(
    vapply(regmatches(found, regexec(pattern, found)), function(parts) parts[-1][seq_len(group)], character(group)),
    nrow = group
  )
}
match_of <- function(html, pattern) {
  found <- regmatches(html, regexec(pattern, html))[[1]]
  stopifnot(length(found) > 1)
  found[-1]
}
spanner_label <- function(tab, prefix) {
  spanners <- gt:::dt_spanners_get(tab)
  hit <- vapply(spanners$vars, function(vars) startsWith(vars[1], prefix), logical(1))
  stopifnot(sum(hit) == 1)
  list(label = spanners$spanner_label[[which(hit)]], columns = spanners$vars[[which(hit)]])
}
spanner_note <- function(key, spanner) note(key, fields(paste0("label ", spanner$label), paste0("columns ", list_of(spanner$columns))))

tab <- report_table(data$guide)
html <- gt::as_raw_html(tab, inline_css = FALSE)

boxhead <- gt:::dt_boxhead_get(tab)
vars <- boxhead$var[boxhead$type == "default"]
ids <- ifelse(grepl("^[0-9]", vars), paste0("a", vars), vars)

body <- cells_of(html, '<td headers="([^"]*)" class="([^"]*)"[^>]*>([^<]*)</td>', 3)
stubs <- cells_of(html, '<th id="(stub_[^"]*)" scope="row" class="[^"]*">([^<]*)</th>', 2)
header_parts <- strsplit(body[1, ], " ", fixed = TRUE)
stopifnot(all(lengths(header_parts) == 3))
stub_ids <- vapply(header_parts, `[`, "", 2)
column_ids <- vapply(header_parts, `[`, "", 3)
columns <- vars[match(column_ids, ids)]
stopifnot(!anyNA(columns), all(stub_ids %in% stubs[1, ]))
row_names <- stubs[2, match(stub_ids, stubs[1, ])]
row_numbers <- as.integer(sub("^stub_[0-9]+_", "", stub_ids))

styles <- gt:::dt_styles_get(gt:::build_data(tab, context = "html"))
stopifnot(all(styles$locname == "data"))
style_of <- function(row, column, part) {
  hit <- which(styles$rownum == row & styles$colname == column)
  if (length(hit) == 0) return(NA_character_)
  stopifnot(length(hit) == 1)
  styles$styles[[hit]][[part]]$color
}
report_rows <- data.frame(
  row = row_names, column = columns, text = body[3, ],
  background = mapply(style_of, row_numbers, columns, "cell_fill", USE.NAMES = FALSE),
  color = mapply(style_of, row_numbers, columns, "cell_text", USE.NAMES = FALSE),
  stringsAsFactors = FALSE
)

title_cell <- match_of(html, '<td colspan="([0-9]+)" class="gt_heading gt_title[^"]*"[^>]*>([^<]*)</td>')
subtitle_cell <- match_of(html, '<td colspan="([0-9]+)" class="gt_heading gt_subtitle[^"]*"[^>]*>([^<]*)</td>')
group_cell <- match_of(html, '<th colspan="[0-9]+" class="gt_group_heading"[^>]*>([^<]*)</th>')
stopifnot(title_cell[1] == subtitle_cell[1])

one_cell_table <- function(table) {
  match_of(gt::as_raw_html(table, inline_css = FALSE), '<td[^>]*>([^<]*)</td>')
}
number_cell <- one_cell_table(
  gt::fmt_number(gt::gt(data.frame(a = 12345)), columns = "a", decimals = 0)
)
percent_cell <- one_cell_table(
  gt::fmt_percent(gt::gt(data.frame(b = 0.125)), columns = "b", decimals = 2, scale_values = FALSE)
)
# The minus sign of gt is U+2212. A note writes it as the text U+2212, thus the file stays ASCII.
ascii_minus <- function(text) gsub(intToUtf8(0x2212), "U+2212", text, fixed = TRUE)
special_number_cell <- function(value) {
  ascii_minus(one_cell_table(
    gt::fmt_number(gt::gt(data.frame(a = value)), columns = "a", decimals = 0)
  ))
}
special_percent_cell <- function(value) {
  ascii_minus(one_cell_table(
    gt::fmt_percent(gt::gt(data.frame(b = value)), columns = "b", decimals = 2, scale_values = FALSE)
  ))
}
special_numbers <- c(`Inf` = Inf, `-Inf` = -Inf, `-0.4` = -0.4, `-5` = -5, `-1234.4` = -1234.4, `-0` = -0, `NaN` = NaN, `NA` = NA_real_)
special_percents <- c(
  `Inf` = Inf, `-Inf` = -Inf, `-0.4` = -0.4, `-5` = -5, `-0.004` = -0.004, `1e-9` = 1e-9,
  `12345.678` = 12345.678, `-0` = -0, `NaN` = NaN, `NA` = NA_real_
)
special_notes <- c(
  vapply(
    names(special_numbers),
    function(label) note(paste("fmt_number", label), special_number_cell(special_numbers[[label]])),
    character(1), USE.NAMES = FALSE
  ),
  vapply(
    names(special_percents),
    function(label) note(paste("fmt_percent", label), special_percent_cell(special_percents[[label]])),
    character(1), USE.NAMES = FALSE
  )
)
short_guide_tab <- report_table(substr(data$guide, 1, 21))

write_plot(
  report_rows, "report_table.tsv",
  c(
    "the cells of gt::as_raw_html(tab, inline_css = FALSE) for tab <- SeedMatchReport(res, seqs, guide)$table, in the order of the rows and then of the columns",
    "row is the stub text of the row, and column is the column name without the prefix a that gt adds to an ID that starts with a digit",
    "background and color are the cell_fill and the cell_text of gt:::dt_styles_get(gt:::build_data(tab, context = \"html\")), or NA without a style",
    note("title", title_cell[2]),
    note("subtitle", subtitle_cell[2]),
    note("group", group_cell),
    note("columns", list_of(columns[row_names == row_names[1]])),
    note("colspan", title_cell[1]),
    spanner_note("spanner Full", spanner_label(tab, "Full")),
    spanner_note("spanner 18mer", spanner_label(tab, "18mer")),
    spanner_note("spanner 15mer", spanner_label(tab, "15mer")),
    note("id", match_of(html, '^<div id="([^"]*)"')),
    note("css selector", trimws(match_of(html, "<style>([^{]*)\\{"))),
    note("striped rows", list_of(sort(unique(row_numbers[grepl("gt_striped", body[2, ])])))),
    note("spanner Full of a 21-base guide", spanner_label(short_guide_tab, "Full")$label),
    note("fmt_number 12345", number_cell),
    note("fmt_percent 0.125", percent_cell),
    special_notes
  )
)
