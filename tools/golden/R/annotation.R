source("tools/golden/R/common.R")
source("tools/golden/R/annotation_common.R")

# tools/fixtures/annotation.py writes the fixture.
FIXTURE <- file.path(DATA_DIR, "annotation")
SCRIPT <- "annotation.R"

txdb <- ensembldb::EnsDb(file.path(FIXTURE, "synthetic.sqlite"))
dna <- rtracklayer::TwoBitFile(file.path(FIXTURE, "synthetic.2bit"))
utils::assignInNamespace("load_species_anno_db", local_anno_db(txdb, dna), ns = "SeedMatchR")

version_notes <- vapply(
  c("ensembldb", "GenomicFeatures", "rtracklayer", "Rsamtools"),
  function(p) paste0(p, ": ", packageVersion(p)),
  character(1)
)
patch_path <- function(name) file.path("tools/golden/patches", name)

configs <- list(
  utr3 = list(feature = "3UTR"),
  utr5 = list(feature = "5UTR"),
  exons = list(feature = "exons"),
  cds = list(feature = "cds"),
  utr3_all = list(feature = "3UTR", args = list(canonical = FALSE)),
  utr3_min0 = list(feature = "3UTR", args = list(canonical = FALSE, min.feature.width = 0)),
  exons_all = list(feature = "exons", args = list(canonical = FALSE)),
  exons_everything = list(
    feature = "exons",
    args = list(canonical = FALSE, protein.coding = FALSE, standard.chroms = FALSE)
  ),
  utr3_longest = list(feature = "3UTR", args = list(canonical = FALSE, longest.utr = TRUE)),
  utr3_longest_noncoding = list(
    feature = "3UTR",
    args = list(canonical = FALSE, protein.coding = FALSE, longest.utr = TRUE),
    patch = "get_longest_utr.R",
    unpatched_error = TRUE
  ),
  utr3_union_r = list(feature = "3UTR", args = list(canonical = FALSE, reduce.features = TRUE)),
  utr3_union = list(
    feature = "3UTR",
    args = list(canonical = FALSE, reduce.features = TRUE),
    patch = "reduce_features.R"
  ),
  utr3_union_min0 = list(
    feature = "3UTR",
    args = list(canonical = FALSE, reduce.features = TRUE, min.feature.width = 0),
    patch = "reduce_features.R"
  ),
  exons_longest_r = list(feature = "exons", args = list(canonical = FALSE, longest.utr = TRUE)),
  tsl1 = list(feature = "3UTR", args = list(canonical = FALSE, transcript.support = 1)),
  gene_ids = list(
    feature = "3UTR",
    args = list(canonical = FALSE, gene_id = c("ENSSYNG00000000001", "ENSSYNG00000000002"))
  ),
  tx_ids = list(
    feature = "3UTR",
    args = list(canonical = FALSE, tx_id = c("ENSSYNT00000000002", "ENSSYNT00000000004"))
  ),
  symbols = list(feature = "3UTR", args = list(canonical = FALSE, symbol = c("Syn1", "Syn5"))),
  entrez = list(feature = "3UTR", args = list(canonical = FALSE, entrez_id = c(1006, 1007, 1001))),
  nmd = list(
    feature = "exons",
    args = list(
      canonical = FALSE, protein.coding = FALSE,
      add.filter = AnnotationFilter::AnnotationFilterList(
        AnnotationFilter::TxBiotypeFilter("nonsense_mediated_decay")
      )
    )
  ),
  lnc = list(
    feature = "exons",
    args = list(
      canonical = FALSE, protein.coding = FALSE,
      add.filter = AnnotationFilter::AnnotationFilterList(AnnotationFilter::GeneBiotypeFilter("lncRNA"))
    )
  )
)

load_config <- function(config) {
  call <- function(return_gene_name) {
    args <- c(
      list(reference.name = "rnor7", feature.type = config$feature),
      config$args,
      list(return_gene_name = return_gene_name)
    )
    suppressMessages(do.call(SeedMatchR::load_annotations, args))
  }
  list(by_tx = call(FALSE), by_gene = call(TRUE))
}

# One row for each range with a width above 0. R gives a range of width 0 to a transcript without
# the feature, and to an exon that ends at the CDS boundary.
features_table <- function(features) {
  ranges <- unlist(features, use.names = FALSE)
  element <- rep(seq_along(features), S4Vectors::elementNROWS(features))
  keep <- BiocGenerics::width(ranges) > 0
  ranges <- ranges[keep]
  element <- element[keep]
  data.frame(
    name = names(features)[element],
    rank = as.integer(stats::ave(element, element, FUN = seq_along)),
    seqnames = as.character(GenomicRanges::seqnames(ranges)),
    start = BiocGenerics::start(ranges),
    end = BiocGenerics::end(ranges),
    strand = as.character(BiocGenerics::strand(ranges))
  )
}

sequences_table <- function(result) {
  seqs <- result$by_tx$seqs
  if (!identical(unname(as.character(seqs)), unname(as.character(result$by_gene$seqs)))) {
    stop("return_gene_name changed the sequences")
  }
  data.frame(
    name_tx = names(seqs),
    name_gene = names(result$by_gene$seqs),
    width = BiocGenerics::width(seqs),
    sequence = unname(as.character(seqs))
  )
}

transcripts_table <- function(filter) {
  tx <- as.data.frame(ensembldb::transcripts(txdb, filter = filter))
  columns <- c(
    "tx_id", "gene_id", "seqnames", "start", "end", "strand", "tx_biotype", "tx_cds_seq_start",
    "tx_cds_seq_end", "tx_support_level", "tx_is_canonical"
  )
  tx <- tx[, columns]
  tx$seqnames <- as.character(tx$seqnames)
  tx$strand <- as.character(tx$strand)
  rownames(tx) <- NULL
  tx
}

tx_tables <- list()
for (name in names(configs)) {
  config <- configs[[name]]
  notes <- version_notes
  if (isTRUE(config$unpatched_error)) {
    error <- tryCatch(
      {
        load_config(config)
        NULL
      },
      error = conditionMessage
    )
    if (is.null(error)) stop(name, ": the unpatched function gave no error")
    notes <- c(notes, paste0("unpatched ", sub("\\.R$", "", config$patch), " stops: ", error))
  }
  if (is.null(config$patch)) {
    result <- load_config(config)
  } else {
    result <- with_seedmatchr_patch(patch_path(config$patch), load_config(config))
    notes <- c(notes, paste0("patched with tools/golden/patches/", config$patch))
  }
  prefix <- paste0("annotation/", name)
  write_golden(features_table(result$by_tx$features), paste0(prefix, "_features.tsv"), SCRIPT, notes)
  write_golden(sequences_table(result), paste0(prefix, "_sequences.tsv"), SCRIPT, notes)
  tx_tables[[name]] <- transcripts_table(result$by_tx$filter)
  write_golden(tx_tables[[name]], paste0("annotation/transcripts_", name, ".tsv"), SCRIPT, notes)
}

# Genome ranges. The soft-masked run and the 3' UTR of T10 come from the fixture files, thus the
# ranges follow the specification of tools/fixtures/annotation.py.
lengths <- GenomeInfoDb::seqlengths(dna)
fasta <- readLines(file.path(FIXTURE, "synthetic.fa"))
headers <- grep("^>", fasta)
contig1 <- paste(fasta[(headers[1] + 1):(headers[2] - 1)], collapse = "")
masked <- regexpr("[acgtn]+", contig1)
if (masked < 0) stop("contig 1 of synthetic.fa has no soft-masked run")
t10 <- range(unlist(ensembldb::threeUTRsByTranscript(
  txdb,
  filter = AnnotationFilter::TxIdFilter("ENSSYNT00000000010")
)))

ranges <- rbind(
  data.frame(label = "contig_X", contig = "X", start = 1L, end = lengths[["X"]]),
  data.frame(label = "contig_SCAF1.1", contig = "SCAF1.1", start = 1L, end = lengths[["SCAF1.1"]]),
  data.frame(label = "contig_1_head", contig = "1", start = 1L, end = 120L),
  data.frame(
    label = paste0("byte_offset_", 0:3), contig = "1",
    start = 1001L + 0:3, end = c(1012L, 1010L, 1017L, 1021L)
  ),
  data.frame(
    label = "t10_utr3_block_boundary", contig = "1",
    start = BiocGenerics::start(t10), end = BiocGenerics::end(t10)
  ),
  data.frame(
    label = "soft_masked_run", contig = "1",
    start = as.integer(masked) - 2L, end = as.integer(masked) + attr(masked, "match.length") + 1L
  ),
  data.frame(label = "width_1", contig = "1", start = 5002L, end = 5002L),
  data.frame(
    label = paste0("last_7_", names(lengths)), contig = names(lengths),
    start = unname(lengths) - 6L, end = unname(lengths)
  )
)
rownames(ranges) <- NULL
gr <- GenomicRanges::GRanges(ranges$contig, IRanges::IRanges(ranges$start, ranges$end), strand = "+")

twobit <- as.character(Biostrings::getSeq(dna, gr))
fa <- Rsamtools::FaFile(
  file.path(FIXTURE, "synthetic.fa.gz"),
  index = file.path(FIXTURE, "synthetic.fa.gz.fai"),
  gzindex = file.path(FIXTURE, "synthetic.fa.gz.gzi")
)
bgzf <- as.character(Rsamtools::scanFa(fa, gr))

# The UCSC library of rtracklayer gives the reason as a warning, then fails. The read goes last,
# after the other reads, because the failure leaves the C stack of R unbalanced.
v1_note <- character(0)
v1_warnings <- character(0)
twobit_v1 <- tryCatch(
  withCallingHandlers(
    as.character(Biostrings::getSeq(rtracklayer::TwoBitFile(file.path(FIXTURE, "synthetic.v1.2bit")), gr)),
    warning = function(w) {
      v1_warnings <<- c(v1_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ),
  error = function(e) {
    v1_note <<- paste0(
      "rtracklayer cannot read synthetic.v1.2bit: ",
      paste(c(v1_warnings, conditionMessage(e)), collapse = "; ")
    )
    rep("unsupported", nrow(ranges))
  }
)
ranges$twobit <- unname(twobit)
ranges$twobit_v1 <- unname(twobit_v1)
ranges$bgzf <- unname(bgzf)

if (!identical(ranges$twobit, ranges$bgzf)) stop("the 2bit and the BGZF sequences differ")
if (any(grepl("[a-z]", ranges$twobit))) stop("a sequence has lowercase letters")
write_golden(ranges, "annotation/genome_ranges.tsv", SCRIPT, c(version_notes, v1_note))

order <- tx_tables$exons_everything$tx_id
tie <- order[order %in% c("ENSSYNT00000000001", "ENSSYNT00000000002", "ENSSYNT00000000010")]
cat("tie order in transcripts_exons_everything: ", paste(tie, collapse = ", "), "\n", sep = "")
cat("twobit_v1: ", if (length(v1_note)) v1_note else "read", "\n", sep = "")
