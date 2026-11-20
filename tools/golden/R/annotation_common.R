# Helpers for the annotation golden scripts.

# load_annotations gets its EnsDb and 2bit file from load_species_anno_db, which reads AnnotationHub.
# This function returns a replacement that gives the local files. Put it into the SeedMatchR
# namespace with assignInNamespace, so that load_annotations runs unchanged.
local_anno_db <- function(txdb, dna) {
  function(reference.name, filter = NULL) {
    list(gtf = ensembldb::transcripts(txdb, filter = filter), dna = dna, txdb = txdb)
  }
}

ANNOTATION_VERSION_NOTES <- vapply(
  c("ensembldb", "GenomicFeatures", "rtracklayer", "Rsamtools"),
  function(p) paste0(p, ": ", packageVersion(p)),
  character(1)
)
patch_path <- function(name) file.path("tools/golden/patches", name)

ANNOTATION_CONFIGS <- list(
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

# Run load_annotations for the config, with the patch that the config names. Return the result of
# load_config and the notes for the golden files.
run_config <- function(name, config) {
  notes <- ANNOTATION_VERSION_NOTES
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
  list(result = result, notes = notes)
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
