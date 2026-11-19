# Helpers for the annotation golden scripts.

# load_annotations gets its EnsDb and 2bit file from load_species_anno_db, which reads AnnotationHub.
# This function returns a replacement that gives the local files. Put it into the SeedMatchR
# namespace with assignInNamespace, so that load_annotations runs unchanged.
local_anno_db <- function(txdb, dna) {
  function(reference.name, filter = NULL) {
    list(gtf = ensembldb::transcripts(txdb, filter = filter), dna = dna, txdb = txdb)
  }
}

# Put each function of the patch file into the SeedMatchR namespace, evaluate code, and then put
# back the original functions. The patched functions get the namespace as the parent environment,
# thus they find the imports of SeedMatchR, for example %>%.
with_seedmatchr_patch <- function(patch, code) {
  ns <- asNamespace("SeedMatchR")
  patched <- new.env(parent = ns)
  sys.source(patch, envir = patched)
  names <- ls(patched, all.names = TRUE)
  originals <- mget(names, envir = ns)
  for (name in names) utils::assignInNamespace(name, patched[[name]], ns = "SeedMatchR")
  on.exit(for (name in names) utils::assignInNamespace(name, originals[[name]], ns = "SeedMatchR"))
  code
}
