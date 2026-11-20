source("tools/golden/R/common.R")

lib <- "tools/golden/.rlib"
src <- "../SeedMatchR"

git <- function(...) {
  out <- suppressWarnings(system2("git", c("-C", src, ...), stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(out, "status"))) {
    stop("git ", paste(c(...), collapse = " "), " failed in ", src, ": ", paste(out, collapse = "\n"))
  }
  out
}

commit <- git("rev-parse", "HEAD")
if (!startsWith(commit, SEEDMATCHR_COMMIT)) {
  stop("SeedMatchR is at ", commit, ", but the goldens need a commit that starts with ", SEEDMATCHR_COMMIT)
}
changes <- git("status", "--porcelain", "--untracked-files=no")
if (length(changes) > 0) {
  stop("SeedMatchR has changes to tracked files:\n", paste(changes, collapse = "\n"))
}

dir.create(lib, recursive = TRUE, showWarnings = FALSE)
lib <- normalizePath(lib)
.libPaths(c(lib, .libPaths()))
if (isTRUE(getOption("repos")["CRAN"] == "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = lib)
}

imports <- read.dcf(file.path(src, "DESCRIPTION"), fields = "Imports")[1, 1]
imports <- trimws(sub("\\(.*\\)", "", strsplit(imports, ",")[[1]]))
missing <- setdiff(imports, rownames(installed.packages()))
cat("Imports missing from the library path:", if (length(missing)) missing else "none", "\n")
if (length(missing) > 0) {
  BiocManager::install(missing, lib = lib, update = FALSE, ask = FALSE)
}

still_missing <- setdiff(imports, rownames(installed.packages()))
if (length(still_missing) > 0) {
  stop("Imports not installed: ", paste(still_missing, collapse = ", "))
}

# igraph is a Suggests of SeedMatchR, and the network goldens call it, thus the Imports list does not hold it.
# check_env.R stops unless the version is 2.3.x.
igraph_found <- "igraph" %in% rownames(installed.packages()) &&
  startsWith(as.character(packageVersion("igraph")), "2.3.")
cat("igraph 2.3.x on the library path:", if (igraph_found) "yes" else "no", "\n")
if (!igraph_found) {
  install.packages("igraph", lib = lib)
}

# PRROC is a Suggests of SeedMatchR, and benchmark_preds calls it. The goldens pin PRROC 1.4 by its digest.
prroc_md5 <- "229140aa270cdad059daf08783b083bb"
if (!("PRROC" %in% rownames(installed.packages(lib.loc = lib))) ||
    as.character(packageVersion("PRROC", lib.loc = lib)) != "1.4") {
  tarball <- file.path(tempdir(), "PRROC_1.4.tar.gz")
  urls <- c(
    "https://cran.r-project.org/src/contrib/PRROC_1.4.tar.gz",
    "https://cran.r-project.org/src/contrib/Archive/PRROC/PRROC_1.4.tar.gz"
  )
  got <- FALSE
  for (url in urls) {
    got <- tryCatch(utils::download.file(url, tarball, mode = "wb", quiet = TRUE) == 0, error = function(e) FALSE)
    if (got) break
  }
  if (!got) stop("PRROC 1.4 cannot be downloaded from ", paste(urls, collapse = " or "))
  if (unname(tools::md5sum(tarball)) != prroc_md5) {
    stop("PRROC_1.4.tar.gz has the MD5 ", tools::md5sum(tarball), ", expected ", prroc_md5)
  }
  install.packages(tarball, repos = NULL, type = "source", lib = lib)
}

install.packages(src, repos = NULL, type = "source", lib = lib, dependencies = FALSE)
