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

install.packages(src, repos = NULL, type = "source", lib = lib, dependencies = FALSE)
