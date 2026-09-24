# Restore a separate package library for each requested figure module.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || any(!args %in% as.character(1:5))) stop("Usage: Rscript scripts/restore.R 1 [2 3 4 5]")
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", repos = "https://cloud.r-project.org")
for (n in args) {
  project <- normalizePath(file.path("figures", paste0("figure", n)), mustWork = TRUE)
  renv::restore(project = project, lockfile = file.path(project, "renv.lock"), prompt = FALSE)
  renv::activate(project = project)
}
