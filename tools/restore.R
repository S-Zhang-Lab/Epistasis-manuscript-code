# Restore a separate package library for each requested figure module.
#
# Each module is restored in its own subprocess whose working directory is that
# module, because renv resolves the project library from the working directory.
# The subprocess activates the project and then loads it: activate writes the
# module's renv/activate.R, and load puts the project library on .libPaths().
# Without the load step renv compares the lockfile against the caller's library
# instead, which on a machine that already has the packages reports "already
# synchronized" and installs nothing. Writing renv/activate.R also means the
# module's .Rprofile selects the pinned library afterwards, for anyone running
# an individual script.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || any(!args %in% as.character(1:5)))
  stop("Usage: Rscript tools/restore.R 1 [2 3 4 5]")

if (!requireNamespace("renv", quietly = TRUE))
  install.packages("renv", repos = "https://cloud.r-project.org")

root <- normalizePath(getwd())
rscript <- file.path(R.home("bin"), "Rscript")
failed <- character()

for (n in args) {
  project <- normalizePath(file.path(root, "figures", paste0("figure", n)),
                           mustWork = TRUE)
  message("[restore] figure ", n, " -> ", project)
  setwd(project)
  status <- tryCatch(
    system2(rscript, c("--vanilla", "-e",
                       shQuote(paste("renv::activate(project = '.');",
                                     "renv::load(project = '.');",
                                     "renv::restore(project = '.', prompt = FALSE)")))),
    error = function(e) 1L,
    finally = setwd(root)
  )
  if (!identical(as.integer(status), 0L)) {
    failed <- c(failed, paste0("figure", n))
    message("[restore] FAILED: figure ", n)
  } else {
    message("[restore] OK: figure ", n)
  }
}

if (length(failed)) {
  stop("renv::restore failed for: ", paste(failed, collapse = ", "), "\n",
       "  Packages that build from source need a working toolchain. On macOS a\n",
       "  mismatched gfortran is the usual cause; check that the paths reported\n",
       "  by `R CMD config FLIBS` exist. Compatible prebuilt binaries can avoid\n",
       "  compilation where available for the requested R version and platform.",
       call. = FALSE)
}

message("[restore] All requested modules restored.")
