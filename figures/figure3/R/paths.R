suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("R/paths.R requires the 'here' package. Install it or run renv::restore().")
  }
})

PROJECT_ROOT <- here::here()

.path <- function(...) {
  p <- file.path(PROJECT_ROOT, ...)
  normalizePath(p, mustWork = FALSE)
}

raw <- function(...) .path("data", "raw", ...)

derived <- function(...) .path("data", "derived", ...)

ref <- function(...) .path("data", "reference", ...)

fig <- function(...) .path("output", "figures", ...)

fig_main <- function(...) .path("output", "figures", "main", ...)

fig_supp <- function(...) .path("output", "figures", "supplementary", ...)

fig_archive <- function(...) .path("output", "figures", "archive", ...)

tbl <- function(...) .path("output", "tables", ...)

log_path <- function(...) .path("output", "logs", ...)

interm <- function(...) .path("output", "intermediates", ...)

cache <- function(...) .path("output", "cache", ...)

assets <- function(...) .path("assets", ...)

ensure_output_dirs <- function() {
  for (d in c(fig_main(), fig_supp(), fig_archive(), tbl(), log_path(), interm(), cache())) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(TRUE)
}

invisible(NULL)
