if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here", repos = "https://cloud.r-project.org")
}

REPO_ROOT     <- here::here()
DATA_DIR      <- here::here("data")
DATA_RAW      <- here::here("data", "raw")
DATA_DERIVED  <- here::here("data", "derived")
DATA_REF      <- here::here("data", "reference")
OUT_DIR       <- here::here("output")
OUT_FIGS      <- here::here("output", "figures")
OUT_TABLES    <- here::here("output", "tables")
OUT_LOGS      <- here::here("output", "logs")
DOCS_DIR      <- here::here("docs")

for (d in c(DATA_RAW, DATA_DERIVED, DATA_REF, OUT_FIGS, OUT_TABLES, OUT_LOGS)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

raw      <- function(...) file.path(DATA_RAW, ...)
derived  <- function(...) file.path(DATA_DERIVED, ...)
ref      <- function(...) file.path(DATA_REF, ...)
fig      <- function(...) file.path(OUT_FIGS, ...)
tbl      <- function(...) file.path(OUT_TABLES, ...)
log_path <- function(...) file.path(OUT_LOGS, ...)
