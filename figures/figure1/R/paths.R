if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here", repos = "https://cloud.r-project.org")
}

BASE_DIR <- here::here()
REPO_ROOT <- BASE_DIR

DATA_DIR     <- here::here("data")
DATA_RAW     <- here::here("data", "raw")
DATA_DERIVED <- here::here("data", "derived")
DATA_REF     <- here::here("data", "reference")

INPUT_DIR <- DATA_RAW

OUTPUT_DIR     <- here::here("output")
OUTPUT_FIG_DIR <- here::here("output", "figures", "main")

OUTPUT_SUP_DIR     <- here::here("output", "figures", "supplementary")
OUTPUT_SUP_VAR_DIR <- here::here("output", "figures", "supplementary", "variants")
OUTPUT_SUP_DIA_DIR <- here::here("output", "figures", "supplementary", "diagnostics")
OUTPUT_TBL_DIR <- here::here("output", "tables")
OUTPUT_INT_DIR <- here::here("output", "intermediates")
OUTPUT_LOG_DIR <- here::here("output", "logs")

DOCS_DIR       <- here::here("docs")
GENE_LISTS_DIR <- here::here("docs", "gene_lists")

for (d in c(DATA_RAW, DATA_DERIVED, DATA_REF,
            OUTPUT_FIG_DIR, OUTPUT_SUP_DIR, OUTPUT_SUP_VAR_DIR, OUTPUT_SUP_DIA_DIR,
            OUTPUT_TBL_DIR, OUTPUT_INT_DIR, OUTPUT_LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

raw      <- function(...) file.path(DATA_RAW, ...)
derived  <- function(...) file.path(DATA_DERIVED, ...)
ref      <- function(...) file.path(DATA_REF, ...)
fig_main <- function(...) file.path(OUTPUT_FIG_DIR, ...)
fig_supp <- function(...) file.path(OUTPUT_SUP_DIR, ...)
tbl      <- function(...) file.path(OUTPUT_TBL_DIR, ...)
log_path <- function(...) file.path(OUTPUT_LOG_DIR, ...)
