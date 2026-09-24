DATA_GUIDE <- "../../DATA.md"

.get_script_dir <- function() {
  if (interactive()) {
    return(tryCatch(dirname(rstudioapi::getActiveDocumentContext()$path),
                    error = function(e) getwd()))
  }
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(fa) > 0 && nzchar(fa[1])) return(dirname(normalizePath(fa[1])))
  this_frame <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(this_frame) && nzchar(this_frame)) return(dirname(normalizePath(this_frame)))
  getwd()
}

SCRIPT_DIR   <- .get_script_dir()
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."))

REQUIRED_INPUTS <- c(
  "data/raw/Perturb_InVivo_merged.rds",
  "data/raw/Perturb_InVitro_HTO.rds",
  "data/derived/screens/Day3.gene_summary.txt",
  "data/derived/screens/Day18.gene_summary.txt",
  "data/derived/screens/2nd_DC_counts.csv",
  "data/derived/screens/plasmid_freq.txt",
  "data/derived/screens/cas12a_cutting_efficiency.csv",
  "data/derived/hallmark_pathways_v26.1.0_MM.csv"
)

REQUIRED_META <- list(
  "data/raw/Perturb_InVivo_merged.rds"  = c("perturbation", "lane", "hash.ID"),
  "data/raw/Perturb_InVitro_HTO.rds"    = c("perturbation", "pool")
)

message("Figure 3 (v7) — input data check")
message(sprintf("Project root: %s", PROJECT_ROOT))
message("")

missing <- character(0)
for (rel in REQUIRED_INPUTS) {
  abs_path <- file.path(PROJECT_ROOT, rel)
  if (!file.exists(abs_path)) {
    missing <- c(missing, rel); message(sprintf("[MISSING] %s", rel))
  } else {
    message(sprintf("[ok]      %s", rel))
  }
}

if (length(missing) > 0L) {
  stop(sprintf(
    paste0("00_check_input_data.R: %d / %d required input(s) missing.\n",
           "First missing: %s\n",
           "The two data/raw/*.rds must be supplied separately.\n",
           "See %s for required files; the derived/ tables are included."),
    length(missing), length(REQUIRED_INPUTS),
    file.path(PROJECT_ROOT, missing[1]), DATA_GUIDE))
}

message("")
message("All required inputs present. md5 checksums:")
hashes <- tools::md5sum(file.path(PROJECT_ROOT, REQUIRED_INPUTS))
for (i in seq_along(hashes)) {
  message(sprintf("  %s  %s", unname(hashes[i]), REQUIRED_INPUTS[i]))
}

manifest_path <- file.path(PROJECT_ROOT, "data", "raw_checksums.tsv")
if (file.exists(manifest_path)) {
  man <- utils::read.table(manifest_path, header = TRUE, sep = "\t",
                           comment.char = "#", stringsAsFactors = FALSE)
  observed <- stats::setNames(unname(hashes), REQUIRED_INPUTS)
  bad <- character(0)
  for (i in seq_len(nrow(man))) {
    p <- man$path[i]
    if (p %in% names(observed) && !is.na(observed[[p]]) &&
        !identical(observed[[p]], man$md5[i])) {
      bad <- c(bad, sprintf("  %s\n    expected %s\n    observed %s",
                            p, man$md5[i], observed[[p]]))
    }
  }
  if (length(bad) > 0L) {
    stop(sprintf(paste0(
      "Input checksum mismatch - the object(s) below are not the expected version:\n%s\n",
      "See %s for input provenance. Update data/raw_checksums.tsv only if the input ",
      "was intentionally revised."), paste(bad, collapse = "\n"), DATA_GUIDE))
  }
  message("  -> all listed checksums match data/raw_checksums.tsv")
} else {
  warning("data/raw_checksums.tsv not found; input identity was NOT verified.")
}

if (identical(Sys.getenv("FIG3_CHECK_COLUMNS"), "1")) {
  message("")
  message("FIG3_CHECK_COLUMNS=1 -> loading objects to validate metadata columns ...")
  for (rel in names(REQUIRED_META)) {
    obj <- readRDS(file.path(PROJECT_ROOT, rel))
    md  <- tryCatch(obj@meta.data, error = function(e) as.data.frame(obj[[]]))
    have <- colnames(md)
    need <- REQUIRED_META[[rel]]
    absent <- setdiff(need, have)
    reds <- tryCatch(names(obj@reductions), error = function(e) character(0))
    if (length(absent) > 0L) {
      stop(sprintf("[%s] missing metadata column(s): %s\n(present: %s)",
                   rel, paste(absent, collapse = ", "), paste(have, collapse = ", ")))
    }
    if (!"pca" %in% reds) {
      stop(sprintf("[%s] missing 'pca' reduction (present: %s)",
                   rel, paste(reds, collapse = ", ")))
    }
    message(sprintf("  [ok] %s  columns + pca reduction present (reductions: %s)",
                    rel, paste(reds, collapse = ", ")))
    rm(obj, md); gc(verbose = FALSE)
  }
}

message("")
message("00_check_input_data.R: PASS")
