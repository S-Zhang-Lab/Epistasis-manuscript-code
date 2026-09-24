source(here::here("R", "paths.R"))

REQUIRED_RAW <- c(
  "Cx3cr1_raw_counts.csv"
)

REQUIRED_DERIVED <- c(

  "pathwayGI_PTEN_CX3CL1.csv",
  "pathwayGI_CHEK2_CX3CL1.csv",

  "network_scaffold_Cx3cr1_nodes.csv",
  "network_scaffold_Cx3cr1_edges.csv"
)

OPTIONAL_DERIVED <- character(0)

REQUIRED_REFERENCE <- character(0)

status_table <- function(root, files) {
  data.frame(
    file = files,
    path = file.path(root, files),
    present = file.exists(file.path(root, files)),
    stringsAsFactors = FALSE
  )
}

report_missing <- function(label, tbl) {
  if (nrow(tbl) == 0L) {
    message(sprintf("[00_check_input_data] OK %s: (none required at this stage)", label))
    return(invisible(character()))
  }
  missing <- tbl[!tbl$present, , drop = FALSE]
  if (nrow(missing) == 0L) {
    message(sprintf("[00_check_input_data] OK %s: %d file(s)", label, nrow(tbl)))
    return(invisible(character()))
  }
  message(sprintf("[00_check_input_data] MISSING %s:", label))
  for (p in missing$path) message("  - ", p)
  missing$path
}

check_duplicate_generated_tables <- function() {
  derived_csv <- list.files(DATA_DERIVED, pattern = "\\.csv$", full.names = FALSE)
  output_csv  <- list.files(OUT_TABLES,   pattern = "\\.csv$", full.names = FALSE)
  common <- intersect(derived_csv, output_csv)
  if (length(common) == 0L) {
    message("[00_check_input_data] OK duplicate generated tables: none")
    return(invisible(TRUE))
  }
  drifted <- character()
  for (f in common) {
    same <- isTRUE(tools::md5sum(derived(f)) == tools::md5sum(tbl(f)))
    if (!same) drifted <- c(drifted, f)
  }
  if (length(drifted)) {
    stop("[00_check_input_data] Duplicate data/derived and output/tables CSVs disagree:\n  ",
         paste(drifted, collapse = "\n  "),
         "\nMove generated outputs out of data/derived or update the frozen snapshot metadata.",
         call. = FALSE)
  }
  message("[00_check_input_data] OK duplicate generated tables: byte-identical")
  invisible(TRUE)
}

verify_manifest <- function() {
  manifest <- raw("MANIFEST.sha256")
  if (!file.exists(manifest)) {
    message("[00_check_input_data] No data/raw/MANIFEST.sha256 found; checksum verification skipped.")
    return(invisible(TRUE))
  }
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("[00_check_input_data] Package 'digest' is required to verify MANIFEST.sha256.",
         call. = FALSE)
  }
  expected <- read.table(manifest, header = FALSE, col.names = c("sha", "file"),
                         stringsAsFactors = FALSE)
  bad <- character(); missing <- character()
  for (i in seq_len(nrow(expected))) {

    f <- if (grepl("/", expected$file[i], fixed = TRUE)) {
      here::here(expected$file[i])
    } else {
      raw(expected$file[i])
    }
    if (!file.exists(f)) { missing <- c(missing, expected$file[i]); next }
    got <- digest::digest(file = f, algo = "sha256")
    if (!identical(tolower(got), tolower(expected$sha[i])))
      bad <- c(bad, expected$file[i])
  }
  if (length(missing) || length(bad)) {
    pieces <- c(
      if (length(missing)) paste0("Missing manifest files:\n  ", paste(missing, collapse = "\n  ")),
      if (length(bad))     paste0("SHA-256 mismatch:\n  ",    paste(bad,     collapse = "\n  "))
    )
    stop("[00_check_input_data] MANIFEST.sha256 verification failed.\n",
         paste(pieces, collapse = "\n"), call. = FALSE)
  }
  message("[00_check_input_data] OK MANIFEST.sha256")
  invisible(TRUE)
}

raw_status       <- status_table(DATA_RAW,     REQUIRED_RAW)
derived_status   <- status_table(DATA_DERIVED, REQUIRED_DERIVED)
optional_status  <- status_table(DATA_DERIVED, OPTIONAL_DERIVED)
reference_status <- status_table(DATA_REF,     REQUIRED_REFERENCE)

for (i in seq_len(nrow(optional_status))) {
  if (optional_status$present[i]) {
    message(sprintf("[00_check_input_data] OK optional: %s",
                    optional_status$file[i]))
  } else {
    message(sprintf("[00_check_input_data] optional MISSING: %s (panel will skip)",
                    optional_status$file[i]))
  }
}

missing <- c(
  report_missing("data/raw inputs",                   raw_status),
  report_missing("data/derived curated inputs",       derived_status),
  report_missing("data/reference external references", reference_status)
)

check_duplicate_generated_tables()
verify_manifest()

if (length(missing)) {
  stop("[00_check_input_data] Required Figure 4 inputs are missing. ",
       "See data/README.md for provenance and where each file should live.",
       call. = FALSE)
}

message("[00_check_input_data] OK -- all required Figure 4 inputs are present.")
