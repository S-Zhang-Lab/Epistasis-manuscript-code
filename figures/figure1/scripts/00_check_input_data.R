source(here::here("R", "paths.R"))

PANCAN_STUDIES <- c("brca", "luad", "lusc", "coadread", "stad",
                    "hnsc", "blca", "prad", "ov", "lihc")

REQUIRED_RAW <- c(

  file.path("GSE110590", c(
    "GSE110590-GPL11154_series_matrix.txt.gz",
    "GSE110590-GPL16791_series_matrix.txt.gz",
    "GSE110590_RAP_A16_log2.sne.tsv.gz")),

  file.path("MSK_MET_2021_breast", c(
    "clinical.csv.gz", "mutations.csv.gz", "cna.csv.gz")),

  file.path("TCGA_PanCan",
            paste0(PANCAN_STUDIES, "_tcga_pan_can_atlas_2018"),
            "data_mutations.txt.gz")
)

accepts_plain <- function(files) grepl("TCGA_PanCan", files, fixed = TRUE)

REQUIRED_DERIVED <- c("Sample_Metadata_Table.csv", "dawnrank_panel_subset.csv")

status_table <- function(root, files) {
  path <- file.path(root, files)
  present <- file.exists(path)

  alt <- accepts_plain(files) & grepl("\\.gz$", files)
  if (any(alt)) {
    plain <- sub("\\.gz$", "", path[alt])
    present[alt] <- present[alt] | file.exists(plain)
  }
  data.frame(file = files, path = path, present = present,
             stringsAsFactors = FALSE)
}

report_missing <- function(label, st) {
  missing <- st[!st$present, , drop = FALSE]
  if (nrow(missing) == 0L) {
    message(sprintf("[00_check_input_data] OK %s: %d file(s)", label, nrow(st)))
    return(invisible(character()))
  }
  message(sprintf("[00_check_input_data] MISSING %s:", label))
  for (p in missing$path) message("  - ", p)
  missing$path
}

verify_manifest <- function() {
  manifest <- NULL
  for (cand in c(raw("MANIFEST.sha256"), raw(".checksums.txt"))) {
    if (file.exists(cand)) { manifest <- cand; break }
  }
  if (is.null(manifest)) {
    message("[00_check_input_data] No checksum manifest found; verification skipped.")
    return(invisible(TRUE))
  }
  lines <- readLines(manifest, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines) & nzchar(trimws(lines))]
  if (!length(lines)) {
    message("[00_check_input_data] Checksum manifest is empty; verification skipped.")
    return(invisible(TRUE))
  }
  bad <- character()
  for (ln in lines) {
    parts <- strsplit(trimws(ln), "\\s+")[[1]]
    if (length(parts) < 2) next
    sha <- parts[1]
    rel <- paste(parts[-1], collapse = " ")
    f   <- file.path(REPO_ROOT, rel)
    if (!file.exists(f)) f <- raw(basename(rel))
    if (!file.exists(f)) next
    got <- tools::sha256sum(f)
    if (!identical(tolower(unname(got)), tolower(sha))) bad <- c(bad, rel)
  }
  if (length(bad)) {
    stop("[00_check_input_data] Checksum mismatch against ", basename(manifest),
         ":\n  ", paste(bad, collapse = "\n  "), call. = FALSE)
  }
  message("[00_check_input_data] OK checksum manifest (", basename(manifest), ")")
  invisible(TRUE)
}

SCHEMA_CONTRACT <- c(
  list(
    list(path = file.path("MSK_MET_2021_breast", "clinical.csv.gz"),
         cols = c("sampleId","CANCER_TYPE","SAMPLE_TYPE","METASTATIC_SITE","MET_SITE_COUNT",
                  "SUBTYPE_ABBREVIATION","DMETS_DX_BONE","DMETS_DX_LIVER","DMETS_DX_LUNG",
                  "DMETS_DX_CNS_BRAIN"),
         min_rows = 2000L),
    list(path = file.path("MSK_MET_2021_breast", "mutations.csv.gz"),
         cols = c("sampleId","hugo","entrez","mutationType"), min_rows = 5000L),
    list(path = file.path("MSK_MET_2021_breast", "cna.csv.gz"),
         cols = c("sampleId","hugo","entrez","alteration"), min_rows = 500000L)
  ),

  lapply(paste0(PANCAN_STUDIES, "_tcga_pan_can_atlas_2018"), function(st)
    list(path = file.path("TCGA_PanCan", st, "data_mutations.txt.gz"),
         cols = c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode",
                  "Chromosome", "Start_Position", "Variant_Type"),
         min_rows = 500L, sep = "\t"))
)

check_schema <- function() {
  problems <- character()
  for (s in SCHEMA_CONTRACT) {
    f <- file.path(DATA_RAW, s$path)

    if (!file.exists(f) && accepts_plain(s$path)) f <- sub("\\.gz$", "", f)
    if (!file.exists(f)) next
    sep <- if (is.null(s$sep)) "," else s$sep
    open_con <- function() if (grepl("\\.gz$", f)) gzfile(f) else file(f)

    head_lines <- tryCatch(readLines(open_con(), n = 200L, warn = FALSE),
                           error = function(e) character())
    head_body <- head_lines[!grepl("^#", head_lines)]
    hdr <- if (length(head_body))
      trimws(gsub('"', "", strsplit(head_body[1], sep, fixed = TRUE)[[1]])) else character()
    miss <- setdiff(s$cols, hdr)
    if (length(miss))
      problems <- c(problems, sprintf("%s: missing column(s) %s", s$path, paste(miss, collapse = ", ")))

    if (file.size(f) <= 20e6) {
      body <- tryCatch(readLines(open_con(), warn = FALSE), error = function(e) character())
      body <- body[!grepl("^#", body)]
      n <- length(body) - 1L
      if (n >= 0L && n < s$min_rows)
        problems <- c(problems, sprintf("%s: %d data rows, expected at least %d (wrong release or truncated file?)",
                                        s$path, n, s$min_rows))
    }
  }
  if (length(problems)) {
    message("[00_check_input_data] SCHEMA problems:")
    for (p in problems) message("  - ", p)
  } else {
    message("[00_check_input_data] OK schema contract: ", length(SCHEMA_CONTRACT), " file(s)")
  }
  problems
}

raw_status     <- status_table(DATA_RAW, REQUIRED_RAW)
derived_status <- status_table(DATA_DERIVED, REQUIRED_DERIVED)

missing <- c(
  report_missing("data/raw Tier-2 public inputs", raw_status),
  report_missing("data/derived Tier-1 curated inputs", derived_status)
)

verify_manifest()
schema_problems <- check_schema()

if (length(schema_problems) && !length(missing)) {
  stop("[00_check_input_data] Input files are present but do not match the expected schema (",
       length(schema_problems), " problem(s) listed above).\n",
       "  Every input listed here ships with the repository, so this normally means\n",
       "  a local edit or a partial checkout. Restore them with:\n",
       "      git checkout -- data/raw/\n",
       "  Provenance for how each was derived: scripts/etl_msk_met_breast.R (MSK-MET),\n",
       "  scripts/subset_tcga_pancan.sh (Panel B MAFs).",
       call. = FALSE)
}

if (length(missing)) {
  stop("[00_check_input_data] Required Figure 1 inputs are missing (",
       length(missing), " file(s)).\n",
       "These all ship with the repository, so a complete clone should never reach\n",
       "this message. Restore them with:\n",
       "  (a) git checkout -- data/raw/           (recover a deleted or edited input), or\n",
       "  (b) re-clone the repository.\n",
       "See data/README.md for the full provenance manifest.",
       call. = FALSE)
}

message("[00_check_input_data] OK — all required Figure 1 inputs are present.")
