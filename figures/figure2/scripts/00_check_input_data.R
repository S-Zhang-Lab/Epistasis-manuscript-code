source(here::here("R", "paths.R"))

REQUIRED_RAW <- c(

  "R1_construct_counts.csv",
  "R1_plasmid_counts.csv"
)

REQUIRED_DERIVED <- c(
  "R2_GI_scores.xlsx",
  "R2_Day18.gene_summary.txt",
  "R2_Day18_legacy.gene_summary.txt",
  "R2_construct_counts.csv",
  "R2_Day3.gene_summary.txt",
  "R2_plasmid_counts.txt"
)

REQUIRED_REFERENCE <- c(
  "R1_library_gene_list.xlsx",
  "R2_STRING_network.tsv"
)

status_table <- function(root, files) {
  data.frame(
    file    = files,
    path    = file.path(root, files),
    present = file.exists(file.path(root, files)),
    stringsAsFactors = FALSE
  )
}

report_missing <- function(label, tbl) {
  missing <- tbl[!tbl$present, , drop = FALSE]
  if (nrow(missing) == 0L) {
    message(sprintf("[00_check_input_data] OK %s: %d file(s)", label, nrow(tbl)))
    return(invisible(character()))
  }
  message(sprintf("[00_check_input_data] MISSING %s:", label))
  for (p in missing$path) message("  - ", p)
  missing$path
}

missing <- c(
  report_missing("data/raw primary inputs", status_table(DATA_RAW, REQUIRED_RAW)),
  report_missing("data/derived curated inputs", status_table(DATA_DERIVED, REQUIRED_DERIVED)),
  report_missing("data/reference external references", status_table(DATA_REF, REQUIRED_REFERENCE))
)

if (length(missing)) {
  stop("[00_check_input_data] Required Figure 2 inputs are missing.\n",
       "Every input is tracked in git, so a complete clone should already have\n",
       "them. If files are absent, the checkout is incomplete: try\n",
       "  git status --short data/   and   git checkout -- data/\n",
       "See data/README.md for the full inventory.",
       call. = FALSE)
}

message("[00_check_input_data] OK -- all required Figure 2 inputs are present.")
