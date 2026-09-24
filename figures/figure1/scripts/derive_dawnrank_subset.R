# Rebuild the DawnRank panel subset used by Panel J from its published source.
#
# Panel J reads data/derived/dawnrank_panel_subset.csv, which ships with the
# repository. That CSV holds the DawnRank scores for the panel genes across the
# eight lung-metastasis specimens, extracted from Supplementary Table 5 of
# Siegel et al., J Clin Invest 2018;128(4):1371-1383 (doi:10.1172/JCI96153).
# The published workbook is copyright American Society for Clinical
# Investigation and is not redistributed here.
#
# This script is provenance only; the panel workflow never calls it. To rebuild
# or audit the shipped CSV, download the supplementary workbook and place it at
# data/raw/GSE110590/JCI96153.sdt1-8.xlsx:
#
#     curl -fsSL -o data/raw/GSE110590/JCI96153.sdt1-8.xlsx \
#       https://dm5migu4zj3pb.cloudfront.net/manuscripts/96000/96153/JCI96153.sdt1-8.xlsx
#
# Then:
#     Rscript scripts/derive_dawnrank_subset.R            # rewrite the CSV
#     Rscript scripts/derive_dawnrank_subset.R --verify   # compare, do not write

source(here::here("R", "paths.R"))
suppressPackageStartupMessages({ library(readxl); library(dplyr) })

VERIFY <- "--verify" %in% commandArgs(trailingOnly = TRUE)

SOURCE_WORKBOOK <- raw("GSE110590", "JCI96153.sdt1-8.xlsx")
SOURCE_SHEET    <- "Supplementary Table 5"
SOURCE_SHA256   <- "1039544685fe3ef35a10716287890192372740343f298ef7d935383863eb9add"
TARGET_CSV      <- derived("dawnrank_panel_subset.csv")

LUNG_MET_COLS <- c(
  "A2.LungMet", "A11.LUNG.MET", "A20.LungMet.L", "A20.LungMet.R",
  "A26.Lung.Met.1", "A26.Lung.Met.2", "A28.LungMet.LLL", "A28.LungMet.RUL"
)

if (!file.exists(SOURCE_WORKBOOK)) {
  stop("Source workbook not found: ", SOURCE_WORKBOOK, "\n",
       "  It is not redistributed with this repository. See the header of this\n",
       "  script for the download command, or DATA.md for provenance.",
       call. = FALSE)
}

got_sha <- unname(tools::sha256sum(SOURCE_WORKBOOK))
if (!identical(tolower(got_sha), SOURCE_SHA256)) {
  message("[dawnrank] NOTE: workbook SHA-256 differs from the release used here.")
  message("  expected: ", SOURCE_SHA256)
  message("  found:    ", got_sha)
}

gl <- readLines(file.path(GENE_LISTS_DIR, "panel_G_dawnrank_35.txt"))
intrinsic_idx <- grep("^## Intrinsic", gl)[1]
niche_idx     <- grep("^## Niche", gl)[1]
stopifnot(!is.na(intrinsic_idx) & !is.na(niche_idx))

extract_block <- function(start, end_excl) {
  block <- trimws(gl[(start + 1):(end_excl - 1)])
  block[nzchar(block) & !startsWith(block, "#")]
}
panel_genes <- c(extract_block(intrinsic_idx, niche_idx),
                 extract_block(niche_idx, length(gl) + 1))

dawnrank_raw <- read_excel(SOURCE_WORKBOOK, sheet = SOURCE_SHEET, skip = 1)
names(dawnrank_raw)[1] <- "Gene"
dawnrank_raw$Gene <- sub("\\|.*$", "", dawnrank_raw$Gene)

missing_cols <- setdiff(LUNG_MET_COLS, colnames(dawnrank_raw))
if (length(missing_cols))
  stop("Missing lung-met columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)

subset_tbl <- dawnrank_raw %>%
  filter(Gene %in% panel_genes) %>%
  select(Gene, all_of(LUNG_MET_COLS))
found <- intersect(panel_genes, subset_tbl$Gene)
subset_tbl <- as.data.frame(subset_tbl[match(found, subset_tbl$Gene), ])

not_found <- setdiff(panel_genes, found)
message("[dawnrank] panel genes: ", length(panel_genes),
        "; present in Supplementary Table 5: ", nrow(subset_tbl))
if (length(not_found))
  message("[dawnrank] absent: ", paste(not_found, collapse = ", "))

if (VERIFY) {
  if (!file.exists(TARGET_CSV))
    stop("No shipped CSV to verify against: ", TARGET_CSV, call. = FALSE)
  shipped <- read.csv(TARGET_CSV, check.names = FALSE, stringsAsFactors = FALSE)
  same_genes <- identical(shipped$Gene, subset_tbl$Gene)
  a <- as.matrix(shipped[, LUNG_MET_COLS]); storage.mode(a) <- "numeric"
  b <- as.matrix(subset_tbl[, LUNG_MET_COLS]); storage.mode(b) <- "numeric"
  same_values <- isTRUE(all.equal(a, b, tolerance = 0))
  if (same_genes && same_values) {
    message("[dawnrank] VERIFY OK: the shipped CSV matches a fresh rebuild.")
  } else {
    stop("[dawnrank] VERIFY FAILED: gene order match = ", same_genes,
         "; value match = ", same_values, call. = FALSE)
  }
} else {
  write.csv(subset_tbl, TARGET_CSV, row.names = FALSE)
  message("[dawnrank] wrote ", TARGET_CSV)
  message("[dawnrank] Next: re-run Panel J and confirm ",
          "output/tables/Fig1G_DawnRank_correlation_matrix.csv is unchanged.")
}
