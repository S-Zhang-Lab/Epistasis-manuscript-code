suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
})

source(here::here("R", "paths.R"))
source(here::here("R", "helpers.R"))

GENE_SUMMARY   <- derived("R2_Day18.gene_summary.txt")
LEGACY_SUMMARY <- derived("R2_Day18_legacy.gene_summary.txt")
REFERENCE_WB   <- derived("R2_GI_scores.xlsx")

WB_TOL <- 1e-04

is_nt <- function(x) !is.na(x) & tolower(x) %in% c("nt", "ntc")

derive_r2_gi <- function(path) {
  raw <- read.delim(path, header = TRUE, sep = "\t",
                    check.names = FALSE, quote = "", comment.char = "")

  required <- c("id", "num", "neg|p-value", "neg|lfc", "pos|p-value", "pos|lfc")
  miss <- setdiff(required, names(raw))
  if (length(miss))
    stop("gene summary is missing required columns: ", paste(miss, collapse = ", "))
  if (anyDuplicated(tolower(raw$id)))
    stop("the gene summary contains duplicated pair ids (case-insensitive)")

  parts <- strsplit(raw$id, "*", fixed = TRUE)
  raw$Gene1 <- vapply(parts, `[`, character(1), 1)
  raw$Gene2 <- vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_,
                      character(1))

  id_index <- setNames(seq_len(nrow(raw)), tolower(raw$id))
  control_index <- function(gene) {
    if (is.na(gene) || is_nt(gene)) return(NA_integer_)
    hit <- id_index[tolower(c(paste0(gene, "*Nt"), paste0("Nt*", gene)))]
    hit <- hit[!is.na(hit)]
    if (!length(hit)) return(NA_integer_)
    as.integer(hit[1])
  }
  c1 <- vapply(raw$Gene1, control_index, integer(1))
  c2 <- vapply(raw$Gene2, control_index, integer(1))

  orphans <- unique(c(raw$Gene1[!is_nt(raw$Gene1) & !is.na(raw$Gene1) & is.na(c1)],
                      raw$Gene2[!is_nt(raw$Gene2) & !is.na(raw$Gene2) & is.na(c2)]))
  if (length(orphans))
    stop("no Gene*Nt or Nt*Gene control found for: ", paste(orphans, collapse = ", "))

  pick <- function(i, col) {
    o <- rep(NA_real_, length(i)); k <- !is.na(i)
    o[k] <- as.numeric(raw[[col]][i[k]]); o
  }
  pick_chr <- function(i, col) {
    o <- rep(NA_character_, length(i)); k <- !is.na(i)
    o[k] <- as.character(raw[[col]][i[k]]); o
  }

  g1 <- pick(c1, "neg|lfc"); g2 <- pick(c2, "neg|lfc")

  expected <- ifelse(is.na(g1), 0, g1) + ifelse(is.na(g2), 0, g2)
  gi_value <- as.numeric(raw$`neg|lfc`) - expected

  include <- !is_nt(raw$Gene1) & !is_nt(raw$Gene2) & !is.na(raw$Gene2)
  gi_mean <- mean(gi_value[include], na.rm = TRUE)
  gi_sd   <- sd(gi_value[include], na.rm = TRUE)
  if (!is.finite(gi_sd) || gi_sd == 0)
    stop("GI standard deviation is zero or not finite")

  out <- data.frame(
    Gene = raw$id, Count = as.integer(raw$num),
    LFC_Mean = as.numeric(raw$`neg|lfc`),
    Gene1 = raw$Gene1, Gene2 = raw$Gene2,
    Gene1_NT = pick_chr(c1, "id"), Gene1_NT_Count = pick(c1, "num"),
    Gene1_NT_LFC_Mean = g1,
    Gene2_NT = pick_chr(c2, "id"), Gene2_NT_Count = pick(c2, "num"),
    Gene2_NT_LFC_Mean = g2,
    Expected = expected, GI = gi_value,
    GI_Z_Include = as.integer(include),
    GI_Z_Score = ifelse(include, (gi_value - gi_mean) / gi_sd, NA_real_),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  out <- out[order(out$GI), ]
  rownames(out) <- NULL
  attr(out, "gi_mean") <- gi_mean
  attr(out, "gi_sd")   <- gi_sd
  out
}

message("--- Deriving the R2 GI table from its MAGeCK summary ---")

gi_tab <- derive_r2_gi(GENE_SUMMARY)
write.csv(gi_tab, tbl("R2_GI_derived.csv"), row.names = FALSE, na = "")
message(sprintf("  %s: %d rows, %d true gene pairs | GI mean %.6f, SD %.6f",
                basename(GENE_SUMMARY), nrow(gi_tab), sum(gi_tab$GI_Z_Include == 1),
                attr(gi_tab, "gi_mean"), attr(gi_tab, "gi_sd")))

if (file.exists(LEGACY_SUMMARY) && file.exists(REFERENCE_WB)) {
  legacy <- derive_r2_gi(LEGACY_SUMMARY)
  wb <- as.data.frame(read_excel(REFERENCE_WB, sheet = "GI"))
  cmp <- merge(
    transform(legacy, .k = toupper(Gene))[, c(".k", "LFC_Mean", "Expected", "GI")],
    transform(wb,     .k = toupper(Gene))[, c(".k", "LFC_Mean", "Expected", "GI")],
    by = ".k", suffixes = c("_new", "_wb")
  )
  devs <- vapply(c("LFC_Mean", "Expected", "GI"), function(v) {
    a <- cmp[[paste0(v, "_new")]]; b <- cmp[[paste0(v, "_wb")]]
    ok <- is.finite(a) & is.finite(b)
    if (!any(ok)) return(NA_real_)
    max(abs(a[ok] - b[ok]))
  }, numeric(1))
  if (any(devs > WB_TOL, na.rm = TRUE))
    stop(sprintf("construction no longer reproduces the legacy workbook (max |delta| %.3g)",
                 max(devs, na.rm = TRUE)))
  message(sprintf("  regression guard: legacy input still reproduces the legacy workbook (max |delta| %.2g)",
                  max(devs, na.rm = TRUE)))
} else {
  warning("legacy pair absent; construction not regression-checked")
}

log_session("04_r2_derive_gi")
