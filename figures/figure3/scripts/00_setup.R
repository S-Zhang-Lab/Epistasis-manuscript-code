.required_pkgs <- c("tidyverse", "Matrix", "ggrepel")

.optional_pkgs <- c("Seurat", "ComplexHeatmap", "circlize", "fgsea",
                    "msigdbr", "ggExtra", "ggalt", "patchwork")

suppressPackageStartupMessages({
  for (p in .required_pkgs) library(p, character.only = TRUE)
  for (p in .optional_pkgs) {
    if (requireNamespace(p, quietly = TRUE)) {
      library(p, character.only = TRUE)
    } else {
      message(sprintf("[00_setup.R] optional package not installed: %s (skipping)", p))
    }
  }
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

.get_setup_script_dir <- function() {
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

source(here::here("R", "paths.R"), local = FALSE)

RAW_DIR       <- raw()
DERIVED_DIR   <- derived()
TABLES_DIR    <- tbl()
INTERMED_DIR  <- interm()
MAIN_DIR      <- fig_main()
SUPP_DIR      <- fig_supp()
LOG_DIR       <- log_path()

INPUT_DIR     <- raw()

ensure_output_dirs()
dir.create(interm("legacy_rerun"), recursive = TRUE, showWarnings = FALSE)

GLOBAL_SEED          <- 42L
ROW_KMEANS_SEED      <- 123L
set.seed(GLOBAL_SEED)

PARTNERS_SYNERGY <- c("Cdh1", "Cx3cl1")
PARTNERS_BUFFER  <- c("Cxcr5", "Tlr7")
PARTNERS_ALL     <- c(PARTNERS_SYNERGY, PARTNERS_BUFFER)

INVIVO_RDS  <- raw("Perturb_InVivo_merged.rds")
INVITRO_RDS <- raw("Perturb_InVitro_HTO.rds")

GI_BLUE        <- "#6FA3D2"
GI_RED         <- "#B2182B"
GI_WHITE       <- "#FFFFFF"

DNES_LIMITS_SYMMETRIC      <- c(-2, 2)
SIGNATURE_LIMITS_SYMMETRIC <- c(-2, 2)

GI_Z_LIMITS_ASYMMETRIC <- c(-3.5, 1.4)

PALETTE_CONTEXT <- c(
  InVitro            = "#999999",
  InVivoBuffer       = "#2166AC",
  InVivoSynergistic  = "#B2182B"
)

PALETTE_PHENO <- c(
  Synergistic = "#B2182B",
  Buffer      = "#2166AC"
)

PALETTE_ANNOT <- c(
  InVivo_enriched  = "#7B3294",
  InVitro_enriched = "#008837"
)

source(file.path(PROJECT_ROOT, "R", "palettes.R"),      local = FALSE)
source(file.path(PROJECT_ROOT, "R", "seurat_compat.R"), local = FALSE)

source(file.path(PROJECT_ROOT, "R", "fig3_theme.R"),    local = FALSE)

HALLMARK_FROZEN_PATH <- derived("hallmark_pathways_v26.1.0_MM.csv")

load_hallmark_genesets <- function(force_refresh = FALSE) {
  if (force_refresh) {
    if (!requireNamespace("msigdbr", quietly = TRUE)) {
      stop("load_hallmark_genesets(force_refresh = TRUE) requires `msigdbr`.")
    }
    df <- msigdbr::msigdbr(species = "mouse",
                           db_species = "MM",
                           collection = "MH")
    out <- unique(data.frame(gs_name     = df$gs_name,
                             gene_symbol = df$gene_symbol,
                             stringsAsFactors = FALSE))
    return(out[order(out$gs_name, out$gene_symbol), ])
  }
  if (!file.exists(HALLMARK_FROZEN_PATH)) {
    stop("Frozen Hallmark snapshot missing at: ", HALLMARK_FROZEN_PATH,
         "\nRegenerate with `load_hallmark_genesets(force_refresh = TRUE)` ",
         "and commit the result.")
  }

  df <- utils::read.csv(HALLMARK_FROZEN_PATH, comment.char = "#",
                        stringsAsFactors = FALSE)
  if (!all(c("pathway", "gene_symbol") %in% colnames(df))) {
    stop("Frozen Hallmark snapshot has unexpected columns: ",
         paste(colnames(df), collapse = ", "))
  }
  data.frame(gs_name     = df$pathway,
             gene_symbol = df$gene_symbol,
             stringsAsFactors = FALSE)
}

invisible(NULL)
