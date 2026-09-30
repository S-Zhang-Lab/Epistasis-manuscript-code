#!/usr/bin/env Rscript
# =============================================================================
# Standard group-level perturb-seq analysis BY PERTURBATION LABEL (lane-agnostic;
# lanes pooled, modest batch treated as noise per the user's direction).
#   (1) DEGs: each perturbation vs NT*NT (all NT*NT pooled)  + synergy-group vs
#       buffering-group (aggressive vs indolent).
#   (2) Hallmark UCell (50) per group: means + enrichment (Wilcoxon vs NT*NT, Cohen d).
# In vitro + in vivo. Caches per-cell UCell + writes tidy tables for figures.
# v7 repo overhaul 2026-07-15: relocated from XL_code/fig3_panels_new_build/group_DE_pathway.R;
#   cache builder for per-cell UCell Hallmark scores (ucell_hallmark_invitro.rds /
#   ucell_hallmark_invivo.rds via cache(), basenames preserved) — reads both raw
#   objects; Hallmark set now from load_hallmark_genesets(); DEG-count / top-gene /
#   hallmark-by-perturbation / pathway-enrichment tables dropped (superseded by
#   dedicated panel scripts).
# =============================================================================
if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
set.seed(1234)
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(presto); library(UCell); library(matrixStats) })
HALL <- load_hallmark_genesets(); HALL <- split(HALL$gene_symbol, HALL$gs_name)

run <- function(file, tag){
  cat("\n#########", tag, "#########\n")
  o <- readRDS(raw(file)); D <- GetAssayData(o, assay="SCT", layer="data")
  pert <- as.character(o$perturbation); names(pert) <- colnames(o)
  ## Hallmark UCell per cell
  U <- t(ScoreSignatures_UCell(D, features=HALL, name="", ncores=1))     # pathways x cells
  saveRDS(list(U=U, pert=pert), cache(paste0("ucell_hallmark_", tolower(gsub(" ","",tag)),".rds")))
  rm(o,D,U); gc(verbose=FALSE)
}
run("Perturb_InVitro_HTO.rds","In vitro"); run("Perturb_InVivo_merged.rds","In vivo")
cat("\nucell hallmark cache done.\n")
