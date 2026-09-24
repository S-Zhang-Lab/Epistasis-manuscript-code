if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
set.seed(1234)
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(presto); library(UCell); library(matrixStats) })
HALL <- load_hallmark_genesets(); HALL <- split(HALL$gene_symbol, HALL$gs_name)

run <- function(file, tag){
  cat("\n#########", tag, "#########\n")
  o <- readRDS(raw(file)); D <- GetAssayData(o, assay="SCT", layer="data")
  pert <- as.character(o$perturbation); names(pert) <- colnames(o)

  U <- t(ScoreSignatures_UCell(D, features=HALL, name="", BPPARAM=BiocParallel::SerialParam()))
  saveRDS(list(U=U, pert=pert), cache(paste0("ucell_hallmark_", tolower(gsub(" ","",tag)),".rds")))
  rm(o,D,U); gc(verbose=FALSE)
}
run("Perturb_InVitro_HTO.rds","In vitro"); run("Perturb_InVivo_merged.rds","In vivo")
cat("\nucell hallmark cache done.\n")
