suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  library(SeuratObject); library(AUCell); library(msigdbr); library(dplyr)
})
set.seed(123L)

RDS_PTEN  <- local({
  a <- raw("RDS/Perturb_UTSW37_GSE322809_InVivo_merged.rds")
  b <- raw("RDS/Perturb_InVivo_merged.rds")
  if (file.exists(a)) a else b
})
RDS_CHEK2 <- raw("RDS/Perturb_UTSW33_GSE322808.rds")

mdf <- tryCatch(msigdbr(species = "Mus musculus", collection = "H"),
                error = function(e) msigdbr(species = "Mus musculus", category = "H"))
gs_col <- if ("gene_symbol" %in% names(mdf)) "gene_symbol" else "db_gene_symbol"
HALL <- lapply(split(mdf[[gs_col]], mdf$gs_name), unique)
message(sprintf("[recompute] Hallmark sets: %d", length(HALL)))

auc_matrix <- function(counts) {
  rk <- AUCell_buildRankings(as.matrix(counts), plotStats = FALSE, verbose = FALSE)
  A  <- getAUC(AUCell_calcAUC(HALL, rk, verbose = FALSE))
  rownames(A) <- paste0("HALLMARK-", gsub("_", "-", sub("^HALLMARK_", "", rownames(A))))
  A
}

gi_table <- function(A, cond, dual, a, b, ntnt = "NT*NT") {
  keep <- cond %in% c(dual, a, b, ntnt)
  A <- A[, keep, drop = FALSE]; cond <- droplevels(factor(cond[keep]))
  m <- sapply(levels(cond), function(g) rowMeans(A[, cond == g, drop = FALSE]))
  gi <- m[, dual] - m[, a] - m[, b] + m[, ntnt]
  tibble::tibble(Pathway = rownames(A),
                 dual = m[, dual], singleA = m[, a], singleB = m[, b], NTNT = m[, ntnt],
                 Expected = m[, a] + m[, b] - m[, ntnt], GI = gi)
}

message("[recompute] PTEN x CX3CL1  (UTSW37, Synergy lane) ...")
o <- readRDS(RDS_PTEN)
o[["RNA"]] <- suppressWarnings(SeuratObject::JoinLayers(o[["RNA"]]))
sel  <- o$lane == "Synergy" &
        o$crRNA_classification.global %in% c("Pten*Cx3cl1","Pten*NT","NT*Cx3cl1","NT*NT")
A    <- auc_matrix(SeuratObject::LayerData(o, assay = "RNA", layer = "counts")[, sel])
pten <- gi_table(A, o$crRNA_classification.global[sel], "Pten*Cx3cl1", "Pten*NT", "NT*Cx3cl1")
readr::write_csv(pten, derived("pathwayGI_PTEN_CX3CL1.csv"))
message(sprintf("    wrote pathwayGI_PTEN_CX3CL1.csv (GI range %.4f..%.4f)", min(pten$GI), max(pten$GI)))
rm(o, A); gc(verbose = FALSE)

message("[recompute] CHEK2 x CX3CL1  (UTSW33) ...")
o <- readRDS(RDS_CHEK2)
gp  <- o$Gene_Pair
sel <- gp %in% c("CHEK2*CX3CL1","CHEK2*NT","NT*CX3CL1","NT*NT","APC*CX3CL1","APC*NT")
A   <- auc_matrix(SeuratObject::LayerData(o, assay = "RNA", layer = "counts")[, sel])
gpk <- droplevels(factor(gp[sel]))

chek2 <- gi_table(A, gpk, "CHEK2*CX3CL1", "CHEK2*NT", "NT*CX3CL1")
readr::write_csv(chek2, derived("pathwayGI_CHEK2_CX3CL1.csv"))
message(sprintf("    wrote pathwayGI_CHEK2_CX3CL1.csv (GI range %.4f..%.4f)", min(chek2$GI), max(chek2$GI)))

cat(sprintf("\n[recompute] cross-pair cor(PTEN, CHEK2) = %.3f  (Panel H)\n",
            cor(pten$GI, chek2$GI[match(pten$Pathway, chek2$Pathway)])))
if (!interactive()) tryCatch(log_session("recompute_pathway_GI_from_RDS"), error = function(e) NULL)
