#!/usr/bin/env Rscript
# =============================================================================
# AUCell Hallmark pathway enrichment in THREE self-contained analyses, each vs
# its OWN NT*NT (avoids any cross-10x-lane effect):
#   1) In vitro           (one HTO object)      vs in-vitro NT*NT
#   2) Synergistic         (in-vivo Synergy lane) vs Synergy-lane NT*NT
#   3) Buffering           (in-vivo Buffering lane) vs Buffering-lane NT*NT
# Per perturbation x pathway: delta AUCell (mean pert - mean NT*NT), Cohen d, BH p.
# AUCell (AUCell_buildRankings + calcAUC) on the 50 Hallmark sets, RNA counts.
#
# v7 repo overhaul 2026-07 : relocated from XL_code/fig3_panels_new_build/pathway_by_context.R;
#   this is a cache builder (no figure panel) -> group_cache/aucell_hallmark_{invitro,invivo}.rds
#   become cache(aucell_hallmark_{invitro,invivo}.rds); the hard-coded hallmark_common.rds read
#   is replaced by the version-locked load_hallmark_genesets() from 00_setup.R.
# =============================================================================
if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
set.seed(1234)
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(AUCell); library(matrixStats) })
.hall_df <- load_hallmark_genesets()
HALL <- split(.hall_df$gene_symbol, .hall_df$gs_name)   # 50 named HALLMARK_* gene-symbol vectors

score_aucell <- function(file, tag, lanecol){
  f <- cache(paste0("aucell_hallmark_", tag, ".rds"))
  if (file.exists(f)) { cat("[",tag,"] cached\n"); return(readRDS(f)) }
  cat("[",tag,"] scoring AUCell ...\n"); o <- readRDS(raw(file))
  if (length(SeuratObject::Layers(o[["RNA"]]))>1) o <- JoinLayers(o, assay="RNA")
  expr <- GetAssayData(o, assay="RNA", layer="counts")
  rk <- AUCell_buildRankings(expr, plotStats=FALSE, verbose=FALSE)
  au <- getAUC(AUCell_calcAUC(HALL, rk, aucMaxRank=ceiling(0.05*nrow(rk)), verbose=FALSE))  # 50 x cells
  res <- list(auc=au, pert=as.character(o$perturbation),
              lane=if(!is.na(lanecol)) as.character(o@meta.data[[lanecol]]) else rep("all", ncol(o)))
  saveRDS(res, f); rm(o,expr,rk); gc(verbose=FALSE); res
}
IV  <- score_aucell("Perturb_InVitro_HTO.rds",   "invitro", NA)
IVV <- score_aucell("Perturb_InVivo_merged.rds", "invivo",  "lane")

cohen <- function(a,b) (mean(a)-mean(b))/sqrt(((length(a)-1)*var(a)+(length(b)-1)*var(b))/(length(a)+length(b)-2))
CTX <- list(
  "In vitro"    = list(d=IV,  lane="all",       perts=c("Pten*NT","NT*Cdh1","Pten*Cdh1","NT*Cx3cl1","Pten*Cx3cl1","NT*Cxcr5","Pten*Cxcr5","NT*Tlr7","Pten*Tlr7")),
  "Synergistic" = list(d=IVV, lane="Synergy",   perts=c("Pten*NT","NT*Cdh1","Pten*Cdh1","NT*Cx3cl1","Pten*Cx3cl1")),
  "Buffering"   = list(d=IVV, lane="Buffering", perts=c("Pten*NT","NT*Cxcr5","Pten*Cxcr5","NT*Tlr7","Pten*Tlr7")))

rows <- list()
for (cx in names(CTX)) {
  C <- CTX[[cx]]; au <- C$d$auc; pert <- C$d$pert; lane <- C$d$lane
  inlane <- lane == C$lane
  ctrl <- which(inlane & pert=="NT*NT")
  for (p in C$perts) {
    pc <- which(inlane & pert==p); if (length(pc)<15) next
    delta <- rowMeans(au[,pc,drop=FALSE]) - rowMeans(au[,ctrl,drop=FALSE])
    d  <- apply(au,1,function(r) cohen(r[pc], r[ctrl]))
    pv <- apply(au,1,function(r) wilcox.test(r[pc], r[ctrl])$p.value); pa <- p.adjust(pv,"BH")
    rows[[paste(cx,p)]] <- data.frame(context=cx, perturbation=p, n_pert=length(pc), n_ctrl=length(ctrl),
      pathway=sub("HALLMARK_","",names(delta)), delta=delta, cohen_d=d, padj=pa, row.names=NULL)
  }
}
EN <- do.call(rbind, rows)
write.csv(EN, tbl("pathway_by_context_enrichment.csv"), row.names=FALSE)

cat("\n===== AUCell Hallmark enrichment vs each context's OWN NT*NT =====\n")
cat("# pathways enriched (padj<0.05 & |Cohen d|>0.2) per perturbation:\n\n")
sig <- EN[EN$padj<0.05 & abs(EN$cohen_d)>0.2,]
tab <- aggregate(pathway~context+perturbation, sig, length); colnames(tab)[3] <- "n_sig_pathways"
tab <- tab[order(match(tab$context,names(CTX)), -tab$n_sig_pathways),]
print(tab, row.names=FALSE)
cat("\n--- top enriched pathway per perturbation (by |Cohen d|) ---\n")
for (cx in names(CTX)) for (p in CTX[[cx]]$perts) {
  s <- sig[sig$context==cx & sig$perturbation==p,]; if(!nrow(s)) next
  s <- s[order(-abs(s$cohen_d)),]
  cat(sprintf("[%-11s] %-12s: %s\n", cx, p, paste(sprintf("%s(%+.2f)", s$pathway[1:min(4,nrow(s))], s$cohen_d[1:min(4,nrow(s))]), collapse=", ")))
}
