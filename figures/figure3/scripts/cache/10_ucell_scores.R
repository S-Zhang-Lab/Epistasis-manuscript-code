#!/usr/bin/env Rscript
# =============================================================================
# Figure 3 — SINGLE-CELL epistasis (cells = biological replicates)
# Data-first: no presupposed convergence/divergence. Each context analyzed on
# its own; results laid out for side-by-side in-vitro | in-vivo reading.
#
# Method (per S. Zhang directives):
#   1. Per-cell UCell score for each of the 50 mouse Hallmark programs.
#   2. CELL-LEVEL epistasis: within each partner's pool(vitro)/lane(vivo), fit
#      score ~ Pten * gene on per-cell UCell scores; the Pten:gene interaction
#      term IS the non-additive epistatic residual, with effect size + p
#      (cells are replicates). Controls come from the same pool/lane.
#   3. Cluster/state composition per perturbation + compositional epistasis.
# Pseudobulk is NOT used. Seed 1234.
#
# v7 repo overhaul 2026-07-15: relocated from
#   XL_code/fig3_singlecell_epistasis/sc_epistasis.R; cache builder (no panel) —
#   heavy caches objects/{hallmark_common,ucell_scores}.rds -> cache() and
#   tables/*.csv -> cache() (all basenames preserved); reads both raw objects;
#   live msigdbr() Hallmark build -> load_hallmark_genesets() (frozen snapshot).
# =============================================================================
if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
set.seed(1234)
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(UCell); library(msigdbr); library(Matrix)
})
say <- function(...) cat(sprintf("[%s] ",format(Sys.time(),"%H:%M:%S")), ..., "\n")

PARTNERS <- list(Cdh1="Synergy", Cx3cl1="Synergy", Cxcr5="Buffering", Tlr7="Buffering")

# ---- Hallmark, restricted to genes present in BOTH contexts ------------------
say("loading objects")
vivo  <- readRDS(raw("Perturb_InVivo_merged.rds"))
vitro <- readRDS(raw("Perturb_InVitro_HTO.rds"))
g_both <- intersect(rownames(vivo[["RNA"]]), rownames(vitro[["RNA"]]))
say(sprintf("common RNA genes = %d", length(g_both)))
mdf <- load_hallmark_genesets()
H <- lapply(split(mdf$gene_symbol, mdf$gs_name), unique)
H <- lapply(H, function(s) intersect(s, g_both)); H <- H[sapply(H,length) >= 5]
say(sprintf("Hallmark sets = %d (scored on genes common to both contexts)", length(H)))
saveRDS(H, cache("hallmark_common.rds"))

# ---- per-cell UCell -----------------------------------------------------------
score_ctx <- function(obj, tag){
  say(sprintf("UCell scoring %s (%d cells)", tag, ncol(obj)))
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])                   # merged objects split counts across layers
  m <- LayerData(obj, assay="RNA", layer="counts")
  u <- ScoreSignatures_UCell(m, features=H, maxRank=1500,
                             BPPARAM=BiocParallel::SerialParam(), name="")
  colnames(u) <- sub("_$","",colnames(u))
  u
}
U_vivo  <- score_ctx(vivo,  "in vivo")
U_vitro <- score_ctx(vitro, "in vitro")
saveRDS(list(vivo=U_vivo, vitro=U_vitro), cache("ucell_scores.rds"))

# ---- cell-level factorial epistasis: score ~ Pten * gene ---------------------
# returns per Hallmark: interaction coef (non-additive residual), se, t, p,
# plus the four condition mean UCell scores.
factorial_epistasis <- function(U, meta, grpcol, ctx){
  rows <- list()
  for (partner in names(PARTNERS)){
    poolval <- PARTNERS[[partner]]
    in_pool <- meta[[grpcol]] == poolval
    pert    <- meta$perturbation
    keep <- in_pool & pert %in% c("NT*NT","Pten*NT",paste0("NT*",partner),paste0("Pten*",partner))
    idx  <- which(keep)
    Pten <- grepl("Pten", pert[idx])
    gene <- grepl(partner, pert[idx])
    for (hm in colnames(U)){
      y  <- U[idx, hm]
      fit <- tryCatch(summary(lm(y ~ Pten * gene)), error=function(e) NULL)
      if (is.null(fit)) next
      co <- fit$coefficients
      itn <- grep(":", rownames(co))
      if (length(itn)!=1) next
      # condition means
      mm <- tapply(y, interaction(Pten,gene,drop=TRUE), mean)
      rows[[length(rows)+1]] <- data.frame(
        context=ctx, partner=partner, pool=poolval,
        pathway=sub("HALLMARK_","",hm),
        interaction=round(co[itn,1],4), se=round(co[itn,2],4),
        t=round(co[itn,3],3), p=signif(co[itn,4],3),
        mean_NTNT=round(mm["FALSE.FALSE"],4), mean_Pten=round(mm["TRUE.FALSE"],4),
        mean_gene=round(mm["FALSE.TRUE"],4), mean_double=round(mm["TRUE.TRUE"],4),
        n_double=sum(Pten&gene), row.names=NULL)
    }
  }
  res <- do.call(rbind, rows)
  res$padj <- p.adjust(res$p, "BH")
  res
}
say("cell-level factorial epistasis (in vivo)")
epi_vivo  <- factorial_epistasis(U_vivo,  vivo@meta.data,  "lane", "InVivo")
say("cell-level factorial epistasis (in vitro)")
epi_vitro <- factorial_epistasis(U_vitro, vitro@meta.data, "pool", "InVitro")
epi <- rbind(epi_vivo, epi_vitro)
write.csv(epi, cache("sc_epistasis_interaction.csv"), row.names=FALSE)

# ---- cluster composition per perturbation + compositional epistasis ----------
composition <- function(meta, ctx){
  tb <- as.data.frame.matrix(table(meta$perturbation, meta$joint_clusters))
  fr <- sweep(tb, 1, rowSums(tb), "/")                       # perturbation x cluster fractions
  write.csv(cbind(perturbation=rownames(fr), round(fr,4)),
            cache(sprintf("cluster_fractions_%s.csv", ctx)), row.names=FALSE)
  # compositional epistasis per partner: frac_double - [frac_Pten + frac_gene - frac_NTNT]
  out <- list()
  for (partner in names(PARTNERS)){
    dd<-paste0("Pten*",partner); gg<-paste0("NT*",partner)
    if (all(c("NT*NT","Pten*NT",gg,dd) %in% rownames(fr))){
      resid <- fr[dd,] - (fr["Pten*NT",] + fr[gg,] - fr["NT*NT",])
      out[[partner]] <- data.frame(context=ctx, partner=partner,
        cluster=colnames(fr), comp_epistasis=round(as.numeric(resid),4), row.names=NULL)
    }
  }
  do.call(rbind, out)
}
comp <- rbind(composition(vivo@meta.data,"InVivo"), composition(vitro@meta.data,"InVitro"))
write.csv(comp, cache("sc_compositional_epistasis.csv"), row.names=FALSE)

# ---- neutral console summary (top |interaction|, FDR<0.05) -------------------
cat("\n=========== TOP CELL-LEVEL EPISTATIC PROGRAMS (|interaction|, FDR<0.05) ===========\n")
for (ctx in c("InVitro","InVivo")) for (partner in names(PARTNERS)){
  sub <- epi[epi$context==ctx & epi$partner==partner & epi$padj<0.05,]
  sub <- sub[order(-abs(sub$interaction)),]
  cat(sprintf("\n%-8s | %-7s (n_double=%s):\n", ctx, partner, if(nrow(sub)) sub$n_double[1] else NA))
  if (nrow(sub)) for (i in seq_len(min(5,nrow(sub))))
    cat(sprintf("   %+0.3f  %-32s  p=%.1e\n", sub$interaction[i], sub$pathway[i], sub$p[i]))
  else cat("   (none FDR<0.05)\n")
}
say("done. tables in", cache())
