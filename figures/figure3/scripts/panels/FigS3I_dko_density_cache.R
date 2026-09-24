if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
set.seed(1)
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(Matrix); library(matrixStats) })
K <- 40; NPG <- 150
PARTNERS <- c("Cdh1","Cx3cl1","Cxcr5","Tlr7")
POOL <- c(Cdh1="Synergy", Cx3cl1="Synergy", Cxcr5="Buffering", Tlr7="Buffering")
CTX  <- list("In vitro"=list(file="Perturb_InVitro_HTO.rds", col="pool"),
             "In vivo" =list(file="Perturb_InVivo_merged.rds", col="lane"))

dko_cvscore <- function(Dg, grpv, dko, K){
  allc <- colnames(Dg); grps <- unique(grpv)
  fold <- sample(rep(1:2, length.out=length(allc))); names(fold) <- allc
  sc <- setNames(rep(NA_real_, length(allc)), allc)
  for (f in 1:2){
    tr <- allc[fold==f]; te <- allc[fold!=f]; gtr <- grpv[tr]
    mu <- sapply(grps, function(g) rowMeans(Dg[, tr[gtr==g], drop=FALSE]))
    stat <- mu[,dko] - rowMeans(mu[, setdiff(grps,dko), drop=FALSE])
    sig  <- names(sort(stat, decreasing=TRUE))[1:K]
    Zte  <- Dg[sig, te, drop=FALSE]; Zte <- (Zte - rowMeans(Zte))/(rowSds(Zte)+1e-9)
    sc[te] <- colMeans(Zte)
  }
  sc
}

OUT <- list()
for (cx in names(CTX)){
  o <- readRDS(raw(CTX[[cx]]$file))
  pert <- as.character(o$perturbation); names(pert) <- colnames(o)
  meta <- as.character(o@meta.data[[CTX[[cx]]$col]]); names(meta) <- colnames(o)
  Dfull <- GetAssayData(o, assay="SCT", layer="data")
  for (pn in PARTNERS){
    grps <- c("NT*NT","Pten*NT",paste0("NT*",pn),paste0("Pten*",pn)); dko <- paste0("Pten*",pn)
    cells <- lapply(grps, function(g){ cc <- names(pert)[pert==g & meta==POOL[pn]]; sample(cc, min(NPG,length(cc))) })
    names(cells) <- grps
    allc <- unlist(cells, use.names=FALSE); grpv <- rep(grps, lengths(cells)); names(grpv) <- allc
    Dg <- as.matrix(Dfull[, allc, drop=FALSE]); Dg <- Dg[rowSums(Dg>0)>=5,,drop=FALSE]
    sc <- dko_cvscore(Dg, grpv, dko, K)
    v  <- rowVars(Dg); hvg <- rownames(Dg)[order(-v)][1:min(2000,sum(v>1e-8))]
    pc <- prcomp(t(Dg[hvg,]), center=TRUE, scale.=TRUE, rank.=2)$x
    OUT[[paste(cx,pn)]] <- data.frame(context=cx, partner=pn, genotype=grpv[allc],
                                      dko_score=sc[allc], PC1=pc[allc,1], PC2=pc[allc,2],
                                      row.names=NULL)
    cat(sprintf("[%s %s] cells=%d genes=%d\n", cx, pn, length(allc), nrow(Dg)))
  }
  rm(o, Dfull); gc(verbose=FALSE)
}
DF <- do.call(rbind, OUT)
write.csv(DF, cache("compression_concrete_cells.csv"), row.names=FALSE)
cat("\nwrote compression_concrete_cells.csv  (", nrow(DF), "cells )\n")
