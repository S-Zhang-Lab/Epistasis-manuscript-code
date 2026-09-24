if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
set.seed(1234)
suppressPackageStartupMessages({ library(Seurat); library(data.table); library(ggplot2) })

GENO <- c("Pten*NT","NT*Cdh1","Pten*Cdh1","NT*Cx3cl1","Pten*Cx3cl1","NT*Cxcr5","Pten*Cxcr5","NT*Tlr7","Pten*Tlr7")
catg <- function(g) fifelse(g=="Pten*NT","PTEN single", fifelse(startsWith(g,"NT*"),"partner single","DKO"))
NPC <- 30; B <- 10000; RN <- 400

ediD <- function(D, ia, ib) 2*mean(D[ia,ib]) - mean(D[ia,ia]) - mean(D[ib,ib])
analyze <- function(file, grpcol, ctx){
  o <- readRDS(raw(file)); pc <- Embeddings(o,"pca"); if (ncol(pc)>NPC) pc <- pc[,1:NPC]
  md <- data.table(cell=colnames(o), pert=as.character(o$perturbation), grp=as.character(o@meta.data[[grpcol]])); rm(o); gc(verbose=FALSE)
  cellsOf <- function(g) md[pert==g, cell]; grpsOf <- function(g) unique(md[pert==g, grp])
  refOf   <- function(g) md[pert=="NT*NT" & grp %in% grpsOf(g), cell]
  sizes <- sapply(GENO, function(g) length(cellsOf(g)))
  ntg <- split(md[pert=="NT*NT", cell], md[pert=="NT*NT", grp])
  n0 <- min(min(sizes), floor(min(sapply(ntg, length))/2)); ntg <- ntg[sapply(ntg, length) >= 2*n0]

  cn <- unlist(lapply(ntg, function(cc) replicate(RN, { s <- sample(cc, 2*n0); Dp <- as.matrix(dist(pc[s,,drop=FALSE])); ediD(Dp, 1:n0, (n0+1):(2*n0)) })))
  thr95 <- as.numeric(quantile(cn, 0.95))
  rows <- rbindlist(lapply(GENO, function(g){
    gc_ <- sample(cellsOf(g), n0); rc <- sample(refOf(g), n0)
    Dp <- as.matrix(dist(pc[c(gc_,rc),,drop=FALSE]))
    obs <- ediD(Dp, 1:n0, (n0+1):(2*n0))
    cnt <- sum(replicate(B, { idx <- sample(2*n0); ediD(Dp, idx[1:n0], idx[(n0+1):(2*n0)]) }) >= obs)
    data.table(context=ctx, geno=g, cat=catg(g), edist=obs, p=(1+cnt)/(B+1)) }))
  cat(sprintf("[%s] n0=%d control-95pct=%.3f\n", ctx, n0, thr95))
  list(rows=rows, thr95=thr95)
}
IV <- analyze("Perturb_InVitro_HTO.rds","pool","In vitro")
VV <- analyze("Perturb_InVivo_merged.rds","lane","In vivo")
D <- rbind(IV$rows, VV$rows); band <- max(IV$thr95, VV$thr95)
D[, context := factor(context, levels=c("In vitro","In vivo"))]
D[, q := p.adjust(p, "BH"), by=context]
fq <- function(q){ s <- ifelse(q>=0.001, sprintf("q=%.2f",q), sprintf("q=%.1e",q)); gsub("e-0","e-",s) }
D[, qlab := fq(q)]
D[, band := band]
fwrite(D, tbl("Fig3F_edist9_stats.csv"))
fwrite(data.table(band = band), tbl("Fig3F_edist9_band.csv"))
print(D[order(context,-edist), .(context,geno,cat,edist=round(edist,3),p=signif(p,3),q=signif(q,3),qlab)])
cat(sprintf("persisted control-baseline band=%.3f -> Fig3F_edist9_band.csv\n", band))
