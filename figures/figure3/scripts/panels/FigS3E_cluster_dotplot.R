if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
set.seed(1234)
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(ggplot2) })

U <- readRDS(cache("ucell_scores.rds"))
short <- function(v) sub("HALLMARK_","",v)

build <- function(uc, file, tag){
  o  <- readRDS(raw(file))
  cl <- setNames(as.character(o$joint_clusters), colnames(o))
  uc <- uc[names(cl)[names(cl) %in% rownames(uc)], , drop=FALSE]
  clv <- cl[rownames(uc)]
  colnames(uc) <- short(colnames(uc))
  clord <- as.character(sort(as.integer(unique(clv))))

  mu <- t(apply(t(uc), 1, function(x) tapply(x, clv, mean)))
  z  <- t(scale(t(mu)))[, clord, drop=FALSE]

  med <- apply(uc, 2, median)
  hi  <- sweep(uc, 2, med, ">")
  pct <- t(apply(t(hi), 1, function(x) tapply(x, clv, mean)))[, clord, drop=FALSE] * 100

  vv   <- apply(z,1,function(x) max(x)-min(x)); topv <- names(sort(vv,decreasing=TRUE))[1:22]
  zt <- z[topv, ]; pt <- pct[topv, ]
  featord <- rownames(zt)[order(apply(zt,1,which.max))]

  dd <- expand.grid(program=topv, cluster=clord, KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
  dd$z   <- as.vector(zt); dd$pct <- as.vector(pt)
  dd$program <- factor(dd$program, levels=featord)
  dd$cluster <- factor(dd$cluster, levels=rev(clord))

  p <- ggplot(dd, aes(program, cluster)) +
    geom_point(aes(size=pct, colour=z)) +
    scale_colour_gradient2(low="#2166ac", mid="grey92", high="#b2182b", midpoint=0,
                           name="mean\nz-score",
                           guide=guide_colourbar(raster=FALSE)) +
    scale_size_continuous(range=c(0.2,2.6), limits=c(0,100),
                          breaks=c(25,50,75,100), name="% cells\n> median") +
    labs(x=NULL, y="cluster",
         title=paste0("Cluster identity: ", tag, " (UCell Hallmark)")) +
    fig_theme +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=FS-0.5),
          axis.line=element_blank(), axis.ticks=element_line(linewidth=0.25),
          panel.grid.major=element_line(linewidth=0.15, colour="grey90"),
          legend.key.size=unit(7,"pt"))
  save_panel(p, paste0("FigS3E_cluster_identity_dotplot_", gsub("_","",tag)), 4.6, 2.9, fig_supp())
  cat(sprintf("dotplot %s written: %d programs x %d clusters\n", tag, length(topv), length(clord)))
  rm(o); gc(verbose=FALSE)
}
build(U$vitro, "Perturb_InVitro_HTO.rds", "in_vitro")
build(U$vivo,  "Perturb_InVivo_merged.rds","in_vivo")
cat("\ncluster identity dotplots done.\n")
