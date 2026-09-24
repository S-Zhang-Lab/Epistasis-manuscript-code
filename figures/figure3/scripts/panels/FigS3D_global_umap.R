if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
set.seed(1234)
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(Matrix)
  library(ggplot2); library(patchwork); library(scales) })
say <- function(...) cat(sprintf("[%s] ",format(Sys.time(),"%H:%M:%S")), ..., "\n")
clpal <- function(n) setNames(hue_pal(l=58, c=105)(n), as.character(0:(n-1)))
PERTPAL <- c("NT*NT"="#9aa0a6","Pten*NT"="#111111",
             "NT*Cdh1"="#f4a3a0","Pten*Cdh1"="#d6604d","NT*Cx3cl1"="#f6b56b","Pten*Cx3cl1"="#e08214",
             "NT*Cxcr5"="#a6cee3","Pten*Cxcr5"="#2166ac","NT*Tlr7"="#b7e4c7","Pten*Tlr7"="#1b7837")
GRPPAL <- c(Synergy="#B2182B", Buffering="#2166AC")

build_ctx <- function(file, tag, groupcol){
  say(paste("loading", tag)); o <- readRDS(raw(file)); DefaultAssay(o) <- "RNA"
  cl <- sort(as.integer(unique(as.character(o$joint_clusters)))); nclu <- length(cl)
  df <- data.frame(UMAP1=Embeddings(o,"umap")[,1], UMAP2=Embeddings(o,"umap")[,2],
                   cluster=factor(as.character(o$joint_clusters), levels=as.character(cl)),
                   perturbation=factor(o$perturbation, levels=names(PERTPAL)),
                   grp=o@meta.data[[groupcol]],
                   nFeature=o$nFeature_RNA, nCount=o$nCount_RNA, mt=o$percent.mt)
  pal <- clpal(nclu)

  cen <- aggregate(cbind(UMAP1,UMAP2)~cluster, df, median)
  p_cl <- umap_scaffold(df) + pts("cluster") + scale_colour_manual(values=pal, guide="none") +
    geom_text(data=cen, aes(UMAP1,UMAP2,label=cluster), size=FS/.MM, fontface="bold", colour="black") +
    labs(title=paste0(tag, "  |  ", nclu, " clusters"))
  save_panel(p_cl, paste0("FigS3D_umap_clusters_", gsub("_", "", tag)), 2.05, 2.05, fig_supp())

  say(sprintf("%s: %d cells, %d clusters — panels written", tag, ncol(o), nclu))

}

invisible(build_ctx("Perturb_InVitro_HTO.rds", "in_vitro", "pool"))
invisible(build_ctx("Perturb_InVivo_merged.rds","in_vivo",  "lane"))
say("global QC/PCA-elbow/UMAP panels complete. (Condition-level PCA: build_condition_pca.R)")
