#!/usr/bin/env Rscript
# =============================================================================
# Fig 3 STEP 1 (rebuild): global QC / PCA / clustering panels + condition-level PCA.
# In vitro & in vivo analyzed/clustered SEPARATELY. Convention: fig3_theme.R
# (120x120pt square UMAPs, circles, ~5pt Helvetica, editable). Scripts live with panels.
# v7 repo overhaul 2026-07-15: relocated from XL_code/fig3_panels_new_build/build_global_panels.R;
#   was umap_clusters_in_vitro / umap_clusters_in_vivo -> FigS3D_umap_clusters_invitro / FigS3D_umap_clusters_invivo (fig_supp).
#   Supplementary Fig. S3D keeps ONLY the two global UMAP-by-cluster panels; the umap-by-perturbation,
#   umap-by-pool/lane, QC violins, PCA elbow, and condition-PCA pseudobulk outputs are dropped
#   (covered elsewhere) and commented out below.
# =============================================================================
if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
set.seed(1234)
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(Matrix)
  library(ggplot2); library(patchwork); library(scales); library(readr); library(dplyr) })
say <- function(...) cat(sprintf("[%s] ",format(Sys.time(),"%H:%M:%S")), ..., "\n")
clpal <- function(n) setNames(hue_pal(l=58, c=105)(n), as.character(0:(n-1)))
PERTPAL <- c("NT*NT"="#9aa0a6","Pten*NT"="#111111",
             "NT*Cdh1"="#f4a3a0","Pten*Cdh1"="#d6604d","NT*Cx3cl1"="#f6b56b","Pten*Cx3cl1"="#e08214",
             "NT*Cxcr5"="#a6cee3","Pten*Cxcr5"="#2166ac","NT*Tlr7"="#b7e4c7","Pten*Tlr7"="#1b7837")
GRPPAL <- c(Synergy="#B2182B", Buffering="#2166AC")

# ---- per-context global panels ----------------------------------------------
build_ctx <- function(file, tag, groupcol){
  say(paste("loading", tag)); o <- readRDS(raw(file)); DefaultAssay(o) <- "RNA"
  cl <- sort(as.integer(unique(as.character(o$joint_clusters)))); nclu <- length(cl)
  emb <- Embeddings(o, "umap")
  df <- data.frame(cell_barcode=rownames(emb), UMAP1=emb[,1], UMAP2=emb[,2],
                   cluster=factor(as.character(o$joint_clusters), levels=as.character(cl)),
                   perturbation=factor(o$perturbation, levels=names(PERTPAL)),
                   grp=o@meta.data[[groupcol]],
                   nFeature=o$nFeature_RNA, nCount=o$nCount_RNA, mt=o$percent.mt)
  pal <- clpal(nclu)

  # UMAP by cluster (on-data labels; MAIN) -- 120pt square
  cen <- aggregate(cbind(UMAP1,UMAP2)~cluster, df, median)
  p_cl <- umap_scaffold(df) + pts("cluster") + scale_colour_manual(values=pal, guide="none") +
    geom_text(data=cen, aes(UMAP1,UMAP2,label=cluster), size=FS/.MM, fontface="bold", colour="black") +
    labs(title=paste0(tag, "  |  ", nclu, " clusters"))
  save_panel(p_cl, paste0("FigS3D_umap_clusters_", gsub("_", "", tag)), 2.05, 2.05, fig_supp())

  source_data <- df |>
    transmute(
      context = tag,
      cell_barcode,
      umap_1 = UMAP1,
      umap_2 = UMAP2,
      cluster = as.character(cluster),
      cluster_colour = unname(pal[as.character(cluster)])
    ) |>
    left_join(
      cen |>
        transmute(
          cluster = as.character(cluster),
          cluster_label_umap_1 = UMAP1,
          cluster_label_umap_2 = UMAP2
        ),
      by = "cluster"
    )

  # DROPPED (v7): UMAP by perturbation (supp) -- covered elsewhere.
  # p_pt <- umap_scaffold(df) + pts("perturbation", size=0.24) +
  #   scale_colour_manual(values=PERTPAL, name=NULL) +
  #   guides(colour=guide_legend(override.aes=list(size=1.1), ncol=1, keyheight=unit(5,"pt"))) +
  #   labs(title=paste0(tag, "  |  perturbation"))
  # save_panel(p_pt, paste0("umap_perturbation_", tag), 3.0, 2.05, fig_supp())

  # DROPPED (v7): UMAP by pool/lane (supp) -- covered elsewhere.
  # p_g <- umap_scaffold(df) + pts("grp", size=0.24) + scale_colour_manual(values=GRPPAL, name=NULL) +
  #   guides(colour=guide_legend(override.aes=list(size=1.2))) +
  #   labs(title=paste0(tag, "  |  ", groupcol))
  # save_panel(p_g, paste0("umap_", groupcol, "_", tag), 2.7, 2.05, fig_supp())

  # DROPPED (v7): QC violins by cluster (supp) -- covered elsewhere.
  # qv <- function(y, ylab) ggplot(df, aes(cluster, .data[[y]], fill=cluster)) +
  #   geom_violin(scale="width", linewidth=0.15, colour="grey35") +
  #   geom_boxplot(width=0.12, outlier.size=0.05, linewidth=0.15, fill="white", colour="grey20") +
  #   scale_fill_manual(values=pal, guide="none") + labs(x="cluster", y=ylab) + fig_theme
  # qc <- (qv("nFeature","genes/cell")/qv("nCount","UMIs/cell")/qv("mt","% mito")) +
  #   plot_annotation(title=paste0("QC by cluster — ", tag), theme=theme(plot.title=element_text(size=FS,face="bold",family="Helvetica")))
  # save_panel(qc, paste0("qc_violins_", tag), 4.0, 3.6, fig_supp())

  # DROPPED (v7): PCA elbow (supp) -- covered elsewhere.
  # sdv <- Stdev(o,reduction="pca")[1:min(50,ncol(Embeddings(o,"pca")))]
  # eb <- ggplot(data.frame(PC=seq_along(sdv), SD=sdv), aes(PC, SD)) +
  #   geom_vline(xintercept=30, linetype=2, colour="grey60", linewidth=0.3) +
  #   geom_point(size=0.5, colour="#2166ac") +
  #   labs(x="PC", y="SD", title=paste0("PCA elbow — ", tag)) + fig_theme
  # save_panel(eb, paste0("pca_elbow_", tag), 2.3, 1.8, fig_supp())

  say(sprintf("%s: %d cells, %d clusters — panels written", tag, ncol(o), nclu))
  # DROPPED (v7): pseudobulk return for the condition PCA -- covered by the condition-PCA script.
  # o[["RNA"]] <- JoinLayers(o[["RNA"]]); cts <- LayerData(o,assay="RNA",layer="counts")
  # pb <- sapply(split(colnames(o), o$perturbation), function(cc) Matrix::rowSums(cts[,cc,drop=FALSE]))
  # colnames(pb) <- paste(tag, colnames(pb), sep="|"); rm(o); gc(verbose=FALSE); pb
  rm(o, emb, df, p_cl); gc(verbose = FALSE)
  source_data
}

source_data <- bind_rows(
  build_ctx("Perturb_InVitro_HTO.rds", "in_vitro", "pool"),
  build_ctx("Perturb_InVivo_merged.rds", "in_vivo", "lane")
)
write_csv(source_data, tbl("FigS3D_global_umap_source_data.csv"))
say("global QC/PCA-elbow/UMAP panels complete. (Condition-level PCA: build_condition_pca.R)")
