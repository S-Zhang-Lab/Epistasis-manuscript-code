#!/usr/bin/env Rscript
# =============================================================================
# Fig. 3C — condition-level PCA.
# 8 points = 4 partners x 2 contexts (Pten x partner double-KO pseudobulk).
# shape = PTEN-loss partner; fill = grey (in vitro) / blue (in-vivo Buffer) /
# red (in-vivo Synergistic); in-vitro & in-vivo hulls; point labels; dashed
# origin; PC1/PC2 % axes; two legends. Base pdf(), Helvetica, ~5pt.
# v7 repo overhaul 2026-07 : relocated from XL_code/fig3_panels_new_build/build_condition_pca.R; was panelB_condition_pca -> Fig3C_condition_pca.
# =============================================================================
if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
set.seed(1234)
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(Matrix)
  library(ggplot2); library(ggforce); library(ggrepel) })
PARTNERS <- c("Cdh1","Cx3cl1","Cxcr5","Tlr7")

dko_pb <- function(file, ctx){          # pseudobulk (summed RNA counts) of each Pten x partner double-KO
  o <- readRDS(raw(file)); o[["RNA"]] <- JoinLayers(o[["RNA"]])
  cts <- LayerData(o, assay="RNA", layer="counts")
  out <- sapply(PARTNERS, function(p){
    cc <- colnames(o)[o$perturbation == paste0("Pten*",p)]; Matrix::rowSums(cts[,cc,drop=FALSE]) })
  colnames(out) <- paste(ctx, PARTNERS, sep="|"); rm(o); gc(verbose=FALSE); out
}
pv <- dko_pb("Perturb_InVitro_HTO.rds","in vitro")
po <- dko_pb("Perturb_InVivo_merged.rds","in vivo")
g  <- intersect(rownames(pv), rownames(po)); M <- cbind(pv[g,], po[g,])
lg <- log2(sweep(M,2,colSums(M),"/")*1e6 + 1)
vv <- apply(lg,1,var); top <- names(sort(vv,decreasing=TRUE))[1:2000]
pc <- prcomp(t(lg[top,]), scale.=TRUE); ve <- round(100*pc$sdev^2/sum(pc$sdev^2),1)

pd <- data.frame(PC1=pc$x[,1], PC2=pc$x[,2],
                 context=sub("\\|.*","",colnames(M)), partner=sub(".*\\|","",colnames(M)))
pd$partner <- factor(pd$partner, levels=PARTNERS)
pd$eclass  <- ifelse(pd$partner %in% c("Cdh1","Cx3cl1"), "Synergistic", "Buffer")
pd$label   <- ifelse(pd$context=="in vitro", paste0("InVitro\n",pd$partner),
                     paste0(pd$eclass,"\n",pd$partner))
vit <- subset(pd, context=="in vitro"); viv <- subset(pd, context=="in vivo")

SHP <- c(Cdh1=21, Cx3cl1=22, Cxcr5=24, Tlr7=23)
p <- ggplot(pd, aes(PC1, PC2)) +
  geom_mark_hull(aes(fill=context, group=context), concavity=2.6, expand=unit(3,"mm"),
                 radius=unit(3,"mm"), alpha=0.22, colour=NA) +
  scale_fill_manual(values=c("in vitro"="#bcd6ea","in vivo"="#eecccc"), guide="none") +
  geom_hline(yintercept=0, linetype=2, linewidth=0.25, colour="grey55") +
  geom_vline(xintercept=0, linetype=2, linewidth=0.25, colour="grey55") +
  ggnewscale::new_scale_fill() +
  geom_point(data=vit, aes(shape=partner), fill="#9e9e9e", colour="black", size=2.4, stroke=0.3) +
  geom_point(data=viv, aes(shape=partner, fill=eclass), colour="black", size=2.4, stroke=0.3) +
  scale_shape_manual(values=SHP, name="PTEN-loss partner") +
  scale_fill_manual(values=c("Buffer"="#3b6fb0","Synergistic"="#b2182b"),
                    name="Epistatic with PTEN-loss") +
  geom_text_repel(aes(label=label), size=FS/.MM, lineheight=0.9, seed=1234,
                  segment.size=0.2, min.segment.length=0, box.padding=0.35, max.overlaps=Inf) +
  guides(shape=guide_legend(override.aes=list(fill="#9e9e9e"), order=1),
         fill=guide_legend(override.aes=list(shape=21), order=2)) +
  labs(x=sprintf("PC1 (%.1f%%)", ve[1]), y=sprintf("PC2 (%.1f%%)", ve[2])) +
  fig_theme + sq_panel(130) +
  theme(legend.key=element_blank())
save_panel(p, "Fig3C_condition_pca", 3.6, 2.7, fig_main())
write.csv(pd, tbl("Fig3C_condition_pca_coords.csv"), row.names=FALSE)
cat(sprintf("condition PCA (draft style) written. PC1=%.1f%% PC2=%.1f%%\n", ve[1], ve[2]))
