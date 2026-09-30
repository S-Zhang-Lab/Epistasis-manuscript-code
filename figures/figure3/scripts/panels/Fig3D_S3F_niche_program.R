#!/usr/bin/env Rscript
# =============================================================================
# Fig. 3D / Supplementary Fig. S3F — the NICHE-ADAPTATION PROGRAM (what the metastatic niche imposes).
# The pure environmental effect: NT*NT cells (the common control genotype) in vitro
# vs in vivo. HARMONIZED RNA -> logCPM (SCT models are dataset-specific, not
# comparable across objects). Two parts:
#   (C1) volcano of niche-adaptation DEGs (in vivo vs in vitro)
#   (C2) Hallmark pathways the niche turns on/off (UCell, Cohen d in vivo vs in vitro)
#
# v7 repo overhaul 2026-07 : relocated from XL_code/fig3_panels_new_build/build_panelC_niche_program.R;
#   was panelC_niche_pathways -> Fig3D_niche_program_hallmark (main), panelC_niche_volcano -> FigS3F_niche_volcano (supp).
# =============================================================================
if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
set.seed(1234)
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(presto); library(UCell); library(matrixStats); library(ggplot2) })
HALL <- readRDS(cache("hallmark_common.rds"))

get_ntnt <- function(file){
  o <- readRDS(raw(file))
  if (length(SeuratObject::Layers(o[["RNA"]]))>1) o <- JoinLayers(o, assay="RNA")
  nt <- colnames(o)[as.character(o$perturbation)=="NT*NT"]
  m <- GetAssayData(o, assay="RNA", layer="counts")[, nt, drop=FALSE]; rm(o); gc(verbose=FALSE); m
}
Cv <- get_ntnt("Perturb_InVitro_HTO.rds"); Co <- get_ntnt("Perturb_InVivo_merged.rds")
g  <- intersect(rownames(Cv), rownames(Co)); C <- cbind(Cv[g,,drop=FALSE], Co[g,,drop=FALSE])
env <- factor(c(rep("In vitro",ncol(Cv)), rep("In vivo",ncol(Co))), levels=c("In vitro","In vivo"))
cat(sprintf("NT*NT cells: in vitro %d, in vivo %d; common genes %d\n", ncol(Cv), ncol(Co), length(g)))

## ---- harmonized logCPM + DE (in vivo vs in vitro) ---------------------------
cpm  <- t(t(C)/colSums(C))*1e6; lcpm <- log1p(cpm)
w <- presto::wilcoxauc(lcpm, as.character(env)); w <- w[w$group=="In vivo",]
mvivo <- rowMeans(cpm[, env=="In vivo"]); mvitro <- rowMeans(cpm[, env=="In vitro"])
DE <- data.frame(gene=rownames(C), log2FC=log2((mvivo+1)/(mvitro+1)),
                 padj=w$padj[match(rownames(C), w$feature)],
                 pct_vivo=w$pct_in[match(rownames(C), w$feature)], pct_vitro=w$pct_out[match(rownames(C), w$feature)])
DE <- DE[is.finite(DE$log2FC) & !is.na(DE$padj) & (DE$pct_vivo>10 | DE$pct_vitro>10), ]
DE$sig <- DE$padj<0.05 & abs(DE$log2FC)>1
write.csv(DE, tbl("FigS3F_niche_DE.csv"), row.names=FALSE)
cat(sprintf("niche-adaptation DEGs (padj<0.05,|log2FC|>1): %d up in vivo, %d up in vitro\n",
            sum(DE$sig & DE$log2FC>0), sum(DE$sig & DE$log2FC<0)))

## ---- Hallmark pathways (UCell, Cohen d in vivo vs in vitro) -----------------
U <- t(ScoreSignatures_UCell(C, features=HALL, name="", ncores=1))
cohen <- function(a,b) (mean(a)-mean(b))/sqrt(((length(a)-1)*var(a)+(length(b)-1)*var(b))/(length(a)+length(b)-2))
PA <- data.frame(pathway=sub("HALLMARK_","",rownames(U)),
                 cohen_d=apply(U,1,function(r) cohen(r[env=="In vivo"], r[env=="In vitro"])),
                 padj=p.adjust(apply(U,1,function(r) wilcox.test(r[env=="In vivo"], r[env=="In vitro"])$p.value),"BH"))
write.csv(PA, tbl("Fig3D_niche_pathways.csv"), row.names=FALSE)

## ---- C1 volcano (compact ~115pt, visible dots, 5pt italic gene labels) ------
DE$y <- pmin(-log10(DE$padj), 300)
up <- DE[DE$sig & DE$log2FC>0,]; dn <- DE[DE$sig & DE$log2FC<0,]
lab <- rbind(head(up[order(-up$log2FC),],6), head(dn[order(dn$log2FC),],3))   # 9 genes only
pC1 <- ggplot(DE, aes(log2FC, y)) +
  geom_point(data=subset(DE,!sig), colour="grey82", size=0.35, stroke=0) +
  geom_point(data=subset(DE, sig), aes(colour=log2FC>0), size=0.9, stroke=0) +   # sig dots BIG
  geom_vline(xintercept=c(-1,1), linetype="dashed", colour="grey75", linewidth=0.25) +
  ggrepel::geom_text_repel(data=lab, aes(label=gene), size=FS/.MM, fontface="italic",
                           max.overlaps=Inf, min.segment.length=0, segment.size=0.2, box.padding=0.35, seed=1) +
  scale_colour_manual(values=c(`TRUE`="#b2182b",`FALSE`="#2166ac"), guide="none") +
  labs(x="log2FC (vivo/vitro)", y="-log10 padj", title="Niche-adaptation program") +
  fig_theme
save_panel(pC1, "FigS3F_niche_volcano", 1.6, 1.6, fig_supp())

## ---- C2 Hallmark bar (compact; 13 pathways at 5pt) -------------------------
ps <- PA[PA$padj<0.05,]; ps <- ps[order(-ps$cohen_d),]; ps <- rbind(head(ps,7), tail(ps,6))
ps$pathway <- factor(ps$pathway, levels=rev(ps$pathway))
pC2 <- ggplot(ps, aes(cohen_d, pathway, fill=cohen_d>0)) +
  geom_col(colour="black", linewidth=0.2, width=0.72) + geom_vline(xintercept=0, linewidth=0.3) +
  scale_fill_manual(values=c(`TRUE`="#b2182b",`FALSE`="#2166ac"), guide="none") +
  labs(x="Cohen d (vivo - vitro)", y=NULL, title="Hallmark: niche on / off") +
  fig_theme + theme(axis.text.y=element_text(size=FS))
save_panel(pC2, "Fig3D_niche_program_hallmark", 2.9, 1.7, fig_main())

cat("\nsaved FigS3F_niche_volcano, Fig3D_niche_program_hallmark\n")
cat("\ntop UP in vivo (niche-induced):\n"); print(head(DE[DE$sig&DE$log2FC>0,][order(-DE$log2FC[DE$sig&DE$log2FC>0]),c("gene","log2FC")],12), row.names=FALSE)
cat("\ntop Hallmark UP in vivo:\n"); print(head(PA[order(-PA$cohen_d),c("pathway","cohen_d","padj")],10), row.names=FALSE)
cat("\ntop Hallmark DOWN in vivo (up in vitro):\n"); print(head(PA[order(PA$cohen_d),c("pathway","cohen_d","padj")],6), row.names=FALSE)
