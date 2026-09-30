#!/usr/bin/env Rscript
# =============================================================================
# Supp Fig 3H — SQUARED CENTROID DISTANCE, 2*||mu_pert - mu_ctrl||^2, for each
# perturbation vs the annotated NT*NT in z-scored 30-PC space; significance by
# label permutation. In vitro: large & significant. In vivo: collapses.
#
# NOT the same statistic as main Fig 3F. This is the energy distance under a
# SQUARED-Euclidean cost, where the within-group terms cancel and only the
# centroid separation survives; Fig 3F uses a Euclidean cost (dist()) and keeps
# them. This panel also differs deliberately in its null: it pools NT*NT ACROSS
# lanes, uses R = 1,000 permutations, and applies ONE BH family across both
# contexts, whereas Fig 3F uses same-lane controls, B = 10,000 and per-context
# BH. Describe it that way in the legend; do not call it an E-test.
# v7 repo overhaul 2026-07-15: relocated from XL_code/fig3_panels_new_build/build_edistance_merged.R; was panelEdist_merged -> FigS3H_edist_bars.
# =============================================================================
if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
set.seed(1234)
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(ggplot2) })
NPC <- 30; R <- 1000
PERTS <- c("Pten*NT","Pten*Cdh1","Pten*Cx3cl1","Pten*Cxcr5","Pten*Tlr7")

edist_context <- function(file, tag){
  o <- readRDS(raw(file))
  if (!"pca" %in% Reductions(o)) o <- RunPCA(o, npcs=NPC, verbose=FALSE)
  PCz <- scale(Embeddings(o,"pca")[,1:NPC]); rownames(PCz) <- colnames(o)
  pert <- setNames(as.character(o$perturbation), colnames(o))
  ctrl <- names(pert)[pert=="NT*NT"]
  ed <- function(a,b) 2*sum((colMeans(PCz[a,,drop=FALSE]) - colMeans(PCz[b,,drop=FALSE]))^2)
  rows <- lapply(PERTS, function(p){
    pc <- names(pert)[pert==p]; if (length(pc)<15) return(NULL)
    obs <- ed(pc, ctrl); pool <- c(pc, ctrl); n1 <- length(pc)
    Z <- PCz[pool,,drop=FALSE]
    nul <- replicate(R, { i <- sample.int(nrow(Z)); 2*sum((colMeans(Z[i[1:n1],,drop=FALSE]) - colMeans(Z[i[-(1:n1)],,drop=FALSE]))^2) })
    data.frame(context=tag, perturbation=p, edist=obs, p=(1+sum(nul>=obs))/(R+1)) })
  rm(o,PCz); gc(verbose=FALSE); do.call(rbind, rows)
}
E <- rbind(edist_context("Perturb_InVitro_HTO.rds","In vitro"),
           edist_context("Perturb_InVivo_merged.rds","In vivo"))
E$padj <- p.adjust(E$p,"BH")
E$star <- ifelse(E$padj<0.001,"***", ifelse(E$padj<0.01,"**", ifelse(E$padj<0.05,"*","n.s.")))
write.csv(E, tbl("edistance_merged.csv"), row.names=FALSE)
cat("\n=== E-distance (energy distance vs merged NT*NT; E-test) ===\n"); print(E, row.names=FALSE, digits=3)

E$perturbation <- factor(pair_lab(E$perturbation), levels=pair_lab(PERTS))   # display "_"; data keep "*"
E$context <- factor(E$context, levels=c("In vitro","In vivo"))
p <- ggplot(E, aes(perturbation, edist, fill=context)) +
  geom_col(position=position_dodge(width=0.72), width=0.66, colour="black", linewidth=0.2) +
  geom_text(aes(label=star), position=position_dodge(width=0.72), vjust=-0.3, size=(FS-1)/.MM) +
  scale_fill_manual(values=c("In vitro"="#2166ac","In vivo"="#d62728"), name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,0.14))) +
  # NOTE: this panel plots 2*||mu_pert - mu_ctrl||^2 (see `ed()` above), i.e. a squared
  # CENTROID separation. That equals the energy distance only under a SQUARED-Euclidean
  # cost; main Fig 3F uses a Euclidean cost (dist()) and so computes a different
  # statistic. Do not label this "E-distance" - the two panels are not comparable.
  # The exact definition and this panel's own null are stated in the supplementary legend.
  labs(x=NULL, y="squared centroid distance", title="Perturbation effect size collapses in vivo") +
  fig_theme + theme(axis.text.x=element_text(angle=35, hjust=1, face="italic"),
                    legend.position=c(0.98,0.98), legend.justification=c(1,1), legend.key.size=unit(6,"pt"))
save_panel(p, "FigS3H_edist_bars", 2.4, 1.9, fig_supp())
cat("\nsaved FigS3H_edist_bars\n")
