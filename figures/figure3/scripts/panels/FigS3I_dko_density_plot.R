#!/usr/bin/env Rscript
# =============================================================================
# Concrete "what niche compression means" prototypes (Supplementary Fig. S3I), from
# compression_concrete_cells.csv. Three variations to review:
#   V1 embedding   : local PCA of each 2x2 -> genotype clouds SEPARATE in vitro,
#                    MERGE in vivo (free axes; PC coords are arbitrary per panel).
#   V2 score-dens  : distribution of the DOUBLE-KO signature score in DKO cells vs
#                    NT*NT cells. In vitro two humps apart (signature IDs its cells);
#                    in vivo they overlap (signature no longer separates) -> graded.
#   V3 score-violin: same score, all four genotypes as violins.
# v7 repo overhaul 2026-07-15: relocated from XL_code/fig3_panels_new_build/build_compression_concrete_panels.R; was panelE_concrete_scoredens -> FigS3I_dko_signature_density. KEEP ONLY the V2 density panel; V1 embedding and V3 violin variants dropped.
# =============================================================================
if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(ggplot2) })
D <- read.csv(cache("compression_concrete_cells.csv"), check.names=FALSE)

PARTNERS <- c("Cdh1","Cx3cl1","Cxcr5","Tlr7")
typ <- function(g,pn) ifelse(g=="NT*NT","ctrl", ifelse(g=="Pten*NT","pten", ifelse(g==paste0("NT*",pn),"partner","double")))
D$type <- mapply(typ, D$genotype, D$partner)
D$type <- factor(D$type, levels=c("ctrl","pten","partner","double"))
D$partner <- factor(D$partner, levels=PARTNERS)
D$context <- factor(D$context, levels=c("In vitro","In vivo"))
TC <- c(ctrl="#8f97c9", pten="#e69138", partner="#7cbaa8", double="#b0463a")
TLAB <- c(ctrl="NT_NT", pten="Pten_NT", partner="NT_par", double="Pten_par")   # display "_" (house style)

## ---- V1 embedding (DROPPED in v7 overhaul; kept only V2 density panel) --------

## ---- V2 score density: DKO signature IDs its cells in vitro, not in vivo ------
D2 <- subset(D, type %in% c("ctrl","double")); D2$type <- droplevels(D2$type)
p2 <- ggplot(D2, aes(dko_score, colour=type)) +
  geom_density(fill=NA, linewidth=0.6) +
  facet_grid(context ~ partner) +
  scale_colour_manual(values=TC[c("ctrl","double")], labels=TLAB[c("ctrl","double")], name=NULL) +
  labs(x="double-KO signature score (z)", y="cell density",
       title="The double-KO signature separates its cells in vitro, not in vivo") +
  fig_theme +
  theme(strip.background=element_blank(), strip.text=element_text(size=FS,face="bold"),
        panel.spacing=unit(3,"pt"), legend.position="top",
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour=guide_legend(nrow=1))
save_panel(p2, "FigS3I_dko_signature_density", 5.0, 2.5, fig_supp())

## ---- V3 score violin (DROPPED in v7 overhaul; kept only V2 density panel) -----

cat("wrote FigS3I_dko_signature_density\n")
