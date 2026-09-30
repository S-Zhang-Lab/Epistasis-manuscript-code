#!/usr/bin/env Rscript
# =============================================================================
# Fig S3C -- enAsCas12a cutting efficiency (CRISPResso2 % editing) per target gene.
# Each point = one dual-guide construct/sample measurement (amplicon-NGS/CRISPResso2);
# genes ordered PTEN anchor -> synergy partners (Cdh1, Cx3cl1) -> buffer partners
# (Cxcr5, Tlr7); crossbar = median. Confirms the perturbations were edited before
# injection. v7 overhaul 2026-07-15: built from XL_code/CE/Perturb_Cutting.xlsx.
# =============================================================================
if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

D <- fread(derived("screens", "cas12a_cutting_efficiency.csv"))
GENE_ORD <- c("Pten", "Cdh1", "Cx3cl1", "Cxcr5", "Tlr7")
D[, target_gene := factor(target_gene, levels = GENE_ORD)]
D[, eclass := fifelse(target_gene == "Pten", "PTEN anchor",
              fifelse(target_gene %chin% c("Cdh1", "Cx3cl1"), "Synergistic", "Buffering"))]
D[, eclass := factor(eclass, levels = c("PTEN anchor", "Synergistic", "Buffering"))]

set.seed(1234)
p <- ggplot(D, aes(target_gene, pct_indel, colour = eclass)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.55, linewidth = 0.3,
               colour = "grey35", fatten = 0) +
  geom_jitter(width = 0.14, height = 0, size = 1.1, alpha = 0.9) +
  scale_colour_manual(values = c("PTEN anchor" = "#6e6e6e", "Synergistic" = "#b2182b",
                                 "Buffering" = "#3b6fb0"), name = NULL) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = NULL, y = "% editing (indel, CRISPResso2)",
       title = "Cas12a cutting efficiency") +
  fig_theme + sq_panel(110) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "italic"),
        legend.position = "top", legend.justification = "right",
        panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.2))

save_panel(p, "FigS3C_cutting_efficiency", 2.2, 2.0, fig_supp())
fwrite(D[, .(n = .N, median = round(median(pct_indel), 1),
             min = round(min(pct_indel), 1), max = round(max(pct_indel), 1)), by = target_gene],
       tbl("FigS3C_cutting_efficiency_summary.csv"))
cat("Fig S3C cutting efficiency written (", nrow(D), "measurements )\n")
