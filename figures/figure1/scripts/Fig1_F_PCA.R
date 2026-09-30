# =============================================================================
#  Fig1_F_PCA.R — Panel F
#  Paired-cohort PCA (cohort QC + tissue-type separation)
# -----------------------------------------------------------------------------
#  Design (locked in plan): cosmetic-only changes from the v2 draft.
#    - Drop the stale "Figure 1C:" title remnant.
#    - Compress sample labels to short patient + tissue codes.
#    - Emphasize PC1 variance (% on axis).
#    - Group-level 95% normal ellipses (Primary vs LungMet).
#    - Intra-patient connecting lines retained (paired trajectory cue).
#    - Declarative on-panel title:
#        "Paired primary and lung metastasis transcriptomes separate along PC1"
#
#  Requires shared infrastructure from 01_load_GSE110590.R:
#    exprs_sub, group_sub, patient_sub, custom_meta
#
#  Outputs:
#    output/figures/main/Fig1_F_PCA_paired_subset.pdf
#
#  Legacy reference: scripts/legacy/Figure1_publication_monolithic.R  lines 310-357
# =============================================================================

message("[Panel F]  Paired-subset PCA — building...")

# -- 1. PCA on the paired subset ---------------------------------------------
pca_res <- prcomp(t(exprs_sub), scale. = TRUE)
var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

# -- 2. Compress sample labels -----------------------------------------------
# Maps the verbose GSE110590 metadata sample IDs to compact patient + tissue
# codes that fit on the panel without overlap (e.g. "A28-LungMet-RUL-RNA" →
# "A28-RUL"). The transformations are explicit per-pattern rather than
# relying on a generic regex, because the source metadata uses three different
# naming conventions across the 5 patients.
short_label <- function(x) {
  out <- as.character(x)
  out <- sub("-Primary-FFPE-RNA$", "-P",   out)
  out <- sub("-PT-FFPE-RNA$",      "-P",   out)
  out <- sub("-LungMet-RNA$",      "-Met", out)
  out <- sub("-LUNG-MET-RNA$",     "-Met", out)
  out <- sub("-LungMet-L-RNA$",    "-L",   out)
  out <- sub("-LungMet-R-RNA$",    "-R",   out)
  out <- sub("-Lung-Met-1-RNA$",   "-M1",  out)
  out <- sub("-Lung-Met-2-RNA$",   "-M2",  out)
  out <- sub("-LungMet-LLL-RNA$",  "-LLL", out)
  out <- sub("-LungMet-RUL-RNA$",  "-RUL", out)
  out
}

# -- 3. Build tidy plotting frame --------------------------------------------
pca_df <- data.frame(
  Sample  = rownames(pca_res$x),
  PC1     = pca_res$x[, 1],
  PC2     = pca_res$x[, 2],
  Group   = group_sub,
  Patient = patient_sub,
  Label   = short_label(custom_meta$sample_id),
  stringsAsFactors = FALSE
)

# Numerical source data for Figure 1F: exact sample coordinates and grouping
# variables supplied to ggplot below. Variance percentages reproduce the axis
# labels and are repeated per row to keep the CSV self-contained.
pca_source_data <- pca_df
pca_source_data$PC1_variance_percent <- var_explained[1]
pca_source_data$PC2_variance_percent <- var_explained[2]
write.csv(
  pca_source_data[, c("Sample", "Label", "Patient", "Group", "PC1", "PC2",
                      "PC1_variance_percent", "PC2_variance_percent")],
  file.path(OUTPUT_TBL_DIR, "Fig1F_PCA_paired_subset_source_data.csv"),
  row.names = FALSE,
  quote = FALSE
)

# -- 4. Plot -----------------------------------------------------------------
p_D <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  # zero-axis guides
  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray70") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "gray70") +
  # group-level 95% normal ellipses (Primary vs LungMet)
  stat_ellipse(aes(group = Group, linetype = Group),
               colour = "gray30", type = "norm", level = 0.95, size = 0.5) +
  # intra-patient connection lines
  geom_line(aes(group = Patient, colour = Patient), alpha = 0.4, size = 0.8) +
  # points: shape = tissue type, fill = patient
  geom_point(aes(fill = Patient, shape = Group),
             size = 5, colour = "black", stroke = 0.6) +
  scale_shape_manual(values = c("Primary" = 22, "LungMet" = 24),
                     labels = c("Primary" = "Primary tumor",
                                "LungMet" = "Lung metastasis"),
                     name   = "Tissue (95% ellipse)") +
  scale_linetype_manual(values = c("Primary" = "solid", "LungMet" = "dashed"),
                        guide  = "none") +
  geom_text_repel(aes(label = Label), size = 3, max.overlaps = 50,
                  segment.color = "grey40", segment.size = 0.2,
                  box.padding = unit(0.4, "lines")) +
  labs(
    title = "Paired primary and lung metastasis transcriptomes separate along PC1",
    subtitle = "GSE110590 paired subset \u2014 5 patients, 5 primary + 8 lung metastasis samples",
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC2 (", var_explained[2], "% variance)")
  ) +
  guides(
    fill   = guide_legend(override.aes = list(shape = 21, size = 5),
                          title = "Patient"),
    colour = "none"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title          = element_text(face = "bold"),
    plot.title.position = "plot",
    panel.grid.minor    = element_blank(),
    legend.position     = "right",
    legend.box          = "vertical"
  )

ggsave(file.path(OUTPUT_FIG_DIR, "Fig1_F_PCA_paired_subset.pdf"),
       p_D, width = 8, height = 6, device = pdf_device)

message("[Panel F]  Done. \u2192 ",
        file.path(OUTPUT_FIG_DIR, "Fig1_F_PCA_paired_subset.pdf"))
