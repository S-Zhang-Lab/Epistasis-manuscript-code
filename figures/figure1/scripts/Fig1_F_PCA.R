message("[Panel D]  Paired-subset PCA — building...")

pca_res <- prcomp(t(exprs_sub), scale. = TRUE)
var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

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

pca_df <- data.frame(
  PC1     = pca_res$x[, 1],
  PC2     = pca_res$x[, 2],
  Group   = group_sub,
  Patient = patient_sub,
  Label   = short_label(custom_meta$sample_id),
  stringsAsFactors = FALSE
)

p_D <- ggplot(pca_df, aes(x = PC1, y = PC2)) +

  geom_hline(yintercept = 0, linetype = "dotted", colour = "gray70") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "gray70") +

  stat_ellipse(aes(group = Group, linetype = Group),
               colour = "gray30", type = "norm", level = 0.95, size = 0.5) +

  geom_line(aes(group = Patient, colour = Patient), alpha = 0.4, size = 0.8) +

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

message("[Panel D]  Done. \u2192 ",
        file.path(OUTPUT_FIG_DIR, "Fig1_F_PCA_paired_subset.pdf"))
