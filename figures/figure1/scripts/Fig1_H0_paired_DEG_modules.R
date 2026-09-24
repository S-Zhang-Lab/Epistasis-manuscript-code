message("[Panel E (new)]  Paired DEG volcano + Hallmark module architecture — building...")

MODULE_LEVELS <- c(
  "Immune_Inflammation",
  "Interferon",
  "Developmental_Signaling",
  "EMT_Adhesion",
  "Hormone_Response",
  "Hypoxia_Metabolism",
  "Proliferation",
  "Stress_Response",
  "Other",
  "None"
)

MODULE_COLORS <- c(
  "Immune_Inflammation"     = "#E08214",
  "Interferon"              = "#FDB863",
  "Developmental_Signaling" = "#5AAE61",
  "EMT_Adhesion"            = "#1B7837",
  "Hormone_Response"        = "#9970AB",
  "Hypoxia_Metabolism"      = "#B2182B",
  "Proliferation"           = "#D6604D",
  "Stress_Response"         = "#542788",
  "Other"                   = "#878787",
  "None"                    = "#E0E0E0"
)

pretty_module <- function(x) gsub("_", " ", x)

SIG_FDR  <- 0.05
SIG_LFC  <- 1.0

message("  [1/4] Paired limma::duplicateCorrelation + eBayes...")

corfit <- duplicateCorrelation(exprs_sub, design_sub, block = patient_sub)
fit_paired <- eBayes(
  lmFit(exprs_sub, design_sub,
        block = patient_sub, correlation = corfit$consensus),
  trend = TRUE
)

deg_paired <- topTable(fit_paired, coef = 2, number = Inf) %>%
  rownames_to_column("Gene") %>%
  arrange(P.Value)

write.csv(deg_paired,
          file.path(OUTPUT_TBL_DIR, "DEG_results_paired_subset.csv"),
          row.names = FALSE)
message("    DEG rows: ", nrow(deg_paired))

message("  [2/4] Joining curated Hallmark_Module_Collapsed annotation...")

mapping_file <- file.path(GENE_LISTS_DIR, "panel_F_hallmark_module_mapping.csv")
if (!file.exists(mapping_file)) {
  stop("Required Panel F curated mapping not found:\n  ", mapping_file,
       "\n\nThis file is the source of truth for Gene -> Hallmark_Module_Collapsed\n",
       "and must exist before Panel F can render. It is committed to the repo;\n",
       "if it is missing, restore it from version control.")
}

prior <- read_csv(mapping_file, show_col_types = FALSE) %>%
  select(Gene, Hallmark_Module, Hallmark_Module_Collapsed) %>%
  distinct(Gene, .keep_all = TRUE)

deg_annot <- deg_paired %>%
  left_join(prior, by = "Gene")

n_matched <- sum(!is.na(deg_annot$Hallmark_Module_Collapsed))
message("    matched mapping for ", n_matched, " of ", nrow(deg_annot),
        " DEG rows (", sprintf("%.1f%%", 100 * n_matched / nrow(deg_annot)), ")")

write.csv(
  deg_annot,
  file.path(OUTPUT_TBL_DIR, "DEG_results_paired_subset_with_Hallmark_modules.csv"),
  row.names = FALSE
)

message("  [3/4] Building volcano + below-volcano bar...")

df_plot <- deg_annot %>%
  mutate(
    Module     = ifelse(is.na(Hallmark_Module_Collapsed),
                        "None", Hallmark_Module_Collapsed),
    Module     = factor(Module, levels = MODULE_LEVELS),
    ModuleNice = pretty_module(as.character(Module)),
    is_sig     = !is.na(adj.P.Val) & adj.P.Val < SIG_FDR & abs(logFC) > SIG_LFC,
    direction  = ifelse(logFC > 0, "Up", "Down")
  ) %>%

  arrange(desc(Module == "None"), Module)

y_max <- max(-log10(df_plot$P.Value), na.rm = TRUE)

p_main <- ggplot(df_plot,
                 aes(x = logFC, y = -log10(P.Value), color = ModuleNice)) +
  geom_point(alpha = 0.65, size = 1.0) +
  scale_color_manual(
    values = setNames(MODULE_COLORS[as.character(MODULE_LEVELS)],
                      pretty_module(MODULE_LEVELS)),
    breaks = pretty_module(setdiff(MODULE_LEVELS, "None")),
    name   = "Hallmark functional module"
  ) +
  geom_vline(xintercept = c(-SIG_LFC, SIG_LFC),
             linetype = "dashed", color = "grey50", size = 0.3) +
  geom_hline(yintercept = -log10(SIG_FDR),
             linetype = "dashed", color = "grey50", size = 0.3) +
  labs(
    title    = "Hallmark module architecture of paired primary \u2192 lung metastasis DEGs",
    subtitle = "GSE110590 paired subset; limma duplicateCorrelation block = patient",
    x        = expression(log[2]~"fold change (LungMet vs Primary)"),
    y        = expression(-log[10]~"(P value)")
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title          = element_text(face = "bold"),
    plot.title.position = "plot",
    panel.grid.minor    = element_blank(),
    legend.position     = "right",
    legend.title        = element_text(face = "bold")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

df_bar <- df_plot %>%
  filter(is_sig, Module != "None") %>%
  count(ModuleNice, direction, name = "Count") %>%
  mutate(
    ModuleNice = factor(ModuleNice,
                        levels = pretty_module(setdiff(MODULE_LEVELS, "None"))),
    direction  = factor(direction, levels = c("Up", "Down"))
  )

p_bar <- ggplot(df_bar, aes(x = ModuleNice, y = Count,
                            fill = ModuleNice, alpha = direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           color = "grey30", size = 0.2) +
  geom_text(aes(label = Count, group = direction),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 2.7) +
  scale_fill_manual(
    values = setNames(MODULE_COLORS[as.character(MODULE_LEVELS)],
                      pretty_module(MODULE_LEVELS)),
    guide  = "none"
  ) +
  scale_alpha_manual(values = c("Up" = 1.0, "Down" = 0.5),
                     name   = "Direction in lung metastasis") +
  labs(
    x = NULL,
    y = "Significant DEG count",
    caption = sprintf("Significant = adj.P.Val < %.2f & |log2 FC| > %g  \u2014  n = 5 paired patients (paired limma); effects directional, hypothesis-generating \u2014 see Methods.",
                      SIG_FDR, SIG_LFC)
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x        = element_text(angle = 35, hjust = 1),
    plot.caption       = element_text(hjust = 0, size = 8,
                                       face = "italic", color = "grey30"),
    legend.position    = "right"
  )

p_F <- (p_main / p_bar) + patchwork::plot_layout(heights = c(3, 1.2))

ggsave(file.path(OUTPUT_SUP_VAR_DIR, "FigS1_paired_volcano_hallmark_modules.pdf"),
       p_F, width = 11, height = 10, device = pdf_device)

message("[Panel F]  Done. \u2192 ",
        file.path(OUTPUT_SUP_VAR_DIR, "FigS1_paired_volcano_hallmark_modules.pdf"))
