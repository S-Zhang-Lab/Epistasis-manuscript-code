suppressPackageStartupMessages({library(dplyr); library(ggplot2)})

if (exists("REPO_ROOT")) {
  ROOT <- REPO_ROOT
} else {
  .args <- commandArgs(trailingOnly = FALSE)
  .fa   <- sub("^--file=", "", .args[grepl("^--file=", .args)])
  ROOT  <- if (length(.fa) > 0) normalizePath(file.path(dirname(.fa), "..")) else normalizePath(getwd())
}
TBL  <- file.path(ROOT, "output/tables/DEG_results_paired_subset_with_Hallmark_modules.csv")
OUT  <- file.path(ROOT, "output/figures/main")

MODULE_LEVELS <- c("Immune_Inflammation","Interferon","Developmental_Signaling",
                   "EMT_Adhesion","Hormone_Response","Hypoxia_Metabolism",
                   "Proliferation","Stress_Response","Other")
pretty <- function(x) gsub("_", " ", x)
SIG_FDR <- 0.05; SIG_LFC <- 1.0

d  <- read.csv(TBL)
db <- d %>%
  filter(!is.na(adj.P.Val), adj.P.Val < SIG_FDR, abs(logFC) > SIG_LFC,
         Hallmark_Module_Collapsed %in% MODULE_LEVELS) %>%
  mutate(direction = ifelse(logFC > 0, "Up", "Down")) %>%
  count(Hallmark_Module_Collapsed, direction, name = "Count") %>%
  mutate(Module    = factor(pretty(Hallmark_Module_Collapsed), levels = pretty(MODULE_LEVELS)),
         direction = factor(direction, levels = c("Up","Down")))

p <- ggplot(db, aes(Module, Count, fill = direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72,
           color = "grey25", linewidth = 0.25) +
  geom_text(aes(label = Count, color = direction),
            position = position_dodge(width = 0.8),
            vjust = -0.35, size = 3, fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c(Up = "#C0392B", Down = "#2166AC"),
                    labels = c(Up = "Up in Lung Met", Down = "Down in Lung Met"),
                    name = NULL) +
  scale_color_manual(values = c(Up = "#1a1a1a", Down = "#8a8a8a"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(x = NULL, y = "Significant DEG count",
       caption = paste0("Significant DEG = adj.P.Val < 0.05 and |log2FC| > 1 (paired limma, n = 5 paired patients).\n",
                        "Up / Down = higher / lower in lung metastasis; directional, hypothesis-generating (see Methods).")) +
  theme_classic(base_size = 12) +
  theme(axis.text.x     = element_text(angle = 40, hjust = 1, color = "black"),
        axis.text.y     = element_text(color = "black"),
        axis.line       = element_line(color = "black"),
        legend.position = c(0.15, 0.94),
        legend.background = element_blank(),
        legend.key.size = unit(12, "pt"),
        plot.caption    = element_text(hjust = 0, size = 8, color = "grey35"))

ggsave(file.path(OUT, "Fig1_H_module_DEG_bar.pdf"), p,
       width = 9, height = 5, device = "pdf")
if (requireNamespace("svglite", quietly = TRUE))
  ggsave(file.path(OUT, "Fig1_H_module_DEG_bar.svg"), p, width = 9, height = 5)

cat("Saved Fig1_H_module_DEG_bar.pdf",
    if (requireNamespace("svglite", quietly=TRUE)) "+ .svg" else "(svglite absent -> pdf only)",
    "to output/figures/main/\n")
print(db %>% tidyr::pivot_wider(id_cols = Module, names_from = direction, values_from = Count))
