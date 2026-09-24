suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "theme_pub.R"))
  source(here::here("R", "palettes.R"))
  source(here::here("R", "helpers.R"))
  for (pkg in c("ggplot2", "dplyr", "ggrepel")) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
})
SEED  <- 123L
N_LAB <- 8L

md <- load_merged_data() |>
  dplyr::filter(!grepl("NT", gene)) |>
  dplyr::mutate(gene_label = gsub("\\*", "×", gene))

lab <- dplyr::bind_rows(
  dplyr::slice_max(md, diff_niche_dep, n = N_LAB, with_ties = FALSE),
  dplyr::slice_min(md, diff_niche_dep, n = N_LAB, with_ties = FALSE))

lim <- max(abs(c(md$dLFC_Cx, md$dLFC_Mrc)), na.rm = TRUE) * 1.05

p <- ggplot(md, aes(dLFC_Cx, dLFC_Mrc)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "gray55", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  geom_point(aes(color = diff_niche_dep), size = 1.5, alpha = 0.85, stroke = 0) +
  scale_colour_gradient2(low = HOST_COLORS[["Cx3cr1"]], mid = "gray85",
                         high = HOST_COLORS[["Mrc1"]], midpoint = 0,
                         name = "Differential\nniche-dependence\n(Mrc1 − Cx3cr1)") +
  geom_point(data = lab, shape = 21, fill = NA, color = "gray25",
             size = 2.2, stroke = 0.4, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = lab, aes(label = gene_label),
                           size = 1.9, fontface = "bold",
                           min.segment.length = 0, max.overlaps = Inf,
                           box.padding = 0.4, segment.size = 0.25,
                           segment.alpha = 0.55, seed = SEED,
                           show.legend = FALSE) +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(title = "First-order: which niche sets each tumor-cell dependency",
       subtitle = "per-pair ΔLFC (LFC_Dep − LFC_WT); above y=x = Mrc1-biased (WT-filtered, n=269)",
       x = "ΔLFC under Cx3cr1+ depletion", y = "ΔLFC under Mrc1+ depletion") +
  theme_pub() +
  theme(legend.position = "right",
        legend.text = element_text(size = 5), legend.title = element_text(size = 5.2),
        legend.key.height = unit(14, "pt"), legend.key.width = unit(6, "pt"),
        plot.subtitle = element_text(size = 6, color = "gray30"))

save_pdf(p, fig("Fig5_first_order_scatter.pdf"), 121, 121)
message(sprintf("[Fig5_F] first-order per-pair scatter; n=%d pairs", nrow(md)))
log_session("Fig5_first_order_scatter")
