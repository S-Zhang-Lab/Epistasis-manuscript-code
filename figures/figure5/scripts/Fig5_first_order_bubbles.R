suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "theme_pub.R"))
  source(here::here("R", "palettes.R"))
  source(here::here("R", "helpers.R"))
  for (pkg in c("ggplot2", "dplyr", "readr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
})

a <- readr::read_csv(tbl("niche_anchor_analysis.csv"), show_col_types = FALSE) |>
  dplyr::filter(layer == "first_order") |>
  dplyr::mutate(
    pos_label = factor(dplyr::recode(pos,
                         driver  = "Tumor driver (geneA)",
                         partner = "Niche / secreted partner (geneB)"),
                       levels = c("Tumor driver (geneA)",
                                  "Niche / secreted partner (geneB)")),
    sig = ifelse(q < 0.05, "q<0.05", "n.s."))

p <- ggplot(a, aes(mean, reorder(gene, mean))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55",
             linewidth = 0.4) +
  geom_point(aes(size = n, fill = mean, color = sig), shape = 21, stroke = 0.5) +
  scale_fill_gradient2(low = HOST_COLORS[["Cx3cr1"]], mid = "gray90",
                       high = HOST_COLORS[["Mrc1"]], midpoint = 0,
                       name = "mean diff.\nniche-dependence") +
  scale_colour_manual(values = c("q<0.05" = "black", "n.s." = "gray70"),
                      name = NULL) +
  scale_size_continuous(range = c(1.2, 4), name = "pairs (n)") +
  scale_x_continuous(breaks = seq(-0.5, 1.5, 0.5), expand = expansion(mult = 0.09)) +
  facet_wrap(~pos_label, scales = "free_y") +
  labs(title = "First-order anchors by screen position",
       subtitle = "per-gene mean ΔLFC_Mrc − ΔLFC_Cx (+ = Mrc1-biased) · CDH1 = strongest partner-side",
       x = "Differential niche-dependence (Mrc1 − Cx3cr1)", y = NULL) +
  guides(color = guide_legend(override.aes = list(size = 3), order = 3),
         fill = guide_colourbar(order = 1, barwidth = unit(42, "pt"), barheight = unit(4, "pt")),
         size = guide_legend(order = 2)) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 4.6),
        legend.position = "bottom", legend.box = "horizontal",
        legend.text = element_text(size = 4.6), legend.title = element_text(size = 5),
        legend.key.size = unit(7, "pt"), legend.box.spacing = unit(2, "pt"),
        legend.spacing.x = unit(4, "pt"),
        plot.subtitle = element_text(size = 5.4, color = "gray30"),
        panel.spacing = unit(4, "pt"))

save_pdf(p, fig("Fig5_first_order_bubbles.pdf"), 121, 118)
message("[Fig5_G] first-order position-split anchor bubbles")
log_session("Fig5_first_order_bubbles")
