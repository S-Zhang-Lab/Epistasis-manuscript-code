#!/usr/bin/env Rscript
# scripts/Fig5S_S2_sign_distribution.R ----------------------------------------
# FINAL SUBMISSION PANEL: Supplementary Fig. 5H (the internal S2 asset name
# predates the final supplementary-panel lettering). Sign-distribution QC,
# companion to the WT-consistency filter:
# after filtering, the WT->Depleted sign categories per host. Most pairs are
# WT(+),Dep(+) (consistent positive selection); the minority quadrants are the
# residual cross-host disagreement the filter is built to contain. Ported +
# restyled from archive Fig5B_pie_chart.
#
# Inputs:  data/derived/geneLFC_{Cx3cr1,Mrc1}.csv  (via load_merged_data())
# Outputs: output/figures/FigS5H_sign_distribution.pdf
#          output/tables/FigS5H_sign_distribution_source_data.csv

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

md <- load_merged_data() |> dplyr::filter(!grepl("NT", gene))  # dual pairs only

quad <- function(wt, dep) dplyr::case_when(
  wt >  0 & dep >  0 ~ "WT(+) Dep(+)",
  wt >  0 & dep <= 0 ~ "WT(+) Dep(-)",
  wt <= 0 & dep >  0 ~ "WT(-) Dep(+)",
  TRUE               ~ "WT(-) Dep(-)")

pie <- dplyr::bind_rows(
  dplyr::tibble(host = "Cx3cr1+ depletion", cat = quad(md$WT_Cx,  md$Dep_Cx)),
  dplyr::tibble(host = "Mrc1+ depletion",   cat = quad(md$WT_Mrc, md$Dep_Mrc))) |>
  dplyr::count(host, cat) |>
  dplyr::group_by(host) |>
  dplyr::mutate(total_pairs = sum(n), proportion = n / total_pairs,
                percent = 100 * proportion) |>
  dplyr::ungroup() |>
  dplyr::mutate(cat = factor(cat, levels = names(QUADRANT_COLORS)))

# Numerical source data for the exact host-by-quadrant wedges plotted in
# Supplementary Fig. 5H.
readr::write_csv(
  pie |>
    dplyr::transmute(
      host,
      sign_quadrant = as.character(cat),
      n_pairs = n,
      total_pairs,
      proportion,
      percent
    ),
  tbl("FigS5H_sign_distribution_source_data.csv")
)

p <- ggplot(pie, aes(x = "", y = n, fill = cat)) +
  geom_col(width = 1, color = "white", linewidth = 0.3) +
  coord_polar("y", start = 0) +
  facet_wrap(~host) +
  scale_fill_manual(values = QUADRANT_COLORS, name = NULL, drop = FALSE) +
  labs(title = "Sign distribution: WT vs depleted (post-filter)",
       subtitle = sprintf("WT->Depleted sign categories per host; n=%d pairs", nrow(md))) +
  theme_void(base_family = PUB_FONT_FAMILY) +
  theme(plot.title    = element_text(size = PUB_FONT_PANEL, hjust = 0),
        plot.subtitle = element_text(size = PUB_FONT_BODY, color = "gray30", hjust = 0),
        strip.text    = element_text(size = PUB_FONT_TITLE),
        legend.text   = element_text(size = PUB_FONT_LEGEND),
        legend.position = "right",
        plot.margin   = margin(2, 2, 2, 2, "mm"))

save_pdf(p, fig("FigS5H_sign_distribution.pdf"), 121, 65)
message("[FigS5H] sign-distribution pies")
log_session("FigS5H_sign_distribution")
