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

comb <- readr::read_csv(tbl("ddGI_combined.csv"), show_col_types = FALSE,
                        progress = FALSE)

st <- comb |>
  dplyr::group_by(geneA) |>
  dplyr::summarize(mean_dd = mean(ddGI, na.rm = TRUE),
                   sd_dd   = sd(ddGI,  na.rm = TRUE),
                   n       = dplyr::n(), .groups = "drop") |>
  dplyr::filter(n >= 3)

p <- ggplot(st, aes(mean_dd, reorder(geneA, mean_dd))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55",
             linewidth = 0.4) +
  geom_errorbar(aes(xmin = mean_dd - sd_dd, xmax = mean_dd + sd_dd),
                orientation = "y", width = 0, color = "gray65", linewidth = 0.4) +
  geom_point(aes(fill = mean_dd), size = 2.8, shape = 21, stroke = 0.4,
             color = "gray25") +
  scale_fill_gradient2(low = HOST_COLORS[["Mrc1"]], mid = "gray90",
                       high = HOST_COLORS[["Cx3cr1"]], midpoint = 0,
                       name = "mean ΔΔGI") +
  labs(title = "Per-driver ΔΔGI across partners (mean ± SD)",
       subtitle = "ΔΔGI = ΔGI_Cx − ΔGI_Mrc (+ = stronger in Cx3cr1⁺ cell depleted host) · narrow bar = consistent across partners",
       x = "mean ΔΔGI across partners", y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 5.5),
        legend.position = "right",
        plot.subtitle = element_text(size = 5.6, color = "gray30"))

save_pdf(p, fig("Fig5S_S7_stability.pdf"), 121, 100)
message("[Fig5S_S7] driver-anchor ddGI stability")
log_session("Fig5S_S7_stability")
