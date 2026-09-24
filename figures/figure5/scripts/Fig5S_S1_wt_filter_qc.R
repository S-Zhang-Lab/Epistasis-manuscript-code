suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "theme_pub.R"))
  source(here::here("R", "palettes.R"))
  source(here::here("R", "helpers.R"))
  for (pkg in c("ggplot2", "dplyr", "readr")) library(pkg, character.only = TRUE)
})

mk <- function(path, nm) readr::read_csv(path, show_col_types = FALSE) |>
  dplyr::filter(type == "dual") |>
  dplyr::transmute(gene, !!nm := mean_top_LFC_WT_vs_Cell_corrected)
m <- dplyr::inner_join(mk(tbl("geneLFC_Cx3cr1.csv"), "WT_Cx"),
                       mk(tbl("geneLFC_Mrc1.csv"),   "WT_Mrc"), by = "gene") |>
  dplyr::mutate(removed = (WT_Cx * WT_Mrc < 0) & abs(WT_Mrc - WT_Cx) > 2,
                status = ifelse(removed, "removed (discordant WT)", "retained"))
n_rm <- sum(m$removed); n_tot <- nrow(m)
lim <- 5

p <- ggplot(m, aes(WT_Cx, WT_Mrc, color = status)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray55", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "gray75", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "gray75", linewidth = 0.3) +
  geom_point(aes(size = removed), alpha = 0.8, stroke = 0) +
  scale_colour_manual(values = c("retained" = "gray65", "removed (discordant WT)" = "#B2182B"), name = NULL) +
  scale_size_manual(values = c(`FALSE` = 1.1, `TRUE` = 1.9), guide = "none") +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(title = "WT-consistency filter (first-order layer)",
       subtitle = sprintf("WT differs across screens; %d of %d dual pairs removed (opposite-signed WT, |ΔWT|>2)", n_rm, n_tot),
       x = "WT LFC (Cx3cr1 screen)", y = "WT LFC (Mrc1 screen)") +
  theme_pub() +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.text = element_text(size = 5.2), legend.key.size = unit(7, "pt"),
        legend.background = element_rect(fill = alpha("white", 0.85), color = "gray80", linewidth = 0.3),
        plot.subtitle = element_text(size = 5.4, color = "gray25"))

save_pdf(p, fig("Fig5S_S1_wt_filter_qc.pdf"), 121, 110)
message(sprintf("[Fig5S_S1_wt_filter_qc] removed %d of %d dual pairs", n_rm, n_tot))
log_session("Fig5S_S1_wt_filter_qc")
