suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "theme_pub.R"))
  source(here::here("R", "palettes.R"))
  source(here::here("R", "helpers.R"))
  for (pkg in c("ggplot2", "dplyr", "readr", "tibble", "patchwork")) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
})

HERO <- "Pten*Cdh1"
syn  <- compute_synergy_table() |> dplyr::filter(gene == HERO)
h    <- readr::read_csv(tbl("ddGI_combined.csv"), show_col_types = FALSE, progress = FALSE) |>
  dplyr::filter(gene == HERO)
stopifnot(nrow(syn) == 1, nrow(h) == 1)
HC <- c("Cx3cr1+ depletion" = HOST_COLORS[["Cx3cr1"]], "Mrc1+ depletion" = HOST_COLORS[["Mrc1"]])

ip <- dplyr::bind_rows(
  tibble::tibble(host = "Cx3cr1+ depletion", pten = rep(c("Pten WT", "Pten KO"), each = 2),
                 cdh1 = rep(c("Cdh1 WT", "Cdh1 KO"), 2), y = c(0, syn$SingleB_Cx,  syn$SingleA_Cx,  syn$Effect_Cx)),
  tibble::tibble(host = "Mrc1+ depletion",   pten = rep(c("Pten WT", "Pten KO"), each = 2),
                 cdh1 = rep(c("Cdh1 WT", "Cdh1 KO"), 2), y = c(0, syn$SingleB_Mrc, syn$SingleA_Mrc, syn$Effect_Mrc))) |>
  dplyr::mutate(host = factor(host, levels = names(HC)), cdh1 = factor(cdh1, levels = c("Cdh1 WT", "Cdh1 KO")))
expd <- tibble::tibble(host = factor(names(HC), levels = names(HC)),
                       cdh1 = factor("Cdh1 KO", levels = c("Cdh1 WT", "Cdh1 KO")),
                       y = c(syn$SingleA_Cx + syn$SingleB_Cx, syn$SingleA_Mrc + syn$SingleB_Mrc))
obsd  <- ip |> dplyr::filter(pten == "Pten KO", cdh1 == "Cdh1 KO")
gapdf <- data.frame(host = expd$host, yexp = expd$y, yobs = obsd$y,
                    lab = c(sprintf("synergy %+.2f", syn$Synergy_Cx), sprintf("synergy %+.2f", syn$Synergy_Mrc)))

p_ip <- ggplot(ip, aes(cdh1, y, group = pten, color = pten)) +
  geom_hline(yintercept = 0, color = "gray85", linewidth = 0.3) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  geom_point(data = expd, aes(cdh1, y), inherit.aes = FALSE, shape = 21, size = 2.4,
             fill = "white", color = "gray45", stroke = 0.5) +
  geom_segment(data = gapdf, aes(x = 2.1, xend = 2.1, y = yexp, yend = yobs), inherit.aes = FALSE,
               color = "gray15", linewidth = 0.5, arrow = grid::arrow(ends = "both", length = unit(2, "pt"), type = "closed")) +
  geom_text(data = gapdf, aes(x = 2.16, y = (yexp + yobs) / 2, label = lab), inherit.aes = FALSE,
            hjust = 0, size = 1.8, fontface = "bold", color = "gray15") +
  scale_colour_manual(values = c("Pten WT" = "gray55", "Pten KO" = "#C51B7D"), name = NULL) +
  facet_wrap(~host) + scale_x_discrete(expand = expansion(add = c(0.25, 0.95))) +
  labs(title = "Genetic-interaction plot (non-parallel = epistasis)",
       subtitle = "lines parallel = additive; Pten-KO overshoots its additive expectation (○) only under Mrc1 loss",
       x = NULL, y = "effect on niche depletion (Dep − WT)") +
  theme_pub() +
  theme(legend.position = "top", plot.subtitle = element_text(size = 5.0, color = "gray25"),
        panel.spacing = unit(6, "pt"))

forest <- tibble::tibble(
  quantity = factor(c("ΔGI_Cx", "ΔGI_Mrc", "ΔΔGI (Cx−Mrc)"),
                    levels = rev(c("ΔGI_Cx", "ΔGI_Mrc", "ΔΔGI (Cx−Mrc)"))),
  est = c(h$dGI_Cx, h$dGI_Mrc, h$ddGI),
  lo  = c(h$dGI_Cx_ci_lo, h$dGI_Mrc_ci_lo, h$ddGI_ci_lo),
  hi  = c(h$dGI_Cx_ci_hi, h$dGI_Mrc_ci_hi, h$ddGI_ci_hi),
  q   = c(h$q_dGI_Cx, h$q_dGI_Mrc, h$ddGI_qvalue),
  kind = c("screen", "screen", "ddGI"))
p_for <- ggplot(forest, aes(est, quantity, color = kind)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.18, linewidth = 0.5) +
  geom_point(size = 2.4) +
  geom_text(aes(label = sprintf("q=%.3f", q)), vjust = -1.0, size = 1.9, color = "gray25") +
  scale_colour_manual(values = c(screen = "gray35", ddGI = "#E08214"), guide = "none") +
  labs(title = "ΔGI / ΔΔGI with 95% bootstrap CIs", subtitle = "dashed = null (0)",
       x = "interaction effect", y = NULL) +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 5.0, color = "gray25"), plot.margin = unit(c(2, 2, 2, 9), "pt"))

p <- p_ip + p_for + patchwork::plot_layout(widths = c(1.25, 1))
save_pdf(p, fig("Fig5S_S12_PtenCdh1_detail.pdf"), 183, 80)
message("[Fig5S_S12] Pten×Cdh1 detail: interaction plot + ΔΔGI forest (95% CI)")
log_session("Fig5S_S12_PtenCdh1_detail")
