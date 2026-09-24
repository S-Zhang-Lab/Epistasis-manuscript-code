suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "theme_pub.R"))
  source(here::here("R", "palettes.R"))
  source(here::here("R", "helpers.R"))
  for (pkg in c("ggplot2", "dplyr", "tidyr", "readr", "ggrepel")) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
})
SEED <- 123L

gi <- compute_synergy_table() |>
  dplyr::transmute(gene,
    exp_Cx  = SingleA_Cx + SingleB_Cx,  obs_Cx  = Effect_Cx,  dGI_Cx  = Synergy_Cx,
    exp_Mrc = SingleA_Mrc + SingleB_Mrc, obs_Mrc = Effect_Mrc, dGI_Mrc = Synergy_Mrc)

qtab <- readr::read_csv(tbl("ddGI_combined.csv"), show_col_types = FALSE, progress = FALSE) |>
  dplyr::select(gene, q_Cx = q_dGI_Cx, q_Mrc = q_dGI_Mrc)

long <- dplyr::bind_rows(
  dplyr::transmute(gi, gene, host = "Cx3cr1+ depletion", expected = exp_Cx,  observed = obs_Cx,  dGI = dGI_Cx),
  dplyr::transmute(gi, gene, host = "Mrc1+ depletion",   expected = exp_Mrc, observed = obs_Mrc, dGI = dGI_Mrc)) |>
  dplyr::left_join(dplyr::bind_rows(
      dplyr::transmute(qtab, gene, host = "Cx3cr1+ depletion", q = q_Cx),
      dplyr::transmute(qtab, gene, host = "Mrc1+ depletion",   q = q_Mrc)),
    by = c("gene", "host")) |>
  dplyr::mutate(sig = !is.na(q) & q < 0.05,
                pair = gsub("\\*", "×", gene))

NTOP <- 4L
hero <- "Pten*Cdh1"
sg   <- dplyr::bind_rows(
  long |> dplyr::filter(sig, dGI > 0) |> dplyr::group_by(host) |>
    dplyr::slice_max(dGI, n = NTOP, with_ties = FALSE) |> dplyr::ungroup(),
  long |> dplyr::filter(sig, dGI < 0) |> dplyr::group_by(host) |>
    dplyr::slice_min(dGI, n = NTOP, with_ties = FALSE) |> dplyr::ungroup(),
  long |> dplyr::filter(gene == hero, sig)) |>
  dplyr::distinct(gene, host, .keep_all = TRUE)
ns   <- dplyr::anti_join(long, sg, by = c("gene", "host"))
labs <- sg
lim  <- max(abs(c(long$expected, long$observed)), na.rm = TRUE) * 1.02

p <- ggplot(mapping = aes(expected, observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray55", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "gray85", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "gray85", linewidth = 0.25) +
  geom_point(data = ns, aes(fill = dGI), shape = 21, size = 0.9, stroke = 0,   alpha = 0.7) +
  geom_point(data = sg, aes(fill = dGI), shape = 21, size = 1.9, stroke = 0.4, color = "gray15", alpha = 0.95) +
  scale_fill_gradient2(low = "#1B7837", mid = "gray85", high = "#762A83", midpoint = 0,
                       name = "ΔGI\n(obs − exp)") +
  ggrepel::geom_text_repel(data = labs, aes(label = pair), size = 2.0, fontface = "bold",
                           color = "gray15", box.padding = 0.45, point.padding = 0.3,
                           min.segment.length = 0, segment.size = 0.22,
                           max.overlaps = Inf, seed = SEED) +
  facet_wrap(~host) +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(title = "Epistasis: observed double vs additive single-gene expectation",
       subtitle = "diagonal = additive; off-diagonal = ΔGI (purple +, green −). Labeled = strongest each direction + Pten×Cdh1 exemplar",
       x = "expected if additive  (single A + single B)", y = "observed (double KO)") +
  theme_pub() +
  theme(legend.position = "right",
        legend.text = element_text(size = 5), legend.title = element_text(size = 5.2),
        legend.key.height = unit(12, "pt"), legend.key.width = unit(5, "pt"),
        plot.subtitle = element_text(size = 5.3, color = "gray30"),
        panel.spacing = unit(5, "pt"))

save_pdf(p, fig("Fig5_additivity.pdf"), 121, 80)
message("[Fig5_additivity] observed-vs-expected epistasis, faceted by host")
log_session("Fig5_additivity")
