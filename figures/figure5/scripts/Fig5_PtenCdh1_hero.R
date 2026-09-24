suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "theme_pub.R"))
  source(here::here("R", "palettes.R"))
  source(here::here("R", "helpers.R"))
  for (pkg in c("ggplot2", "dplyr", "readr", "tibble")) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
})

HERO <- "Pten*Cdh1"
syn  <- compute_synergy_table() |> dplyr::filter(gene == HERO)
comb <- readr::read_csv(tbl("ddGI_combined.csv"), show_col_types = FALSE, progress = FALSE) |>
  dplyr::filter(gene == HERO)
stopifnot(nrow(syn) == 1, nrow(comb) == 1)
qlab <- function(q) if (q < 0.05) sprintf("q=%.3f", q) else sprintf("n.s. (q=%.2f)", q)
HC <- c("Cx3cr1+ depletion" = HOST_COLORS[["Cx3cr1"]], "Mrc1+ depletion" = HOST_COLORS[["Mrc1"]])

mk <- function(host, pten, cdh1, obs, synv, q) tibble::tibble(
  host = host, part = c("Pten", "Cdh1", "expected", "observed"),
  val = c(pten, cdh1, pten + cdh1, obs), synv = synv, q = q)
dec <- dplyr::bind_rows(
  mk("Cx3cr1+ depletion", syn$SingleA_Cx,  syn$SingleB_Cx,  syn$Effect_Cx,  syn$Synergy_Cx,  comb$q_dGI_Cx),
  mk("Mrc1+ depletion",   syn$SingleA_Mrc, syn$SingleB_Mrc, syn$Effect_Mrc, syn$Synergy_Mrc, comb$q_dGI_Mrc)) |>
  dplyr::mutate(
    host = factor(host, levels = names(HC)),
    part = factor(part, levels = c("Pten", "Cdh1", "expected", "observed")),
    fillc = dplyr::case_when(part == "observed" ~ HC[as.character(host)],
                             part == "expected" ~ "gray55", TRUE ~ "gray80"))
gap <- dec |> dplyr::group_by(host) |>
  dplyr::summarize(exp = val[part == "expected"], obs = val[part == "observed"],
                   synv = synv[1], q = q[1], .groups = "drop")

p <- ggplot(dec, aes(part, val, fill = fillc)) +
  geom_hline(yintercept = 0, color = "gray55", linewidth = 0.4) +
  geom_col(width = 0.68, color = "white", linewidth = 0.3) +
  scale_fill_identity() +
  geom_segment(data = gap, aes(x = 2.6, xend = 4.42, y = exp, yend = exp), inherit.aes = FALSE,
               linetype = "dashed", color = "gray40", linewidth = 0.35) +
  geom_segment(data = gap, aes(x = 4.42, xend = 4.42, y = exp, yend = obs), inherit.aes = FALSE,
               color = "gray15", linewidth = 0.5,
               arrow = grid::arrow(ends = "both", length = unit(2.2, "pt"), type = "closed")) +
  geom_text(data = gap, aes(x = 4.55, y = (exp + obs) / 2,
            label = sprintf("ΔGI %+.2f\n%s", synv, vapply(q, qlab, character(1)))),
            inherit.aes = FALSE, hjust = 0, size = 1.8, fontface = "bold", lineheight = 0.9, color = "gray15") +
  facet_wrap(~host) +
  scale_x_discrete(expand = expansion(add = c(0.6, 1.9))) +
  labs(title = "Pten×Cdh1: interaction is much stronger in the Mrc1⁺ cell depleted host",
       subtitle = sprintf("ΔLFC = LFC_Dep − LFC_WT · expected = Pten + Cdh1 · gap = ΔGI · between-host ΔΔGI = %.2f (q=%.3f)",
                          comb$ddGI, comb$ddGI_qvalue),
       x = NULL, y = "ΔLFC (LFC_Dep − LFC_WT)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 32, hjust = 1, size = 4.8),
        plot.subtitle = element_text(size = 4.8, color = "gray25"),
        panel.spacing = unit(6, "pt"))

save_pdf(p, fig("Fig5_PtenCdh1_hero.pdf"), 121, 74)
message(sprintf("[Fig5_PtenCdh1_hero] synergy Cx=%.2f (q=%.3f) Mrc=%.2f (q=%.3f) ΔΔGI=%.2f (q=%.3f)",
                syn$Synergy_Cx, comb$q_dGI_Cx, syn$Synergy_Mrc, comb$q_dGI_Mrc, comb$ddGI, comb$ddGI_qvalue))
log_session("Fig5_PtenCdh1_hero")
