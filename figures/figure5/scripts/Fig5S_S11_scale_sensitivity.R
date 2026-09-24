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

mk <- function(df, suf) df |> dplyr::transmute(gene, type,
        !!paste0("WT_",  suf) := mean_top_LFC_WT_vs_Cell_corrected,
        !!paste0("Dep_", suf) := mean_top_LFC_Depleted_vs_Cell_corrected)
glc <- readr::read_csv(tbl("geneLFC_Cx3cr1.csv"), show_col_types = FALSE)
glm <- readr::read_csv(tbl("geneLFC_Mrc1.csv"),   show_col_types = FALSE)

dual <- dplyr::inner_join(mk(glc, "Cx"), mk(glm, "Mrc") |> dplyr::select(-type),
                          by = "gene") |>
  dplyr::mutate(sign_flip = (WT_Cx * WT_Mrc) < 0, abs_wt_diff = abs(WT_Mrc - WT_Cx),
                dLFC_Cx = Dep_Cx - WT_Cx, dLFC_Mrc = Dep_Mrc - WT_Mrc) |>
  dplyr::filter(!(sign_flip & abs_wt_diff > 2), type == "dual") |>
  tidyr::separate(gene, into = c("A", "B"), sep = "\\*", remove = FALSE) |>
  dplyr::mutate(diff_niche_dep = dLFC_Mrc - dLFC_Cx,
                dnd_z = scale(dLFC_Mrc)[, 1] - scale(dLFC_Cx)[, 1])

anc <- function(d, v) dplyr::bind_rows(
  d |> dplyr::group_by(gene = A) |> dplyr::summarize(pos = "driver",  m = mean(.data[[v]]), .groups = "drop"),
  d |> dplyr::group_by(gene = B) |> dplyr::summarize(pos = "partner", m = mean(.data[[v]]), .groups = "drop"))
cmp <- dplyr::inner_join(anc(dual, "diff_niche_dep") |> dplyr::rename(raw = m),
                         anc(dual, "dnd_z") |> dplyr::rename(z = m),
                         by = c("gene", "pos"))
rho <- cor(cmp$raw, cmp$z, method = "spearman")

KEY <- c("Cdh1", "Brca2", "Sox10", "Pten", "Smad4", "Tgfbr2")
cmp <- dplyr::mutate(cmp,
  pos_lab = dplyr::recode(pos, driver = "driver (geneA)", partner = "partner (geneB)"),
  lab = ifelse(gene %in% KEY, gene, ""))

p <- ggplot(cmp, aes(raw, z)) +
  geom_hline(yintercept = 0, color = "gray80", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "gray80", linewidth = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              color = "gray70", linetype = "dashed", linewidth = 0.4) +
  geom_point(aes(color = pos_lab), size = 1.5, alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = lab, color = pos_lab), size = 1.8,
                           fontface = "bold", min.segment.length = 0,
                           box.padding = 0.4, max.overlaps = Inf,
                           segment.size = 0.22, seed = SEED, show.legend = FALSE) +
  scale_colour_manual(values = c("driver (geneA)" = "#4D4D4D",
                                 "partner (geneB)" = "#1B7837"), name = NULL) +
  labs(title = "Robust to per-screen scale",
       subtitle = sprintf("raw vs per-screen-standardized · ρ = %.2f", rho),
       x = "raw (LFC units)", y = "standardized (SD units)") +
  theme_pub() +
  theme(legend.position = "top", legend.text = element_text(size = 5.4),
        legend.key.size = unit(7, "pt"),
        plot.subtitle = element_text(size = 5.7, color = "gray30"))

save_pdf(p, fig("Fig5S_S11_scale_sensitivity.pdf"), 57, 62)
message(sprintf("[Fig5S_S11] scale sensitivity; Spearman rho(raw,z)=%.3f", rho))
log_session("Fig5S_S11_scale_sensitivity")
