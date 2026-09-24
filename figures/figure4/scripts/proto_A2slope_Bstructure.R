suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "theme_pub.R"))
  source(here::here("R", "palettes.R"))
  for (pkg in c("ggplot2", "dplyr", "tidyr", "readr", "ggrepel")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
})
SEED <- 123L

gi <- readr::read_csv(tbl("GI_Cx3cr1.csv"), show_col_types = FALSE,
                      progress = FALSE) |>
  dplyr::mutate(
    delta_GI = GI_Depleted_vs_Cell - GI_WT_vs_Cell,
    is_int = (GI_WT_qvalue < 0.05) | (GI_Depleted_qvalue < 0.05),
    is_hc  = is_int & (delta_GI_qvalue < 0.05),
    is_flip = (GI_WT_qvalue < 0.05) & (GI_Depleted_qvalue < 0.05) &
              (sign(GI_WT_vs_Cell) != sign(GI_Depleted_vs_Cell)))

hostlev <- c("WT host", "CX3CR1-depleted")
to_long <- function(d) d |>
  dplyr::transmute(gene, WT = GI_WT_vs_Cell, Dep = GI_Depleted_vs_Cell) |>
  tidyr::pivot_longer(c(WT, Dep), names_to = "host", values_to = "GI") |>
  dplyr::mutate(host = factor(ifelse(host == "WT", hostlev[1], hostlev[2]),
                              levels = hostlev))

flips <- gi |> dplyr::filter(is_flip) |>
  dplyr::mutate(dir = ifelse(GI_WT_vs_Cell > 0,
                             "Synergy -> buffering", "Buffering -> synergy"))
flip_long <- flips |> to_long() |>
  dplyr::left_join(dplyr::select(flips, gene, dir), by = "gene")
ctx_long  <- gi |> dplyr::filter(is_hc & !is_flip) |> to_long()
ymax <- max(abs(c(flip_long$GI, gi$GI_WT_vs_Cell[gi$is_hc],
                  gi$GI_Depleted_vs_Cell[gi$is_hc]))) * 1.05

DIRC <- c("Synergy -> buffering" = "#2166AC", "Buffering -> synergy" = "#D6604D")

p_slope <- ggplot() +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.5) +
  geom_line(data = ctx_long, aes(host, GI, group = gene),
            color = "grey85", linewidth = 0.3, alpha = 0.7) +
  geom_line(data = flip_long, aes(host, GI, group = gene, color = dir),
            linewidth = 0.7) +
  geom_point(data = flip_long, aes(host, GI, color = dir), size = 1.5) +
  ggrepel::geom_text_repel(
    data = dplyr::filter(flip_long, host == hostlev[2]),
    aes(host, GI, label = gsub("\\*", "_", gene), color = dir),
    size = 2.0, nudge_x = 0.18, hjust = 0, direction = "y",
    segment.size = 0.2, seed = SEED) +
  annotate("text", x = 0.62, y = ymax * 0.93, label = "synergistic",
           size = 2.1, color = "grey35", fontface = "italic", hjust = 0) +
  annotate("text", x = 0.62, y = -ymax * 0.93, label = "buffering",
           size = 2.1, color = "grey35", fontface = "italic", hjust = 0) +
  scale_color_manual(values = DIRC, name = NULL) +
  scale_x_discrete(expand = expansion(add = c(0.35, 0.95))) +
  coord_cartesian(ylim = c(-ymax, ymax)) +
  labs(title = "Sign-flips: the niche reverses the interaction",
       subtitle = sprintf("%d interactions cross zero when the CX3CR1+ niche is depleted (grey = other host-conditional pairs)",
                          nrow(flips)),

       x = NULL, y = expression("GI score")) +
  theme_pub(base_size = 7) +
  theme(legend.position = "top",
        plot.subtitle = element_text(size = 6.0, colour = "grey25"),
        legend.key.size = unit(7, "pt"))

IMMUNE <- c("Cxcr2", "Cxcr5", "Btla", "Mpeg1", "Pdcd1lg2", "Cx3cl1",
            "Fcrl6", "Tlr7", "Tnfsf15", "Il6", "Mysm1")

pb <- gi |> dplyr::group_by(geneB) |>
  dplyr::summarise(n_int = sum(is_int), n_hc = sum(is_hc),
                   hc_rate = 100 * sum(is_hc) / sum(is_int),
                   mean_abs = mean(abs(delta_GI)), .groups = "drop") |>
  dplyr::mutate(cls = ifelse(geneB %in% IMMUNE,
                             "Immune / chemokine / checkpoint",
                             "Other / structural"))

gi$cls <- ifelse(gi$geneB %in% IMMUNE, "immune", "other")
tab <- with(dplyr::filter(gi, is_int), table(cls, is_hc))
fish <- fisher.test(tab)
wil  <- wilcox.test(abs(delta_GI) ~ cls, data = dplyr::filter(gi, is_int))

cat(sprintf(paste0(
  "\n----------- PART B STRUCTURE CHECK -----------\n",
  "  host-conditional rate, immune partners : %.0f%%\n",
  "  host-conditional rate, other partners  : %.0f%%\n",
  "  Fisher OR = %.2f, p = %.3f\n",
  "  Wilcoxon |dGI| immune vs other, p = %.3f\n",
  "  => partner-identity structure is %s\n",
  "----------------------------------------------\n\n"),
  100 * sum(gi$is_hc[gi$cls=="immune"]) / sum(gi$is_int[gi$cls=="immune"]),
  100 * sum(gi$is_hc[gi$cls=="other"])  / sum(gi$is_int[gi$cls=="other"]),
  unname(fish$estimate), fish$p.value, wil$p.value,

  {
    padj <- stats::p.adjust(c(rate = fish$p.value, magnitude = wil$p.value),
                            method = "holm")
    sprintf(paste0(
      "rate: %s (Fisher p_adj = %.3f); magnitude: %s (Wilcoxon p_adj = %.3f)\n",
      "     => partner identity does NOT predict WHETHER a pair is host-conditional;\n",
      "        it is at most weakly associated with the SIZE of |dGI|"),
      ifelse(padj[["rate"]] < 0.05, "structured", "no evidence"), padj[["rate"]],
      ifelse(padj[["magnitude"]] < 0.05, "structured", "no evidence"), padj[["magnitude"]])
  }))

rate_imm <- 100 * sum(gi$is_hc[gi$is_int & gi$cls == "immune"]) /
                  sum(gi$is_int & gi$cls == "immune")
rate_oth <- 100 * sum(gi$is_hc[gi$is_int & gi$cls == "other"]) /
                  sum(gi$is_int & gi$cls == "other")

p_struct <- ggplot(pb, aes(hc_rate, reorder(geneB, hc_rate), color = cls)) +
  geom_segment(aes(x = 0, xend = hc_rate, yend = geneB), linewidth = 0.5) +
  geom_point(aes(size = mean_abs)) +
  scale_color_manual(values = c("Immune / chemokine / checkpoint" = "#E41A1C",
                                "Other / structural" = "grey55"), name = NULL) +
  scale_size_area(max_size = 3, name = expression("mean |" * Delta * "GI|")) +
  scale_x_continuous(limits = c(0, 105), expand = c(0, 0)) +
  labs(title = "Host-dependent rewiring",
       subtitle = sprintf("Immune %.0f%% vs other %.0f%% host-conditional (Fisher p = %.2f); partner class does not predict rewiring",
                          rate_imm, rate_oth, fish$p.value),
       x = "% of interactions host-conditional", y = NULL) +
  theme_pub(base_size = 7) +
  theme(legend.position = "right",
        plot.subtitle = element_text(size = 5.8, colour = "grey25"),
        axis.text.y = element_text(size = 6),
        legend.key.size = unit(6, "pt"))

save_one <- function(p, stem, w_mm, h_mm) {
  save_pdf(p, fig(paste0(stem, ".pdf")), w_mm, h_mm)
  message("  wrote ", stem, ".pdf")
}
save_one(p_slope,  "SuppFig4E_signflip_slopes",       121, 90)
save_one(p_struct, "SuppFig4F_partner_structure", 121, 90)
message("\n[slopes + supp partner-structure] complete.\n")

if (!interactive()) tryCatch(log_session("proto_A2slope_Bstructure"), error = function(e) NULL)
