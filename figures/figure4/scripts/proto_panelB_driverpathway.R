suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "theme_pub.R"))
  for (pkg in c("ggplot2", "dplyr", "readr")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
})

GROWTH <- c("Pten", "Tsc2", "Stk11", "Apc")
GENOME <- c("Brca1", "Brca2", "Chek2", "Rb1", "Smad4", "Tgfbr2")
GRP_COLORS <- c("Growth (PI3K / Wnt)" = "#377EB8",
                "Genome-integrity / TGFb" = "#E41A1C",
                "Other / unclassified" = "grey62")

gi <- readr::read_csv(tbl("GI_Cx3cr1.csv"), show_col_types = FALSE,
                      progress = FALSE) |>
  dplyr::mutate(
    delta_GI = GI_Depleted_vs_Cell - GI_WT_vs_Cell,
    is_int = (GI_WT_qvalue < 0.05) | (GI_Depleted_qvalue < 0.05),
    is_hc  = is_int & (delta_GI_qvalue < 0.05),
    grp = dplyr::case_when(geneA %in% GROWTH ~ "Growth (PI3K / Wnt)",
                           geneA %in% GENOME ~ "Genome-integrity / TGFb",
                           TRUE ~ "Other / unclassified"))

pa <- gi |> dplyr::group_by(geneA, grp) |>
  dplyr::summarise(n_int = sum(is_int), n_hc = sum(is_hc),
                   rate = 100 * sum(is_hc) / sum(is_int),
                   mean_abs = mean(abs(delta_GI)), .groups = "drop")

g2 <- gi |> dplyr::filter(is_int, grp != "Other / unclassified") |>
  dplyr::mutate(labe = grp == "Genome-integrity / TGFb")
ft <- fisher.test(table(g2$labe, g2$is_hc))
rate_growth <- 100 * sum(g2$is_hc[!g2$labe]) / sum(!g2$labe)
rate_genome <- 100 * sum(g2$is_hc[ g2$labe]) / sum( g2$labe)

cat(sprintf("\n[partB] Growth %.0f%% vs Genome/TGFb %.0f%% | OR=%.1f p=%.2g\n\n",
            rate_growth, rate_genome, unname(ft$estimate), ft$p.value))

p <- ggplot(pa, aes(rate, reorder(geneA, rate), color = grp)) +
  geom_segment(aes(x = 0, xend = rate, yend = geneA), linewidth = 0.5) +
  geom_point(aes(size = mean_abs)) +
  scale_color_manual(values = GRP_COLORS, name = NULL) +
  scale_size_area(max_size = 3.2, name = expression("mean |" * Delta * "GI|")) +
  scale_x_continuous(limits = c(0, 105), expand = c(0, 0)) +
  labs(
    title = "Host-dependent rewiring",
    subtitle = sprintf(
      "Genome-integrity / TGFb %.0f%% vs growth %.0f%% (OR=%.1f, p<1e-4) — aggregate trend, not a per-pair predictor",
      rate_genome, rate_growth, unname(ft$estimate)),
    x = "% of interactions host-conditional", y = NULL) +
  theme_pub(base_size = 7) +
  theme(legend.position = "right",
        plot.subtitle = element_text(size = 5.8, colour = "grey25"),
        axis.text.y = element_text(size = 6.5),
        legend.key.size = unit(7, "pt"),
        legend.spacing.y = unit(1, "pt"))

pten_y <- which(levels(reorder(pa$geneA, pa$rate)) == "Pten")
p <- p + annotate("text", x = 48, y = 2.6, hjust = 0, size = 1.8,
                  fontface = "italic", color = "grey30",
                  label = "Pten: niche-independent on average,\nyet Pten_Cx3cl1 has the largest dGI\nof any PTEN interaction")

save_pdf(p, fig("Fig4E_driver_pathway.pdf"), 121, 90)
message("  wrote Fig4E_driver_pathway.pdf\n")

if (!interactive()) tryCatch(log_session("proto_panelB_driverpathway"), error = function(e) NULL)
