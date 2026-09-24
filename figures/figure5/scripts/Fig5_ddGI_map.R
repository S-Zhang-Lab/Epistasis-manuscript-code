suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "theme_pub.R"))
  source(here::here("R", "palettes.R"))
  source(here::here("R", "helpers.R"))
  for (pkg in c("ggplot2", "dplyr", "tidyr", "tibble", "scales", "readr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
})

HEATMAP_Q <- 0.97
sym_lim <- function(vals, q = 0.97) stats::quantile(abs(vals), q, na.rm = TRUE)
q_to_star <- function(q) ifelse(is.na(q), "",
                          ifelse(q < 0.05, "**", ifelse(q < 0.10, "*", "")))

comb <- readr::read_csv(tbl("ddGI_combined.csv"), show_col_types = FALSE,
                        progress = FALSE) |>
  dplyr::mutate(star_dd  = q_to_star(ddGI_qvalue),
                star_Cx  = q_to_star(q_dGI_Cx),
                star_Mrc = q_to_star(q_dGI_Mrc))

ddmat <- comb |>
  dplyr::select(geneA, geneB, ddGI) |>
  tidyr::pivot_wider(names_from = geneB, values_from = ddGI) |>
  tibble::column_to_rownames("geneA") |>
  as.matrix()
m0 <- ddmat; m0[is.na(m0)] <- 0
geneA_order <- rownames(ddmat)[hclust(dist(m0),    "complete")$order]
geneB_order <- colnames(ddmat)[hclust(dist(t(m0)), "complete")$order]

comb <- comb |>
  dplyr::mutate(geneA = factor(geneA, levels = geneA_order),
                geneB = factor(geneB, levels = geneB_order))

lim_pb <- sym_lim(comb$ddGI, HEATMAP_Q)
lim_rb <- sym_lim(c(comb$dGI_Cx, comb$dGI_Mrc), HEATMAP_Q)

scale_pb <- function(lim) scale_fill_gradientn(
  colors = DIV_PB_COLORS, values = scales::rescale(DIV_ANCHOR_PROPS * lim),
  limits = c(-lim, lim), oob = scales::squish, name = "ΔΔGI",
  guide = guide_colorbar(title.position = "top", title.hjust = 0.5))
scale_rb <- function(lim, nm) scale_fill_gradientn(
  colors = DIV_RB_COLORS, values = scales::rescale(DIV_ANCHOR_PROPS * lim),
  limits = c(-lim, lim), oob = scales::squish, name = nm,
  guide = guide_colorbar(title.position = "top", title.hjust = 0.5))

theme_hm <- function() theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(axis.text.x = element_text(size = 5.5, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 5.5, color = "black"),
        axis.title  = element_text(size = 7, face = "bold"),
        plot.title  = element_text(size = 8, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 6, color = "gray30", hjust = 0.5),
        legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
        legend.key.width = unit(0.7, "cm"), legend.key.height = unit(0.22, "cm"),
        panel.grid = element_blank(), plot.margin = margin(3, 3, 3, 3))

make_hm <- function(value_col, star_col, title_str, subtitle_str, scale_obj) {
  ggplot(comb, aes(x = geneA, y = geneB, fill = .data[[value_col]])) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = .data[[star_col]]), size = 1.9, color = "black",
              lineheight = 0.7, vjust = 0.7) +
    scale_obj +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(title = title_str, subtitle = subtitle_str,
         x = "Tumor driver (geneA)", y = "Niche / immune partner (geneB)") +
    coord_fixed(ratio = 1) + theme_hm()
}

p_main <- make_hm("ddGI", "star_dd",
                  "Host-specific interaction landscape (ΔΔGI = Cx − Mrc)",
                  "pink = stronger in Cx3cr1⁺ cell depleted host · blue = stronger in Mrc1⁺ cell depleted host · ** ΔΔGI q<0.05",
                  scale_pb(lim_pb))

p_main <- p_main +
  geom_tile(data = dplyr::filter(comb, gene == "Pten*Cdh1"),
            aes(x = geneA, y = geneB), fill = NA, color = "black",
            linewidth = 0.8, inherit.aes = FALSE)
p_cx <- make_hm("dGI_Cx", "star_Cx",
                "ΔGI under Cx3cr1+ depletion (Supp.)",
                "per-screen niche-conditional interaction", scale_rb(lim_rb, "ΔGI_Cx"))
p_mrc <- make_hm("dGI_Mrc", "star_Mrc",
                 "ΔGI under Mrc1+ depletion (Supp.)",
                 "per-screen niche-conditional interaction", scale_rb(lim_rb, "ΔGI_Mrc"))

save_pdf(p_main, fig("Fig5_ddGI_map.pdf"),         121, 130)
save_pdf(p_cx,   fig("Fig5S_S4_dGI_Cx.pdf"),  121, 130)
save_pdf(p_mrc,  fig("Fig5S_S5_dGI_Mrc.pdf"), 121, 130)
message("[proto_C] wrote ΔΔGI heatmap + 2 per-screen supplements")
log_session("Fig5_ddGI_map")
