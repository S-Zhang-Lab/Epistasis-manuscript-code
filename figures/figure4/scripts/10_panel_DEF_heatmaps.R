suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "theme_pub.R"))
  source(here::here("R", "palettes.R"))
  source(here::here("R", "helpers.R"))
  for (pkg in c("ggplot2", "dplyr", "tidyr", "tibble", "patchwork",
                "scales", "readr")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
})

HEATMAP_Q <- 0.97

message("\n=== 10_panel_DEF_heatmaps ===\n")

gi_raw <- read_csv(tbl("GI_Cx3cr1.csv"), show_col_types = FALSE,
                   progress = FALSE)

has_stats <- all(c("GI_WT_qvalue", "GI_Depleted_qvalue", "delta_GI_qvalue")
                 %in% names(gi_raw))
if (has_stats) {
  message("[10_panel_DEF] q-value columns present -- overlaying FDR asterisks")
} else {
  message("[10_panel_DEF] q-value columns absent -- skipping FDR overlay ",
          "(run scripts/05_stats.R for significance markers)")
}

stat_q <- function(col_q) if (has_stats) col_q else NA_character_

gi <- gi_raw %>%
  select(geneA, geneB,
         GI_WT  = GI_WT_vs_Cell,
         GI_Dep = GI_Depleted_vs_Cell,
         dplyr::any_of(c("GI_WT_qvalue",
                         "GI_Depleted_qvalue",
                         "delta_GI_qvalue"))) %>%
  mutate(delta_GI = GI_Dep - GI_WT)

q_to_star <- function(q) {
  ifelse(is.na(q), "",
         ifelse(q < 0.05, "**",
                ifelse(q < 0.10, "*", "")))
}
if (has_stats) {
  gi <- gi %>%
    mutate(star_WT    = q_to_star(GI_WT_qvalue),
           star_Dep   = q_to_star(GI_Depleted_qvalue),
           star_delta = q_to_star(delta_GI_qvalue))
} else {
  gi <- gi %>%
    mutate(star_WT = "", star_Dep = "", star_delta = "")
}

delta_mat <- gi %>%
  select(geneA, geneB, delta_GI) %>%
  pivot_wider(names_from = geneB, values_from = delta_GI) %>%
  column_to_rownames("geneA") %>%
  as.matrix()

row_clust <- hclust(dist(delta_mat,    method = "euclidean"), method = "complete")
col_clust <- hclust(dist(t(delta_mat), method = "euclidean"), method = "complete")
geneA_order <- rownames(delta_mat)[row_clust$order]
geneB_order <- colnames(delta_mat)[col_clust$order]

gi <- gi %>%
  mutate(geneA = factor(geneA, levels = geneA_order),
         geneB = factor(geneB, levels = geneB_order))

lim_rb <- sym_lim(c(gi$GI_WT, gi$GI_Dep), q = HEATMAP_Q)
lim_pg <- sym_lim(gi$delta_GI,            q = HEATMAP_Q)

theme_hm <- function(axis_text_size = 7) {
  theme_minimal(base_size = 9) +
    theme(
      axis.text.x  = element_text(size = axis_text_size, angle = 45,
                                  hjust = 1, color = "black"),
      axis.text.y  = element_text(size = axis_text_size, color = "black"),
      axis.title.x = element_text(size = 9, face = "bold",
                                  margin = margin(t = 4)),
      axis.title.y = element_text(size = 9, face = "bold",
                                  margin = margin(r = 4)),
      plot.title   = element_text(size = 9.5, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 8),
      legend.text  = element_text(size = 7.5),
      legend.key.width  = unit(0.9,  "cm"),
      legend.key.height = unit(0.28, "cm"),
      panel.grid   = element_blank(),
      plot.margin  = margin(4, 4, 4, 4)
    )
}
scale_rb <- function(lim) scale_fill_gradientn(
  colors = DIV_RB_COLORS,
  values = rescale(DIV_ANCHOR_PROPS * lim),
  limits = c(-lim, lim), oob = squish, name = "GI",
  guide  = guide_colorbar(title.position = "top", title.hjust = 0.5))

scale_pg <- function(lim) scale_fill_gradientn(
  colors = DIV_PG_COLORS,
  values = rescale(DIV_ANCHOR_PROPS * lim),
  limits = c(-lim, lim), oob = squish, name = expression(Delta*"GI"),
  guide  = guide_colorbar(title.position = "top", title.hjust = 0.5))

make_hm <- function(df, value_col, star_col, title_str, scale_fn, lim,
                    show_x = TRUE, show_y = TRUE, swap_xy = FALSE) {
  if (swap_xy) {
    p <- ggplot(df, aes(x = geneA, y = geneB, fill = .data[[value_col]])) +
      geom_tile(color = "white", linewidth = 0.25) +
      scale_fn(lim) +
      labs(title = title_str,
           x = if (show_x) "geneA (before *)" else NULL,
           y = if (show_y) "geneB (after *)"  else NULL) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      theme_hm()
  } else {
    p <- ggplot(df, aes(x = geneB, y = geneA, fill = .data[[value_col]])) +
      geom_tile(color = "white", linewidth = 0.25) +
      scale_fn(lim) +
      labs(title = title_str,
           x = if (show_x) "geneB (after *)"  else NULL,
           y = if (show_y) "geneA (before *)" else NULL) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      theme_hm()
  }
  if (any(nzchar(df[[star_col]]))) {
    p <- p + geom_text(aes(label = .data[[star_col]]),
                       size = 2.2, color = "black", lineheight = 0.7,
                       vjust = 0.65)
  }
  if (!show_x) p <- p + theme(axis.text.x  = element_blank(),
                              axis.ticks.x = element_blank())
  if (!show_y) p <- p + theme(axis.text.y  = element_blank(),
                              axis.ticks.y = element_blank())

  p + coord_fixed(ratio = 1)
}

p_C <- make_hm(gi, "delta_GI", "star_delta",
               expression(Delta*"GI (CX3CR1 depleted - WT)"),
               scale_pg, lim_pg,
               show_x = TRUE, show_y = TRUE, swap_xy = TRUE)

p_SuppA1 <- make_hm(gi, "GI_WT",  "star_WT",
                    "WT",
                    scale_rb, lim_rb,
                    show_x = TRUE, show_y = TRUE, swap_xy = TRUE)
p_SuppA2 <- make_hm(gi, "GI_Dep", "star_Dep",
                    "CX3CR1 Depleted",
                    scale_rb, lim_rb,
                    show_x = TRUE, show_y = TRUE, swap_xy = TRUE)

PANEL_W_MM <- 140
PANEL_H_MM <- 140

save_pdf(p_C,      fig("Fig4B_deltaGI_heatmap.pdf"),       PANEL_W_MM, PANEL_H_MM)
save_pdf(p_SuppA1, fig("SuppFig4C_GI_heatmap_Cx3cr1_WT.pdf"),    PANEL_W_MM, PANEL_H_MM)
save_pdf(p_SuppA2, fig("SuppFig4D_GI_heatmap_Cx3cr1_Dep.pdf"),   PANEL_W_MM, PANEL_H_MM)

if (!interactive()) log_session("10_panel_DEF_heatmaps")

invisible(list(C = p_C, SuppA1 = p_SuppA1, SuppA2 = p_SuppA2, gi = gi))
