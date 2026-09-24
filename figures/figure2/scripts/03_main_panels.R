required_pkgs <- c(
  "ggplot2", "ggrepel", "dplyr",
  "readxl", "pheatmap",
  "igraph", "ggraph", "tidygraph",
  "stringr", "tidyr"
)
missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

library(ggplot2)
library(ggrepel)
library(dplyr)
library(readxl)
library(pheatmap)
library(igraph)
library(ggraph)
library(tidygraph)
library(stringr)
library(tidyr)

source(here::here("R", "paths.R"))
source(here::here("R", "palettes.R"))
source(here::here("R", "helpers.R"))

source(here::here("R", "theme_pub.R"))

OUTPUT_DIR   <- OUT_TABLES
FIG_DIR      <- OUT_FIGS

COL_ENR <- "#d62728"
COL_DEP <- "#1f77b4"
COL_NT  <- "#ff7f0e"
COL_BG  <- "#aaaaaa"

MODULE_COLORS <- c("0" = "#E41A1C", "1" = "#377EB8", "2" = "#4DAF4A")

GI_STRONG_Z_THRESHOLD <- 1.5

message("\n=== Figure 2 Consolidated R Script ===")
message("Repo root  : ", REPO_ROOT)
message("Output dir : ", FIG_DIR, "\n")

message("--- Fig2C: Sample PCA ---")

PCA_DAY0_MIN <- 5

group_order <- c("Cell_Injection", "In_Vitro_Day14",
                 "In_Vivo_Day3", "In_Vivo_Day18")

group_info <- c(
  setNames(rep("Cell_Injection", 3),
           c("Cell_Day0", "Cell_Day0_2", "Cell_Day0_3")),
  setNames(rep("In_Vitro_Day14", 5), paste0("Cell_Day14_", 1:5)),
  setNames(rep("In_Vivo_Day3",   6), paste0("Mice_Day3_",  1:6)),
  setNames(rep("In_Vivo_Day18",  6), paste0("Mice_Day18_", 1:6))
)
DAY0_REPS <- c("Cell_Day0", "Cell_Day0_2", "Cell_Day0_3")

pca_raw <- read.csv(
  raw("R1_construct_counts.csv"),
  header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
)

sample_cols <- intersect(names(group_info), colnames(pca_raw))
missing_s   <- setdiff(names(group_info), colnames(pca_raw))
if (length(missing_s) > 0)
  stop("Samples defined but not found in the count matrix: ",
       paste(missing_s, collapse = ", "))
stopifnot(length(sample_cols) == 20L)

mat <- as.matrix(pca_raw[, sample_cols, drop = FALSE])
storage.mode(mat) <- "double"
mat[is.na(mat)] <- 0

mat <- sweep(mat, 2, colSums(mat), "/") * 1e7

gate <- rowSums(pca_raw[, DAY0_REPS, drop = FALSE]) >= PCA_DAY0_MIN
mat  <- mat[gate, , drop = FALSE]
message(sprintf("  constructs with pooled Day-0 >= %d: %d of %d",
                PCA_DAY0_MIN, nrow(mat), nrow(pca_raw)))

mat <- log2(mat + 1)
keep <- apply(mat, 1, function(x) sd(x) > 0)
mat  <- mat[keep, , drop = FALSE]

pca_input <- t(mat)

pca_res  <- prcomp(pca_input, center = TRUE, scale. = FALSE)

.sc <- pca_res$x
if (.sc["Cell_Day0", "PC1"] < 0) { .sc[, "PC1"] <- -.sc[, "PC1"] }
.d18 <- paste0("Mice_Day18_", 1:6)
if (mean(.sc[.d18, "PC2"]) < 0) { .sc[, "PC2"] <- -.sc[, "PC2"] }
pca_res$x <- .sc

pca_df   <- as.data.frame(pca_res$x)
pca_df$Sample <- rownames(pca_df)
pca_df$Group  <- factor(group_info[pca_df$Sample], levels = group_order)
pct_var  <- round((pca_res$sdev^2) / sum(pca_res$sdev^2) * 100, 2)

plot_pca <- function(df, pc_x = 1, pc_y = 2) {
  xc <- paste0("PC", pc_x); yc <- paste0("PC", pc_y)
  ellipse_df <- df |>
    group_by(Group) |>
    filter(n() >= 3, qr(cbind(.data[[xc]], .data[[yc]]))$rank >= 2) |>
    ungroup()
  ggplot(df, aes(x = .data[[xc]], y = .data[[yc]])) +
    geom_point(aes(color = Group, shape = Group), size = SZ_POINT) +

    stat_ellipse(data = ellipse_df, aes(fill = Group), geom = "polygon",
                 type = "norm", level = 0.80, alpha = 0.15,
                 show.legend = FALSE) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    scale_shape_manual(values = c(
      "Cell_Injection" = 16, "In_Vitro_Day14" = 16,
      "In_Vivo_Day3" = 17, "In_Vivo_Day18" = 15
    ), drop = FALSE) +
    xlab(paste0(xc, " (", pct_var[pc_x], "%)")) +
    ylab(paste0(yc, " (", pct_var[pc_y], "%)")) +
    ggtitle("PCA \u2014 sgRNA library representation") +
    theme_panel() +
    theme(panel.grid.major.x = element_line(colour = "#eeeeee", linewidth = 0.2))
}

p2b <- plot_pca(pca_df)

save_panel(p2b, file.path(FIG_DIR, "Fig2C_PCA.pdf"))

message("  Saved Fig2C_PCA.pdf/.png")

message("--- Fig2D: Rank Plot ---")

FDR_THRESH <- 0.05
N_LABEL    <- 10

make_rank_plot <- function(tbl, y_lab, title, subtitle, out_dir, out_stem) {
  gp <- read.csv(file.path(OUTPUT_DIR, tbl), stringsAsFactors = FALSE) |>
    arrange(desc(weighted_mean_lfc)) |>
    mutate(
      rank        = seq_len(n()),
      label_clean = pair_label(gene_pair),
      is_nt       = pair_type == "nt_nt",
      category    = case_when(
        is_nt                                     ~ "NT control",
        fdr < FDR_THRESH & weighted_mean_lfc > 0 ~ "Enriched",
        fdr < FDR_THRESH & weighted_mean_lfc < 0 ~ "Depleted",
        TRUE                                      ~ "Not significant"
      ),
      category = factor(category, levels = c("Enriched", "Depleted",
                                             "NT control", "Not significant"))
    )
  n_enr <- sum(gp$category == "Enriched")
  n_dep <- sum(gp$category == "Depleted")

  top_enr <- gp |> filter(category == "Enriched", pair_type == "dual_KO") |>
    slice_head(n = N_LABEL)
  nt_row  <- gp |> filter(is_nt) |> slice_head(n = 1)
  n_total <- nrow(gp)

  col_map   <- c("Enriched" = COL_ENR, "Depleted" = COL_DEP,
                 "NT control" = COL_NT, "Not significant" = COL_BG)
  label_map <- c(
    "Enriched"        = paste0("Enriched  FDR<", FDR_THRESH, "  (n=", n_enr, ")"),
    "Depleted"        = paste0("Depleted  FDR<", FDR_THRESH, "  (n=", n_dep, ")"),
    "NT control"      = paste0("Nt", PAIR_SEP, "Nt control"),
    "Not significant" = "Not significant")

  p <- ggplot(gp, aes(rank, weighted_mean_lfc, color = category)) +
    geom_hline(yintercept = 0, color = "black", linewidth = LW_RULE) +
    geom_point(data = filter(gp, category == "Not significant"),
               size = SZ_POINT_SM, alpha = 0.45, shape = 16) +
    geom_point(data = filter(gp, category != "Not significant"),
               size = SZ_POINT_HL, alpha = 0.8, shape = 16) +
    {if (nrow(top_enr) > 0)
      geom_text_repel(
        data = top_enr, aes(label = label_clean),
        color = COL_ENR, size = pt(FS_LABEL), fontface = "plain",
        nudge_x = n_total * 0.30, direction = "y", hjust = 0,
        segment.color = COL_ENR, segment.size = 0.2, segment.alpha = 0.6,
        min.segment.length = 0, box.padding = 0.14, point.padding = 0.06,
        max.overlaps = Inf, force = 2, seed = 42, show.legend = FALSE,

        ylim = c(max(gp$weighted_mean_lfc, na.rm = TRUE) * 0.06, NA))} +
    {if (nrow(nt_row) > 0)
      geom_text_repel(
        data = nt_row, aes(label = label_clean),
        color = COL_NT, size = pt(FS_LABEL), fontface = "plain",
        nudge_x = n_total * 0.10,
        nudge_y = max(gp$weighted_mean_lfc, na.rm = TRUE) * 0.20,
        segment.color = COL_NT, segment.size = 0.25, segment.alpha = 0.85,
        arrow = arrow(length = unit(2, "pt"), type = "closed"),
        min.segment.length = 0, max.overlaps = Inf, seed = 42,
        show.legend = FALSE)} +
    scale_color_manual(values = col_map, labels = label_map, name = NULL,
                       drop = FALSE,
                       guide = guide_legend(override.aes = list(
                         size = 1.1, alpha = 1, shape = 15))) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +

    labs(x = expression("Rank  (by weighted mean log"[2]*"FC)"),
         y = y_lab, title = title, subtitle = subtitle) +
    theme_panel() +
    theme(
      legend.position  = "inside",
      legend.position.inside = c(0.02, 0.02),
      legend.justification = c("left", "bottom"),
      legend.background = element_rect(fill = "white", color = "#cccccc",
                                       linewidth = LW_THIN),
      legend.key.size  = unit(4.5, "pt")
    )
  save_panel(p, file.path(out_dir, paste0(out_stem, ".pdf")))
  message(sprintf("  Saved %s.pdf/.png | %d enriched, %d depleted (FDR<%.2f)",
                  out_stem, n_enr, n_dep, FDR_THRESH))
  invisible(p)
}

p2c <- make_rank_plot(
  "R1_enrichment_day0_invivo.csv",
  expression("Weighted mean log"[2]*"FC  (Day 18 / Day 0)"),
  "Rank Plot \u2014 All Gene Pairs, in vivo",
  paste0("Nt", PAIR_SEP, "Nt-corrected; FDR<0.05 (BH, random-effects Wald)"),
  FIG_DIR, "Fig2D_rank_plot")

if (file.exists(file.path(OUTPUT_DIR, "R1_enrichment_day0_invitro.csv"))) {
  make_rank_plot(
    "R1_enrichment_day0_invitro.csv",
    expression("Weighted mean log"[2]*"FC  (Day 14 / Day 0)"),
    "Rank Plot \u2014 All Gene Pairs, in vitro control",
    "Same library and baseline, cultured rather than injected",
    OUT_SUPP, "SuppFig2N_invitro_rank_plot")
}

message("--- Fig2G: Volcano Plot ---")

gs_df <- read.table(
  derived("R2_Day18.gene_summary.txt"),
  header = TRUE, sep = "\t", check.names = FALSE
)

gs_df <- gs_df |>
  rename(
    label    = id,
    lfc      = `pos|lfc`,
    pos_pval = `pos|p-value`,
    neg_pval = `neg|p-value`
  ) |>
  mutate(
    pval        = ifelse(lfc >= 0, pos_pval, neg_pval),
    neg_log10_p = -log10(pval)
  ) |>
  filter(!is.na(pval))

sig_thr <- -log10(0.05)

gs_df <- gs_df |>
  mutate(color_group = case_when(
    neg_log10_p >= sig_thr & lfc > 0  ~ "enriched",
    neg_log10_p >= sig_thr & lfc < 0  ~ "depleted",
    TRUE                               ~ "ns"
  ))

to_label_vol <- rbind(
  head(gs_df[gs_df$color_group == "depleted", ][
    order(gs_df$pval[gs_df$color_group == "depleted"]), , drop = FALSE], 3),
  head(gs_df[gs_df$color_group == "enriched", ][
    order(gs_df$pval[gs_df$color_group == "enriched"]), , drop = FALSE], 8)
)
to_label_vol <- to_label_vol[!duplicated(to_label_vol$label), , drop = FALSE]

p2i <- ggplot(gs_df, aes(lfc, neg_log10_p, color = color_group)) +
  geom_point(size = SZ_POINT * 0.8, alpha = 0.85) +
  geom_hline(yintercept = sig_thr, linetype = "dashed",
             color = "gray50", linewidth = LW_THIN) +

  {set.seed(42); geom_text_repel(
    data = to_label_vol, aes(label = pair_label(label)),
    size = pt(FS_LABEL), fontface = "plain", seed = 42,
    segment.size = 0.2, max.overlaps = 20, box.padding = 0.16
  )} +
  annotate("text", x = Inf, y = sig_thr + 0.12,
           label = "p = 0.05", hjust = 1.05, color = "gray50", size = pt(FS_LABEL)) +
  annotate("segment",
           x = min(gs_df$lfc) * 0.3, xend = max(gs_df$lfc) * 1.05,
           y = -0.5, yend = -0.5,
           arrow = arrow(length = unit(2.5, "pt"), type = "closed"),
           color = "#CC3300", linewidth = LW_RULE) +
  annotate("text",
           x    = (min(gs_df$lfc) * 0.3 + max(gs_df$lfc) * 1.05) / 2,
           y    = -0.85,
           label = "enriched crRNA pairs in overt lung mets",
           color = "#CC3300", size = pt(FS_LABEL), fontface = "plain") +
  scale_color_manual(
    values = c("enriched" = "#D84A1B", "depleted" = "#3157D5", "ns" = "#A6A6A6"),
    guide = "none"
  ) +
  labs(
    title = "Volcano \u2014 crRNA frequency change, 2nd round",
    x     = expression("Log"[2]*" Fold Change (LFC)"),
    y     = expression("\u2212Log"[10]*"(p-value)")
  ) +
  theme_panel() +
  theme(axis.title.x = element_text(margin = margin(t = 9)),
        plot.margin  = margin(3, 3, 3, 3)) +
  coord_cartesian(clip = "off")

save_panel(p2i, file.path(FIG_DIR, "Fig2G_volcano.pdf"))

message(sprintf("  Saved Fig2G_volcano.pdf/.png | enriched=%d, depleted=%d, ns=%d",
                sum(gs_df$color_group == "enriched"),
                sum(gs_df$color_group == "depleted"),
                sum(gs_df$color_group == "ns")))

.is_nt <- function(x) !is.na(x) & tolower(x) %in% c("nt", "ntc")

message("--- Fig2H: GI Z-score Heatmap ---")

gi_df <- read.csv(file.path(OUTPUT_DIR, "R2_GI_derived.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
gi_df <- as.data.frame(gi_df)
gi_df <- gi_df[!.is_nt(gi_df$Gene1) & !.is_nt(gi_df$Gene2) & !is.na(gi_df$Gene2), ]

gi_df$GI_Zscore <- gi_df$GI_Z_Score
global_sd      <- sd(gi_df$GI, na.rm = TRUE)

gi_out <- gi_df[!is.na(gi_df$GI_Zscore), c("Gene", "Gene1", "Gene2",
                                              "Count", "GI", "GI_Zscore")]
write.csv(gi_out[order(gi_out$GI_Zscore), ],
          file.path(OUTPUT_DIR, "GI_Zscore_values.csv"), row.names = FALSE)

mat_wide <- reshape(gi_df[, c("Gene1", "Gene2", "GI_Zscore")],
                    idvar = "Gene2", timevar = "Gene1", direction = "wide")
rownames(mat_wide) <- mat_wide$Gene2
mat_wide$Gene2 <- NULL
colnames(mat_wide) <- sub("GI_Zscore\\.", "", colnames(mat_wide))
gi_mat <- as.matrix(mat_wide)
storage.mode(gi_mat) <- "numeric"

HEAT_LIM <- 3
gi_mat_plot <- pmax(pmin(gi_mat, HEAT_LIM), -HEAT_LIM)
breaks_heat <- c(
  seq(-HEAT_LIM, -0.3, length.out = 61),
  seq(-0.3, 0.3, length.out = 13)[-1],
  seq(0.3, HEAT_LIM, length.out = 61)[-1]
)
colors_heat <- c(
  colorRampPalette(c("#2166AC", "#FFFFFF"))(60),
  rep("#FFFFFF", 12),
  colorRampPalette(c("#FFFFFF", "#B2182B"))(60)
)
leg_breaks  <- c(-3, -2, -1, 0, 1, 2, 3)

gi_mat_impute <- gi_mat
gi_mat_impute[is.na(gi_mat_impute)] <- 0

hc_row <- hclust(dist(gi_mat_impute),      method = "average")
hc_col <- hclust(dist(t(gi_mat_impute)),   method = "average")

heat_args <- list(
  gi_mat_plot,
  color             = colors_heat,
  breaks            = breaks_heat,
  na_col            = "#D6D6D6",
  cluster_rows      = hc_row,
  cluster_cols      = hc_col,
  border_color      = "#C7C7C7",

  cellwidth         = PANEL_PT / ncol(gi_mat_plot),
  cellheight        = PANEL_PT / nrow(gi_mat_plot),
  fontsize          = MIN_FONT_PT,
  fontsize_row      = MIN_FONT_PT,
  fontsize_col      = MIN_FONT_PT,
  fontfamily        = PANEL_FAMILY,
  angle_col         = 270,
  treeheight_row    = 14, treeheight_col = 14,
  legend_breaks     = leg_breaks,
  legend_labels     = as.character(leg_breaks)
)

pdf(nullfile(), width = 8, height = 8, family = PANEL_FAMILY,
    encoding = PDF_ENCODING, useDingbats = FALSE)
ph  <- do.call(pheatmap, c(heat_args, list(filename = NA, silent = TRUE)))
.hw <- sum(grid::convertWidth(ph$gtable$widths,   "in", valueOnly = TRUE))
.hh <- sum(grid::convertHeight(ph$gtable$heights, "in", valueOnly = TRUE))
invisible(dev.off())

pdf(file.path(FIG_DIR, "Fig2H_GI_heatmap.pdf"), width = .hw, height = .hh,
    family = PANEL_FAMILY, encoding = PDF_ENCODING, useDingbats = FALSE)
grid::grid.newpage(); grid::grid.draw(ph$gtable)
invisible(dev.off())

png(file.path(FIG_DIR, "Fig2H_GI_heatmap.png"), width = .hw, height = .hh,
    units = "in", res = 600, type = "cairo", family = PANEL_FAMILY)
grid::grid.newpage(); grid::grid.draw(ph$gtable)
invisible(dev.off())

.record_panel_size(file.path(FIG_DIR, "Fig2H_GI_heatmap"),
                   PANEL_PT, PANEL_PT, .hw, .hh, 1L)

message("  Saved Fig2H_GI_heatmap.pdf/.png")

message("--- Section 5: Fig2H-derived R2 GI graph prep ---")

r2_gene1_modules <- c(
  Pten = "Growth_Control", Apc = "Growth_Control", Rb1 = "Growth_Control",
  Smad4 = "Growth_Control", Stk11 = "Growth_Control", Tsc2 = "Growth_Control",
  Brca1 = "DNA_Repair", Brca2 = "DNA_Repair", Chek2 = "DNA_Repair",
  Calml3 = "Barrier_ECM", Serpinb5 = "Barrier_ECM",
  Sox10 = "Cell_Identity", Krt15 = "Cell_Identity", Clca2 = "Cell_Identity",
  Tgfbr2 = "Unknown", Kit = "Unknown", Chit1 = "Unknown"
)
r2_gene2_modules <- c(
  Btla = "Immune_Evasion", Cxcr2 = "Immune_Evasion", Cxcr5 = "Immune_Evasion",
  Fcrl6 = "Immune_Evasion", Il6 = "Immune_Evasion", Mpeg1 = "Immune_Evasion",
  Pdcd1lg2 = "Immune_Evasion", Tlr7 = "Immune_Evasion",
  Tnfsf15 = "Immune_Evasion", Cx3cl1 = "Barrier_ECM", Mfap5 = "Barrier_ECM",
  Muc16 = "Barrier_ECM", Pla2g2a = "Barrier_ECM", Pla2r1 = "Barrier_ECM",
  Cdh1 = "Cell_Identity", Hgf = "Growth_Signaling", Mysm1 = "Unknown"
)
r2_modules <- c(r2_gene1_modules, r2_gene2_modules)
r2_modules <- setNames(unname(r2_modules), toupper(names(r2_modules)))

gi2_net <- read.csv(file.path(OUTPUT_DIR, "R2_GI_derived.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
gi2_gene_rows <- gi2_net |>
  filter(!.is_nt(Gene1), !.is_nt(Gene2), !is.na(Gene2)) |>
  mutate(gene_A = toupper(Gene1), gene_B = toupper(Gene2))
r2_genes <- sort(unique(c(gi2_gene_rows$gene_A, gi2_gene_rows$gene_B)))

gi2_global_sd <- sd(gi2_gene_rows$GI, na.rm = TRUE)
gi2_global_mean <- mean(gi2_gene_rows$GI, na.rm = TRUE)

gi2_edges <- gi2_gene_rows |>

  mutate(
    GI_Zscore = (GI - gi2_global_mean) / gi2_global_sd,
    abs_GI_Z = abs(GI_Zscore)
  )

strong_gi_partner_counts <- gi2_edges |>
  filter(abs_GI_Z > GI_STRONG_Z_THRESHOLD) |>
  select(gene_A, gene_B) |>
  pivot_longer(cols = everything(), names_to = "endpoint", values_to = "name") |>
  count(name, name = "strong_gi_partners")

r2_single_lfc <- bind_rows(
  gi2_net |> filter(!.is_nt(Gene1)) |>
    transmute(gene = toupper(Gene1), single_lfc = Gene1_NT_LFC_Mean),
  gi2_net |> filter(!.is_nt(Gene2)) |>
    transmute(gene = toupper(Gene2), single_lfc = Gene2_NT_LFC_Mean)
) |>
  filter(!is.na(single_lfc)) |>
  group_by(gene) |>
  summarise(single_lfc = first(single_lfc), .groups = "drop")

r2_vertices <- data.frame(name = r2_genes, stringsAsFactors = FALSE) |>
  left_join(r2_single_lfc, by = c("name" = "gene")) |>
  left_join(strong_gi_partner_counts, by = "name") |>
  mutate(
    strong_gi_partners = coalesce(strong_gi_partners, 0L),
    module_f = factor(ifelse(name %in% names(r2_modules), r2_modules[name], "Unknown"),
                      levels = names(R2_MODULE_COLORS))
  )

write.csv(
  r2_vertices |>
    transmute(
      gene = name,
      module = as.character(module_f),
      strong_gi_partners = strong_gi_partners,
      single_lfc = single_lfc
    ),
  file.path(OUTPUT_DIR, "Fig2I_node_strong_GI_partner_counts.csv"),
  row.names = FALSE
)

message(sprintf(
  "  Fig2H-derived R2 GI graph: %d genes, %d GI edges, %d strong edges (|GI Z| > %.1f)",
  length(r2_genes), nrow(gi2_edges), sum(gi2_edges$abs_GI_Z > GI_STRONG_Z_THRESHOLD),
  GI_STRONG_Z_THRESHOLD
))

message("--- Fig2I: Fig2H-derived strong-edge GI Network ---")

to_mouse_symbol <- function(x) {
  ifelse(is.na(x) | x == "",
         x,
         paste0(toupper(substr(x, 1, 1)),
                tolower(substring(x, 2))))
}

gi2_strong_edges <- gi2_edges |>
  filter(abs_GI_Z > GI_STRONG_Z_THRESHOLD)

strong_vertices <- r2_vertices |>
  filter(strong_gi_partners > 0)

write.csv(
  gi2_strong_edges |>
    transmute(gene_A, gene_B, GI, GI_Zscore, abs_GI_Z, Count),
  file.path(OUTPUT_DIR, "Fig2I_strong_GI_edges.csv"),
  row.names = FALSE
)

g_strong <- graph_from_data_frame(
  gi2_strong_edges[, c("gene_A", "gene_B", "GI", "GI_Zscore", "abs_GI_Z")],
  vertices = strong_vertices,
  directed = FALSE
)

tg_strong <- as_tbl_graph(g_strong) |>
  activate(nodes) |>
  mutate(
    module_f = factor(module_f, levels = names(R2_MODULE_COLORS))
  ) |>
  activate(edges) |>
  mutate(
    zscore = GI_Zscore,
    edge_weight_plot = abs_GI_Z
  )

set.seed(123)
p2l <- ggraph(tg_strong, layout = "fr", weights = abs_GI_Z) +
  geom_edge_link(aes(color = zscore, width = edge_weight_plot), alpha = 0.78) +
  scale_edge_color_gradient2(
    low = GI_BLUE, mid = "white", high = GI_RED, midpoint = 0, name = "GI Z"
  ) +
  scale_edge_width(range = c(0.25, 2.2), guide = "none") +

  geom_node_point(
    aes(fill = module_f),
    shape = 21, color = "white", size = 1.9, stroke = 0.2
  ) +

  geom_node_point(
    data = function(x) dplyr::filter(x, name == "PTEN"),
    shape = 21, fill = NA, color = "black", size = 2.4, stroke = 0.45,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = R2_MODULE_COLORS, name = "Functional\nmodule") +
  geom_node_text(
    aes(label = to_mouse_symbol(name),
        fontface = "plain"),
    repel = TRUE, size = pt(FS_LABEL), max.overlaps = Inf,
    segment.alpha = 0.45, segment.size = 0.2,
    box.padding = 0.12, point.padding = 0.07
  ) +
  theme_void(base_size = MIN_FONT_PT) +
  labs(
    title = "Strong GI network of Fig 2H genes",
    subtitle = paste0(
      ecount(g_strong), " strong GI edges (|GI Z| > ", GI_STRONG_Z_THRESHOLD,
      "); ", vcount(g_strong), " connected genes\n",
      "FR layout; node color = functional module; edge color/width = GI Z"
    )
  ) +
  theme(
    plot.title      = element_text(face = "plain", hjust = 0, size = FS_TITLE),
    plot.subtitle   = element_text(hjust = 0, color = "gray40", size = FS_SUBTITLE),
    legend.position = "right",
    legend.title    = element_text(size = FS_LEGEND, face = "plain"),
    legend.text     = element_text(size = FS_LEGEND),
    legend.key.size = unit(6, "pt"),
    plot.margin     = margin(3, 3, 3, 3)
  )

save_panel(p2l, unplaced_fig("supp_extra_fig2I_fr_layout.pdf"))

message(sprintf(
  "  Saved Fig2I strong-edge FR network | %d connected genes, %d edges with |GI Z| > %.1f",
  vcount(g_strong), ecount(g_strong), GI_STRONG_Z_THRESHOLD
))

message("--- SuppFig2Q: Functional Module Epistasis Pattern ---")

gene1_modules <- c(
  Pten   = "Growth_Control", Apc    = "Growth_Control",
  Rb1    = "Growth_Control", Smad4  = "Growth_Control",
  Stk11  = "Growth_Control", Tsc2   = "Growth_Control",
  Brca1  = "DNA_Repair",     Brca2  = "DNA_Repair",
  Chek2  = "DNA_Repair",
  Calml3 = "Barrier_ECM",    Serpinb5 = "Barrier_ECM",
  Sox10  = "Cell_Identity",  Krt15  = "Cell_Identity",
  Clca2  = "Cell_Identity",
  Tgfbr2 = "Unknown",        Kit    = "Unknown",
  Chit1  = "Unknown"
)

gene2_modules <- c(
  Btla     = "Immune_Evasion", Cxcr2    = "Immune_Evasion",
  Cxcr5    = "Immune_Evasion", Fcrl6    = "Immune_Evasion",
  Il6      = "Immune_Evasion", Mpeg1    = "Immune_Evasion",
  Pdcd1lg2 = "Immune_Evasion", Tlr7    = "Immune_Evasion",
  Tnfsf15  = "Immune_Evasion",
  Cx3cl1   = "Barrier_ECM",   Mfap5    = "Barrier_ECM",
  Muc16    = "Barrier_ECM",   Pla2g2a  = "Barrier_ECM",
  Pla2r1   = "Barrier_ECM",
  Cdh1     = "Cell_Identity",
  Hgf      = "Growth_Signaling",
  Mysm1    = "Unknown"
)

gi2_df <- read.csv(file.path(OUTPUT_DIR, "R2_GI_derived.csv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
gi2_df <- as.data.frame(gi2_df)
gi2_df <- gi2_df[!.is_nt(gi2_df$Gene1) & !.is_nt(gi2_df$Gene2) &
                   !is.na(gi2_df$Gene2), ]

gi2_df$Gene1_Pathway <- gene1_modules[gi2_df$Gene1]
gi2_df$Gene2_Pathway <- gene2_modules[gi2_df$Gene2]
gi2_df <- gi2_df[!is.na(gi2_df$Gene1_Pathway) & !is.na(gi2_df$Gene2_Pathway), ]
gi2_df$interaction_class <- ifelse(gi2_df$GI > 0, "synergy", "buffering")

gi2_df$Module_Interaction <- mapply(function(a, b) {
  pair <- sort(c(a, b))
  if (pair[1] == pair[2]) paste0(pair[1], "_Redundancy")
  else paste(pair[1], pair[2], sep = "_x_")
}, gi2_df$Gene1_Pathway, gi2_df$Gene2_Pathway)

counts_m <- as.data.frame(table(gi2_df$Module_Interaction,
                                 gi2_df$interaction_class))
names(counts_m) <- c("Module_Interaction", "interaction_class", "n")

wide_m <- reshape(counts_m,
                  idvar    = "Module_Interaction",
                  timevar  = "interaction_class",
                  direction = "wide")
names(wide_m) <- gsub("n\\.", "", names(wide_m))
for (col in c("buffering", "synergy")) {
  if (!col %in% names(wide_m)) wide_m[[col]] <- 0
}
wide_m$total        <- wide_m$buffering + wide_m$synergy
wide_m$synergy_pct  <- wide_m$synergy   / wide_m$total * 100
wide_m$buffering_pct<- wide_m$buffering / wide_m$total * 100
summary_m <- wide_m[order(wide_m$synergy_pct), ]

write.csv(summary_m, file.path(OUTPUT_DIR, "module_epistasis_summary.csv"),
          row.names = FALSE)

make_epistasis_barplot <- function(smry) {

  plot_long <- rbind(
    data.frame(Module_Interaction = smry$Module_Interaction,
               class = "Buffering (GI<0)", pct = smry$buffering_pct),
    data.frame(Module_Interaction = smry$Module_Interaction,
               class = "Synergy (GI>0)",   pct = smry$synergy_pct)
  )
  plot_long$Module_Interaction <- factor(plot_long$Module_Interaction,
                                          levels = smry$Module_Interaction)

  lbl <- smry
  lbl$Module_Interaction <- factor(lbl$Module_Interaction,
                                    levels = smry$Module_Interaction)

  lbl$syn_x  <- lbl$synergy_pct / 2
  lbl$buff_x <- lbl$synergy_pct + lbl$buffering_pct / 2

  ggplot(plot_long, aes(pct, Module_Interaction, fill = class)) +
    geom_col(alpha = 0.85, color = "black", linewidth = LW_THIN) +
    geom_text(data = subset(lbl, synergy_pct > 5),
              aes(x = syn_x, y = Module_Interaction,
                  label = paste0(round(synergy_pct), "%")),
              inherit.aes = FALSE, color = "white",
              fontface = "plain", size = pt(FS_LABEL)) +
    geom_text(data = subset(lbl, buffering_pct > 5),
              aes(x = buff_x, y = Module_Interaction,
                  label = paste0(round(buffering_pct), "%")),
              inherit.aes = FALSE, color = "white",
              fontface = "plain", size = pt(FS_LABEL)) +

    scale_fill_manual(values = c("Buffering (GI<0)" = "#3498DB",
                                  "Synergy (GI>0)"  = "#E74C3C")) +
    scale_x_continuous(limits = c(0, 100), expand = c(0, 0),
                       breaks = c(0, 25, 50, 75, 100)) +
    labs(
      title = "Epistasis by functional module",
      x = "Percentage (%)", y = NULL, fill = NULL
    ) +
    theme_panel(grid_y = FALSE) +
    theme(
      axis.text.y     = element_text(size = FS_AXIS_TEXT, face = "plain"),
      legend.position = "bottom",
      legend.key.size = unit(5, "pt")
    )
}

summary_nounknown <- summary_m[!grepl("Unknown", summary_m$Module_Interaction), ]
p2l <- make_epistasis_barplot(summary_nounknown)
save_panel(p2l, supp_fig("SuppFig2Q_module_epistasis.pdf"))

message(sprintf("  Saved SuppFig2Q_module_epistasis.pdf/.png -> _supplementary/ | %d combinations, %d pairs",
                nrow(summary_nounknown), sum(summary_nounknown$total)))

message("\n=== Main-figure R panels saved to: ", FIG_DIR, " ===\n")
message("This script (Step 3) produces main Fig2C (PCA), Fig2D (rank), Fig2G (volcano),")
message("  Fig2H (GI heatmap), Fig2I (GI network), plus SuppFig2Q (module epistasis).")
message("Step 1 (01_prepare_r1_gi_tables.R) produces Fig2E (guide LFC) + R1 tables;")
message("Step 2 (02_r1_modules.R) produces SuppFig2J (spectral) + SuppFig2K (forest).")
message("Full workflow: Rscript scripts/make_figure_2.R")

log_session("03_main_panels")
