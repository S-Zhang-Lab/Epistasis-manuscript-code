suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "helpers.R"))
  source(here::here("R", "gi_pipeline.R"))
  source(here::here("R", "theme_pub.R"))
  source(here::here("R", "palettes.R"))
  for (pkg in c("ggplot2", "dplyr", "readr", "patchwork", "ggrepel")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
})

SEED    <- 123L
N_PERM  <- 1000L
set.seed(SEED)

message("\n=== proto_panelA_hostconditional ===\n")

gi <- readr::read_csv(tbl("GI_Cx3cr1.csv"), show_col_types = FALSE,
                      progress = FALSE) |>
  dplyr::mutate(delta_GI = GI_Depleted_vs_Cell - GI_WT_vs_Cell)

norm <- load_and_normalize(raw("Cx3cr1_raw_counts.csv"))
log2cpm     <- norm$log2cpm
sample_cols <- norm$sample_cols
gene_map    <- norm$gene_map
groups      <- classify_samples(sample_cols)

cell_cols <- groups$cell
host_cols <- c(groups$wt, groups$depleted)
n_wt      <- length(groups$wt)
stopifnot(length(host_cols) == n_wt + length(groups$depleted))

mean_cell <- rowMeans(log2cpm[, cell_cols, drop = FALSE])

guide_gene <- vapply(unname(gene_map[rownames(log2cpm)]),
                     convert_gene_pair, character(1L), USE.NAMES = FALSE)
ntnt_mask  <- guide_gene == "NT*NT"

gene_groups <- split(which(!ntnt_mask), guide_gene[!ntnt_mask])
gene_names  <- names(gene_groups)
gene_type   <- vapply(gene_names, pair_type, character(1L))
top_n_vec   <- ifelse(gene_type == "dual", TOP_DUAL, TOP_NT)
names(top_n_vec) <- gene_names

dual_tok <- gi$gene
sA_tok   <- gi$singleA_pair
sB_tok   <- gi$singleB_pair
stopifnot(all(c(dual_tok, sA_tok, sB_tok) %in% gene_names))

delta_gi_for_split <- function(wt_idx, dep_idx) {
  lfc_wt  <- rowMeans(log2cpm[, wt_idx,  drop = FALSE]) - mean_cell
  lfc_dep <- rowMeans(log2cpm[, dep_idx, drop = FALSE]) - mean_cell
  ntnt_wt  <- mean(lfc_wt[ntnt_mask])
  ntnt_dep <- mean(lfc_dep[ntnt_mask])

  corr_wt <- vapply(gene_names, function(g) {
    idx <- gene_groups[[g]]; n <- min(top_n_vec[[g]], length(idx))
    mean(sort(lfc_wt[idx], decreasing = TRUE)[seq_len(n)]) - ntnt_wt
  }, numeric(1L))
  corr_dep <- vapply(gene_names, function(g) {
    idx <- gene_groups[[g]]; n <- min(top_n_vec[[g]], length(idx))
    mean(sort(lfc_dep[idx], decreasing = TRUE)[seq_len(n)]) - ntnt_dep
  }, numeric(1L))

  gi_wt  <- corr_wt[dual_tok]  - corr_wt[sA_tok]  - corr_wt[sB_tok]
  gi_dep <- corr_dep[dual_tok] - corr_dep[sA_tok] - corr_dep[sB_tok]
  unname(gi_dep - gi_wt)
}

dgi_recomputed <- delta_gi_for_split(groups$wt, groups$depleted)
max_dev <- max(abs(dgi_recomputed - gi$delta_GI))
message(sprintf("[validate] max |recomputed - CSV delta_GI| = %.4g", max_dev))
if (max_dev > 1e-3)
  stop("Fast permutation path diverges from canonical GI math; aborting.")

message(sprintf("[perm] running %d host-label permutations ...", N_PERM))
t0 <- Sys.time()
null_mat <- matrix(NA_real_, nrow = nrow(gi), ncol = N_PERM)
for (b in seq_len(N_PERM)) {
  perm <- sample(host_cols)
  null_mat[, b] <- delta_gi_for_split(perm[seq_len(n_wt)],
                                      perm[(n_wt + 1L):length(host_cols)])
}
message(sprintf("[perm] done in %.1f s", as.numeric(difftime(Sys.time(), t0,
                                                             units = "secs"))))

abs_obs   <- abs(gi$delta_GI)
perm_p    <- (rowSums(abs(null_mat) >= abs_obs) + 1) / (N_PERM + 1)
gi$delta_GI_perm_p <- perm_p
gi$delta_GI_perm_q <- p.adjust(perm_p, method = "BH")
null_band <- as.numeric(quantile(abs(as.vector(null_mat)), 0.95))

qW <- gi$GI_WT_qvalue; qD <- gi$GI_Depleted_qvalue; qd <- gi$delta_GI_qvalue
gi <- gi |>
  dplyr::mutate(
    is_interaction = (qW < 0.05) | (qD < 0.05),
    is_hostcond    = is_interaction & (qd < 0.05),
    is_signflip    = (qW < 0.05) & (qD < 0.05) &
                     (sign(GI_WT_vs_Cell) != sign(GI_Depleted_vs_Cell)),
    rewire_class = dplyr::case_when(
      is_signflip                ~ "GI flip",
      is_hostcond                ~ "Niche-reprogrammed",
      is_interaction             ~ "Niche-independent",
      TRUE                       ~ "No interaction")
  )

n_int   <- sum(gi$is_interaction)
n_hc    <- sum(gi$is_hostcond)
n_flip  <- sum(gi$is_signflip)
pct_hc  <- 100 * n_hc / n_int
pct_all <- 100 * mean(qd < 0.05)

n_perm_sig   <- sum(gi$delta_GI_perm_q < 0.05)
n_concord    <- sum(gi$is_hostcond & gi$delta_GI_perm_q < 0.05)

null_display  <- 5

cat(sprintf(paste0(
  "\n---------------- PROTOTYPE NUMBERS ----------------\n",
  "  dual pairs total              : %d\n",
  "  genuine interactions (q<.05)  : %d\n",
  "  host-conditional (bootstrap)  : %d  (%.1f%% of interactions)\n",
  "  host-conditional (perm FDR)   : %d  (concordant with bootstrap: %d)\n",
  "  strict sign-flips             : %d\n",
  "  host-conditional, all pairs   : %.1f%%\n",
  "  global 95%% null band (|dGI|)  : %.3f\n",
  "---------------------------------------------------\n\n"),
  nrow(gi), n_int, n_hc, pct_hc, n_perm_sig, n_concord, n_flip, pct_all,
  null_band))

gi <- gi |>
  dplyr::mutate(
    wpos = GI_WT_vs_Cell >= 0, dpos = GI_Depleted_vs_Cell >= 0,
    category = dplyr::case_when(
      !wpos & !dpos ~ "Conserved Buffering",
       wpos &  dpos ~ "Conserved Synergistic",
      !wpos &  dpos ~ "Niche-masked GI",
       wpos & !dpos ~ "Niche-enabled GI"))

CAT_COLORS <- c("Conserved Buffering"   = "#4DAF4A",
                "Conserved Synergistic" = "#FF7F00",
                "Niche-masked GI"       = "#E41A1C",
                "Niche-enabled GI"      = "#377EB8")

ZL <- 2.9
band_df <- data.frame(x = c(-ZL, ZL))

off <- gi |> dplyr::filter(abs(GI_WT_vs_Cell) > ZL |
                           abs(GI_Depleted_vs_Cell) > ZL)

hero <- gi |> dplyr::filter(gene == "Pten*Cx3cl1")

build_scatter <- function(s = 1, n_lab = 10) {
  top_rep <- gi |> dplyr::filter(is_hostcond) |>
    dplyr::slice_max(abs(delta_GI), n = n_lab)
  lab_df <- gi |>
    dplyr::filter(is_signflip | gene %in% top_rep$gene) |>
    dplyr::filter(abs(GI_WT_vs_Cell) <= ZL, abs(GI_Depleted_vs_Cell) <= ZL)

  p_scatter <- ggplot(gi, aes(GI_WT_vs_Cell, GI_Depleted_vs_Cell)) +

  geom_ribbon(data = band_df,
              aes(x = x, ymin = x - null_band, ymax = x + null_band),
              inherit.aes = FALSE, fill = "grey80", alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey65",
             linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey65",
             linewidth = 0.3) +
  geom_point(aes(fill = category, size = abs(delta_GI), alpha = is_hostcond),
             shape = 21, color = "white", stroke = 0.25) +
  ggrepel::geom_text_repel(data = lab_df,
             aes(label = gsub("\\*", "_", gene)),
             size = 2.0 * s, color = "grey15", max.overlaps = Inf,
             segment.size = 0.2, min.segment.length = 0,
             box.padding = 0.3, seed = SEED) +
  scale_fill_manual(values = CAT_COLORS, name = NULL,
                    breaks = names(CAT_COLORS)) +
  scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.22),
                     guide = "none") +
  scale_size_area(max_size = 3.4 * s, name = expression("|" * Delta * "GI|")) +
  coord_equal(xlim = c(-ZL, ZL), ylim = c(-ZL, ZL), clip = "on") +
  labs(
    title = "The niche reprograms epistasis",
    subtitle = sprintf(
      "%.0f%% of interactions host-conditional vs ~5%% by chance (FDR<0.05); opaque = rewired, faint = conserved",
      pct_hc),

    x = expression("GI score    WT host  (niche intact)"),
    y = expression("GI score    CX3CR1-depleted host")) +
  theme_pub(base_size = 7 * s) +
  theme(legend.position = "right",
        plot.title    = element_text(size = 9 * s, hjust = 0),
        plot.subtitle = element_text(size = 5.5 * s, colour = "grey25"),
        axis.title    = element_text(size = 8 * s),
        axis.text     = element_text(size = 7 * s),
        legend.text   = element_text(size = 5.5 * s),
        legend.title  = element_text(size = 6 * s),
        legend.key.height = unit(7 * s, "pt"))

if (nrow(hero) == 1) {
  p_scatter <- p_scatter +
    annotate("segment",
             x = hero$GI_WT_vs_Cell, xend = hero$GI_WT_vs_Cell,
             y = hero$GI_WT_vs_Cell, yend = hero$GI_Depleted_vs_Cell,
             arrow = arrow(length = unit(3, "pt"), type = "closed"),
             color = "#B2182B", linewidth = 0.5)
}

catlab <- function(x, y, hjust, vjust, txt, col)
  annotate("text", x = x, y = y, hjust = hjust, vjust = vjust, size = 1.8 * s,
           lineheight = 0.85, fontface = "bold", colour = col, label = txt)
p_scatter <- p_scatter +
  catlab(-ZL * 0.97,  ZL * 0.97, 0, 1, "Niche-masked\nGI", "#E41A1C") +
  catlab(-ZL * 0.97, -ZL * 0.97, 0, 0, "Conserved\nBuffering",     "#4DAF4A") +
  catlab( ZL * 0.97, -ZL * 0.97, 1, 0, "Niche-enabled\nGI",    "#377EB8")

if (nrow(off) > 0) {
  p_scatter <- p_scatter +
    annotate("point", x = ZL * 0.985, y = ZL * 0.985, shape = 24,
             fill = "#FF7F00", colour = "white", size = 2.1 * s) +
    annotate("text", x = ZL * 0.97, y = ZL * 0.72, hjust = 1, vjust = 1,
             size = 1.75 * s, lineheight = 0.85, colour = "#B25900",
             label = sprintf("Conserved Synergistic\n%s off-scale\n(%.1f, %.1f) Mrc1-gated (Fig 5)",
                             gsub("\\*", "_", off$gene[1]),
                             off$GI_WT_vs_Cell[1], off$GI_Depleted_vs_Cell[1]))
} else {
  p_scatter <- p_scatter +
    catlab(ZL * 0.97, ZL * 0.97, 1, 1, "Conserved\nSynergistic", "#FF7F00")
}

  p_scatter
}

RC_COLORS <- c("Niche-independent"          = "grey62",
               "Niche-reprogrammed" = "#377EB8",
               "GI flip"          = "#E41A1C")
stk <- dplyr::bind_rows(
  data.frame(bar = "Observed",
             class = c("Niche-independent", "Niche-reprogrammed", "GI flip"),
             pct = 100 * c(sum(gi$rewire_class == "Niche-independent"),
                           sum(gi$rewire_class == "Niche-reprogrammed"),
                           sum(gi$rewire_class == "GI flip")) / n_int),
  data.frame(bar = "Chance\n(null)",
             class = c("Niche-independent", "Niche-reprogrammed", "GI flip"),
             pct = c(100 - null_display, null_display, 0))
)
stk$class <- factor(stk$class,
                    levels = c("Niche-independent", "Niche-reprogrammed", "GI flip"))
stk$bar <- factor(stk$bar, levels = c("Observed", "Chance\n(null)"))

p_stack <- ggplot(stk, aes(bar, pct, fill = class)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(data = subset(stk, bar == "Observed" & class != "Niche-independent"),
            aes(label = sprintf("%.0f%%", pct)),
            position = position_stack(vjust = 0.5),
            size = 2.1, color = "white", fontface = "bold") +
  scale_fill_manual(values = RC_COLORS, name = NULL) +
  labs(title = "Pervasiveness", x = NULL, y = "% of interactions") +
  theme_pub(base_size = 7) +
  theme(legend.position = "bottom",
        legend.key.size = unit(7, "pt"),
        plot.title = element_text(size = 8))

save_one <- function(p, stem, w_mm, h_mm) {
  pdf <- fig(paste0(stem, ".pdf"))
  save_pdf(p, pdf, w_mm, h_mm)
  message("  wrote ", pdf)
}

save_one(build_scatter(s = 1, n_lab = 10), "Fig4D_GxE_scatter",   121, 110)
save_one(p_stack,                          "Fig4C_pervasiveness",  57,  85)

local({
  f <- fig("Fig4D_GxE_scatter_large.pdf")
  grDevices::pdf(f, width = 280 / 25.4, height = 250 / 25.4,
                 family = PUB_FONT_FAMILY, bg = "white", useDingbats = FALSE)
  on.exit(grDevices::dev.off())
  print(build_scatter(s = 1.9, n_lab = 16))
  message("  wrote ", f, " (large review render)")
})

message("\n[proto] complete.\n")

if (!interactive()) tryCatch(log_session("proto_panelA_hostconditional"), error = function(e) NULL)
