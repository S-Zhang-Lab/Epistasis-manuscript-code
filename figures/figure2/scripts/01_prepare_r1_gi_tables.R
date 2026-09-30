suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

source(here::here("R", "paths.R"))
source(here::here("R", "helpers.R"))

source(here::here("R", "theme_pub.R"))

INPUT_CSV <- raw("R1_construct_counts.csv")
OUT       <- OUT_TABLES
FIG_DIR   <- OUT_FIGS

PSEUDO <- 0.5
N_CAP <- 25
DAY3_THR <- 3.0
MIN_GUIDES <- 5
FDR_THRESH <- 0.05
N_LABEL <- 10
SE_FLOOR <- 0.05

day18_cols <- paste0("Mice_Day18_", 1:6)
day3_cols <- paste0("Mice_Day3_", 1:6)

day0_reps <- c("Cell_Day0", "Cell_Day0_2", "Cell_Day0_3")
day0_gate <- "Cell_Day0_pooled"

DAY0_MIN <- 5
C_ENR <- "#d62728"
C_DEP <- "#1f77b4"
C_BG <- "#aaaaaa"
C_NT <- "#ff7f0e"

message("Loading R1 counts: ", INPUT_CSV)
df <- read.csv(INPUT_CSV, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("sgRNA", "gene", day18_cols, day3_cols, day0_reps) %in% colnames(df)))
df[[day0_gate]] <- rowSums(df[, day0_reps, drop = FALSE])

NTNT_LABEL <- "NTC"
gene_parts <- strsplit(df$gene, "\\*", fixed = FALSE)
df$gene_A <- vapply(gene_parts, `[`, character(1), 1)
df$gene_B <- vapply(gene_parts, function(x) if (length(x) >= 2) x[2] else x[1],
                    character(1))
df$gene_A[df$gene == NTNT_LABEL] <- "Nt"
df$gene_B[df$gene == NTNT_LABEL] <- "Nt"

message("Normalizing counts with median-ratio size factors")
counts <- as.matrix(df[, c(day18_cols, day3_cols)])
storage.mode(counts) <- "double"
counts <- counts + PSEUDO
geo <- exp(rowMeans(log(counts)))
ratio <- sweep(counts, 1, geo, "/")
size_factors <- apply(ratio, 2, median, na.rm = TRUE)
norm_counts <- sweep(counts, 2, size_factors, "/")

norm_d18 <- norm_counts[, day18_cols, drop = FALSE]
norm_d3 <- norm_counts[, day3_cols, drop = FALSE]
log2_d18 <- log2(norm_d18)
log2_d3 <- log2(norm_d3)

df$lfc_guide <- rowMeans(log2_d18) - rowMeans(log2_d3)
df$mean_norm_count <- (rowMeans(norm_d18) + rowMeans(norm_d3)) / 2
df$se_lfc <- sqrt(apply(log2_d18, 1, var) / 6 + apply(log2_d3, 1, var) / 6)
df$mean_d3 <- rowMeans(norm_d3)

df_filt <- df[df$mean_d3 >= DAY3_THR, , drop = FALSE]
message(sprintf("Guides passing Day3 >= %.1f: %d / %d", DAY3_THR, nrow(df_filt), nrow(df)))

message("Aggregating per gene pair")
aggregate_one_pair <- function(grp) {
  n_orig <- nrow(grp)
  if (nrow(grp) > N_CAP) {
    grp <- grp[order(grp$mean_norm_count, decreasing = TRUE), , drop = FALSE]
    grp <- grp[seq_len(N_CAP), , drop = FALSE]
  }
  ses <- pmax(grp$se_lfc, SE_FLOOR)
  w <- 1 / ses^2
  ws <- sum(w)
  wm <- sum(w * grp$lfc_guide) / ws
  se_c <- 1 / sqrt(ws)
  data.frame(
    gene_pair = grp$gene[1],
    gene_A = grp$gene_A[1],
    gene_B = grp$gene_B[1],
    n_guides_original = n_orig,
    n_guides_total = nrow(grp),
    weighted_mean_lfc = wm,
    se_combined = se_c,
    z_score = ifelse(se_c > 0, wm / se_c, 0),
    stringsAsFactors = FALSE
  )
}

aggregate_one_pair_re <- function(grp) {
  n_orig <- nrow(grp)
  if (nrow(grp) > N_CAP) {
    grp <- grp[order(grp$mean_norm_count, decreasing = TRUE), , drop = FALSE]
    grp <- grp[seq_len(N_CAP), , drop = FALSE]
  }
  ses <- pmax(grp$se_lfc, SE_FLOOR)
  w <- 1 / ses^2
  wm_fe <- sum(w * grp$lfc_guide) / sum(w)
  k <- nrow(grp)
  tau2 <- if (k > 1) {
    Q <- sum(w * (grp$lfc_guide - wm_fe)^2)
    max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  } else 0
  w <- 1 / (ses^2 + tau2)
  wm <- sum(w * grp$lfc_guide) / sum(w)
  se_c <- 1 / sqrt(sum(w))
  data.frame(
    gene_pair = grp$gene[1], gene_A = grp$gene_A[1], gene_B = grp$gene_B[1],
    n_guides_original = n_orig, n_guides_total = k,
    weighted_mean_lfc = wm, se_combined = se_c, tau2 = tau2,
    z_score = ifelse(se_c > 0, wm / se_c, 0),
    stringsAsFactors = FALSE
  )
}

results <- df_filt |>
  group_by(gene) |>
  group_modify(~ aggregate_one_pair(.x), .keep = TRUE) |>
  ungroup() |>
  select(-gene) |>
  filter(n_guides_total >= MIN_GUIDES)

nt_row <- results[toupper(results$gene_pair) %in% c("NT*NT", toupper(NTNT_LABEL)),
                  , drop = FALSE]
if (nrow(nt_row) == 0) stop("Nt x Nt reference not found after filtering")
NT_LFC <- nt_row$weighted_mean_lfc[1]
message(sprintf("NT*NT raw LFC: %.4f", NT_LFC))

results <- results |>
  mutate(
    weighted_mean_lfc_raw = weighted_mean_lfc,
    weighted_mean_lfc = weighted_mean_lfc - NT_LFC,
    z_score = weighted_mean_lfc / se_combined,
    pval = 2 * pnorm(-abs(z_score)),
    fdr = p.adjust(pval, method = "BH"),
    pair_type = case_when(
      toupper(gene_pair) == toupper(NTNT_LABEL)          ~ "nt_nt",
      toupper(gene_A) == "NT" & toupper(gene_B) == "NT" ~ "nt_nt",
      toupper(gene_A) == "NT" | toupper(gene_B) == "NT" ~ "single_KO",
      TRUE ~ "dual_KO"
    )
  )

write.csv(results, file.path(OUT, "all_gene_pairs_full.csv"),
          row.names = FALSE, quote = FALSE)

df_guides <- df |>
  select(sgRNA, gene, gene_A, gene_B, lfc_guide, se_lfc,
         mean_norm_count, mean_d3) |>
  mutate(passed_day3_filter = mean_d3 >= DAY3_THR)
write.csv(df_guides, file.path(OUT, "all_guides_lfc.csv"),
          row.names = FALSE, quote = FALSE)

message("Building single-KO references and GI scores")
single_ko <- results |> filter(pair_type == "single_KO")
single_long <- single_ko |>
  mutate(gene = ifelse(toupper(gene_A) == "NT", toupper(gene_B), toupper(gene_A)),
         se_use = pmax(se_combined, SE_FLOOR),
         weight = 1 / se_use^2)

singleKO_df <- single_long |>
  group_by(gene) |>
  summarise(
    meta_lfc = sum(weight * weighted_mean_lfc) / sum(weight),
    meta_se = 1 / sqrt(sum(weight)),
    n_guides_singleKO = sum(n_guides_total),
    .groups = "drop"
  ) |>
  arrange(gene)

write.csv(singleKO_df |> select(gene, meta_lfc, meta_se),
          file.path(OUT, "singleKO_reference_table.csv"),
          row.names = FALSE, quote = FALSE)

meta_lfc <- setNames(singleKO_df$meta_lfc, singleKO_df$gene)
meta_se <- setNames(singleKO_df$meta_se, singleKO_df$gene)
n_guides_single <- setNames(singleKO_df$n_guides_singleKO, singleKO_df$gene)

dual_ko <- results |> filter(pair_type == "dual_KO") |>
  mutate(gene_A_u = toupper(gene_A), gene_B_u = toupper(gene_B)) |>
  filter(gene_A_u %in% names(meta_lfc), gene_B_u %in% names(meta_lfc))

gi_df <- dual_ko |>
  transmute(
    gene_pair, gene_A, gene_B,
    lfc_AB = weighted_mean_lfc,
    se_AB = pmax(se_combined, SE_FLOOR),
    lfc_A_meta = meta_lfc[gene_A_u],
    se_A_meta = meta_se[gene_A_u],
    lfc_B_meta = meta_lfc[gene_B_u],
    se_B_meta = meta_se[gene_B_u],
    n_guides_AB = n_guides_total,
    n_guides_A = as.integer(n_guides_single[gene_A_u]),
    n_guides_B = as.integer(n_guides_single[gene_B_u])
  ) |>
  mutate(
    GI_score = lfc_AB - lfc_A_meta - lfc_B_meta,
    SE_GI = sqrt(se_AB^2 + se_A_meta^2 + se_B_meta^2),
    GI_Z = ifelse(SE_GI > 0, GI_score / SE_GI, 0),
    pval_GI = 2 * pnorm(-abs(GI_Z)),
    fdr_GI = p.adjust(pval_GI, method = "BH")
  ) |>
  select(gene_pair, gene_A, gene_B, lfc_AB, se_AB, lfc_A_meta, se_A_meta,
         lfc_B_meta, se_B_meta, GI_score, SE_GI, GI_Z, n_guides_AB,
         n_guides_A, n_guides_B, pval_GI, fdr_GI)

write.csv(gi_df, file.path(OUT, "GI_scores_all_dualKO.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(gi_df |> filter(fdr_GI < FDR_THRESH, GI_score > 0),
          file.path(OUT, "GI_scores_synergy_sig.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(gi_df |> filter(fdr_GI < FDR_THRESH, GI_score < 0),
          file.path(OUT, "GI_scores_buffering_sig.csv"),
          row.names = FALSE, quote = FALSE)

message("Computing Day-0-referenced enrichment")

enrichment_vs_day0 <- function(treat_cols, tag) {
  cnt <- as.matrix(df[, c(treat_cols, day0_reps)])
  storage.mode(cnt) <- "double"
  cnt <- cnt + PSEUDO

  gate <- df[[day0_gate]] >= DAY0_MIN
  geo <- exp(rowMeans(log(cnt[gate, , drop = FALSE])))
  sf <- apply(sweep(cnt[gate, , drop = FALSE], 1, geo, "/"), 2,
              median, na.rm = TRUE)
  nc <- sweep(cnt, 2, sf, "/")
  lt <- log2(nc[, treat_cols, drop = FALSE])
  l0 <- log2(nc[, day0_reps, drop = FALSE])
  n0 <- length(day0_reps)

  g <- data.frame(
    sgRNA = df$sgRNA, gene = df$gene,
    gene_A = df$gene_A, gene_B = df$gene_B,
    raw_day0 = df[[day0_gate]],
    lfc_guide = rowMeans(lt) - rowMeans(l0),
    mean_norm_count = (rowMeans(nc[, treat_cols, drop = FALSE]) +
                         rowMeans(nc[, day0_reps, drop = FALSE])) / 2,

    se_lfc = sqrt(apply(lt, 1, var) / length(treat_cols) +
                    apply(l0, 1, var) / n0),
    stringsAsFactors = FALSE
  )
  g <- g[g$raw_day0 >= DAY0_MIN, , drop = FALSE]
  message(sprintf("  [%s] guides with pooled Day-0 >= %d: %d / %d",
                  tag, DAY0_MIN, nrow(g), nrow(df)))

  res <- g |>
    group_by(gene) |>
    group_modify(~ aggregate_one_pair_re(.x), .keep = TRUE) |>
    ungroup() |> select(-gene) |>
    filter(n_guides_total >= MIN_GUIDES)

  ntr <- res[toupper(res$gene_pair) %in% c("NT*NT", toupper(NTNT_LABEL)), ]
  if (nrow(ntr) == 0) stop("Nt x Nt reference absent from the ", tag, " contrast")
  pairs <- res |>
    mutate(
      weighted_mean_lfc_raw = weighted_mean_lfc,
      weighted_mean_lfc = weighted_mean_lfc - ntr$weighted_mean_lfc[1],
      z_score = weighted_mean_lfc / se_combined,
      pval = 2 * pnorm(-abs(z_score)),
      fdr = p.adjust(pval, method = "BH"),
      pair_type = case_when(
        toupper(gene_pair) == toupper(NTNT_LABEL)          ~ "nt_nt",
        toupper(gene_A) == "NT" & toupper(gene_B) == "NT" ~ "nt_nt",
        toupper(gene_A) == "NT" | toupper(gene_B) == "NT" ~ "single_KO",
        TRUE                                              ~ "dual_KO"
      ),
      contrast = tag
    )

  g$lfc_guide <- g$lfc_guide - ntr$weighted_mean_lfc[1]
  list(pairs = pairs, guides = g)
}

invivo_cols  <- day18_cols
invitro_cols <- grep("^Cell_Day14_", colnames(df), value = TRUE)

enr_invivo <- enrichment_vs_day0(invivo_cols, "Day18_vs_Day0")
write.csv(enr_invivo$pairs, file.path(OUT, "R1_enrichment_day0_invivo.csv"),
          row.names = FALSE, quote = FALSE)
message(sprintf("  in vivo: %d pairs | enriched %d, depleted %d at FDR < %.2f",
                nrow(enr_invivo$pairs),
                sum(enr_invivo$pairs$fdr < FDR_THRESH & enr_invivo$pairs$weighted_mean_lfc > 0),
                sum(enr_invivo$pairs$fdr < FDR_THRESH & enr_invivo$pairs$weighted_mean_lfc < 0),
                FDR_THRESH))

if (length(invitro_cols) > 0) {
  enr_invitro <- enrichment_vs_day0(invitro_cols, "Day14cells_vs_Day0")
  write.csv(enr_invitro$pairs, file.path(OUT, "R1_enrichment_day0_invitro.csv"),
            row.names = FALSE, quote = FALSE)
  message(sprintf("  in vitro: %d pairs | enriched %d, depleted %d at FDR < %.2f",
                  nrow(enr_invitro$pairs),
                  sum(enr_invitro$pairs$fdr < FDR_THRESH & enr_invitro$pairs$weighted_mean_lfc > 0),
                  sum(enr_invitro$pairs$fdr < FDR_THRESH & enr_invitro$pairs$weighted_mean_lfc < 0),
                  FDR_THRESH))
} else {
  warning("No Cell_Day14_* columns found; skipping the in vitro control contrast")
}

message("Drawing Fig2E guide LFC distribution")

ranked <- enr_invivo$pairs |>
  filter(pair_type == "dual_KO", fdr < FDR_THRESH, weighted_mean_lfc > 0) |>
  arrange(desc(weighted_mean_lfc))
if (nrow(ranked) == 0) stop("No dual-knockout pair is enriched at FDR < ", FDR_THRESH)
top10_pairs <- ranked$gene_pair[seq_len(min(N_LABEL, nrow(ranked)))]
message(sprintf("  Fig2E: top %d of %d enriched dual-KO pairs (FDR < %.2f)",
                length(top10_pairs), nrow(ranked), FDR_THRESH))

plot_guides <- enr_invivo$guides |>
  filter(gene %in% top10_pairs) |>
  mutate(
    gene = factor(gene, levels = rev(top10_pairs)),
    color_group = if_else(lfc_guide > 0, "LFC > 0", "LFC <= 0")
  )

# Numerical source data for Figure 2E: this is the exact data frame supplied to
# ggplot below, after Day-0 gating, NT x NT centering, and top-pair selection.
write.csv(
  plot_guides |>
    select(sgRNA, gene, gene_A, gene_B, raw_day0, lfc_guide,
           mean_norm_count, se_lfc, color_group),
  file.path(OUT, "Fig2E_guide_lfc_source_data.csv"),
  row.names = FALSE,
  quote = FALSE
)

x_range <- range(plot_guides$lfc_guide, na.rm = TRUE)
p2d <- ggplot(plot_guides, aes(x = lfc_guide, y = gene, color = color_group,
                               alpha = color_group)) +
  geom_vline(xintercept = 0, linewidth = LW_RULE, color = "black") +
  geom_segment(aes(xend = lfc_guide,
                   y = as.numeric(gene) - 0.34,
                   yend = as.numeric(gene) + 0.34),
               linewidth = 0.55) +
  scale_y_continuous(
    breaks = seq_along(levels(plot_guides$gene)),
    labels = pair_label(levels(plot_guides$gene))
  ) +
  scale_color_manual(values = c("LFC > 0" = C_ENR, "LFC <= 0" = C_DEP),
                     name = NULL) +
  scale_alpha_manual(values = c("LFC > 0" = 1, "LFC <= 0" = 1), guide = "none") +
  coord_cartesian(xlim = x_range + c(-0.3, 0.3)) +
  labs(
    title = paste0("Per-guide LFC, top ", length(top10_pairs), " enriched pairs"),
    subtitle = paste0("Guides with pooled Day-0 >= ", DAY0_MIN,
                      " raw reads; pairs at FDR < ", FDR_THRESH),

    x = expression("Log"[2]*"(Fold Change)  (Day 18 / Day 0)"),
    y = NULL
  ) +
  theme_panel(grid_y = FALSE) +
  theme(
    axis.text.y     = element_text(size = FS_AXIS_TEXT, face = "plain"),
    legend.position = "bottom",
    legend.key.size = unit(5, "pt")
  )

save_panel(p2d, file.path(FIG_DIR, "Fig2E_guide_lfc.pdf"))

message(sprintf("Wrote R1 outputs: %d gene pairs, %d GI pairs",
                nrow(results), nrow(gi_df)))

log_session("01_prepare_r1_gi_tables")
