pkgs <- c("readr", "dplyr", "ggplot2", "scales", "patchwork")
inst <- rownames(installed.packages())
if (length(setdiff(pkgs, inst)))
  install.packages(setdiff(pkgs, inst), repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
  library(scales); library(patchwork)
})

source(here::here("R", "paths.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "theme_pub.R"))
source(here::here("R", "qc_library.R"))

accent <- "#2C7FB8"
FILTER <- 3L
REPS   <- c("Cell_Day0", "Cell_Day0_2", "Cell_Day0_3")

raw_df <- read.csv(raw("R1_construct_counts.csv"), check.names = FALSE,
                   stringsAsFactors = FALSE)
stopifnot(all(REPS %in% colnames(raw_df)))

rep_mat <- as.matrix(raw_df[, REPS, drop = FALSE])
storage.mode(rep_mat) <- "double"
depths  <- colSums(rep_mat)

keep <- apply(rep_mat, 1, max) >= FILTER
con <- data.frame(
  construct = raw_df$sgRNA[keep],
  reads     = rowSums(rep_mat[keep, , drop = FALSE]),
  stringsAsFactors = FALSE
)
con$pair <- qc_gene_pair(con$construct)

m <- qc_metrics(con$reads, con$pair)

cpm <- sweep(rep_mat[keep, , drop = FALSE], 2, depths, "/") * 1e6
lg  <- log10(cpm + 1)
pairs_idx <- combn(seq_along(REPS), 2, simplify = FALSE)
r_vals <- vapply(pairs_idx, function(ij) cor(lg[, ij[1]], lg[, ij[2]]), numeric(1))
names(r_vals) <- vapply(pairs_idx,
                        function(ij) paste(REPS[ij[1]], REPS[ij[2]], sep = " vs "),
                        character(1))
message(sprintf(
  "[SuppFig2I] cell: %s constructs (%.1f%% of %s), %s dual-KO pairs, median %s, Gini %.3f, mean r %.3f",
  comma(m$n_constructs), m$pct_detected, comma(R1_WHITELIST),
  comma(m$n_dual_pairs), comma(m$median_reads), m$gini, mean(r_vals)))

base_t <- theme_panel()

pA <- ggplot(con, aes(log10(reads))) +
  geom_histogram(aes(y = after_stat(density)), bins = 45,
                 fill = accent, colour = "white", alpha = 0.9) +
  geom_density(colour = "#3F007D", linewidth = 0.45) +
  geom_vline(xintercept = log10(m$median_reads), linetype = "dashed",
             colour = "grey30") +
  labs(x = expression("log"[10]*" (summed reads per construct)"), y = "Density",
       title = "Read-count distribution") +
  base_t

lz <- qc_lorenz(con$reads)
pB <- ggplot(lz, aes(p_constructs, p_reads)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(colour = "#1B7837", linewidth = 0.55) +
  annotate("text", x = 0.04, y = 0.93, hjust = 0, size = pt(FS_LABEL),
           label = sprintf("Gini = %.2f", m$gini)) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Cumulative % of constructs", y = "Cumulative % of reads",
       title = "Library evenness") +
  base_t

rep_df <- data.frame(x = lg[, 1], y = lg[, 2])
pC <- ggplot(rep_df, aes(x, y)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(colour = accent, alpha = 0.18, size = 0.25, stroke = 0) +
  annotate("text", x = min(rep_df$x), y = max(rep_df$y), hjust = 0,
           vjust = 1, size = pt(FS_LABEL),
           label = sprintf("r = %.2f", r_vals[1])) +
  labs(x = expression("Replicate 1, log"[10]*"(CPM + 1)"),
       y = expression("Replicate 2, log"[10]*"(CPM + 1)"),
       title = "Replicate reproducibility") +
  base_t

tbl_df <- tibble::tribble(
  ~Metric,                      ~Value,
  "Whitelist (designed)",       comma(R1_WHITELIST),
  "Constructs (min 3)",         comma(m$n_constructs),
  "Detected (% of whitelist)",  sprintf("%.1f%%", m$pct_detected),
  "Dual-KO pairs",              comma(m$n_dual_pairs),
  "Single-KO controls",         comma(m$n_single_ctrl),
  "Replicate depths (M)",       paste(sprintf("%.1f", depths / 1e6), collapse = " / "),
  "Median reads/construct",     comma(m$median_reads),
  "Gini coefficient",           sprintf("%.3f", m$gini),
  "Reads in top 10%",           sprintf("%.1f%%", m$top10_share),
  "Replicate r (log CPM)",      sprintf("%.2f-%.2f", min(r_vals), max(r_vals))
) |> mutate(y = rev(seq_len(dplyr::n())))

pD <- ggplot(tbl_df, aes(y = y)) +
  geom_text(aes(x = 0, label = Metric), hjust = 0, size = pt(FS_LABEL)) +
  geom_text(aes(x = 1, label = Value),  hjust = 1, size = pt(FS_LABEL),
            fontface = "plain") +
  scale_x_continuous(limits = c(-0.04, 1.04)) +
  labs(title = "QC metrics summary", x = NULL, y = NULL) +
  theme_void(base_size = MIN_FONT_PT) +
  theme(plot.title = element_text(size = FS_TITLE, hjust = 0))

for (nm in c("a", "b", "c", "d")) {
  p <- switch(nm, a = pA, b = pB, c = pC, d = pD)
  save_panel(p, candidate_fig(sprintf("SuppFig2I_cell_QC_%s.pdf", nm)))
}
save_panel_grid(list(pA, pB, pC, pD),
                candidate_fig("SuppFig2I_cell_QC_grid.pdf"), ncol = 2)

write.csv(
  data.frame(Metric = tbl_df$Metric, Value = tbl_df$Value),
  tbl("SuppFig2I_r1_cell_QC_metrics.csv"), row.names = FALSE)

log_session("supp_10_r1_cell_qc")
message("[SuppFig2I] wrote 4 candidate sub-panels + grid -> _supplementary/_candidates/")
