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

design <- qc_design_genes()

all_con <- read_csv(raw("R1_plasmid_counts.csv"), show_col_types = FALSE) |>
  filter(!is.na(freq), freq >= FILTER) |>
  mutate(class = qc_classify(crRNA, design),
         pair  = qc_gene_pair(crRNA))

brk <- all_con |> count(class, wt = freq, name = "reads") |>
  left_join(count(all_con, class, name = "constructs"), by = "class")
for (i in seq_len(nrow(brk))) {
  message(sprintf("  [SuppFig2H] %-15s %7s constructs  %13s reads (%5.2f%%)",
                  brk$class[i], comma(brk$constructs[i]), comma(brk$reads[i]),
                  100 * brk$reads[i] / sum(brk$reads)))
}

con <- filter(all_con, class == "canonical")
m <- qc_metrics(con$freq, con$pair)
message(sprintf(
  "[SuppFig2H] plasmid (cross-program only): %s constructs (%.1f%% of %s), %s dual-KO pairs + %s single-KO controls, median %s, Gini %.3f",
  comma(m$n_constructs), m$pct_detected, comma(R1_WHITELIST),
  comma(m$n_dual_pairs), comma(m$n_single_ctrl), comma(m$median_reads), m$gini))

R1_DESIGNED_PAIRS <- 7098L
stopifnot(m$n_dual_pairs <= R1_DESIGNED_PAIRS, m$n_single_ctrl == 169L)
message(sprintf("[SuppFig2H] design coverage: %s of %s dual-KO pairs (%.2f%%)",
                comma(m$n_dual_pairs), comma(R1_DESIGNED_PAIRS),
                100 * m$n_dual_pairs / R1_DESIGNED_PAIRS))

base_t <- theme_panel()

pA <- ggplot(con, aes(log10(freq))) +
  geom_histogram(aes(y = after_stat(density)), bins = 45,
                 fill = accent, colour = "white", alpha = 0.9) +
  geom_density(colour = "#3F007D", linewidth = 0.45) +
  geom_vline(xintercept = log10(m$median_reads), linetype = "dashed",
             colour = "grey30") +
  labs(x = expression("log"[10]*" (reads per construct)"), y = "Density",
       title = "Read-count distribution") +
  base_t

lz <- qc_lorenz(con$freq)
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

rk <- con |> arrange(desc(freq)) |> mutate(rank = row_number())
pC <- ggplot(rk, aes(rank, freq)) +
  geom_line(colour = accent, linewidth = 0.45) +
  geom_hline(yintercept = m$median_reads, linetype = "dashed", colour = "grey30") +
  scale_y_log10(labels = comma) +
  scale_x_continuous(labels = comma) +
  labs(x = "Constructs ranked by abundance", y = "Reads (log10)",
       title = "Ranked abundance") +
  base_t

tbl_df <- tibble::tribble(
  ~Metric,                      ~Value,
  "Whitelist (designed)",       comma(R1_WHITELIST),
  "Constructs (freq >= 3)",     comma(m$n_constructs),
  "Detected (% of whitelist)",  sprintf("%.1f%%", m$pct_detected),
  "Dual-KO pairs detected",     sprintf("%s / %s", comma(m$n_dual_pairs),
                                        comma(R1_DESIGNED_PAIRS)),
  "Single-KO controls",         comma(m$n_single_ctrl),
  "Total reads",                comma(m$total_reads),
  "Median reads/construct",     comma(m$median_reads),
  "Gini coefficient",           sprintf("%.3f", m$gini),
  "Skew ratio (90th/10th)",     sprintf("%.0fx", m$skew_ratio),
  "Reads in top 10%",           sprintf("%.1f%%", m$top10_share)
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
  save_panel(p, candidate_fig(sprintf("SuppFig2H_plasmid_QC_%s.pdf", nm)))
}
save_panel_grid(list(pA, pB, pC, pD),
                candidate_fig("SuppFig2H_plasmid_QC_grid.pdf"), ncol = 2)

write.csv(
  data.frame(Metric = tbl_df$Metric, Value = tbl_df$Value),
  tbl("SuppFig2H_r1_plasmid_QC_metrics.csv"), row.names = FALSE)

log_session("supp_09_r1_plasmid_qc")
message("[SuppFig2H] wrote 4 candidate sub-panels + grid -> _supplementary/_candidates/")
