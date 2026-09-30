# =====================================================================
# Day 0 cell library QC  (in vivo dual-guide library)
# Unit = sgRNA-pair construct (1248); baseline = Cell_1/2/3
# Panels:  A read-count distribution   B Lorenz curve + Gini
#          C replicate reproducibility D QC metrics summary
# =====================================================================

if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
need <- c("readr","dplyr","tidyr","ggplot2","scales","patchwork")
invisible(lapply(need, library, character.only = TRUE))

REPS <- c("Cell_1","Cell_2","Cell_3")
accent <- "#2C7FB8"

raw <- read_csv(derived("screens","2nd_DC_counts.csv"), show_col_types = FALSE)

# ---- per-construct: raw reads (sum reps) + CPM per replicate ----
cpm <- sweep(as.matrix(raw[REPS]), 2, colSums(raw[REPS]), "/") * 1e6
con <- tibble(sgRNA = raw$sgRNA,
              reads = rowSums(as.matrix(raw[REPS])),
              mean_cpm = rowMeans(cpm))
cpm_df <- as_tibble(cpm)

# ---- metrics ----
gini <- function(x){ x <- sort(x); n <- length(x)
  sum((2*seq_len(n) - n - 1) * x) / (n * sum(x)) }
N        <- nrow(con)
n_zero   <- sum(con$reads == 0)
pct_det  <- 100 * mean(con$reads > 0)
G        <- gini(con$reads[con$reads > 0])
q        <- quantile(con$reads[con$reads>0], c(.1,.5,.9))
skew9010 <- q[3] / q[1]
top10sh  <- 100 * sum(sort(con$reads, decreasing=TRUE)[1:ceiling(0.1*N)]) / sum(con$reads)
# replicate correlations (Pearson on log10 cpm; Spearman on raw cpm)
pr <- function(a,b) cor(log10(cpm_df[[a]]+1), log10(cpm_df[[b]]+1))
rP <- c(pr("Cell_1","Cell_2"), pr("Cell_1","Cell_3"), pr("Cell_2","Cell_3"))

base_t <- theme_bw(base_size = 11) +
  theme(plot.title = element_text(face="bold", size=11),
        panel.grid.minor = element_blank())

# ---- A: read-count distribution ----
det <- con |> filter(reads > 0) |> mutate(l10 = log10(reads))
pA <- ggplot(det, aes(l10)) +
  geom_histogram(aes(y = after_stat(density)), bins = 45,
                 fill = accent, colour = "white", alpha = 0.9) +
  geom_density(colour = "#08519C", linewidth = 0.9) +
  geom_vline(xintercept = log10(median(det$reads)), linetype = "dashed", colour="grey30") +
  labs(x = expression(log[10]~(Day~0~reads~per~construct)), y = "Density",
       title = "A  Read-count distribution") + base_t

# ---- B: Lorenz curve + Gini ----
lor <- con |> filter(reads>0) |> arrange(reads) |>
  mutate(cp = row_number()/n(), cr = cumsum(reads)/sum(reads))
pB <- ggplot(lor, aes(cp, cr)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey60") +
  geom_line(colour = "#1B7837", linewidth = 1.1) +
  annotate("text", x=.03, y=.9, hjust=0, fontface="bold", size=4,
           label = sprintf("Gini = %.2f", G)) +
  scale_x_continuous(labels=percent_format()) + scale_y_continuous(labels=percent_format()) +
  labs(x="Cumulative % of constructs", y="Cumulative % of reads",
       title="B  Library evenness (Lorenz curve)") + base_t

# ---- C: replicate reproducibility ----
pC <- ggplot(cpm_df, aes(Cell_2, Cell_3)) +
  geom_abline(slope=1, intercept=0, colour="grey70") +
  geom_point(alpha=0.25, size=0.9, colour=accent) +
  scale_x_log10(labels=comma) + scale_y_log10(labels=comma) +
  annotate("text", x=Inf, y=-Inf, hjust=1.05, vjust=-0.6, size=3.4,
           label = sprintf("Pearson r (log10 CPM):\nrep1-2 %.3f | rep1-3 %.3f | rep2-3 %.3f",
                           rP[1], rP[2], rP[3])) +
  labs(x="Cell_2 (CPM)", y="Cell_3 (CPM)",
       title="C  Replicate reproducibility") + base_t

# ---- D: QC metrics summary ----
metrics <- tribble(~Metric, ~Value,
  "Constructs (sgRNA pairs)", comma(N),
  "Detected (>=1 read)",      sprintf("%d (%.1f%%)", N - n_zero, pct_det),
  "Zero-count (dropout)",     as.character(n_zero),
  "Median reads/construct",   comma(round(q[2])),
  "Gini coefficient",         sprintf("%.3f", G),
  "Skew ratio (90th/10th)",   sprintf("%.1fx", skew9010),
  "Reads in top 10%",         sprintf("%.1f%%", top10sh),
  "Mean replicate r (log)",   sprintf("%.3f", mean(rP)))
write_csv(metrics, tbl("day0_library_QC_metrics.csv"))

pD <- ggplot(metrics, aes(x=0, y=factor(Metric, levels=rev(Metric)))) +
  geom_text(aes(label=Metric), hjust=0, size=3.6) +
  geom_text(aes(x=1, label=Value), hjust=1, size=3.6, fontface="bold") +
  scale_x_continuous(limits=c(0,1), expand=expansion(mult=c(0.02,0.02))) +
  labs(title="D  QC metrics summary", x=NULL, y=NULL) +
  theme_void(base_size=11) +
  theme(plot.title=element_text(face="bold", size=11, hjust=0),
        plot.margin=margin(10,10,10,10))

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(title = "Day 0 cell library QC (in vivo dual-guide screen)",
                  theme = theme(plot.title = element_text(face="bold", size=14)))
ggsave(fig_supp("fig_day0_library_QC.pdf"), fig, width = 10, height = 8)
ggsave(fig_supp("fig_day0_library_QC.png"), fig, width = 10, height = 8, dpi = 300)
message("Saved fig_day0_library_QC.pdf / .png and day0_library_QC_metrics.csv")
print(metrics)
