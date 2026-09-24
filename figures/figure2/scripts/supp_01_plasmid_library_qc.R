pkgs <- c("readr", "dplyr", "ggplot2", "scales", "patchwork", "stringr")
inst <- rownames(installed.packages())
if (length(setdiff(pkgs, inst)))
  install.packages(setdiff(pkgs, inst), repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
  library(scales); library(patchwork); library(stringr)
})

source(here::here("R", "paths.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "theme_pub.R"))
dir.create(OUT_SUPP, recursive = TRUE, showWarnings = FALSE)

accent <- "#6A51A3"

con <- read_delim(derived("R2_plasmid_counts.txt"), delim = ":", trim_ws = TRUE,
                  col_names = c("construct","reads"), show_col_types = FALSE) |>
  filter(construct != "Name1*Name2") |>
  mutate(reads = as.numeric(reads))

gini <- function(x){ x <- sort(x); n <- length(x)
  sum((2*seq_len(n) - n - 1) * x) / (n * sum(x)) }
N        <- nrow(con)
n_zero   <- sum(con$reads == 0)
pct_det  <- 100 * mean(con$reads > 0)
G        <- gini(con$reads[con$reads > 0])
q        <- quantile(con$reads[con$reads>0], c(.01,.1,.5,.9,.99))
skew9010 <- q["90%"] / q["10%"]
top10sh  <- 100 * sum(sort(con$reads, decreasing=TRUE)[1:ceiling(0.1*N)]) / sum(con$reads)
dyn      <- q["99%"] / q["1%"]

base_t <- theme_panel()

det <- con |> filter(reads>0) |> mutate(l10 = log10(reads))
pA <- ggplot(det, aes(l10)) +
  geom_histogram(aes(y=after_stat(density)), bins=45, fill=accent, colour="white", alpha=0.9) +
  geom_density(colour="#3F007D", linewidth=0.45) +
  geom_vline(xintercept=log10(median(det$reads)), linetype="dashed", colour="grey30") +
  labs(x=expression("log"[10]*" (reads per construct)"), y="Density",
       title="Read-count distribution") + base_t

lor <- con |> filter(reads>0) |> arrange(reads) |>
  mutate(cp=row_number()/n(), cr=cumsum(reads)/sum(reads))
pB <- ggplot(lor, aes(cp, cr)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey60") +
  geom_line(colour="#1B7837", linewidth=0.55) +
  annotate("text", x=.03, y=.9, hjust=0, fontface = "plain", size=pt(FS_LABEL),
           label=sprintf("Gini = %.2f", G)) +
  scale_x_continuous(labels=percent_format()) + scale_y_continuous(labels=percent_format()) +
  labs(x="Cumulative % of constructs", y="Cumulative % of reads",
       title="Library evenness") + base_t

ranked <- con |> arrange(desc(reads)) |> mutate(rank=row_number())
pC <- ggplot(ranked, aes(rank, reads)) +
  geom_line(colour=accent, linewidth=0.45) +
  geom_hline(yintercept=median(con$reads), linetype="dashed", colour="grey30") +
  annotate("text", x=N, y=median(con$reads), hjust=1, vjust=-0.6, size=pt(FS_LABEL), colour="grey30",
           label=sprintf("median = %d", round(median(con$reads)))) +
  scale_y_log10(labels=comma) +
  labs(x="Constructs ranked by abundance", y="Reads (log10)",
       title="Ranked abundance (dynamic range)") + base_t

metrics <- tibble::tribble(~Metric, ~Value,
  "Constructs",              comma(N),
  "Gene pairs",              "324",
  "Detected (>=1 read)",     sprintf("%d (%.1f%%)", N - n_zero, pct_det),
  "Zero-count (dropout)",    as.character(n_zero),
  "Median reads/construct",  comma(round(q["50%"])),
  "Gini coefficient",        sprintf("%.3f", G),
  "Skew ratio (90th/10th)",  sprintf("%.1fx", skew9010),
  "Dynamic range (99th/1st)",sprintf("%.1fx", dyn),
  "Reads in top 10%",        sprintf("%.1f%%", top10sh))
write_csv(metrics, tbl("SuppFig2_plasmid_library_QC_metrics.csv"))

pD <- ggplot(metrics, aes(x=0, y=factor(Metric, levels=rev(Metric)))) +
  geom_text(aes(label=Metric), hjust=0, size=pt(FS_LABEL)) +
  geom_text(aes(x=1, label=Value), hjust=1, size=pt(FS_LABEL), fontface = "plain") +
  scale_x_continuous(limits=c(0,1), expand=expansion(mult=c(0.02,0.02))) +
  labs(title="QC metrics summary", x=NULL, y=NULL) +
  theme_void(base_size=MIN_FONT_PT) +
  theme(plot.title=element_text(face = "plain", size=FS_TITLE, hjust=0),
        plot.margin=margin(10,10,10,10))

save_panel_grid(list(pA, pB, pC, pD),
                supp_fig("SuppFig2J_plasmid_library_QC.pdf"), ncol = 2)

log_session("supp_01_plasmid_library_qc")
message("[supp_01] Saved SuppFig2J_plasmid_library_QC.{pdf,png} and metrics CSV")
print(metrics)
