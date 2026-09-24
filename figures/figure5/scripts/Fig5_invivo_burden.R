suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "theme_pub.R"))
  source(here::here("R", "palettes.R"))
  source(here::here("R", "helpers.R"))
  for (pkg in c("ggplot2", "dplyr", "tidyr", "readr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
})

WT_LAB  <- "Mrc1-B6 (WT)"
DTR_LAB <- "Mrc1-DTR"

counts_csv <- file.path(DATA_DIR, "Tumor_Burden_Quantification",
                        "metastasis_counts_summary.csv")
quant <- readr::read_csv(counts_csv, show_col_types = FALSE) |>
  tidyr::pivot_longer(c("Depleted_count", "Control_WT_count"),
                      names_to = "arm", values_to = "foci") |>
  dplyr::mutate(
    host = factor(ifelse(arm == "Depleted_count", DTR_LAB, WT_LAB),
                  levels = c(WT_LAB, DTR_LAB))
  ) |>
  dplyr::select(host, mouse, foci) |>
  dplyr::arrange(host, mouse)

readr::write_csv(quant, tbl("tumor_burden_foci.csv"))

wt  <- quant$foci[quant$host == WT_LAB]
dtr <- quant$foci[quant$host == DTR_LAB]
wtest <- wilcox.test(dtr, wt, alternative = "two.sided", exact = TRUE)
pval  <- wtest$p.value
fold  <- mean(dtr) / mean(wt)
message(sprintf("[Fig5I] WT median=%.0f  DTR median=%.0f  fold(mean)=%.1fx  Mann-Whitney exact p=%.4f",
                median(wt), median(dtr), fold, pval))

summ <- quant |>
  dplyr::group_by(host) |>
  dplyr::summarize(mean = mean(foci), sem = sd(foci) / sqrt(dplyr::n()), .groups = "drop")

fills <- c(HOST_COLORS[["WT"]], HOST_COLORS[["Mrc1"]]); names(fills) <- c(WT_LAB, DTR_LAB)
ymax  <- max(quant$foci)
ytop  <- ymax * 1.20
brk_y <- ymax * 1.06

p <- ggplot(quant, aes(host, foci)) +
  geom_errorbar(data = summ, aes(x = host, ymin = mean - sem, ymax = mean + sem, color = host),
                width = 0.16, linewidth = 0.45, inherit.aes = FALSE) +
  geom_segment(data = summ, aes(x = as.integer(host) - 0.24, xend = as.integer(host) + 0.24,
                                y = mean, yend = mean, color = host),
               linewidth = 0.6, inherit.aes = FALSE) +
  geom_jitter(aes(fill = host), width = 0.11, height = 0, size = 1.7,
              shape = 21, stroke = 0.35, color = "gray20") +
  annotate("segment", x = 1, xend = 2, y = brk_y, yend = brk_y, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1, y = brk_y, yend = brk_y - ymax * 0.02, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2, y = brk_y, yend = brk_y - ymax * 0.02, linewidth = 0.4) +
  annotate("text", x = 1.5, y = ymax * 1.13,
           label = sprintf("Mann–Whitney P = %.4f", pval), size = 2.35) +
  scale_fill_manual(values = fills, guide = "none") +
  scale_colour_manual(values = fills, guide = "none") +
  scale_x_discrete(labels = c(setNames("Mrc1-B6\n(WT)", WT_LAB),
                              setNames("Mrc1-DTR", DTR_LAB))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0)),
                     breaks = seq(0, 3000, 1000)) +
  coord_cartesian(ylim = c(0, ytop), clip = "off") +
  labs(title = "GFP+ lung metastatic burden",
       subtitle = sprintf("Pten×Cdh1 · %.1f× · n = 6/arm", fold),
       x = NULL, y = "GFP+ metastatic foci per lung") +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 5.8, color = "gray25"))

save_pdf(p, fig("Fig5_invivo_burden.pdf"), 57, 75)
log_session("Fig5_invivo_burden")
