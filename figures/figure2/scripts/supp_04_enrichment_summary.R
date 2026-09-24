pkgs <- c("readr", "dplyr", "tidyr", "ggplot2", "ggridges", "patchwork", "scales")
inst <- rownames(installed.packages())
if (length(setdiff(pkgs, inst)))
  install.packages(setdiff(pkgs, inst), repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(ggridges); library(patchwork); library(scales)
})

source(here::here("R", "paths.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "theme_pub.R"))
dir.create(OUT_SUPP, recursive = TRUE, showWarnings = FALSE)

day_levels <- c("Day 0", "Day 3", "Day 18")
fdr_cut <- 0.05

read_tp <- function(path, day) {
  read_tsv(path, show_col_types = FALSE) |>
    transmute(construct = id, lfc = `neg|lfc`, fdr = `pos|fdr`, timepoint = day)
}

d3  <- read_tp(derived("R2_Day3.gene_summary.txt"),  "Day 3")
d18 <- read_tp(derived("R2_Day18.gene_summary.txt"), "Day 18")
d0  <- d3 |> mutate(lfc = 0, fdr = NA_real_, timepoint = "Day 0")

dat <- bind_rows(d0, d3, d18) |>
  mutate(timepoint = factor(timepoint, levels = day_levels))

cat_levels <- c("Depleted (LFC < -1)", "Neutral (-1 to 0.5)",
                "Mildly enriched (0.5-1)", "Strongly enriched (LFC > 1)")
dat <- dat |>
  mutate(category = case_when(
    lfc <  -1   ~ cat_levels[1],
    lfc <   0.5 ~ cat_levels[2],
    lfc <=  1   ~ cat_levels[3],
    TRUE        ~ cat_levels[4]) |> factor(levels = cat_levels))

pal <- c("Depleted (LFC < -1)"         = "#3B6FB6",
         "Neutral (-1 to 0.5)"         = "#D9D9D9",
         "Mildly enriched (0.5-1)"     = "#F6B26B",
         "Strongly enriched (LFC > 1)" = "#CC3311")

prop <- dat |> count(timepoint, category, name = "n") |>
  group_by(timepoint) |> mutate(prop = n / sum(n)) |> ungroup()

panel_I <- ggplot(prop, aes(timepoint, prop, fill = category)) +
  geom_col(width = 0.7, colour = "white", linewidth = LW_THIN) +
  scale_y_continuous(labels = percent_format(), expand = expansion(c(0, 0.02))) +
  scale_fill_manual(values = pal, name = NULL) +
  labs(x = NULL, y = "% of crRNA constructs",
       title = "Enrichment categories redistribute between day 3 and day 18") +
  theme_panel() +
  theme(legend.position = "right")

save_panel(panel_I, supp_fig("SuppFig2L_enrichment_proportion.pdf"))

dat_rl <- dat |> mutate(timepoint = factor(timepoint, levels = rev(day_levels)))
rl_stats <- dat_rl |>
  group_by(timepoint) |>
  summarise(med   = median(lfc),
            p_enr = mean(lfc > 1),
            p90   = quantile(lfc, 0.90),
            n_sig = sum(fdr < fdr_cut, na.rm = TRUE),
            .groups = "drop") |>
  mutate(med_lab  = sprintf("median = %.2f", med),
         tail_lab = sprintf("%.1f%% with LFC > 1\n90th pct LFC = %.2f\n%d pairs FDR < %.2f",
                            100 * p_enr, p90, n_sig, fdr_cut))

ridge <- ggplot(dat_rl, aes(x = lfc, y = timepoint, fill = timepoint)) +
  geom_density_ridges(scale = 1.5, alpha = 0.85, colour = "white",
                      quantile_lines = TRUE, quantiles = 2, vline_width = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_text(data = rl_stats, inherit.aes = FALSE,
            aes(x = med, y = timepoint, label = med_lab),
            nudge_y = 0.10, hjust = -0.05, vjust = 0, size = pt(FS_LABEL), colour = "grey20") +
  geom_text(data = rl_stats, inherit.aes = FALSE,
            aes(x = 2.3, y = timepoint, label = tail_lab),
            nudge_y = 0.30, hjust = 0, vjust = 0, size = pt(FS_LABEL),
            lineheight = 0.95, fontface = "plain") +
  coord_cartesian(xlim = c(-3, 5)) +
  scale_fill_manual(values = c("Day 0"  = "#BBBBBB",
                               "Day 3"  = "#88A0C8",
                               "Day 18" = "#CC3311"), guide = "none") +
  labs(x = "Log2 fold change vs Day 0", y = NULL,
       title = "Bulk shifts down while the positive tail extends") +
  theme_panel() +
  theme(plot.title = element_text(face = "plain", size = 12))

save_panel(ridge, supp_fig("SuppFig2M_enrichment_ridgeline.pdf"))

log_session("supp_04_enrichment_summary")
message("[supp_04] Saved SuppFig2L_enrichment_proportion.{pdf,png} (+ ridgeline extra)")
print(rl_stats |> select(timepoint, med, p_enr, p90, n_sig))
