# =====================================================================
# Time-course CRISPR screen summary - VERSION 2
# Panel A: 100% stacked proportion bar
# Panel B: ridgeline with richer annotations
#          (median LFC + % LFC>1 + 90th-pct LFC + # FDR-sig pairs)
# No ECDF panel in this version.
# =====================================================================

if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
need <- c("readr", "dplyr", "tidyr", "ggplot2", "ggridges", "patchwork", "scales")
invisible(lapply(need, library, character.only = TRUE))

day_levels <- c("Day 0", "Day 3", "Day 18")
fdr_cut <- 0.05

# ---- 1. read (keep LFC and pos|fdr) ----
read_tp <- function(path, day) {
  read_tsv(path, show_col_types = FALSE) |>
    transmute(construct = id, lfc = `neg|lfc`, fdr = `pos|fdr`, timepoint = day)
}
d3  <- read_tp(derived("screens","Day3.gene_summary.txt"),  "Day 3")
d18 <- read_tp(derived("screens","Day18.gene_summary.txt"), "Day 18")
d0  <- d3 |> mutate(lfc = 0, fdr = NA_real_, timepoint = "Day 0")   # Day0 = baseline

dat <- bind_rows(d0, d3, d18) |>
  mutate(timepoint = factor(timepoint, levels = day_levels))

# ---- 2. enrichment categories (Panel A) ----
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

# ---- 3. Panel A ----
prop <- dat |> dplyr::count(timepoint, category, name = "n") |>
  group_by(timepoint) |> mutate(prop = n / sum(n)) |> ungroup()

pA <- ggplot(prop, aes(timepoint, prop, fill = category)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  scale_y_continuous(labels = percent_format(), expand = expansion(c(0, 0.02))) +
  scale_fill_manual(values = pal, name = NULL) +
  labs(x = NULL, y = "% of crRNA constructs",
       title = "A. crRNA proportion shifts toward enrichment over time") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 12))

# ---- 4. Panel B with richer annotations ----
dat_rl <- dat |> mutate(timepoint = factor(timepoint, levels = rev(day_levels)))

rl_stats <- dat_rl |>
  group_by(timepoint) |>
  summarise(med   = median(lfc),
            p_enr = mean(lfc > 1),
            p90   = quantile(lfc, 0.90),
            n_sig = sum(fdr < fdr_cut, na.rm = TRUE),  # Day0 fdr=NA -> 0
            .groups = "drop") |>
  mutate(med_lab  = sprintf("median = %.2f", med),
         tail_lab = sprintf("%.1f%% with LFC > 1\n90th pct LFC = %.2f\n%d pairs FDR < %.2f",
                            100 * p_enr, p90, n_sig, fdr_cut))

pB <- ggplot(dat_rl, aes(x = lfc, y = timepoint, fill = timepoint)) +
  geom_density_ridges(scale = 1.5, alpha = 0.85, colour = "white",
                      quantile_lines = TRUE, quantiles = 2, vline_width = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_text(data = rl_stats, inherit.aes = FALSE,
            aes(x = med, y = timepoint, label = med_lab),
            nudge_y = 0.10, hjust = -0.05, vjust = 0, size = 3.0, colour = "grey20") +
  geom_text(data = rl_stats, inherit.aes = FALSE,
            aes(x = 2.3, y = timepoint, label = tail_lab),
            nudge_y = 0.30, hjust = 0, vjust = 0, size = 3.0,
            lineheight = 0.95, fontface = "bold") +
  coord_cartesian(xlim = c(-3, 5)) +    # clip extreme Day3 outlier (LFC~10)
  scale_fill_manual(values = c("Day 0"  = "#BBBBBB",
                               "Day 3"  = "#88A0C8",
                               "Day 18" = "#CC3311"), guide = "none") +
  labs(x = "Log2 fold change vs Day 0", y = NULL,
       title = "B. LFC distribution widens / shifts positive over time") +
  theme_ridges(center_axis_labels = TRUE) +
  theme(plot.title = element_text(face = "bold", size = 12))

# ---- 5. combine (A + B only) & save ----
fig <- pA / pB + plot_layout(heights = c(1, 1.2))
ggsave(fig_supp("crRNA_enrichment_summary_v2.pdf"), fig, width = 8, height = 8)
ggsave(fig_supp("crRNA_enrichment_summary_v2.png"), fig, width = 8, height = 8, dpi = 300)

message("Saved crRNA_enrichment_summary_v2.pdf / .png")
print(rl_stats |> select(timepoint, med, p_enr, p90, n_sig))
