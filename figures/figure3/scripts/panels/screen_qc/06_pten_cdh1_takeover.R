# =====================================================================
# PUBLICATION FIGURE PANEL
# PTEN*CDH1: rare at Day 0 -> dominant by Day 18 (in vivo pool takeover)
# Panel A: abundance-rank trajectory
# Panel B: pool-share takeover
# Source: 2nd_DC_counts.csv (Cell=Day0, D3, D18)
# NOTE: all numbers quoted in titles/captions below are computed from the
# data at run time (sprintf'd from hp/pair_tp) -- no hard-coded literals --
# so they stay correct if the input counts table changes on a future rerun.
# =====================================================================

if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
need <- c("readr","dplyr","tidyr","ggplot2","scales","patchwork","stringr")
invisible(lapply(need, library, character.only = TRUE))

HILITE <- "PTEN*CDH1"
crimson <- "#C0392B"; greyln <- "grey80"

raw <- read_csv(derived("screens","2nd_DC_counts.csv"), show_col_types = FALSE)
grp <- list(`Day 0` = grep("^Cell", names(raw), value = TRUE),
            `Day 3` = grep("^D3_",  names(raw), value = TRUE),
            `Day 18`   = grep("^D18_", names(raw), value = TRUE))

# pair-level pool fraction (%) and rank per timepoint
pair_tp <- lapply(names(grp), function(tp){
  cols <- grp[[tp]]
  raw |> transmute(gene, reads = rowSums(across(all_of(cols), as.numeric))) |>
    group_by(gene) |> summarise(reads = sum(reads), .groups = "drop") |>
    mutate(timepoint = tp,
           pct  = 100 * reads / sum(reads),
           rank = rank(-reads, ties.method = "min"))
}) |> bind_rows() |>
  mutate(timepoint = factor(timepoint, levels = c("Day 0","Day 3","Day 18")),
         hl = toupper(gene) == toupper(HILITE))

hp <- pair_tp |> filter(hl) |> arrange(timepoint)   # highlight path
ends <- hp |> filter(timepoint %in% c("Day 0","Day 18"))

# ---- computed values used to build self-consistent titles/captions ----
n_pairs   <- n_distinct(pair_tp$gene)
day0_rank <- hp$rank[hp$timepoint == "Day 0"]
d14_rank  <- hp$rank[hp$timepoint == "Day 18"]
d14_share <- hp$pct[hp$timepoint == "Day 18"]

base_theme <- theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 13),
        axis.title = element_text(size = 12),
        plot.margin = margin(8, 16, 8, 8))

# ---- Panel A: rank trajectory (1 = top) ----
pA <- ggplot() +
  geom_line(data = filter(pair_tp, !hl),
            aes(timepoint, rank, group = gene), colour = greyln,
            alpha = 0.35, linewidth = 0.4) +
  geom_line(data = hp, aes(timepoint, rank, group = gene),
            colour = crimson, linewidth = 1.6) +
  geom_point(data = hp, aes(timepoint, rank), colour = crimson, size = 3.4) +
  geom_text(data = ends, aes(timepoint, rank,
            label = paste0("rank ", rank, "/", n_pairs)),
            vjust = ifelse(ends$timepoint == "Day 0", 2.1, -1.3),
            colour = crimson, fontface = "bold", size = 4) +
  annotate("text", x = 1, y = 1, label = HILITE, hjust = -0.05, vjust = -0.9,
           fontface = "bold.italic", colour = crimson, size = 4.4) +
  scale_y_reverse(breaks = c(1, 50, 100, 200, 300, n_pairs),
                  expand = expansion(mult = c(0.07, 0.07))) +
  labs(x = NULL, y = "Abundance rank  (1 = most abundant)",
       title = sprintf("A   Rank trajectory: %d/%d to #%d", day0_rank, n_pairs, d14_rank)) +
  base_theme

# ---- Panel B: pool-share takeover ----
pB <- ggplot() +
  geom_line(data = filter(pair_tp, !hl),
            aes(timepoint, pct, group = gene), colour = greyln,
            alpha = 0.35, linewidth = 0.4) +
  geom_line(data = hp, aes(timepoint, pct, group = gene),
            colour = crimson, linewidth = 1.6) +
  geom_point(data = hp, aes(timepoint, pct), colour = crimson, size = 3.4) +
  geom_text(data = hp, aes(timepoint, pct,
            label = ifelse(pct < 1, sprintf("%.2f%%", pct), sprintf("%.1f%%", pct))),
            vjust = -1.2, hjust = c(0.1, 0.5, 0.6),
            colour = crimson, fontface = "bold", size = 4) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0.02, 0.13))) +
  labs(x = NULL, y = "Share of the screened pool",
       title = sprintf("B   One pair seizes ~%.1f%% of the population", d14_share)) +
  base_theme

fig <- pA + pB +
  plot_annotation(
    title = expression(paste(italic("Pten*Cdh1"), " expands from a rare clone to population-dominant in vivo")),
    caption = sprintf("In vivo combinatorial screen (n=%d pairs). Pool share = reads per pair / total reads per timepoint; ranks by abundance. Grey = all other pairs.", n_pairs),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.caption = element_text(size = 9, colour = "grey35", hjust = 0)))

ggsave(fig_supp("fig_pten_cdh1_takeover.pdf"), fig, width = 10, height = 5.2)
ggsave(fig_supp("fig_pten_cdh1_takeover.png"), fig, width = 10, height = 5.2, dpi = 600)
message("Saved fig_pten_cdh1_takeover.pdf / .png")
print(hp |> select(gene, timepoint, reads, pct, rank))
