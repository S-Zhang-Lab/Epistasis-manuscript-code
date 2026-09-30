# =====================================================================
# Representative gene-pair enrichment across the time course (Day0/3/18)
# Narrative: candidates are restricted to pairs that actually cross the
# FDR<0.1 significance line at some timepoint in the current MAGeCK
# tables (verified across all 286 Day-18 pairs: only Pten*Cdh1, Pten*Nt,
# Pten*Muc16 do so -- no non-Pten pair reaches significance at Day 18).
# Per user request, the figure shows only 2 of those 3 (Pten*Nt dropped
# to keep a clean 2-line comparison; it remains a real Day18 hit at
# fdr=0.009, just not plotted here):
#   - Pten*Cdh1:  significant Day 3 (fdr=0.0016) AND Day 18 (fdr=0.0007)
#     -> persistent hit.
#   - Pten*Muc16: NOT significant Day 3 (fdr=0.188) but significant
#     Day 18 (fdr=0.0007) -> true emergent hit (rank 7 -> rank 4/286).
#
# Panel A - dot plot : pair x timepoint; colour = LFC, size = -log10(FDR)
# Panel B - significance trajectory : -log10(pos FDR) over time
#
# NOTE: Day3 and Day18 are separate MAGeCK runs (each median-normalised),
# so raw LFC magnitude is NOT strictly comparable across timepoints.
# Rank / FDR / goodsgrna are the robust cross-timepoint metrics.
#
# Day 0 baseline (redesigned): there is no MAGeCK test at Day 0 (no
# enrichment vs. itself), so Day 0 is NOT drawn as a fabricated
# lfc=0/fdr=1 point. Instead:
#   - colour (LFC-like) = log2(observed Day-0 pool share / uniform
#     expectation 1/n_pairs), i.e. how over/under-represented each pair
#     already was in the starting plasmid->cell pool, from real
#     2nd_DC_counts.csv reads.
#   - size (significance) is left undefined at Day 0: Panel A draws the
#     Day-0 dot as a fixed-size open diamond (shape 18) instead of an
#     FDR-scaled circle, and Panel B's significance-trajectory lines
#     start at Day 3 (no Day-0 point plotted at all).
# =====================================================================

if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
need <- c("readr", "dplyr", "tidyr", "ggplot2", "patchwork", "scales")
invisible(lapply(need, library, character.only = TRUE))

# Robust to being launched either from the parent directory or from inside
# final_rerun/ itself (e.g. via Rscript run from within the folder).

# ---- 1. read both files; key on uppercase id (files differ in casing) ----
read_tp <- function(path, day) {
  read_tsv(path, show_col_types = FALSE) |>
    transmute(key = toupper(id),
              lfc = `neg|lfc`,
              fdr = `pos|fdr`,
              rank = `pos|rank`,
              good = `pos|goodsgrna`,
              timepoint = day)
}
d3  <- read_tp(derived("screens","Day3.gene_summary.txt"),  "Day 3")
d18 <- read_tp(derived("screens","Day18.gene_summary.txt"), "Day 18")

# ---- 2. featured pairs (edit here to swap) ----
# NOTE (this rerun, revised): earlier drafts of this figure featured
# several pairs (Apc*Nt, Kit*Nt, Clca2*Mpeg1, Krt15*Cdh1) chosen only by
# Day-18 rank improvement, without a significance filter -- none of them
# ever reach pos|fdr < 0.1 (all sit at fdr 0.77-0.96 at Day 18), so their
# trajectory lines in Panel B never cross the FDR=0.1 reference line.
# Per user request, the featured set is now restricted to pairs that
# ACTUALLY cross FDR<0.1 at some timepoint (checked against all 286
# pairs in the current Day18.gene_summary.txt -- exactly 3 qualify, all
# Pten pairs; no non-Pten pair reaches Day-18 significance):
sel <- tribble(
  ~key,            ~pair,          ~class,
  "PTEN*CDH1",     "Pten*Cdh1",    "Persistent (Day 3 + Day 18 hit)",
  "PTEN*MUC16",    "Pten*Muc16",   "Emerges by Day 18 (not sig. Day 3)"
)

tp <- bind_rows(d3, d18) |> inner_join(sel, by = "key")

# ---- 2b. Day 0 baseline (redesigned): real pool-share, no FDR, no LFC ----
# There is no MAGeCK enrichment test at Day 0 (nothing to compare it to
# yet), so we do NOT fabricate lfc=0/fdr=1. We DO compute, from the real
# Day-0 (Cell_*) read counts in 2nd_DC_counts.csv, how over/under-
# represented each featured pair already was in the starting pool
# relative to a perfectly uniform library (1/n_pairs share each) -- but
# this "pool_lfc" quantity is NOT the same thing as the Day3/Day18
# MAGeCK lfc (one compares real reads to a uniform-library assumption;
# the other compares two real timepoints via MAGeCK's own normalisation)
# and must NOT be plotted on the same colour scale/legend as those, or
# it reads as if all timepoints show one consistent "LFC vs Day 0"
# quantity. So Day 0 keeps its own separate column (pool_lfc, pool_pct)
# used only for the point label text below -- the Day-0 point itself is
# drawn as a neutral grey diamond with NO colour-scale mapping at all.
counts <- read_csv(derived("screens","2nd_DC_counts.csv"), show_col_types = FALSE)
day0_cols <- grep("^Cell", names(counts), value = TRUE)
n_pairs_lib <- n_distinct(counts$gene)

day0_share <- counts |>
  transmute(key = toupper(gene), reads = rowSums(across(all_of(day0_cols), as.numeric))) |>
  group_by(key) |> summarise(reads = sum(reads), .groups = "drop") |>
  mutate(pool_pct = 100 * reads / sum(reads),
         pool_lfc = log2(pool_pct/100 / (1 / n_pairs_lib)))   # vs uniform expectation; NOT comparable to MAGeCK lfc

day0 <- sel |>
  left_join(day0_share, by = "key") |>
  mutate(lfc = NA_real_, fdr = NA_real_, rank = NA, good = NA, timepoint = "Day 0") |>
  select(key, pair, class, lfc, fdr, rank, good, timepoint, pool_pct, pool_lfc)

dat <- bind_rows(day0, tp |> mutate(pool_pct = NA_real_, pool_lfc = NA_real_)) |>
  mutate(timepoint = factor(timepoint, levels = c("Day 0", "Day 3", "Day 18")),
         class     = factor(class, levels = unique(sel$class)),
         neglog10fdr = -log10(fdr),   # NA at Day 0 -- intentional, not plotted as a sized point
         # y-axis order: by Day-18 significance, grouped by class
         pair = factor(pair, levels = sel$pair[order(sel$class, decreasing = TRUE)]))

# ---- 3. Panel A: dot plot ----
# Day 0 is drawn as a FIXED-SIZE, UNCOLOURED (neutral grey) open diamond
# (shape 18): it carries no colour-scale value at all, because there is
# no MAGeCK test at Day 0 and its "pool_lfc" (vs. a uniform-library
# assumption) is not the same quantity as the Day3/Day18 MAGeCK lfc --
# putting it on the same colour legend would wrongly imply one
# consistent "LFC vs Day 0" scale across timepoints. Instead,
# the real Day-0 read-share signal is shown as a text label (% of pool)
# next to the grey diamond. Day 3 / Day 18 keep the original FDR-sized,
# lfc-coloured circles -- colour there is a real, well-defined MAGeCK lfc.
pA <- ggplot(dat, aes(timepoint, pair)) +
  geom_point(data = dat |> filter(timepoint != "Day 0"),
             aes(colour = lfc, size = neglog10fdr)) +
  geom_point(data = dat |> filter(timepoint == "Day 0"),
             size = 3.6, shape = 18, colour = "grey55", fill = "grey55") +
  geom_text(data = dat |> filter(timepoint == "Day 0"),
            aes(label = sprintf("%.2f%%", pool_pct)),
            nudge_y = 0.32, size = 2.9, colour = "grey35") +
  facet_grid(class ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_gradient2(low = "#3B6FB6", mid = "grey92", high = "#CC3311",
                         midpoint = 0, limits = c(-2, 6), oob = squish,
                         name = "LFC\n(MAGeCK, Day3/Day18)") +
  scale_size_continuous(range = c(1.5, 10), name = expression(-log[10]~FDR),
                        breaks = c(0, 1, 2, 3)) +
  scale_x_discrete(limits = levels(dat$timepoint)) +
  labs(x = NULL, y = NULL,
       title = "A. Two pairs cross FDR < 0.1",
       caption = paste0(
         paste(strwrap("One persistent hit (Day 3 + Day 18) and one Day-18-only emergent hit. Colour capped at LFC=6; LFC not strictly comparable across the Day3/Day18 MAGeCK runs -- size (FDR) is the robust metric.", width = 78), collapse = "\n"),
         "\n",
         paste(strwrap("Day 0 (grey diamond) carries no MAGeCK lfc/FDR -- no test exists at Day 0. Label = pair's real share (%) of total Day-0 pool reads (2nd_DC_counts.csv), for reference only, not on the colour scale.", width = 78), collapse = "\n"),
         "\n",
         paste(strwrap("All 286 Day-18 pairs were screened for pos|fdr<0.1; only 3 (all Pten pairs) qualify -- 2 are shown here (Pten*Nt, also significant, omitted for clarity).", width = 78), collapse = "\n"))) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        strip.placement = "outside",
        panel.grid.minor = element_blank())

# ---- 4. Panel B: significance trajectory ----
# Day 0 has no FDR by construction, so it is intentionally omitted here
# (lines start at Day 3) rather than plotting a fabricated significance
# value.
pB <- ggplot(dat |> filter(timepoint != "Day 0"),
             aes(timepoint, neglog10fdr, group = pair, colour = class)) +
  geom_hline(yintercept = -log10(0.10), linetype = "dashed", colour = "grey50") +
  annotate("text", x = 0.6, y = -log10(0.10) + 0.12, label = "FDR = 0.1",
           size = 3, colour = "grey40", hjust = 0) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.4) +
  geom_label(
    data = dat |> filter(timepoint == "Day 18") |>
      mutate(ylab = neglog10fdr + c(0.12, -0.12)[match(pair, c("Pten*Cdh1", "Pten*Muc16"))]),
    aes(y = ylab, label = pair, colour = class),
    hjust = 0, nudge_x = 0.12, size = 3.2, label.size = 0,
    label.padding = unit(0.1, "lines"), fill = "white", show.legend = FALSE) +
  scale_colour_manual(values = c("Persistent (Day 3 + Day 18 hit)" = "#CC3311",
                                 "Emerges by Day 18 (not sig. Day 3)" = "#1B7837"), name = NULL) +
  scale_x_discrete(limits = c("Day 3", "Day 18"), expand = expansion(mult = c(0.08, 0.75))) +
  labs(x = NULL, y = expression(-log[10]~(pos~FDR)),
       title = "B. Significance trajectory: both cross FDR = 0.1",
       caption = paste(strwrap("Day 0 omitted: no MAGeCK significance test exists at Day 0 (see Panel A for the real Day-0 pool-share signal).", width = 78), collapse = "\n")) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top",
        panel.grid.minor = element_blank())

fig <- pA / pB + plot_layout(heights = c(1, 1))
ggsave(fig_supp("representative_pairs.pdf"), fig, width = 8.5, height = 9)
ggsave(fig_supp("representative_pairs.png"), fig, width = 8.5, height = 9, dpi = 300)

message("Saved representative_pairs.pdf / .png")
dat |> arrange(class, pair, timepoint) |>
  select(pair, class, timepoint, lfc, fdr, rank, good) |> print(n = Inf)
