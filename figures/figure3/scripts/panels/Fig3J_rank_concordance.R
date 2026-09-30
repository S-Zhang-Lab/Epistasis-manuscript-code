# =====================================================================
# Rank concordance panel: Day 3 vs Day 18 enrichment rank (in vivo)
#
# Data facet forced to EXACTLY 120 x 120 pt (excl. legend); Helvetica; 5 pt
# font floor; dots sized proportional to the panel. Built at FINAL placed
# size so no rescaling is needed in Illustrator (lab convention, R/fig3_theme.R).
#   x / y   = Day 18 / Day 3 positive-selection rank (MAGeCK), top 100
#   colour  = rank log2FC  ->  purple = stable, blue = early, orange = late
#   size    = -log10(min FDR across the two timepoints)  [significance]
#   labels  = up to 10 pairs with the strongest FDR (9 in the current top-100 overlap)
#
# Inputs : Day3.gene_summary.txt, Day18.gene_summary.txt  (MAGeCK gene_summary)
# Outputs: Fig3J_rank_concordance.pdf  (vector, panel = 150x150 pt)
#          Fig3J_rank_concordance.png  (600 dpi preview)
#          output/tables/Fig3J_rank_concordance_source_data.csv
#
# v7 repo overhaul 2026-07 : relocated from XL_code/rank_concordance_panel/.
# 2026-08 (XL screen-data update): reverted to the true Day 18 endpoint
#   (the earlier "Day 14" labelling was corrected to Day 18, per PI); the
#   underlying analysis is unchanged. Case-insensitive pair key already
#   tolerates the PTEN*CDH1 -> Pten*Cdh1 rename in the updated tables.
# =====================================================================

if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))

need <- c("readr","dplyr","ggplot2","ggrepel","ragg")
invisible(lapply(need, library, character.only = TRUE))

FONT <- "Helvetica"
FDR  <- 0.15     # highlight threshold (colour/label pairs below this)
TOP  <- 100      # keep only pairs ranked within the top 100 at both timepoints (zoom on the hits)

# ---- helpers ----------------------------------------------------------
# unordered, lowercase gene-pair key so "A*B" and "b*a" match across files
norm_key <- function(x) vapply(strsplit(tolower(x), "\\*"), function(g)
  paste(sort(g), collapse = "*"), character(1))
# display name: mouse-style Title case around the "*"
titlecase_pair <- function(x) sapply(strsplit(x, "\\*"), function(g)
  paste(vapply(g, function(t) paste0(toupper(substr(t,1,1)), tolower(substr(t,2,nchar(t)))),
               character(1)), collapse = "_"))   # display "_" (house style); the join key keeps "*"
read_rank <- function(f) read_tsv(f, show_col_types = FALSE) |>
  transmute(key = norm_key(id), rank = `pos|rank`, fdr = `pos|fdr`)

# ---- data -------------------------------------------------------------
d3  <- read_rank(derived("screens", "Day3.gene_summary.txt"))  |> rename(rank3 = rank,  fdr3 = fdr)
d18 <- read_rank(derived("screens", "Day18.gene_summary.txt")) |> rename(rank18 = rank, fdr18 = fdr)

cmp <- inner_join(d3, d18, by = "key") |>
  mutate(pair    = titlecase_pair(key),
         min_fdr = pmin(fdr3, fdr18),          # strongest evidence at either timepoint
         sig     = -log10(min_fdr),            # -> point size
         rlfc    = log2(rank3 / rank18)) |>    # rank log2FC: + gained by D18, - higher at D3
  arrange(sig)
sig_max <- max(cmp$sig)

hits  <- filter(cmp, min_fdr <  FDR, rank3 <= TOP, rank18 <= TOP)   # enriched, coloured
noise <- filter(cmp, min_fdr >= FDR, rank3 <= TOP, rank18 <= TOP)   # grey background
labels <- slice_max(hits, sig, n = 10)
rlfc_lim <- max(abs(hits$rlfc))                                     # symmetric colour limits

# Numerical source data for every point and label actually shown in Fig. 3J.
plot_data <- bind_rows(
  mutate(noise, point_class = "not significant"),
  mutate(hits,  point_class = "FDR < 0.15")
) |>
  mutate(labelled = key %in% labels$key) |>
  arrange(rank18, rank3)

write_csv(
  plot_data |>
    transmute(
      pair,
      pair_key = key,
      day3_positive_selection_rank = rank3,
      day18_positive_selection_rank = rank18,
      day3_positive_selection_fdr = fdr3,
      day18_positive_selection_fdr = fdr18,
      minimum_fdr = min_fdr,
      negative_log10_minimum_fdr = sig,
      rank_log2_fold_change = rlfc,
      point_class,
      labelled
    ),
  tbl("Fig3J_rank_concordance_source_data.csv")
)

# ---- theme: Helvetica, smallest text = 5 pt (lab floor) ---------------
compact <- theme_bw(base_size = 5, base_family = FONT) +
  theme(text             = element_text(family = FONT),
        plot.title       = element_text(face = "bold", size = 6),
        plot.subtitle    = element_text(size = 5),
        plot.caption     = element_text(hjust = 0, size = 5, colour = "grey35"),
        axis.title       = element_text(size = 5),
        axis.text        = element_text(size = 5),
        legend.title     = element_text(size = 5),
        legend.text      = element_text(size = 5),
        legend.key.size  = unit(2.8, "mm"),
        legend.margin    = margin(1,1,1,1),
        panel.grid.minor = element_blank(),
        plot.margin      = margin(2,2,2,2))
LAB_PT <- 5/.pt   # geom_text 'size' that renders at 5 pt

# ---- plot -------------------------------------------------------------
p <- ggplot(mapping = aes(rank18, rank3)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_point(data = noise, fill = "grey85", shape = 21, colour = "grey65",
             size = 1.1, alpha = 0.55, stroke = 0.2) +
  geom_point(data = hits, aes(fill = rlfc, size = sig), shape = 21, colour = "grey25",
             stroke = 0.35, alpha = 0.95) +
  geom_text_repel(data = labels, aes(label = pair, colour = rlfc),
                  size = LAB_PT, family = FONT, fontface = "italic", box.padding = 0.55,
                  point.padding = 0.35, force = 14, force_pull = 0.15, max.overlaps = Inf,
                  min.segment.length = 0, max.iter = 1e5, max.time = 2, segment.size = 0.12,
                  segment.colour = "grey65", seed = 1) +
  scale_x_reverse(breaks = c(1,25,50,75,100)) +   # rank 1 (most enriched) at top-right
  scale_y_reverse(breaks = c(1,25,50,75,100)) +
  scale_fill_gradient2(low = "#3690C0", mid = "#6A1B9A", high = "#E08214", midpoint = 0,
                       limits = c(-rlfc_lim, rlfc_lim), name = "Rank log2FC\n(0 = stable)") +
  scale_colour_gradient2(low = "#2A6F9E", mid = "#511578", high = "#B5670F", midpoint = 0,
                         limits = c(-rlfc_lim, rlfc_lim), guide = "none") +
  scale_size_continuous(range = c(1.8, 5), limits = c(0, sig_max),
                        name = expression(-log[10]~FDR), breaks = c(1,2,3)) +
  labs(x = "Day 18 enrichment rank", y = "Day 3 enrichment rank",
       title = "Rank concordance: D3 vs D18 (top 100)",
       subtitle = "FDR<0.15. Size = significance; colour: stable / early / late",
       caption = "purple = stable, blue = early, orange = late; grey = not significant.") +
  compact

# ---- force PANEL (data area) to exactly 120 x 120 pt ------------------
# (legend / axes / titles are laid out AROUND it, so the full page is larger)
pdf(tempfile(fileext = ".pdf"), width = 15, height = 15)   # measuring device (open before grob)
gt <- ggplotGrob(p)
pidx <- gt$layout[gt$layout$name == "panel", , drop = FALSE]
gt$widths[pidx$l]  <- unit(120, "pt")
gt$heights[pidx$t] <- unit(120, "pt")
tw <- grid::convertWidth(sum(gt$widths),   "in", valueOnly = TRUE)
th <- grid::convertHeight(sum(gt$heights), "in", valueOnly = TRUE)
invisible(dev.off())
pdf(fig_main("Fig3J_rank_concordance.pdf"), width = tw, height = th, family = FONT)
grid::grid.draw(gt); invisible(dev.off())
ragg::agg_png(fig_main("Fig3J_rank_concordance.png"), width = tw, height = th, units = "in", res = 600)
grid::grid.draw(gt); invisible(dev.off())
if (file.exists("Rplots.pdf")) invisible(file.remove("Rplots.pdf"))   # tidy any stray default device
message(sprintf("panel = 120x120 pt; full figure = %.0f x %.0f pt", tw*72, th*72))
