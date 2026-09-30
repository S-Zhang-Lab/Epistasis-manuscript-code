# =====================================================================
# Five gene-pair abundance-rank trajectory panel (in vivo)
#   Pten*Cdh1, Pten*Nt, Pten*Cx3cl1, Pten*Cxcr5, Pten*Tlr7
#
# No endpoint labels; the lines fill a SQUARE panel (aspect.ratio = 1);
# Helvetica. All 324 pairs shown in grey; the five highlighted in colour
# with a colour legend. y = pool-abundance rank (1 = most abundant).
#
# Input : 2nd_DC_counts.csv  (Cell_1/2/3 = Day 0; D3_* = Day 3; D18_* = Day 18)
# Outputs: Fig3K_seed_trajectories.pdf / .png
#          output/tables/Fig3K_seed_trajectories_source_data.csv
#
# v7 repo overhaul 2026-07-15: relocated from XL_code/four_pair_trajectory_panel/.
# 2026-08 (XL screen-data update): (1) endpoint reverted to the true Day 18
#   (count columns relabelled D14_ -> D18_ per PI); (2) pair matching is now
#   case-insensitive via toupper(gene) -> tolerates the updated tables' mixed
#   case (Pten*Cdh1) where an exact "PTEN*CDH1" join would have drawn no lines;
#   (3) Pten*Nt added as a 5th reference trajectory.
# =====================================================================

if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))

need <- c("readr","dplyr","tidyr","ggplot2","scales","tibble","ragg")
invisible(lapply(need, library, character.only = TRUE))

FONT <- "Helvetica"

# highlighted pairs (key is upper-cased for matching) -> display -> colour (Okabe-Ito)
sel <- tribble(
  ~id,           ~disp,          ~col,
  "PTEN*CDH1",   "Pten_Cdh1",    "#D55E00",
  "PTEN*NT",     "Pten_Nt",      "#E69F00",
  "PTEN*CX3CL1", "Pten_Cx3cl1",  "#0072B2",
  "PTEN*CXCR5",  "Pten_Cxcr5",   "#009E73",
  "PTEN*TLR7",   "Pten_Tlr7",    "#CC79A7"
)

# ---- data: pool-abundance rank per gene pair at each timepoint --------
raw <- read_csv(derived("screens", "2nd_DC_counts.csv"), show_col_types = FALSE)
grp <- list(`Day 0`  = grep("^Cell", names(raw), value = TRUE),
            `Day 3`  = grep("^D3_",  names(raw), value = TRUE),
            `Day 18` = grep("^D18_", names(raw), value = TRUE))

pair_tp <- lapply(names(grp), function(tp){
  raw |> transmute(gene, reads = rowSums(across(all_of(grp[[tp]]), as.numeric))) |>
    group_by(gene) |> summarise(reads = sum(reads), .groups = "drop") |>
    mutate(timepoint = tp, rank = rank(-reads, ties.method = "min"))
}) |> bind_rows() |>
  mutate(timepoint = factor(timepoint, levels = c("Day 0","Day 3","Day 18")))
N <- n_distinct(pair_tp$gene)

hp  <- pair_tp |> mutate(.key = toupper(gene)) |>
  inner_join(sel, by = c(".key" = "id")) |>
  mutate(disp = factor(disp, levels = sel$disp))
pal <- setNames(sel$col, sel$disp)

# Numerical source data for every trajectory plotted in Fig. 3K. `reads` is
# the timepoint-specific replicate sum from which the plotted rank is derived.
write_csv(
  pair_tp |>
    mutate(.key = toupper(gene)) |>
    left_join(sel, by = c(".key" = "id")) |>
    transmute(
      gene_pair = gene,
      timepoint = as.character(timepoint),
      summed_reads = reads,
      abundance_rank = rank,
      trajectory_class = if_else(is.na(disp), "other pair", "highlighted pair"),
      display_label = disp,
      plot_colour = if_else(is.na(col), "grey82", col)
    ) |>
    arrange(gene_pair, factor(timepoint, levels = c("Day 0", "Day 3", "Day 18"))),
  tbl("Fig3K_seed_trajectories_source_data.csv")
)

# ---- plot: square 120pt data facet, 5 pt font floor, proportional geoms
p <- ggplot() +
  geom_line(data = filter(pair_tp, !toupper(gene) %in% sel$id),
            aes(timepoint, rank, group = gene), colour = "grey82",
            alpha = 0.35, linewidth = 0.25) +
  geom_line(data = hp, aes(timepoint, rank, group = disp, colour = disp), linewidth = 0.8) +
  geom_point(data = hp, aes(timepoint, rank, colour = disp), size = 1.5) +
  scale_y_reverse(breaks = c(1,50,100,200,300,N), expand = expansion(mult = c(.04,.04))) +
  scale_x_discrete(expand = expansion(mult = c(0.06, 0.06))) +   # minimal padding -> lines fill panel
  scale_colour_manual(values = pal, name = NULL) +
  labs(x = NULL, y = "Abundance rank  (1 = most abundant)",
       title = "Gene-pair rank trajectories (in vivo)",
       caption = "n = 324 pairs; grey = other pairs.") +
  theme_classic(base_size = 5, base_family = FONT) +
  theme(text            = element_text(family = FONT, colour = "black"),
        plot.title      = element_text(face = "bold", size = 6),
        axis.text       = element_text(size = 5, colour = "black"),
        axis.title      = element_text(size = 5),
        axis.line       = element_line(linewidth = 0.3, colour = "black"),
        axis.ticks      = element_line(linewidth = 0.3, colour = "black"),
        legend.position = "right",
        legend.title    = element_blank(),
        legend.text     = element_text(size = 5),
        legend.key.size = unit(3, "mm"),
        legend.margin   = margin(0,0,0,0),
        plot.caption    = element_text(size = 5, colour = "grey35", hjust = 0),
        plot.margin     = margin(3,3,3,3))

# ---- force the data facet to exactly 120 x 120 pt (excl. legend/axes) -
pdf(tempfile(fileext = ".pdf"), width = 12, height = 12)   # measuring device
gt <- ggplotGrob(p)
pidx <- gt$layout[gt$layout$name == "panel", , drop = FALSE]
gt$widths[pidx$l]  <- unit(120, "pt")
gt$heights[pidx$t] <- unit(120, "pt")
tw <- grid::convertWidth(sum(gt$widths),   "in", valueOnly = TRUE)
th <- grid::convertHeight(sum(gt$heights), "in", valueOnly = TRUE)
invisible(dev.off())
pdf(fig_main("Fig3K_seed_trajectories.pdf"), width = tw, height = th, family = FONT)
grid::grid.draw(gt); invisible(dev.off())
ragg::agg_png(fig_main("Fig3K_seed_trajectories.png"), width = tw, height = th, units = "in", res = 600)
grid::grid.draw(gt); invisible(dev.off())
if (file.exists("Rplots.pdf")) invisible(file.remove("Rplots.pdf"))
message(sprintf("panel = 120x120 pt; full figure = %.0f x %.0f pt", tw*72, th*72))
print(hp |> arrange(disp, timepoint) |> select(disp, timepoint, rank))
