if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))

need <- c("readr","dplyr","ggplot2","ggrepel","ragg")
invisible(lapply(need, library, character.only = TRUE))

FONT <- "Helvetica"
FDR  <- 0.15
TOP  <- 100

norm_key <- function(x) vapply(strsplit(tolower(x), "\\*"), function(g)
  paste(sort(g), collapse = "*"), character(1))

titlecase_pair <- function(x) sapply(strsplit(x, "\\*"), function(g)
  paste(vapply(g, function(t) paste0(toupper(substr(t,1,1)), tolower(substr(t,2,nchar(t)))),
               character(1)), collapse = "_"))
read_rank <- function(f) read_tsv(f, show_col_types = FALSE) |>
  transmute(key = norm_key(id), rank = `pos|rank`, fdr = `pos|fdr`)

d3  <- read_rank(derived("screens", "Day3.gene_summary.txt"))  |> rename(rank3 = rank,  fdr3 = fdr)
d18 <- read_rank(derived("screens", "Day18.gene_summary.txt")) |> rename(rank18 = rank, fdr18 = fdr)

cmp <- inner_join(d3, d18, by = "key") |>
  mutate(pair    = titlecase_pair(key),
         min_fdr = pmin(fdr3, fdr18),
         sig     = -log10(min_fdr),
         rlfc    = log2(rank3 / rank18)) |>
  arrange(sig)
sig_max <- max(cmp$sig)

hits  <- filter(cmp, min_fdr <  FDR, rank3 <= TOP, rank18 <= TOP)
noise <- filter(cmp, min_fdr >= FDR, rank3 <= TOP, rank18 <= TOP)
rlfc_lim <- max(abs(hits$rlfc))

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
LAB_PT <- 5/.pt

p <- ggplot(mapping = aes(rank18, rank3)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_point(data = noise, fill = "grey85", shape = 21, colour = "grey65",
             size = 1.1, alpha = 0.55, stroke = 0.2) +
  geom_point(data = hits, aes(fill = rlfc, size = sig), shape = 21, colour = "grey25",
             stroke = 0.35, alpha = 0.95) +
  geom_text_repel(data = slice_max(hits, sig, n = 10), aes(label = pair, colour = rlfc),
                  size = LAB_PT, family = FONT, fontface = "italic", box.padding = 0.55,
                  point.padding = 0.35, force = 14, force_pull = 0.15, max.overlaps = Inf,
                  min.segment.length = 0, max.iter = 1e5, max.time = 2, segment.size = 0.12,
                  segment.colour = "grey65", seed = 1) +
  scale_x_reverse(breaks = c(1,25,50,75,100)) +
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

pdf(tempfile(fileext = ".pdf"), width = 15, height = 15)
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
if (file.exists("Rplots.pdf")) invisible(file.remove("Rplots.pdf"))
message(sprintf("panel = 120x120 pt; full figure = %.0f x %.0f pt", tw*72, th*72))
