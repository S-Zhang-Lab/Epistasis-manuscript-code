suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

source(here::here("R", "paths.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "theme_pub.R"))

e <- read.csv(tbl("Fig2I_strong_GI_edges.csv"), stringsAsFactors = FALSE)
n <- read.csv(tbl("Fig2I_node_strong_GI_partner_counts.csv"), stringsAsFactors = FALSE)

tc <- function(x) paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))
e$A <- tc(e$gene_A); e$B <- tc(e$gene_B); n$g <- tc(n$gene)

left  <- unique(e$A); right <- unique(e$B)
deg   <- setNames(n$strong_gi_partners, n$g)
un    <- setdiff(n$g, c(left, right))
side  <- function(gs) tibble(g = gs, d = deg[gs]) |> arrange(desc(d), g)

L <- side(left)  |> mutate(x = 0, y = -seq_len(n()))
R <- side(right) |> mutate(x = 1, y = -seq_len(n()))
gi <- read.csv(tbl("R2_GI_derived.csv"), check.names = FALSE)
p1 <- unique(tc(gi$Gene1[!is.na(gi$Gene1)])); p1 <- p1[p1 != "Nt"]
UL <- tibble(g = un[un %in% p1], d = 0) |> arrange(g) |> mutate(x = 0, y = -nrow(L) - seq_len(n()) - 0.6)
UR <- tibble(g = un[!un %in% p1], d = 0) |> arrange(g) |> mutate(x = 1, y = -nrow(R) - seq_len(n()) - 0.6)
V  <- bind_rows(L, R, UL, UR)

ed <- e |>
  left_join(V |> select(g, xa = x, ya = y), by = c("A" = "g")) |>
  left_join(V |> select(g, xb = x, yb = y), by = c("B" = "g"))

p <- ggplot() +
  geom_segment(data = ed,
               aes(x = xa, y = ya, xend = xb, yend = yb,
                   colour = GI_Zscore, linewidth = abs_GI_Z), alpha = 0.75) +
  scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "GI Z-score") +
  scale_linewidth(range = c(0.2, 1.6), guide = "none") +
  geom_point(data = V, aes(x, y, size = pmax(d, 0.8)),
             shape = 21, fill = ifelse(V$d > 0, "#4D4D4D", "#DDDDDD"),
             colour = "white", stroke = 0.25) +
  scale_size(range = c(1.0, 3.2), guide = "none") +
  geom_text(data = filter(V, x == 0), aes(x - 0.045, y, label = g),
            hjust = 1, size = pt(FS_LABEL)) +
  geom_text(data = filter(V, x == 1), aes(x + 0.045, y, label = g),
            hjust = 0, size = pt(FS_LABEL)) +
  annotate("text", x = 0, y = 0.9, label = "Tumor–intrinsic",
           hjust = 0.5, size = pt(FS_AXIS_TITLE)) +
  annotate("text", x = 1, y = 0.9, label = "Niche–associated",
           hjust = 0.5, size = pt(FS_AXIS_TITLE)) +
  coord_cartesian(xlim = c(-0.42, 1.42), ylim = c(min(V$y) - 0.8, 1.6), clip = "off") +
  theme_void(base_size = MIN_FONT_PT) +
  theme(legend.title  = element_text(size = FS_LEGEND),
        legend.text   = element_text(size = FS_LEGEND),
        legend.key.size = unit(5, "pt"))

save_panel(p, supp_fig("SuppFig2O_bipartite_GI_network.pdf"))
log_session("supp_08_bipartite_gi")
message(sprintf("[SuppFig2O] bipartite GI network: %d edges | left %d, right %d, unconnected %d",
                nrow(e), nrow(L), nrow(R), length(un)))
