pkgs <- c("ggplot2", "dplyr", "igraph", "ggrepel")
inst <- rownames(installed.packages())
if (length(setdiff(pkgs, inst)))
  install.packages(setdiff(pkgs, inst), repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(igraph)
})

source(here::here("R", "paths.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "theme_pub.R"))
source(here::here("R", "palettes.R"))

HUB <- "PTEN"
FR_SEED <- 42L
e   <- read.csv(tbl("Fig2I_strong_GI_edges.csv"), stringsAsFactors = FALSE)
nd  <- read.csv(tbl("Fig2I_node_strong_GI_partner_counts.csv"),
                stringsAsFactors = FALSE)
tc  <- function(x) paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))

all_genes <- toupper(nd$gene)
stopifnot(length(all_genes) == 34L)

mod <- setNames(nd$module, toupper(nd$gene))

g <- graph_from_data_frame(
  e[, c("gene_A", "gene_B")], directed = FALSE,
  vertices = data.frame(name = all_genes)
)
set.seed(FR_SEED)
lay <- layout_with_fr(g, weights = e$abs_GI_Z)
rownames(lay) <- V(g)$name

deg_vec <- degree(g)
V_ <- tibble(
  name   = V(g)$name,
  x      = lay[, 1],
  y      = lay[, 2],
  deg    = deg_vec[name],
  lab    = tc(name),
  module = droplevels(factor(mod[name], levels = names(R2_MODULE_COLORS)))
)

ed <- e |>
  left_join(V_ |> select(name, x, y) |> rename(xa = x, ya = y),
            by = c("gene_A" = "name")) |>
  left_join(V_ |> select(name, x, y) |> rename(xb = x, yb = y),
            by = c("gene_B" = "name")) |>
  mutate(hub = gene_A == HUB | gene_B == HUB)

n_connected <- sum(V_$deg > 0)
n_isolated  <- sum(V_$deg == 0)

p <- ggplot() +
  geom_segment(data = filter(ed, !hub),
               aes(xa, ya, xend = xb, yend = yb,
                   colour = GI_Zscore, linewidth = abs_GI_Z),
               alpha = 0.45) +
  geom_segment(data = filter(ed, hub),
               aes(xa, ya, xend = xb, yend = yb,
                   colour = GI_Zscore, linewidth = abs_GI_Z),
               alpha = 0.90) +
  scale_colour_gradient2(low = GI_BLUE, mid = "grey92", high = GI_RED,
                         midpoint = 0, name = "GI Z-score") +
  scale_linewidth(range = c(0.22, 2.4), guide = "none") +
  geom_point(data = V_, aes(x, y, fill = module, size = pmax(deg, 0.5)),
             shape = 21, colour = "white", stroke = 0.3) +
  geom_point(data = filter(V_, name == HUB), aes(x, y),
             shape = 21, fill = NA, colour = "black",
             size = 5.4, stroke = 0.5) +
  scale_fill_manual(values = R2_MODULE_COLORS, name = "Functional\nmodule") +
  scale_size(range = c(1.0, 4.4), guide = "none") +
  ggrepel::geom_text_repel(
    data = V_,
    aes(x, y, label = lab),
    size = pt(FS_LABEL), seed = FR_SEED, max.overlaps = Inf,
    box.padding = 0.18, point.padding = 0.12,
    segment.size = 0.15, segment.alpha = 0.5,
    min.segment.length = 0.25
  ) +
  coord_equal(clip = "off") +
  labs(title = "GI network among 34 round-2 genes",
       subtitle = sprintf(
         "%d strong interactions (|GI Z| > 1.5); %d connected, %d isolated; FR layout",
         nrow(e), n_connected, n_isolated)) +
  theme_void(base_size = MIN_FONT_PT) +
  theme(plot.title    = element_text(size = FS_TITLE, hjust = 0),
        plot.subtitle = element_text(size = FS_SUBTITLE, colour = "grey40",
                                     hjust = 0),
        legend.title  = element_text(size = FS_LEGEND),
        legend.text   = element_text(size = FS_LEGEND),
        legend.key.size = unit(5, "pt"))

save_panel(p, fig("Fig2I_GI_network.pdf"))
log_session("06_fig2I_pten_ego")
message(sprintf(
  "[Fig2I] FR layout: %d genes (%d connected, %d isolated) | %d edges | seed %d",
  nrow(V_), n_connected, n_isolated, nrow(e), FR_SEED))
