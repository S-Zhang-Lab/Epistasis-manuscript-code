suppressPackageStartupMessages({
  library(ggplot2); library(ggh4x); library(grid)
})
FS   <- 5
.MM  <- 2.845276
PANEL_PT <- 100

fig_theme <- theme_classic(base_size = FS, base_family = "Helvetica") +
  theme(text        = element_text(family = "Helvetica", colour = "black"),
        axis.text   = element_text(size = FS, colour = "black"),
        axis.title  = element_text(size = FS),
        axis.line   = element_line(linewidth = 0.3, colour = "black"),
        axis.ticks  = element_line(linewidth = 0.3, colour = "black"),
        plot.title  = element_text(size = FS, face = "bold", hjust = 0, margin = margin(b = 2)),
        legend.text = element_text(size = FS), legend.title = element_text(size = FS),
        legend.key.size = unit(6, "pt"), legend.margin = margin(0,0,0,0),
        legend.background = element_blank(),
        plot.margin = margin(3,3,3,3))

pair_lab <- function(x) gsub("*", "_", as.character(x), fixed = TRUE)

sq_panel <- function(pt = PANEL_PT) ggh4x::force_panelsizes(rows = unit(pt,"pt"), cols = unit(pt,"pt"))

pts <- function(colvar, size = 0.28)
  geom_point(aes(colour = .data[[colvar]]), shape = 16, size = size, stroke = 0)

umap_scaffold <- function(df)
  ggplot(df, aes(UMAP1, UMAP2)) + fig_theme + sq_panel() +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "UMAP1", y = "UMAP2")

save_panel <- function(p, stem, w, h, dir){
  ggsave(file.path(dir, paste0(stem, ".pdf")), p, device = "pdf", width = w, height = h,
         units = "in", useDingbats = FALSE)
  ggsave(file.path(dir, paste0(stem, ".png")), p, width = w, height = h, units = "in", dpi = 300, bg = "white")
}
