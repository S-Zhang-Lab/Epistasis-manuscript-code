suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  library(ggplot2)
})

SCIENCE_WIDTH_MM <- c(single = 57, one_and_half = 121, double = 183)

PUB_FONT_FAMILY  <- "Helvetica"
PUB_FONT_BODY    <- 7
PUB_FONT_AXIS    <- 7
PUB_FONT_TITLE   <- 8
PUB_FONT_PANEL   <- 9
PUB_FONT_LEGEND  <- 6
PUB_FONT_MIN     <- 5

theme_pub <- function(base_size = PUB_FONT_BODY, base_family = PUB_FONT_FAMILY) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line       = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks      = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(2, "pt"),
      axis.text       = element_text(size = PUB_FONT_AXIS, colour = "black"),
      axis.title      = element_text(size = PUB_FONT_TITLE, colour = "black"),
      plot.title      = element_text(size = PUB_FONT_PANEL, face = "plain", hjust = 0),
      plot.subtitle   = element_text(size = PUB_FONT_BODY, colour = "grey25"),
      legend.text     = element_text(size = PUB_FONT_LEGEND),
      legend.title    = element_text(size = PUB_FONT_LEGEND),
      legend.key.size = unit(8, "pt"),
      legend.background = element_blank(),
      panel.grid      = element_blank(),
      strip.background  = element_blank(),
      strip.text      = element_text(size = PUB_FONT_TITLE),
      plot.margin     = margin(2, 2, 2, 2, "mm")
    )
}

save_pdf <- function(plot, file, width_mm, height_mm, snap_width = TRUE) {
  stopifnot(grepl("\\.pdf$", file, ignore.case = TRUE))

  if (snap_width && !width_mm %in% SCIENCE_WIDTH_MM) {
    snapped <- SCIENCE_WIDTH_MM[which.min(abs(SCIENCE_WIDTH_MM - width_mm))]
    warning(sprintf(
      "save_pdf(): width %.1f mm is not a Science column width; snapping to %.0f mm (%s).",
      width_mm, snapped, names(snapped)
    ))
    width_mm <- as.numeric(snapped)
  }

  is_heatmap <- inherits(plot, c("Heatmap", "HeatmapList"))

  grDevices::pdf(
    file        = file,
    width       = width_mm  / 25.4,
    height      = height_mm / 25.4,
    family      = PUB_FONT_FAMILY,
    bg          = if (is_heatmap) "white" else "transparent",
    useDingbats = FALSE
  )
  on.exit(grDevices::dev.off(), add = TRUE)

  if (inherits(plot, "ggplot")) {
    print(plot)
  } else if (is_heatmap) {
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      stop("ComplexHeatmap is required to save Heatmap objects.")
    }
    ComplexHeatmap::draw(plot, newpage = FALSE)
  } else {
    grid::grid.draw(plot)
  }

  invisible(file)
}
