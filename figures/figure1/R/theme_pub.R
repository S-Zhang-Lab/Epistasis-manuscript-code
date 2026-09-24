suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2", repos = "https://cloud.r-project.org")
  }
  library(ggplot2)
})

SCIENCE_WIDTH_MM <- c(single = 57, one_and_half = 121, double = 183)

PUB_FONT_FAMILY <- "Helvetica"
PUB_FONT_MIN    <- 5
PUB_FONT_LEGEND <- 6
PUB_FONT_BODY   <- 7
PUB_FONT_AXIS   <- 7
PUB_FONT_TITLE  <- 8
PUB_FONT_PANEL  <- 9

theme_pub <- function(base_size = PUB_FONT_BODY,
                      base_family = PUB_FONT_FAMILY) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title    = element_text(size = PUB_FONT_TITLE, face = "bold"),
      plot.subtitle = element_text(size = PUB_FONT_BODY,  color = "grey25"),
      plot.caption  = element_text(size = PUB_FONT_MIN,   color = "grey35",
                                   hjust = 0),
      axis.title    = element_text(size = PUB_FONT_AXIS),
      axis.text     = element_text(size = PUB_FONT_AXIS, color = "black"),
      legend.title  = element_text(size = PUB_FONT_LEGEND),
      legend.text   = element_text(size = PUB_FONT_LEGEND),
      legend.key.size  = unit(0.7, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

save_pub_pdf <- function(plot, path, width_mm = 121, height_mm = 90) {
  if (!any(abs(width_mm - SCIENCE_WIDTH_MM) < 1)) {
    warning(sprintf(
      "save_pub_pdf: width %.0f mm is not a Science column width (%s).",
      width_mm, paste(SCIENCE_WIDTH_MM, collapse = " / ")))
  }
  ggplot2::ggsave(
    filename = path, plot = plot,
    width = width_mm, height = height_mm, units = "mm",
    device = pdf_device
  )
  invisible(path)
}
