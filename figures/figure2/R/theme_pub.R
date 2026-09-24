suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gtable)
})

PANEL_PT     <- 120
MIN_FONT_PT  <- 5

PANEL_FAMILY <- "Helvetica"

PAIR_SEP <- "_"

pair_label <- function(x) {
  vapply(strsplit(as.character(x), "*", fixed = TRUE), function(parts) {
    parts <- parts[nzchar(parts)]
    parts <- paste0(toupper(substring(parts, 1, 1)),
                    tolower(substring(parts, 2)))
    paste(parts, collapse = PAIR_SEP)
  }, character(1))
}
PDF_ENCODING <- "WinAnsi"

FS_AXIS_TEXT  <- 5
FS_AXIS_TITLE <- 6.5
FS_TITLE      <- 7
FS_SUBTITLE   <- 5.5
FS_LEGEND     <- 5
FS_LABEL      <- 5
FS_STRIP      <- 5

SZ_POINT      <- 1.3
SZ_POINT_SM   <- 0.30
SZ_POINT_HL   <- 0.55
LW_THIN       <- 0.25
LW_RULE       <- 0.4

pt <- function(x) x / .pt

theme_panel <- function(base_size = MIN_FONT_PT, grid_y = TRUE,
                        base_family = PANEL_FAMILY) {
  th <- theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text              = element_text(size = base_size, family = base_family,
                                       face = "plain"),
      plot.title        = element_text(face = "plain", size = FS_TITLE,
                                       hjust = 0, margin = margin(b = 1.5)),
      plot.subtitle     = element_text(size = FS_SUBTITLE, colour = "#555555",
                                       hjust = 0, margin = margin(b = 3)),

      plot.caption      = element_text(size = FS_SUBTITLE, colour = "#555555",
                                       hjust = 0, margin = margin(t = 3)),
      axis.title.x      = element_text(size = FS_AXIS_TITLE, margin = margin(t = 2)),
      axis.title.y      = element_text(size = FS_AXIS_TITLE, margin = margin(r = 2)),
      axis.text         = element_text(size = FS_AXIS_TEXT, colour = "#222222"),
      axis.line         = element_line(colour = "#333333", linewidth = LW_THIN),
      axis.ticks        = element_line(colour = "#333333", linewidth = LW_THIN),
      axis.ticks.length = unit(1.5, "pt"),
      legend.title      = element_blank(),
      legend.text       = element_text(size = FS_LEGEND),
      legend.key.size   = unit(6, "pt"),
      legend.margin     = margin(0, 0, 0, 2),
      legend.background = element_blank(),
      legend.spacing    = unit(2, "pt"),
      strip.text        = element_text(size = FS_STRIP),
      strip.background  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin       = margin(3, 3, 3, 3)
    )
  if (grid_y) {
    th <- th + theme(panel.grid.major.y = element_line(colour = "#eeeeee",
                                                       linewidth = 0.2))
  } else {
    th <- th + theme(panel.grid.major.y = element_blank())
  }
  th
}

.fix_panel <- function(plot, panel_w = PANEL_PT, panel_h = PANEL_PT) {

  .dev_open <- FALSE
  try({
    grDevices::pdf(nullfile(), width = 8, height = 8, family = PANEL_FAMILY,
                   encoding = PDF_ENCODING, useDingbats = FALSE)
    .dev_open <- TRUE
  }, silent = TRUE)
  on.exit(if (.dev_open) grDevices::dev.off(), add = TRUE)

  g <- ggplotGrob(plot)
  lay <- g$layout[grepl("^panel", g$layout$name), , drop = FALSE]
  if (nrow(lay) == 0) stop("no panel found in this plot's gtable")
  g$widths[unique(lay$l)]  <- unit(panel_w, "pt")
  g$heights[unique(lay$t)] <- unit(panel_h, "pt")

  zap <- function(u) {
    v <- vapply(seq_along(u), function(i) {
      ui <- u[i]
      if (identical(as.character(unitType(ui)), "null")) 0
      else convertWidth(ui, "pt", valueOnly = TRUE)
    }, numeric(1))
    sum(v)
  }
  w_pt <- zap(g$widths)
  h_pt <- zap(g$heights)

  deficit <- 0
  txt_rows <- character(0)
  for (nm in c("title", "subtitle", "caption")) {
    i <- which(g$layout$name == nm)
    if (!length(i)) next
    gr <- g$grobs[[i]]
    if (inherits(gr, "zeroGrob")) next
    wpt <- try(convertWidth(grobWidth(gr), "pt", valueOnly = TRUE), silent = TRUE)
    if (inherits(wpt, "try-error") || !is.finite(wpt) || wpt <= 0) next
    from <- g$layout$l[i]
    avail <- sum(vapply(seq(from, length(g$widths)), function(j)
      convertWidth(g$widths[j], "pt", valueOnly = TRUE), numeric(1)))
    deficit <- max(deficit, wpt - avail)
    txt_rows <- c(txt_rows, nm)
  }
  if (deficit > 0) {
    g <- gtable::gtable_add_cols(g, unit(deficit + 6, "pt"), pos = -1)
    last <- ncol(g)
    for (nm in txt_rows) {
      k <- which(g$layout$name == nm)
      if (length(k)) g$layout$r[k] <- last
    }
    w_pt <- w_pt + deficit + 6
  }

  iy <- which(g$layout$name == "ylab-l")
  if (length(iy)) {
    gy <- g$grobs[[iy]]
    if (!inherits(gy, "zeroGrob")) {
      hy <- try(convertHeight(grobHeight(gy), "pt", valueOnly = TRUE), silent = TRUE)
      if (!inherits(hy, "try-error") && is.finite(hy) && hy > 0) {
        from <- g$layout$t[iy]
        avail <- sum(vapply(seq(from, length(g$heights)), function(j)
          convertHeight(g$heights[j], "pt", valueOnly = TRUE), numeric(1)))
        if (hy > avail) {
          g <- gtable::gtable_add_rows(g, unit(hy - avail + 6, "pt"), pos = -1)
          g$layout$b[iy] <- nrow(g)
          h_pt <- h_pt + (hy - avail + 6)
        }
      }
    }
  }

  list(grob = g, width_in = w_pt / 72.27, height_in = h_pt / 72.27)
}

save_panel <- function(plot, path, panel_w = PANEL_PT, panel_h = PANEL_PT,
                       png_dpi = 600) {
  fx <- .fix_panel(plot, panel_w, panel_h)
  stem <- sub("\\.pdf$", "", path)
  grDevices::pdf(paste0(stem, ".pdf"), width = fx$width_in,
                 height = fx$height_in, family = PANEL_FAMILY,
                 encoding = PDF_ENCODING, useDingbats = FALSE)
  grid.newpage(); grid.draw(fx$grob); grDevices::dev.off()

  grDevices::png(paste0(stem, ".png"), width = fx$width_in, height = fx$height_in,
                 units = "in", res = png_dpi, type = "cairo",
                 family = PANEL_FAMILY)
  grid.newpage(); grid.draw(fx$grob); grDevices::dev.off()
  .record_panel_size(stem, panel_w, panel_h, fx$width_in, fx$height_in, 1L)
  invisible(c(pdf = paste0(stem, ".pdf"), png = paste0(stem, ".png")))
}

.record_panel_size <- function(stem, panel_w, panel_h, w_in, h_in, n_panels) {
  message(sprintf("    [size] %-38s %s%gx%g pt | page %.1fx%.1f pt",
                  basename(paste0(stem, ".pdf")),
                  if (n_panels > 1) sprintf("%d panels ", n_panels) else "panel ",
                  panel_w, panel_h, w_in * 72.27, h_in * 72.27))
  f <- file.path(OUT_DIR, "panel_sizes.csv")
  row <- data.frame(panel = basename(stem), n_panels = n_panels,
                    panel_w_pt = panel_w, panel_h_pt = panel_h,
                    page_w_pt = round(w_in * 72.27, 2),
                    page_h_pt = round(h_in * 72.27, 2))
  old <- if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) else NULL
  if (!is.null(old)) old <- old[old$panel != row$panel, , drop = FALSE]
  utils::write.csv(rbind(old, row), f, row.names = FALSE, quote = FALSE)
  invisible(row)
}

panel_size_pt <- function(plot, panel_w = PANEL_PT, panel_h = PANEL_PT) {
  g <- .fix_panel(plot, panel_w, panel_h)$grob
  lay <- g$layout[grepl("^panel", g$layout$name), , drop = FALSE]
  c(width  = convertWidth(g$widths[unique(lay$l)][1], "pt", valueOnly = TRUE),
    height = convertHeight(g$heights[unique(lay$t)][1], "pt", valueOnly = TRUE))
}

save_panel_grid <- function(plots, path, ncol = length(plots),
                            panel_w = PANEL_PT, panel_h = PANEL_PT,
                            png_dpi = 600) {
  if (!requireNamespace("gridExtra", quietly = TRUE))
    stop("gridExtra is required to compose multi-panel figures")
  fx <- lapply(plots, .fix_panel, panel_w = panel_w, panel_h = panel_h)
  grobs <- lapply(fx, `[[`, "grob")
  nrow_ <- ceiling(length(plots) / ncol)
  ws <- vapply(fx, `[[`, numeric(1), "width_in")
  hs <- vapply(fx, `[[`, numeric(1), "height_in")

  idx <- seq_along(plots)
  w_in <- sum(vapply(split(ws, (idx - 1) %% ncol), max, numeric(1)))
  h_in <- sum(vapply(split(hs, (idx - 1) %/% ncol), max, numeric(1)))
  arranged <- gridExtra::arrangeGrob(grobs = grobs, ncol = ncol)
  stem <- sub("\\.pdf$", "", path)
  grDevices::pdf(paste0(stem, ".pdf"), width = w_in, height = h_in,
                 family = PANEL_FAMILY, encoding = PDF_ENCODING,
                 useDingbats = FALSE)
  grid.newpage(); grid.draw(arranged); grDevices::dev.off()
  grDevices::png(paste0(stem, ".png"), width = w_in, height = h_in,
                 units = "in", res = png_dpi, type = "cairo",
                 family = PANEL_FAMILY)
  grid.newpage(); grid.draw(arranged); grDevices::dev.off()
  .record_panel_size(stem, panel_w, panel_h, w_in, h_in, length(plots))
  invisible(paste0(stem, ".pdf"))
}
