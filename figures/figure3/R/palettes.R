build_symmetric_scale <- function(limits = c(-2, 2),
                                  low    = "#6FA3D2",
                                  mid    = "#FFFFFF",
                                  high   = "#B2182B",
                                  n      = 121L) {
  stopifnot(length(limits) == 2, limits[1] < 0, limits[2] > 0)
  stopifnot(abs(limits[1] + limits[2]) < 1e-9)
  if (n %% 2 == 0) n <- n + 1L
  half <- (n - 1L) %/% 2L
  breaks <- c(seq(limits[1], 0,         length.out = half + 1L),
              seq(0,         limits[2], length.out = half + 1L)[-1])
  colors <- c(colorRampPalette(c(low,  mid))(half + 1L),
              colorRampPalette(c(mid,  high))(half + 1L)[-1])
  list(breaks = breaks, colors = colors)
}

build_asymmetric_scale <- function(limits = c(-3.5, 1.4),
                                   low    = "#6FA3D2",
                                   mid    = "#FFFFFF",
                                   high   = "#B2182B",
                                   n      = 121L) {
  stopifnot(length(limits) == 2, limits[1] < 0, limits[2] > 0)
  range_total <- limits[2] - limits[1]
  n_neg <- round(n * (-limits[1]) / range_total)
  n_pos <- n - n_neg
  breaks <- c(seq(limits[1], 0,         length.out = n_neg + 1L),
              seq(0,         limits[2], length.out = n_pos + 1L)[-1])
  colors <- c(colorRampPalette(c(low,  mid))(n_neg + 1L),
              colorRampPalette(c(mid,  high))(n_pos + 1L)[-1])
  list(breaks = breaks, colors = colors)
}

to_colorRamp2 <- function(scale) {
  circlize::colorRamp2(scale$breaks, scale$colors)
}

if (sys.nframe() == 0L) {
  s <- build_asymmetric_scale()
  zero_idx <- which.min(abs(s$breaks))
  stopifnot(s$colors[zero_idx] == "#FFFFFF")
  message("utils_colorscale.R: asymmetric scale anchors zero at white.")
}
