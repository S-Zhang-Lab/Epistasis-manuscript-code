suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "theme_pub.R"))
  source(here::here("R", "palettes.R"))
  source(here::here("R", "helpers.R"))
  for (pkg in c("ggplot2", "dplyr", "tidyr", "readr", "tibble")) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
})

load_counts <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  if (!"sgRNA" %in% names(df)) stop("Expected 'sgRNA' column in ", path)
  num_cols <- setdiff(names(df), c("sgRNA", "gene"))
  df |>
    dplyr::select(sgRNA, dplyr::all_of(num_cols)) |>
    dplyr::group_by(sgRNA) |>
    dplyr::summarize(dplyr::across(dplyr::everything(), sum), .groups = "drop") |>
    tibble::column_to_rownames("sgRNA") |>
    as.matrix()
}

log_cpm <- function(m) {
  cpm <- sweep(m, 2, colSums(m), "/") * 1e6
  log2(cpm + 1)
}

cx_counts  <- load_counts(raw("Cx3cr1_raw_counts.csv"))
mrc_counts <- load_counts(raw("Mrc1_raw_counts.csv"))

shared <- intersect(rownames(cx_counts), rownames(mrc_counts))
if (length(shared) == 0L) stop("No shared sgRNAs between the two screens.")

cx_l  <- log_cpm(cx_counts [shared, , drop = FALSE])
mrc_l <- log_cpm(mrc_counts[shared, , drop = FALSE])
colnames(cx_l)  <- paste0("Cx3cr1_", colnames(cx_l))
colnames(mrc_l) <- paste0("Mrc1_",   colnames(mrc_l))

run_pca <- function(combined) {
  X  <- t(combined)
  X  <- X[, !apply(X, 2, function(v) any(!is.finite(v))), drop = FALSE]
  Xc <- scale(X, center = TRUE, scale = FALSE)
  sd_col <- apply(Xc, 2, stats::sd)
  Xs <- sweep(Xc, 2, sd_col + 1e-12, "/")
  s  <- svd(Xs)
  coords <- sweep(s$u[, 1:2, drop = FALSE], 2, s$d[1:2], "*")
  rownames(coords) <- rownames(X)

  is_cx_dep <- grepl("Cx3cr1_Depleted", rownames(coords))
  is_cx_all <- grepl("Cx3cr1",          rownames(coords))
  if (mean(coords[is_cx_dep, 2]) < 0) coords[, 2] <- -coords[, 2]
  if (mean(coords[is_cx_all, 1]) < 0) coords[, 1] <- -coords[, 1]

  explained <- (s$d^2) / sum(s$d^2) * 100
  list(coords = coords, pc1v = explained[1], pc2v = explained[2])
}

cond_colors  <- c(Cell = "#E15759", WT = "#4E79A7", Depleted = "#59A14F")
shape_values <- c(Cx3cr1 = 21, Mrc1 = 22, Cell = 24)

variants <- list(
  list(label = "nocell", file = "Fig5S_S3_pca.pdf",    keep_cell = FALSE)
)

for (v in variants) {
  combined <- cbind(cx_l, mrc_l)
  if (!v$keep_cell) combined <- combined[, !grepl("_Cell_", colnames(combined)), drop = FALSE]
  pca  <- run_pca(combined)

  meta <- as.data.frame(pca$coords)
  names(meta) <- c("PC1", "PC2")
  meta$sample <- rownames(meta)
  meta$exp    <- sub("^([^_]+)_.*$", "\\1", meta$sample)
  meta$cond   <- sub("^[^_]+_([^_]+)_[0-9]+$", "\\1", meta$sample)
  meta$cond   <- factor(meta$cond, levels = names(cond_colors))
  meta$shape_class <- ifelse(meta$cond == "Cell", "Cell",
                             ifelse(meta$exp == "Cx3cr1", "Cx3cr1", "Mrc1"))
  meta$shape_class <- factor(meta$shape_class, levels = names(shape_values))

  p <- ggplot(meta, aes(x = PC1, y = PC2, fill = cond, shape = shape_class)) +
    geom_hline(yintercept = 0, color = "gray85", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "gray85", linewidth = 0.3) +
    geom_point(size = 2.6, stroke = 0.4, color = "white", alpha = 0.95) +
    scale_fill_manual(values = cond_colors,  name = "Condition",
                      guide = guide_legend(override.aes = list(shape = 21))) +
    scale_shape_manual(values = shape_values, name = "Sample") +
    labs(title = "PCA of Cas12a dual-gene library screen",
         subtitle = "Cx3cr1 and Mrc1 depletion experiments",
         x = sprintf("PC1 (%.1f%% variance)", pca$pc1v),
         y = sprintf("PC2 (%.1f%% variance)", pca$pc2v)) +
    theme_pub() +
    theme(legend.position = "right")

  save_pdf(p, fig(v$file), width_mm = 121, height_mm = 95)
  message(sprintf("[10_panel_A_pca] %s — n=%d samples, PC1=%.1f%%, PC2=%.1f%%",
                  v$file, nrow(meta), pca$pc1v, pca$pc2v))
}

log_session("10_panel_A_pca")
