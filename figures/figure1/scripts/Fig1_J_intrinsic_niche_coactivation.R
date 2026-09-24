message("[Panel J]  DawnRank intrinsic \u00d7 niche correlation matrix — building...")

DAWNRANK_FILE <- file.path(INPUT_DIR, "GSE110590", "JCI96153.sdt1-8.xlsx")
DAWNRANK_SHEET <- "Supplementary Table 5"

LUNG_MET_COLS <- c(
  "A2.LungMet",
  "A11.LUNG.MET",
  "A20.LungMet.L",
  "A20.LungMet.R",
  "A26.Lung.Met.1",
  "A26.Lung.Met.2",
  "A28.LungMet.LLL",
  "A28.LungMet.RUL"
)

R_ANNOT_FLOOR <- 0.55

message("  [1/4] Loading 35-gene panel from snapshot...")

gl <- readLines(file.path(GENE_LISTS_DIR, "panel_G_dawnrank_35.txt"))

intrinsic_idx <- grep("^## Intrinsic", gl)[1]
niche_idx     <- grep("^## Niche", gl)[1]
stopifnot(!is.na(intrinsic_idx) & !is.na(niche_idx))

extract_block <- function(start, end_excl) {
  block <- gl[(start + 1):(end_excl - 1)]
  block <- trimws(block)
  block <- block[nzchar(block) & !startsWith(block, "#")]
  block
}

intrinsic_genes <- extract_block(intrinsic_idx, niche_idx)
niche_genes     <- extract_block(niche_idx, length(gl) + 1)

panel_genes <- c(intrinsic_genes, niche_genes)
message("    intrinsic (", length(intrinsic_genes), "): ",
        paste(intrinsic_genes, collapse = ", "))
message("    niche-associated (", length(niche_genes), "): ",
        paste(niche_genes, collapse = ", "))
message("    total panel: ", length(panel_genes))

message("  [2/4] Loading DawnRank Supp Table 5 + subsetting...")

dawnrank_raw <- read_excel(DAWNRANK_FILE, sheet = DAWNRANK_SHEET, skip = 1)

names(dawnrank_raw)[1] <- "Gene"

dawnrank_raw$Gene <- sub("\\|.*$", "", dawnrank_raw$Gene)

missing_cols <- setdiff(LUNG_MET_COLS, colnames(dawnrank_raw))
if (length(missing_cols) > 0) {
  stop("Missing lung-met columns: ", paste(missing_cols, collapse = ", "))
}

dawnrank_sub <- dawnrank_raw %>%
  filter(Gene %in% panel_genes) %>%
  select(Gene, all_of(LUNG_MET_COLS))

found    <- intersect(panel_genes, dawnrank_sub$Gene)
not_found <- setdiff(panel_genes, dawnrank_sub$Gene)
if (length(not_found) > 0) {
  message("    NOTE: ", length(not_found),
          " panel genes not in DawnRank Supp Table 5: ",
          paste(not_found, collapse = ", "))
}

dawnrank_sub <- dawnrank_sub[match(found, dawnrank_sub$Gene), ]
mat_G <- as.matrix(dawnrank_sub[, LUNG_MET_COLS])
rownames(mat_G) <- dawnrank_sub$Gene
storage.mode(mat_G) <- "numeric"

message("  [3/4] Pearson correlation (n = ", ncol(mat_G), " lung mets)...")
cor_G <- cor(t(mat_G), method = "pearson")

write.csv(cor_G,
          file.path(OUTPUT_TBL_DIR, "Fig1G_DawnRank_correlation_matrix.csv"))

gene_groups <- ifelse(rownames(cor_G) %in% intrinsic_genes,
                      "Intrinsic", "Niche-associated")
gene_groups <- factor(gene_groups, levels = c("Intrinsic", "Niche-associated"))

message("  [4/4] Rendering correlation heatmap...")

cor_palette <- colorRamp2(c(-1, 0, 1),
                          c("#2166AC", "white", "#B2182B"))

top_anno <- HeatmapAnnotation(
  Type    = as.character(gene_groups),
  col     = list(Type = LIBRARY_COLORS),
  annotation_name_side  = "left",
  annotation_name_gp    = gpar(fontsize = 8),
  simple_anno_size      = unit(3.5, "mm"),
  show_legend           = TRUE,
  annotation_legend_param = list(
    Type = list(title          = "Driver class",
                title_gp       = gpar(fontsize = 9, fontface = "bold"),
                labels_gp      = gpar(fontsize = 9))
  )
)

pdf_device(file.path(OUTPUT_FIG_DIR, "Fig1_J_intrinsic_niche_coactivation.pdf"),
          width = 10, height = 9)

ht_G <- Heatmap(
  cor_G,
  name              = "Pearson r",
  col               = cor_palette,
  row_split         = gene_groups,
  column_split      = gene_groups,
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  cluster_row_slices    = FALSE,
  cluster_column_slices = FALSE,
  row_title         = NULL,
  column_title      = NULL,
  row_names_gp      = gpar(fontsize = 8),
  column_names_gp   = gpar(fontsize = 8),
  row_names_side    = "right",
  column_names_side = "bottom",
  rect_gp           = gpar(type = "none"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (i <= j) {
      grid.rect(x, y, width, height,
                gp = gpar(fill = fill, col = "gray90"))
      if (i != j && abs(cor_G[i, j]) >= R_ANNOT_FLOOR) {
        grid.text(sprintf("%.2f", cor_G[i, j]), x, y,
                  gp = gpar(fontsize = 6, fontface = "bold"))
      }
    }
  },
  top_annotation    = top_anno,
  heatmap_legend_param = list(
    title       = "Pearson r",
    title_gp    = gpar(fontsize = 9, fontface = "bold"),
    at          = c(-1, -0.55, 0, 0.55, 1),
    labels_gp   = gpar(fontsize = 9),
    legend_height = unit(3.5, "cm")
  ),
  show_heatmap_legend = TRUE
)

draw(
  ht_G,
  column_title    = "Intrinsic \u00d7 niche driver co-activation in clinical lung metastases (n = 8 lung mets)",
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  padding         = unit(c(12, 4, 4, 4), "mm")
)

grid.text(
  "Gene-level co-activation; the pathway-level view is Panel G",
  x = unit(0.5, "npc"),
  y = unit(4, "mm"),
  just = "centre",
  gp = gpar(fontsize = 10, fontface = "italic", col = "grey25")
)

dev.off()

message("[Panel J]  Done. \u2192 ",
        file.path(OUTPUT_FIG_DIR, "Fig1_J_intrinsic_niche_coactivation.pdf"))
