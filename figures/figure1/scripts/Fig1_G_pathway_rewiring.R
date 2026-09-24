message("[Panel H]  Hallmark module co-activation rewiring — building...")

MODULE_ORDER <- c(

  "MYC_TARGETS_V1", "MYC_TARGETS_V2", "E2F_TARGETS", "G2M_CHECKPOINT",
  "MITOTIC_SPINDLE", "DNA_REPAIR", "P53_PATHWAY", "APOPTOSIS",

  "OXIDATIVE_PHOSPHORYLATION", "GLYCOLYSIS", "FATTY_ACID_METABOLISM",
  "ADIPOGENESIS", "BILE_ACID_METABOLISM", "CHOLESTEROL_HOMEOSTASIS",
  "XENOBIOTIC_METABOLISM", "PEROXISOME", "HEME_METABOLISM",

  "HYPOXIA", "UNFOLDED_PROTEIN_RESPONSE", "REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "UV_RESPONSE_UP", "UV_RESPONSE_DN", "PROTEIN_SECRETION",

  "PI3K_AKT_MTOR_SIGNALING", "MTORC1_SIGNALING", "KRAS_SIGNALING_UP",
  "KRAS_SIGNALING_DN", "TGF_BETA_SIGNALING", "NOTCH_SIGNALING",
  "HEDGEHOG_SIGNALING", "WNT_BETA_CATENIN_SIGNALING", "ANGIOGENESIS",

  "ANDROGEN_RESPONSE", "ESTROGEN_RESPONSE_EARLY", "ESTROGEN_RESPONSE_LATE",
  "MYOGENESIS", "SPERMATOGENESIS", "PANCREAS_BETA_CELLS",

  "IL6_JAK_STAT3_SIGNALING", "IL2_STAT5_SIGNALING", "INFLAMMATORY_RESPONSE",
  "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE",
  "TNFA_SIGNALING_VIA_NFKB", "ALLOGRAFT_REJECTION", "COMPLEMENT", "COAGULATION",

  "EPITHELIAL_MESENCHYMAL_TRANSITION", "APICAL_JUNCTION", "APICAL_SURFACE"
)

col_corr  <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))
col_delta <- colorRamp2(c(-1, 0, 1), c("#762A83", "white", "#1B7837"))

SUP_DELTA_THRESHOLD <- 1.0

REWIRING_CALLOUTS <- list(
  "Estrogen-Early decoupling" = list(
    anchor  = "ESTROGEN_RESPONSE_EARLY",
    members = c("EPITHELIAL_MESENCHYMAL_TRANSITION",
                "WNT_BETA_CATENIN_SIGNALING",
                "IL6_JAK_STAT3_SIGNALING",
                "HYPOXIA"),
    color   = "#1B7837"
  ),
  "MYC \u00d7 Cell-Cycle \u00d7 WNT emergence" = list(
    anchor  = NULL,
    members = c("MYC_TARGETS_V1", "MYC_TARGETS_V2", "E2F_TARGETS",
                "G2M_CHECKPOINT", "WNT_BETA_CATENIN_SIGNALING"),
    color   = "#762A83"
  )
)

FIG2_PAIRS <- list(
  list(intrinsic = "PI3K_AKT_MTOR_SIGNALING",   niche = "INFLAMMATORY_RESPONSE",
       label = "PTEN \u00d7 CX3CL1/TLR7"),
  list(intrinsic = "MTORC1_SIGNALING",          niche = "INFLAMMATORY_RESPONSE",
       label = "TSC2/STK11 \u00d7 CX3CL1/TLR7"),
  list(intrinsic = "TGF_BETA_SIGNALING",        niche = "EPITHELIAL_MESENCHYMAL_TRANSITION",
       label = "TGFBR2/SMAD4 \u00d7 CDH1"),
  list(intrinsic = "WNT_BETA_CATENIN_SIGNALING", niche = "EPITHELIAL_MESENCHYMAL_TRANSITION",
       label = "APC \u00d7 CDH1"),
  list(intrinsic = "PI3K_AKT_MTOR_SIGNALING",   niche = "TNFA_SIGNALING_VIA_NFKB",
       label = "PTEN \u00d7 TNFSF15")
)

message("  [1/6] Per-sample Hallmark module scores (mean-z)...")

hallmark_msig <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_msig$gene_symbol, hallmark_msig$gs_name)

exprs_z <- t(scale(t(exprs_mat)))

hallmark_scores <- vapply(
  hallmark_list,
  function(g) {
    g_in <- intersect(g, rownames(exprs_z))
    if (length(g_in) < 10) return(rep(NA_real_, ncol(exprs_z)))
    colMeans(exprs_z[g_in, , drop = FALSE], na.rm = TRUE)
  },
  numeric(ncol(exprs_z))
)
hallmark_scores <- as.data.frame(hallmark_scores)
rownames(hallmark_scores) <- colnames(exprs_z)
colnames(hallmark_scores) <- sub("HALLMARK_", "", colnames(hallmark_scores))

scores_sub  <- hallmark_scores[custom_meta$sample_id, , drop = FALSE]
primary_ids <- custom_meta$sample_id[group_sub == "Primary"]
met_ids     <- custom_meta$sample_id[group_sub == "LungMet"]

valid_modules <- colnames(scores_sub)[!is.na(colSums(scores_sub))]
message("    measurable modules: ", length(valid_modules), " / 50 Hallmarks")

ordered_modules <- intersect(MODULE_ORDER, valid_modules)
extra_modules   <- setdiff(valid_modules, MODULE_ORDER)
if (length(extra_modules) > 0) {
  message("    appending non-themed modules: ",
          paste(extra_modules, collapse = ", "))
  ordered_modules <- c(ordered_modules, extra_modules)
}

message("  [2/6] Spearman correlation matrices + Δρ...")

cor_primary <- cor(scores_sub[primary_ids, ordered_modules], method = "spearman")
cor_lungmet <- cor(scores_sub[met_ids,     ordered_modules], method = "spearman")
cor_delta   <- cor_lungmet - cor_primary

write.csv(cor_primary,
          file.path(OUTPUT_TBL_DIR, "Hallmark_corr_primary_paired_subset.csv"))
write.csv(cor_lungmet,
          file.path(OUTPUT_TBL_DIR, "Hallmark_corr_lungmet_paired_subset.csv"))
write.csv(cor_delta,
          file.path(OUTPUT_TBL_DIR, "Hallmark_delta_correlation_paired_subset.csv"))

pretty_mod <- function(x) gsub("_", " ", x)

message("  [3/6] Building decorator helpers (callout boxes + Fig 2 glyphs)...")

decorate_callouts <- function(modules) {
  n <- length(modules)
  for (callout in REWIRING_CALLOUTS) {
    members <- intersect(c(callout$anchor, callout$members), modules)
    if (length(members) < 2) next
    col_idx <- which(modules %in% members)
    row_idx <- which(modules %in% members)

    x_lo <- (min(col_idx) - 1) / n
    x_hi <-  max(col_idx)      / n
    y_lo <- 1 - max(row_idx)      / n
    y_hi <- 1 - (min(row_idx) - 1) / n
    grid.rect(x = unit((x_lo + x_hi) / 2, "npc"),
              y = unit((y_lo + y_hi) / 2, "npc"),
              width  = unit(x_hi - x_lo, "npc"),
              height = unit(y_hi - y_lo, "npc"),
              gp = gpar(col = callout$color, lwd = 2.2, fill = NA, lty = "solid"))
  }
}

decorate_fig2_glyphs <- function(modules) {
  n <- length(modules)
  for (pair in FIG2_PAIRS) {
    if (!(pair$intrinsic %in% modules) || !(pair$niche %in% modules)) next
    i <- which(modules == pair$intrinsic)
    j <- which(modules == pair$niche)

    for (cell in list(c(i, j), c(j, i))) {
      r <- cell[1]; c <- cell[2]
      x <- (c - 0.5) / n
      y <- 1 - (r - 0.5) / n
      grid.points(x = unit(x, "npc"), y = unit(y, "npc"),
                  pch = 21, size = unit(2.4, "mm"),
                  gp = gpar(col = "black", fill = "white", lwd = 1.0))
    }
  }
}

DELTA_RHO      <- expression(Delta * rho)
DELTA_RHO_DIFF <- expression(paste(Delta, rho, "  (LungMet - Primary)"))
REWIRE_TITLE   <- expression(paste("Hallmark co-activation rewiring  (",
                                   Delta, rho, " = LungMet - Primary)"))
SPEARMAN_RHO   <- expression(paste("Spearman ", rho))
HM_DELTA_RHO   <- "delta_rho"

message("  [4/6] Rendering Variant A (Δρ only)...")

ht_delta <- Heatmap(
  cor_delta,
  name             = HM_DELTA_RHO,
  col              = col_delta,
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  row_order        = ordered_modules,
  column_order     = ordered_modules,
  row_labels       = pretty_mod(ordered_modules),
  column_labels    = pretty_mod(ordered_modules),
  row_names_gp     = gpar(fontsize = 7),
  column_names_gp  = gpar(fontsize = 7),
  border           = TRUE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title          = DELTA_RHO_DIFF,
    title_gp       = gpar(fontsize = 9, fontface = "bold"),
    legend_height  = unit(3.5, "cm"),
    at             = c(-1, -0.5, 0, 0.5, 1)
  )
)

pdf_device(panel_output_path("H", "variant_A",
                            "Fig1_G_pathway_rewiring.pdf"),
          width = 10, height = 9)
draw(ht_delta,
     column_title    = REWIRE_TITLE,
     column_title_gp = gpar(fontsize = 13, fontface = "bold"),
     padding         = unit(c(2, 4, 4, 4), "mm"))
decorate_heatmap_body(HM_DELTA_RHO, code = {
  decorate_callouts(ordered_modules)
  decorate_fig2_glyphs(ordered_modules)
})

grid.text("\u2022  pathway pair interrogated in Fig 2",
          x = unit(0.99, "npc"), y = unit(2, "mm"),
          just = c("right", "bottom"),
          gp = gpar(fontsize = 8, fontface = "italic", col = "grey25"))
dev.off()

message("  [5/6] Rendering Variant B (Supp. Fig. 1 panel E: two source matrices)...")

make_corr_heatmap <- function(mat, title, palette, name_id, legend_title,
                              show_legend = TRUE) {
  Heatmap(
    mat,
    name             = name_id,
    col              = palette,
    cluster_rows     = FALSE,
    cluster_columns  = FALSE,
    row_order        = ordered_modules,
    column_order     = ordered_modules,
    row_labels       = pretty_mod(ordered_modules),
    column_labels    = pretty_mod(ordered_modules),
    row_names_gp     = gpar(fontsize = 6),
    column_names_gp  = gpar(fontsize = 6),
    border           = TRUE,
    column_title     = title,
    column_title_gp  = gpar(fontsize = 11, fontface = "bold"),
    show_heatmap_legend = show_legend,
    heatmap_legend_param = list(
      title         = legend_title,
      title_gp      = gpar(fontsize = 8, fontface = "bold"),
      legend_height = unit(2.5, "cm")
    )
  )
}

ht_p <- make_corr_heatmap(cor_primary, "Primary Tumor Hallmark Pathway Co-activation",
                          col_corr, "rho_primary", SPEARMAN_RHO, show_legend = TRUE)
ht_l <- make_corr_heatmap(cor_lungmet, "LungMet Hallmark Pathway Co-activation",
                          col_corr, "rho_lungmet", SPEARMAN_RHO, show_legend = FALSE)

pdf_device(panel_output_path("H", "variant_B",
                            "FigS1_G_coactivation_matrices.pdf",
                            placed = TRUE),
          width = 13, height = 6.5)
draw(ht_p + ht_l, ht_gap = unit(6, "mm"))
dev.off()

message("  [6/6] Rendering Variant C (Δρ dominant + small baselines)...")

ht_d_big <- Heatmap(
  cor_delta,
  name             = HM_DELTA_RHO,
  col              = col_delta,
  cluster_rows     = FALSE, cluster_columns = FALSE,
  row_order        = ordered_modules,
  column_order     = ordered_modules,
  row_labels       = pretty_mod(ordered_modules),
  column_labels    = pretty_mod(ordered_modules),
  row_names_gp     = gpar(fontsize = 7),
  column_names_gp  = gpar(fontsize = 7),
  border           = TRUE,
  column_title     = REWIRE_TITLE,
  column_title_gp  = gpar(fontsize = 12, fontface = "bold"),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title         = DELTA_RHO,
    title_gp      = gpar(fontsize = 9, fontface = "bold"),
    legend_height = unit(3, "cm")
  )
)

ht_p_small <- Heatmap(
  cor_primary,
  name             = "rho_primary",
  col              = col_corr,
  cluster_rows     = FALSE, cluster_columns = FALSE,
  row_order        = ordered_modules,
  column_order     = ordered_modules,
  show_row_names   = FALSE, show_column_names = FALSE,
  border           = TRUE,
  column_title     = "Primary tumor",
  column_title_gp  = gpar(fontsize = 9, fontface = "bold"),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(title = expression(rho), legend_height = unit(2, "cm"))
)
ht_l_small <- Heatmap(
  cor_lungmet,
  name             = "rho_lungmet",
  col              = col_corr,
  cluster_rows     = FALSE, cluster_columns = FALSE,
  row_order        = ordered_modules,
  column_order     = ordered_modules,
  show_row_names   = FALSE, show_column_names = FALSE,
  border           = TRUE,
  column_title     = "Lung metastasis",
  column_title_gp  = gpar(fontsize = 9, fontface = "bold"),
  show_heatmap_legend = FALSE
)

pdf_device(panel_output_path("H", "variant_C",
                            "Fig1.H_pathway_rewiring_variant_C_delta_dominant.pdf"),
          width = 12, height = 12)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(2.6, 1), "null"))))

pushViewport(viewport(layout.pos.row = 1))
draw(ht_d_big, newpage = FALSE)
decorate_heatmap_body(HM_DELTA_RHO, code = {
  decorate_callouts(ordered_modules)
  decorate_fig2_glyphs(ordered_modules)
})
grid.text("\u2022  pathway pair interrogated in Fig 2",
          x = unit(0.99, "npc"), y = unit(2, "mm"),
          just = c("right", "bottom"),
          gp = gpar(fontsize = 8, fontface = "italic", col = "grey25"))
popViewport()

pushViewport(viewport(layout.pos.row = 2,
                      layout = grid.layout(1, 2, widths = unit(c(1, 1), "null"))))
pushViewport(viewport(layout.pos.col = 1))
draw(ht_p_small, newpage = FALSE)
popViewport()
pushViewport(viewport(layout.pos.col = 2))
draw(ht_l_small, newpage = FALSE)
popViewport()
popViewport()

popViewport()
dev.off()

message("  Supplementary barplot: top rewiring events...")

delta_long <- as.data.frame(as.table(cor_delta)) %>%
  rename(Module1 = Var1, Module2 = Var2, DeltaRho = Freq) %>%
  filter(as.character(Module1) < as.character(Module2)) %>%
  arrange(DeltaRho)

top_rewire <- delta_long %>%
  filter(abs(DeltaRho) >= SUP_DELTA_THRESHOLD) %>%
  mutate(pair = paste(pretty_mod(Module1), "\u00d7", pretty_mod(Module2)))

write.csv(top_rewire,
          file.path(OUTPUT_TBL_DIR, "Fig1H_top_rewiring_events.csv"),
          row.names = FALSE)

if (nrow(top_rewire) > 0) {
  p_H_supp <- ggplot(top_rewire,
                     aes(x = reorder(pair, DeltaRho), y = DeltaRho,
                         fill = DeltaRho > 0)) +
    geom_col(alpha = 0.85) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE"  = "#1B7837",
                                 "FALSE" = "#762A83"),
                      labels = c("FALSE" = "Lost in metastasis",
                                 "TRUE"  = "Gained in metastasis"),
                      name   = NULL) +
    labs(x = NULL,
         y = expression(paste(Delta, rho, "  (Spearman: LungMet - Primary)")),
         title = bquote("Top Hallmark co-activation rewiring events  (|" *
                        Delta * rho * "|" >= .(sprintf("%.1f", SUP_DELTA_THRESHOLD)) * ")")) +
    theme_classic(base_size = 10) +
    theme(axis.text.y     = element_text(size = 7),
          plot.title      = element_text(face = "bold"),
          legend.position = "top")

  ggsave(file.path(OUTPUT_SUP_VAR_DIR, "Fig1.Sup_Network_delta_rewiring.pdf"),
         p_H_supp, width = 9, height = 9, device = pdf_device)
  message("    supp barplot: ", nrow(top_rewire), " rewiring events at |\u0394\u03c1| \u2265 ",
          SUP_DELTA_THRESHOLD)
}

message("[Panel H]  Done. Three variants + supp barplot rendered.")
