# =============================================================================
#  Fig1_G_pathway_rewiring.R — Panel G
#  Hallmark module co-activation rewiring (Δρ = ρ_LungMet − ρ_Primary)
# -----------------------------------------------------------------------------
#  Design (locked in plan) — the crescendo of Figure 1.
#    - Per-sample module activity = mean z-score over each Hallmark gene set.
#    - Spearman correlation per group: Primary (n = 5) vs LungMet (n = 8).
#    - Δρ matrix = correlation difference.
#    - GENERATE all three layout variants for side-by-side review:
#        A — single Δρ panel
#        B — three-panel Primary | LungMet | Δρ
#        C — Δρ dominant + small Primary/LungMet insets below
#    - Rewiring-event callout boxes (2 + optional 3rd):
#        * Estrogen-Early axis decoupling
#        * MYC × Cell-Cycle × WNT emergence
#    - Fig 2 library-pair glyphs (●) on Δρ cells corresponding to Hallmark
#      module pairs perturbed by the in vivo Cas12a screen.
#    - Supp barplot of top |Δρ| ≥ 0.5 events.
#    - No on-panel small-n footer; caveats live in figure legend + Methods.
#
#  Requires shared infrastructure from 01_load_GSE110590.R:
#    exprs_mat, custom_meta, group_sub
#
#  Outputs (routed by docs/selected_main_variants.yml; G selects variant_A):
#    output/figures/main/Fig1_G_pathway_rewiring.pdf                variant_A, main panel G
#    output/figures/supplementary/FigS1_G_coactivation_matrices.pdf variant_B, supp panel E
#    output/figures/supplementary/variants/
#      Fig1.H_pathway_rewiring_variant_C_delta_dominant.pdf         not placed
#      Fig1.Sup_Network_delta_rewiring.pdf                          not placed
#    output/tables/Hallmark_corr_primary_paired_subset.csv
#    output/tables/Hallmark_corr_lungmet_paired_subset.csv
#    output/tables/Hallmark_delta_correlation_paired_subset.csv
#    output/tables/Fig1H_top_rewiring_events.csv
#
#  Legacy reference: scripts/legacy/Figure1_publication_monolithic.R  lines 703-818
# =============================================================================

message("[Panel G]  Hallmark module co-activation rewiring — building...")

# -- 0. Configuration --------------------------------------------------------

# Thematic module ordering (groups related Hallmarks adjacent for readability).
MODULE_ORDER <- c(
  # Proliferation / cell cycle / DDR
  "MYC_TARGETS_V1", "MYC_TARGETS_V2", "E2F_TARGETS", "G2M_CHECKPOINT",
  "MITOTIC_SPINDLE", "DNA_REPAIR", "P53_PATHWAY", "APOPTOSIS",
  # Metabolism
  "OXIDATIVE_PHOSPHORYLATION", "GLYCOLYSIS", "FATTY_ACID_METABOLISM",
  "ADIPOGENESIS", "BILE_ACID_METABOLISM", "CHOLESTEROL_HOMEOSTASIS",
  "XENOBIOTIC_METABOLISM", "PEROXISOME", "HEME_METABOLISM",
  # Stress / hypoxia
  "HYPOXIA", "UNFOLDED_PROTEIN_RESPONSE", "REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "UV_RESPONSE_UP", "UV_RESPONSE_DN", "PROTEIN_SECRETION",
  # Signaling
  "PI3K_AKT_MTOR_SIGNALING", "MTORC1_SIGNALING", "KRAS_SIGNALING_UP",
  "KRAS_SIGNALING_DN", "TGF_BETA_SIGNALING", "NOTCH_SIGNALING",
  "HEDGEHOG_SIGNALING", "WNT_BETA_CATENIN_SIGNALING", "ANGIOGENESIS",
  # Hormone / reproductive / developmental
  "ANDROGEN_RESPONSE", "ESTROGEN_RESPONSE_EARLY", "ESTROGEN_RESPONSE_LATE",
  "MYOGENESIS", "SPERMATOGENESIS", "PANCREAS_BETA_CELLS",
  # Immune / inflammation
  "IL6_JAK_STAT3_SIGNALING", "IL2_STAT5_SIGNALING", "INFLAMMATORY_RESPONSE",
  "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE",
  "TNFA_SIGNALING_VIA_NFKB", "ALLOGRAFT_REJECTION", "COMPLEMENT", "COAGULATION",
  # EMT / structure
  "EPITHELIAL_MESENCHYMAL_TRANSITION", "APICAL_JUNCTION", "APICAL_SURFACE"
)

# Color scales.
col_corr  <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))
col_delta <- colorRamp2(c(-1, 0, 1), c("#762A83", "white", "#1B7837"))

# Δρ rewiring threshold for the supplementary barplot.
# Tightened to |Δρ| ≥ 1.0 because the paired n = 5 / n = 8 cohort produces
# noisy correlation differences; |Δρ| ≥ 0.5 yields ~430 pairs (unreadable);
# |Δρ| ≥ 1.0 retains the most striking ~50–80 rewiring events.
SUP_DELTA_THRESHOLD <- 1.0

# Rewiring-event callout boxes (2 — optional third left as TODO).
# Each entry maps to a list of modules whose pairwise rewiring forms the
# "block" to outline on the Δρ heatmap.
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
    anchor  = NULL,  # block of mutual co-activation, no single anchor
    members = c("MYC_TARGETS_V1", "MYC_TARGETS_V2", "E2F_TARGETS",
                "G2M_CHECKPOINT", "WNT_BETA_CATENIN_SIGNALING"),
    color   = "#762A83"
  )
)

# Fig 2 library-pair glyphs (intrinsic Hallmark × niche Hallmark).
# Maps the perturbation pairs in Fig 2 onto Hallmark module pairs interrogated
# at the pathway level by Panel G.
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

# -- 1. Per-sample Hallmark module activity (mean-z, ssGSEA-like) ------------
message("  [1/6] Per-sample Hallmark module scores (mean-z)...")

hallmark_msig <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_msig$gene_symbol, hallmark_msig$gs_name)

# Gene-level z-score across full 83-sample matrix to stabilize per-sample
# module score at small-n.
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

# Restrict to paired-subset samples.
scores_sub  <- hallmark_scores[custom_meta$sample_id, , drop = FALSE]
primary_ids <- custom_meta$sample_id[group_sub == "Primary"]
met_ids     <- custom_meta$sample_id[group_sub == "LungMet"]

# Filter to modules measurable in BOTH groups (≥10 detected genes).
valid_modules <- colnames(scores_sub)[!is.na(colSums(scores_sub))]
message("    measurable modules: ", length(valid_modules), " / 50 Hallmarks")

# Reorder to thematic order; drop modules not in MODULE_ORDER (none expected).
ordered_modules <- intersect(MODULE_ORDER, valid_modules)
extra_modules   <- setdiff(valid_modules, MODULE_ORDER)
if (length(extra_modules) > 0) {
  message("    appending non-themed modules: ",
          paste(extra_modules, collapse = ", "))
  ordered_modules <- c(ordered_modules, extra_modules)
}

# -- 2. Spearman correlation matrices per group + Δρ ------------------------
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

# Numerical source data for Supplementary Figure S1E (Variant B). This tidy
# table contains every cell of the two full, unclustered matrices displayed side
# by side, in the same thematic row/column order used by ComplexHeatmap below.
s1e_source <- expand.grid(
  row_module = ordered_modules,
  column_module = ordered_modules,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
s1e_source$primary_spearman_rho <- as.vector(cor_primary)
s1e_source$lungmet_spearman_rho <- as.vector(cor_lungmet)
s1e_source$n_primary <- length(primary_ids)
s1e_source$n_lungmet <- length(met_ids)
write.csv(
  s1e_source,
  file.path(OUTPUT_TBL_DIR,
            "FigS1E_Hallmark_coactivation_matrices_source_data.csv"),
  row.names = FALSE,
  quote = FALSE
)

# Pretty module names: drop underscore.
pretty_mod <- function(x) gsub("_", " ", x)

# -- 3. Decorator helpers ---------------------------------------------------
message("  [3/6] Building decorator helpers (callout boxes + Fig 2 glyphs)...")

# After Heatmap draws, decorate_heatmap_body lets us overlay grid graphics.
# Cell coordinates inside the body viewport: x,y in [0,1]; column index = 1..n_col,
# row index = 1..n_row; cell width/height = 1/n_col, 1/n_row.

decorate_callouts <- function(modules) {
  n <- length(modules)
  for (callout in REWIRING_CALLOUTS) {
    members <- intersect(c(callout$anchor, callout$members), modules)
    if (length(members) < 2) next
    col_idx <- which(modules %in% members)
    row_idx <- which(modules %in% members)
    # Bounding box of the smallest contiguous rectangle covering all member
    # cells (member rows × member cols). Outline only.
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
    # Mark BOTH (i, j) and (j, i) cells (matrix is symmetric Δρ).
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

# -- 4. Variant A: single Δρ heatmap with callouts + glyphs ------------------
# -- Text under the base pdf() device ----------------------------------------
# Panels are written with pdf_device() so Illustrator receives each label as its
# own editable text object (see R/devices.R). That device maps text through a
# single-byte encoding, which cannot carry Greek. Delta and rho are therefore
# drawn as plotmath expressions, which R renders from the Symbol font and
# Illustrator imports as ordinary editable text. Heatmap `name=` stays an ASCII
# identifier, because decorate_heatmap_body() looks the heatmap up by it; only
# the displayed legend title becomes an expression.
DELTA_RHO      <- expression(Delta * rho)
DELTA_RHO_DIFF <- expression(paste(Delta, rho, "  (LungMet - Primary)"))
REWIRE_TITLE   <- expression(paste("Hallmark co-activation rewiring  (",
                                   Delta, rho, " = LungMet - Primary)"))
SPEARMAN_RHO   <- expression(paste("Spearman ", rho))
HM_DELTA_RHO   <- "delta_rho"   # ASCII heatmap id, never displayed

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

pdf_device(panel_output_path("G", "variant_A",
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
# Corner legend for Fig 2 glyph
grid.text("\u2022  pathway pair interrogated in Fig 2",
          x = unit(0.99, "npc"), y = unit(2, "mm"),
          just = c("right", "bottom"),
          gp = gpar(fontsize = 8, fontface = "italic", col = "grey25"))
dev.off()

# -- 5. Variant B: the two source matrices — Supplementary Figure 1 panel E --
# Primary and lung-metastasis co-activation side by side on one shared colour
# scale. Their difference is main panel G and is deliberately NOT repeated here.
# An earlier three-matrix version put the Delta-rho in both figures and carried
# three legends, two of them identical, so the assembled panel had to be cropped
# by hand and left an orphaned legend behind. This renders what is placed.
message("  [5/6] Rendering Variant B (Supp. Fig. 1 panel E: two source matrices)...")

# `name_id` is ComplexHeatmap's internal identifier and must stay a plain
# character string; `legend_title` is what the reader sees and may be a plotmath
# expression, which is how the Greek rho survives the base pdf() encoding.
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

# Both matrices are Spearman rho over the same range, so one legend serves both.
ht_p <- make_corr_heatmap(cor_primary, "Primary Tumor Hallmark Pathway Co-activation",
                          col_corr, "rho_primary", SPEARMAN_RHO, show_legend = TRUE)
ht_l <- make_corr_heatmap(cor_lungmet, "LungMet Hallmark Pathway Co-activation",
                          col_corr, "rho_lungmet", SPEARMAN_RHO, show_legend = FALSE)

# placed = TRUE: this non-selected variant IS Supplementary Figure 1 panel E, so
# it belongs in supplementary/ proper rather than supplementary/variants/.
pdf_device(panel_output_path("G", "variant_B",
                            "FigS1_G_coactivation_matrices.pdf",
                            placed = TRUE),                 # Supp Fig 1 panel E
          width = 13, height = 6.5)
draw(ht_p + ht_l, ht_gap = unit(6, "mm"))
dev.off()

# -- 6. Variant C: Δρ-dominant + small Primary/LungMet baselines below ------
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

pdf_device(panel_output_path("G", "variant_C",
                            "Fig1.H_pathway_rewiring_variant_C_delta_dominant.pdf"),
          width = 12, height = 12)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(2.6, 1), "null"))))

# Top: large Δρ
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

# Bottom: two small baseline heatmaps side-by-side
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

# -- 7. Supplementary barplot: top |Δρ| ≥ 0.5 events ------------------------
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
  # Not placed in either assembled figure; kept for comparison.
  ggsave(file.path(OUTPUT_SUP_VAR_DIR, "Fig1.Sup_Network_delta_rewiring.pdf"),
         p_H_supp, width = 9, height = 9, device = pdf_device)
  message("    supp barplot: ", nrow(top_rewire), " rewiring events at |\u0394\u03c1| \u2265 ",
          SUP_DELTA_THRESHOLD)
}

message("[Panel G]  Done. Three variants + supp barplot rendered.")
