message("[Panel F]  Tumor cell-microenvironment interface gene downregulation — building...")

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(msigdbr)
})

NICHE_HUGO_FILE <- file.path(GENE_LISTS_DIR, "panel_F_niche_engagement_HUGO_curated.csv")

N_LABEL <- 14

HUGO_PALETTE <- c(
  Cytokine            = "#C6322A",
  Chemokine_ligand    = "#E58978",
  Chemokine_receptor  = "#E5A638",
  MHC_class_I         = "#7E4FA0",
  MHC_class_II        = "#A78BC1",
  Adhesion            = "#2A6FB0",
  ECM                 = "#3FA67D",
  Checkpoint_ligand   = "#B83D8C"
)

MSIG_PALETTE <- c(
  HALLMARK_INFLAMMATORY_RESPONSE              = "#C6322A",
  HALLMARK_TNFA_SIGNALING_VIA_NFKB            = "#E58978",
  HALLMARK_INTERFERON_ALPHA_RESPONSE          = "#E5A638",
  HALLMARK_INTERFERON_GAMMA_RESPONSE          = "#7E4FA0",
  HALLMARK_ALLOGRAFT_REJECTION                = "#A78BC1",
  KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION = "#2A6FB0",
  KEGG_CELL_ADHESION_MOLECULES_CAMS                = "#3FA67D",
  KEGG_ECM_RECEPTOR_INTERACTION               = "#B83D8C"
)

.parse_series_matrix <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con))
  lines <- readLines(con)
  get_row <- function(tag) {
    r <- grep(paste0("^!", tag), lines, value = TRUE)
    if (!length(r)) return(NULL)
    unlist(strsplit(sub(paste0("^!", tag, "\\s*"), "", r[1]), "\t"))
  }
  title <- gsub("\"", "", get_row("Sample_title"))
  ch_rows <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
  ch <- lapply(ch_rows, function(r)
    gsub("\"", "", unlist(strsplit(sub("^!Sample_characteristics_ch1\\s*", "", r), "\t"))))
  ch_mat <- do.call(rbind, ch)
  get_char <- function(prefix) {
    row <- ch_mat[grep(paste0("^", prefix), ch_mat[, 1])[1], ]
    if (is.null(row) || all(is.na(row))) return(rep(NA, length(title)))
    sub(paste0("^", prefix, ":\\s*"), "", row)
  }
  data.frame(
    title                = title,
    `patient id:ch1`     = get_char("patient id"),
    `tumor location:ch1` = get_char("tumor location"),
    check.names          = FALSE,
    stringsAsFactors     = FALSE
  )
}

message("  [1/7] Parsing GEO series matrices for full-cohort sample annotation...")
pdata_gpl11154 <- .parse_series_matrix(file.path(INPUT_DIR, "GSE110590",
                                                 "GSE110590-GPL11154_series_matrix.txt.gz"))
pdata_gpl11154$Platform <- "GPL11154"
pdata_gpl16791 <- .parse_series_matrix(file.path(INPUT_DIR, "GSE110590",
                                                 "GSE110590-GPL16791_series_matrix.txt.gz"))
pdata_gpl16791$Platform <- "GPL16791"
pdata_full <- rbind(pdata_gpl11154, pdata_gpl16791)

expr_samples      <- colnames(exprs_mat)
expr_samples_norm <- normalize_sample_name(expr_samples)

sample_anno_full <- pdata_full %>%
  select(title, `patient id:ch1`, `tumor location:ch1`, Platform) %>%
  filter(`tumor location:ch1` %in% c("PRIMARY", "LUNG")) %>%
  mutate(NormalizedTitle = normalize_sample_name(title),
         MatchedSample   = expr_samples[match(NormalizedTitle, expr_samples_norm)]) %>%
  filter(!is.na(MatchedSample),
         !(title %in% "A12-LIV-MET-RNA")) %>%
  transmute(Sample   = MatchedSample,
            Patient  = `patient id:ch1`,
            Group    = ifelse(`tumor location:ch1` == "PRIMARY",
                              "Primary", "LungMet"),
            Platform = Platform)

n_primary <- sum(sample_anno_full$Group == "Primary")
n_lungmet <- sum(sample_anno_full$Group == "LungMet")
message("    full-cohort samples: ", n_primary, " primary + ",
        n_lungmet, " lung met (n = ", nrow(sample_anno_full), ")")

message("  [2/7] Fitting unpaired limma on full GSE110590 cohort...")
exprs_full <- exprs_mat[, sample_anno_full$Sample, drop = FALSE]
exprs_full <- apply(exprs_full, 2, as.numeric)
rownames(exprs_full) <- rownames(exprs_mat)
exprs_full <- exprs_full[complete.cases(exprs_full), ]
exprs_full <- exprs_full[apply(exprs_full, 1, var) > 1e-6, ]

group_full  <- factor(sample_anno_full$Group, levels = c("Primary", "LungMet"))
design_full <- model.matrix(~ group_full)
fit_full    <- eBayes(lmFit(exprs_full, design_full), trend = TRUE)
deg_full    <- topTable(fit_full, coef = 2, number = Inf) %>%
                tibble::rownames_to_column("Gene") %>%
                arrange(P.Value)

deg_full <- deg_full[!grepl("^ENTREZ_|\\.[0-9]+$", deg_full$Gene), ]

write.csv(deg_full,
          file.path(OUTPUT_TBL_DIR, "Fig1F_full_cohort_DEG.csv"),
          row.names = FALSE)
message("    DEG genes: ", nrow(deg_full),
        "  | adj.P.Val < 0.05 down: ",
        sum(deg_full$adj.P.Val < 0.05 & deg_full$logFC < 0))

message("  [3/7] Loading HUGO + curated niche-engagement gene list...")
.read_hugo_csv <- function(path) {
  if (!file.exists(path)) stop("Missing HUGO niche-engagement CSV: ", path)
  raw <- readLines(path, warn = FALSE)
  raw <- raw[!grepl("^#", raw)]
  read.csv(text = paste(raw, collapse = "\n"),
           stringsAsFactors = FALSE)
}
niche_hugo <- .read_hugo_csv(NICHE_HUGO_FILE)
niche_hugo$Category <- factor(niche_hugo$Category, levels = names(HUGO_PALETTE))
message("    HUGO niche-engagement genes: ", nrow(niche_hugo),
        "  (categories: ", paste(unique(niche_hugo$Category), collapse = ", "), ")")

MSIG_SETS <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "KEGG_CELL_ADHESION_MOLECULES_CAMS",
  "KEGG_ECM_RECEPTOR_INTERACTION"
)

message("  [4/7] Pulling MSigDB Hallmark + KEGG niche-engagement gene sets...")

hallmark_df <- msigdbr(species = "Homo sapiens", category = "H") %>%
  filter(gs_name %in% MSIG_SETS) %>%
  select(Source = gs_name, Gene = gene_symbol)

.try_kegg <- function() {
  out <- tryCatch(
    msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY"),
    error = function(e) NULL
  )
  if (is.null(out) || nrow(out) == 0) {
    out <- tryCatch(
      msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG"),
      error = function(e) NULL
    )
  }
  out
}
kegg_all <- .try_kegg()
if (is.null(kegg_all)) {
  warning("[Panel F] MSigDB KEGG subcollection not retrievable in this msigdbr build; ",
          "MSigDB variant will fall back to Hallmark sets only.")
  kegg_df <- data.frame(Source = character(0), Gene = character(0))
} else {
  kegg_df <- kegg_all %>%
    filter(gs_name %in% MSIG_SETS) %>%
    select(Source = gs_name, Gene = gene_symbol)
}

niche_msig <- bind_rows(hallmark_df, kegg_df) %>%
  distinct(Gene, Source) %>%
  mutate(Source = factor(Source, levels = MSIG_SETS))

message("    MSigDB union genes (with source): ",
        length(unique(niche_msig$Gene)), " unique (",
        nrow(niche_msig), " gene-source rows across ",
        length(unique(niche_msig$Source)), " source pathways)")

write.csv(niche_msig,
          file.path(OUTPUT_TBL_DIR, "Fig1F_niche_engagement_MSigDB_membership.csv"),
          row.names = FALSE)

message("  [5/7] Joining niche-engagement gene sets onto the DEG table...")

deg_hugo <- deg_full %>%
  inner_join(niche_hugo, by = "Gene") %>%
  arrange(P.Value)
write.csv(deg_hugo,
          file.path(OUTPUT_TBL_DIR, "Fig1F_niche_engagement_HUGO_DEG.csv"),
          row.names = FALSE)
message("    HUGO niche-engagement genes detected in DEG matrix: ",
        nrow(deg_hugo))

niche_msig_first <- niche_msig %>%
  group_by(Gene) %>%
  summarise(Source = first(as.character(Source)), .groups = "drop")
deg_msig <- deg_full %>%
  inner_join(niche_msig_first, by = "Gene") %>%
  mutate(Source = factor(Source, levels = MSIG_SETS)) %>%
  arrange(P.Value)
write.csv(deg_msig,
          file.path(OUTPUT_TBL_DIR, "Fig1F_niche_engagement_MSigDB_DEG.csv"),
          row.names = FALSE)
message("    MSigDB niche-engagement genes detected in DEG matrix: ",
        nrow(deg_msig))

.render_niche_panel <- function(deg_full, niche_df, palette,
                                group_col, title, subtitle,
                                variant_legend_title,
                                variant_legend_break_relabel = NULL) {
  group_sym <- rlang::sym(group_col)
  niche_genes <- niche_df$Gene

  bg_df    <- deg_full %>% filter(!Gene %in% niche_genes)
  fg_df    <- deg_full %>%
                inner_join(niche_df, by = "Gene") %>%
                mutate(grp_factor = factor(as.character(!!group_sym),
                                            levels = names(palette)))

  ks_res <- ks.test(fg_df$logFC, bg_df$logFC,
                    alternative = "greater")

  wc_res <- wilcox.test(fg_df$logFC, bg_df$logFC,
                        alternative = "less", conf.int = FALSE)

  med_bg    <- median(bg_df$logFC, na.rm = TRUE)
  med_niche <- median(fg_df$logFC, na.rm = TRUE)

  to_label <- fg_df %>%
    filter(logFC < 0) %>%
    arrange(P.Value) %>%
    head(N_LABEL) %>%
    mutate(label_text = Gene)

  p_density <- ggplot() +
    geom_density(data = bg_df,
                 aes(x = logFC), fill = "grey80", color = "grey60",
                 alpha = 0.7, linewidth = 0.3) +
    geom_density(data = fg_df,
                 aes(x = logFC), fill = "#C6322A", color = "#7E1A12",
                 alpha = 0.45, linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey30", linewidth = 0.4) +
    geom_vline(xintercept = c(med_bg, med_niche),
               linetype = "dotted",
               color = c("grey25", "#7E1A12"),
               linewidth = 0.4) +
    annotate("text",
             x = -Inf, y = Inf,
             hjust = -0.05, vjust = 1.5,

             label = sprintf(
               "Median log2FC: niche-engagement %+.3f vs all genes %+.3f\n(descriptive shift; genes are co-regulated, not independent)",
               med_niche, med_bg),
             size = 2.7, color = "grey20", lineheight = 0.95) +
    scale_x_continuous(limits = NULL, expand = expansion(mult = 0.02)) +
    labs(y = "density",
         x = NULL,
         title    = title,
         subtitle = subtitle) +
    theme_classic(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "grey25"),
          axis.text.x   = element_blank(),
          axis.ticks.x  = element_blank(),
          axis.title.x  = element_blank(),
          axis.title.y  = element_text(size = 8),
          axis.text.y   = element_text(size = 8))

  legend_breaks <- levels(droplevels(fg_df$grp_factor))
  legend_labels <- legend_breaks
  if (!is.null(variant_legend_break_relabel)) {
    legend_labels <- variant_legend_break_relabel[legend_breaks]
  }

  p_vol <- ggplot() +
    geom_vline(xintercept = c(-1, 0, 1),
               linetype = c("dashed", "solid", "dashed"),
               color    = c("grey55", "grey30", "grey55"),
               linewidth = c(0.4, 0.5, 0.4)) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dotted", color = "grey55",
               linewidth = 0.4) +
    geom_point(data = bg_df,
               aes(x = logFC, y = -log10(P.Value)),
               color = "grey85", alpha = 0.35, size = 0.55) +
    geom_point(data = fg_df,
               aes(x = logFC, y = -log10(P.Value),
                   color = grp_factor,
                   size  = pmin(8, -log10(P.Value))),
               alpha = 0.92) +
    geom_text_repel(
      data = to_label,
      aes(x = logFC, y = -log10(P.Value), label = label_text),
      size                = 2.9,
      color               = "black",
      box.padding         = 0.45,
      point.padding       = 0.30,
      segment.color       = "grey50",
      segment.size        = 0.25,
      max.overlaps        = Inf,
      min.segment.length  = 0,
      show.legend         = FALSE
    ) +
    scale_color_manual(name   = variant_legend_title,
                       values = palette,
                       breaks = legend_breaks,
                       labels = legend_labels,
                       drop   = FALSE) +
    scale_size_continuous(name = "-log10(P)",
                          range = c(1.0, 4.0),
                          guide = guide_legend(override.aes = list(alpha = 1))) +
    labs(x = expression(log[2]*"(LungMet / Primary)"),
         y = expression(-log[10]*"(P)"),
         caption = paste0(
           "Unpaired moderated t-test (limma::eBayes, trend=TRUE) on full GSE110590 cohort (",
           n_primary + n_lungmet,
           " samples; ", n_primary, " primary + ", n_lungmet,
           " lung met).\nDashed verticals: 2-fold change. Dotted horizontal: P = 0.05.\n",
           "Background grey = all expressed genes; colored points = niche-engagement subset."
         )) +
    theme_classic(base_size = 11) +
    theme(plot.caption     = element_text(size = 8, color = "grey35", hjust = 0),
          legend.key.size  = unit(0.7, "lines"),
          legend.position  = "right",
          axis.title       = element_text(size = 10),
          axis.text        = element_text(color = "black"))

  composite <- (p_density / p_vol) +
    plot_layout(heights = c(1.0, 3.0))
  composite
}

message("  [7/7] Rendering both variants and routing per docs/selected_main_variants.yml...")

selected <- read_main_variant("F")
message("    selected main variant for Panel F: ", selected)

title_hugo <- "Tumor cell–microenvironment interface genes are broadly shifted downward in lung metastasis"
subtitle_hugo <- paste0(
  "HUGO-curated cytokines, chemokines, MHC, adhesion, ECM, and inhibitory ligands ",
  "(n = ", nrow(deg_hugo), " genes detected). Discovery model: specimen-level, deliberately ",
  "permissive, to nominate a wide candidate set for the combinatorial library (Fig. 2); ",
  "it is not an independence-corrected significance test. Patient- and platform-aware ",
  "sensitivity is reported in the supplement."
)
hugo_relabel <- c(
  Cytokine            = "Cytokine",
  Chemokine_ligand    = "Chemokine ligand",
  Chemokine_receptor  = "Chemokine receptor",
  MHC_class_I         = "MHC class I",
  MHC_class_II        = "MHC class II",
  Adhesion            = "Adhesion (CAM / cadherin / integrin)",
  ECM                 = "ECM",
  Checkpoint_ligand   = "Checkpoint ligand"
)
plot_hugo <- .render_niche_panel(
  deg_full = deg_full,
  niche_df = niche_hugo,
  palette  = HUGO_PALETTE,
  group_col = "Category",
  title    = title_hugo,
  subtitle = subtitle_hugo,
  variant_legend_title = "HUGO category",
  variant_legend_break_relabel = hugo_relabel
)

ggsave(
  filename = panel_output_path("F", "hugo",
                               "Fig1_I_niche_engagement_volcano.pdf"),
  plot     = plot_hugo,
  width    = 9.0,
  height   = 8.0,
  device = pdf_device
)

title_msig <- "Niche-engagement gene sets (Hallmark + KEGG) are broadly shifted downward in lung metastasis"
subtitle_msig <- paste0(
  "Union of 5 immune-axis Hallmarks and 3 niche-interaction KEGG sets ",
  "(n = ", nrow(deg_msig), " unique genes detected). The niche-engagement distribution ",
  "sits below the background distribution. Discovery model: specimen-level, deliberately ",
  "permissive, to nominate a wide candidate set for the combinatorial library (Fig. 2); ",
  "it is not an independence-corrected significance test. Patient- and platform-aware ",
  "sensitivity is reported in the supplement."
)
msig_relabel <- c(
  HALLMARK_INFLAMMATORY_RESPONSE              = "Hallmark: Inflammatory response",
  HALLMARK_TNFA_SIGNALING_VIA_NFKB            = "Hallmark: TNFA via NFKB",
  HALLMARK_INTERFERON_ALPHA_RESPONSE          = "Hallmark: IFN-alpha response",
  HALLMARK_INTERFERON_GAMMA_RESPONSE          = "Hallmark: IFN-gamma response",
  HALLMARK_ALLOGRAFT_REJECTION                = "Hallmark: Allograft rejection",
  KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION = "KEGG: Cytokine-receptor interaction",
  KEGG_CELL_ADHESION_MOLECULES_CAMS                = "KEGG: Cell adhesion molecules",
  KEGG_ECM_RECEPTOR_INTERACTION               = "KEGG: ECM-receptor interaction"
)
plot_msig <- .render_niche_panel(
  deg_full = deg_full,
  niche_df = niche_msig_first,
  palette  = MSIG_PALETTE,
  group_col = "Source",
  title    = title_msig,
  subtitle = subtitle_msig,
  variant_legend_title = "MSigDB source",
  variant_legend_break_relabel = msig_relabel
)

ggsave(
  filename = panel_output_path("F", "msigdb",
                               "Fig1.F_niche_engagement_MSigDB_union.pdf"),
  plot     = plot_msig,
  width    = 9.5,
  height   = 8.0,
  device = pdf_device
)

message("[Panel F]  Done. Selected variant routed to main/, the other to supplementary/.")
