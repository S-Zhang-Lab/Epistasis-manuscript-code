message("[Panel A]  Pan-cancer co-occurrence chord — building...")

PANCAN_STUDIES <- c(
  BRCA     = "brca_tcga_pan_can_atlas_2018",
  LUAD     = "luad_tcga_pan_can_atlas_2018",
  LUSC     = "lusc_tcga_pan_can_atlas_2018",
  COADREAD = "coadread_tcga_pan_can_atlas_2018",
  STAD     = "stad_tcga_pan_can_atlas_2018",
  HNSC     = "hnsc_tcga_pan_can_atlas_2018",
  BLCA     = "blca_tcga_pan_can_atlas_2018",
  PRAD     = "prad_tcga_pan_can_atlas_2018",
  OV       = "ov_tcga_pan_can_atlas_2018",
  LIHC     = "lihc_tcga_pan_can_atlas_2018"
)

PANCAN_DRIVERS <- c(

  "TP53", "KRAS", "PIK3CA", "PTEN", "APC", "BRAF", "EGFR", "NRAS",
  "IDH1", "IDH2", "RB1", "CDKN2A", "CDKN2B", "MYC",

  "BRCA1", "BRCA2", "ATM", "CHEK2", "BAP1", "FANCA",

  "MLH1", "MSH2", "MSH6", "PMS2", "POLE",

  "ARID1A", "ARID2", "KMT2C", "KMT2D", "SETD2", "EP300", "CREBBP",
  "SPOP", "FOXA1", "GATA3",

  "ERBB2", "ERBB3", "FGFR2", "FGFR3", "MET", "AKT1", "STK11", "NF1",

  "CTNNB1", "SMAD4", "TGFBR2", "ACVR2A", "NOTCH1", "FBXW7",

  "VHL", "NF2", "CDH1", "CDKN1B"
)

PANCAN_DRIVER_CONTRACT <- file.path(GENE_LISTS_DIR,
                                    "panel_B_pancan_drivers_53.txt")
if (file.exists(PANCAN_DRIVER_CONTRACT)) {
  contract_genes <- readLines(PANCAN_DRIVER_CONTRACT, warn = FALSE)
  contract_genes <- trimws(sub("#.*$", "", contract_genes))
  contract_genes <- contract_genes[nzchar(contract_genes)]
  if (!setequal(PANCAN_DRIVERS, contract_genes)) {
    stop("Panel B driver panel disagrees with its data contract.\n",
         "  in the script only: ",
         paste(setdiff(PANCAN_DRIVERS, contract_genes), collapse = ", "), "\n",
         "  in the gene list only: ",
         paste(setdiff(contract_genes, PANCAN_DRIVERS), collapse = ", "), "\n",
         "  Reconcile PANCAN_DRIVERS with docs/gene_lists/",
         "panel_B_pancan_drivers_53.txt, then rebuild the committed subsets.\n",
         "  Running scripts/subset_tcga_pancan.sh with no full MAF present\n",
         "  prints the commands that fetch them:\n",
         "    bash scripts/subset_tcga_pancan.sh",
         call. = FALSE)
  }
} else {
  warning("Panel B driver contract missing (", PANCAN_DRIVER_CONTRACT,
          "); cannot confirm the committed MAF subsets cover this gene panel.")
}

PANCAN_TIER_COLORS <- c(
  `1`  = "#E0E0E0",
  `2`  = "#FDAE6B",
  `3`  = "#E6550D",
  `4+` = "#7F2704"
)

PANCAN_TIER_ALPHA <- c(`1` = 0.12, `2` = 0.25, `3` = 0.85, `4+` = 1.0)

RECURRENCE_FLOOR <- 2

recurrence_tier <- function(n) {
  ifelse(n >= 4, "4+", as.character(n))
}

alpha_hex <- function(hex, a) {
  rgb_mat <- col2rgb(hex) / 255
  rgb(rgb_mat[1, ], rgb_mat[2, ], rgb_mat[3, ], alpha = a)
}

message("  [1/6] Loading ", length(PANCAN_STUDIES), " per-lineage MAFs...")

maf_path <- function(study) {
  plain <- file.path(INPUT_DIR, "TCGA_PanCan", study, "data_mutations.txt")
  if (file.exists(plain)) plain else paste0(plain, ".gz")
}

mafs <- lapply(PANCAN_STUDIES, function(study) {
  f <- maf_path(study)
  if (!file.exists(f)) stop("Missing MAF: ", f)
  read.maf(maf = f, verbose = FALSE)
})

message("  [2/6] Pooled pan-cancer somaticInteractions on ",
        length(PANCAN_DRIVERS), " curated drivers...")
pooled_maf <- merge_mafs(mafs, verbose = FALSE)

pooled_si <- somaticInteractions(
  maf    = pooled_maf,
  genes  = PANCAN_DRIVERS,
  pvalue = c(0.05, 0.1)
)

pooled_df <- as.data.frame(pooled_si)
pooled_pairs <- pooled_df[
  pooled_df$Event == "Co_Occurence" & pooled_df$pAdj < 0.05,
  c("gene1", "gene2", "pAdj")
]
stopifnot(nrow(pooled_pairs) > 0)
message("    pooled significant positive co-occurrences: ", nrow(pooled_pairs))

message("  [3/6] Per-lineage somaticInteractions on the curated driver set...")

per_lineage_status <- vapply(names(PANCAN_STUDIES), function(ctype) {
  tryCatch({
    suppressWarnings(somaticInteractions(
      maf    = mafs[[ctype]],
      genes  = PANCAN_DRIVERS,
      pvalue = c(0.05, 0.1)
    ))
    "ok"
  }, error = function(e) {
    paste0("FAIL: ", conditionMessage(e))
  })
}, character(1))

write.csv(
  data.frame(cancer_type = names(per_lineage_status),
             status      = unname(per_lineage_status)),
  file.path(OUTPUT_TBL_DIR, "Fig1A_per_lineage_status.csv"),
  row.names = FALSE
)

failed <- per_lineage_status[per_lineage_status != "ok"]
if (length(failed) > 0) {
  stop("Per-lineage somaticInteractions failed for ", length(failed),
       " cohort(s):\n  ",
       paste(names(failed), unname(failed), sep = " — ", collapse = "\n  "),
       "\n\nFix the failing cohort(s) or remove from PANCAN_STUDIES.\n",
       "Status table: output/tables/Fig1A_per_lineage_status.csv")
}

per_lineage_si <- lapply(names(PANCAN_STUDIES), function(ctype) {
  d <- as.data.frame(somaticInteractions(
    maf    = mafs[[ctype]],
    genes  = PANCAN_DRIVERS,
    pvalue = c(0.05, 0.1)
  ))
  d$CancerType <- ctype
  d
})
per_lineage_df <- do.call(rbind, per_lineage_si)

per_lineage_sig <- per_lineage_df[
  !is.na(per_lineage_df$pAdj) &
    per_lineage_df$Event == "Co_Occurence" &
    per_lineage_df$pAdj < 0.1,
  c("gene1", "gene2", "CancerType", "pAdj")
]

message("  [4/6] Counting pan-cancer recurrence per pair...")
key <- function(a, b) apply(cbind(a, b), 1, function(x) paste(sort(x), collapse = "|"))
pooled_pairs$pair_key <- key(pooled_pairs$gene1, pooled_pairs$gene2)

if (nrow(per_lineage_sig) > 0) {
  per_lineage_sig$pair_key <- key(per_lineage_sig$gene1, per_lineage_sig$gene2)
  recur_tbl <- tapply(
    per_lineage_sig$CancerType,
    per_lineage_sig$pair_key,
    function(x) length(unique(x))
  )
} else {
  recur_tbl <- integer(0)
}

pooled_pairs$n_cancers <- as.integer(recur_tbl[pooled_pairs$pair_key])
pooled_pairs$n_cancers[is.na(pooled_pairs$n_cancers)] <- 0L

pooled_pairs$pooled_only <- pooled_pairs$n_cancers == 0L

pooled_pairs$tier <- ifelse(pooled_pairs$pooled_only,
                            "0", recurrence_tier(pooled_pairs$n_cancers))
pooled_pairs$ribbon <- character(nrow(pooled_pairs))
is_po  <- pooled_pairs$pooled_only
pooled_pairs$ribbon[is_po]  <- alpha_hex("#CCCCCC", 0.10)
pooled_pairs$ribbon[!is_po] <- alpha_hex(
  PANCAN_TIER_COLORS[pooled_pairs$tier[!is_po]],
  PANCAN_TIER_ALPHA[pooled_pairs$tier[!is_po]]
)

message("    recurrence distribution (pooled-only = no per-lineage support):")
print(table(pooled_pairs$tier))
message("    pooled-only pairs (excluded from chord by RECURRENCE_FLOOR=", RECURRENCE_FLOOR, "): ",
        sum(pooled_pairs$pooled_only))

message("  [5/6] Writing tables...")
write.csv(
  pooled_pairs[, c("gene1", "gene2", "pAdj", "n_cancers", "pooled_only", "tier")],
  file      = file.path(OUTPUT_TBL_DIR, "Fig1A_pancancer_cooccurrence_pairs.csv"),
  row.names = FALSE
)
write.csv(per_lineage_sig,
          file      = file.path(OUTPUT_TBL_DIR, "Fig1A_pancancer_cooccurrence_per_lineage.csv"),
          row.names = FALSE)

message("  [6/6] Drawing chord diagram...")

chord_pairs <- pooled_pairs[pooled_pairs$n_cancers >= RECURRENCE_FLOOR, ]

chord_pairs <- chord_pairs[order(chord_pairs$n_cancers), ]
message("    pairs displayed in chord (recurrence ≥ ", RECURRENCE_FLOOR, "): ",
        nrow(chord_pairs))

chord_df <- chord_pairs[, c("gene1", "gene2")]
genes    <- unique(c(chord_df$gene1, chord_df$gene2))

set.seed(42)
grid_colors <- structure(rand_color(length(genes), transparency = 0.2),
                         names = genes)

pdf_device(file.path(OUTPUT_FIG_DIR, "Fig1_B_pancancer_cooccurrence_chord.pdf"),
          width = 11, height = 11)
circos.clear()

circos.par(start.degree     = 90,
           gap.degree       = 0.5,
           track.margin     = c(0.005, 0.005),
           cell.padding     = c(0, 0, 0, 0),
           canvas.xlim      = c(-1.55, 1.55),
           canvas.ylim      = c(-1.55, 1.55),
           points.overflow.warning = FALSE)

gene_degree_chord <- table(c(chord_df$gene1, chord_df$gene2))
gene_degree_sorted <- names(sort(gene_degree_chord, decreasing = TRUE))
n_g <- length(gene_degree_sorted)
interleaved_genes <- character(n_g)
hi_idx <- 1; lo_idx <- n_g
for (i in seq_len(n_g)) {
  if (i %% 2L == 1L) {
    interleaved_genes[i] <- gene_degree_sorted[hi_idx]; hi_idx <- hi_idx + 1L
  } else {
    interleaved_genes[i] <- gene_degree_sorted[lo_idx]; lo_idx <- lo_idx - 1L
  }
}

chordDiagram(
  x                 = chord_df,
  grid.col          = grid_colors,
  col               = chord_pairs$ribbon,
  annotationTrack   = "grid",
  order             = interleaved_genes,
  preAllocateTracks = list(track.height = 0.20)
)

ALWAYS_LABEL_GENES <- c(
  "TP53", "BRCA1", "BRCA2", "PIK3CA", "PTEN", "MYC",
  "KRAS", "EGFR", "RB1", "APC", "ERBB2", "AKT1",
  "CDKN2A", "STK11", "IDH1", "IDH2", "NF1", "FOXA1", "GATA3"
)
LABEL_DEGREE_FLOOR <- 3
gene_degree_chord <- table(c(chord_df$gene1, chord_df$gene2))

labeled_genes <- union(
  intersect(ALWAYS_LABEL_GENES, names(gene_degree_chord)),
  names(gene_degree_chord)[gene_degree_chord >= LABEL_DEGREE_FLOOR]
)
unlabeled_genes <- setdiff(names(gene_degree_chord), labeled_genes)

message("    labeled genes (canonical override + degree >= ", LABEL_DEGREE_FLOOR,
        "): ", length(labeled_genes), " of ", length(gene_degree_chord))
if (length(unlabeled_genes) > 0) {
  message("    unlabeled (small-arc, non-canonical): ",
          paste(unlabeled_genes, collapse = ", "))
}

circos.track(
  track.index = 1,
  bg.border   = NA,
  panel.fun = function(x, y) {
    if (CELL_META$sector.index %in% labeled_genes) {
      circos.text(
        CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing     = "clockwise",
        niceFacing = TRUE,
        adj        = c(0, 0.5),
        cex        = 1.30,
        font       = 2
      )
    }
  }
)

title(main      = "Driver co-occurrence is pervasive across metastatic cancers",
      line      = -0.5,
      cex.main  = 1.5,
      font.main = 2)

legend("bottomleft",
       title       = "Recurrence\n(# cancer types)",
       legend      = names(PANCAN_TIER_COLORS),
       fill        = PANCAN_TIER_COLORS,
       border      = NA,
       bty         = "n",
       cex         = 1.30,
       title.cex   = 1.40,
       text.font   = 2)

legend("bottomright",
       legend = c(paste0("Cancer types pooled: ", length(PANCAN_STUDIES)),
                  paste0("Curated drivers (pool): ", length(PANCAN_DRIVERS)),
                  paste0("Significant pairs (pAdj<0.05): ", nrow(pooled_pairs)),
                  paste0("Pairs displayed (recurrence ≥ ", RECURRENCE_FLOOR, "): ",
                         nrow(chord_pairs))),
       bty = "n",
       cex = 1.05)

dev.off()
circos.clear()

message("[Panel A]  Done. → ",
        file.path(OUTPUT_FIG_DIR, "Fig1_B_pancancer_cooccurrence_chord.pdf"))
