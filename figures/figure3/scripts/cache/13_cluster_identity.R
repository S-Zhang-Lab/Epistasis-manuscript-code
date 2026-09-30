#!/usr/bin/env Rscript
# =============================================================================
# Fig 3 — cluster IDENTITY by UCell Hallmark enrichment (each context separate).
# Per-cluster mean UCell (reused from fig3_singlecell_epistasis/objects/ucell_scores.rds),
# z-scored across clusters -> identity heatmap + PROPOSED biologist-readable labels
# (curate manually). Lets in-vitro & in-vivo cluster phenotypes be cross-read by
# biology, WITHOUT forcing cluster-to-cluster mapping. Convention: fig3_theme.R.
# v7 repo overhaul 2026-07-15: relocated from XL_code/fig3_panels_new_build/build_cluster_identity.R;
#   cache builder (no panel; the identity DotPlot is FigS3C's job) — reads
#   ucell_scores.rds via cache() + both raw objects; emits label + z-score CSVs
#   under cache() with basenames preserved.
# =============================================================================
if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
set.seed(1234)
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(ggplot2) })

U <- readRDS(cache("ucell_scores.rds"))  # $vivo, $vitro : cells x 50 Hallmark
# Hallmark -> broad biological process (for PROPOSED labels; curate)
CAT <- c(
  E2F_TARGETS="Proliferative", G2M_CHECKPOINT="Proliferative", MYC_TARGETS_V1="Proliferative",
  MYC_TARGETS_V2="Proliferative", MITOTIC_SPINDLE="Proliferative", DNA_REPAIR="Proliferative",
  INTERFERON_ALPHA_RESPONSE="Inflammatory/IFN", INTERFERON_GAMMA_RESPONSE="Inflammatory/IFN",
  INFLAMMATORY_RESPONSE="Inflammatory/IFN", TNFA_SIGNALING_VIA_NFKB="Inflammatory/IFN",
  IL6_JAK_STAT3_SIGNALING="Inflammatory/IFN", IL2_STAT5_SIGNALING="Inflammatory/IFN",
  COMPLEMENT="Inflammatory/IFN", ALLOGRAFT_REJECTION="Inflammatory/IFN",
  EPITHELIAL_MESENCHYMAL_TRANSITION="EMT/matrix", ANGIOGENESIS="EMT/matrix",
  APICAL_JUNCTION="EMT/matrix", COAGULATION="EMT/matrix", APICAL_SURFACE="EMT/matrix",
  HYPOXIA="Hypoxia/glycolysis", GLYCOLYSIS="Hypoxia/glycolysis",
  OXIDATIVE_PHOSPHORYLATION="OXPHOS/metabolic", FATTY_ACID_METABOLISM="OXPHOS/metabolic",
  ADIPOGENESIS="OXPHOS/metabolic", BILE_ACID_METABOLISM="OXPHOS/metabolic",
  XENOBIOTIC_METABOLISM="OXPHOS/metabolic", HEME_METABOLISM="OXPHOS/metabolic",
  PEROXISOME="OXPHOS/metabolic", CHOLESTEROL_HOMEOSTASIS="OXPHOS/metabolic",
  UNFOLDED_PROTEIN_RESPONSE="Stress/proteostasis", REACTIVE_OXYGEN_SPECIES_PATHWAY="Stress/proteostasis",
  P53_PATHWAY="Stress/proteostasis", APOPTOSIS="Stress/proteostasis", PROTEIN_SECRETION="Stress/proteostasis",
  MTORC1_SIGNALING="mTORC1/anabolic", PI3K_AKT_MTOR_SIGNALING="mTORC1/anabolic",
  WNT_BETA_CATENIN_SIGNALING="Signaling/dev", NOTCH_SIGNALING="Signaling/dev",
  HEDGEHOG_SIGNALING="Signaling/dev", TGF_BETA_SIGNALING="Signaling/dev",
  KRAS_SIGNALING_UP="Signaling/dev", KRAS_SIGNALING_DN="Signaling/dev",
  ANDROGEN_RESPONSE="Hormone/other", ESTROGEN_RESPONSE_EARLY="Hormone/other",
  ESTROGEN_RESPONSE_LATE="Hormone/other", UV_RESPONSE_UP="Hormone/other", UV_RESPONSE_DN="Hormone/other",
  MYOGENESIS="Hormone/other", SPERMATOGENESIS="Hormone/other", PANCREAS_BETA_CELLS="Hormone/other")
short <- function(v) sub("HALLMARK_","",v)

build <- function(uc, file, tag){
  o <- readRDS(raw(file))
  cl <- setNames(as.character(o$joint_clusters), colnames(o))
  uc <- uc[names(cl)[names(cl) %in% rownames(uc)], , drop=FALSE]
  clv <- cl[rownames(uc)]
  colnames(uc) <- short(colnames(uc))
  # per-cluster mean, z across clusters
  mu <- t(apply(t(uc), 1, function(x) tapply(x, clv, mean)))          # hallmark x cluster
  z  <- t(scale(t(mu)))                                               # z per hallmark across clusters
  clord <- as.character(sort(as.integer(colnames(z))))
  z <- z[, clord]
  # top-variable hallmarks for the heatmap
  vv <- apply(z,1,function(x) max(x)-min(x)); topv <- names(sort(vv,decreasing=TRUE))[1:22]
  # per-cluster top enriched hallmarks + proposed label
  lab <- do.call(rbind, lapply(clord, function(cc){
    top <- names(sort(z[,cc], decreasing=TRUE))[1:3]
    propose <- CAT[top[1]]; if (is.na(propose)) propose <- CAT[top[2]]   # key off the top hit
    data.frame(context=tag, cluster=cc, n=sum(clv==cc),
               top1=top[1], top2=top[2], top3=top[3],
               proposed_label=ifelse(is.na(propose),"(curate)",propose), row.names=NULL)
  }))
  write.csv(lab, cache(paste0("cluster_identity_labels_", tag, ".csv")), row.names=FALSE)
  write.csv(cbind(hallmark=rownames(z), round(z,3)),
            cache(paste0("cluster_identity_zscore_", tag, ".csv")), row.names=FALSE)
  # NOTE (2026-07-09): the identity HEATMAP is RETIRED. The Seurat-convention
  # DotPlot (build_cluster_dotplot.R) is the identity panel (supplement). This
  # script is kept only to emit the label + z-score CSVs that feed the dotplot
  # interpretation; `topv` is still used by the dotplot's program selection.
  cat("\n===", tag, "proposed cluster labels ===\n"); print(lab[,c("cluster","n","top1","proposed_label")], row.names=FALSE)
  rm(o); gc(verbose=FALSE)
}
build(U$vitro, "Perturb_InVitro_HTO.rds", "in_vitro")
build(U$vivo,  "Perturb_InVivo_merged.rds","in_vivo")
cat("\ncluster identity done.\n")
