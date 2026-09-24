suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

if (!exists("OUTPUT_SUP_DIR")) {
  .args <- commandArgs(trailingOnly = FALSE)
  .fa   <- sub("^--file=", "", .args[grepl("^--file=", .args)])
  REPO_ROOT <- if (length(.fa) > 0)
    normalizePath(file.path(dirname(.fa), "..", "..")) else normalizePath(getwd())
  OUTPUT_SUP_DIR <- file.path(REPO_ROOT, "output", "figures", "supplementary")
  OUTPUT_TBL_DIR <- file.path(REPO_ROOT, "output", "tables")
  DATA_RAW       <- file.path(REPO_ROOT, "data", "raw")
  dir.create(OUTPUT_SUP_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(OUTPUT_TBL_DIR, recursive = TRUE, showWarnings = FALSE)
}

if (!exists("pdf_device")) source(file.path(REPO_ROOT, "R", "devices.R"))
if (!exists("OUTPUT_SUP_DIA_DIR"))
  OUTPUT_SUP_DIA_DIR <- file.path(OUTPUT_SUP_DIR, "diagnostics")
dir.create(OUTPUT_SUP_DIA_DIR, recursive = TRUE, showWarnings = FALSE)
IN  <- file.path(DATA_RAW, "MSK_MET_2021_breast"); SUP <- OUTPUT_SUP_DIR

msk <- function(base) { g <- file.path(IN, paste0(base, ".gz"))
                        if (file.exists(g)) g else file.path(IN, base) }

sv <- function(p,f,w,h,dir=SUP){ ggsave(file.path(dir,paste0(f,".pdf")),p,width=w,height=h,device = pdf_device)
                                 ggsave(file.path(dir,paste0(f,".png")),p,width=w,height=h,dpi=200) }
th <- theme_classic(base_size=10) +
  theme(plot.title=element_text(face="bold",size=10.5),
        plot.subtitle=element_text(size=7.6,color="grey35"), legend.position="right")

SILENT <- c("Silent","Intron","3'UTR","5'UTR","3'Flank","5'Flank","IGR","RNA","Splice_Region")
SC     <- c("DMETS_DX_BONE","DMETS_DX_LIVER","DMETS_DX_LUNG","DMETS_DX_CNS_BRAIN")
HIST   <- list(Lung="DMETS_DX_LUNG", Liver="DMETS_DX_LIVER", Bone="DMETS_DX_BONE", Brain="DMETS_DX_CNS_BRAIN")
BIOP   <- list(Lung="Lung", Liver="Liver", Bone="Bone", Brain="CNS/Brain")
SITECOL<- c(Lung="#16a085", Liver="#c0392b", Bone="#e67e22", Brain="#8e44ad")

clin <- fread(msk("clinical.csv"))
mut  <- fread(msk("mutations.csv"))[!(mutationType %in% SILENT) & !is.na(hugo)]
cna  <- fread(msk("cna.csv"))[alteration %in% c(2L,-2L) & !is.na(hugo)]
pm <- rowSums(as.matrix(clin[,..SC])=="Yes")>0
cp <- clin[pm]
cm <- cp[grepl("Metas", SAMPLE_TYPE, ignore.case=TRUE)]
nmet <- nrow(cm)
yb <- lapply(BIOP, function(v) grepl(v, cm$METASTATIC_SITE, ignore.case=TRUE))
amt <- mut[,.N,by=sampleId]$N[match(cm$sampleId, mut[,.N,by=sampleId]$sampleId)]; amt[is.na(amt)] <- 0
act <- cna[,.N,by=sampleId]$N[match(cm$sampleId, cna[,.N,by=sampleId]$sampleId)]; act[is.na(act)] <- 0
altload <- amt + act

cat("[FigS1_cohort] exposure mismatch (diagnostic) ...\n")
mm <- rbindlist(lapply(names(HIST), function(nm) data.table(
  site = nm,
  `patients with history of involvement` = sum(cm[[HIST[[nm]]]]=="Yes"),
  `specimens actually biopsied there`    = sum(yb[[nm]]))))
fwrite(mm, file.path(OUTPUT_TBL_DIR, "FigS1_exposure_mismatch.csv"))
mml <- melt(mm, id.vars="site", variable.name="measure", value.name="n")
mml[, site := factor(site, levels=names(HIST))]
pMM <- ggplot(mml, aes(site, n, fill=measure)) +
  geom_col(position=position_dodge(width=.75), width=.68) +
  geom_text(aes(label=n), position=position_dodge(width=.75), vjust=-0.3, size=2.9) +
  scale_fill_manual(values=c("grey72","#16a085"), name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,.16))) +
  labs(title="The site flags record patient history, not where the specimen came from",
       subtitle=paste0("Metastasis-derived specimens (n = ",nmet,"). Only 83 of the 408 specimens from patients with a history of\n",
                       "lung involvement are lung biopsies (20%); the rest are liver, lymph node, bone, pleura, skin and brain.\n",
                       "The divergence analysis therefore uses biopsy site, and patient history is kept only as a covariate."),
       x=NULL, y="number of metastatic specimens") + th +
  theme(legend.position="top", legend.text=element_text(size=7.6))
sv(pMM, "FigS1_exposure_mismatch", 7.4, 4.4, dir = OUTPUT_SUP_DIA_DIR)

cat("[FigS1_cohort] biopsy-site composition (panel A) ...\n")
prev <- data.table(site=factor(names(BIOP), levels=names(BIOP)),
                   n=vapply(names(BIOP), function(nm) sum(yb[[nm]]), 0L))
prev[, pct := 100*n/nmet]
fwrite(prev, file.path(OUTPUT_TBL_DIR, "FigS1_cohort_by_biopsy_site.csv"))
pA <- ggplot(prev, aes(site, pct, fill=site)) + geom_col(width=.7) +
  geom_text(aes(label=sprintf("%.0f%%\n(n=%d)", pct, n)), vjust=-0.2, size=2.9) +
  scale_fill_manual(values=SITECOL, guide="none") +
  scale_y_continuous(limits=c(0,42), expand=expansion(mult=c(0,.1))) +
  labs(title="Biopsy-site composition of the metastatic specimens",
       subtitle=paste0("Site of the sequenced specimen; metastasis-derived specimens (n = ", nmet,
                       "). The remainder are lymph node,\npleura, skin, ovary and unspecified sites."),
       x=NULL, y="% of metastatic specimens") + th
sv(pA, "FigS1_cohort_by_biopsy_site", 6.2, 4.0)

cat("[FigS1_cohort] alteration load (panel B) ...\n")
cdt <- rbindlist(lapply(names(BIOP), function(nm)
  data.table(site=nm, member=ifelse(yb[[nm]],"this site","other sites"), load=altload)))
cdt[, site:=factor(site, levels=names(BIOP))][, member:=factor(member, levels=c("this site","other sites"))]
pv <- vapply(names(BIOP), function(nm)
  suppressWarnings(wilcox.test(altload[yb[[nm]]], altload[!yb[[nm]]])$p.value), 0)
fwrite(data.table(site=names(BIOP), MWU_P=pv), file.path(OUTPUT_TBL_DIR,"FigS1_alteration_load.csv"))
lab <- data.table(site=factor(names(BIOP),levels=names(BIOP)), p=pv, y=max(altload)*0.97)
pC <- ggplot(cdt, aes(member, load, fill=member)) +
  geom_boxplot(outlier.size=.35, width=.6) + facet_wrap(~site, nrow=1) +
  geom_text(data=lab, aes(x=1.5, y=y, label=sprintf("MWU P=%.2g", p)), inherit.aes=FALSE, size=2.5) +
  scale_fill_manual(values=c("this site"="#16a085","other sites"="grey75"), guide="none") +
  labs(title="Lung-derived specimens are not enriched for alteration load",
       subtitle=paste0("Per-specimen count of nonsilent mutations plus high-level copy-number alterations on a targeted panel.\n",
                       "This is not mutations per megabase. Full range shown, no observations clipped. Lung (P = 0.14) and bone\n",
                       "(P = 0.097) do not differ from the rest; liver (P = 0.023) and brain (P = 1.5e-05) do, which does not bear on\n",
                       "the lung claim. A non-significant test does not by itself establish equivalence."),
       x=NULL, y="alteration load per specimen") + th +
  theme(axis.text.x=element_text(size=7.6))
sv(pC, "FigS1_alteration_load", 9, 3.9)

cat(sprintf("[FigS1_cohort] DONE. metastatic specimens=%d | biopsy-site n: %s | lung load MWU P=%.2g\n",
            nmet, paste(sprintf("%s=%d", prev$site, prev$n), collapse=" "), pv[["Lung"]]))
