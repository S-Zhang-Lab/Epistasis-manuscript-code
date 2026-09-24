suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(circlize); library(scales)
})

if (!exists("OUTPUT_FIG_DIR")) {
  .args <- commandArgs(trailingOnly = FALSE)
  .fa   <- sub("^--file=", "", .args[grepl("^--file=", .args)])
  REPO_ROOT <- if (length(.fa) > 0)
    normalizePath(file.path(dirname(.fa), "..")) else normalizePath(getwd())
  OUTPUT_FIG_DIR <- file.path(REPO_ROOT, "output", "figures", "main")
  OUTPUT_TBL_DIR <- file.path(REPO_ROOT, "output", "tables")
  DATA_RAW       <- file.path(REPO_ROOT, "data", "raw")
  dir.create(OUTPUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(OUTPUT_TBL_DIR, recursive = TRUE, showWarnings = FALSE)
}

if (!exists("pdf_device")) source(file.path(REPO_ROOT, "R", "devices.R"))
IN  <- file.path(DATA_RAW, "MSK_MET_2021_breast")

msk <- function(base) { g <- file.path(IN, paste0(base, ".gz"))
                        if (file.exists(g)) g else file.path(IN, base) }

FIG <- OUTPUT_FIG_DIR; TBL <- OUTPUT_TBL_DIR

GENES <- c("TP53","PIK3CA","CDH1","GATA3","MAP3K1","KMT2C","ESR1","ARID1A","AKT1","NF1",
           "NCOR1","FOXA1","TBX3","RB1","ERBB2","CCND1","MYC","FGFR1","PTEN","CDKN2A")
SILENT <- c("Silent","Intron","3'UTR","5'UTR","3'Flank","5'Flank","IGR","RNA","Splice_Region")
SC     <- c("DMETS_DX_BONE","DMETS_DX_LIVER","DMETS_DX_LUNG","DMETS_DX_CNS_BRAIN")

BIOPSY_SITES <- list(Lung="Lung", Liver="Liver", Bone="Bone", Brain="CNS/Brain")
ng <- length(GENES)

clin <- fread(msk("clinical.csv"))
mut  <- fread(msk("mutations.csv"))[!(mutationType %in% SILENT) & !is.na(hugo)]
cna  <- fread(msk("cna.csv"))[alteration %in% c(2L, -2L) & !is.na(hugo)]
samples <- clin$sampleId
alt <- unique(rbind(mut[,.(sampleId,hugo)], cna[,.(sampleId,hugo)]))[sampleId %in% samples & hugo %in% GENES]
M <- matrix(0L, length(samples), ng, dimnames = list(samples, GENES))
M[cbind(match(alt$sampleId, samples), match(alt$hugo, GENES))] <- 1L
poolmask <- rowSums(as.matrix(clin[, ..SC]) == "Yes") > 0
cp <- clin[poolmask]; Mp <- M[samples[poolmask], , drop = FALSE]

lnor <- function(X){ k<-ncol(X); b<-crossprod(X); ct<-diag(b); a<-b
  bb<-matrix(ct,k,k)-a; cc<-matrix(ct,k,k,byrow=TRUE)-a; d<-nrow(X)-a-bb-cc
  log(((a+.5)*(d+.5))/((bb+.5)*(cc+.5))) }
mad_dl <- function(A,B){ L<-lnor(A)-lnor(B); mean(abs(L)[upper.tri(L)]) }

cat("[Fig1_msk] co-occurrence chord ...\n")
both <- crossprod(Mp); ct <- diag(both); N <- nrow(Mp)
A <- both; B <- matrix(ct,ng,ng)-A; Cc <- matrix(ct,ng,ng,byrow=TRUE)-A; D <- N-A-B-Cc
L <- log(((A+.5)*(D+.5))/((B+.5)*(Cc+.5)))
P <- matrix(1, ng, ng)
for(i in 1:(ng-1)) for(j in (i+1):ng){
  P[i,j] <- fisher.test(matrix(c(A[i,j],B[i,j],Cc[i,j],D[i,j]), 2))$p.value; P[j,i] <- P[i,j] }
q <- p.adjust(P[upper.tri(P)], "BH"); Q <- matrix(1,ng,ng)
Q[upper.tri(Q)] <- q; Q[lower.tri(Q)] <- t(Q)[lower.tri(Q)]
sig <- which(upper.tri(P) & Q < 0.10, arr.ind = TRUE)
ed  <- data.table(from=GENES[sig[,1]], to=GENES[sig[,2]], lnor=L[sig], q=Q[sig])
ed[, value := pmin(abs(lnor), 2)][, col := ifelse(lnor > 0, "#b2182b", "#2166ac")]
fwrite(ed[order(q)], file.path(TBL, "Fig1C_MBC_cooccurrence_edges.csv"))
cat(sprintf("   edges BH<0.10: %d (%d co-occur, %d exclusive)\n",
            nrow(ed), sum(ed$lnor>0), sum(ed$lnor<0)))

gcol <- setNames(colorRampPalette(c("#16a085","#8e44ad","#e67e22"))(ng), GENES)
draw_circos <- function(){
  circos.clear(); circos.par(gap.degree=2, start.degree=90, canvas.ylim=c(-1.58,1.0))
  chordDiagram(ed[,.(from,to,value)], grid.col=gcol, col=ed$col, transparency=0.25,
               annotationTrack="grid", preAllocateTracks=list(track.height=0.10),
               directional=0, link.lwd=0.5, link.border=NA)
  circos.trackPlotRegion(track.index=1, panel.fun=function(x,y){
    circos.text(CELL_META$xcenter, CELL_META$ylim[1]+0.15, CELL_META$sector.index,
                facing="clockwise", niceFacing=TRUE, adj=c(0,0.5), cex=0.65)}, bg.border=NA)
  lx <- -0.98
  text(lx,-1.075,"Chord = significant gene-pair association:", adj=c(0,0.5), cex=0.62, font=2)
  rect(lx,-1.205, lx+0.055,-1.155, col="#b2182b", border=NA)
  text(lx+0.085,-1.18,"co-occurrence  (OR > 1, altered together)", adj=c(0,0.5), cex=0.6)
  rect(lx,-1.315, lx+0.055,-1.265, col="#2166ac", border=NA)
  text(lx+0.085,-1.29,"mutual exclusivity  (OR < 1, rarely together)", adj=c(0,0.5), cex=0.6)
  text(lx,-1.415,"Chord width = log odds ratio magnitude:", adj=c(0,0.5), cex=0.6)
  segments(lx+0.44,-1.415, lx+0.55,-1.415, lwd=1, col="#555"); text(lx+0.575,-1.415,"0.5", adj=c(0,0.5), cex=0.56)
  segments(lx+0.70,-1.415, lx+0.81,-1.415, lwd=5, col="#555"); text(lx+0.835,-1.415,"2 or more", adj=c(0,0.5), cex=0.56)
  text(lx,-1.505,"Arc = gene; arc color = gene identity (no quantitative meaning); arc length = total connectivity.",
       adj=c(0,0.5), cex=0.52, col="#555")
  text(lx,-1.565,"Fisher exact test, Benjamini-Hochberg FDR < 0.10;  MSK-MET breast, n = 1,512.",
       adj=c(0,0.5), cex=0.52, col="#555")
}
for(dev in c("pdf","png")){
  fn <- file.path(FIG, paste0("Fig1_C_MBC_cooccurrence_chord.", dev))
  if(dev=="pdf") pdf_device(fn, width=6.6, height=8.5) else png(fn, width=6.6, height=8.5, units="in", res=200)
  draw_circos(); title(main="Co-occurrence architecture of metastatic breast cancer",
                       cex.main=1.0, font.main=2, adj=0); dev.off()
}
circos.clear()

cat("[Fig1_msk] niche divergence by biopsy site ...\n")
ismet <- grepl("Metas", cp$SAMPLE_TYPE, ignore.case = TRUE)
im <- which(ismet); Xm <- Mp[im, , drop = FALSE]; cm <- cp[im]
subtype <- ifelse(is.na(cm$SUBTYPE_ABBREVIATION), "NA", cm$SUBTYPE_ABBREVIATION)
burden  <- suppressWarnings(as.numeric(cm$MET_SITE_COUNT)); burden[is.na(burden)] <- median(burden, na.rm=TRUE)
bt      <- as.character(as.integer(cut(burden, quantile(burden, c(0,1/3,2/3,1)), include.lowest=TRUE)))
CTRLS <- list("none"=NULL, "+ burden"=bt, "+ subtype"=subtype, "+ both"=paste(subtype, bt, sep="|"))

set.seed(1); N_PERM <- 2500
res <- list(); PERM_BOTH <- list()
for(nm in names(BIOPSY_SITES)){
  y <- grepl(BIOPSY_SITES[[nm]], cm$METASTATIC_SITE, ignore.case = TRUE)
  for(cn in names(CTRLS)){
    obs <- mad_dl(Xm[y,,drop=FALSE], Xm[!y,,drop=FALSE]); st <- CTRLS[[cn]]; pp <- numeric(N_PERM)
    if(is.null(st)){ n<-length(y); nL<-sum(y)
      for(b in 1:N_PERM){ s2<-sample.int(n); pp[b]<-mad_dl(Xm[s2[1:nL],,drop=FALSE], Xm[s2[(nL+1):n],,drop=FALSE]) }
    } else { lv<-split(seq_along(y), st); nv<-vapply(lv, function(k) sum(y[k]), integer(1))
      for(b in 1:N_PERM){ yp<-logical(length(y))
        for(j in seq_along(lv)){ k<-lv[[j]]; if(nv[j]>0) yp[k[sample(length(k), nv[j])]]<-TRUE }
        pp[b]<-mad_dl(Xm[yp,,drop=FALSE], Xm[!yp,,drop=FALSE]) } }
    if(cn == "+ both") PERM_BOTH[[nm]] <- list(obs=obs, perm=pp, n=sum(y))
    res[[length(res)+1]] <- data.table(site=nm, n=sum(y), ctrl=cn, obs=obs,
                                       P=(sum(pp>=obs)+1)/(N_PERM+1), z=(obs-mean(pp))/sd(pp))
  }}
R <- rbindlist(res)
R[, q := p.adjust(P, "BH"), by = ctrl]
R[, site := factor(site, levels = rev(names(BIOPSY_SITES)))]
R[, ctrl := factor(ctrl, levels = names(CTRLS))]
R[, lab := sprintf("%.3f\n(q=%.2f)", P, q)]
fwrite(R[order(ctrl, site)], file.path(TBL, "Fig1_niche_divergence_biopsy_site.csv"))

pE <- ggplot(R, aes(ctrl, site, fill = z)) +
  geom_tile(color="white", linewidth=1.2) +
  geom_text(aes(label=lab), size=2.6, lineheight=0.9,
            fontface=ifelse(R$q<0.05,"bold","plain"),
            color=ifelse(R$z>2.2,"white","grey15")) +
  scale_fill_gradient2(low="#f7fbff", mid="#9ecae1", high="#08519c", midpoint=1,
                       limits=c(-1,3), oob=scales::squish, name="SD above\nchance") +
  scale_y_discrete(labels=function(x) paste0(x, "\n(n=", R$n[match(x, as.character(R$site))], ")")) +
  labs(title="Only lung-derived metastases diverge in co-occurrence structure",
       subtitle=paste0("Exposure = biopsy site of the sequenced specimen; metastasis-derived biopsies only (n = ",
                       length(im), ").\nCell = permutation P, with Benjamini-Hochberg q across the four sites.",
                       "  Bold = q < 0.05.  Exploratory."),
       x="confounders controlled for", y=NULL) +
  coord_equal() + theme_minimal(base_size=10) +
  theme(panel.grid=element_blank(), plot.title=element_text(face="bold", size=11),
        plot.subtitle=element_text(size=7.4, color="grey35"),
        axis.text=element_text(size=9), legend.key.width=unit(8,"pt"))
ggsave(file.path(FIG,"Fig1_E_niche_divergence_heatmap.pdf"), pE, width=7.2, height=4.0, device = pdf_device)
ggsave(file.path(FIG,"Fig1_E_niche_divergence_heatmap.png"), pE, width=7.2, height=4.0, dpi=200)

cat("[Fig1_msk] matched-null histograms (same computation as panel E) ...\n")
Rboth <- R[ctrl == "+ both"]
hd <- rbindlist(lapply(names(BIOPSY_SITES), function(nm){
  r  <- PERM_BOTH[[nm]]; rr <- Rboth[site == nm]
  data.table(site = nm,
             lab  = sprintf("%s (n=%d)   P=%.3f, q=%.2f", nm, r$n, rr$P, rr$q),
             obs  = r$obs, x = r$perm) }))
hd[, site := factor(site, levels = names(BIOPSY_SITES))]
hd[, lab  := factor(lab,  levels = unique(lab))]
ol <- hd[, .(obs = obs[1], site = site[1]), by = lab]
SITECOL <- c(Lung="#16a085", Liver="#c0392b", Bone="#e67e22", Brain="#8e44ad")
pD <- ggplot(hd, aes(x)) +
  geom_histogram(bins = 38, fill = "grey82") +
  geom_vline(data = ol, aes(xintercept = obs, color = site), linewidth = 1.05, show.legend = FALSE) +
  geom_text(data = ol, aes(x = obs, y = Inf, label = " observed", color = site),
            hjust = 0, vjust = 1.6, size = 2.5, show.legend = FALSE) +
  scale_color_manual(values = SITECOL) +
  facet_wrap(~lab, nrow = 2, scales = "free") +
  labs(title = "Only lung-derived metastases fall outside their matched null",
       subtitle = paste0("Grey, distribution of the divergence statistic over ", N_PERM,
                         " permutations of the site label within\nintrinsic-subtype and metastatic-burden strata; coloured line, observed value.\n",
                         "Same computation as panel E, fully controlled model. Benjamini-Hochberg q is across the four sites."),
       x = "difference in driver co-occurrence structure\n(mean absolute change in log odds ratio, 190 gene pairs)", y = "permutations") +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 7.0, color = "grey35"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 7.8, hjust = 0),
        panel.spacing = unit(9, "pt"))
ggsave(file.path(FIG,"Fig1_D_matched_null_histograms.pdf"), pD, width=6.8, height=5.2, device = pdf_device)
ggsave(file.path(FIG,"Fig1_D_matched_null_histograms.png"), pD, width=6.8, height=5.2, dpi=200)

cat(sprintf("[Fig1_msk] DONE. lung P (none..both) = %s | lung q (none..both) = %s\n",
            paste(sprintf("%.3f", R[site=="Lung"][order(ctrl)]$P), collapse="/"),
            paste(sprintf("%.3f", R[site=="Lung"][order(ctrl)]$q), collapse="/")))
cat("[Fig1_msk] panels D and E share one computation; values agree by construction:\n")
print(Rboth[, .(site, n, P=signif(P,3), q=signif(q,3))])
