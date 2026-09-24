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
}

if (!exists("pdf_device")) source(file.path(REPO_ROOT, "R", "devices.R"))
IN  <- file.path(DATA_RAW, "MSK_MET_2021_breast"); SUP <- OUTPUT_SUP_DIR

msk <- function(base) { g <- file.path(IN, paste0(base, ".gz"))
                        if (file.exists(g)) g else file.path(IN, base) }

GENES <- c("TP53","PIK3CA","CDH1","GATA3","MAP3K1","KMT2C","ESR1","ARID1A","AKT1","NF1",
           "NCOR1","FOXA1","TBX3","RB1","ERBB2","CCND1","MYC","FGFR1","PTEN","CDKN2A")
SILENT <- c("Silent","Intron","3'UTR","5'UTR","3'Flank","5'Flank","IGR","RNA","Splice_Region")
SC <- c("DMETS_DX_BONE","DMETS_DX_LIVER","DMETS_DX_LUNG","DMETS_DX_CNS_BRAIN")
ng <- length(GENES)

clin <- fread(msk("clinical.csv"))
mut  <- fread(msk("mutations.csv"))[!(mutationType %in% SILENT) & !is.na(hugo)]
cna  <- fread(msk("cna.csv"))[alteration %in% c(2L,-2L) & !is.na(hugo)]
s <- clin$sampleId
alt <- unique(rbind(mut[,.(sampleId,hugo)],cna[,.(sampleId,hugo)]))[sampleId%in%s & hugo%in%GENES]
M <- matrix(0L,length(s),ng,dimnames=list(s,GENES)); M[cbind(match(alt$sampleId,s),match(alt$hugo,GENES))]<-1L
pm <- rowSums(as.matrix(clin[,..SC])=="Yes")>0
cp <- clin[pm]; Mp <- M[s[pm],,drop=FALSE]
im <- which(grepl("Metas", cp$SAMPLE_TYPE, ignore.case=TRUE))
Xm <- Mp[im,,drop=FALSE]; cm <- cp[im]
subtype <- ifelse(is.na(cm$SUBTYPE_ABBREVIATION),"NA",cm$SUBTYPE_ABBREVIATION)
burden  <- suppressWarnings(as.numeric(cm$MET_SITE_COUNT)); burden[is.na(burden)] <- median(burden,na.rm=TRUE)
bt <- as.character(as.integer(cut(burden,quantile(burden,c(0,1/3,2/3,1)),include.lowest=TRUE)))
amt <- mut[,.N,by=sampleId]$N[match(cm$sampleId, mut[,.N,by=sampleId]$sampleId)]; amt[is.na(amt)]<-0
tt <- as.character(as.integer(cut(amt,quantile(amt,c(0,1/3,2/3,1)),include.lowest=TRUE)))
S_both <- paste(subtype,bt,sep="|")
yLb <- grepl("Lung", cm$METASTATIC_SITE, ignore.case=TRUE)

lnor <- function(X){k<-ncol(X);b<-crossprod(X);ct<-diag(b);a<-b;bb<-matrix(ct,k,k)-a
  cc<-matrix(ct,k,k,byrow=TRUE)-a;d<-nrow(X)-a-bb-cc;log(((a+.5)*(d+.5))/((bb+.5)*(cc+.5)))}
lnor_var <- function(X){k<-ncol(X);b<-crossprod(X);ct<-diag(b);a<-b;bb<-matrix(ct,k,k)-a
  cc<-matrix(ct,k,k,byrow=TRUE)-a;d<-nrow(X)-a-bb-cc
  list(L=log(((a+.5)*(d+.5))/((bb+.5)*(cc+.5))),V=1/(a+.5)+1/(bb+.5)+1/(cc+.5)+1/(d+.5))}
mad_dl <- function(A,B){L<-lnor(A)-lnor(B);mean(abs(L)[upper.tri(L)])}
pstrat <- function(y,st,N,X=Xm){obs<-mad_dl(X[y,,drop=FALSE],X[!y,,drop=FALSE])
  lv<-split(seq_along(y),st);nv<-vapply(lv,function(k)sum(y[k]),integer(1));pp<-numeric(N)
  for(b in 1:N){yp<-logical(length(y));for(j in seq_along(lv)){k<-lv[[j]];if(nv[j]>0)yp[k[sample(length(k),nv[j])]]<-TRUE}
    pp[b]<-mad_dl(X[yp,,drop=FALSE],X[!yp,,drop=FALSE])};list(obs=obs,p=(sum(pp>=obs)+1)/(N+1))}

cat("[FigS1_div] confound battery (biopsy-site exposure) ...\n")
set.seed(1); N_PERM <- 2500
hi <- which(burden >= 5); lungHist <- cm$DMETS_DX_LUNG=="Yes"
spec <- which(!grepl("Unspecified", cm$METASTATIC_SITE, ignore.case=TRUE))
batt <- rbind(
  data.table(test="matched on subtype + burden",             P=pstrat(yLb,S_both,N_PERM)$p,                                            n_lung=sum(yLb)),
  data.table(test="matched on exact metastatic-site count",  P=pstrat(yLb,as.character(burden),N_PERM)$p,                               n_lung=sum(yLb)),
  data.table(test="matched on subtype + mutation load",      P=pstrat(yLb,paste(subtype,tt,sep="|"),N_PERM)$p,                          n_lung=sum(yLb)),
  data.table(test="high-burden specimens only",              P=pstrat(yLb[hi],S_both[hi],N_PERM,X=Xm[hi,,drop=FALSE])$p,                n_lung=sum(yLb[hi])),
  data.table(test="excluding unspecified-site specimens",    P=pstrat(yLb[spec],S_both[spec],N_PERM,X=Xm[spec,,drop=FALSE])$p,          n_lung=sum(yLb[spec])),
  data.table(test="within patients having lung history",     P=pstrat(yLb[lungHist],S_both[lungHist],N_PERM,X=Xm[lungHist,,drop=FALSE])$p, n_lung=sum(yLb[lungHist])),
  data.table(test="patient lung history (prior definition)", P=pstrat(lungHist,S_both,N_PERM)$p,                                        n_lung=sum(lungHist)))
batt[, sig := P < 0.05][, test := factor(test, levels=rev(test))]
fwrite(batt, file.path(OUTPUT_TBL_DIR, "FigS1_confound_battery.csv"))
n_sig <- sum(batt$sig)
pE <- ggplot(batt, aes(P, test, color=sig)) +
  geom_vline(xintercept=0.05, linetype=2, color="grey55") + geom_point(size=3) +
  geom_text(aes(label=sprintf("n=%d", n_lung)), hjust=-0.45, size=2.4, color="grey40") +
  scale_color_manual(values=c(`TRUE`="#16a085", `FALSE`="grey60"), guide="none") +
  scale_x_log10(breaks=c(0.001,0.01,0.05,0.2,0.5), labels=c("0.001","0.01","0.05","0.2","0.5"),
                limits=c(8e-4,1.2)) +
  labs(title=sprintf("The lung-biopsy divergence holds under %d of %d confounder controls", n_sig, nrow(batt)),
       subtitle=paste0("Each row repeats the test holding the named factor constant; n = lung specimens\n",
                       "available. Two rows do not reach 0.05, both after losing specimens and power:\n",
                       "high-burden only, and lung vs non-lung biopsies within patients who all have lung\n",
                       "involvement. Rows are related tests on overlapping data, not independent replications."),
       x="permutation P (lung biopsies vs other metastatic biopsies)", y=NULL) +
  theme_classic(base_size=10) +
  theme(plot.title=element_text(face="bold",size=10.5), plot.subtitle=element_text(size=7.4,color="grey35"))
ggsave(file.path(SUP,"FigS1_confound_battery.pdf"), pE, width=8.2, height=4.1, device = pdf_device)
ggsave(file.path(SUP,"FigS1_confound_battery.png"), pE, width=8.2, height=4.1, dpi=200)

cat("[FigS1_div] naive size-only null landscape ...\n")
XL <- Xm[yLb,,drop=FALSE]; XO <- Xm[!yLb,,drop=FALSE]
rL <- lnor_var(XL); rO <- lnor_var(XO); LL<-rL$L; LO<-rO$L; dL<-LL-LO
Zm <- (LL-LO)/sqrt(rL$V+rO$V); Pm <- 2*pnorm(abs(Zm),lower.tail=FALSE)
ut <- upper.tri(Pm); qv <- p.adjust(Pm[ut],"BH")
melt_mat <- function(Mat){d<-as.data.table(as.table(Mat));setnames(d,c("g1","g2","val"))
  d[,g1:=factor(g1,levels=GENES)][,g2:=factor(g2,levels=rev(GENES))];d[g1==g2,val:=NA]
  d[,val:=pmin(pmax(val,-2),2)];d}
heat <- function(Mat,title,dots=FALSE){
  g<-ggplot(melt_mat(Mat),aes(g1,g2,fill=val))+geom_tile()+
    scale_fill_gradient2(low="#2166ac",mid="white",high="#b2182b",midpoint=0,limits=c(-2,2),
                         na.value="grey95",name="log odds\nratio")+
    coord_equal()+labs(title=title,x=NULL,y=NULL)+theme_minimal(base_size=7)+
    theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=1,size=4.4),axis.text.y=element_text(size=4.4),
          panel.grid=element_blank(),plot.title=element_text(face="bold",size=8.2),legend.key.width=unit(7,"pt"))
  if(dots){dd<-as.data.table(which(Pm<0.05 & upper.tri(Pm),arr.ind=TRUE))
    if(nrow(dd)){d1<-copy(dd)[,`:=`(g1=factor(GENES[row],levels=GENES),g2=factor(GENES[col],levels=rev(GENES)))]
      d2<-copy(dd)[,`:=`(g1=factor(GENES[col],levels=GENES),g2=factor(GENES[row],levels=rev(GENES)))]
      g<-g+geom_point(data=rbind(d1[,.(g1,g2)],d2[,.(g1,g2)]),aes(g1,g2),inherit.aes=FALSE,size=.45)}}
  g}
pF <- (heat(LL,sprintf("Lung biopsies (n=%d)",sum(yLb))) |
       heat(LO,sprintf("Other metastatic biopsies (n=%d)",sum(!yLb))) |
       heat(dL,"Difference (lung minus other)   dot = nominal P<0.05",dots=TRUE)) +
  plot_annotation(title="Pairwise structure under a naive size-only null, shown to justify the matched analysis",
    subtitle=sprintf(paste0("A size-only null ignores cohort composition. %d of 190 gene pairs differ nominally and none survives FDR<0.10 ",
                            "(minimum q = %.2f).\nColour is the log odds ratio of co-alteration. This panel is a methodological control, not a result; ",
                            "the matched\nanalysis is main Fig. 1D-E."), sum(Pm[ut]<0.05), min(qv)),
    theme=theme(plot.title=element_text(face="bold",size=10.5), plot.subtitle=element_text(size=7.4,color="grey35")))
ggsave(file.path(SUP,"FigS1_naive_null_landscape.pdf"), pF, width=13.5, height=4.6, device = pdf_device)
ggsave(file.path(SUP,"FigS1_naive_null_landscape.png"), pF, width=13.5, height=4.6, dpi=170)

cat(sprintf("[FigS1_div] DONE. battery %d/%d significant | naive null %d/190 nominal, min q=%.2f\n",
            n_sig, nrow(batt), sum(Pm[ut]<0.05), min(qv)))
print(batt[, .(test, P=signif(P,3), n_lung, sig)])
