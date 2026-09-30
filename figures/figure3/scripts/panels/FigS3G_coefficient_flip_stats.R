#!/usr/bin/env Rscript
# =============================================================================
# Epistasis residual (Norman additive model) BY PERTURBATION LABEL (lanes pooled).
#   delta_X = pseudobulk(X) - pseudobulk(NT*NT)   [all by label, both lanes]
#   fit  delta_DKO = c1*delta_PTEN + c2*delta_partner + eps   (through origin, over genes)
#   c1,c2 = additive coefficients; R2 = additivity; eps = neomorphic (epistatic) residual.
# The residual is a difference-of-differences, so the modest lane offset cancels
# (partner arms share it; balanced baselines cancel) -> label-pooled is clean here.
# Permutation null: shuffle the 4 labels within the 2x2 (preserve sizes), RE-SELECT
# responsive genes each shuffle (avoids selection circularity) -> p for residual & DKO.
# Gene level (SCT) + Hallmark AUCell level (reuse group_cache). In vitro + in vivo.
# v7 repo overhaul 2026-07-15: relocated from XL_code/fig3_panels_new_build/label_epistasis_residual.R; was panelJ_interaction_coeffs (stats half) -> FigS3G_coefficient_flip_stats.
# =============================================================================
if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
set.seed(1234)
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(matrixStats) })
PARTNERS <- c("Cdh1","Cx3cl1","Cxcr5","Tlr7"); CAP <- 700; R <- 400; SEL <- 0.10

# through-origin 2-predictor regression on responsive genes -> list(c1,c2,R2,resid,sel)
fit_gi <- function(dPT,dPA,dDK){
  sel <- abs(dPT)>SEL | abs(dPA)>SEL | abs(dDK)>SEL
  if(sum(sel)<20) return(NULL)
  X<-cbind(dPT[sel],dPA[sel]); y<-dDK[sel]
  b<-tryCatch(solve(crossprod(X), crossprod(X,y)), error=function(e) c(NA,NA))
  fitv<-as.numeric(X%*%b); r<-y-fitv
  list(c1=b[1], c2=b[2], R2=1-sum(r^2)/sum(y^2), resid_norm=sqrt(sum(r^2)),
       dko_norm=sqrt(sum(y^2)), nsel=sum(sel), sel=which(sel), resid=r) }

analyze <- function(M, cellsets, tag, level){
  out<-list(); neom<-list()
  for(p in PARTNERS){
    g<-cellsets[[paste0("Pten*",p)]]; pa<-cellsets[[paste0("NT*",p)]]
    pt<-cellsets[["Pten*NT"]]; nt<-cellsets[["NT*NT"]]
    if(any(sapply(list(nt,pt,pa,g),length)<15)) next
    cap<-function(x) if(length(x)>CAP) sample(x,CAP) else x
    nt<-cap(nt); pt<-cap(pt); pa<-cap(pa); g<-cap(g)
    pb<-function(cells) rowMeans(M[,cells,drop=FALSE])
    obs<-fit_gi(pb(pt)-pb(nt), pb(pa)-pb(nt), pb(g)-pb(nt)); if(is.null(obs)) next
    # permutation: pool the 4 groups, shuffle labels preserving sizes
    pool<-c(nt,pt,pa,g); nk<-c(length(nt),length(pt),length(pa),length(g)); cut<-cumsum(nk)
    Sp<-M[,pool,drop=FALSE]
    nulls<-replicate(R,{ pr<-sample.int(sum(nk))
      i1<-pr[1:cut[1]]; i2<-pr[(cut[1]+1):cut[2]]; i3<-pr[(cut[2]+1):cut[3]]; i4<-pr[(cut[3]+1):cut[4]]
      m1<-rowMeans(Sp[,i1,drop=FALSE]); m2<-rowMeans(Sp[,i2,drop=FALSE]); m3<-rowMeans(Sp[,i3,drop=FALSE]); m4<-rowMeans(Sp[,i4,drop=FALSE])
      f<-fit_gi(m2-m1, m3-m1, m4-m1); if(is.null(f)) c(NA,NA) else c(f$resid_norm, f$dko_norm) })
    p_resid<-(1+sum(nulls[1,]>=obs$resid_norm,na.rm=TRUE))/(R+1)
    p_dko  <-(1+sum(nulls[2,]>=obs$dko_norm, na.rm=TRUE))/(R+1)
    out[[p]]<-data.frame(context=tag,level=level,partner=p,c1=obs$c1,c2=obs$c2,R2=obs$R2,
      resid_norm=obs$resid_norm,resid_p=p_resid,dko_norm=obs$dko_norm,dko_p=p_dko,
      n_sel=obs$nsel,n_DKO=length(g))
    # neomorphic features (top |residual|)
    feats<-rownames(M)[obs$sel]; rr<-obs$resid
    ord<-order(-abs(rr)); k<-min(30,length(rr))
    neom[[p]]<-data.frame(context=tag,level=level,partner=p,feature=feats[ord][1:k],residual=rr[ord][1:k])
  }
  list(coef=do.call(rbind,out), neom=do.call(rbind,neom))
}

COEF<-list(); NEOM<-list()
for(cc in list(c("Perturb_InVitro_HTO.rds","In vitro","invitro"), c("Perturb_InVivo_merged.rds","In vivo","invivo"))){
  cat("\n#########",cc[2],"#########\n")
  o<-readRDS(raw(cc[1])); pert<-as.character(o$perturbation); names(pert)<-colnames(o)
  cs<-split(colnames(o),pert)
  # gene level (SCT data, expressed genes)
  D<-GetAssayData(o,assay="SCT",layer="data"); det<-rowMeans(D[,unlist(cs)]>0); D<-D[det>0.10,]
  rg<-analyze(D,cs,cc[2],"gene"); COEF[[paste0("g",cc[2])]]<-rg$coef; NEOM[[cc[2]]]<-rg$neom
  # pathway level (Hallmark AUCell from group_cache)
  U<-readRDS(cache(paste0("ucell_hallmark_",cc[3],".rds")))$U
  rp<-analyze(U,cs,cc[2],"pathway"); COEF[[paste0("p",cc[2])]]<-rp$coef
  # also save per-pathway residual (all 50) for a heatmap
  pbU<-function(cells) rowMeans(U[,cells,drop=FALSE])
  for(p in PARTNERS){ g<-cs[[paste0("Pten*",p)]]; if(is.null(g)) next
    dPT<-pbU(cs[["Pten*NT"]])-pbU(cs[["NT*NT"]]); dPA<-pbU(cs[[paste0("NT*",p)]])-pbU(cs[["NT*NT"]]); dDK<-pbU(g)-pbU(cs[["NT*NT"]])
    b<-solve(crossprod(cbind(dPT,dPA)),crossprod(cbind(dPT,dPA),dDK)); res<-dDK-cbind(dPT,dPA)%*%b
    NEOM[[paste0("path_",cc[2],"_",p)]]<-data.frame(context=cc[2],partner=p,pathway=sub("HALLMARK_","",names(dDK)),residual=as.numeric(res)) }
  rm(o,D,U); gc(verbose=FALSE)
}
CO<-do.call(rbind,COEF[grep("^[gp]",names(COEF))])
write.csv(CO, tbl("label_epistasis_coeffs.csv"), row.names=FALSE)
write.csv(do.call(rbind,NEOM[c("In vitro","In vivo")]), tbl("label_neomorphic_genes.csv"), row.names=FALSE)
write.csv(do.call(rbind,NEOM[grep("^path_",names(NEOM))]), tbl("label_pathway_residual.csv"), row.names=FALSE)
cat("\n===== EPISTASIS (Norman additive model) by label — coefficients =====\n")
print(CO[,c("context","level","partner","c1","c2","R2","resid_p","dko_p")], row.names=FALSE, digits=3)
cat("\nresid_p = is the neomorphic residual above the label-shuffle null?  dko_p = is the DKO effect itself real?\n")
