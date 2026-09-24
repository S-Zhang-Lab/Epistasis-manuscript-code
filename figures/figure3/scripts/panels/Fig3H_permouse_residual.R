if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(Seurat); library(data.table); library(ggplot2) })
FS6 <- 6
th6 <- fig_theme + theme(
  text=element_text(size=FS6), axis.text=element_text(size=FS6), axis.title=element_text(size=FS6),
  legend.text=element_text(size=FS6), legend.title=element_text(size=FS6),
  strip.text=element_text(size=FS6, face="bold"), plot.title=element_text(size=FS6+0.5, face="bold"),
  plot.subtitle=element_text(size=FS6-0.5, colour="grey35"))

au_raw <- readRDS(cache("aucell_hallmark_invivo.rds"))
cat("cache class:", class(au_raw)[1], "\n")
au <- if (is.list(au_raw) && !is.data.frame(au_raw)) {
  ix <- which(sapply(au_raw, function(x) (is.matrix(x)||inherits(x,"Matrix")) && length(dim(x))==2)); as.matrix(au_raw[[ix[1]]])
} else as.matrix(au_raw)
if (nrow(au) > ncol(au)) au <- t(au)
cat("AUCell (path x cell):", paste(dim(au), collapse=" x "), " | eg rows:", paste(head(rownames(au),2),collapse=", "), "\n")

o <- readRDS(raw("Perturb_InVivo_merged.rds"))
md <- data.table(cell=colnames(o), pert=as.character(o$perturbation),
                 lane=as.character(o$lane),
                 mouse=paste(as.character(o$lane), as.character(o$hash.ID), sep="|"))
rm(o); gc(verbose=FALSE)
cat("\nlane x mouse:\n"); print(table(md$lane, md$mouse))

common <- intersect(colnames(au), md$cell)
cat("\ncells AUCell:", ncol(au), " meta:", nrow(md), " common:", length(common), "\n")
au <- au[, common, drop=FALSE]; md <- md[match(common, cell)]

niche_up <- c("INTERFERON_GAMMA_RESPONSE","INTERFERON_ALPHA_RESPONSE","INFLAMMATORY_RESPONSE",
              "IL6_JAK_STAT3_SIGNALING","TNFA_SIGNALING_VIA_NFKB","EPITHELIAL_MESENCHYMAL_TRANSITION","ANGIOGENESIS")
niche_dn <- c("MYC_TARGETS_V1","E2F_TARGETS","G2M_CHECKPOINT","MTORC1_SIGNALING","OXIDATIVE_PHOSPHORYLATION","GLYCOLYSIS")
niche <- c(niche_up, niche_dn); rn <- paste0("HALLMARK_", niche)
stopifnot(all(rn %in% rownames(au)))
disp <- c(INTERFERON_GAMMA_RESPONSE="IFN-gamma", INTERFERON_ALPHA_RESPONSE="IFN-alpha", INFLAMMATORY_RESPONSE="Inflammatory",
          IL6_JAK_STAT3_SIGNALING="IL6-JAK-STAT3", TNFA_SIGNALING_VIA_NFKB="TNFa-NFkB", EPITHELIAL_MESENCHYMAL_TRANSITION="EMT",
          ANGIOGENESIS="Angiogenesis", MYC_TARGETS_V1="MYC targets", E2F_TARGETS="E2F targets", G2M_CHECKPOINT="G2M checkpoint",
          MTORC1_SIGNALING="mTORC1", OXIDATIVE_PHOSPHORYLATION="OxPhos", GLYCOLYSIS="Glycolysis")

S <- as.data.table(t(au[rn, , drop=FALSE])); setnames(S, niche)
S[, `:=`(pert=md$pert, lane=md$lane, mouse=md$mouse)]
long <- melt(S, id.vars=c("pert","lane","mouse"), variable.name="pathway", value.name="score")
agg <- long[, .(mean_score=mean(score), n=.N), by=.(mouse, lane, pert, pathway)]

DKO <- c(Cdh1="Pten*Cdh1", Cx3cl1="Pten*Cx3cl1", Cxcr5="Pten*Cxcr5", Tlr7="Pten*Tlr7"); MINCELL <- 20
base <- agg[pert=="NT*NT", .(mouse, pathway, base=mean_score, n_base=n)]
dev  <- merge(agg[pert %chin% DKO], base, by=c("mouse","pathway"))
dev[, partner := names(DKO)[match(pert, DKO)]]
dev[, eclass  := ifelse(partner %chin% c("Cdh1","Cx3cl1"),"Synergistic","Buffering")]
dev[, dev := mean_score - base]
cat("\ncells-per (mouse x DKO), min:", min(dev$n), " NT*NT min:", min(dev$n_base), "\n")
dev <- dev[n>=MINCELL & n_base>=MINCELL]
cat("mouse-genotype observations kept:", nrow(dev), " (of", length(niche), "paths x", nrow(unique(dev[,.(mouse,partner)])), "mouse-genotypes)\n")

summ_p <- dev[, .(m=mean(dev), se=sd(dev)/sqrt(.N), nmice=.N), by=.(pathway, partner, eclass)]
summ_e <- dev[, .(m=mean(dev), se=sd(dev)/sqrt(.N), nmice=.N,
                  sign_consistent=(all(dev>0)|all(dev<0))), by=.(pathway, eclass)]
fwrite(dev,   cache("P2_permouse_deviations.csv"))
fwrite(summ_e, tbl("P2_permouse_summary_by_eclass.csv"))
cat("\n--- syn - buf (mouse-mean dev) per pathway ---\n")
wide <- dcast(summ_e, pathway ~ eclass, value.var="m"); wide[, syn_minus_buf := Synergistic - Buffering]
print(wide[order(-abs(syn_minus_buf))][, lapply(.SD, function(x) if(is.numeric(x)) round(x,4) else x)])

ford <- rev(disp[niche])
summ_p[, plab := factor(disp[as.character(pathway)], levels=ford)]
summ_e[, plab := factor(disp[as.character(pathway)], levels=ford)]
summ_p[, partner := factor(partner, levels=c("Cdh1","Cx3cl1","Cxcr5","Tlr7"))]
summ_p[, eclass := factor(eclass, levels=c("Synergistic","Buffering"))]
summ_e[, eclass := factor(eclass, levels=c("Synergistic","Buffering"))]

lim <- max(abs(summ_p$m))
pA <- ggplot(summ_p, aes(partner, plab)) +
  geom_point(aes(fill=m, size=abs(m)), shape=21, colour="grey40", stroke=0.2) +
  facet_grid(. ~ eclass, scales="free_x", space="free_x") +
  scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b", midpoint=0, limits=c(-lim,lim), name="DKO - own-mouse\nNT*NT (AUCell)") +
  scale_size(range=c(0.6,3.6), guide="none") +
  labs(x=NULL, y=NULL, title="P2  Niche program residual (per-mouse)",
       subtitle="each DKO minus its own mouse's NT*NT; mean over mice") +
  th6 + theme(panel.grid.major=element_line(colour="grey92", linewidth=0.2), axis.line=element_blank(),
              axis.ticks=element_blank(), strip.background=element_rect(fill="grey92", colour=NA),
              legend.position="right", legend.key.height=unit(14,"pt"))
ggsave(fig_main("Fig3H_permouse_residual_dotplot.pdf"), pA, device="pdf", width=3.1, height=2.25, units="in", useDingbats=FALSE)
ggsave(fig_main("Fig3H_permouse_residual_dotplot.png"), pA, width=3.1, height=2.25, units="in", dpi=300, bg="white")

cat("\nDONE. wrote Fig3H_permouse_residual_dotplot (pdf/png) to", fig_main(),
    "and P2_permouse_deviations.csv (per-mouse deviations) to", cache(), "\n")
