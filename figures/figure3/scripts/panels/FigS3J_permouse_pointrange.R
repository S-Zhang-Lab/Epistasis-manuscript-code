if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
FS6 <- 6
th6 <- fig_theme + theme(text=element_text(size=FS6), axis.text=element_text(size=FS6), axis.title=element_text(size=FS6),
  legend.text=element_text(size=FS6), legend.title=element_text(size=FS6),
  plot.title=element_text(size=FS6+0.5, face="bold"), plot.subtitle=element_text(size=FS6-0.5, colour="grey35"))

niche <- c("INTERFERON_GAMMA_RESPONSE","INTERFERON_ALPHA_RESPONSE","INFLAMMATORY_RESPONSE","IL6_JAK_STAT3_SIGNALING",
           "TNFA_SIGNALING_VIA_NFKB","EPITHELIAL_MESENCHYMAL_TRANSITION","ANGIOGENESIS","MYC_TARGETS_V1","E2F_TARGETS",
           "G2M_CHECKPOINT","MTORC1_SIGNALING","OXIDATIVE_PHOSPHORYLATION","GLYCOLYSIS")
disp <- c(INTERFERON_GAMMA_RESPONSE="IFN-gamma", INTERFERON_ALPHA_RESPONSE="IFN-alpha", INFLAMMATORY_RESPONSE="Inflammatory",
          IL6_JAK_STAT3_SIGNALING="IL6-JAK-STAT3", TNFA_SIGNALING_VIA_NFKB="TNFa-NFkB", EPITHELIAL_MESENCHYMAL_TRANSITION="EMT",
          ANGIOGENESIS="Angiogenesis", MYC_TARGETS_V1="MYC targets", E2F_TARGETS="E2F targets", G2M_CHECKPOINT="G2M checkpoint",
          MTORC1_SIGNALING="mTORC1", OXIDATIVE_PHOSPHORYLATION="OxPhos", GLYCOLYSIS="Glycolysis")

dev <- fread(cache("P2_permouse_deviations.csv"))

mouseDT <- dev[, .(dev_m=mean(dev)), by=.(mouse, lane, eclass, pathway)]
cat("mice per class:\n"); print(unique(mouseDT[,.(mouse,eclass)])[, .N, by=eclass])

stat <- rbindlist(lapply(niche, function(p){
  s <- mouseDT[eclass=="Synergistic" & pathway==p, dev_m]
  b <- mouseDT[eclass=="Buffering"   & pathway==p, dev_m]
  tt <- tryCatch(t.test(s, b), error=function(e) list(p.value=NA_real_))
  ww <- tryCatch(suppressWarnings(wilcox.test(s, b, exact=FALSE)), error=function(e) list(p.value=NA_real_))
  data.table(pathway=p, n_syn=length(s), n_buf=length(b), mean_syn=mean(s), mean_buf=mean(b),
             diff=mean(s)-mean(b), p_t=tt$p.value, p_wilcox=ww$p.value)
}))
stat[, q_t := p.adjust(p_t, "BH")]
setorder(stat, p_t)
fwrite(stat, tbl("P2_permouse_stats.csv"))
cat("\n=== Synergistic vs Buffering, per pathway (mouse-level, n=3 vs 4) ===\n")
print(stat[, .(pathway=disp[pathway], mean_syn=round(mean_syn,4), mean_buf=round(mean_buf,4),
               diff=round(diff,4), p_t=round(p_t,3), q_FDR=round(q_t,3), p_wilcox=round(p_wilcox,3))])

summ <- mouseDT[, .(m=mean(dev_m), se=sd(dev_m)/sqrt(.N), n=.N), by=.(pathway, eclass)]
ford <- rev(disp[niche])
summ[, plab := factor(disp[pathway], levels=ford)]
summ[, eclass := factor(eclass, levels=c("Synergistic","Buffering"))]
summ <- merge(summ, stat[, .(pathway, p_t, q_t)], by="pathway")
ptxt <- unique(summ[, .(plab, p_t, q_t)])
ptxt[, lab := sprintf("p=%.2f%s", p_t, ifelse(!is.na(q_t) & q_t<0.1, "*", ""))]

xr   <- max(abs(c(summ$m + summ$se, summ$m - summ$se)), na.rm=TRUE)
xtxt <- xr * 1.75
pB <- ggplot(summ, aes(m, plab, colour=eclass)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey60") +
  geom_errorbarh(aes(xmin=m-se, xmax=m+se), height=0, linewidth=0.4, position=position_dodge(width=0.62)) +
  geom_point(size=1.5, position=position_dodge(width=0.62)) +
  geom_text(data=ptxt, aes(x=xtxt, y=plab, label=lab), inherit.aes=FALSE, hjust=1, size=FS6/.MM, colour="grey25") +
  scale_colour_manual(values=c(Synergistic="#b2182b", Buffering="#3b6fb0"), name=NULL) +
  scale_x_continuous(limits=c(-xr*1.15, xtxt*1.02)) +
  labs(x="DKO - own-mouse NT*NT  (AUCell)", y=NULL,
       title="P2  Niche program residual (per-mouse)",
       subtitle="mean +/- SE over mice (n=3 syn, 4 buf); Welch t per pathway; * FDR<0.1") +
  th6 + theme(panel.grid.major.y=element_line(colour="grey93", linewidth=0.2), legend.position="bottom")
ggsave(fig_supp("FigS3J_permouse_pointrange.pdf"), pB, device="pdf", width=3.0, height=2.55, units="in", useDingbats=FALSE)
ggsave(fig_supp("FigS3J_permouse_pointrange.png"), pB, width=3.0, height=2.55, units="in", dpi=320, bg="white")
cat("\nDONE -> FigS3J_permouse_pointrange + P2_permouse_stats.csv\n")
