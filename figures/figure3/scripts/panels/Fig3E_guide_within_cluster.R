if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
set.seed(1234)
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(ggplot2) })

cramersV <- function(tab){
  chi <- suppressWarnings(chisq.test(tab)$statistic); n <- sum(tab)
  sqrt(as.numeric(chi)/(n*(min(dim(tab))-1)))
}
jsd <- function(p,q){ m<-(p+q)/2
  kl<-function(a,b){ idx<-a>0; sum(a[idx]*log2(a[idx]/b[idx])) }
  0.5*kl(p,m)+0.5*kl(q,m) }
mean_pair_js <- function(P){
  g <- rownames(P); v <- c()
  for (a in 1:(length(g)-1)) for (b in (a+1):length(g)) v <- c(v, jsd(P[a,],P[b,]))
  mean(v) }

meta <- function(file, poolcol, tag){
  o <- readRDS(raw(file))
  data.frame(context=tag, pool=as.character(o@meta.data[[poolcol]]),
             cluster=as.character(o$joint_clusters),
             guide=as.character(o$perturbation), stringsAsFactors=FALSE)
}
M <- rbind(meta("Perturb_InVitro_HTO.rds","pool","In vitro"),
           meta("Perturb_InVivo_merged.rds","lane","In vivo"))
M$guide <- sub("\\*","_",M$guide)
M$context <- factor(M$context, levels=c("In vitro","In vivo"))
M$pool    <- factor(M$pool,    levels=c("Synergy","Buffering"))

stat <- do.call(rbind, lapply(split(M, list(M$context,M$pool), drop=TRUE), function(d){
  tab <- table(d$guide, d$cluster)
  P <- sweep(tab,1,rowSums(tab),"/")
  data.frame(context=d$context[1], pool=d$pool[1],
             cramersV=round(cramersV(tab),3),
             mean_pairwise_JSD=round(mean_pair_js(P),4),
             n_clusters=ncol(tab), n_guides=nrow(tab))
}))
stat <- stat[order(stat$pool, stat$context),]
write.csv(stat, tbl("composition_flatness_stats.csv"), row.names=FALSE)
cat("\n=== composition flatness (lower V / JSD = flatter = more compressed) ===\n")
print(stat, row.names=FALSE)
cat("\nin-vivo / in-vitro Cramer's V ratio by pool:\n")
for (pl in levels(M$pool)){
  vv <- stat$cramersV[stat$pool==pl & stat$context=="In vitro"]
  vo <- stat$cramersV[stat$pool==pl & stat$context=="In vivo"]
  cat(sprintf("  %-9s  in vitro %.3f -> in vivo %.3f   (ratio %.2f)\n", pl, vv, vo, vo/vv))
}

prop <- as.data.frame(with(M, table(context, pool, cluster, guide)))
prop <- subset(prop, Freq>0)
prop$cluster <- factor(prop$cluster, levels=as.character(0:20))

prop <- do.call(rbind, lapply(split(prop, list(prop$context,prop$pool,prop$cluster), drop=TRUE),
  function(d){ d$frac <- d$Freq/sum(d$Freq); d }))

GORDER <- c("NT_NT","Pten_NT","NT_Cdh1","Pten_Cdh1","NT_Cx3cl1","Pten_Cx3cl1",
            "NT_Cxcr5","Pten_Cxcr5","NT_Tlr7","Pten_Tlr7")
prop$guide <- factor(prop$guide, levels=GORDER)

PAL <- c(NT_NT="#7f7f7f", Pten_NT="#f28e2b",
         NT_Cdh1="#59a14f", Pten_Cdh1="#4e79a7", NT_Cx3cl1="#edc948", Pten_Cx3cl1="#b07aa1",
         NT_Cxcr5="#e15759", Pten_Cxcr5="#76b7b2", NT_Tlr7="#9c755f", Pten_Tlr7="#ff9da7")
p <- ggplot(prop, aes(cluster, frac, fill=guide)) +
  geom_col(width=0.9, colour="white", linewidth=0.08) +
  facet_grid(pool ~ context, scales="free_x", space="free_x") +
  scale_fill_manual(values=PAL, name=NULL, breaks=GORDER) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x="Cluster", y="Proportion", title="Guide distribution within each cell cluster") +
  fig_theme +
  guides(fill=guide_legend(nrow=2, byrow=TRUE, keywidth=unit(5,"pt"), keyheight=unit(5,"pt"))) +
  theme(panel.spacing=unit(4,"pt"), strip.background=element_blank(),
        strip.text=element_text(size=FS, face="bold"),
        legend.position="bottom", legend.box.spacing=unit(2,"pt"),
        legend.text=element_text(size=FS-0.5), legend.key.size=unit(5,"pt"))
save_panel(p, "Fig3E_guide_within_cluster", 3.5, 3.1, fig_main())
cat("\npanelC guide-within-cluster written.\n")
