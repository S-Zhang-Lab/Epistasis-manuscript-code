if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(ggplot2) })
CO <- read.csv(tbl("label_epistasis_coeffs.csv"), check.names=FALSE)
G  <- CO[CO$level=="gene",]; G$context <- factor(G$context, levels=c("In vitro","In vivo"))
shp <- c(Cdh1=16, Cx3cl1=17, Cxcr5=15, Tlr7=18)

pA <- ggplot(G, aes(c1, c2, colour=context, shape=partner)) +
  geom_abline(slope=1,intercept=0,colour="grey80",linewidth=0.3,linetype="dashed") +
  geom_hline(yintercept=0,colour="grey90",linewidth=0.3)+geom_vline(xintercept=0,colour="grey90",linewidth=0.3) +
  geom_point(size=1.5, stroke=0.4) +
  ggrepel::geom_text_repel(aes(label=partner), size=FS/.MM, segment.size=0.15, max.overlaps=20,
                           min.segment.length=0, box.padding=0.2, show.legend=FALSE) +
  annotate("text", x=Inf, y=-Inf, label="in vitro", colour="#2166ac", hjust=1.05, vjust=-0.7, size=FS/.MM, fontface="bold") +
  annotate("text", x=-Inf, y=Inf, label="in vivo", colour="#d62728", hjust=-0.05, vjust=1.5, size=FS/.MM, fontface="bold") +
  scale_colour_manual(values=c("In vitro"="#2166ac","In vivo"="#d62728"), guide="none") +
  scale_shape_manual(values=shp, guide="none") +
  labs(x="c1 (PTEN)", y="c2 (partner)", title="Interaction coefficients") +
  fig_theme + theme(legend.position="none")
save_panel(pA, "FigS3G_coefficient_flip", 1.45, 1.45, fig_supp())

pB <- ggplot(G, aes(factor(partner,levels=c("Cdh1","Cx3cl1","Cxcr5","Tlr7")), R2, fill=context)) +
  geom_col(position=position_dodge(width=0.72), width=0.66, colour="black", linewidth=0.2) +
  geom_text(aes(label=sprintf("%.2f",R2)), position=position_dodge(width=0.72), vjust=-0.3, size=FS/.MM) +
  scale_fill_manual(values=c("In vitro"="#2166ac","In vivo"="#d62728"), name=NULL) +
  scale_y_continuous(limits=c(0,1), expand=expansion(mult=c(0,0.1))) +
  labs(x=NULL, y="additivity R²", title="Additivity of the double-KO (gene level)") +
  fig_theme + theme(axis.text.x=element_text(face="italic"), legend.position="top", legend.key.size=unit(7,"pt"))
save_panel(pB, "FigS3G_additivity_R2", 3.4, 2.1, fig_archive())

PR <- read.csv(tbl("label_pathway_residual.csv"), check.names=FALSE)
PR$context <- factor(PR$context, levels=c("In vitro","In vivo"))
agg <- aggregate(abs(residual)~pathway, PR, max); top <- agg$pathway[order(-agg$`abs(residual)`)][1:22]
PRt <- PR[PR$pathway %in% top,]; PRt$pathway <- factor(PRt$pathway, levels=rev(top))
PRt$partner <- factor(PRt$partner, levels=c("Cdh1","Cx3cl1","Cxcr5","Tlr7"))
lim <- max(abs(PRt$residual))
pC <- ggplot(PRt, aes(partner, pathway, fill=residual)) + geom_tile() +
  facet_wrap(~context, nrow=1) +
  scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b", midpoint=0, limits=c(-lim,lim), name="residual\n(DKO - additive)") +
  labs(x=NULL, y=NULL, title="Neomorphic pathway residual (Hallmark AUCell)") +
  fig_theme + theme(axis.text.x=element_text(angle=45,hjust=1,face="italic"), axis.text.y=element_text(size=FS-0.5),
                    axis.line=element_blank(), axis.ticks=element_blank(), legend.key.size=unit(6,"pt"))
save_panel(pC, "FigS3G_pathway_residual", 4.6, 3.4, fig_archive())
cat("saved: FigS3G_coefficient_flip (supplementary/);",
    "FigS3G_additivity_R2, FigS3G_pathway_residual (archive/, not placed)\n")
