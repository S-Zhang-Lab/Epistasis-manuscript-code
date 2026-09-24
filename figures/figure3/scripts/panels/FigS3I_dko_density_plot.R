if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(ggplot2) })
D <- read.csv(cache("compression_concrete_cells.csv"), check.names=FALSE)

PARTNERS <- c("Cdh1","Cx3cl1","Cxcr5","Tlr7")
typ <- function(g,pn) ifelse(g=="NT*NT","ctrl", ifelse(g=="Pten*NT","pten", ifelse(g==paste0("NT*",pn),"partner","double")))
D$type <- mapply(typ, D$genotype, D$partner)
D$type <- factor(D$type, levels=c("ctrl","pten","partner","double"))
D$partner <- factor(D$partner, levels=PARTNERS)
D$context <- factor(D$context, levels=c("In vitro","In vivo"))
TC <- c(ctrl="#8f97c9", pten="#e69138", partner="#7cbaa8", double="#b0463a")
TLAB <- c(ctrl="NT_NT", pten="Pten_NT", partner="NT_par", double="Pten_par")

D2 <- subset(D, type %in% c("ctrl","double")); D2$type <- droplevels(D2$type)
p2 <- ggplot(D2, aes(dko_score, colour=type)) +
  geom_density(fill=NA, linewidth=0.6) +
  facet_grid(context ~ partner) +
  scale_colour_manual(values=TC[c("ctrl","double")], labels=TLAB[c("ctrl","double")], name=NULL) +
  labs(x="double-KO signature score (z)", y="cell density",
       title="The double-KO signature separates its cells in vitro, not in vivo") +
  fig_theme +
  theme(strip.background=element_blank(), strip.text=element_text(size=FS,face="bold"),
        panel.spacing=unit(3,"pt"), legend.position="top",
        axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour=guide_legend(nrow=1))
save_panel(p2, "FigS3I_dko_signature_density", 5.0, 2.5, fig_supp())

cat("wrote FigS3I_dko_signature_density\n")
