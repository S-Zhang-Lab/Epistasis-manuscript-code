if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
PT <- 1/72
base5 <- theme_classic(base_size=5, base_family="Helvetica") +
  theme(text=element_text(colour="black", size=5), axis.text=element_text(colour="black", size=5),
        axis.title=element_text(size=5), axis.line=element_line(linewidth=0.25), axis.ticks=element_line(linewidth=0.25),
        legend.text=element_text(size=5), legend.title=element_blank(), legend.key.size=unit(5,"pt"),
        legend.margin=margin(0,0,0,0), legend.box.spacing=unit(1,"pt"), plot.margin=margin(1,3,1,1), plot.title=element_blank())

D <- fread(tbl("Fig3F_edist9_stats.csv"))
BAND <- D$band[1]
D[, context := factor(context, levels=c("In vitro","In vivo"))]
ord <- c("Pten*Cdh1","Pten*Cx3cl1","Pten*Cxcr5","Pten*Tlr7","NT*Cdh1","NT*Cx3cl1","NT*Cxcr5","NT*Tlr7","Pten*NT")
Dw <- dcast(D, geno+cat ~ context, value.var="edist"); setnames(Dw, c("geno","cat","vr","vo"))
ql <- dcast(D, geno ~ context, value.var="qlab"); setnames(ql, c("geno","qv","qo")); Dw <- merge(Dw, ql, by="geno")
Dw[, geno := factor(pair_lab(geno), levels=pair_lab(rev(ord)))]
xmax <- max(Dw$vr)*1.12
p9 <- ggplot(Dw) +
  annotate("rect", xmin=-Inf, xmax=BAND, ymin=-Inf, ymax=Inf, fill="grey92") +
  geom_vline(xintercept=0, linewidth=0.25, colour="grey70") +
  geom_segment(aes(x=vr, xend=vo, y=geno, yend=geno), colour="grey72", linewidth=0.4) +
  geom_point(aes(x=vo, y=geno, colour="in vivo",  shape=cat), size=1.5) +
  geom_point(aes(x=vr, y=geno, colour="in vitro", shape=cat), size=1.5) +
  geom_text(aes(x=vr, y=geno, label=qv), nudge_y= 0.34, size=4/.MM, colour="#b2182b") +
  geom_text(aes(x=vo, y=geno, label=qo), nudge_y=-0.34, size=4/.MM, colour="#5d6b7a") +
  annotate("text", x=BAND, y=9.95, label="control baseline", size=3.6/.MM, colour="grey45", hjust=1.03) +
  scale_colour_manual(values=c("in vitro"="#b2182b","in vivo"="#8aa0b6")) +
  scale_shape_manual(values=c("DKO"=16,"partner single"=17,"PTEN single"=15)) +
  scale_x_continuous(limits=c(0,xmax), expand=expansion(mult=c(0.01,0.03))) +
  scale_y_discrete(expand=expansion(add=c(0.65,0.75))) +
  labs(x="E-distance from control (energy distance)", y=NULL) +
  base5 + theme(panel.grid.major.y=element_line(colour="grey94", linewidth=0.2),
                legend.position="top", legend.justification="right", legend.spacing.x=unit(2,"pt"))
ggsave(fig_main("Fig3F_edist_dumbbell9.pdf"), p9, device="pdf", width=190*PT, height=185*PT, units="in")
ggsave(fig_main("Fig3F_edist_dumbbell9.png"), p9, width=190*PT, height=185*PT, units="in", dpi=600, bg="white")
cat("replotted Fig3F_edist_dumbbell9 (staggered q labels)\n")
