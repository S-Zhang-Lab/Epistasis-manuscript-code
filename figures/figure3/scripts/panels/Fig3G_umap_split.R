if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(Seurat); library(data.table); library(ggplot2); library(ggh4x) })
OUT <- fig_main(); dir.create(OUT, showWarnings=FALSE)
FS6 <- 6

o  <- readRDS(raw("Perturb_InVivo_merged.rds"))
um <- Embeddings(o, "umap")
df <- data.table(UMAP1=um[,1], UMAP2=um[,2],
                 cluster=as.character(o$seurat_clusters), lane=as.character(o$lane))
rm(o); gc(verbose=FALSE)

lab <- fread(cache("cluster_identity_labels_in_vivo.csv"))[context=="in_vivo"]
cmap <- setNames(lab$proposed_label, as.character(lab$cluster))
df[, state := cmap[cluster]]
df[, lane := factor(ifelse(lane=="Synergy","Synergistic","Buffering"), levels=c("Synergistic","Buffering"))]
states <- c("Proliferative","Hypoxia/glycolysis","EMT/matrix","Inflammatory/IFN","OXPHOS/metabolic","Signaling/dev","Hormone/other")
df[, state := factor(state, levels=states)]
cat("cells:", nrow(df), " | state counts:\n"); print(df[, .N, by=state][order(-N)])
cat("\nlane x state proportions (convergence check):\n"); print(round(prop.table(table(df$lane, df$state),1),3))

bg <- rbindlist(lapply(levels(df$lane), function(L){ data.table(UMAP1=df$UMAP1, UMAP2=df$UMAP2, lane=factor(L, levels=levels(df$lane))) }))

pal <- c("Proliferative"="#E8724C","Hypoxia/glycolysis"="#2E8B57","EMT/matrix"="#38A6C9",
         "Inflammatory/IFN"="#3B6FB0","OXPHOS/metabolic"="#7E5AA2","Signaling/dev"="#C98BB9","Hormone/other"="#D64A8B")

p <- ggplot() +
  geom_point(data=bg, aes(UMAP1,UMAP2), colour="grey88", size=0.28, stroke=0) +
  geom_point(data=df, aes(UMAP1,UMAP2, colour=state), size=0.40, stroke=0) +
  facet_wrap(~lane, nrow=1) +
  scale_colour_manual(values=pal, name="Cell state (in vivo)") +
  guides(colour=guide_legend(override.aes=list(size=1.7))) +
  labs(x="UMAP1", y="UMAP2") +
  fig_theme +
  theme(text=element_text(size=FS6), axis.title=element_text(size=FS6),
        axis.text=element_blank(), axis.ticks=element_blank(),
        legend.text=element_text(size=FS6), legend.title=element_text(size=FS6),
        legend.key.size=unit(7,"pt"), legend.position="right",
        strip.background=element_rect(fill="grey92", colour=NA),
        strip.text=element_text(size=FS6, face="bold"),
        panel.spacing=unit(4,"pt")) +
  force_panelsizes(rows=unit(120,"pt"), cols=unit(118,"pt"))

ggsave(file.path(OUT,"Fig3G_umap_split_invivo.pdf"), p, device="pdf", width=4.6, height=2.05, units="in", useDingbats=FALSE)
ggsave(file.path(OUT,"Fig3G_umap_split_invivo.png"), p, width=4.6, height=2.05, units="in", dpi=320, bg="white")
cat("\nDONE. Fig3G_umap_split_invivo written to", OUT, "\n")
