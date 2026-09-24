suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "helpers.R"))
  library(dplyr); library(tidyr); library(readr)
})
set.seed(1)

WT_DIFF_CUT <- 2.0

DRV_CLASS <- c(Apc="Prolif", Pten="Prolif", Kit="Prolif", Tsc2="Prolif",
               Stk11="Prolif", Rb1="Prolif", Brca1="DNArep", Brca2="DNArep",
               Chek2="DNArep", Smad4="TGFb", Tgfbr2="TGFb", Krt15="EMT",
               Serpinb5="EMT", Sox10="EMT", Calml3="Other", Chit1="Other", Clca2="Other")
PRT_CLASS <- c(Btla="Immune", Cx3cl1="Immune", Cxcr2="Immune", Cxcr5="Immune",
               Fcrl6="Immune", Il6="Immune", Mpeg1="Immune", Mysm1="Immune",
               Pdcd1lg2="Immune", Tlr7="Immune", Tnfsf15="Immune", Cdh1="Adhesion",
               Mfap5="ECM", Muc16="ECM", Pla2g2a="ECM", Pla2r1="ECM", Hgf="GF")

mk <- function(df, suf) df %>% transmute(gene, type,
        !!paste0("WT_",  suf) := mean_top_LFC_WT_vs_Cell_corrected,
        !!paste0("Dep_", suf) := mean_top_LFC_Depleted_vs_Cell_corrected)
glc <- read_csv(tbl("geneLFC_Cx3cr1.csv"), show_col_types = FALSE)
glm <- read_csv(tbl("geneLFC_Mrc1.csv"),   show_col_types = FALSE)

merged <- inner_join(mk(glc, "Cx"), mk(glm, "Mrc") %>% select(-type), by = "gene") %>%
  mutate(sign_flip   = (WT_Cx * WT_Mrc) < 0,
         abs_wt_diff = abs(WT_Mrc - WT_Cx),
         dLFC_Cx = Dep_Cx - WT_Cx,
         dLFC_Mrc= Dep_Mrc - WT_Mrc,
         diff_niche_dep = dLFC_Mrc - dLFC_Cx)
merged_f <- merged %>% filter(!(sign_flip & abs_wt_diff > WT_DIFF_CUT))

dual <- merged_f %>% filter(type == "dual") %>%
  separate(gene, into = c("A", "B"), sep = "\\*", remove = FALSE) %>%
  mutate(dclass = DRV_CLASS[A], pclass = PRT_CLASS[B])
n_dual_all <- sum(merged$type == "dual")
cat(sprintf("\n[filter] first-order dual pairs: %d of %d survive (%d removed by WT filter)\n",
            nrow(dual), n_dual_all, n_dual_all - nrow(dual)))

anchor_layer <- function(d, val, layer) {
  bind_rows(
    d %>% group_by(gene = A) %>% summarize(pos = "driver", n = n(),
          mean = mean(.data[[val]]),
          p = tryCatch(t.test(.data[[val]])$p.value, error = function(e) NA_real_),
          .groups = "drop"),
    d %>% group_by(gene = B) %>% summarize(pos = "partner", n = n(),
          mean = mean(.data[[val]]),
          p = tryCatch(t.test(.data[[val]])$p.value, error = function(e) NA_real_),
          .groups = "drop")) %>%
    mutate(layer = layer, q = p.adjust(p, "BH"))
}

a1 <- anchor_layer(dual, "diff_niche_dep", "first_order")

perm_max <- replicate(2000, {
  iv <- sample(dual$diff_niche_dep)
  mA <- tapply(iv, dual$A, mean); mB <- tapply(iv, dual$B, mean)
  max(abs(c(mA, mB)))
})
obs_max <- max(abs(a1$mean))

dd <- read_csv(tbl("ddGI_combined.csv"), show_col_types = FALSE) %>%
  separate(gene, into = c("A", "B"), sep = "\\*", remove = FALSE)
a2 <- anchor_layer(dd, "ddGI", "second_order")

rpt <- function(a, lab, dirpos, dirneg) {
  a <- arrange(a, mean)
  cat(sprintf("\n===== %s =====\n", lab))
  cat(sprintf("  significant anchors (BH q<0.05): %d\n", sum(a$q < 0.05, na.rm = TRUE)))
  cat(sprintf("  strongest %s:\n", dirneg)); print(as.data.frame(head(a %>% transmute(gene,pos,n,mean=round(mean,2),q=round(q,3)), 5)))
  cat(sprintf("  strongest %s:\n", dirpos)); print(as.data.frame(tail(a %>% transmute(gene,pos,n,mean=round(mean,2),q=round(q,3)), 5)))
}
rpt(a1, "FIRST-ORDER  diff_niche_dep (+ Mrc1 / - Cx3cr1), filtered", "Mrc1-biased", "Cx3cr1-biased")
cat(sprintf("  permutation: top anchor |mean|=%.2f ; P(>= by chance)=%.3f\n", obs_max, mean(perm_max >= obs_max)))
rpt(a2, "SECOND-ORDER ddGI (+ Cx3cr1 / - Mrc1), unfiltered 289", "stronger in Cx3cr1", "stronger in Mrc1")

j <- inner_join(dual %>% select(gene, diff_niche_dep),
                dd %>% transmute(gene, ddGI_mrc = -ddGI), by = "gene")
cat(sprintf("\n[concordance] first- vs second-order, sign-aligned to Mrc1+: r=%.3f (n=%d)\n",
            cor(j$diff_niche_dep, j$ddGI_mrc), nrow(j)))

cat("\n[class] first-order mean diff_niche_dep by driver / partner class:\n")
print(dual %>% group_by(dclass) %>% summarize(n=n(), meanInt=round(mean(diff_niche_dep),2), .groups="drop") %>% arrange(desc(meanInt)) %>% as.data.frame())
print(dual %>% group_by(pclass) %>% summarize(n=n(), meanInt=round(mean(diff_niche_dep),2), .groups="drop") %>% arrange(desc(meanInt)) %>% as.data.frame())

out <- bind_rows(a1, a2) %>%
  mutate(class = ifelse(pos == "driver", DRV_CLASS[gene], PRT_CLASS[gene])) %>%
  select(layer, gene, pos, class, n, mean, p, q) %>%
  arrange(layer, mean)
write_csv(out, tbl("niche_anchor_analysis.csv"))
cat(sprintf("\nwrote %s (%d anchor rows)\n", tbl("niche_anchor_analysis.csv"), nrow(out)))

if (!interactive()) log_session("09_filtered_objective_analysis")
