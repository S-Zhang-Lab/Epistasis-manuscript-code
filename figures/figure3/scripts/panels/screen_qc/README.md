# Screen-QC panels — retained for provenance, NOT part of Figure 3

These six scripts reproduce quality-control and enrichment views of the in vivo
combinatorial screen from the same committed tables the Figure 3 panels use
(`data/derived/screens/`). They were consolidated here when Xiyu Liu's updated
screen data was incorporated (2026-08).

**None of them is a Figure 3 panel, and `scripts/make_figure_3.R` does not run
them.** Status (set by S. Zhang, 2026-08-28):

| Script | Produces | Status |
|---|---|---|
| `01_plasmid_library_QC.R` | plasmid library read-count distribution, Lorenz/Gini | already in the **Figure 2** supplement |
| `02_day0_cell_library_QC.R` | Day-0 library evenness, replicate reproducibility | already in the **Figure 2** supplement |
| `03_day0_top_bottom_pairs.R` | 10 most / 10 least abundant pairs at baseline | already in the **Figure 2** supplement |
| `04_invivo_enrichment_summary.R` | enrichment-category bar + ridgeline (Day 3/18) | already in the **Figure 2** supplement |
| `05_representative_pairs.R` | LFC × FDR dot plot + significance trajectory | **dropped** from Supplementary Fig 3 (low information) |
| `06_pten_cdh1_takeover.R` | Pten\*Cdh1 rank trajectory + pool share | **dropped** from Supplementary Fig 3 (superseded by main Fig 3K) |

Run any of them individually if needed:

```bash
Rscript scripts/panels/screen_qc/06_pten_cdh1_takeover.R
```

They are tidyverse-first and must run in a clean R session: sourcing them into
the umbrella driver's long-lived environment lets `matrixStats::count` mask
`dplyr::count`.

Note the endpoint is **Day 18** throughout (the count file's `D14_*` columns were
relabelled `D18_*` when the data were incorporated).
