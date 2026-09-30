# Data inputs and provenance

Included analytical inputs are enumerated with byte sizes and SHA-256 hashes in [inputs.tsv](../inputs.tsv). `python3 tools/validate.py` verifies every listed file. Data values are retained from the source figure workflows. Machine-specific workbook paths and personal edit metadata have been removed; worksheet and shared-string XML are unchanged. The input manifest also records the original source-file hashes. Paths below are relative to the relevant `figures/figureN/` module.

## Included inputs

| Figure | Files / directory | Provenance and role |
|---|---|---|
| 1 | `data/raw/TCGA_PanCan/` | Ten TCGA PanCancer Atlas 2018 mutation subsets, from cBioPortal DataHub; somatic co-occurrence analysis. |
| 1 | `data/raw/MSK_MET_2021_breast/` | Breast-cancer clinical, mutation and copy-number subsets from MSK-MET 2021 via cBioPortal DataHub. |
| 1 | `data/raw/GSE110590/` | GEO series matrices and the expression table for the paired primary/metastatic cohort. |
| 1 | `data/derived/dawnrank_panel_subset.csv` | DawnRank scores for the 27 panel genes present in Supplementary Table 5 of Siegel et al., across the eight lung-metastasis specimens; the only values Panel J reads from that source. |
| 1 | `data/derived/Sample_Metadata_Table.csv`, `docs/gene_lists/` | Curated paired-subset metadata and gene sets used by the analyses. |
| 2 | `data/raw/R1_construct_counts.csv`, `R1_plasmid_counts.csv` | Round-1 construct and plasmid-library count matrices. |
| 2 | `data/derived/R2_*` | Round-2 counts, MAGeCK summaries, supplied GI workbook and a regression input used to check the GI calculation. |
| 2 | `data/reference/` | Library gene list, cached STRING network and a reference clustered GI matrix. |
| 3 | `data/derived/screens/` | Gene-level screen summaries, construct counts, plasmid frequencies and cutting-efficiency measurements. |
| 3 | `data/derived/hallmark_pathways_v26.1.0_MM.csv` | Frozen mouse Hallmark gene-set snapshot used for single-cell scoring. |
| 4 | `data/raw/Cx3cr1_raw_counts.csv` | Construct counts for the CX3CR1-positive-cell depletion screen. |
| 4 | `data/derived/pathwayGI_*.csv` | Supplied PTEN/CX3CL1 and CHEK2/CX3CL1 pathway-interaction summaries. |
| 4 | `data/derived/network_scaffold_Cx3cr1_*.csv` | Fixed network geometry for the pathway plot. |
| 5 | `data/raw/`, `data/derived/` | CX3CR1 and MRC1 screen counts and supplied gene-level reference summaries. |
| 5 | `data/Tumor_Burden_Quantification/metastasis_counts_summary.csv` | Per-mouse counts for targeted Pten/Cdh1 validation. |
| 5 | `data/Pooled_Burden_Quantification/` | Per-mouse nodule counts for pooled-screen burden controls. |

The round-1 enrichment contrasts use the injected Day-0 pool, whereas round-1 genetic-interaction scoring uses Day 3 as its baseline. Round-2 genetic-interaction code starts from the supplied gene summaries; raw-read processing and MAGeCK execution are upstream of this repository. The file named `R2_Day18_legacy.gene_summary.txt` is retained only because the calculation explicitly uses it as a regression fixture.

Figure 1 public-source references:

- [cBioPortal DataHub](https://github.com/cBioPortal/datahub/tree/db2f8008a119f6008fba5a99102d085a863995a1): pinned source revision for the clinical genomic datasets.
- [GEO GSE110590](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110590): paired primary/metastatic expression data.
- [Siegel et al., JCI article 96153](https://www.jci.org/articles/view/96153): source of the DawnRank driver-activity scores. The published workbook is copyright American Society for Clinical Investigation and is not redistributed here; `data/derived/dawnrank_panel_subset.csv` carries only the 27 genes by 8 specimens that Panel J uses. `scripts/derive_dawnrank_subset.R` rebuilds or verifies that extract against the workbook, which its header explains how to download.
- [HGNC gene groups](https://www.genenames.org/tools/search/#!/groups): provenance categories are retained in the curated niche-engagement CSV. Bulk expression can include immune and stromal contributions.
- [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/): Hallmark gene-set annotations; their source terms apply.
- [STRING](https://string-db.org/): interaction-network reference used in supplementary Figure 2P.

Figure 1 includes scripts to rederive the cBioPortal subsets from upstream files. Those larger upstream files are optional and are not needed for the panel workflow. Figures 1 and 4 may consult the installed `msigdbr` package/cache; installing or fetching annotations can require network access.

## External inputs for Figure 3

Place these exact files in `figures/figure3/data/raw/`. They are excluded from Git. Both are deposited at Zenodo under [doi:10.5281/zenodo.21384699](https://doi.org/10.5281/zenodo.21384699), which always resolves to the current version of the record; access opens at publication. The sequencing data underlying them are in GEO under accessions GSE322808 and GSE322809.

| Filename | Bytes | MD5 |
|---|---:|---|
| `Perturb_InVivo_merged.rds` | 882515458 | `2bf87b061c0368f8d2927fe7f155dcc5` |
| `Perturb_InVitro_HTO.rds` | 516211443 | `c8a3728cf3a88d0d036c9d46762ab755` |

`--check-inputs` checks both file identities. Set `FIG3_CHECK_COLUMNS=1` to additionally load the objects and inspect required metadata and PCA reductions. Single-cell alignment, demultiplexing and object construction are upstream steps.

## Optional upstream reconstruction for Figure 4

The Figure 4 panel workflow uses the included pathway-GI tables. `scripts/recompute_pathway_GI_from_RDS.R` can recalculate them when the following processed objects are supplied in `figures/figure4/data/raw/RDS/`:

- `Perturb_InVivo_merged.rds`, the PTEN in-vivo object, corresponding to GEO accession GSE322809. Some earlier notes call this file `Perturb_UTSW37_GSE322809_InVivo_merged.rds`; it is the same object, and the Zenodo record carries it under the shorter name. Either filename works, because the script accepts both.
- `Perturb_UTSW33_GSE322808.rds`, the CHEK2/CX3CL1 object, corresponding to GEO accession GSE322808.

These optional upstream objects are not included in Git. Both are in the same Zenodo record as the Figure 3 objects, [doi:10.5281/zenodo.21384699](https://doi.org/10.5281/zenodo.21384699). This script obtains Hallmark annotations through `msigdbr` and writes the pathway CSVs in `data/derived/`; a different annotation release may change those tables. The regular panel workflow starts from the supplied, checksummed summaries.

## Third-party data and terms

Code in this repository is MIT-licensed. That license does not extend to the third-party data included here, which remains under the terms of its original source.

The TCGA PanCancer Atlas and MSK-MET files are not verbatim copies of the upstream releases. Both are derived subsets produced by the scripts in this repository: `scripts/subset_tcga_pancan.sh` cuts the ten mutation files down to the driver panel, and `scripts/etl_msk_met_breast.R` filters MSK-MET to breast specimens and to the columns the analyses read. Each script has a `--verify` mode that rebuilds its subset from the upstream release and compares. The upstream files are open access through the [cBioPortal DataHub](https://github.com/cBioPortal/datahub/tree/db2f8008a119f6008fba5a99102d085a863995a1) at the pinned revision recorded above; cite the original studies when reusing them.

GEO series matrices and expression tables for GSE110590 are redistributed from NCBI GEO. The DawnRank scores are not: the published workbook is copyright American Society for Clinical Investigation, and only the derived extract described above travels with this repository.

STRING interaction data and MSigDB Hallmark gene sets are redistributed under their respective Creative Commons terms, with attribution to those projects.

## Externally produced panels

The panel guide identifies experimental images, flow-cytometry outputs, CRISPResso composites and schematics that are not regenerated by these workflows. Raw images, FASTQs and flow-cytometry files are not bundled. A panel listed as external has no corresponding automated reconstruction claim.
