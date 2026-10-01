# Validation

The numerical and rendering comparisons below were performed on 2026-09-24 with R 4.5.1 on macOS arm64 using the installed package library.

On 2026-10-01, all five full figure workflows and Figure 3's light workflow completed from an isolated source snapshot using the installed R 4.5.1 library. The full Figure 3 run required restoring the explicit serial BiocParallel backend in two UCell calls; `ncores=1` alone selected a backend that attempted to open worker sockets. Scoring formulas and inputs were unchanged. This repeat check establishes successful execution of the corrected workflows; the earlier comprehensive numerical and rendering comparisons were not repeated. Fresh dependency restoration remains unverified, as detailed below.

## Execution and numerical checks

All five selected panel workflows completed in isolated test directories. Figures 1, 2, 4 and 5 used the included inputs. Figure 3 used the two separately supplied processed RDS files with matching checksums; its light workflow also completed. UCell scoring uses an explicit serial BiocParallel backend for portability.

Comparisons against existing source-workflow outputs found 82 CSV/TSV tables identical after line-ending normalization. The Figure 5 combined interaction table reproduced point estimates and classifications. Most text differences in that table were floating-point serialization differences below 5 × 10⁻¹⁶. One upper confidence bound differed: `Brca2*Hgf`, `ddGI_ci_hi`, 2.4768900000000005 in the stored output versus 2.476987500000001 in the fresh run (absolute difference 0.0000975). The current analysis code and input data were retained. This bound difference remains recorded rather than replaced with the stored value.

The Figure 2 panel-size manifest contains only the selected workflow's panels; all 27 generated rows match their corresponding source rows.

## Rendering checks

At 72 dpi, 94 of 97 generated PDFs matched their stored source counterparts pixel for pixel. Three exceptions were inspected:

- Figure 1J: the current source script contains an updated footer. The heatmap data are unchanged.
- Supplementary Figure 3G: minor label/connector placement differences from `geom_text_repel`; extracted text and numerical tables match.
- Figure 5 targeted burden plot: horizontal point jitter differs; counts, summary statistics and extracted text match.

These checks concern individual generated plots. Final manuscript page composition includes external artwork and is not regenerated here.

## Repository checks and limits

On 2026-10-01, the repository validator passed for 229 tracked files, 61 checksummed inputs and all 99 main/supplementary panel labels, and all 126 R source files parsed. The validator checks script paths, local documentation links, file sizes, nested histories and common credential/machine-path patterns, including workbook XML and compressed inputs. Earlier negative tests verified detection of altered inputs, missing scripts and plain-text/hidden-workbook machine paths.

Worksheet and shared-string XML in the two included workbooks are unchanged from their source files; only machine-specific metadata was sanitized. Staged input bytes match the distribution's SHA-256 manifest.

Panel 1J now reads `data/derived/dawnrank_panel_subset.csv` rather than the published supplementary workbook, which is not redistributed here. The extract was checked against the workbook in two ways: `scripts/derive_dawnrank_subset.R --verify` rebuilds it and compares, and a run of the panel from the extract reproduces `Fig1G_DawnRank_correlation_matrix.csv` byte for byte against the stored source-workflow output. The input matrix and the correlation matrix are identical under `identical()`.

Restoration was tested from a clean checkout with an isolated package cache. `tools/restore.R` previously activated each module after restoring instead of before, and did not load the project, so renv compared the lockfile against the caller's library; on a machine that already held the packages it reported "already synchronized" and installed nothing. The script now activates and loads each module in a subprocess rooted at that module, and a clean checkout downloads and installs as expected. Note that `renv/activate.R` is not distributed, so each module's `.Rprofile` only selects the pinned library after that module has been restored once.

A complete restoration of all five lockfiles has not finished on this machine: three packages that build from source (`RcppArmadillo`, `graphlayouts`, `igraph`) previously failed to link because the host's R reports `FLIBS` paths under `/opt/gfortran` that do not exist. Those configured directories were still absent on 2026-10-01. This toolchain gap prevents validation of a complete fresh installation on this host. Compatible prebuilt binaries can avoid compilation where available. The installed libraries also differ from some lockfiles, so execution with those libraries does not verify restoration of every pinned version.

The [GitHub Actions run for commit ed05284](https://github.com/S-Zhang-Lab/Epistasis-manuscript-code/actions/runs/36769389229) passed. That workflow checks structure, input integrity and syntax only; full figure execution and dependency restoration on GitHub Actions have not been tested. Some installed packages report that they were built under R 4.5.2. These checks establish the reported computational scope and do not constitute a separate statistical or experimental-design review. External-data access and non-scripted panels are described in [DATA.md](DATA.md) and [PANELS.md](PANELS.md).
