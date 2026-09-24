# Validation

Validated on 2026-09-24 with R 4.5.1 on macOS arm64 using the installed package library.

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

The repository validator passes for 61 checksummed inputs and all 99 main/supplementary panel labels. It checks script paths, local documentation links, file sizes, nested histories and common credential/machine-path patterns, including workbook XML and compressed inputs. Negative tests verified detection of altered inputs, missing scripts and plain-text/hidden-workbook machine paths. All 119 R source files parse.

Worksheet and shared-string XML in the three included workbooks are unchanged from their source files; only machine-specific metadata was sanitized. Staged input bytes match the distribution's SHA-256 manifest.

A fresh restoration of all five package lockfiles and execution on GitHub Actions have not been tested. Some installed packages report that they were built under R 4.5.2. The CI workflow checks structure, input integrity and syntax only. These checks establish the reported computational scope and do not constitute a separate statistical or experimental-design review. External-data access and non-scripted panels are described in [DATA.md](DATA.md) and [PANELS.md](PANELS.md).
