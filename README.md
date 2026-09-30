# Epistasis manuscript code

Analysis and plotting code accompanying **“Niche-dependent epistasis within the cancer genome shapes metastatic fitness”** by Xiyu Liu and colleagues.

The repository contains five figure modules, their manuscript-related supplementary panels, and small analytical inputs. Each module retains its own R package lockfile. Start with the [panel guide](docs/PANELS.md) to locate an analysis and [data guide](docs/DATA.md) for input provenance and external-data requirements.

## Repository layout

```
README.md            this file
CITATION.cff         software citation metadata
LICENSE              MIT, covering the code only
inputs.tsv           every included input file with its size and SHA-256
panels.tsv           machine-readable panel-to-script map; the authority for panel assignment
docs/                DATA.md (input provenance), PANELS.md (panel guide), VALIDATION.md (tested scope)
tools/               repository-wide commands: validate.py, check_syntax.R, restore.R, reproduce.py
figures/figureN/     one self-contained module per manuscript figure
    .here, .Rprofile, renv.lock   project root marker and pinned package library
    R/                            module-local helpers, themes and palettes
    data/                         raw/, derived/ and reference/ inputs
    scripts/                      the analyses; run.R is the module entry point
    output/                       created on run, not tracked
```

Two directories are named `scripts`, and they differ. `tools/` at the root holds commands that act on the
whole repository. `figures/figureN/scripts/` holds the analyses for one figure.

Script naming is not uniform across modules, because the five modules began as five separate repositories.
Figure 3 is organised into `cache/` and `panels/` subdirectories; the others are flat and mix numbered
prefixes (`01_`, `03_`), figure-prefixed names (`Fig1_G_`, `Fig5S_S11_`) and prototype prefixes (`proto_`).
Filenames also carry earlier panel lettering in places. Resolve any panel through
[panels.tsv](panels.tsv) or the [panel guide](docs/PANELS.md) rather than by reading filenames.

## Figure modules

| Module | Analyses | Entry point |
|---|---|---|
| [Figure 1](figures/figure1) | Clinical genomic co-occurrence and paired primary/metastatic transcriptomes | `python3 tools/reproduce.py --figures 1` |
| [Figure 2](figures/figure2) | Dual-gene screening, enrichment, genetic interactions and library QC | `python3 tools/reproduce.py --figures 2` |
| [Figure 3](figures/figure3) | Perturb-seq, cell-state programs, editing efficiency and screen time courses | `python3 tools/reproduce.py --figures 3` |
| [Figure 4](figures/figure4) | Genetic interactions following CX3CR1-positive-cell depletion | `python3 tools/reproduce.py --figures 4` |
| [Figure 5](figures/figure5) | Comparison of CX3CR1-positive and MRC1-positive depletion contexts | `python3 tools/reproduce.py --figures 5` |

## Reproduction

Use R 4.5.1 and Python 3.9 or newer. Python uses only its standard library. Run these commands from the repository root. Package installation requires network access; Bioconductor packages and their system dependencies are needed for several modules. Several packages compile from source and need a working Fortran toolchain, which is the usual cause of a failed restore on macOS; configuring a binary repository such as `options(repos = c(P3M = "https://packagemanager.posit.co/cran/latest"))` avoids the compilers entirely.

```bash
# Check the included files and the panel-to-code mapping.
python3 tools/validate.py
Rscript --vanilla tools/check_syntax.R

# Restore separate R libraries for the modules you intend to run.
Rscript tools/restore.R 1 2 3 4 5

# These modules use included inputs.
python3 tools/reproduce.py --figures 1 2 4 5

# Figure 3 screen-based panels also use included inputs.
python3 tools/reproduce.py --figures 3 --light

# After obtaining the two external RDS files described in docs/DATA.md:
python3 tools/reproduce.py --figures 3 --check-inputs
python3 tools/reproduce.py --figures 3
```

Figure 3's full analysis requires about 1.4 GB of processed input objects and substantial additional memory; at least 32 GB RAM is recommended. Bootstrap and permutation analyses can take considerably longer than plot rendering. Outputs are written beneath each module's `output/` directory. Scripts retain their established output filenames; the panel guide maps these to manuscript labels.

The default runner executes every selected module. `--light` is limited to Figure 3 panels J–K and supplementary panel C; it does not run the transcriptomic analyses. To run an individual R script, first change into its figure module so `here::here()` resolves the module's `.here` file. Run upstream steps in `scripts/run.R` before downstream plotting scripts.

## Scope

The workflows start from the count matrices, gene-level summaries and processed single-cell objects specified in the data guide. They do not perform sequencing alignment, initial single-cell processing, microscopy segmentation, flow-cytometry gating, or final figure composition. Schematics and externally assembled panels are identified individually in the panel guide. Some source scripts emit additional diagnostic plots; manuscript panel assignments are defined by `panels.tsv`.

`renv.lock` files record the source environments. The automated GitHub workflow checks file integrity, panel coverage and R syntax; it does not claim to reproduce the complete figures. See [validation notes](docs/VALIDATION.md) for the tested scope.

## Citation and license

Please cite the accompanying manuscript and the software version used. Software citation metadata is provided in [CITATION.cff](CITATION.cff). Code is distributed under the [MIT license](LICENSE). Third-party datasets and annotations retain their original attribution and terms; see [DATA.md](docs/DATA.md).
