# Figure 2

This module contains the analysis code and small inputs for Figure 2 and its supplementary figure.

From the repository root:

```bash
Rscript scripts/restore.R 2
python3 scripts/reproduce.py --figures 2
```

See the [panel guide](../../docs/PANELS.md) for the manuscript panel mapping, the [data guide](../../docs/DATA.md) for input provenance, and [validation notes](../../docs/VALIDATION.md) for the tested scope.

`scripts/run.R` records execution order. Local helper functions are in `R/`; generated plots and tables are written under `output/`. Run individual R scripts from this module directory. The module has an independent `renv.lock` and `.here` marker.
