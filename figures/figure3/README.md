# Figure 3

This module contains the analysis code and small inputs for Figure 3 and its supplementary figure.

From the repository root:

```bash
Rscript tools/restore.R 3
python3 tools/reproduce.py --figures 3
```

See the [panel guide](../../docs/PANELS.md) for the manuscript panel mapping, the [data guide](../../docs/DATA.md) for input provenance, and [validation notes](../../docs/VALIDATION.md) for the tested scope.

`scripts/run.R` records execution order. Local helper functions are in `R/`; generated plots and tables are written under `output/`. Run individual R scripts from this module directory. The module has an independent `renv.lock` and `.here` marker.

The full workflow requires two external processed RDS files. Use `python3 tools/reproduce.py --figures 3 --light` from the root for the three panels that use included inputs.
