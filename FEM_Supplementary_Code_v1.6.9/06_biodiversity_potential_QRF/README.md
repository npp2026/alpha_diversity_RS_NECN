# QRF potential and static q95 export

Run `run_potential_workflow.R` after setting `FEM_POTENTIAL_DATA_DIR`. Use `--help` for options and `--check` for input/dependency preflight. Full commands are in [the run guide](../README.md#run-guide).

## Input contract

The directory contains `train4pot.csv`, `ENV.tif` and `templ_1km.tif`.

- CSV columns: `long`, `lat`, `Forest_age`, `Rich_tree`, `Shannon_wiener`, `DEMc`, `bio6_wc`, `bio10_wc`, `bio12_wc`, `bio17_wc`, `bio15_wc`, `bio4_wc`.
- Raster bands: `DEM`, `BIO6`, `BIO10`, `BIO12`, `BIO17`, `BIO15`, `BIO4`.
- Template: one layer, projected metre CRS, 1000 × 1000 m; only finite nonzero cells are in the domain.

## Execution

Default steps are `predict,oldage`, followed by logged static export in manuscript mode. `--steps=tune,predict,oldage` adds tuning diagnostics. Scripts share finite/sentinel cleaning, content-based cache fingerprints and fresh-output checks through `R/qrf_contracts.R`. `FEM_QRF_THREADS` sets model threads.

The manuscript profile uses 1,500 trees: richness mtry 4/node 10/forced age; Shannon mtry 2/node 10/no forced age. Surfaces retain ages 80/100/120 and q90/q95. Absolute age holdouts are 100/120/148/152; within-old-age validation uses age >=100 for both training and testing, ten folds × three repeats. Small groups can be explicitly skipped by minimum-size guards.

`FEM_QRF_PARAMETER_MODE=retuned` uses the tuning snippet and prevents automatic manuscript static export. Old-age validation and static export require a completed prediction run, matching input hashes and matching parameters. Reusing nonempty requested output directories is rejected.

`static_q95_age100/` contains two positive age-100 q95 rasters and a provenance manifest for module 05. Rebuild current Figure S8/Table S2 from the saved diagnostics and model outputs. Methodological qualifications are in [IMPLEMENTATION_CHOICES.md](../IMPLEMENTATION_CHOICES.md).
