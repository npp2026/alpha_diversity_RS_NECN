# Dependencies

Use R >= 4.3 and Python >= 3.9. Exact historical package versions were not available; record versions from the final validated run.

| Component | Dependencies |
|---|---|
| R analytical workflows | `ranger`, `VSURF`, `xgboost`, `CAST`, `terra`, `sf`, `dplyr`, `tidyr`, `yaml`, `digest`, `jsonlite` |
| Optional trend and nearest-neighbour paths | `trend`, `modifiedmk`, `mutoss`, `FNN` (the package has deterministic base-R fallbacks where documented) |
| Optional trend diagnostics | `trend`, `modifiedmk`, `mutoss` |
| Python Moran workflows | `numpy`, `scipy`, `rasterio`, `tqdm` |
| Launchers / release checks | Bash; Node.js for local JavaScript syntax checks |
| Remote predictor extraction | Google Earth Engine access and the study plot FeatureCollection |

For an R installation with the required GDAL/GEOS/PROJ system libraries:

```r
install.packages(c("ranger", "VSURF", "xgboost", "CAST", "terra", "sf",
                   "dplyr", "tidyr", "yaml", "digest", "jsonlite"))

# Optional package-backed diagnostics. Omit these when using the documented
# project-local fallbacks.
install.packages(c("trend", "modifiedmk", "mutoss", "FNN"))
```

```bash
python3 -m pip install -r 04_multiscale_Moran_block_diagnostics/requirements.txt
Rscript scripts/check_environment.R
```

The check covers analytical dependencies. Package availability alone does not verify an analysis.

Record the actual full-data environment with `sessionInfo()` and `python3 -m pip freeze` after validation. This submission package does not include an unpinned environment template.
