# Stage 1 — Google Earth Engine predictor extraction

Run the JavaScript files in the Google Earth Engine Code Editor after setting `CONFIG.assetId` to the private plot FeatureCollection.

Expected plot attributes include `plot_id`, `year`, `lat`, and `long`. Plot records and coordinates are intentionally not included in this code-only archive.

## Files

- `Code_S1_Landsat_predictors.js` — Landsat 5 TM / 7 ETM+ / 8 OLI preprocessing and predictor extraction, including OLI harmonization and seasonal metrics.
- `Code_S2_MODIS_predictors.js` — MODIS fPAR, LAI, GPP, phenology, and related predictor extraction.

Both scripts use the target-year ±1 temporal matching design. These supplied GEE scripts extract predictors at plot locations; the uploaded code set did not include a separate full-domain annual image-export script.

## Extraction contracts

LAI uses 0.1; phenology amplitude uses 0.0001. Dates become per-image-year DOY before averaging. Native/analysis projection is fixed before neighborhood operations, missing composites are masked, and date-window ends are exclusive correctly. TEST_MODE restricts the Landsat years. Local JavaScript syntax checks do not verify Earth Engine execution. Use consistent extraction definitions for plot and annual raster predictors.
