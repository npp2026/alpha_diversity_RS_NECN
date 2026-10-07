# Stage 4 — Multi-scale Moran's I and spatial-block diagnostics

This folder contains the publication-relevant Moran/block-bootstrap subset from the uploaded trend package. It **does not** compute the six-class trend states; Stage 5 contains the supplied R implementation of Theil–Sen + Hamed–Rao modified Mann–Kendall + BH-FDR.

## Inputs

Provide the Sen-slope rasters for the requested periods/metrics using the directory convention configured in `run_multiscale_moran_only_dual_periods.sh`, or edit the environment variables in `config/moran.env.example`.

## Main files

- `spatial_autocorrelation.py` — global Moran's I over valid slope pixels and optional spatial-block bootstrap.
- `spatial_block_bootstrap.py` — block-resampling utilities.
- `run_multiscale_moran_only_dual_periods.sh` — 30 m to 100 km multi-scale driver; publication bootstrap default is 2,000.

The driver produces Moran and block diagnostic results; figure assembly is outside this analysis-code release.


Input handoff: module 08 writes native 30-m Sen rasters for both trend windows. Use those for this scale diagnostic; module 05 remains the separate 1-km six-class workflow.


## Native regional slope handoff

Module 08 exports unscaled floating-point slopes. Stage the four files in the driver's expected layout, then use `SLOPE_INPUT_SCALE=1`. The legacy driver's default 100 applies only to older integer-encoded inputs.

```bash
export REGIONAL_DIR=/absolute/results/regional30m
export MORAN_ROOT=/absolute/results/moran
for period in 2001_2020 2005_2020; do
  mkdir -p "$MORAN_ROOT/period_$period/trend/diversity_rich_BASIC"
  mkdir -p "$MORAN_ROOT/period_$period/trend/diversity_shannon_BASIC"
  cp "$REGIONAL_DIR/Rich_tree_${period}_Sen_30m.tif" \
    "$MORAN_ROOT/period_$period/trend/diversity_rich_BASIC/Rich_tree_${period}_BASIC_TrendSlope.tif"
  cp "$REGIONAL_DIR/Shannon_wiener_${period}_Sen_30m.tif" \
    "$MORAN_ROOT/period_$period/trend/diversity_shannon_BASIC/Shannon_wiener_${period}_BASIC_TrendSlope.tif"
done
OUTPUT_ROOT="$MORAN_ROOT" SLOPE_INPUT_SCALE=1 USE_SIGNIFICANCE_MASK=0 \
  bash 04_multiscale_Moran_block_diagnostics/run_multiscale_moran_only_dual_periods.sh
```

Run this from the package root. A single positive scalar would cancel in Moran's I, but preserving the actual slope units keeps intermediate rasters and metadata correct. Configure block sizes and review the reported block-quality diagnostics using [config/moran.env.example](config/moran.env.example).
