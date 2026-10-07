# v5.6 real-data subset audit

Run from a disposable copy of the engine directory after the synthetic E2E has passed. Set `DATA_ROOT` to the real input root and `BREAK_CONFIG` to the intended analysis configuration. These diagnostics generate intermediate files within that copy.

```r
Sys.setenv(MAX_CORES = "4")
source("09_run_real_subset_audit.R")
```

Default mode is `standard`. It keeps the publication SupF MC and response-specific rho calibration settings, but reduces expensive bootstrap repetitions ; no figures are rendered.

Optional controls:

```r
Sys.setenv(REAL_SUBSET_MODE = "quick")    # quick | standard | full
Sys.setenv(REAL_SUBSET_SIDE = "256")      # crop side length in raster cells
Sys.setenv(REAL_SUBSET_REUSE = "true")    # reuse a crop for the same input fingerprint/window
Sys.setenv(MAX_CORES = "4")
source("09_run_real_subset_audit.R")
```

The script prints the final report path. The key files in the run's `validation/` directory are:

- `REAL_DATA_SUBSET_AUDIT_REPORT.md`
- `real_data_subset_audit_metrics.csv`
- `real_data_subset_audit_summary.json`
- `full_input_geometry.csv`
- `full_input_geometry_signature.csv`
- `full_input_layer_stats.csv`
- `full_input_complete_series_coverage.csv`
- `subset_selection_candidates.csv`
- `subset_selection_chosen.csv`

Interpretation rule: the crop is for implementation/data auditing. Its BH/BY family is smaller than the full study family, and clipped regional trajectories do not represent full regions. Do not quote subset significance fractions or regional turning points as final scientific results.


## Region identifier configuration

The real-data default now separates feature-level IDs from display labels and broader mountain-system grouping:

```yaml
data:
  region_id_field: "Reg_EN"
  region_label_field: "Reg_CN"
  region_group_field: "L2_code"
  region_group_label_field: "L2_name"
```

`Reg_EN` is the machine-stable ID for the five province-by-mountain-system polygons. `L2_code` is intentionally **not** used as the feature ID because the Changbaishan-Qianshan group (`I2`) spans multiple province features. Regional output tables retain `region` for backward compatibility and add `region_label`, `region_group`, and `region_group_label`.

If a configured ID field is missing, preflight now reports `region_ids_complete`/`region_ids_unique` as `NOT_TESTED` rather than emitting cascading FAILs, and writes `region_id_candidates.csv` with complete nonblank unique character-field candidates.

## Exact rho-model stress gate

The real-data configuration includes `linear_diagnostic` in `sensitivity.rho.values`. This token is resolved dynamically from the current run's `rho_model_comparison.csv`, separately for each response. It is a scientific stress test only; the primary estimator remains the prevalidated hinge-residual corrected rho. The audit report prints the exact scenario and records recovery-count retention versus primary as an informational diagnostic. Subset recovery counts remain audit-only because the FDR family is cropped.


## Run and review the rho diagnostic

After inspecting the subset outputs, set `RHO_STRESS_SOURCE_RUN` to the completed subset run directory and invoke `10_run_rho_stress_gate.R`. Set `RHO_STRESS_GATE_RUN` to a completed gate directory and invoke `11_review_rho_stress_gate.R` to review its stored comparisons without rerunning inference.

GO indicates comparatively robust recovery retention; CONDITIONAL_GO records material model sensitivity; NO_GO indicates an integrity failure requiring repair. Review all conditions and Monte Carlo floor diagnostics before interpreting a full-data run. These labels do not confer a release certificate. The archive preserves historical publication workflows, but the current analysis package does not invoke them.
