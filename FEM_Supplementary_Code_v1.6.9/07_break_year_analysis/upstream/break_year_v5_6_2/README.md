# Breakpoint analysis engine

This is the FEM 1.6.9 analysis derivative of the supplied v5.6.2 code. The directory name identifies its source lineage. Numerical algorithms and retained statistical settings are preserved; publication-freeze workflows, historical certifications and old result logs are excluded.

Use the [FEM adapter](../../README.md) for manuscript runs. Its response-specific Monte Carlo settings and [method guard](../../method_profile.R) enforce the manuscript profile. Direct engine launchers use `config/v5_6.yml` unless `BREAK_CONFIG` selects another configuration; that inherited baseline is not the FEM manuscript profile. Set `DATA_ROOT` to real inputs and `OUTPUT_ROOT` outside the source package. Paths in the baseline configuration are relative to `DATA_ROOT`.

| Retained entry point | Purpose |
|---|---|
| `00_run_all.R` | All cumulative analysis and validation stages |
| `01_run_primary.R`–`05_build_tables.R` | Cumulative primary, sensitivity, uncertainty, spatial and table stages |
| `07_run_validation.R` | Algorithm unit and smoke checks |
| `08_run_synthetic_e2e.R` | Synthetic end-to-end diagnostic |
| `09_run_real_subset_audit.R` | Geometry, data and implementation checks on a real-data crop |
| `10_run_rho_stress_gate.R` | Exact linear-residual rho stress diagnostic from an existing subset run |
| `11_review_rho_stress_gate.R` | Re-evaluate an existing diagnostic using its recorded comparison |

The upstream synthetic/subset tools can generate intermediate files beneath their project directory. Run those tools in a disposable copy; use the FEM adapter with an external output directory for ordinary analysis.

From this directory, available checks are:

```bash
python3 validation/python_reference_checks.py
python3 validation/static_audit.py
Rscript validation/run_all_tests.R
Rscript validation/run_all_tests.R --synthetic --reuse-data
```

The Python reference check is independent of R execution. The static audit checks source patterns and delimiters; it is not an R parser. `validation/synthetic_reference_test.py` additionally requires NumPy, pandas, SciPy, PyYAML and rasterio, and writes its own result CSVs. Existing synthetic inputs and expected truth are retained; old test-result CSVs are not current evidence.

See [REAL_DATA_SUBSET_AUDIT_GUIDE.md](REAL_DATA_SUBSET_AUDIT_GUIDE.md) for subset and rho diagnostics. `validation/risk_register.csv` lists current method limitations; inherited resolved-run claims have been removed. Diagnostic GO/CONDITIONAL_GO/NO_GO labels do not certify publication readiness.

The pixel bootstrap still uses `frozen_primary_fdr_p_cutoff`: this fixes the primary rejection threshold across bootstrap draws and is part of the numerical method. It must not be removed as release metadata.
