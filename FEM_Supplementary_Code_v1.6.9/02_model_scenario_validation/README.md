# Nested model and scenario validation

Use `model_fit/run_nested_MS.R` with matched numeric plot data and a separate shared kNNDM fold CSV. An optional fourth argument supplies a feature/group dictionary. [The run guide](../README.md#run-guide) gives exact commands and [the input contract](../README.md#required-inputs-and-outputs) defines the schemas.

Nine scenarios share samples and paired outer folds across RF/XGBoost and responses. Screening and tuning are repeated within training data. Figure2 scenario R2, source contributions and their bootstrap use **raw outer-test RF OOF predictions**. Fold-local RF_QM predictions are exported separately for calibration diagnostics; they are not inputs to these statistics. Review the recorded choices in [IMPLEMENTATION_CHOICES.md](../IMPLEMENTATION_CHOICES.md).

VSURF screening belongs to this validation workflow. Its selected features and fold audits do not automatically replace the independent fixed 15-variable production list in module 03.

`run_postprocess_submission.sh RESULTS_DIR` validates paired OOF identities, computes 50 km block-bootstrap/BH statistics and source contrasts. `postprocess/Make_consolidated_Table_S1.R` retains a legacy filename for Tables S3/S4 contrasts. The current algorithm Table S1 comes from `model_fit/run_nested_MS.R`, as `Table_S1_algorithm_means.csv`.

All three postprocessing scripts require an explicit results directory. Interactive calls to the two statistical scripts can still use `FIG23_OUTPUT_DIR` and `FIG4_OUTPUT_DIR`.

Each entry explicitly selects `raw_predictions` before pairing checks and bootstrap. This also supports original v1.6.9 results where `predictions` contained QM and `raw_predictions` contained raw RF. Alternatively, raw-only objects must label `predictions` with `prediction_type="RF_raw"`; missing raw inputs fail rather than falling back to QM or an unlabeled vector. The input RDS is not rewritten. Output CSVs add `Prediction_Type=RF_raw` while retaining their filenames and `Model=RF`.

The model launcher requires explicit `DATA_CSV`, `OUTPUT_DIR` and `SHARED_KNNDM_FOLD_CSV`; paths are resolved from the caller's working directory. Older model implementations and default development result paths are excluded from this release.
