# Final RF / QM / annual mapping

01 fits paired, survey-year-matched data using 2,000 trees, mtry=5, min.node.size=5. Metadata plot_id/plot_year are required. Nonfinite/sentinel inputs are excluded; no response-tail trimming. Richness QM applies the same deterministic continuity correction as nested validation through ../R/qm_core.R. The sampling fraction remains an explicit implementation choice.

Production uses the author-confirmed **fixed 15-variable list**, unchanged for both responses. VSURF screening in module 02 is separate: final fitting does not read its outputs, and those outputs never automatically replace `CONFIGS$production$predictors`. Training metadata records `predictor_selection=fixed_15`, `vsurf_screening_used=FALSE` and `qm_calibration_source=training_OOB`. The optional predefined sensitivity configurations remain separate from the default production configuration. The exact list is documented in [INPUTS_OUTPUTS.md](../README.md#required-inputs-and-outputs).

Annual mapping continues to use **OOB-fitted empirical QM** on RF means and QRF quantiles. Figure2 and source-contribution statistics use raw RF OOF predictions in module 02; this change does not disable annual QM.

The training script writes saved OOB/QM calibration metrics and paired predictions. Script 03 selects model bands before prediction, checks ranger quantile support, writes mean/q05/q50/q95/width rasters through staged task files, reports failures with nonzero exit, and cleans only its own temporary files. Output directories must be fresh. annual_manifest.csv is directly usable by module 08; RUN_COMPLETE.txt is present only after all tasks pass. Manifest VRT hashes cover the definition, not every external source pixel file; archive original inputs separately for scientific provenance.

See ../README.md#run-guide for commands. Quantile widths are not confidence intervals for a mean.
