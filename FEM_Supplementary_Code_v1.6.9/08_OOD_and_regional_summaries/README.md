# OOD and native-resolution regional summaries
run_OOD_prepare_1km.R takes a 40-row annual manifest, a positive/NA forest mask, five region polygons, the exact QRF template and an empty output directory. It fits a training reference from the saved final ranger predictor names and training table; applies range first, then ridge Mahalanobis; writes area-weighted gate-specific rates with their correct denominators and 2000 spatial-cluster percentile intervals. It masks OOD before native summaries and overlap-weighted common-grid preparation.

run_regional_30m.R takes that output directory plus a new empty output directory. It writes native 30-m pixel Sen/HR/BH, regional mean pixel slopes, early-to-late relative changes, paired regional contrasts, annual pixel-weighted means with pointwise region-nested50-km CIs, and balanced200000-per-region common-support contemporary samples. Exact native BH can require substantial memory.

See ../README.md#run-guide for the full schema and commands and ../IMPLEMENTATION_CHOICES.md for decisions not explicit in MS. The block table, saved bootstrap replicates, manifests and sessionInfo make each numerical summary auditable.

Diagnostics separately report range failure / valid domain, MD failure / range-passing domain, and total OOD / valid domain. Support diagnostics live below prepared_1km/<response>/support/. Forest values outside the positive finite mask never enter aggregation. Regional Table S6 uses paired draws for the two trend windows.

Period means explicitly use terra::app; base mean dispatch is not used for SpatRaster objects. The common template excludes zero and nonfinite cells. The native regional entrypoint supports explicit --profile=diagnostic --bootstrap=B --sample-per-region=N for low-cost testing. Counts cannot be overridden in manuscript mode; the diagnostic profile is recorded in outputs. Default manuscript calculations and counts remain unchanged.
