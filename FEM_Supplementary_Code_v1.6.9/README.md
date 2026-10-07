# FEM supplementary code v1.6.9

Code for *Integrating state and trajectory for monitoring tree alpha diversity in Northeast China's forests*.

This analysis-only package contains eight analysis modules, shared R helpers, the supplied breakpoint-engine source, synthetic fixtures required by that engine, and source-provenance records. It produces statistical tables and raster outputs; figure rendering is outside this release. The breakpoint engine is derived from the supplied v5.6.2 source, while historical publication workflows are excluded.

The v1.6.9 method update retains the fixed 15-variable production model and separates it from validation-only VSURF screening. Figure 2 R², source contributions and bootstrap use raw RF OOF predictions. Annual mapping retains training-OOB quantile mapping. The package is organized for a journal supplementary-code submission: the former root `tests/`, `validation/`, `environment/` and `docs/` directories are not included. Submission-facing guidance and the MS/SI alignment record are consolidated here. The `07_break_year_analysis/upstream/break_year_v5_6_2/validation/` directory is retained because it is part of the upstream engine and its synthetic fixtures.

## Start here

1. Read [DEPENDENCIES.md](DEPENDENCIES.md) and install the analytical dependencies. Exact historical package versions were not supplied; record the versions from the validated study run.
2. Run the package integrity and syntax check:

   ```bash
   python3 scripts/validate_release.py
   ```

   If R is available, run the dependency check and the package self-check:

   ```bash
   Rscript scripts/check_environment.R
   Rscript scripts/run_MS_validation.R . ../FEM_validation_v169
   ```

   The self-check writes its evidence outside this source tree. These checks inspect package structure and selected computational contracts; they do not reproduce the study without the original inputs.
3. Prepare the [required inputs and outputs](#required-inputs-and-outputs), then follow the [run guide](#run-guide). Use fresh output directories outside this source tree.

## Dependencies and package checks

Use R >= 4.3 and Python >= 3.9. For R installations, the required GDAL/GEOS/PROJ system libraries must also be available.

| Component | Dependencies |
|---|---|
| R analytical workflows | `ranger`, `VSURF`, `xgboost`, `CAST`, `terra`, `sf`, `dplyr`, `tidyr`, `yaml`, `digest`, `jsonlite` |
| Optional trend and nearest-neighbour paths | `trend`, `modifiedmk`, `mutoss`, `FNN` (deterministic base-R fallbacks are used where documented) |
| Python Moran workflows | `numpy`, `scipy`, `rasterio`, `tqdm` |
| Launchers and release checks | Bash; Node.js for local JavaScript syntax checks |
| Remote predictor extraction | Google Earth Engine access and the study plot FeatureCollection |

```r
install.packages(c("ranger", "VSURF", "xgboost", "CAST", "terra", "sf",
                   "dplyr", "tidyr", "yaml", "digest", "jsonlite"))
# Optional package-backed diagnostics; omit when using the documented fallbacks.
install.packages(c("trend", "modifiedmk", "mutoss", "FNN"))
```

```bash
python3 -m pip install -r 04_multiscale_Moran_block_diagnostics/requirements.txt
Rscript scripts/check_environment.R
```

Package availability alone does not verify an analysis. Record `sessionInfo()` and `python3 -m pip freeze` from the final study run.

## Run guide

The commands below were reconstructed from the source interfaces and have not been executed on original study inputs. Replace every `/absolute/...` placeholder. Run from the package root and use separate fresh output directories.

### Dependency order

1. Run module 01 in Earth Engine and externally assemble matched plot predictors and annual raster predictors. Supply the shared kNNDM folds. These external steps are not completely automated by this archive.
2. Run module 02 for model validation and contrasts. Run module 03 final fitting/mapping from the matched data. Run module 06 to obtain the common template and potential surfaces.
3. Run module 08 OOD/common-grid preparation using outputs of 03 and the exact template used by 06.
4. From 08, run regional summaries → 04 Moran; prepared annual 1-km maps plus 06 static potential → 05 classification; prepared annual 1-km maps → 07 break years.

### Model validation and final mapping

```bash
Rscript 02_model_scenario_validation/model_fit/run_nested_MS.R \
  /absolute/matched_plots.csv /absolute/results/nested /absolute/shared_knndm.csv
# Optional fourth argument: /absolute/feature_groups.csv
bash 02_model_scenario_validation/run_postprocess_submission.sh /absolute/results/nested

DATA_FILE=/absolute/matched_plots.csv \
MODEL_OUTPUT_DIR=/absolute/results/models \
TRAIN_DATA_DIR=/absolute/results/train \
QM_DIAG_DIR=/absolute/results/qm_diag \
CALIB_DIR=/absolute/results/qm \
Rscript 03_final_RFM_QM_annual_mapping/01_fit_final_RFM_and_QM.R

VRT_DIR=/absolute/annual_predictor_vrts \
MODEL_DIR=/absolute/results/models \
TRAIN_DATA_DIR=/absolute/results/train \
CALIB_DIR=/absolute/results/qm \
PREDICTION_OUTPUT_DIR=/absolute/results/annual \
N_WORKERS=1 N_THREADS=1 \
Rscript 03_final_RFM_QM_annual_mapping/03_predict_annual_30m_maps.R
```

Figure 2, paired source contrasts and their bootstrap use raw outer-test RF OOF predictions. New `all_results.rds` objects expose them as `predictions` and `raw_predictions`, with `prediction_type=RF_raw`; QM is stored separately. VSURF screening remains local to validation and never automatically updates the fixed 15 production variables. Annual mapping still uses OOB-fitted QM.

For existing v1.6.9 outputs, copy `all_results.rds` into a separate results directory and run the postprocess command there. Raw RF vectors must exist as `raw_predictions`; the postprocess scripts do not infer raw predictions from QM. Retain the previous QM-based tables separately. See [UPDATE_NOTES.md](UPDATE_NOTES.md).

### Potential and common-grid preparation

```bash
export FEM_POTENTIAL_DATA_DIR=/absolute/potential_data
Rscript 06_biodiversity_potential_QRF/run_potential_workflow.R --check
Rscript 06_biodiversity_potential_QRF/run_potential_workflow.R

Rscript 08_OOD_and_regional_summaries/run_OOD_prepare_1km.R \
  /absolute/results/annual/annual_manifest.csv \
  /absolute/forest_mask_30m.tif /absolute/regions.gpkg \
  /absolute/potential_data/templ_1km.tif /absolute/results/ood

Rscript 08_OOD_and_regional_summaries/run_regional_30m.R \
  /absolute/results/ood /absolute/results/regional30m
```

The default QRF steps are `predict,oldage`, followed by static export in manuscript mode. Tuning is explicitly requested with `--steps=tune,predict,oldage`; it does not override manuscript parameters unless `FEM_QRF_PARAMETER_MODE=retuned`. Retuned output is ineligible for the manuscript static export. Module 04 needs staged native Sen rasters; use `SLOPE_INPUT_SCALE=1` for unscaled floating-point slopes from module 08 and follow its module README for the remaining handoff.

### State–trajectory and break years

```bash
export DATA_DIR=/absolute/results/ood/prepared_1km
export OBS_ANNUAL_DIRS="$DATA_DIR/Rich_tree;$DATA_DIR/Shannon_wiener"
export POTENTIAL_DIR=/absolute/potential_data/static_q95_age100
export OUT_ROOT_NAME=/absolute/results/state_trajectory
bash 05_state_trajectory_1km/run_all_submission.sh

export FEM_BREAK_DATA_DIR=/absolute/results/ood/prepared_1km
export FEM_BREAK_OUTPUT_DIR=/absolute/results/break_year
Rscript 07_break_year_analysis/run_break_year.R --check
Rscript 07_break_year_analysis/run_break_year.R
```

The state workflow defaults to `trend_2005_2020` consistently across calculation, preflight and downstream summaries. To calculate both windows, export `RUN_TREND_PERIODS_R='c("trend_2005_2020","trend_2001_2020")'`. Publication-data assembly and final acceptance use `ANALYSIS_GROUP`, default `trend_2005_2020`; use a separate output root for the other group because `fig6_final_data` is replaced during assembly.

The break workflow defaults to cumulative `sensitivity`; `--stage=primary` runs primary only. A small supplied-fixture run can be requested with `--synthetic` and a separate `FEM_BREAK_OUTPUT_DIR`; it is diagnostic and does not certify manuscript reproduction. Historical publication workflows are excluded from this package.

## Required inputs and outputs

The following contracts were written from the supplied source code. Manuscript documents and original study data are not included.

| Module | Required inputs | Main outputs / next consumer |
|---|---|---|
| 01 | Private Earth Engine FeatureCollection with `plot_id`, `year`, `lat`, `long` | Landsat/MODIS plot CSVs; join on `plot_id + year + aug_offset`, retain `plot_year/image_year` |
| 02 | Numeric matched plot CSV; shared kNNDM assignments; optional `feature,group` dictionary | `all_results.rds`, `nested_OOF_predictions.csv`, `nested_metrics.csv`, `Table_S1_algorithm_means.csv`, temporal validation, fold audits; postprocess reads `all_results.rds` |
| 03 fitting | Matched plot CSV and fixed production predictors | Final ranger models, training tables, OOB diagnostics, shared empirical-QM calibrators |
| 03 mapping | `Combined_2001.vrt` … `Combined_2020.vrt`, fitted models, training tables, QM objects | Annual mean/quantiles/width rasters and `annual_manifest.csv` for 08 |
| 06 | `train4pot.csv`, `ENV.tif`, `templ_1km.tif` | Age-specific potential/AOA diagnostics and two static q95 age-100 rasters for 05 |
| 08 preparation | Annual manifest, 30-m forest mask, five regions, exact QRF 1-km template | OOD-masked 30-m maps, 1-km annual maps, region vectors, block statistics and OOD rates |
| 08 regional | Completed OOD output directory | Native 30-m Sen/HR/BH maps, regional summaries, paired bootstrap outputs and contemporary samples |
| 04 | Native slope rasters from 08, staged in documented period/response folders | Multiscale Moran and block diagnostics |
| 05 | Prepared annual 1-km maps from 08, static potential from 06, region vector | State–trajectory classes, area CIs, metric contrasts, sensitivities and publication data tables |
| 07 | Complete annual 1-km maps from 08, region vector and protected inference profile | Break-year maps, distributions, primary/sensitivity diagnostics |

### Plot/model data

Module 02 requires `plot_id`, `plot_year`, `Lon_Export`, `Lat_Export`, `Rich_tree`, `Shannon_wiener` and all selected predictors. Module 03 requires the same plot identity/year and both responses, plus its fixed predictor set; coordinates are not used for final fitting. Keep identifiers as strings. Survey years are 2008–2012 and 2014–2017. When `image_year` is present it must equal `plot_year`. After matching, independent plot IDs must be nonempty and unique. Model columns must be numeric; nonfinite values and `-9999` are missing. A common complete sample is selected for both responses.

The kNNDM file has `plot_id,knndm_fold`, one row per plot, with complete fold labels 1–10. It is an external required input, not generated by module 02. Group dictionaries use `feature,group`; allowed groups are `Env,Prod,Het,Temp`, all nonempty and disjoint.

The production list in module 03 is:

`NIR_Std_GS`, `GRVI_Med_Aut`, `NDMI_CV_900m_Sum`, `NIR_TSD_Spr`, `EVI_mean_w7_Spr`, `Prod_GPP_CV`, `Prod_GPP_Mean`, `Prod_LAI_Max`, `Pheno_Amp`, `Pheno_GSL`, `BIO12_baseline_mean`, `BIO4_baseline_mean`, `ELEV`, `SPEI_12_annual_baseline_mean`, `STN`.

This fixed 15-variable list is author-confirmed and independent of module 02 VSURF screening. Selected fold features are validation outputs only; final fitting does not import them. Training metadata records this separation. Annual RF/QRF map predictions still use training-OOB QM.

### Scenario OOF results and postprocessing

New module 02 `all_results.rds` entries retain matched `observed`, `plot_id`, coordinates and shared folds. Each scenario stores `predictions=raw_predictions` (raw outer-test RF OOF), `prediction_type="RF_raw"`, and separate `calibrated_predictions`/`calibrated_prediction_type="RF_QM"`. The long `nested_OOF_predictions.csv` distinguishes `RF_raw`, `RF_QM` and `XGBoost_raw` via `algorithm`; calibrated outputs are diagnostics.

All three postprocessing entries select `raw_predictions` explicitly before validation and resampling. If that field is absent, only `predictions` explicitly labeled `prediction_type="RF_raw"` is accepted. Original mixed objects with QM `predictions` and raw RF `raw_predictions` are supported. QM-only or unlabeled objects fail; raw values are never inferred from QM. R², source removals, marginal contributions, substitutions, CIs and BH values all use the selected raw vectors. Tables add `Prediction_Type=RF_raw` while preserving their filenames and existing `Model=RF` column.

**Unresolved handoff:** module 01 emits `Elev`, while production expects `ELEV`. Climate, SPEI and soil predictors also require external preparation; the archive does not implement their complete merge/export workflow. Verify the elevation source and transformation before creating a matching column. Annual raster predictors must use the same definitions and named bands as training.

### Raster/manifest data

The module 03 annual manifest has `response,year,prediction,predictors,train_rds,model_rds`; module 08 requires exactly 40 unique rows, both `Rich_tree` and `Shannon_wiener` for 2001–2020. Use absolute paths. Each response uses one final model and one training table. Prediction and forest rasters are single layers; predictor stacks expose model feature names. Native rasters must align in a projected metre CRS with 30 × 30 m pixels. Positive finite forest-mask cells define forest support.

The common QRF template is one layer, 1000 × 1000 m in the same projected metre CRS. Finite nonzero cells define its domain. The 30→1000 m ratio is not an integer; module 08 uses overlap-weighted averaging and requires ≥80% valid forest support. Module 05 consumes these pre-aligned rasters rather than aggregating raw 30-m inputs.

Five multipart region features are allowed, one per code: `HLJDXAL`, `HLJXXAL`, `HLJCBS`, `JLSCBS`, `LNCBS` (I–V). Module 08 writes the shared region fields and `prepared_1km/NE_Mountain_Output/NE_Mountain_Regions_All.shp` for 05/07.

Module 06 CSV columns: `long,lat,Forest_age,Rich_tree,Shannon_wiener,DEMc,bio6_wc,bio10_wc,bio12_wc,bio17_wc,bio15_wc,bio4_wc`. `ENV.tif` bands: `DEM,BIO6,BIO10,BIO12,BIO17,BIO15,BIO4`. Module 06's README describes the manuscript and retuned profiles.

Module 05 needs one static positive potential raster per response: `Q95_Rich_tree_age100_1km.tif` and `Q95_Shannon_wiener_age100_1km.tif`. It rejects `POT_REF_YEARS_R`. RF is the 2016–2020 observed mean / age-100 q95, with high RF ≥0.80. Primary trends are 2005–2020. Module 07 instead requires all 20 annual observations, 2001–2020, and maps `Rich_tree→SR`, `Shannon_wiener→Shannon` in its configuration.

Keep new output directories outside this source tree and use fresh locations. `RUN_COMPLETE.txt` and workflow-specific acceptance tables are completion evidence only for the run that wrote them. Historical logs and freeze records are not evidence for a new run.

## Analysis modules

Module numbers identify components; use the dependency order in the run guide.

| Module | Purpose | Main entry point |
|---|---|---|
| [01](01_GEE_predictor_extraction/README.md) | Plot predictor extraction in Earth Engine | `Code_S1_Landsat_predictors.js`, `Code_S2_MODIS_predictors.js` |
| [02](02_model_scenario_validation/README.md) | Nested scenario/algorithm validation and paired contrasts | `model_fit/run_nested_MS.R` |
| [03](03_final_RFM_QM_annual_mapping/README.md) | Final RF, quantile mapping and annual 30 m predictions | `01_fit_final_RFM_and_QM.R`, `03_predict_annual_30m_maps.R` |
| [04](04_multiscale_Moran_block_diagnostics/README.md) | Multiscale Moran diagnostics | `run_multiscale_moran_only_dual_periods.sh` |
| [05](05_state_trajectory_1km/README.md) | Common-grid state–trajectory classification | `run_all_submission.sh` |
| [06](06_biodiversity_potential_QRF/README.md) | QRF potential, validation and static age-100 q95 export | `run_potential_workflow.R` |
| [07](07_break_year_analysis/README.md) | SupF/BH break-year analysis | `run_break_year.R` |
| [08](08_OOD_and_regional_summaries/README.md) | OOD masks, 1 km preparation and regional summaries | `run_OOD_prepare_1km.R`, `run_regional_30m.R` |

## Methods map

| Module | Analysis | Entry point | Code contract | Verification scope |
|---|---|---|---|---|
| 01 | Plot predictor extraction | `01_GEE_predictor_extraction/Code_S1_Landsat_predictors.js` | Target-year ±1; Landsat and MODIS predictors | Reviewed against MS 2.3 and SI S1–S2; remote GEE execution requires the study FeatureCollection |
| 02 | Nested validation | `02_model_scenario_validation/model_fit/run_nested_MS.R` | Shared outer folds; fold-local VSURF screening and tuning; Figure 2/source/bootstrap raw RF OOF; QM diagnostics separate | Reviewed against MS 2.4 and SI S3.1–S3.2; raw RF convention aligned |
| 03 | Final RF and annual mapping | `03_final_RFM_QM_annual_mapping/01_fit_final_RFM_and_QM.R` | 2000 trees; mtry 5; node 5; fixed 15 independent of VSURF; annual training-OOB QM | Reviewed against MS 2.5 and author-confirmed fixed-15 convention |
| 04 | Spatial autocorrelation | `04_multiscale_Moran_block_diagnostics/run_multiscale_moran_only_dual_periods.sh` | Native Sen maps; multiscale Moran and block diagnostics | Reviewed against MS 2.6 and SI S5; Python execution requires study rasters |
| 05 | State and trajectory | `05_state_trajectory_1km/run_all_submission.sh` | 2016–2020 mean/static age100 q95; RF≥0.80; 2005–2020 Sen/HR/BH | Reviewed against MS 2.6 and SI S6; static potential and annual inputs are required |
| 06 | Potential QRF | `06_biodiversity_potential_QRF/run_potential_workflow.R` | Manuscript response-specific parameters; age100 q95; AOA and old-age validation | Reviewed against SI S5.4; manuscript profile is default and retuned mode is diagnostic |
| 07 | Break-year analysis | `07_break_year_analysis/run_break_year.R` | Unconstrained minimum SSE; AR1 supF; per-response BH; direction after selection | Reviewed against MS 2.6 and SI S5.5; original annual series are required |
| 08 | OOD and regional summaries | `08_OOD_and_regional_summaries/run_OOD_prepare_1km.R` | Range then ridge Mahalanobis; area weights; overlap-weighted 1 km preparation | Reviewed against MS 2.5–2.6 and SI S3.3/S5; study rasters and polygons are required |

## MS/SI alignment

This record documents the method check performed against the supplied manuscript and supplementary-information files on 2026-10-07. It lets a reviewer distinguish author-confirmed conventions from checks that require study inputs.

### Methods implemented by this release

| Manuscript location | Requirement | Release implementation | Status |
|---|---|---|---|
| MS 2.4; SI S3.1 | Sequential scale selection, VIF filtering at 10, VSURF, tuning and calibration inside each outer training fold | `02_model_scenario_validation/model_fit/nested_core_MS.R` performs each operation on the relevant training data; the outer runner writes fold audits | Aligned |
| MS 2.4; SI S3.1 | Nine scenarios with shared samples and shared random and kNNDM folds | `run_nested_MS.R` uses the nine `ms_scenarios`, ten random folds and the supplied ten-fold kNNDM assignment | Aligned |
| SI Table S1; SI S3.1 | RF and XGBoost use the same fold-local screened feature sets and fold assignments for paired algorithm comparisons | `ms_outer()` creates outer and inner training-fold screens before the algorithm loop; both algorithms reuse those screens and partitions | Aligned with the revised SI Table S1 note |
| MS 2.4 and Figure 2 | Scenario R², source removal, marginal contributions, substitutions and bootstrap use raw RF outer-test OOF predictions | `raw_predictions` is selected explicitly and labeled `RF_raw`; QM vectors are retained as separate diagnostics | Aligned |
| Author-confirmed production convention | Final annual model keeps the fixed 15 predictors; VSURF does not replace them | Module 03 checks the fixed list and records `predictor_selection=fixed_15` and `vsurf_screening_used=FALSE` | Aligned |
| MS 2.5; SI S3.2 | Annual mapping uses training-OOB empirical quantile mapping | Module 03 fits the shared QM object from training OOB predictions and module 03.2 applies it to annual outputs | Aligned |
| MS 2.6; SI S5.4 and S6 | Age-100 q95 potential, 1-km common support, RF threshold 0.80, 2005–2020 trend and paired spatial-block bootstrap | Modules 05, 06 and 08 enforce the static age-100 q95 inputs, common grid, threshold and bootstrap contracts | Aligned |
| SI S5.4 | QRF candidates and manuscript profile: 18 candidates, top-five kNNDM diagnostic, richness 4/10/forced age, Shannon 2/10/no forced age, 1,500 trees | Module 06 retains tuning as a diagnostic and uses the manuscript profile by default; retuned mode is explicit and cannot produce the manuscript static export | Aligned |
| SI S5.4 | Age holdouts at 100, 120, 148 and 152 years; repeated ten-fold validation within age ≥100, three repetitions | `oldage_validation_and_aoa_coverage.R` uses those four absolute holdouts and restricts both training and testing of the repeated validation to age ≥100 | Aligned |

### SI Table S1 screening alignment

The Table S1 note in the supplied `SI_f_STYLEPASS_FINAL.docx` now states:

> Feature screening was repeated within each outer training fold. For the algorithm-choice comparison, RF and XGBoost used the same fold-local screened feature set and fold assignments, so paired performance differences reflect algorithm choice rather than feature-set differences.

This wording agrees with MS 2.4, SI Note S3.1 and `ms_outer()` in `02_model_scenario_validation/model_fit/nested_core_MS.R`. Screening uses training rows only, and RF and XGBoost reuse the same selected features and fold assignments. Matching is between algorithms within each fold; selected feature sets may vary across folds. Screening is also repeated within inner training folds for tuning, and both algorithms use the nested Latin hypercube framework described in Note S3.1.

The former SI Table S1 screening conflict is resolved in the supplied SI. This documentation update records that alignment; it does not change the analytical code or certify numerical reproduction of Table S1.

### Claims deliberately not certified by this package

The release does not include the study rasters, shared fold file, original plot data, or the final validated R environment. It therefore does not claim numerical reproduction of manuscript tables or figures. R syntax, R self-checks, model fitting, annual mapping and original-data reproduction must be run in an environment satisfying [DEPENDENCIES.md](DEPENDENCIES.md). Figure rendering remains outside this analysis-only release.

## Reproduction scope and package boundaries

Original plot records, annual predictor rasters, shared spatial folds and study outputs are not included. The GEE scripts extract plot predictors; full-domain annual raster preparation remains external. Current package checks do not establish agreement with manuscript numbers.

[IMPLEMENTATION_CHOICES.md](IMPLEMENTATION_CHOICES.md) records scientific assumptions and unresolved handoffs. [CITATION.cff](CITATION.cff) contains citation metadata. [SOURCE_SPECIFICATION.json](SOURCE_SPECIFICATION.json) identifies the supplied base and inherited breakpoint source; module provenance CSVs record original and current source hashes. No publication DOI or license terms were supplied.

Historical freeze records, old review reports, old checksums and their dedicated workflows are outside this package. The separately delivered `FEM_Historical_Audit_Archive_v1.6.8.zip` preserves the exact supplied base and historical verification records. Analysis and current checks do not require that archive. `FILE_MANIFEST.csv` and `SHA256SUMS.txt` describe only this release.
