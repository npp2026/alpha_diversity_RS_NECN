# Decisions needed beyond the manuscript text

| Item | Current decision and consequence |
|---|---|
| Inner CV | Five grouped random folds inside each outer-training set; MS/SI gives the outer designs but not the inner count/design. Screening is repeated inside inner training as well. |
| Feature-screened RF/XGBoost comparison | Both algorithms share each fold-local screened set and fold assignments. The revised SI Table S1 note now agrees with MS 2.4, SI Note S3.1 and the nested screening implemented here; the former wording conflict is resolved. Screening is repeated within outer and inner training folds. This alignment confirms method wording, not numerical reproduction of Table S1. RF_raw and XGBoost_raw form the algorithm comparison. As confirmed by the authors, Figure2 scenario R2, source contributions and bootstrap use raw outer-test RF OOF; RF_QM is retained separately for calibration diagnostics. No undocumented XGBoost QM is added. |
| Production predictors and QM (author-confirmed) | Module 03 retains its fixed 15-variable list for both responses. VSURF screening is part of module 02 validation and does not automatically replace production predictors. Annual maps retain training-OOB QM; the Figure2 raw-prediction convention does not apply to annual maps. |
| VSURF edge case | All default VSURF steps are requested. When its prediction step is unavailable but interpretation features exist, retain those features, warn, and record selection_stage. A single available feature bypasses an inapplicable multivariable selection step. No usable interpretation/prediction features causes failure. |
| Scale families | Remove `_wN`, `_Nm`, `_Nkm` tokens, retaining variable/season identifiers; choose strongest absolute training Spearman correlation. An explicit feature/group dictionary is accepted. Review feature_groups_used.csv against the extraction nomenclature. |
| Missing plots | Use a single complete paired sample across the union of scenario features, after survey-year matching. Write sample inclusion and shared folds. A count other than 3,066 is flagged; no rows or responses are manufactured to match the MS. Duplicate plot IDs fail rather than pretending repeated measurements are independent. |
| Final RF sampling fraction | Retain 0.65 from supplied code because MS specifies trees/mtry/node but not sampling fraction. Both final models use 2,000/5/5 and permutation importance. |
| QRF parameters | Default manuscript profile fixes the SI-selected response-specific parameters. Optional tuning remains diagnostic. `FEM_QRF_PARAMETER_MODE=retuned` explicitly selects new tuning results, and the MS static export refuses this profile. Validation also uses 1,500 trees for consistency; this changes the supplied 800-tree prediction validation. |
| Old-age validation | “Within age ≥100” is implemented as both training and testing restricted to that subset. Absolute holdouts 100/120/148/152 are fixed, not recalculated sample percentiles. Small groups are reported as skipped by the inherited minimum-size guards. |
| Potential input | One positive static q95 surface at age 100, exactly aligned to the realized 1-km grid. Zero potential is excluded from the ratio. No annual-potential fallback. |
| 30 m → 1 km | Overlap-weighted raster average on the explicit QRF template, with ≥80% valid forest support. This handles the non-integer 1000/30 resolution ratio. All rasters must already share a projected metre CRS; reprojection is not silently chosen. Finite nonzero template cells define the common domain, matching the default QRF mask convention. |
| Multi-year state means | Require ≥80% finite years for 2016–2020 and for regional period means; the manuscript states this threshold for trend windows but does not separately specify a state-mean missing-year rule. |
| Break selection and testing | The supplied v5.6.2-derived engine defines breakpoint inference. Select the unconstrained minimum-SSE hinge, apply response-wise BH q ≤ .05 over all testable pixels, then require negative/positive slopes. A significant omnibus change does not constitute a separate directional test. |
| Break calibration/support | Complete 20-year 1-km series; response-specific bias-corrected hinge-residual AR(1) calibration; MC B=200000 for SR and 1000000 for Shannon. The manuscript settings are explicit. Each new run requires rho sensitivity and finite-MC-resolution review; prior run results do not certify new inputs. |
| OOD | Use squared Mahalanobis distances and their empirical .975 quantile; this yields the same gate as unsquared distances with the corresponding quantile. Range endpoints pass. Missing predictor/prediction values are excluded from the valid prediction domain. OOD rates use pixel-area weights and region-nested 50-km clusters; both denominator weights are output explicitly. |
| Regional uncertainty | Resample 50-km blocks independently within each region; preserve the same block draws across years/statistics. Overall pools block numerators and denominators. These are pointwise percentile intervals, not simultaneous bands or prediction intervals. |
| Contemporary sample | Stream independent uniform priority keys to obtain exactly 200,000 common-support pixels per region without replacement. Pool the five equal samples for Table S5 Overall. Insufficient sample sizes fail explicitly. |
| Native-resolution BH | Exact response/window-specific BH uses a vector proportional to the number of valid 30-m cells. Large-area runs require substantial RAM. No approximate histogram FDR is silently substituted. |

These choices are saved here for review before formal scientific release. They must be carried into the final Methods/SI or adjusted if the authors' original analytical intent differs. None is a claim that MS numerical results have already been regenerated.

## Predictor and statistical conventions

- The MOD15A2H LAI scale is 0.1 and MCD12Q2 EVI amplitude scale is 0.0001. Annual epoch dates are converted to each image year’s DOY before window averaging. Upstream extraction changes require new training predictors and annual raster predictors together.
- Plot year matching is shared by nested and final fitting; no response-dependent outlier trimming is imposed. Richness continuity correction and tied-knot handling use one QM implementation.
- State-class confidence intervals resample eligible blocks independently within each reporting unit as SI S6.3 requires. Paired metrics share draws on common support. Overall in cross-metric agreement covers I–V.
- Hamed–Rao adjustment currently computes autocorrelation over the retained sequence when some years are missing; Sen uses actual calendar spacing. The precise missing-lag autocorrelation rule is not established by MS/SI and remains for author review.
- The authors confirmed that the fixed 15-variable final production list is independent of VSURF screening. It is not inferred from outer/inner-fold selections and is not automatically updated by them.
- Legacy anonymous OOF objects without plot IDs can only be checked by coordinate order; identical coordinates cannot prove unique plot identity. New runs carry independent IDs.

## Domain and validation conventions

- The original QRF default treats template zero as outside, while earlier OOD code treated any non-NA value as inside. The manuscript profile now consistently uses finite nonzero template values; standardize other encodings to 1/NA before all stages.
- Tuning prediction points use the production template domain and jointly valid predictor cells. The inherited tuning sampler is random, while production kNNDM uses environmental stratification. The difference remains explicit; manuscript mode ignores new tuned choices.
- QRF environmental calibration/fill is still fitted before the CV split using environmental covariates, not responses. This inherited transductive preprocessing convention is not changed here; strict fold-local preprocessing would require a separate analytical decision.
- Native regional diagnostic mode changes only explicitly requested bootstrap/sample counts and labels its outputs. Manuscript defaults remain B=2000 and 200000 per region. Low-count tests do not estimate study uncertainty.
- The bundled breakpoint synthetic configuration uses reduced calibration and is labeled diagnostic by the adapter. It is not a production run; the original root integration test was absent from the supplied archive.

## Input and workflow conventions

- Missing `-9999` values are excluded consistently from selected annual predictor bands and OOD eligible domains. This repairs an input-contract inconsistency and can alter predictions/OOD rates when sentinel values occur.
- Module 05 defaults to 2005–2020 throughout calculation, preflight and postprocessing; explicitly request both trend windows when needed. The fixed inferential thresholds and available window definitions are unchanged.
- The GEE `Elev` → production `ELEV` handoff and external climate/soil/annual-raster preparation are not resolved by the supplied code. Authors must establish a consistent variable dictionary and source provenance.
