# Break-year analysis — manuscript method

`run_break_year.R` uses the supplied v5.6.2-derived analysis engine. [SOURCE_PROVENANCE.csv](SOURCE_PROVENANCE.csv) records the inherited source member/hash and the current packaged hash for every retained engine file. Source identity is documented in [SOURCE_SPECIFICATION.json](../SOURCE_SPECIFICATION.json).

## Method and execution

Use complete 2001–2020 annual 1-km inputs prepared by module 08. Set `FEM_BREAK_DATA_DIR` and `FEM_BREAK_OUTPUT_DIR`, run `--check`, then run normally. The default stage is `sensitivity`: primary 5/5 inference, boundary settings 4/4, 3/3, 2/2, temporal-dependence sensitivity and BH/BY comparisons. `--stage=primary` runs only primary inference. Stages are cumulative and start a new run, rather than resume an existing one.

The minimum-SSE continuous hinge is selected over all admissible years before testing slope direction. Primary inference uses stationary Gaussian AR(1) Monte Carlo supF p values, BH q ≤ .05 across all testable pixels within each response, and then negative-pre/positive-post classification. A failed direction filter does not trigger a new candidate search. The 5/5 candidates are 2005–2015; the breakpoint observation belongs to the left segment. All 20 observations are required.

[config/fem_break_year.yml](config/fem_break_year.yml) explicitly sets B=200000 for SR and B=1000000 for Shannon, preserving the supplied manuscript calibration. Hinge-residual rho is estimated separately by response and corrected by simulation inversion. [`method_profile.R`](method_profile.R) rejects changes to the manuscript's primary calibration, boundary/rho sensitivity, RNG and testing settings. `--config=/absolute/config.yml` can customize input locations and operational settings. Reduced calibration requires an explicit `--diagnostic` flag and is labeled in `FEM_RUN_METADATA.json`. `--synthetic` uses the bundled upstream diagnostic fixture and its lower calibration settings.

## Outputs and manuscript figures

Each run has upstream `primary/`, `rho/` and `sensitivity/` directories, configuration/input manifests and logs. Each response/scenario writes `core_multiband.tif` (10 layers), separate named layers, `q_BH.tif`, `q_BY.tif`, `significant_break.tif`, `recovery_break.tif`, `summary.csv`, and `break_year_distribution.csv`. The distribution denominator is the number of decline–recovery pixels, not all significant breaks. Boundary summaries are in `sensitivity/boundary/boundary_sensitivity_summary.csv`; inspect rho and multiplicity summaries alongside them.

Figure S5A uses primary decline–recovery distributions; S6 uses matching distributions for each segment setting. Module 08 supplies S5B/C annual means and region-nested spatial-block intervals. Optional stages `uncertainty`, `spatial`, `tables`, and `validation` retain the analytical outputs, including a different AR(1) regional trajectory bootstrap. Their trajectory panels must not silently replace the MS spatial-bootstrap panels. Default execution does not run these optional stages.

## Validation and scope

Algorithm tests, synthetic fixtures, real-data subset checks and rho stress diagnostics remain available; see the [engine guide](upstream/break_year_v5_6_2/README.md). Their decisions apply to the inputs of that run. Inspect finite-MC resolution and rho-model sensitivity for every new analysis.

Historical freeze/refreeze/finalizer entry points and evidence snapshots are archived separately. Their modes are rejected explicitly even when invoking the retained engine launchers directly. The analytical bootstrap's fixed primary FDR cutoff remains part of the method and is unrelated to historical release freezing.

No study rasters were supplied for this cleanup. Manuscript numbers were not regenerated, and R parsing and execution remain unverified.
