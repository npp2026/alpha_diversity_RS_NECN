from pathlib import Path
import re,sys,yaml
root=Path(__file__).resolve().parents[1];rfiles=list(root.rglob('*.R'));errors=[];warnings=[]
def strip_strings_comments(s):
 out=[];i=0;quote=None;esc=False
 while i<len(s):
  c=s[i]
  if quote:
   if esc:esc=False
   elif c=='\\':esc=True
   elif c==quote:quote=None
   out.append(' ');i+=1;continue
  if c in ('"',"'"):quote=c;out.append(' ');i+=1;continue
  if c=='#':
   while i<len(s) and s[i]!='\n':out.append(' ');i+=1
   continue
  out.append(c);i+=1
 return ''.join(out)
for f in rfiles:
 s=f.read_text(encoding='utf-8');z=strip_strings_comments(s);stack=[];pairs={')':'(',']':'[','}':'{'}
 for j,c in enumerate(z):
  if c in '([{':stack.append((c,j))
  elif c in ')]}':
   if not stack or stack[-1][0]!=pairs[c]:errors.append(f'{f}: mismatched {c} at {j}');break
   stack.pop()
 if stack:errors.append(f'{f}: unclosed delimiters {stack[-3:]}')

# Catch a common generated-R failure class: backslash escapes that R does not recognize.
valid_escape=set("abfnrtv\\'\"01234567xuU")
for f in rfiles:
    raw=f.read_text(encoding='utf-8');i=0;line=1;quote=None
    while i<len(raw):
        c=raw[i]
        if c=='\n': line+=1
        if quote is None and c in ("'", '"'):
            quote=c;i+=1;continue
        if quote is not None:
            if c=='\\':
                if i+1<len(raw) and raw[i+1] not in valid_escape: errors.append(f'{f}: possible invalid R escape at line {line}: \\{raw[i+1]}')
                i+=2;continue
            if c==quote: quote=None
        i+=1

# Catch accidental same-line statement concatenation introduced by text-generation edits.
for f in rfiles:
    for line_no,line in enumerate(f.read_text(encoding='utf-8').splitlines(),1):
        if re.search(r'\)\s{2,}[A-Za-z_.][A-Za-z0-9_.]*\s*<-',line) and ';' not in line:
            errors.append(f'{f}: suspicious concatenated assignment at line {line_no}')

cfg=yaml.safe_load((root/'config/v5_6.yml').read_text(encoding='utf-8')); scfg=yaml.safe_load((root/'config/v5_6_synthetic.yml').read_text(encoding='utf-8'))
text=lambda p:(root/p).read_text(encoding='utf-8')
checks={
 'robust zero-RSS/Inf SupF':'exact_hinge' in text('R/core_break.R') and 'perfect_linear' in text('R/core_break.R'),
 'scale-aware SSE tolerance':'roundoff <-' in text('R/core_break.R') and 'max(1, scale)' not in text('R/core_break.R'),
 'MC ties use >=':'left.open = TRUE' in text('R/ar1_calibration.R') and cfg['supf']['p_rule']=='greater_equal',
 'discrete block FDR':'discrete_fdr_lookup' in text('R/fdr_discrete.R') and 'FLT8S' in text('R/fdr_discrete.R'),
 'live YAML':'yaml::read_yaml' in text('R/utils.R') and 'read_config' in text('R/pipeline.R'),
 'complete-series inference guard':'require_complete_series' in text('R/raster_pipeline.R') and cfg['data']['require_complete_series'],
 'terra app vector/matrix contract':'is.matrix(v)' in text('R/raster_pipeline.R') and 'vapply(seq_len(nrow(v))' in text('R/raster_pipeline.R'),
 'boundary 5/4/3/2':cfg['sensitivity']['boundary']['min_segments']==[5,4,3,2],
 'response-specific corrected rho':'rho_for_calibration' in text('R/rho_estimation.R') and 'rho_by_response' in text('R/raster_pipeline.R'),
 'break-robust primary rho':cfg['rho_estimation']['residual_model']=='hinge' and scfg['rho_estimation']['residual_model']=='hinge',
 'calibration cache includes RNG state':'master_seed' in text('R/ar1_calibration.R') and 'rng_kind' in text('R/ar1_calibration.R') and 'chunk' in text('R/ar1_calibration.R'),
 'rho bias mapping cache':'rho_bias_mapping_cached' in text('R/rho_estimation.R') and '_cache' in text('R/rho_estimation.R'),
 'frozen-FDR bootstrap':'frozen_fdr_p_cutoff' in text('R/bootstrap_pixel.R') and 'bootstrap_rejection_frequency' in text('R/bootstrap_pixel.R'),
 'regional trajectory bootstrap':'bootstrap_regional_trajectory' in text('R/bootstrap_trajectory.R'),
 'cluster IoU full denominator':'primary_n+z$other_n-z$intersection' in text('R/spatial_validation.R'),
 'spatial temporal agreement':'break_year_agreement_rasters' in text('R/spatial_validation.R') and 'within1_year_agreement' in text('R/spatial_validation.R'),
 'synthetic fixture integrated':'generate_synthetic_dataset' in text('R/synthetic_data.R') and (root/'08_run_synthetic_e2e.R').exists() and (root/'validation/synthetic_spec.yml').exists(),
 'area-weighted regional aggregation':'cellSize' in text('R/aggregation.R') and 'valid_area_km2' in text('R/aggregation.R'),
 'cell-center region zones':'touches = FALSE' in text('R/aggregation.R') and ('touches=FALSE' in text('R/tables.R') or 'touches = FALSE' in text('R/tables.R')) and 'touches = TRUE' not in text('R/aggregation.R') and 'touches=TRUE' not in text('R/tables.R') and 'touches = TRUE' not in text('R/tables.R'),
 'pipeline failure manifest':'completed_partial' in text('R/pipeline.R') and 'error_message' in text('R/provenance.R'),
 'PSOCK parallel runtime':'makePSOCKcluster' in text('R/parallel_runtime.R') and 'nested_parallelism' in text('R/parallel_runtime.R'),
 'response parallel':'run_response_scenario_task' in text('R/raster_pipeline.R') and 'stage = "response"' in text('R/raster_pipeline.R'),
 'scenario parallel':'run_sensitivity_scenario_task' in text('R/sensitivity.R') and 'stage = "scenario"' in text('R/sensitivity.R'),
 'response-specific rho diagnostic stress':'linear_diagnostic' in cfg['sensitivity']['rho']['values'] and 'rho_sensitivity_sources' in text('R/sensitivity.R') and 'rho_sensitivity_values.csv' in text('R/sensitivity.R'),
 'targeted rho stress gate':(root/'10_run_rho_stress_gate.R').exists() and (root/'validation/rho_stress_gate.R').exists() and 'RHO_STRESS_GATE_REPORT.md' in text('validation/rho_stress_gate.R'),
 'rho stress diagnostic decision':(root/'11_review_rho_stress_gate.R').exists() and 'rho_stress_gate_decision.csv' in text('validation/rho_stress_gate.R') and 'CONDITIONAL_GO' in text('validation/rho_stress_gate.R') and cfg['validation']['rho_stress_gate']['robust_retention_min']==0.5,
 'analysis-only pipeline': 'publication_freeze_preflight_v56' not in text('R/pipeline.R') and 'publication_freeze_finalize_v56' not in text('R/pipeline.R') and 'publication_freeze' not in cfg and not (root/'12_run_publication_freeze.R').exists(),
 'deterministic MC chunk parallel':'supf_mc_chunk_task' in text('R/ar1_calibration.R') and '"mc_chunk"' in text('R/ar1_calibration.R') and 'stage = "supf_mc"' in text('R/ar1_calibration.R'),
 'pixel batch parallel':'pixel_bootstrap_batch_task' in text('R/bootstrap_pixel.R') and 'stage = "pixel_bootstrap"' in text('R/bootstrap_pixel.R'),
 'parallel reproducibility tests':(root/'validation/parallel_reproducibility.R').exists() and 'mc_serial_vs_psock' in text('validation/parallel_reproducibility.R'),
 'warning regression tests':(root/'validation/warning_regression.R').exists() and 'readStart' in text('validation/warning_regression.R'),
 'region feature/group metadata':cfg['data']['region_id_field']=='Reg_EN' and cfg['data'].get('region_label_field')=='Reg_CN' and cfg['data'].get('region_group_field')=='L2_code' and 'region_metadata' in text('R/utils.R') and 'region_group_label' in text('R/aggregation.R'),
 'region preflight candidate diagnostics':'NOT_TESTED' in text('validation/real_data_subset_audit.R') and 'region_id_candidates.csv' in text('validation/real_data_subset_audit.R') and (root/'validation/region_metadata_tests.R').exists(),
 'real-data subset audit integrated':(root/'09_run_real_subset_audit.R').exists() and (root/'validation/real_data_subset_audit.R').exists() and (root/'validation/real_data_subset_spec.yml').exists() and 'REAL_DATA_SUBSET_AUDIT_REPORT.md' in text('validation/real_data_subset_audit.R'),
}
old=['read_break_config','initialize_run','require_pkgs','set_reproducible_rng','mc_p_value','simulate_ar1_stationary','nrow(bs)']
allr='\n'.join(f.read_text(encoding='utf-8') for f in rfiles)
if any(x in allr for x in ('terra::readStart', 'terra::readStop', 'terra::readValues')):
 errors.append('manual terra read connection API remains; use blockwise terra::values(row=, nrows=)')
for x in old:
 if x in allr:errors.append('obsolete/broken API remains: '+x)
for k,v in checks.items():
 if not v:errors.append('failed feature check: '+k)
print(f'R files checked: {len(rfiles)}')
for k,v in checks.items():print(('PASS' if v else 'FAIL'),k)
if errors:
 print('\nERRORS:');print('\n'.join(errors));sys.exit(1)
print('Static audit passed.')
