#!/usr/bin/env bash
# =============================================================================
# run_multiscale_moran_only_dual_periods.sh
# =============================================================================
# Purpose:
#   Compute ONLY Global Moran's I at multiple aggregation grains for two periods
#   and two response variables, using the existing spatial_autocorrelation.py.
#
# Default behavior:
#   Moran's I is computed from ALL valid slope pixels.
#   It does NOT filter by MK p-value (Mkp < 0.05) or FDR q-value (q < 0.05).
#
# Optional behavior:
#   Set USE_SIGNIFICANCE_MASK=1 to compute Moran's I only on significant pixels.
#   If USE_FDR=1, use *_MannKendall_FDR.tif band 3 as q-value.
#   If USE_FDR=0, use *_MannKendall.tif band 2 as raw MK p-value.
#
# Examples:
#   chmod +x run_multiscale_moran_only_dual_periods.sh
#   FORCE=1 ./run_multiscale_moran_only_dual_periods.sh
#   ONLY="30m,250m,500m,1km,5km,25km" FORCE=1 ./run_multiscale_moran_only_dual_periods.sh
#   USE_SIGNIFICANCE_MASK=1 USE_FDR=1 FORCE=1 ./run_multiscale_moran_only_dual_periods.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-}"
if [[ -n "$CONFIG_FILE" ]]; then
  if [[ ! -r "$CONFIG_FILE" ]]; then
    printf '[FAIL] Configuration file is not readable: %s\n' "$CONFIG_FILE" >&2
    exit 1
  fi
  # shellcheck disable=SC1090
  set -a
  source "$CONFIG_FILE"
  set +a
fi

# Prevent nested BLAS/OpenMP pools from multiplying the explicit Python worker
# count. The block code supplies its own shared-memory parallelism.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

PROJECT_DIR="${PROJECT_DIR:-$SCRIPT_DIR}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
AUTOCORR_PY="${AUTOCORR_PY:-${PROJECT_DIR}/spatial_autocorrelation.py}"
OUTPUT_ROOT="${OUTPUT_ROOT:-${PROJECT_DIR}/output}"

# Periods to run. Use comma-separated values, e.g. PERIODS="2001-2020,2005-2020".
PERIODS="${PERIODS:-2001-2020,2005-2020}"

# Response variable names used by the trend pipeline.
SR_RESPONSE="${SR_RESPONSE:-Rich_tree}"
SH_RESPONSE="${SH_RESPONSE:-Shannon_wiener}"

# Moran mask mode.
# Default = 0: all valid slope pixels, ignoring Mkp/q significance masks.
USE_SIGNIFICANCE_MASK="${USE_SIGNIFICANCE_MASK:-0}"
USE_FDR="${USE_FDR:-0}"
P_THRESHOLD="${P_THRESHOLD:-0.05}"

# Raster scale factors and Moran settings.
SLOPE_INPUT_SCALE="${SLOPE_INPUT_SCALE:-100}"
MK_PVALUE_SCALE="${MK_PVALUE_SCALE:-10000}"
AGGREGATE_EDGE_MODE="${AGGREGATE_EDGE_MODE:-partial}"
AGGREGATE_MIN_VALID_FRACTION="${AGGREGATE_MIN_VALID_FRACTION:-0}"
NEIGHBOR="${NEIGHBOR:-queen}"
KERNEL_RADIUS="${KERNEL_RADIUS:-1}"
N_PERMUTATIONS="${N_PERMUTATIONS:-0}"
SEED="${SEED:-42}"
MIN_MORAN_PIXELS="${MIN_MORAN_PIXELS:-30}"
SMALL_MORAN_SAMPLE_ACTION="${SMALL_MORAN_SAMPLE_ACTION:-point-only}"
MIN_LOCAL_PIXELS="${MIN_LOCAL_PIXELS:-30}"

# Spatial block-bootstrap CI. Publication default: 2,000 replicates.
# Worker counts are intentionally conservative and can be overridden through
# environment variables without changing the statistical settings.
N_SPATIAL_BOOTSTRAP="${N_SPATIAL_BOOTSTRAP:-2000}"
BOOTSTRAP_WORKERS="${BOOTSTRAP_WORKERS:-1}"
BOOTSTRAP_BLOCK_WORKERS="${BOOTSTRAP_BLOCK_WORKERS:-1}"
BOOTSTRAP_CI_LEVEL="${BOOTSTRAP_CI_LEVEL:-0.95}"
BOOTSTRAP_MIN_VALID_PIXELS="${BOOTSTRAP_MIN_VALID_PIXELS:-30}"
BOOTSTRAP_MIN_VALID_FRACTION="${BOOTSTRAP_MIN_VALID_FRACTION:-0}"
BOOTSTRAP_RANDOM_GRID_OFFSET="${BOOTSTRAP_RANDOM_GRID_OFFSET:-0}"
BOOTSTRAP_EDGE_MODE="${BOOTSTRAP_EDGE_MODE:-balanced}"
BOOTSTRAP_INCLUDE_PARTIAL_EDGE_BLOCKS="${BOOTSTRAP_INCLUDE_PARTIAL_EDGE_BLOCKS:-0}"  # deprecated alias
BOOTSTRAP_MAX_PARTITION_DIFFERENCE="${BOOTSTRAP_MAX_PARTITION_DIFFERENCE:-0.05}"
BOOTSTRAP_MAX_VALID_PIXEL_CV="${BOOTSTRAP_MAX_VALID_PIXEL_CV:-0.50}"
BOOTSTRAP_MIN_USABLE_BLOCKS="${BOOTSTRAP_MIN_USABLE_BLOCKS:-20}"
BOOTSTRAP_SAVE_DISTRIBUTION="${BOOTSTRAP_SAVE_DISTRIBUTION:-1}"

# Either set one global block size, e.g. BOOTSTRAP_BLOCK_CELLS=20, or a
# scale-specific map, e.g.
# BOOTSTRAP_BLOCK_CELLS_BY_SCALE="30m:200,100m:100,250m:60,500m:40,1km:25,5km:8,10km:5,25km:3,50km:2,75km:2,100km:2"
# Values are block side lengths in CELLS AFTER aggregation. The same scale map
# is automatically used for both periods, ensuring comparable CIs.
BOOTSTRAP_BLOCK_CELLS="${BOOTSTRAP_BLOCK_CELLS:-}"
BOOTSTRAP_BLOCK_CELLS_BY_SCALE="${BOOTSTRAP_BLOCK_CELLS_BY_SCALE:-}"

FOREST_MASK="${FOREST_MASK:-}"
FORCE="${FORCE:-0}"
FAIL_FAST="${FAIL_FAST:-1}"
DISABLE_RUN_LOCK="${DISABLE_RUN_LOCK:-0}"
LOG_LEVEL="${LOG_LEVEL:-INFO}"

if [[ "$USE_FDR" == "1" ]]; then
  MK_PVALUE_BAND="${MK_PVALUE_BAND:-3}"   # band 3 in *_MannKendall_FDR.tif = q-value
else
  MK_PVALUE_BAND="${MK_PVALUE_BAND:-2}"   # band 2 in *_MannKendall.tif = raw MK p-value
fi

# Override with ONLY="500m,1km" if needed.
DEFAULT_SCALES="30m,100m,250m,500m,1km,5km,10km,25km,50km,75km,100km"
SCALES_TO_RUN="${ONLY:-$DEFAULT_SCALES}"
IFS=',' read -ra SCALES_ARR <<< "$SCALES_TO_RUN"

# Parse optional scale:block_cells map.
declare -A BOOT_BLOCK_MAP=()
if [[ -n "$BOOTSTRAP_BLOCK_CELLS_BY_SCALE" ]]; then
  IFS=',' read -ra _BOOT_PAIRS <<< "$BOOTSTRAP_BLOCK_CELLS_BY_SCALE"
  for _pair in "${_BOOT_PAIRS[@]}"; do
    _pair="$(echo "$_pair" | tr -d '[:space:]')"
    [[ -z "$_pair" ]] && continue
    if [[ ! "$_pair" =~ ^([^:]+):([0-9]+)$ ]]; then
      echo "[FAIL] Invalid BOOTSTRAP_BLOCK_CELLS_BY_SCALE entry: $_pair" >&2
      exit 1
    fi
    if (( 10#${BASH_REMATCH[2]} < 2 )); then
      echo "[FAIL] Bootstrap block size must be >=2 cells: $_pair" >&2
      exit 1
    fi
    BOOT_BLOCK_MAP["${BASH_REMATCH[1]}"]="${BASH_REMATCH[2]}"
  done
fi

get_bootstrap_block_cells() {
  local scale_label="$1"
  if [[ -n "${BOOT_BLOCK_MAP[$scale_label]:-}" ]]; then
    printf '%s' "${BOOT_BLOCK_MAP[$scale_label]}"
  elif [[ -n "$BOOTSTRAP_BLOCK_CELLS" ]]; then
    printf '%s' "$BOOTSTRAP_BLOCK_CELLS"
  else
    printf ''
  fi
}

# Aggregation factor relative to the original 30 m raster.
declare -A LABEL_FACTOR=(
  [30m]=1
  [100m]=3
  [250m]=8
  [500m]=17
  [1km]=33
  [5km]=167
  [10km]=333
  [25km]=833
  [50km]=1667
  [75km]=2500
  [100km]=3333
)

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

ACTIVE_STAGE_DIR=""
ACTIVE_BACKUP_DIR=""
ACTIVE_FINAL_DIR=""

cleanup_on_exit() {
  local status=$?
  trap - EXIT
  if [[ -n "$ACTIVE_BACKUP_DIR" && -d "$ACTIVE_BACKUP_DIR" &&         -n "$ACTIVE_FINAL_DIR" ]]; then
    if [[ ! -d "$ACTIVE_FINAL_DIR" ]]; then
      mv -- "$ACTIVE_BACKUP_DIR" "$ACTIVE_FINAL_DIR" 2>/dev/null || true
    else
      rm -rf -- "$ACTIVE_BACKUP_DIR"
    fi
  fi
  if [[ -n "$ACTIVE_STAGE_DIR" && -d "$ACTIVE_STAGE_DIR" ]]; then
    rm -rf -- "$ACTIVE_STAGE_DIR"
  fi
  exit "$status"
}
trap cleanup_on_exit EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

validate_bool() {
  local name="$1" value="$2"
  if [[ "$value" != "0" && "$value" != "1" ]]; then
    printf '  %b[FAIL]%b %s must be 0 or 1 (got %s)\n' "$RED" "$NC" "$name" "$value" >&2
    return 1
  fi
}

acquire_run_lock() {
  [[ "$DISABLE_RUN_LOCK" == "1" ]] && return 0
  mkdir -p -- "$OUTPUT_ROOT"
  local lock_file="${OUTPUT_ROOT}/.moran_pipeline.lock"
  exec 9>"$lock_file"
  if ! flock -n 9; then
    printf '%b[FAIL]%b Another Moran pipeline process is using %s\n' "$RED" "$NC" "$OUTPUT_ROOT" >&2
    return 1
  fi
}

preflight_check() {
  local errors=0
  echo "Pre-flight checks:"

  if [[ ! -f "$AUTOCORR_PY" ]]; then
    echo -e "  ${RED}[FAIL]${NC} Missing script: $AUTOCORR_PY"
    errors=$((errors + 1))
  else
    echo -e "  ${GREEN}[ OK ]${NC} Found: $AUTOCORR_PY"
  fi

  local boot_module="$(dirname "$AUTOCORR_PY")/spatial_block_bootstrap.py"
  if [[ "$N_SPATIAL_BOOTSTRAP" -gt 0 && ! -f "$boot_module" ]]; then
    echo -e "  ${RED}[FAIL]${NC} Missing bootstrap module: $boot_module"
    errors=$((errors + 1))
  fi

  if ! command -v sha256sum >/dev/null 2>&1 || ! command -v stat >/dev/null 2>&1; then
    echo -e "  ${RED}[FAIL]${NC} GNU sha256sum and stat are required for cache fingerprints"
    errors=$((errors + 1))
  fi

  if [[ "$DISABLE_RUN_LOCK" != "1" ]] && ! command -v flock >/dev/null 2>&1; then
    echo -e "  ${RED}[FAIL]${NC} flock is required unless DISABLE_RUN_LOCK=1"
    errors=$((errors + 1))
  fi

  for _bool_name in USE_SIGNIFICANCE_MASK USE_FDR BOOTSTRAP_RANDOM_GRID_OFFSET \
      BOOTSTRAP_INCLUDE_PARTIAL_EDGE_BLOCKS BOOTSTRAP_SAVE_DISTRIBUTION FORCE \
      FAIL_FAST DISABLE_RUN_LOCK; do
    if ! validate_bool "$_bool_name" "${!_bool_name}"; then
      errors=$((errors + 1))
    fi
  done

  if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
    echo -e "  ${RED}[FAIL]${NC} Python not found: $PYTHON_BIN"
    errors=$((errors + 1))
  else
    echo -e "  ${GREEN}[ OK ]${NC} Python: $($PYTHON_BIN --version 2>&1)"
    if ! "$PYTHON_BIN" -c 'import numpy, scipy, rasterio' >/dev/null 2>&1; then
      echo -e "  ${RED}[FAIL]${NC} Python packages numpy/scipy/rasterio are required"
      errors=$((errors + 1))
    fi
  fi

  if [[ "$AGGREGATE_EDGE_MODE" != "partial" && "$AGGREGATE_EDGE_MODE" != "trim" ]]; then
    echo -e "  ${RED}[FAIL]${NC} AGGREGATE_EDGE_MODE must be partial or trim"
    errors=$((errors + 1))
  fi
  if ! awk -v x="$AGGREGATE_MIN_VALID_FRACTION" 'BEGIN {exit !(x >= 0 && x <= 1)}'; then
    echo -e "  ${RED}[FAIL]${NC} AGGREGATE_MIN_VALID_FRACTION must be between 0 and 1"
    errors=$((errors + 1))
  fi
  if (( MIN_MORAN_PIXELS < 4 )); then
    echo -e "  ${RED}[FAIL]${NC} MIN_MORAN_PIXELS must be >=4"
    errors=$((errors + 1))
  fi
  if [[ "$SMALL_MORAN_SAMPLE_ACTION" != "point-only" && "$SMALL_MORAN_SAMPLE_ACTION" != "error" ]]; then
    echo -e "  ${RED}[FAIL]${NC} SMALL_MORAN_SAMPLE_ACTION must be point-only or error"
    errors=$((errors + 1))
  fi
  if (( MIN_LOCAL_PIXELS < 2 )); then
    echo -e "  ${RED}[FAIL]${NC} MIN_LOCAL_PIXELS must be >=2"
    errors=$((errors + 1))
  fi
  if [[ "$N_SPATIAL_BOOTSTRAP" -gt 0 ]]; then
    if (( N_SPATIAL_BOOTSTRAP < 2 )); then
      echo -e "  ${RED}[FAIL]${NC} N_SPATIAL_BOOTSTRAP must be 0 or >=2"
      errors=$((errors + 1))
    fi
    if (( BOOTSTRAP_WORKERS < 1 || BOOTSTRAP_BLOCK_WORKERS < 1 )); then
      echo -e "  ${RED}[FAIL]${NC} Bootstrap worker counts must be >=1"
      errors=$((errors + 1))
    fi
    if (( BOOTSTRAP_MIN_USABLE_BLOCKS < 2 )); then
      echo -e "  ${RED}[FAIL]${NC} BOOTSTRAP_MIN_USABLE_BLOCKS must be >=2"
      errors=$((errors + 1))
    fi
    if [[ -n "$BOOTSTRAP_BLOCK_CELLS" ]] && (( BOOTSTRAP_BLOCK_CELLS < 2 )); then
      echo -e "  ${RED}[FAIL]${NC} BOOTSTRAP_BLOCK_CELLS must be >=2"
      errors=$((errors + 1))
    fi
    if [[ "$BOOTSTRAP_EDGE_MODE" != "balanced" && "$BOOTSTRAP_EDGE_MODE" != "drop" && "$BOOTSTRAP_EDGE_MODE" != "partial" ]]; then
      echo -e "  ${RED}[FAIL]${NC} BOOTSTRAP_EDGE_MODE must be balanced, drop, or partial"
      errors=$((errors + 1))
    fi
    if [[ "$BOOTSTRAP_RANDOM_GRID_OFFSET" == "1" && "$BOOTSTRAP_EDGE_MODE" == "balanced" && "$BOOTSTRAP_INCLUDE_PARTIAL_EDGE_BLOCKS" != "1" ]]; then
      echo -e "  ${RED}[FAIL]${NC} Random grid offset requires BOOTSTRAP_EDGE_MODE=drop or partial"
      errors=$((errors + 1))
    fi
  fi

  local seen_scales="|"
  local raw_scale clean_scale
  for raw_scale in "${SCALES_ARR[@]}"; do
    clean_scale="$(echo "$raw_scale" | tr -d '[:space:]')"
    if [[ -z "$clean_scale" || -z "${LABEL_FACTOR[$clean_scale]:-}" ]]; then
      echo -e "  ${RED}[FAIL]${NC} Unknown scale in ONLY/SCALES_TO_RUN: '$clean_scale'"
      errors=$((errors + 1))
    elif [[ "$seen_scales" == *"|${clean_scale}|"* ]]; then
      echo -e "  ${RED}[FAIL]${NC} Duplicate scale requested: $clean_scale"
      errors=$((errors + 1))
    else
      seen_scales="${seen_scales}${clean_scale}|"
    fi
  done

  if [[ "$errors" -gt 0 ]]; then
    exit 1
  fi
}

csv_get() {
  local key="$1"
  local csv="$2"
  awk -F',' -v k="$key" '$1==k {print $2; exit}' "$csv"
}

metric_is_finite() {
  local value="$1"
  [[ "$value" =~ ^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$ ]]
}

path_signature() {
  local path="$1"
  if [[ -e "$path" ]]; then
    stat -c '%n|%s|%y' "$path"
  else
    printf '%s|MISSING\n' "$path"
  fi
}

make_run_config_id() {
  local period_tag="$1" response_key="$2" slope_tif="$3" mk_tif="$4"
  local scale_label="$5" agg_factor="$6" block_cells="$7"
  {
    printf 'period=%s\nresponse=%s\nscale=%s\naggregate_factor=%s\n' \
      "$period_tag" "$response_key" "$scale_label" "$agg_factor"
    printf 'slope_scale=%s\naggregate_edge=%s\naggregate_min_valid_fraction=%s\n' \
      "$SLOPE_INPUT_SCALE" "$AGGREGATE_EDGE_MODE" "$AGGREGATE_MIN_VALID_FRACTION"
    printf 'neighbor=%s\nkernel_radius=%s\nmin_moran_pixels=%s\nsmall_sample_action=%s\n' \
      "$NEIGHBOR" "$KERNEL_RADIUS" "$MIN_MORAN_PIXELS" "$SMALL_MORAN_SAMPLE_ACTION"
    printf 'permutations=%s\nseed=%s\nuse_significance=%s\nuse_fdr=%s\np_threshold=%s\n' \
      "$N_PERMUTATIONS" "$SEED" "$USE_SIGNIFICANCE_MASK" "$USE_FDR" "$P_THRESHOLD"
    printf 'mk_band=%s\nmk_scale=%s\nforest_mask=%s\n' \
      "$MK_PVALUE_BAND" "$MK_PVALUE_SCALE" "$FOREST_MASK"
    printf 'bootstrap_n=%s\nblock_cells=%s\nci_level=%s\nmin_valid_pixels=%s\nmin_valid_fraction=%s\n' \
      "$N_SPATIAL_BOOTSTRAP" "$block_cells" "$BOOTSTRAP_CI_LEVEL" \
      "$BOOTSTRAP_MIN_VALID_PIXELS" "$BOOTSTRAP_MIN_VALID_FRACTION"
    printf 'random_offset=%s\nedge_mode=%s\npartial_alias=%s\nmax_partition_diff=%s\nmax_valid_cv=%s\nmin_usable_blocks=%s\nsave_distribution=%s\n' \
      "$BOOTSTRAP_RANDOM_GRID_OFFSET" "$BOOTSTRAP_EDGE_MODE" \
      "$BOOTSTRAP_INCLUDE_PARTIAL_EDGE_BLOCKS" "$BOOTSTRAP_MAX_PARTITION_DIFFERENCE" \
      "$BOOTSTRAP_MAX_VALID_PIXEL_CV" "$BOOTSTRAP_MIN_USABLE_BLOCKS" \
      "$BOOTSTRAP_SAVE_DISTRIBUTION"
    path_signature "$slope_tif"
    [[ "$USE_SIGNIFICANCE_MASK" == "1" ]] && path_signature "$mk_tif"
    [[ -n "$FOREST_MASK" ]] && path_signature "$FOREST_MASK"
    path_signature "$AUTOCORR_PY"
    if [[ "$N_SPATIAL_BOOTSTRAP" -gt 0 ]]; then
      path_signature "$(dirname "$AUTOCORR_PY")/spatial_block_bootstrap.py"
    fi
  } | sha256sum | awk '{print $1}'
}

append_summary_row() {
  # Args: summary_csv period response scale factor source_metric_csv final_metric_path
  local summary_csv="$1"
  local period_tag="$2"
  local response_key="$3"
  local scale_label="$4"
  local agg_factor="$5"
  local source_csv="$6"
  local final_csv_path="${7:-$source_csv}"

  local run_config_id I n_pixels W_total z_rand p_rand mask_mode result_status inference_available
  local boot_se boot_ci_low boot_ci_high boot_blocks boot_block_rows boot_block_cols
  local boot_bias boot_bias_partition boot_partition_diff boot_coverage boot_reliable
  local boot_reason boot_method boot_skipped
  run_config_id="$(csv_get "run_config_id" "$source_csv")"
  I="$(csv_get "morans_I" "$source_csv")"
  n_pixels="$(csv_get "n_pixels" "$source_csv")"
  W_total="$(csv_get "W_total" "$source_csv")"
  z_rand="$(csv_get "z_randomization" "$source_csv")"
  p_rand="$(csv_get "p_randomization_two_tailed" "$source_csv")"
  mask_mode="$(csv_get "moran_mask_mode" "$source_csv")"
  result_status="$(csv_get "moran_result_status" "$source_csv")"
  inference_available="$(csv_get "moran_inference_available" "$source_csv")"
  boot_se="$(csv_get "bootstrap_se" "$source_csv")"
  boot_ci_low="$(csv_get "bootstrap_ci_low" "$source_csv")"
  boot_ci_high="$(csv_get "bootstrap_ci_high" "$source_csv")"
  boot_blocks="$(csv_get "bootstrap_usable_blocks" "$source_csv")"
  boot_block_rows="$(csv_get "bootstrap_block_rows_cells" "$source_csv")"
  boot_block_cols="$(csv_get "bootstrap_block_cols_cells" "$source_csv")"
  boot_bias="$(csv_get "bootstrap_bias_vs_full_I" "$source_csv")"
  boot_bias_partition="$(csv_get "bootstrap_bias_vs_partition_I" "$source_csv")"
  boot_partition_diff="$(csv_get "bootstrap_partition_difference_vs_full_I" "$source_csv")"
  boot_coverage="$(csv_get "bootstrap_valid_pixel_coverage_fraction" "$source_csv")"
  boot_reliable="$(csv_get "bootstrap_ci_reliable" "$source_csv")"
  boot_reason="$(csv_get "bootstrap_ci_reliability_reasons" "$source_csv")"
  boot_method="$(csv_get "bootstrap_ci_method" "$source_csv")"
  [[ -z "$boot_method" ]] && boot_method="$(csv_get "bootstrap_ci_method_default" "$source_csv")"
  boot_skipped="$(csv_get "bootstrap_skipped_reason" "$source_csv")"

  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$period_tag" "$response_key" "$scale_label" "$agg_factor" "$run_config_id" \
    "$I" "$n_pixels" "$W_total" "$z_rand" "$p_rand" "$mask_mode" \
    "$result_status" "$inference_available" "$boot_se" "$boot_ci_low" "$boot_ci_high" \
    "$boot_blocks" "$boot_block_rows" "$boot_block_cols" "$boot_bias" \
    "$boot_bias_partition" "$boot_partition_diff" "$boot_coverage" "$boot_reliable" \
    "$boot_reason" "$boot_method" "$boot_skipped" "$final_csv_path" >> "$summary_csv"
}

run_one() {
  # Args: period response slope mk scale factor staging_dir final_dir summary_csv
  local period_tag="$1"
  local response_key="$2"
  local slope_tif="$3"
  local mk_tif="$4"
  local scale_label="$5"
  local agg_factor="$6"
  local out_dir="$7"
  local final_out_dir="$8"
  local summary_csv="$9"

  local prefix="${response_key}_${period_tag}_${scale_label}"
  local moran_csv="${out_dir}/${prefix}_Moran.csv"
  local final_moran_csv="${final_out_dir}/${prefix}_Moran.csv"
  local boot_npy="${out_dir}/${prefix}_Moran_block_bootstrap.npy"
  local python_log="${out_dir}/${prefix}_log.txt"
  local driver_log="${out_dir}/${prefix}_driver.log"

  if [[ ! -f "$slope_tif" ]]; then
    echo -e "  ${RED}[MISS]${NC} Slope missing: $slope_tif"
    return 1
  fi
  if [[ "$USE_SIGNIFICANCE_MASK" == "1" && ! -f "$mk_tif" ]]; then
    echo -e "  ${RED}[MISS]${NC} MK/FDR mask missing: $mk_tif"
    return 1
  fi

  local block_cells=""
  if [[ "$N_SPATIAL_BOOTSTRAP" -gt 0 ]]; then
    block_cells="$(get_bootstrap_block_cells "$scale_label")"
    if [[ -z "$block_cells" ]]; then
      echo -e "  ${RED}[FAIL]${NC} No bootstrap block size for scale $scale_label"
      echo "         Set BOOTSTRAP_BLOCK_CELLS or BOOTSTRAP_BLOCK_CELLS_BY_SCALE."
      return 1
    fi
  fi
  local run_config_id
  run_config_id="$(make_run_config_id "$period_tag" "$response_key" "$slope_tif" \
    "$mk_tif" "$scale_label" "$agg_factor" "$block_cells")"

  if [[ -f "$moran_csv" && "$FORCE" != "1" ]]; then
    local old_I
    old_I="$(csv_get "morans_I" "$moran_csv")"
    local old_config_id
    old_config_id="$(csv_get "run_config_id" "$moran_csv")"
    if ! metric_is_finite "$old_I"; then
      echo -e "  ${YELLOW}[WARN]${NC} Existing Moran.csv is invalid; rerunning $prefix"
    elif [[ -z "$old_config_id" || "$old_config_id" != "$run_config_id" ]]; then
      echo -e "  ${YELLOW}[STALE]${NC} $prefix configuration/input fingerprint changed; rerunning"
    else
      echo "  [SKIP] $prefix (matching run_config_id; set FORCE=1 to rerun)"
      append_summary_row "$summary_csv" "$period_tag" "$response_key" "$scale_label" \
        "$agg_factor" "$moran_csv" "$final_moran_csv"
      return 0
    fi
  fi

  echo "  [RUN]  $prefix  aggregate_factor=$agg_factor"
  local t0 t1 dt
  t0=$(date +%s)

  local cmd=(
    "$PYTHON_BIN" "$AUTOCORR_PY"
    --slope-input "$slope_tif"
    --slope-input-scale "$SLOPE_INPUT_SCALE"
    --aggregate-factor "$agg_factor"
    --aggregate-edge-mode "$AGGREGATE_EDGE_MODE"
    --aggregate-min-valid-fraction "$AGGREGATE_MIN_VALID_FRACTION"
    --neighbor "$NEIGHBOR"
    --kernel-radius "$KERNEL_RADIUS"
    --min-moran-pixels "$MIN_MORAN_PIXELS"
    --small-moran-sample-action "$SMALL_MORAN_SAMPLE_ACTION"
    --min-local-pixels "$MIN_LOCAL_PIXELS"
    --permutations "$N_PERMUTATIONS"
    --seed "$SEED"
    --global-only
    --output-dir "$out_dir"
    --prefix "$prefix"
    --run-config-id "$run_config_id"
    --log-level "$LOG_LEVEL"
  )

  if [[ "$N_SPATIAL_BOOTSTRAP" -gt 0 ]]; then
    cmd+=(
      --spatial-bootstrap "$N_SPATIAL_BOOTSTRAP"
      --bootstrap-block-size-cells "$block_cells"
      --bootstrap-workers "$BOOTSTRAP_WORKERS"
      --bootstrap-block-workers "$BOOTSTRAP_BLOCK_WORKERS"
      --bootstrap-ci-level "$BOOTSTRAP_CI_LEVEL"
      --bootstrap-seed "$SEED"
      --bootstrap-min-valid-pixels "$BOOTSTRAP_MIN_VALID_PIXELS"
      --bootstrap-min-valid-fraction "$BOOTSTRAP_MIN_VALID_FRACTION"
      --bootstrap-edge-mode "$BOOTSTRAP_EDGE_MODE"
      --bootstrap-max-partition-difference "$BOOTSTRAP_MAX_PARTITION_DIFFERENCE"
      --bootstrap-max-valid-pixel-cv "$BOOTSTRAP_MAX_VALID_PIXEL_CV"
      --bootstrap-min-usable-blocks "$BOOTSTRAP_MIN_USABLE_BLOCKS"
    )
    [[ "$BOOTSTRAP_RANDOM_GRID_OFFSET" == "1" ]] && cmd+=( --bootstrap-random-grid-offset )
    [[ "$BOOTSTRAP_INCLUDE_PARTIAL_EDGE_BLOCKS" == "1" ]] && cmd+=( --bootstrap-include-partial-edge-blocks )
    [[ "$BOOTSTRAP_SAVE_DISTRIBUTION" == "1" ]] && cmd+=( --bootstrap-save-distribution )
  fi
  [[ -n "$FOREST_MASK" ]] && cmd+=( --forest-mask "$FOREST_MASK" )
  if [[ "$USE_SIGNIFICANCE_MASK" == "1" ]]; then
    cmd+=(
      --mk-input "$mk_tif"
      --mk-pvalue-band "$MK_PVALUE_BAND"
      --mk-pvalue-scale "$MK_PVALUE_SCALE"
      --p-threshold "$P_THRESHOLD"
    )
  fi

  # Per-file backup remains useful inside the staging directory, while the
  # period-level transaction guarantees the published directory is unchanged
  # unless every requested task succeeds.
  local backup_csv="${moran_csv}.pre_run_backup"
  local backup_npy="${boot_npy}.pre_run_backup"
  rm -f "$backup_csv" "$backup_npy"
  [[ -f "$moran_csv" ]] && mv "$moran_csv" "$backup_csv"
  [[ -f "$boot_npy" ]] && mv "$boot_npy" "$backup_npy"
  rm -f "$driver_log"

  local cmd_status=0
  set +e
  "${cmd[@]}" > "$driver_log" 2>&1
  cmd_status=$?
  set -e
  t1=$(date +%s)
  dt=$((t1 - t0))

  local I_value=""
  [[ -s "$moran_csv" ]] && I_value="$(csv_get "morans_I" "$moran_csv")"
  if [[ "$cmd_status" -ne 0 || ! -s "$moran_csv" ]] || ! metric_is_finite "$I_value"; then
    rm -f "$moran_csv" "$boot_npy"
    [[ -f "$backup_csv" ]] && mv "$backup_csv" "$moran_csv"
    [[ -f "$backup_npy" ]] && mv "$backup_npy" "$boot_npy"
    echo -e "         ${RED}FAILED${NC} (exit=${cmd_status}); see $driver_log and $python_log"
    tail -n 20 "$driver_log" 2>/dev/null || true
    return 1
  fi

  local result_status
  result_status="$(csv_get "moran_result_status" "$moran_csv")"
  [[ -z "$result_status" ]] && result_status="complete"
  local bootstrap_skipped_reason
  bootstrap_skipped_reason="$(csv_get "bootstrap_skipped_reason" "$moran_csv")"
  if [[ "$N_SPATIAL_BOOTSTRAP" -gt 0 && "$result_status" == "complete" && \
        -z "$bootstrap_skipped_reason" && \
        -z "$(csv_get "bootstrap_n_success" "$moran_csv")" ]]; then
    rm -f "$moran_csv" "$boot_npy"
    [[ -f "$backup_csv" ]] && mv "$backup_csv" "$moran_csv"
    [[ -f "$backup_npy" ]] && mv "$backup_npy" "$boot_npy"
    echo -e "         ${RED}FAILED${NC}; output CSV lacks bootstrap metrics"
    return 1
  fi

  rm -f "$backup_csv" "$backup_npy"
  append_summary_row "$summary_csv" "$period_tag" "$response_key" "$scale_label" \
    "$agg_factor" "$moran_csv" "$final_moran_csv"
  local I_val n_val p_val ci_reliable
  I_val="$(awk -v x="$I_value" 'BEGIN {printf "%.4f", x}')"
  n_val="$(csv_get "n_pixels" "$moran_csv")"
  p_val="$(csv_get "p_randomization_two_tailed" "$moran_csv")"
  ci_reliable="$(csv_get "bootstrap_ci_reliable" "$moran_csv")"
  echo "         done in ${dt}s  Moran I=${I_val}, n=${n_val}, status=${result_status}, p=${p_val}, CI_reliable=${ci_reliable:-NA}"
}

run_period() {
  local start_year="$1"
  local end_year="$2"
  local period_failures=0
  local period_tag="${start_year}_${end_year}"
  local period_dir="${OUTPUT_ROOT}/period_${period_tag}"

  local mask_tag="allvalid"
  if [[ "$USE_SIGNIFICANCE_MASK" == "1" && "$USE_FDR" == "1" ]]; then
    mask_tag="fdr_q${P_THRESHOLD}"
  elif [[ "$USE_SIGNIFICANCE_MASK" == "1" ]]; then
    mask_tag="mkp${P_THRESHOLD}"
  fi
  mask_tag="$(echo "$mask_tag" | tr -d '.')"
  [[ "$N_SPATIAL_BOOTSTRAP" -gt 0 ]] && mask_tag="${mask_tag}_blockboot"

  local final_out_dir="${period_dir}/spatial_moran_multiscale_${mask_tag}"
  local stage_out_dir="${final_out_dir}.staging.$$"
  local backup_out_dir="${final_out_dir}.replace_backup.$$"
  local summary_csv="${stage_out_dir}/Moran_multiscale_summary.csv"

  local sr_trend_dir="${period_dir}/trend/diversity_rich_BASIC"
  local sh_trend_dir="${period_dir}/trend/diversity_shannon_BASIC"
  local sr_prefix="${SR_RESPONSE}_${period_tag}_BASIC"
  local sh_prefix="${SH_RESPONSE}_${period_tag}_BASIC"
  local sr_slope="${sr_trend_dir}/${sr_prefix}_TrendSlope.tif"
  local sh_slope="${sh_trend_dir}/${sh_prefix}_TrendSlope.tif"
  local sr_mk="${sr_trend_dir}/${sr_prefix}_MannKendall.tif"
  local sh_mk="${sh_trend_dir}/${sh_prefix}_MannKendall.tif"
  if [[ "$USE_FDR" == "1" ]]; then
    sr_mk="${sr_trend_dir}/${sr_prefix}_MannKendall_FDR.tif"
    sh_mk="${sh_trend_dir}/${sh_prefix}_MannKendall_FDR.tif"
  fi

  # Build a private staging copy. Published outputs remain untouched until all
  # requested response/scale tasks and the new summary succeed. The EXIT trap
  # removes incomplete staging data and restores a temporarily moved published
  # directory if the process is interrupted during the final rename.
  ACTIVE_STAGE_DIR="$stage_out_dir"
  ACTIVE_BACKUP_DIR="$backup_out_dir"
  ACTIVE_FINAL_DIR="$final_out_dir"
  rm -rf -- "$stage_out_dir" "$backup_out_dir"
  mkdir -p "$stage_out_dir"
  if [[ -d "$final_out_dir" ]]; then
    cp -a "$final_out_dir/." "$stage_out_dir/"
  fi
  echo 'period,response,scale,aggregate_factor,run_config_id,morans_I,n_pixels,W_total,z_randomization,p_randomization_two_tailed,moran_mask_mode,moran_result_status,moran_inference_available,bootstrap_se,bootstrap_ci_low,bootstrap_ci_high,bootstrap_usable_blocks,bootstrap_block_rows_cells,bootstrap_block_cols_cells,bootstrap_bias_vs_full_I,bootstrap_bias_vs_partition_I,bootstrap_partition_difference_vs_full_I,bootstrap_valid_pixel_coverage_fraction,bootstrap_ci_reliable,bootstrap_ci_reliability_reasons,bootstrap_ci_method,bootstrap_skipped_reason,moran_csv' > "$summary_csv"

  echo ""
  echo -e "${BOLD}${CYAN}================================================================${NC}"
  echo -e "${BOLD}${CYAN}Period ${period_tag}: multi-scale Global Moran's I only${NC}"
  echo -e "${BOLD}${CYAN}================================================================${NC}"
  echo "Published output: $final_out_dir"
  echo "Staging output:   $stage_out_dir"
  echo "Scales:           ${SCALES_ARR[*]}"
  echo "Neighbor:         ${NEIGHBOR}, radius=${KERNEL_RADIUS}"
  echo "Slope scale:      ${SLOPE_INPUT_SCALE}"
  echo "Agg edges:        ${AGGREGATE_EDGE_MODE}"
  echo "Agg valid frac:   ${AGGREGATE_MIN_VALID_FRACTION}"
  echo "Moran minimum:    ${MIN_MORAN_PIXELS} (${SMALL_MORAN_SAMPLE_ACTION})"
  if [[ "$N_SPATIAL_BOOTSTRAP" -gt 0 ]]; then
    echo "Block boot:       n=${N_SPATIAL_BOOTSTRAP}, workers=${BOOTSTRAP_WORKERS}, block-workers=${BOOTSTRAP_BLOCK_WORKERS}, min-blocks=${BOOTSTRAP_MIN_USABLE_BLOCKS}, edge=${BOOTSTRAP_EDGE_MODE}"
    echo "Block map:        ${BOOTSTRAP_BLOCK_CELLS_BY_SCALE:-global=${BOOTSTRAP_BLOCK_CELLS:-UNSET}}"
  else
    echo "Block boot:       disabled (set N_SPATIAL_BOOTSTRAP=2000 to match the publication default)"
  fi
  if [[ "$USE_SIGNIFICANCE_MASK" == "1" && "$USE_FDR" == "1" ]]; then
    echo "Moran mask:       FDR q <= ${P_THRESHOLD}"
  elif [[ "$USE_SIGNIFICANCE_MASK" == "1" ]]; then
    echo "Moran mask:       raw MK p <= ${P_THRESHOLD}"
  else
    echo "Moran mask:       all valid slope pixels; Mkp/q masks are NOT used"
  fi
  echo ""

  local response_keys=("SR_Rich_tree" "Shannon_Shannon_wiener")
  for response_key in "${response_keys[@]}"; do
    local slope_tif mk_tif
    case "$response_key" in
      SR_Rich_tree) slope_tif="$sr_slope"; mk_tif="$sr_mk" ;;
      Shannon_Shannon_wiener) slope_tif="$sh_slope"; mk_tif="$sh_mk" ;;
      *) period_failures=$((period_failures + 1)); break ;;
    esac

    echo "=== Response: $response_key ==="
    for scale_label_raw in "${SCALES_ARR[@]}"; do
      local scale_label
      scale_label="$(echo "$scale_label_raw" | tr -d '[:space:]')"
      local factor="${LABEL_FACTOR[$scale_label]}"
      if ! run_one "$period_tag" "$response_key" "$slope_tif" "$mk_tif" \
          "$scale_label" "$factor" "$stage_out_dir" "$final_out_dir" "$summary_csv"; then
        echo -e "  ${YELLOW}[WARN]${NC} Failed: $response_key $period_tag $scale_label"
        period_failures=$((period_failures + 1))
        if [[ "$FAIL_FAST" == "1" ]]; then
          break 2
        fi
      fi
    done
    echo ""
  done

  if (( period_failures > 0 )); then
    rm -rf -- "$stage_out_dir"
    ACTIVE_STAGE_DIR=""
    ACTIVE_BACKUP_DIR=""
    ACTIVE_FINAL_DIR=""
    echo -e "${RED}[PERIOD FAILED]${NC} $period_tag had ${period_failures} failed task(s); published directory was not changed"
    return 1
  fi

  # Transactional publish: keep the prior complete directory until the staging
  # directory is ready, and restore it if the final rename fails.
  if [[ -d "$final_out_dir" ]]; then
    mv -- "$final_out_dir" "$backup_out_dir"
  fi
  if ! mv -- "$stage_out_dir" "$final_out_dir"; then
    [[ -d "$backup_out_dir" ]] && mv -- "$backup_out_dir" "$final_out_dir"
    ACTIVE_STAGE_DIR=""
    ACTIVE_BACKUP_DIR=""
    ACTIVE_FINAL_DIR=""
    echo -e "${RED}[PUBLISH FAILED]${NC} Could not publish $final_out_dir; previous directory restored"
    return 1
  fi
  rm -rf -- "$backup_out_dir"
  ACTIVE_STAGE_DIR=""
  ACTIVE_BACKUP_DIR=""
  ACTIVE_FINAL_DIR=""
  echo -e "${GREEN}[DONE]${NC} Period $period_tag outputs published transactionally: $final_out_dir"
  return 0
}

main() {
  local t0 t1 total
  local failed_periods=0
  t0=$(date +%s)

  echo -e "${BOLD}${CYAN}Multi-scale Global Moran's I only pipeline${NC}"
  preflight_check
  acquire_run_lock

  IFS=',' read -ra PERIOD_ARR <<< "$PERIODS"
  for period in "${PERIOD_ARR[@]}"; do
    period="$(echo "$period" | tr -d '[:space:]')"
    if [[ ! "$period" =~ ^[0-9]{4}-[0-9]{4}$ ]]; then
      echo -e "${RED}[FAIL]${NC} Bad PERIODS entry: $period. Use YYYY-YYYY."
      exit 1
    fi
    local start_year="${period%-*}"
    local end_year="${period#*-}"
    if ! run_period "$start_year" "$end_year"; then
      failed_periods=$((failed_periods + 1))
      if [[ "$FAIL_FAST" == "1" ]]; then
        exit 1
      fi
    fi
  done

  t1=$(date +%s)
  total=$((t1 - t0))
  echo ""
  if (( failed_periods > 0 )); then
    echo -e "${RED}[FAILED]${NC} ${failed_periods} period(s) had failed tasks; elapsed: ${total}s"
    exit 1
  fi
  echo -e "${GREEN}[ALL DONE]${NC} Total elapsed: ${total}s"
}

main "$@"
