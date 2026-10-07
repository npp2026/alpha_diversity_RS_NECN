#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT="${RSCRIPT:-Rscript}"
RESULT_DIR="${1:-${RESULT_DIR:-}}"
: "${RESULT_DIR:?Pass the scenario output directory as the first argument or set RESULT_DIR}"
N_BOOT="${N_BOOT:-2000}"
BLOCK_KM="${BLOCK_KM:-50}"
"$RSCRIPT" "$SCRIPT_DIR/postprocess/Check_postprocess_inputs.R" "$RESULT_DIR"
"$RSCRIPT" "$SCRIPT_DIR/postprocess/Compute_fig23_stats.R" "$RESULT_DIR" "$N_BOOT" "$BLOCK_KM"
"$RSCRIPT" "$SCRIPT_DIR/postprocess/Make_consolidated_Table_S1.R" "$RESULT_DIR" "$N_BOOT" "$BLOCK_KM"
printf '%s\n' "Done. Statistical output mapping:" \
  "  Table_S1_consolidated_block_removal.csv -> current Table S3" \
  "  Table_S1_consolidated_substitution.csv -> current Table S4"
