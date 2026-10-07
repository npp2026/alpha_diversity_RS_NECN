#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT="${RSCRIPT:-Rscript}"
: "${DATA_CSV:?Set DATA_CSV to the matched plot CSV}"
: "${OUTPUT_DIR:?Set OUTPUT_DIR to a fresh analysis directory}"
SHARED_KNNDM_FOLD_CSV="${SHARED_KNNDM_FOLD_CSV:-}"

[[ -f "$DATA_CSV" ]] || { echo "ERROR: Missing DATA_CSV: $DATA_CSV" >&2; exit 1; }
if [[ -n "$SHARED_KNNDM_FOLD_CSV" && ! -f "$SHARED_KNNDM_FOLD_CSV" ]]; then
  echo "ERROR: Missing SHARED_KNNDM_FOLD_CSV: $SHARED_KNNDM_FOLD_CSV" >&2
  exit 1
fi
: "${SHARED_KNNDM_FOLD_CSV:?Set the shared kNNDM fold CSV}"
"$RSCRIPT" "$SCRIPT_DIR/run_nested_MS.R" "$DATA_CSV" "$OUTPUT_DIR" "$SHARED_KNNDM_FOLD_CSV"
