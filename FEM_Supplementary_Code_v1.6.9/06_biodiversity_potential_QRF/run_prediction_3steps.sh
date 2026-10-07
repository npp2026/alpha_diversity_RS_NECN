#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec "${RSCRIPT:-Rscript}" "$HERE/run_potential_workflow.R" "$@"
