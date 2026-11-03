#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
MODE="${1:-figures}"
if [[ $# -gt 0 ]]; then shift; fi
R="${RSCRIPT:-Rscript}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
case "$MODE" in
  figures)
    exec "$R" "$ROOT/scripts/reproduce_results.R" "$@"
    ;;
  validate)
    exec bash "$ROOT/scripts/validate.sh" "$@"
    ;;
  trace|all)
    export SIM_NREP="${SIM_NREP:-50}" SIM_NCORES="${CORES:-48}"
    export SIM_OUTPUT_BASE="${SIM_OUTPUT_BASE:-$ROOT/runs/simulation_$MODE}"
    if [[ "$MODE" == all ]]; then export SIM_FRESH_INPUTS=TRUE; fi
    bash "$ROOT/simulation/run_all_50rep.sh" "$MODE"
    REAL_OUT="${REAL_OUTPUT_DIR:-$ROOT/runs/real_data_$MODE}"
    if [[ "$MODE" == trace ]]; then
      "$R" "$ROOT/real_data/run_real_data.R" --reps "${REAL_NREP:-50}" \
        --workers "${CORES:-48}" --out "$REAL_OUT" "$@"
    else
      "$R" "$ROOT/real_data/full_analysis/run_real_data.R" \
        --data-dir "$ROOT/real_data/data/processed/gtex" \
        --reps "${REAL_NREP:-50}" --workers "${CORES:-48}" --out "$REAL_OUT" "$@"
    fi
    if [[ "$SIM_NREP" != 50 || "${REAL_NREP:-50}" != 50 ]]; then
      echo 'Reduced run completed; publication-only figures require 50 replications.'
      exit 0
    fi
    "$R" "$ROOT/scripts/reproduce_results.R" --sim-root "$SIM_OUTPUT_BASE" \
      --real-root "$REAL_OUT" --out "$ROOT/runs/derived_$MODE"
    ;;
  *)
    echo "Usage: bash run_analysis.sh figures|validate|trace|all" >&2
    exit 2
    ;;
esac
