#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
case "${1:-}" in
  1|01) FOLDER="01_TRACE_identifiable_mixture_contamination" ;;
  2|02) FOLDER="02_TRACE_support_swap_contamination" ;;
  3|03) FOLDER="03_TRACE_mosaic_contamination" ;;
  4|04) FOLDER="04_TRACE_complementary_sources_fixed_noise" ;;
  *) echo "Usage: bash run_setting.sh 01|02|03|04 [trace|all]" >&2; exit 2 ;;
esac
export SIM_MODE="${2:-${SIM_MODE:-trace}}"
case "$SIM_MODE" in trace|all) ;; *) echo "SIM_MODE must be trace or all." >&2; exit 2 ;; esac
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export SIM_NREP="${SIM_NREP:-50}" SIM_NCORES="${SIM_NCORES:-48}" SIM_RESUME="${SIM_RESUME:-TRUE}"
export SIM_SOURCE_N="${SIM_SOURCE_N:-500}"
BASE="${SIM_RESULT_ROOT:-${SIM_OUTPUT_BASE:-$ROOT/results_T0_corrected_$SIM_MODE}/$FOLDER}"
export SIM_RESULT_ROOT="$BASE/trace_refit"
export SIM_COMPARISON_ROOT="$BASE/comparison"
export BORROW_CORE_PATH="$ROOT/_common/borrow_functions_core.R"
exec "${RSCRIPT:-Rscript}" "$ROOT/$FOLDER/run_analysis_50rep.R"
