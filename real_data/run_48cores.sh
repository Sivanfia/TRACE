#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
Rscript run_real_data.R --reps "${REAL_NREP:-50}" --workers "${WORKERS:-48}" --out "${REAL_OUTPUT_DIR:-results_T0_corrected}" "$@"
