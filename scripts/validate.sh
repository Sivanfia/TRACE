#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
R="${RSCRIPT:-Rscript}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
"$R" "$ROOT/scripts/validate_package.R"
"$R" "$ROOT/simulation/tests/test_mode.R"
"$R" "$ROOT/simulation/tests/test_merge.R"
"$R" "$ROOT/simulation/validate_T0.R" "$ROOT/simulation"
"$R" "$ROOT/real_data/tests/validate_T0.R" "$ROOT/real_data"
"$R" "$ROOT/real_data/tests/test_cached_inputs.R" "$ROOT/real_data"
"$R" "$ROOT/real_data/tests/test_merge.R" "$ROOT/real_data"
if [[ "${FULL_VALIDATION:-0}" == 1 ]]; then
  for SETTING in 01 03 04; do "$R" "$ROOT/simulation/tests/test_pipeline.R" "$SETTING"; done
  "$R" "$ROOT/simulation/tests/test_pipeline.R" 02 I
  "$R" "$ROOT/simulation/tests/test_pipeline.R" 02 III
  "$R" "$ROOT/simulation/tests/test_full_methods.R"
  "$R" "$ROOT/real_data/tests/test_full_analysis.R"
  CHECK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/trace_submission_check.XXXXXX")"
  trap 'rm -rf "$CHECK_DIR"' EXIT
  "$R" "$ROOT/real_data/tests/test_pipeline.R" "$ROOT/real_data" "$CHECK_DIR/results"
fi
echo "Validation completed. Rebuild saved-result figures with: bash run_analysis.sh figures"
