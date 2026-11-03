#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
"${RSCRIPT:-Rscript}" "$ROOT/tests/test_mode.R"
"${RSCRIPT:-Rscript}" "$ROOT/tests/test_merge.R"
"${RSCRIPT:-Rscript}" "$ROOT/validate_T0.R" "$ROOT"
for SETTING in 01 03 04; do
  "${RSCRIPT:-Rscript}" "$ROOT/tests/test_pipeline.R" "$SETTING"
done
"${RSCRIPT:-Rscript}" "$ROOT/tests/test_pipeline.R" 02 I
"${RSCRIPT:-Rscript}" "$ROOT/tests/test_pipeline.R" 02 III
