#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
Rscript tests/validate_T0.R .
Rscript tests/test_cached_inputs.R .
Rscript tests/test_merge.R .
validation_dir=$(mktemp -d "${TMPDIR:-/tmp}/trace_real_T0_validation.XXXXXX")
trap 'rm -rf "$validation_dir"' EXIT
Rscript tests/test_pipeline.R . "$validation_dir/results"
