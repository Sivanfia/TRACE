#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${SIM_RESULT_ROOT:-}" ]]; then
  echo "For all four settings, use SIM_OUTPUT_BASE rather than SIM_RESULT_ROOT." >&2
  exit 2
fi
export SIM_MODE="${1:-${SIM_MODE:-trace}}"
for SETTING in 01 02 03 04; do
  echo "Starting Setting $SETTING, mode=$SIM_MODE"
  bash "$ROOT/run_setting.sh" "$SETTING" "$SIM_MODE"
done
BASE="${SIM_OUTPUT_BASE:-$ROOT/results_T0_corrected_$SIM_MODE}"
mkdir -p "$BASE/paper_figures"
for FILE in "$BASE"/0*_TRACE_*/comparison/paper_figures/*_main_f1.pdf "$BASE"/0*_TRACE_*/comparison/paper_figures/*_main_f1.png; do
  if [[ -f "$FILE" ]]; then cp "$FILE" "$BASE/paper_figures/"; fi
done
echo "Combined paper figures: $BASE/paper_figures"
