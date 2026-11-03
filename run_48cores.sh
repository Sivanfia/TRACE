#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export CORES=48
exec bash "$ROOT/run_analysis.sh" "${1:-all}"
