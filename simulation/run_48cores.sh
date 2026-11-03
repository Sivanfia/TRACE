#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export SIM_NCORES=48
exec bash "$ROOT/run_all_50rep.sh" "${1:-trace}"
